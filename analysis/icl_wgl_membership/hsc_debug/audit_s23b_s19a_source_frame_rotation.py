#!/usr/bin/env python3
"""Replay the fixed S19A profile with source-frame angles as the only change."""
import argparse
from collections import defaultdict
import csv
import json
from pathlib import Path

import numpy as np

from analyze_s23b_wl_fixed_common_support_spatial_covariance import (
    DEFAULT_PLAN_DIR, DEFAULT_RESULTS, EXPECTED_FIXED_PROFILE_SUMMARY_SHA256, EXPECTED_NZ_SHA256,
    FIXED_PROFILE_DIR, NZ_PATH, RADIAL_EDGES_MPC, SUPPORT_DIR, Cosmology, cache_geometry_and_efficiency,
    cached_stack, fixed_sample_ids, group_indices, load_approved_pairs, load_hsc_y3_nz_sacc,
    profile_values, select_rows, sha,
)
from ogfkit import wl_y1


def source_angles(pairs, cache):
    phi = cache['phi'].copy()
    mask = cache['active']
    phi[mask] = wl_y1.source_frame_angle(*(pairs[key][mask] for key in
        ('lens_ra_deg','lens_dec_deg','source_ra_deg','source_dec_deg')))
    offset = .5*np.abs(np.angle(np.exp(2j*(phi[mask]-cache['phi'][mask]))))
    summary = dict(zip(('median','q95','max'),map(float,np.rad2deg(np.quantile(offset,[.5,.95,1]))*3600),strict=True))
    return {**cache,'phi':phi},summary


def run(output):
    if output.exists():
        raise FileExistsError(output)
    fixed_file = FIXED_PROFILE_DIR/'summary.json'
    if sha(fixed_file) != EXPECTED_FIXED_PROFILE_SUMMARY_SHA256 or sha(NZ_PATH) != EXPECTED_NZ_SHA256:
        raise ValueError('known profile or n(z) input changed')
    fixed_profile = json.loads(fixed_file.read_text())['radial_profile_msun_mpc2']
    original = np.array([(r['delta_sigma_mean_msun_mpc2'],r['cross_delta_sigma_mean_msun_mpc2']) for r in fixed_profile])
    lens_ids, selection = fixed_sample_ids(SUPPORT_DIR)
    lp, rp, match, provenance = load_approved_pairs(DEFAULT_PLAN_DIR, DEFAULT_RESULTS)
    controls = defaultdict(list)
    for center,lens in match.items():
        if lens in lens_ids:
            controls[lens].append(center)
    if set(controls) != lens_ids or any(len(items)!=5 for items in controls.values()):
        raise ValueError('matched control ledger differs')
    lens_pairs = select_rows(lp,lens_ids)
    random_pairs = select_rows(rp,{c for items in controls.values() for c in items})
    del lp,rp
    nz_z,nz_bins = load_hsc_y3_nz_sacc(NZ_PATH)
    cosmo = Cosmology(H0=70., Om=.3)
    lc = cache_geometry_and_efficiency(lens_pairs,nz_z,nz_bins,cosmo)
    rc = cache_geometry_and_efficiency(random_pairs,nz_z,nz_bins,cosmo)
    lc_source,l_angle = source_angles(lens_pairs,lc)
    rc_source,r_angle = source_angles(random_pairs,rc)
    li = group_indices(lens_pairs['lens_id'])
    ri = group_indices(random_pairs['lens_id'])
    legacy_values, source_values, rows = [], [], []
    for lens in sorted(lens_ids):
        target = li[lens]
        random = np.concatenate([ri[c] for c in sorted(controls[lens]) if c in ri])
        legacy_l = cached_stack(lens_pairs,lc,target,cosmo)
        legacy_r = cached_stack(random_pairs,rc,random,cosmo)
        source_l = cached_stack(lens_pairs,lc_source,target,cosmo)
        source_r = cached_stack(random_pairs,rc_source,random,cosmo)
        for before,after in ((legacy_l,source_l),(legacy_r,source_r)):
            for a,b in zip(before['profile'],after['profile'],strict=True):
                for field in ('n_pairs','n_lenses','n_sources','responsivity','selection_bias'):
                    if a[field]!=b[field]:
                        raise ValueError(f'rotation changed {field}')
        base,n_controls = profile_values(legacy_l,legacy_r)
        source,_ = profile_values(source_l,source_r)
        if np.any(n_controls<4) or not np.isfinite(base).all() or not np.isfinite(source).all():
            raise ValueError('fixed sample lost matched radial support')
        legacy_values.append(base)
        source_values.append(source)
        for k,(a,b) in enumerate(zip(base,source,strict=True)):
            rows.append(dict(matched_lens_id=lens,radial_bin=k,n_controls=int(n_controls[k]),
                original_t=float(a[0]),original_cross=float(a[1]),source_t=float(b[0]),source_cross=float(b[1])))
    legacy = np.mean(legacy_values,axis=0)
    source = np.mean(source_values,axis=0)
    if not np.allclose(legacy,original,rtol=2e-9,atol=1e-3):
        raise ValueError('baseline did not reproduce the pinned fixed-52 profile')
    profile = [dict(r_min=float(lo),r_max=float(hi),legacy_t=float(a[0]),source_t=float(b[0]),
                    delta_t=float(b[0]-a[0]),legacy_cross=float(a[1]),source_cross=float(b[1]),
                    delta_cross=float(b[1]-a[1])) for lo,hi,a,b in
               zip(RADIAL_EDGES_MPC[:-1],RADIAL_EDGES_MPC[1:],legacy,source,strict=True)]
    output.mkdir(parents=True,exist_ok=False)
    table = output/'per_lens_angle_sensitivity.tsv'
    with table.open('x',newline='') as stream:
        writer = csv.DictWriter(stream,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)
    summary = dict(science_measurement=False, selection=selection, row_provenance=provenance, profile=profile,
                   units='Msun per proper Mpc^2', angle_difference_arcsec={'lens':l_angle,'random':r_angle},
                   representative_only='Rotation-only sensitivity of the original fixed-52 diagnostic with the same S19A catalogue, Y3 n(z), radial selection and calibration.',
                   fixed_profile_sha256=EXPECTED_FIXED_PROFILE_SUMMARY_SHA256, nz_sha256=EXPECTED_NZ_SHA256,
                   output_sha256={table.name:sha(table)}, runner_sha256=sha(Path(__file__)), rotation_sha256=sha(Path(wl_y1.__file__)))
    with (output/'summary.json').open('x') as stream:
        json.dump(summary,stream,indent=2,allow_nan=False)
        stream.write('\n')
    print(json.dumps(profile,indent=2))


if __name__=='__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    run(p.parse_args().output.resolve())
