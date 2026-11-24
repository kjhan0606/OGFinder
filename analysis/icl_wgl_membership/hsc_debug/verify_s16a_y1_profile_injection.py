#!/usr/bin/env python3
"""Run a noiseless Y1 profile injection on the actual fixed pilot geometry."""
import argparse
from collections import defaultdict
import csv
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u

from audit_s16a_pointz_background_support import (
    C_KMS, G_MPC, GEOMETRY, GEOMETRY_SHA, OBJECT_QA, PLAN, RADIAL_EDGES, FlatLambdaCDM,
    finite_point_redshifts, load_inputs, sha, write_table,
)
from ogfkit import wl_photoz, wl_y1

SUPPORT = PLAN.parent/'s16a_pointz_background_support_2026-11-18_01'/'summary.json'
SUPPORT_SHA = '9d69cd60216c593bdce2f937297e17e89b833abdc1a3525ddcb8e140a731046b'
GAP = .1
MSEL = .012  # Injected response, not a measured selection correction.


def beta_from_distances(cosmo, lens_z, source_z):
    """Independent angular-diameter-distance reference for the injection."""
    source_z = np.asarray(source_z)
    beta = np.zeros_like(source_z)
    behind = source_z > lens_z
    # Astropy versions before and after the method rename both occur on hosts.
    dls = cosmo.angular_diameter_distance_z1z2(lens_z[behind], source_z[behind]).value
    beta[behind] = (4*np.pi*G_MPC/C_KMS**2 * cosmo.angular_diameter_distance(lens_z[behind]).value
                    * dls/cosmo.angular_diameter_distance(source_z[behind]).value)
    return beta


def profile_bin(records, role, k, cosmo):
    lens_z, z_est, z_ref, shape_w, rms, m, phi, oracle_phi = np.asarray(records, dtype=float).T
    beta_est = beta_from_distances(cosmo, lens_z, z_est)
    beta_ref = beta_from_distances(cosmo, lens_z, z_ref)
    pair_w = shape_w * beta_est**2
    total = pair_w.sum()
    response = 1-np.average(rms**2, weights=pair_w)
    mean_m = np.average(m, weights=pair_w)
    injection = 2e13/(1+.3*k) if role == 'lens' else 0.
    # A constant additive term in sky coordinates. Generate the injected
    # signal from Astropy's independent source position angle.
    additive = np.tile([.003, -.004], (len(records),1))
    tangent_sky = -np.column_stack([np.cos(2*oracle_phi), np.sin(2*oracle_phi)])
    distortion = 2*response*(additive + (beta_ref*injection*(1+mean_m)*(1+MSEL))[:,None]*tangent_sky)
    photoz, observed_weights = [], []
    calibration_cosmology = SimpleNamespace(comoving=lambda z: cosmo.comoving_distance(z).value)
    for lens in np.unique(lens_z):
        use = lens_z == lens
        photoz.append(wl_photoz.point_redshift_response(lens, z_est[use], z_ref[use], shape_w[use],
                                                       np.ones(use.sum()), calibration_cosmology)['response'])
        observed_weights.append(float(pair_w[use].sum()))
    correction = wl_photoz.stack_photoz_response(photoz, observed_weights)['response']
    oracle = float(np.sum(shape_w*beta_est*beta_ref)/total)
    if not np.isclose(correction, oracle, rtol=1e-12):
        raise ValueError('photo-z response differs from independent distance calculation')
    e = wl_y1.rotate_sky_components(*distortion.T, phi)
    c = wl_y1.rotate_sky_components(*additive.T, phi)
    args = dict(lens_weight=np.ones(len(records)), selection_m=MSEL)
    raw = wl_y1.y1_bin_delta_sigma(e,c,shape_w,beta_est,rms,m,photoz_response=1.,**args)
    fitted = wl_y1.y1_bin_delta_sigma(e,c,shape_w,beta_est,rms,m,photoz_response=correction,**args)
    if not np.allclose([fitted['tangential'], fitted['cross']], [injection,0], rtol=1e-9, atol=3e3):
        raise ValueError('noiseless injection failed tangential or cross recovery')
    return dict(role=role, radial_bin=k, n_pairs=len(records), injected_tangential=injection,
                photoz_response=correction, before_photoz=raw['tangential'],
                recovered_tangential=fitted['tangential'], recovered_cross=fitted['cross'])


def run(output):
    if output.exists():
        raise FileExistsError(output)
    centers, pairs, selected, available, _, provenance = load_inputs(PLAN, OBJECT_QA)
    center_by_id = {c['center_id']:c for c in centers}
    point_z = finite_point_redshifts(available)
    table = GEOMETRY/'pair_geometry.tsv'
    geometry_summary = GEOMETRY/'summary.json'
    if sha(SUPPORT) != SUPPORT_SHA or sha(geometry_summary) != GEOMETRY_SHA:
        raise ValueError('changed support or geometry provenance')
    support = json.loads(SUPPORT.read_text())
    geometry = json.loads(geometry_summary.read_text())
    if sha(table) != geometry['output_sha256'][table.name] or any(support[k]!=v or geometry[k]!=v for k,v in provenance.items()):
        raise ValueError('geometry or support lineage differs')
    fixed = set(support['common_0p1_to_3'][str(GAP)])
    by_bin = defaultdict(list)
    with table.open(newline='') as stream:
        for pair, row in zip(pairs, csv.DictReader(stream,delimiter='\t'), strict=True):
            cid,obj = pair['center_id'],pair['object_id']
            if (cid,obj) != (row['center_id'],row['object_id']):
                raise ValueError('geometry pair identity differs')
            c = center_by_id[cid]
            k = int(row['radial_bin'])
            zl = float(c['z_lens'])
            z = point_z.get(obj)
            if k == 0 or c['matched_lens_id'] not in fixed or obj not in selected or z is None or z <= zl+GAP:
                continue
            by_bin[c['role'],k].append((pair,c,float(row['source_angle_rad']),zl,z,.05 if int(obj)%10==0 else z-.05))
    cosmo = FlatLambdaCDM(H0=70., Om0=.3, Tcmb0=0.)
    results = []
    for role in ('lens','random'):
        for k in range(1,len(RADIAL_EDGES)-1):
            rows = by_bin[role,k]
            source = SkyCoord([float(p['source_ra']) for p,*_ in rows]*u.deg,
                              [float(p['source_dec']) for p,*_ in rows]*u.deg)
            center = SkyCoord([float(c['ra_deg']) for _,c,*_ in rows]*u.deg,
                              [float(c['dec_deg']) for _,c,*_ in rows]*u.deg)
            oracle_angles = np.pi/2-source.position_angle(center).rad
            records = [(zl,z,zref,float(p['weight']),float(p['rms_e']),float(p['m']),angle,oracle)
                       for (p,_,angle,zl,z,zref),oracle in zip(rows,oracle_angles,strict=True)]
            results.append(profile_bin(records,role,k,cosmo))
    output.mkdir(parents=True,exist_ok=False)
    result_table = output/'injected_profile.tsv'
    write_table(result_table,results)
    summary = dict(**provenance, science_measurement=False, synthetic_distortion=True, injected_selection_m=MSEL,
                   rows=results, n_pairs=sum(r['n_pairs'] for r in results), fixed_lenses=sorted(fixed), gap=GAP,
                   support_sha256=SUPPORT_SHA, geometry_sha256=GEOMETRY_SHA, artifact_sha256=sha(result_table),
                   runner_sha256=sha(Path(__file__)), y1_kernel_sha256=sha(Path(wl_y1.__file__)),
                   photoz_kernel_sha256=sha(Path(wl_photoz.__file__)),
                   limitation='Noiseless synthetic shape and reference redshifts on existing selected pairs. Within-bin input responsivity and multiplicative response are fixed at the weighted ensemble values.')
    with (output/'summary.json').open('x') as stream:
        json.dump(summary,stream,indent=2,allow_nan=False)
        stream.write('\n')
    print(json.dumps(dict(n_pairs=summary['n_pairs'],bins=len(results),fixed_lenses=sorted(fixed)),indent=2))


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    run(parser.parse_args().output.resolve())
