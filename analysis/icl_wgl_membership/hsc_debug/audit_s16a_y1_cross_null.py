#!/usr/bin/env python3
"""Cross-component null check for the fixed, point-z S16A pilot."""
import argparse
from collections import defaultdict
import csv
import json
from pathlib import Path

import numpy as np
from scipy.stats import chi2

from audit_s16a_y1_reuse_shape_noise import covariance
from verify_s16a_y1_profile_injection import (
    GAP, GEOMETRY, GEOMETRY_SHA, OBJECT_QA, PLAN, RADIAL_EDGES, SUPPORT, SUPPORT_SHA,
    FlatLambdaCDM, beta_from_distances, finite_point_redshifts, load_inputs, sha, wl_y1, write_table,
)


def extract_pairs():
    centers, pairs, selected, available, _, provenance = load_inputs(PLAN, OBJECT_QA)
    center_by_id = {c['center_id']:c for c in centers}
    point_z = finite_point_redshifts(available)
    summary_file, table = GEOMETRY/'summary.json', GEOMETRY/'pair_geometry.tsv'
    if sha(SUPPORT) != SUPPORT_SHA or sha(summary_file) != GEOMETRY_SHA:
        raise ValueError('changed verified support or geometry summary')
    geometry, support = (json.loads(p.read_text()) for p in (summary_file,SUPPORT))
    if (sha(table) != geometry['output_sha256'][table.name]
            or any(geometry[k]!=v or support[k]!=v for k,v in provenance.items())):
        raise ValueError('geometry or support lineage differs')

    fixed = set(support['common_0p1_to_3'][str(GAP)])
    bins = defaultdict(list)
    with table.open(newline='') as stream:
        for pair, geo in zip(pairs,csv.DictReader(stream,delimiter='\t'),strict=True):
            cid,obj = pair['center_id'],pair['object_id']
            if (cid,obj) != (geo['center_id'],geo['object_id']):
                raise ValueError('pair ledgers differ')
            c = center_by_id[cid]
            k,z_lens = int(geo['radial_bin']),float(c['z_lens'])
            z_source = point_z.get(obj)
            if k and c['matched_lens_id'] in fixed and obj in selected and z_source is not None and z_source>z_lens+GAP:
                bins[c['role'],k].append((pair,float(geo['source_angle_rad']),z_lens,z_source))
    return bins, sorted(fixed), provenance


def estimate(bins):
    cosmo = FlatLambdaCDM(H0=70., Om0=.3, Tcmb0=0.)
    n_bins = len(RADIAL_EDGES)-2
    rows, origin, pair_ids, angles, amplitudes, diagnostics = [], [], [], [], [], []
    source_weights = {}
    for j,role in enumerate(('lens','random')):
        for k in range(1,n_bins+1):
            records = bins[role,k]
            pairs, phi, lens_z, source_z = zip(*records,strict=True)
            phi, lens_z, source_z = map(np.asarray,(phi,lens_z,source_z))
            beta = beta_from_distances(cosmo,lens_z,source_z)
            weight,rms,m,e1,e2,c1,c2 = np.array([[float(p[name]) for name in
                ('weight','rms_e','m','e1','e2','c1','c2')] for p in pairs]).T
            for p,w in zip(pairs,weight,strict=True):
                if source_weights.setdefault(p['object_id'],w) != w:
                    raise ValueError('source weight differs between repeated pairs')
            e_cross = wl_y1.rotate_sky_components(e1,e2,phi)[:,1]
            c_cross = wl_y1.rotate_sky_components(c1,c2,phi)[:,1]
            zeros = np.zeros_like(weight)
            result = wl_y1.y1_bin_delta_sigma(np.column_stack((zeros,e_cross)), np.column_stack((zeros,c_cross)),
                     weight,beta,rms,m,lens_weight=np.ones_like(weight),selection_m=0,photoz_response=1)
            gain = 2*result['responsivity']*(1+result['mean_m'])
            linear = np.sqrt(weight)*beta/(result['weight_sum']*gain)
            index = 2*(j*n_bins + k-1)
            rows.extend([index]*len(pairs))
            origin.extend(p['object_id'] for p in pairs)
            pair_ids.extend(f'{p["center_id"]}:{p["object_id"]}' for p in pairs)
            angles.extend(phi)
            amplitudes.extend(linear)
            diagnostics.append(dict(role=role,annulus=k,n_pairs=len(pairs),cross=result['cross'],
                responsivity=result['responsivity'], mean_m=result['mean_m'], independent_variance=1/(result['weight_sum']*gain**2)))
    args = (np.asarray(rows), np.asarray(angles), np.asarray(amplitudes), 4*n_bins)
    actual = covariance(args[0],origin,*args[1:])[1::2,1::2]
    independent = covariance(args[0],pair_ids,*args[1:])[1::2,1::2]
    if not np.allclose(np.diag(independent),[r['independent_variance'] for r in diagnostics],rtol=1e-12,atol=0):
        raise ValueError('independent-pair diagonal does not match the analytic variance')
    transform = np.column_stack((np.eye(n_bins),-np.eye(n_bins)))
    d = np.array([r['cross'] for r in diagnostics])
    tests = {}
    for role, projection in (('lens',np.eye(2*n_bins)[:n_bins]),('random',np.eye(2*n_bins)[n_bins:]),('lens_minus_random',transform)):
        vector, cov = projection @ d, projection @ actual @ projection.T
        # Validate positive-definiteness and solve without forming an inverse.
        lower = np.linalg.cholesky(cov)
        stat = float(np.square(np.linalg.solve(lower,vector)).sum())
        tests[role] = dict(chi2=stat,dof=n_bins,p_value_shape_only=float(chi2.sf(stat,n_bins)))
    return diagnostics, actual, independent, transform@d, transform@actual@transform.T, tests


def run(output):
    if output.exists():
        raise FileExistsError(output)
    bins, fixed, provenance = extract_pairs()
    records, actual, independent, difference, diff_cov, tests = estimate(bins)
    if sum(r['n_pairs'] for r in records) != 48852:
        raise ValueError('expected the sample from the profile injection')
    output.mkdir(parents=True,exist_ok=False)
    table, matrix = output/'cross_null.tsv', output/'cross_null_covariance.npz'
    n_bins = len(difference)
    for k,r in enumerate(records):
        r['shape_sd'] = float(np.sqrt(actual[k,k]))
    records.extend(dict(role='lens_minus_random',annulus=k+1,n_pairs=None,cross=float(difference[k]),
                        responsivity=None,mean_m=None,independent_variance=float(independent[k,k]+independent[k+n_bins,k+n_bins]),
                        shape_sd=float(np.sqrt(diff_cov[k,k]))) for k in range(n_bins))
    write_table(table,records)
    np.savez_compressed(matrix, role_cross_covariance=actual, independent_pair_covariance=independent,
                        difference_cross_covariance=diff_cov)
    summary = dict(**provenance, fixed_lenses=fixed, gap=GAP, units='Msun per proper Mpc^2', tests=tests,
                   science_measurement=False, tangential_measured=False,
                   limitations='Point-z candidate selection with catalog ensemble correction. Shape-only covariance conditioned on fixed support. Selection, redshift-dependent responsivity, COSMOS dilution and spatial survey covariance are pending.',
                   support_sha256=SUPPORT_SHA, geometry_sha256=GEOMETRY_SHA,
                   output_sha256={p.name:sha(p) for p in (table,matrix)},runner_sha256=sha(Path(__file__)))
    with (output/'summary.json').open('x') as stream:
        json.dump(summary,stream,indent=2,allow_nan=False)
        stream.write('\n')
    print(json.dumps(tests,indent=2))


if __name__=='__main__':
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    run(p.parse_args().output.resolve())
