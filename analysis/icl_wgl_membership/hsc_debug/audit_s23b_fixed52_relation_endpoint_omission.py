#!/usr/bin/env python3
"""Refit ICL–WL relations after deleting both WL endpoints in a survey field."""
import argparse
from collections import Counter,defaultdict
import json
import math
from pathlib import Path

import numpy as np

from audit_s23b_fixed52_relation_field_omission import (
    EXPECTED_SHA256, FIELD, FIELD_SHA, MIN_LENSES, load_and_check, make_models,
    pair_block_labels, partial_r, sha256, write_tsv,
)
from adjust_s23b_icl_satellite_wl_fixed52_max_t import read_tsv


def load_profiles(directory):
    legacy_summary = FIELD/'summary.json'
    if sha256(legacy_summary) != FIELD_SHA:
        raise ValueError('historical field deletion summary changed')
    legacy = json.loads(legacy_summary.read_text())
    actual_file = directory/'summary.json'
    actual = json.loads(actual_file.read_text())
    for key in ('status','inputs','sample','full_profiles_msun_mpc2','region_deletions'):
        if legacy[key] != actual[key]:
            raise ValueError(f'replay did not reproduce historical field result: {key}')
    if actual['science_measurement'] is not False:
        raise ValueError('unexpected science flag')
    if sha256(directory/'region_profile_shifts.tsv') != legacy['outputs']['region_profile_shifts_sha256']:
        raise ValueError('aggregate shift table differs from historical result')
    table = directory/'region_lens_profiles.tsv'
    if sha256(table) != actual['outputs']['region_lens_profiles_sha256']:
        raise ValueError('per-lens export changed')
    profiles, radial = defaultdict(dict), defaultdict(list)
    for row in read_tsv(table):
        field, lens, k = row['region_label'],row['matched_lens_id'],int(row['radial_bin'])
        key = (lens,k)
        if key in profiles[field] or int(row['n_controls']) < 4:
            raise ValueError('duplicated or unsupported per-lens profile')
        values = tuple(float(row[c]) for c in ('pooled_t','pooled_cross','equal_t','equal_cross'))
        if not np.all(np.isfinite(values)):
            raise ValueError('nonfinite per-lens profile')
        profiles[field][key] = values
        radial[field,k].append(values)
    sizes = {r['region_label']:r['common_supported_lens_sets_after_deletion'] for r in legacy['region_deletions']}
    if set(sizes) != set(profiles):
        raise ValueError('field sets differ')
    for field, n in sizes.items():
        ids = {lens for lens,_ in profiles[field]}
        if len(ids)!=n or {(lens,k) for lens in ids for k in range(7)}!=set(profiles[field]):
            raise ValueError('per-lens annulus support differs')
    for r in read_tsv(directory/'region_profile_shifts.tsv'):
        v = np.mean(radial[r['region_label'],int(r['radial_bin'])],axis=0)
        a = {'pooled_controls':v[:2], 'equal_control_centers':v[2:], 'equal_minus_pooled':v[2:]-v[:2]}[r['estimator']]
        expected = [float(r[c]) for c in ('leave_region_out_delta_sigma_msun_mpc2','leave_region_out_cross_msun_mpc2')]
        if not np.allclose(a,expected,rtol=2e-9,atol=.01):
            raise ValueError('per-lens export fails to reproduce aggregate profile')
    return profiles, sizes, sha256(actual_file)


def refit(model, field, profiles):
    k = int(model['radial_bin'])
    keep = np.array([(lens,k) in profiles for lens in model['cluster_ids']])
    kept = [i for i,include in zip(model['cluster_ids'],keep,strict=True) if include]
    design = model['design'][keep]
    rank = int(np.linalg.matrix_rank(design)) if kept else 0
    matched_r, endpoint_r = None,None
    status = 'below_original_lens_floor' if len(kept) < MIN_LENSES else 'rank_deficient'
    if len(kept) >= MIN_LENSES and rank == design.shape[1]:
        old = (model['reduced_fit']+model['restricted_resid'])[keep]
        column = 0 if model['component']=='tangential' else 1
        updated = np.array([profiles[lens,k][column] for lens in kept])
        matched_r = partial_r(design[:,-1],old,design[:,:-1])
        endpoint_r = partial_r(design[:,-1],updated,design[:,:-1])
        if math.isfinite(matched_r) and math.isfinite(endpoint_r):
            status = 'complete'
        else:
            status,matched_r,endpoint_r = 'undefined_correlation',None,None
    reference = model['partial_r']
    return dict(**{key:model[key] for key in ('filter','support_group','relationship','predictor','component','radial_bin')},
                omitted_field=field,status=status,n_full=len(model['cluster_ids']),n_retained=len(kept),
                retained_ids='|'.join(kept),full_r=reference,retained_baseline_r=matched_r,endpoint_removed_r=endpoint_r,
                sample_shift=None if matched_r is None else matched_r-reference,
                endpoints_only_shift=None if endpoint_r is None else endpoint_r-matched_r,
                total_shift=None if endpoint_r is None else endpoint_r-reference)


def run(directory,output):
    if output.exists():
        raise FileExistsError(output)
    profile, sizes, field_sha = load_profiles(directory)
    _,metrics,joined,icl = load_and_check()
    ids = sorted(icl)
    target = [{f'{a}_deg':float(icl[i][f'{a}_deg']) for a in ('ra','dec')} for i in ids]
    random = [{f'{a}_deg':float(icl[i][f'control_{a}_deg']) for a in ('ra','dec')} for i in ids]
    _,labels,_,_,_,_ = pair_block_labels(target,random,np.ones(len(ids),dtype=bool))
    models = [m for m in make_models(metrics,joined,icl,dict(zip(ids,map(str,labels),strict=True)))
              if m['relationship'] != 'icl_vs_satellite']
    if len(models)!=112:
        raise ValueError('unexpected WL relation family')
    rows = [refit(m,field,p) for field,p in sorted(profile.items()) for m in models]
    complete = [r for r in rows if r['status']=='complete']
    output.mkdir(parents=True,exist_ok=False)
    table = output/'endpoint_omission_relations.tsv'
    write_tsv(table,rows)
    summary = dict(science_measurement=False, estimator='pooled_controls', all_bin_support_sizes=sizes,
        n_comparisons=len(rows), outcome_counts=dict(Counter(r['status'] for r in rows)),
        endpoint_only_sign_reversals=sum(r['retained_baseline_r']*r['endpoint_removed_r']<0 for r in complete),
        full_to_endpoint_sign_reversals=sum(r['full_r']*r['endpoint_removed_r']<0 for r in complete),
        largest_endpoint_changes=sorted(complete,key=lambda r:abs(r['endpoints_only_shift']),reverse=True)[:5],
        largest_total_changes=sorted(complete,key=lambda r:abs(r['total_shift']),reverse=True)[:5],
        calibration_scope='WL endpoints only. For the same surviving groups, ICL/member measurements and covariates are kept fixed.',
        limitations='Three-field partial-correlation influence diagnostic, 12-lens floor. No significance or covariance inference.',
        input_sha256={**EXPECTED_SHA256,'field_replay_summary':field_sha,'historical_field_summary':FIELD_SHA},
        runner_sha256=sha256(Path(__file__)),output_sha256={table.name:sha256(table)})
    with (output/'summary.json').open('x') as stream:
        json.dump(summary,stream,indent=2,allow_nan=False)
        stream.write('\n')
    print(json.dumps({k:summary[k] for k in ('all_bin_support_sizes','outcome_counts','endpoint_only_sign_reversals','full_to_endpoint_sign_reversals')},indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--field-replay',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    run(args.field_replay.resolve(),args.output.resolve())
