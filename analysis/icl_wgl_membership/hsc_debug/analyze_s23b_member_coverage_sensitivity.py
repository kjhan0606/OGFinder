#!/usr/bin/env python3
"""Compare member outer-light coverage cuts with fixed clusters and candidates.

Uses completed, audited v2 diagnostics only; no network or image access.
Native, common-cluster and common-candidate samples are reported separately.
For fixed samples, coverage nuisance covariates are frozen at the 0.6 reference
so paired slope differences use exactly the same design and spatial blocks.
All results are exploratory sensitivity estimates, not calibrated detections.
"""
import argparse
from collections import defaultdict
import json
import os
from pathlib import Path
import tempfile

import numpy as np
from analyze_s23b_member_outer_icl_relation import (
    TARGET_MANIFEST, CONTROL_MANIFEST, PAIR_PROFILES, EXPECTED_HASHES,
    read_json, read_tsv, write_tsv, sha256, annular_support,
    model_matrix, fit_block_slope, connected_pair_blocks,
)

FLOORS = (.4, .6, .8)
MODES = ('native', 'common_clusters', 'common_candidates')
SEED = 20261117


def candidate_profiles(diagnostics, floors=FLOORS, min_controls=3):
    """Return nested candidate supports while recomputing each control mean."""
    out = {floor: {} for floor in floors}
    seen = set()
    for r in diagnostics:
        if r['radial_bin_kpc'] != '50-75':
            continue
        key = (r['target_id'], r['filter'], r['member_index'])
        if key in seen:
            raise ValueError('duplicate outer-bin candidate')
        seen.add(key)
        if r['sep_matched'] != 'True' or r['target_profile_available'] != 'True':
            continue
        controls = json.loads(r['control_bin_measurements_json'])
        for floor in floors:
            good = [c for c in controls if c['coverage'] >= floor]
            if float(r['target_coverage']) < floor or len(good) < min_controls:
                continue
            target_sb = float(r['target_sb_njy_arcsec2'])
            control_sb = float(np.mean([c['sb_njy_arcsec2'] for c in good]))
            out[floor][key] = dict(delta=target_sb-control_sb,
                target_coverage=float(r['target_coverage']),
                control_coverage=float(np.mean([c['coverage'] for c in good])))
    for lo, hi in zip(floors[:-1], floors[1:]):
        if not set(out[hi]).issubset(out[lo]):
            raise ValueError('coverage supports must be nested')
    return out


def cluster_profiles(candidates, mode):
    common_candidates = set.intersection(*(set(v) for v in candidates.values()))
    common_clusters = set.intersection(*({key[:2] for key in v} for v in candidates.values()))
    result = {}
    for floor, pool in candidates.items():
        grouped = defaultdict(list)
        for key, value in pool.items():
            if mode == 'common_candidates' and key not in common_candidates:
                continue
            if mode == 'common_clusters' and key[:2] not in common_clusters:
                continue
            grouped[key[:2]].append(value)
        result[floor] = {key: dict(
            member_outer_delta_njy_arcsec2=float(np.mean([r['delta'] for r in group])),
            median_member_target_coverage_50_75=float(np.median([r['target_coverage'] for r in group])),
            median_member_control_coverage_50_75=float(np.median([r['control_coverage'] for r in group])),
            n_member_positions_50_75=len(group)) for key,group in grouped.items()}
    if mode != 'native':
        ids = [set(v) for v in result.values()]
        if any(s != ids[0] for s in ids[1:]):
            raise ValueError('fixed sample has different cluster support')
        for pool in result.values():
            for key,row in pool.items():
                for nuisance in ('median_member_target_coverage_50_75', 'median_member_control_coverage_50_75'):
                    row[nuisance] = result[.6][key][nuisance]
    return result


def fit(rows, width, phase, replicas, seed, difference=None):
    record = dict(n_target_clusters=len(rows), status='insufficient_sample',
                  n_connected_blocks=0, slope=None, ci_low=None, ci_high=None, valid_draws=0, reason='')
    if len(rows) < 10:
        return record
    try:
        X,y,_ = model_matrix(rows)
    except ValueError as exc:
        record.update(status='rank_deficient_model', reason=str(exc))
        return record
    controls = {r['target_id']:dict(ra_deg=r['control_ra_deg'],dec_deg=r['control_dec_deg']) for r in rows}
    labels,_ = connected_pair_blocks(rows, controls, width, phase, 84.056247280625)
    if difference is not None:
        y = np.array(difference)
    record['n_connected_blocks'] = len(np.unique(labels))
    try:
        result = fit_block_slope(X,y,labels,replicas,seed)
    except ValueError as exc:
        record.update(status='inadequate_bootstrap_support', reason=str(exc))
        return record
    record.update(result)
    record['n_connected_blocks'] = result['n_blocks']
    return record


def run(run_root, output, replicas, measurement=None, audit_path=None):
    if output.exists():
        raise FileExistsError(f'refusing to overwrite {output}')
    measurement = measurement or run_root/'measurement'
    audit_path = audit_path or run_root/'audit/summary.json'
    audit = read_json(audit_path)
    ms = read_json(measurement/'summary.json')
    if (audit.get('historical_profile_rows_identical') is not True
            or audit.get('measurement_summary_sha256') != sha256(measurement/'summary.json')
            or ms.get('science_measurement') is not False
            or ms.get('schema') != 'astrafek-s23b-camira-member-outer-profile-stack-v2'):
        raise ValueError('completed audited v2 diagnostics required')
    for name,digest in ms['output_sha256'].items():
        if Path(name).name != name or sha256(measurement/name) != digest:
            raise ValueError('measurement output changed')
    for p in (TARGET_MANIFEST,CONTROL_MANIFEST,PAIR_PROFILES):
        if sha256(p) != EXPECTED_HASHES[str(p)]:
            raise ValueError('pinned context changed')
    selected = read_tsv(measurement/'selected_target_frames.tsv')
    targets = {str(r['id']):r for r in read_json(TARGET_MANIFEST)['clusters']}
    controls = {str(r['target_cluster_id']):r for r in read_json(CONTROL_MANIFEST)['clusters']}
    paired = {(r['target_id'],r['filter']):r for r in read_tsv(PAIR_PROFILES)}
    context = {}
    for r in selected:
        key = (r['target_id'],r['filter'])
        t,c,p = targets[key[0]], controls[key[0]], paired[key]
        context[key] = dict(target_id=key[0],filter=key[1],ra_deg=float(t['ra_deg']),dec_deg=float(t['dec_deg']),
            control_ra_deg=float(c['ra_deg']),control_dec_deg=float(c['dec_deg']),
            z_cl=float(r['z']),richness=float(r['richness']),
            icl_contrast_njy_arcsec2=float(r['icl_contrast_50_300_minus_300_800_njy_arcsec2']),
            min_icl_common_support_fraction=min(annular_support(p['profile_json'],50,300),annular_support(p['profile_json'],300,800)))
    raw = candidate_profiles(read_tsv(measurement/'member_profile_bin_diagnostics.tsv'))
    # Validate that the reference flux reconstruction reproduces the saved bin table.
    saved = {(r['target_id'],r['filter'],r['member_index']):float(r['delta_sb_njy_arcsec2'])
             for r in read_tsv(measurement/'member_profile_bins.tsv') if r['radial_bin_kpc']=='50-75'}
    if set(saved) != set(raw[.6]) or any(abs(saved[k]-raw[.6][k]['delta'])>1e-9 for k in saved):
        raise ValueError('reference outer flux does not reconstruct')
    cluster_rows, fits, differences = [], [], []
    for mi,mode in enumerate(MODES):
        pools = cluster_profiles(raw,mode)
        assembled = {}
        for floor,pool in pools.items():
            assembled[floor] = {k:dict(**context[k],**v) for k,v in pool.items()}
            cluster_rows.extend(dict(mode=mode,coverage_floor=floor,**r,science_measurement=False)
                                for r in assembled[floor].values())
        for band in ('i','i2'):
            for width in (2.5,5.,10.):
                for phase in (0.,.5):
                    for floor in FLOORS:
                        rows = [r for k,r in sorted(assembled[floor].items()) if k[1]==band]
                        result = fit(rows,width,phase,replicas,SEED+mi)
                        fits.append(dict(mode=mode,filter=band,coverage_floor=floor,
                            block_width_deg=width,grid_phase_fraction=phase,**result,science_measurement=False))
                        if mode != 'native' and floor != .6:
                            dy = [r['member_outer_delta_njy_arcsec2']-assembled[.6][(r['target_id'],band)]['member_outer_delta_njy_arcsec2'] for r in rows]
                            df = fit(rows,width,phase,replicas,SEED+mi,difference=dy)
                            differences.append(dict(mode=mode,filter=band,coverage_floor=floor,
                                reference_floor=.6,block_width_deg=width,grid_phase_fraction=phase,
                                **df,science_measurement=False))
    output.parent.mkdir(parents=True,exist_ok=True)
    tmp = Path(tempfile.mkdtemp(prefix=output.name+'.partial.',dir=output.parent))
    for name,rows in [('cluster_threshold_profiles.tsv',cluster_rows),('coverage_slope_layouts.tsv',fits),('paired_slope_changes.tsv',differences)]:
        write_tsv(tmp/name,rows,list(rows[0]))
    summary = dict(schema='astrafek-s23b-member-coverage-sensitivity-v1',
        status='complete_exploratory_fixed_support_sensitivity',
        run_root=str(run_root),measurement=str(measurement),audit_summary=str(audit_path),
        n_frames=len(selected),n_bootstrap=replicas,
        modes=list(MODES),coverage_floors=list(FLOORS),
        measurement_summary_sha256=sha256(measurement/'summary.json'),audit_summary_sha256=sha256(audit_path),
        input_context_sha256={str(p):sha256(p) for p in (TARGET_MANIFEST,CONTROL_MANIFEST,PAIR_PROFILES)},
        script_sha256=sha256(Path(__file__)),
        helper_sha256=sha256(Path(__file__).with_name('analyze_s23b_member_outer_icl_relation.py')),
        output_sha256={p.name:sha256(p) for p in tmp.iterdir()},
        primary_layout=[r for r in fits if r['block_width_deg']==5 and r['grid_phase_fraction']==0],
        primary_paired_changes=[r for r in differences if r['block_width_deg']==5 and r['grid_phase_fraction']==0],
        science_measurement=False,images_read=False,network_used=False,
        limitations=['Native samples differ across thresholds.',
            'Fixed modes freeze coverage nuisance covariates at the 0.6 reference.',
            'Common candidates isolates control-mean changes, conditional on the strictest support.',
            'Intervals are marginal sensitivity summaries, not multiple-testing detection claims.',
            'PSF, source Poisson, sky, covariance and photometric membership systematics remain unresolved.'])
    (tmp/'summary.json').write_text(json.dumps(summary,indent=2,allow_nan=False)+'\n')
    os.rename(tmp,output)
    return summary


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run-root',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--measurement',type=Path,help='explicit measurement directory, for a separately stored pilot')
    parser.add_argument('--audit-summary',type=Path,help='explicit completed audit summary')
    parser.add_argument('--replicas',type=int,default=5000)
    args=parser.parse_args()
    if args.replicas < 100:
        parser.error('at least 100 bootstrap draws required')
    result=run(args.run_root.resolve(),args.output.resolve(),args.replicas,
        args.measurement.resolve() if args.measurement else None,
        args.audit_summary.resolve() if args.audit_summary else None)
    print(json.dumps({'status':result['status'],'n_frames':result['n_frames'],'primary_layout':result['primary_layout']},indent=2))
