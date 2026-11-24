#!/usr/bin/env python3
"""Summarize sky-stratified mask counts using clusters and sky-cell bootstrap.

Cells, not individual correlated stamps, are resampled. The equal-cluster
estimand describes the selected sample and is not population weighted.
"""
import argparse
from collections import defaultdict
import csv
import hashlib
import json
from pathlib import Path
import numpy as np

METRICS={
 'quality_bits':'quality_bits_pixels', 'bright_object':'bright_object_pixels',
 'nonfinite_image':'nonfinite_image_pixels', 'invalid_variance':'invalid_variance_pixels',
 'inexact_psf':'MP_INEXACT_PSF_pixels', 'rejected':'MP_REJECTED_pixels',
 'base':'quality_invalid_pixels', 'segmentation':'segmentation_pixels',
 'neighbor_shapes':'neighbor_shapes_pixels', 'unmatched_circles':'unmatched_circles_pixels',
 'base_exclusive':'quality_invalid_exclusive_pixels',
 'neighbor_shapes_exclusive':'neighbor_shapes_exclusive_pixels',
 'bright_object_base_exclusive':'bright_object_base_exclusive_pixels',
 'outside_image':'outside_image_annulus_pixels', 'unmasked':'unmasked_pixels',
}


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def cell_bootstrap(values,cells,replicas,seed):
    unique=sorted(set(cells))
    if not unique:raise ValueError('no sky cells')
    groups={c:np.flatnonzero(np.array([x==c for x in cells])) for c in unique}
    rng=np.random.default_rng(seed)
    draws=[]
    for _ in range(replicas):
        selected=rng.integers(0,len(unique),len(unique))
        indices=np.concatenate([groups[unique[int(i)]] for i in selected])
        draws.append(float(np.mean(values[indices])))
    return np.quantile(draws,[.025,.975])


def run(audit,plan_path,out,replicas):
    if out.exists():raise FileExistsError(out)
    summary=json.loads((audit/'summary.json').read_text())
    plan=json.loads(plan_path.read_text())
    table=audit/'stamp_component_counts.tsv'
    if (summary.get('status')!='complete_exploratory_spatial_mask_sample'
        or summary.get('selection_plan_sha256')!=sha(plan_path)
        or not all(summary.get(k) is True for k in ('mask_union_exact','base_mask_union_exact','discrete_annulus_reconciled'))
        or summary.get('output_sha256')!=sha(table)):
        raise ValueError('predeclared sample or audit integrity failed')
    planned={(r['target_id'],r['filter']):tuple(r['sky_cell']) for r in plan['sample']}
    recorded={tuple(x) for x in summary['frame_ids']}
    if recorded!=set(planned) or len(planned)!=96:raise ValueError('declared and audited frame sets differ')
    grouped=defaultdict(list)
    with table.open(newline='') as f:
        for r in csv.DictReader(f,delimiter='\t'):
            key=(r['target_id'],r['filter'])
            if key not in planned or r['role'] not in ('target','control'):raise ValueError('unexpected stamp')
            grouped[(key[0],key[1],r['role'])].append(r)
    rows=[]
    for (ident,band,role),stamps in sorted(grouped.items()):
        inside=sum(int(r['in_image_ring_pixels']) for r in stamps)
        full=sum(int(r['full_discrete_annulus_pixels']) for r in stamps)
        outside=sum(int(r['outside_image_annulus_pixels']) for r in stamps)
        if not inside or full!=inside+outside:raise ValueError('annulus area mismatch')
        row=dict(target_id=ident,filter=band,role=role,sky_cell_x=planned[(ident,band)][0],
                 sky_cell_y=planned[(ident,band)][1],n_stamps=len(stamps))
        for name,column in METRICS.items():
            row[name]=sum(int(r[column]) for r in stamps)/(full if name=='outside_image' else inside)
        rows.append(row)
    summaries=[]
    for band in ('i','i2'):
        for role in ('target','control'):
            group=[r for r in rows if r['filter']==band and r['role']==role]
            if not group:raise ValueError('no stamps for band and role')
            cells=[(r['sky_cell_x'],r['sky_cell_y']) for r in group]
            for name in METRICS:
                values=np.array([r[name] for r in group],dtype=float)
                interval=cell_bootstrap(values,cells,replicas,20261117)
                summaries.append(dict(filter=band,role=role,metric=name,
                    n_planned=sum(x['filter']==band for x in plan['sample']),
                    n_clusters=len(group),n_cells=len(set(cells)),n_stamps=sum(r['n_stamps'] for r in group),
                    mean=float(values.mean()),median=float(np.median(values)),
                    cell_bootstrap_ci_low=float(interval[0]),cell_bootstrap_ci_high=float(interval[1])))
    out.mkdir(parents=True,exist_ok=False)
    for name,data in [('cluster_fractions.tsv',rows),('summary_by_filter_role.tsv',summaries)]:
        with (out/name).open('x',newline='') as f:
            writer=csv.DictWriter(f,fieldnames=list(data[0]),delimiter='\t',lineterminator='\n');writer.writeheader();writer.writerows(data)
    result=dict(status='complete_exploratory_spatial_mask_summary',science_measurement=False,
        audit_summary_sha256=sha(audit/'summary.json'),sample_plan_sha256=sha(plan_path),
        script_sha256=sha(Path(__file__)),output_sha256={p.name:sha(p) for p in out.iterdir()},
        bootstrap_replicas=replicas,seed=20261117,estimand='equal-cluster mean among sampled clusters with measured stamps',
        limitations=['Stratified sample is not population-weighted.', 'Sky-cell bootstrap intervals are descriptive and reflect only sampled cells.',
         'Masks overlap and exclusive fractions are not causal.', 'Unmatched candidates and clusters without measured stamps are outside the estimand.'])
    (out/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--audit',required=True,type=Path);p.add_argument('--plan',required=True,type=Path)
    p.add_argument('--output',required=True,type=Path);p.add_argument('--replicas',type=int,default=5000)
    a=p.parse_args()
    if a.replicas<100:p.error('at least 100 bootstrap draws required')
    run(a.audit.resolve(),a.plan.resolve(),a.output.resolve(),a.replicas)
