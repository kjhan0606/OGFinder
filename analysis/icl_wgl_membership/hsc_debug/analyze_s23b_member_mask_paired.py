#!/usr/bin/env python3
"""Describe paired target/control mask differences and exploratory associations.

This reads the predeclared sample's saved cluster fractions and pinned ICL
profile context. It performs no image reads or network requests. Intervals use
shared 5-degree sky-cell resampling, and are not population-weighted.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
from scipy.stats import rankdata
from analyze_s23b_member_outer_icl_relation import (
    PAIR_PROFILES, RELATIVE_TABLE, EXPECTED_HASHES, annular_support,
)

METRICS=('quality_bits','bright_object','neighbor_shapes','base_exclusive','inexact_psf','unmasked')
SEED=20261117


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def spearman(x,y):
    a,b=rankdata(x),rankdata(y)
    a=a-a.mean();b=b-b.mean()
    den=np.sqrt((a@a)*(b@b))
    return float(a@b/den) if den>0 else None


def block_indices(cells,replicas,seed):
    unique=sorted(set(cells))
    groups=[np.flatnonzero(np.array([c==u for c in cells])) for u in unique]
    rng=np.random.default_rng(seed)
    return [np.concatenate([groups[int(i)] for i in rng.integers(0,len(unique),len(unique))])
            for _ in range(replicas)]


def run(summary_dir,plan_path,out,replicas):
    if out.exists():raise FileExistsError(out)
    summary=json.loads((summary_dir/'summary.json').read_text())
    plan=json.loads(plan_path.read_text())
    table=summary_dir/'cluster_fractions.tsv'
    if (summary.get('status')!='complete_exploratory_spatial_mask_summary'
        or summary.get('science_measurement') is not False
        or summary.get('sample_plan_sha256')!=sha(plan_path)
        or summary['output_sha256']['cluster_fractions.tsv']!=sha(table)):
        raise ValueError('selection or cluster summary integrity failed')
    for p in (PAIR_PROFILES,RELATIVE_TABLE):
        if sha(p)!=EXPECTED_HASHES[str(p)]:raise ValueError('pinned ICL context changed')
    with table.open(newline='') as f: rows=list(csv.DictReader(f,delimiter='\t'))
    by_key={}
    for r in rows:
        k=(r['target_id'],r['filter'],r['role'])
        if k in by_key:raise ValueError('duplicate cluster-role row')
        by_key[k]=r
    planned={(r['target_id'],r['filter']):tuple(r['sky_cell']) for r in plan['sample']}
    if len(planned)!=96:raise ValueError('planned sample does not contain 96 unique frames')
    with PAIR_PROFILES.open(newline='') as f:
        pairs={(r['target_id'],r['filter']):r for r in csv.DictReader(f,delimiter='\t')}
    with RELATIVE_TABLE.open(newline='') as f:
        contrast={(r['target_id'],r['filter']):float(r['contrast_njy_arcsec2']) for r in csv.DictReader(f,delimiter='\t')
                  if r['contrast']=='50-300_minus_300-800' and r['contrast_njy_arcsec2']}
    paired=[]
    for k,cell in sorted(planned.items()):
        a=by_key.get((*k,'target'));b=by_key.get((*k,'control'))
        if a is None or b is None:continue
        if (int(a['sky_cell_x']),int(a['sky_cell_y']))!=cell or (int(b['sky_cell_x']),int(b['sky_cell_y']))!=cell:
            raise ValueError('sky cell mismatch')
        prof=pairs[k]['profile_json']
        record=dict(target_id=k[0],filter=k[1],sky_cell_x=cell[0],sky_cell_y=cell[1],
            icl_contrast=contrast[k],icl_min_common_support=min(annular_support(prof,50,300),annular_support(prof,300,800)))
        for name in METRICS:
            record['target_'+name]=float(a[name]);record['control_'+name]=float(b[name])
            record['difference_'+name]=float(b[name])-float(a[name])
        if not all(np.isfinite(v) for v in record.values() if isinstance(v,float)):
            raise ValueError('nonfinite paired data')
        paired.append(record)
    outputs=[]
    for band in ('i','i2'):
        group=[r for r in paired if r['filter']==band]
        if len(group)<10:raise ValueError('insufficient paired clusters')
        cells=[(r['sky_cell_x'],r['sky_cell_y']) for r in group]
        draws=block_indices(cells,replicas,SEED)
        for metric in METRICS:
            values=np.array([r['difference_'+metric] for r in group])
            estimates=np.array([values[idx].mean() for idx in draws])
            low,high=np.quantile(estimates,[.025,.975])
            outputs.append(dict(filter=band,kind='paired_mean_control_minus_target',metric=metric,
                context='',n_clusters=len(group),n_cells=len(set(cells)),estimate=float(values.mean()),
                interval_low=float(low),interval_high=float(high),valid_draws=len(estimates),science_measurement=False))
            for context in ('icl_contrast','icl_min_common_support'):
                x=np.array([r[context] for r in group])
                estimate=spearman(values,x)
                corr=[spearman(values[idx],x[idx]) for idx in draws]
                corr=np.array([v for v in corr if v is not None])
                if estimate is None or len(corr)<replicas*.9:raise ValueError('insufficient variable bootstrap correlations')
                low,high=np.quantile(corr,[.025,.975])
                outputs.append(dict(filter=band,kind='spearman_descriptive',metric=metric,context=context,
                    n_clusters=len(group),n_cells=len(set(cells)),estimate=estimate,
                    interval_low=float(low),interval_high=float(high),valid_draws=len(corr),science_measurement=False))
    out.mkdir(parents=True,exist_ok=False)
    for name,data in [('paired_cluster_fractions.tsv',paired),('paired_mask_results.tsv',outputs)]:
        with (out/name).open('x',newline='') as f:
            w=csv.DictWriter(f,fieldnames=list(data[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(data)
    result=dict(status='complete_exploratory_paired_mask_analysis',science_measurement=False,network_used=False,images_read=False,
        bootstrap_replicas=replicas,seed=SEED,cluster_summary_sha256=sha(summary_dir/'summary.json'),
        sample_plan_sha256=sha(plan_path),context_sha256={str(p):sha(p) for p in (PAIR_PROFILES,RELATIVE_TABLE)},
        script_sha256=sha(Path(__file__)),output_sha256={p.name:sha(p) for p in out.iterdir()},
        limitations=['Sample is sky-balanced rather than population weighted.',
            'Rotated controls share the target image; pairs are cluster-level pooled fractions, not one-to-one stamp matches.',
            'Spearman intervals are exploratory, marginal and unadjusted for multiple comparisons.',
            'Mask components overlap; no causal claim or mask relaxation is supported.',
            'PSF, sky, Poisson, covariance and membership gates remain open.'])
    (out/'summary.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--summary-dir',type=Path,required=True);p.add_argument('--plan',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--replicas',type=int,default=5000)
    a=p.parse_args()
    if a.replicas<100:p.error('at least 100 replicates required')
    run(a.summary_dir.resolve(),a.plan.resolve(),a.output.resolve(),a.replicas)
