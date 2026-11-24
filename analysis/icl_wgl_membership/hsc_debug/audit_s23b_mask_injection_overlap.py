#!/usr/bin/env python3
"""Describe existing unconvolved injection recovery in the spatial mask sample.

The one injected member is compared with cluster-pooled mask footprints, not
its own footprint. This is an exploratory overlap audit with no new image
reads, noise injection, network request, causal inference or calibration.
"""
import argparse
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path
import numpy as np
from analyze_s23b_member_mask_paired import spearman, block_indices

DATA=Path('/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13')
INJ=DATA/'s23b_camira_member_outer_injection_recovery_2026-11-07_01/injection_responses.tsv'
MASK=DATA/'s23b_member_mask_paired_2026-11-17_01/paired_cluster_fractions.tsv'
EXPECTED={INJ:'b7cc4deeff36dc95dd0fcbd41cef929b49e6524761f7b561586364d4288bfa30',
 MASK:'bd8646be837d978f8bb6ab88ba1c31bc30c80e4b6a0d94462a55f06a3318624f'}


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):
    with p.open(newline='') as f:return list(csv.DictReader(f,delimiter='\t'))


def run(output,replicas):
    if output.exists():raise FileExistsError(output)
    for p,h in EXPECTED.items():
        if sha(p)!=h:raise ValueError('pinned input changed')
    mask=read(MASK); lookup={(r['target_id'],r['filter']):r for r in mask}
    if len(lookup)!=96:raise ValueError('mask sample is not the reviewed 96 unique clusters')
    joined=[];seen=set()
    for r in read(INJ):
        key=(r['target_id'],r['filter'])
        if key not in lookup:continue
        amplitude=float(r['amplitude_at_50_kpc_njy_arcsec2'])
        unique=key+(amplitude,)
        if unique in seen:raise ValueError('duplicate injection frame/amplitude')
        seen.add(unique);m=lookup[key]
        row=dict(target_id=key[0],filter=key[1],amplitude=amplitude,
            source_row=r['source_row'],injected_status=r['injected_status'],
            sky_cell_x=int(m['sky_cell_x']),sky_cell_y=int(m['sky_cell_y']),
            quality_difference=float(m['difference_quality_bits']),
            unmasked_difference=float(m['difference_unmasked']),
            injected_response=float(r['recovered_increment_over_nominal_input']) if r['recovered_increment_over_nominal_input'] else None,
            masked_model_transfer=float(r['masked_model_over_nominal_input']) if r['masked_model_over_nominal_input'] else None,
            science_measurement=False)
        if (row['injected_status']=='supported') != (row['injected_response'] is not None):
            raise ValueError('response/support mismatch')
        joined.append(row)
    if len(joined)!=34*3:raise ValueError('unexpected overlap: inspect selection provenance')
    records=[]
    for band in ('i','i2'):
        for amplitude in (.5,1.,2.):
            allrows=[r for r in joined if r['filter']==band and r['amplitude']==amplitude]
            supported=[r for r in allrows if r['injected_status']=='supported']
            cells=[(r['sky_cell_x'],r['sky_cell_y']) for r in supported]
            values=np.array([r['injected_response'] for r in supported])
            records.append(dict(filter=band,amplitude=amplitude,n_selected=len(allrows),n_supported=len(supported),
                n_cells=len(set(cells)),statuses=dict(Counter(r['injected_status'] for r in allrows)),
                median_response=float(np.median(values)) if len(values) else None))
    correlations=[]
    for band in ('i','i2'):
        sample=[r for r in joined if r['filter']==band and r['amplitude']==1. and r['injected_status']=='supported']
        y=np.array([r['injected_response'] for r in sample]);cells=[(r['sky_cell_x'],r['sky_cell_y']) for r in sample]
        draws=block_indices(cells,replicas,20261117)
        for metric in ('quality_difference','unmasked_difference'):
            x=np.array([r[metric] for r in sample]);estimate=spearman(x,y)
            valid=[v for idx in draws if (v:=spearman(x[idx],y[idx])) is not None]
            low,high=np.quantile(valid,[.025,.975]) if len(valid)>=.9*replicas else (None,None)
            correlations.append(dict(filter=band,metric=metric,n_clusters=len(sample),n_cells=len(set(cells)),
                rho=estimate,ci_low=None if low is None else float(low),ci_high=None if high is None else float(high),
                valid_draws=len(valid),confidence_status='descriptive' if low is not None else 'insufficient_valid_draws'))
    output.mkdir(parents=True,exist_ok=False)
    with (output/'overlap_ledger.tsv').open('x',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(joined[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(joined)
    result=dict(status='complete_exploratory_unconvolved_overlap',science_measurement=False,network_used=False,images_read=False,
        input_sha256={str(p):sha(p) for p in EXPECTED},script_sha256=sha(Path(__file__)),
        output_sha256={'overlap_ledger.tsv':sha(output/'overlap_ledger.tsv')},
        baseline_reference_amplitude_njy_arcsec2=1.,bootstrap_replicas=replicas,seed=20261117,
        support_by_filter_amplitude=records,exploratory_correlations=correlations,
        limitations=['Post hoc overlapping sky-balanced sample; no population inference.',
          'Mask predictors pool all selected cluster members, not the injected member footprint.',
          'Only previously supported profiles enter recovery: excluded candidates are not completeness.',
          'No injected PSF, source Poisson, coadd covariance or full sky calibration.',
          'Correlation intervals are marginal, unadjusted and noncausal.'])
    (output/'summary.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n')
    print(json.dumps(result,indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True)
    p.add_argument('--replicas',type=int,default=5000);a=p.parse_args()
    if a.replicas<100:p.error('at least 100 replicates required')
    run(a.output.resolve(),a.replicas)
