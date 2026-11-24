#!/usr/bin/env python3
"""Compare the existing injection subset with the remaining ICL pairs."""
import argparse
import json
from pathlib import Path

import numpy as np

from audit_s23b_icl_satellite_field_omission import DATA, RELATION, RELATION_SHA, read_tsv, sha256, write_tsv

INJECTION=DATA/'s23b_injection_support_selection_2026-11-21_01'


def run(output):
    if output.exists():
        raise FileExistsError(output)
    relation_file=RELATION/'summary.json'
    if sha256(relation_file)!=RELATION_SHA:
        raise ValueError('Changed relation summary')
    relation=json.loads(relation_file.read_text())
    relation_table=RELATION/'paired_icl_satellite_profiles.tsv'
    inj_summary=INJECTION/'summary.json'
    inj=json.loads(inj_summary.read_text())
    injection_table=INJECTION/'matched_support_ledger.tsv'
    if sha256(relation_table)!=relation['output_sha256'][relation_table.name] or sha256(injection_table)!=inj['output_sha256'][injection_table.name]:
        raise ValueError('Changed input table')
    matched=read_tsv(injection_table)
    selected={(r['target_id'],r['filter']) for r in matched}
    if len(selected)!=len(matched):
        raise ValueError('Duplicate injection target')
    rows=[r for r in read_tsv(relation_table)
          if r['icl_relative_contrast_njy_arcsec2'] and
          min(float(r['inner_common_support_fraction']),float(r['outer_common_support_fraction']))>=0.1]
    all_ids={(r['target_id'],r['filter']) for r in rows}
    if len(all_ids)!=len(rows) or not selected.issubset(all_ids):
        raise ValueError('Injection membership mismatch')
    results=[]
    fields=('log10_inner_outer_candidate_density_ratio','z_cl','log10_richness')
    for band_index,band in enumerate(('i','i2')):
        sample=[r for r in rows if r['filter']==band]
        cells={}
        for r in sample:
            cell=(int(np.floor(((float(r['ra_deg'])-84.056247280625)%360)/5)),
                  int(np.floor((float(r['dec_deg'])+90)/5)))
            cells.setdefault(cell,[]).append(r)
        keys=sorted(cells)
        for metric_index,metric in enumerate(fields):
            def contrast(group):
                yes=[float(r[metric]) for r in group if (r['target_id'],r['filter']) in selected]
                no=[float(r[metric]) for r in group if (r['target_id'],r['filter']) not in selected]
                return float(np.mean(yes)-np.mean(no)) if yes and no else None
            included=[float(r[metric]) for r in sample if (r['target_id'],r['filter']) in selected]
            excluded=[float(r[metric]) for r in sample if (r['target_id'],r['filter']) not in selected]
            rng=np.random.default_rng(2026112130+band_index*10+metric_index)
            draws=[]
            for _ in range(5000):
                boot=[r for i in rng.integers(0,len(keys),len(keys)) for r in cells[keys[i]]]
                value=contrast(boot)
                if value is not None:
                    draws.append(value)
            if len(draws)<4500:
                raise ValueError('Insufficient bootstrap draws')
            low,high=np.percentile(draws,[2.5,97.5])
            results.append(dict(filter=band,metric=metric,n_injection=len(included),n_other=len(excluded),
                                n_cells=len(keys),injection_mean=float(np.mean(included)),
                                other_mean=float(np.mean(excluded)),mean_difference=contrast(sample),
                                interval_low=float(low),interval_high=float(high),valid_draws=len(draws)))
    output.mkdir(parents=True)
    table=output/'injection_sample_comparison.tsv'
    write_tsv(table,results)
    report=dict(science_measurement=False,relation_summary_sha256=RELATION_SHA,
                injection_summary_sha256=sha256(inj_summary),runner_sha256=sha256(Path(__file__)),
                output_sha256={table.name:sha256(table)},replicas=5000,
                limitation='Post-hoc mean comparisons with marginal 5-degree target-sky-cell bootstrap intervals. '
                           'Comparing observed covariates cannot prove transfer representativeness, '
                           'handle all selection variables or justify a population correction.')
    (output/'summary.json').write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    run(parser.parse_args().output.resolve())
