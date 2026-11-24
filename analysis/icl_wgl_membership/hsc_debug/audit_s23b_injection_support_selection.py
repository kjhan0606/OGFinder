#!/usr/bin/env python3
"""Compare injection support across satellite concentration quartiles."""
import argparse
from collections import Counter
import json
from pathlib import Path

import numpy as np

from audit_s23b_injection_concentration import COMPARE, COMPARE_SHA
from audit_s23b_icl_satellite_field_omission import RELATION, RELATION_SHA, read_tsv, sha256, write_tsv


def run(output, amplitude=1.0):
    if output.exists():
        raise FileExistsError(output)
    rsummary,csummary=RELATION/'summary.json',COMPARE/'summary.json'
    if sha256(rsummary)!=RELATION_SHA or sha256(csummary)!=COMPARE_SHA:
        raise ValueError('Changed pinned summaries')
    rel_meta,inj_meta=(json.loads(p.read_text()) for p in (rsummary,csummary))
    rel_table=RELATION/'paired_icl_satellite_profiles.tsv'
    inj_table=COMPARE/'paired_ledger.tsv'
    if sha256(rel_table)!=rel_meta['output_sha256'][rel_table.name] or sha256(inj_table)!=inj_meta['output_sha256'][inj_table.name]:
        raise ValueError('Changed source table')
    rel=read_tsv(rel_table)
    lookup={(r['target_id'],r['filter']):r for r in rel}
    if len(lookup)!=len(rel):
        raise ValueError('Duplicate relation entries')
    inject=[r for r in read_tsv(inj_table) if float(r['amplitude'])==amplitude]
    if len(inject)!=427 or len({(r['target_id'],r['filter']) for r in inject})!=427:
        raise ValueError('Unexpected injection sample')
    matched=[]
    for item in inject:
        row=lookup.get((item['target_id'],item['filter']))
        if row is None or not row['icl_relative_contrast_njy_arcsec2']:
            continue
        if min(float(row['inner_common_support_fraction']),float(row['outer_common_support_fraction']))<0.1:
            continue
        matched.append(dict(target_id=item['target_id'],filter=item['filter'],
                            concentration=float(row['log10_inner_outer_candidate_density_ratio']),
                            supported=item['both_supported']=='True',psf_status=item['psf_status'],
                            unconvolved_status=item['unconvolved_status'],
                            ra=float(row['ra_deg']),dec=float(row['dec_deg'])))
    results=[]
    for band_index,band in enumerate(('i','i2')):
        rows=[r for r in matched if r['filter']==band]
        values=np.asarray([r['concentration'] for r in rows])
        q1,q3=np.quantile(values,[0.25,0.75])
        low=[r for r in rows if r['concentration']<=q1]
        high=[r for r in rows if r['concentration']>=q3]
        if set(r['target_id'] for r in low)&set(r['target_id'] for r in high):
            raise ValueError('Overlapping quartiles')
        for label,group in (('low',low),('high',high)):
            for r in group:
                r['quartile']=label
        cells={}
        for r in low+high:
            cell=(int(np.floor(((r['ra']-84.056247280625)%360)/5)),
                  int(np.floor((r['dec']+90)/5)))
            cells.setdefault(cell,[]).append(r)
        keys=sorted(cells)
        def difference(sample):
            a=[r['supported'] for r in sample if r['quartile']=='low']
            b=[r['supported'] for r in sample if r['quartile']=='high']
            return float(np.mean(b)-np.mean(a)) if a and b else None
        point=difference(low+high)
        rng=np.random.default_rng(2026112120+band_index)
        draws=[]
        for _ in range(5000):
            sample=[r for i in rng.integers(0,len(keys),len(keys)) for r in cells[keys[i]]]
            value=difference(sample)
            if value is not None:
                draws.append(value)
        if len(draws)<4500:
            raise ValueError('Insufficient valid cell bootstrap draws')
        lo,hi=np.percentile(draws,[2.5,97.5])
        results.append(dict(filter=band,n_matched=len(rows),n_cells=len(keys),
                            n_low=len(low),n_high=len(high),
                            low_supported=sum(r['supported'] for r in low),
                            high_supported=sum(r['supported'] for r in high),
                            low_fraction=float(np.mean([r['supported'] for r in low])),
                            high_fraction=float(np.mean([r['supported'] for r in high])),
                            high_minus_low=point,low95=float(lo),high95=float(hi),
                            valid_draws=len(draws),
                            status_counts=str(dict(Counter((r['unconvolved_status'],r['psf_status']) for r in rows)))))
    output.mkdir(parents=True)
    table=output/'support_by_concentration.tsv'
    ledger=output/'matched_support_ledger.tsv'
    write_tsv(table,results)
    write_tsv(ledger,matched)
    report=dict(science_measurement=False,amplitude_at_50_kpc_njy_arcsec2=amplitude,relation_summary_sha256=RELATION_SHA,
                injection_comparison_sha256=COMPARE_SHA,runner_sha256=sha256(Path(__file__)),
                output_sha256={p.name:sha256(p) for p in (table,ledger)},
                definition='Both PSF and unconvolved injection measurements pass their existing profile-selection criteria.',
                limitation='Post-hoc extreme quartiles of injection sample conditional on 10% ICL support; '
                           'one injected member per cluster, unadjusted association and marginal sky-cell intervals. '
                           'Not a population completeness measurement.')
    (output/'summary.json').write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--amplitude',type=float,choices=(0.5,1.0,2.0),default=1.0)
    args=parser.parse_args()
    run(args.output.resolve(),args.amplitude)
