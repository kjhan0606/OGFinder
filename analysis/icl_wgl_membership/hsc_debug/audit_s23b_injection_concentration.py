#!/usr/bin/env python3
"""Compare PSF-core injection recovery with cluster member concentration."""
import argparse
import json
from pathlib import Path

import numpy as np
from scipy.stats import rankdata

from audit_s23b_icl_satellite_field_omission import (
    DATA, RELATION, RELATION_SHA, read_tsv, sha256, write_tsv,
)

COMPARE=DATA/'s23b_member_psf_injection_run_2026-11-17_01/comparison'
COMPARE_SHA='3ecb1ad6900ac6255b4143b461782a3978e020e3dbed9e22b9209e534bec19cf'


def correlation(x,y):
    a=rankdata(x);b=rankdata(y)
    a-=a.mean();b-=b.mean()
    denom=np.sqrt((a@a)*(b@b))
    return float((a@b)/denom) if denom>0 else None


def run(output,predictor='concentration'):
    if output.exists():
        raise FileExistsError(output)
    paths=(RELATION/'summary.json',COMPARE/'summary.json')
    if sha256(paths[0])!=RELATION_SHA or sha256(paths[1])!=COMPARE_SHA:
        raise ValueError('Changed source summary')
    relation,compare=(json.loads(path.read_text()) for path in paths)
    rel_table=RELATION/'paired_icl_satellite_profiles.tsv'
    inject_table=COMPARE/'paired_ledger.tsv'
    for table,report in ((rel_table,relation),(inject_table,compare)):
        if sha256(table)!=report['output_sha256'][table.name]:
            raise ValueError('Changed source table')
    rel=read_tsv(rel_table)
    lookup={(r['target_id'],r['filter']):r for r in rel}
    if len(lookup)!=len(rel):
        raise ValueError('Duplicate relation target')
    selected=[r for r in read_tsv(inject_table) if float(r['amplitude'])==1.0 and r['both_supported']=='True']
    if len({(r['target_id'],r['filter']) for r in selected})!=len(selected):
        raise ValueError('Multiple injected members per cluster')
    joined=[]
    for item in selected:
        r=lookup.get((item['target_id'],item['filter']))
        if r is None or not r['icl_relative_contrast_njy_arcsec2']:
            continue
        joined.append(dict(target_id=item['target_id'],filter=item['filter'],
                           concentration=float(r['log10_inner_outer_candidate_density_ratio']),
                           richness=float(r['log10_richness']),
                           common_support=min(float(r['inner_common_support_fraction']),float(r['outer_common_support_fraction'])),
                           psf_response=float(item['psf_response']),delta_response=float(item['delta_response']),
                           ra_deg=float(r['ra_deg']),dec_deg=float(r['dec_deg'])))
    results=[]
    for band_index,band in enumerate(('i','i2')):
        for floor_index,floor in enumerate((0.,0.1)):
            sample=[r for r in joined if r['filter']==band and r['common_support']>=floor]
            cells={}
            for r in sample:
                cell=(int(np.floor(((r['ra_deg']-84.056247280625)%360)/5)),
                      int(np.floor((r['dec_deg']+90)/5)))
                cells.setdefault(cell,[]).append(r)
            keys=sorted(cells)
            if len(keys)<5:
                raise ValueError('Too few sky cells')
            for metric_index,metric in enumerate(('psf_response','delta_response')):
                estimate=correlation([r[predictor] for r in sample],[r[metric] for r in sample])
                rng=np.random.default_rng(2026112110+band_index*10+floor_index*2+metric_index)
                draws=[]
                for _ in range(5000):
                    members=[r for index in rng.integers(0,len(keys),len(keys)) for r in cells[keys[index]]]
                    value=correlation([r[predictor] for r in members],[r[metric] for r in members])
                    if value is not None:
                        draws.append(value)
                if len(draws)<4500 or estimate is None:
                    raise ValueError('Insufficient correlation draws')
                lo,hi=np.percentile(draws,[2.5,97.5])
                results.append(dict(filter=band,support_floor=floor,metric=metric,n_clusters=len(sample),
                                    n_cells=len(keys),spearman=estimate,ci_low=float(lo),ci_high=float(hi),valid_draws=len(draws)))
    output.mkdir(parents=True)
    ledger=output/'matched_injection_concentration.tsv'
    result=output/'injection_concentration_correlation.tsv'
    write_tsv(ledger,joined)
    write_tsv(result,results)
    report=dict(science_measurement=False,predictor=predictor,relation_summary_sha256=RELATION_SHA,comparison_summary_sha256=COMPARE_SHA,
                runner_sha256=sha256(Path(__file__)),output_sha256={p.name:sha256(p) for p in (ledger,result)},
                limitation='Post-hoc paired-supported injections, one member per cluster, one halo shape and '
                           'amplitude. Unadjusted Spearman intervals are marginal, no multiplicity correction; '
                           'selection, PSF wings and spatial cosmic variance remain unresolved.')
    (output/'summary.json').write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--predictor',choices=('concentration','richness'),default='concentration')
    args=parser.parse_args()
    run(args.output.resolve(),args.predictor)
