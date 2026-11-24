#!/usr/bin/env python3
"""Individual-cluster influence on richness-standardized injection recovery."""
import argparse
import json
from pathlib import Path
import numpy as np

from audit_s23b_injection_richness_reweight import INJECTION, RELATION, RELATION_SHA, metrics, read_tsv, sha256, write_tsv

PRIOR=RELATION.parent/'s23b_injection_richness_reweight_2026-11-21_02'


def run(output):
    if output.exists():
        raise FileExistsError(output)
    relation_file=RELATION/'summary.json'
    injection_file=INJECTION/'summary.json'
    prior_file=PRIOR/'summary.json'
    relation,injection,prior=(json.loads(p.read_text()) for p in (relation_file,injection_file,prior_file))
    if sha256(relation_file)!=RELATION_SHA or injection['relation_summary_sha256']!=RELATION_SHA:
        raise ValueError('Changed relation lineage')
    if prior['injection_summary_sha256']!=sha256(injection_file):
        raise ValueError('Changed injection lineage')
    for directory,report,name in ((RELATION,relation,'paired_icl_satellite_profiles.tsv'),
                                  (INJECTION,injection,'matched_injection_concentration.tsv'),
                                  (PRIOR,prior,'richness_standardization.tsv')):
        if sha256(directory/name)!=report['output_sha256'][name]:
            raise ValueError('Changed input table')
    injected=read_tsv(INJECTION/'matched_injection_concentration.tsv')
    lookup={(r['target_id'],r['filter']):r for r in injected}
    if len(lookup)!=len(injected):
        raise ValueError('Duplicate injection')
    sample=[r for r in read_tsv(RELATION/'paired_icl_satellite_profiles.tsv') if r['icl_relative_contrast_njy_arcsec2'] and
            min(float(r['inner_common_support_fraction']),float(r['outer_common_support_fraction']))>=0.1]
    prior_values={r['filter']:r for r in read_tsv(PRIOR/'richness_standardization.tsv')}
    summary_rows=[]
    influence_rows=[]
    for band in ('i','i2'):
        rows=[]
        for r in sample:
            if r['filter']!=band:
                continue
            matched=lookup.get((r['target_id'],band))
            rows.append(dict(target_id=r['target_id'],richness=float(r['log10_richness']),
                             response=float(matched['psf_response']) if matched else None))
        values=np.asarray([r['richness'] for r in rows])
        edges=np.r_[-np.inf,np.quantile(values,[.25,.5,.75]),np.inf]
        base,_=metrics(rows,edges)
        if base is None or not np.isclose(base[1],float(prior_values[band]['reweighted_mean']),rtol=1e-12):
            raise ValueError('Baseline mismatch')
        idx=np.clip(np.searchsorted(edges,values,side='right')-1,0,3)
        count=np.bincount(idx,minlength=4)
        counts=np.bincount([i for i,r in zip(idx,rows,strict=True) if r['response'] is not None],minlength=4)
        supported=[(r,i) for r,i in zip(rows,idx,strict=True) if r['response'] is not None]
        responses=np.asarray([r['response'] for r,_ in supported])
        weights=np.asarray([count[i]/counts[i] for _,i in supported])
        order=np.argsort(responses)
        cumulative=np.cumsum(weights[order])
        median=float(responses[order[np.searchsorted(cumulative,cumulative[-1]/2)]])
        for selected,_ in supported:
            # Keep the full target population and remove only this injection response.
            altered=[dict(r, response=None) if r['target_id']==selected['target_id'] else r for r in rows]
            value,missing=metrics(altered,edges)
            influence_rows.append(dict(filter=band,target_id=selected['target_id'],response=selected['response'],
                                       leave_one_out=value[1] if value else None,
                                       delta=value[1]-base[1] if value else None,unsupported_fraction=missing))
        valid=[r for r in influence_rows if r['filter']==band and r['delta'] is not None]
        worst=max(valid,key=lambda r:abs(r['delta']))
        summary_rows.append(dict(filter=band,n_injected=len(supported),unweighted_mean=base[0],reweighted_mean=base[1],
                                 unweighted_median=float(np.median(responses)),reweighted_median=median,
                                 max_abs_leave_one_out=float(abs(worst['delta'])),most_influential_id=worst['target_id'],
                                 min_response=float(responses.min()),max_response=float(responses.max()),
                                 unsupported_deletions=len(supported)-len(valid)))
    output.mkdir(parents=True)
    a,b=output/'richness_influence_summary.tsv',output/'richness_leave_one_out.tsv'
    write_tsv(a,summary_rows)
    write_tsv(b,influence_rows)
    report=dict(science_measurement=False,calibration_applied=False,input_prior_sha256=sha256(prior_file),
                runner_sha256=sha256(Path(__file__)),output_sha256={a.name:sha256(a),b.name:sha256(b)},
                limitation='Fixed full-sample richness bins, remove one supported injection response '
                           'while preserving the target population. Conditional sensitivity only; '
                           'does not provide uncertainties, completeness or population calibration.')
    (output/'summary.json').write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    run(parser.parse_args().output.resolve())
