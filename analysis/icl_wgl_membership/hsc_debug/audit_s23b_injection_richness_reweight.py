#!/usr/bin/env python3
"""Sensitivity of injected-member mean recovery to richness standardization."""
import argparse
import json
from pathlib import Path

import numpy as np

from audit_s23b_icl_satellite_field_omission import DATA, RELATION, RELATION_SHA, read_tsv, sha256, write_tsv

INJECTION=DATA/'s23b_injection_richness_2026-11-21_01'


def metrics(rows, edges):
    x=np.asarray([r['richness'] for r in rows])
    bins=np.searchsorted(edges,x,side='right')-1
    bins=np.clip(bins,0,len(edges)-2)
    target=np.bincount(bins,minlength=len(edges)-1)
    included=np.asarray([r['response'] is not None for r in rows])
    count=np.bincount(bins[included],minlength=len(edges)-1)
    missing=float(target[count==0].sum()/target.sum())
    if missing>0 or not included.any():
        return None,missing
    factors=np.divide(target,count,out=np.zeros_like(target,dtype=float),where=count>0)
    weights=factors[bins[included]]
    response=np.asarray([r['response'] for r in rows if r['response'] is not None])
    raw=float(np.mean(response))
    weighted=float(weights@response/weights.sum())
    return (raw,weighted,float(weights.sum()**2/(weights@weights)),int(count.min())),missing


def run(output):
    if output.exists():
        raise FileExistsError(output)
    relation_path=RELATION/'summary.json'
    if sha256(relation_path)!=RELATION_SHA:
        raise ValueError('Changed relation summary')
    relation=json.loads(relation_path.read_text())
    inj_path=INJECTION/'summary.json'
    inj=json.loads(inj_path.read_text())
    if inj['relation_summary_sha256']!=RELATION_SHA or inj['science_measurement'] is not False:
        raise ValueError('Injection lineage or status differs')
    rt=RELATION/'paired_icl_satellite_profiles.tsv'
    it=INJECTION/'matched_injection_concentration.tsv'
    if sha256(rt)!=relation['output_sha256'][rt.name] or sha256(it)!=inj['output_sha256'][it.name]:
        raise ValueError('Changed input table')
    injections=read_tsv(it)
    lookup={(r['target_id'],r['filter']):r for r in injections}
    if len(lookup)!=len(injections):
        raise ValueError('Duplicate injection target')
    sample=[r for r in read_tsv(rt) if r['icl_relative_contrast_njy_arcsec2'] and
            min(float(r['inner_common_support_fraction']),float(r['outer_common_support_fraction']))>=0.1]
    results=[]
    bins_output=[]
    for band_index,band in enumerate(('i','i2')):
        rows=[]
        for row in sample:
            if row['filter']!=band:
                continue
            injected=lookup.get((row['target_id'],band))
            rows.append(dict(target_id=row['target_id'],richness=float(row['log10_richness']),
                             response=float(injected['psf_response']) if injected else None,
                             ra=float(row['ra_deg']),dec=float(row['dec_deg'])))
        values=np.asarray([r['richness'] for r in rows])
        edges=np.r_[-np.inf,np.quantile(values,[.25,.5,.75]),np.inf]
        if np.any(np.diff(edges)<=0):
            raise ValueError('Duplicate quartile edges')
        original,missing=metrics(rows,edges)
        if original is None:
            raise ValueError('Observed richness bin lacks supported injection')
        cells={}
        for row in rows:
            cell=(int(np.floor(((row['ra']-84.056247280625)%360)/5)),int(np.floor((row['dec']+90)/5)))
            cells.setdefault(cell,[]).append(row)
        keys=sorted(cells)
        rng=np.random.default_rng(2026112140+band_index)
        draws=[]
        for _ in range(5000):
            resampled=[r for i in rng.integers(0,len(keys),len(keys)) for r in cells[keys[i]]]
            value,_=metrics(resampled,edges)
            if value is not None:
                draws.append(value[1]-value[0])
        if len(draws)<4500:
            raise ValueError('Too few supported resamples')
        lo,hi=np.percentile(draws,[2.5,97.5])
        results.append(dict(filter=band,n_target=len(rows),n_injected=sum(r['response'] is not None for r in rows),
                            n_cells=len(keys),unweighted_mean=original[0],reweighted_mean=original[1],
                            difference=original[1]-original[0],interval_low=float(lo),interval_high=float(hi),
                            reweighted_neff=original[2],minimum_injected_per_bin=original[3],
                            unsupported_target_fraction=missing,valid_draws=len(draws)))
        assigned=np.clip(np.searchsorted(edges,values,side='right')-1,0,len(edges)-2)
        for k in range(4):
            group=[r for r,i in zip(rows,assigned,strict=True) if i==k]
            bins_output.append(dict(filter=band,bin=k,n_target=len(group),
                                    n_injected=sum(r['response'] is not None for r in group),
                                    lower=None if not np.isfinite(edges[k]) else float(edges[k]),
                                    upper=None if not np.isfinite(edges[k+1]) else float(edges[k+1])))
    output.mkdir(parents=True)
    a,b=output/'richness_standardization.tsv',output/'richness_bin_support.tsv'
    write_tsv(a,results)
    write_tsv(b,bins_output)
    report=dict(science_measurement=False,calibration_applied=False,relation_summary_sha256=RELATION_SHA,
                injection_summary_sha256=sha256(inj_path),runner_sha256=sha256(Path(__file__)),
                output_sha256={a.name:sha256(a),b.name:sha256(b)},replicas=5000,
                method='Reweight supported 1 nJy PSF-core injection responses to the full 10%-support '
                       'cluster sample in four fixed full-sample richness quantiles. '
                       'Mean contrasts use shared target-sky-cell bootstrap.',
                limitation='Post-hoc richness-only mean sensitivity. Supported one-member injections, '
                           'finite cells, selection, redshift, PSF wings and other covariates prevent '
                           'population calibration or completeness inference.')
    (output/'summary.json').write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    run(parser.parse_args().output.resolve())
