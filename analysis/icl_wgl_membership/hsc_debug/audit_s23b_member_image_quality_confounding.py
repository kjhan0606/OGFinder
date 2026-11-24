#!/usr/bin/env python3
"""Audit image-quality predictability of the satellite concentration."""
import argparse
import json
from pathlib import Path

import numpy as np

from audit_s23b_icl_satellite_image_quality import CONTROL, CONTROL_SHA
from audit_s23b_icl_satellite_field_omission import (
    ORIGIN, RELATION, RELATION_SHA, connected_pair_blocks, model_matrix,
    read_tsv, sha256, write_tsv,
)


def residual_variance(matrix, y, weights):
    use=weights>0
    roots=np.sqrt(weights[use])
    x=matrix[use]*roots[:,None]
    yw=y[use]*roots
    coef,_,rank,singular=np.linalg.lstsq(x,yw,rcond=None)
    if rank!=x.shape[1] or singular[-1]<=singular[0]*1e-5:
        return None
    residual=yw-x@coef
    return float(residual@residual)


def run(output):
    if output.exists():
        raise FileExistsError(output)
    summary=RELATION/'summary.json'
    if sha256(summary)!=RELATION_SHA or sha256(CONTROL)!=CONTROL_SHA:
        raise ValueError('Changed input summary')
    table=RELATION/'paired_icl_satellite_profiles.tsv'
    if sha256(table)!=json.loads(summary.read_text())['output_sha256'][table.name]:
        raise ValueError('Changed relation rows')
    controls={r['id']:r for r in json.loads(CONTROL.read_text())['clusters']}
    rows=[r for r in read_tsv(table) if r['icl_relative_contrast_njy_arcsec2']]
    if len(rows)!=988:
        raise ValueError('Unexpected sample')
    results=[]
    folds=[]
    for band_index,band in enumerate(('i','i2')):
        sample=[r for r in rows if r['filter']==band and min(float(r['inner_common_support_fraction']),float(r['outer_common_support_fraction']))>=0.1]
        X,_,_=model_matrix(sample)
        predictor=X[:,1]
        reduced=np.delete(X,1,axis=1)
        quality=[]
        for row in sample:
            c=controls[row['control_id']]
            if c['target_cluster_id']!=row['target_id']:
                raise ValueError('Changed pairing')
            quality.append([float(c[k]) for k in ('target_depth_mag','delta_depth_mag','target_seeing_arcsec','seeing_ratio')])
        quality=np.asarray(quality)
        if not np.all(np.isfinite(quality)) or np.any(quality.std(axis=0)==0):
            raise ValueError('Invalid quality input')
        expanded=np.column_stack((reduced,(quality-quality.mean(axis=0))/quality.std(axis=0)))
        def partial_r2(weights):
            old=residual_variance(reduced,predictor,weights)
            new=residual_variance(expanded,predictor,weights)
            if old is None or new is None or old<=0:
                return None
            return float(1-new/old)
        value=partial_r2(np.ones(len(sample)))
        control={r['target_id']:dict(ra_deg=r['control_ra_deg'],dec_deg=r['control_dec_deg']) for r in sample}
        labels,blocks=connected_pair_blocks(sample,control,5.,0.,ORIGIN)
        unique,group=np.unique(labels,return_inverse=True)
        if len(unique)<expanded.shape[1]+2:
            raise ValueError('Too few independent blocks')
        total_old=0.
        total_new=0.
        for block in range(len(unique)):
            train=group!=block
            test=~train
            residuals=[]
            for design in (reduced,expanded):
                beta,_,rank,singular=np.linalg.lstsq(design[train],predictor[train],rcond=None)
                if rank!=design.shape[1] or singular[-1]<=singular[0]*1e-5:
                    raise ValueError('Cross-validation training fold is rank deficient')
                err=predictor[test]-design[test]@beta
                residuals.append(float(err@err))
            total_old+=residuals[0]
            total_new+=residuals[1]
            folds.append(dict(filter=band,block=str(unique[block]),n_test=int(test.sum()),
                              baseline_sse=residuals[0],expanded_sse=residuals[1]))
        if total_old<=0:
            raise ValueError('Empty cross-validation residual')
        draws=np.random.default_rng(2026112105+band_index).multinomial(len(unique),np.full(len(unique),1/len(unique)),size=5000)
        estimates=[]
        for d in draws:
            result=partial_r2(d[group])
            if result is not None:
                estimates.append(result)
        if len(estimates)<4500:
            raise ValueError('Too few valid resamples')
        low,high=np.percentile(estimates,[2.5,97.5])
        results.append(dict(filter=band,n_pairs=len(sample),n_blocks=len(blocks),
                            partial_r_squared=value,low=float(low),high=float(high),valid_draws=len(estimates),
                            holdout_baseline_sse=total_old,holdout_expanded_sse=total_new,
                            holdout_fractional_sse_reduction=1-total_new/total_old))
    output.mkdir(parents=True)
    target=output/'member_quality_partial_r2.tsv'
    write_tsv(target,results)
    fold_table=output/'member_quality_holdout_folds.tsv'
    write_tsv(fold_table,folds)
    report=dict(science_measurement=False,control_manifest_sha256=CONTROL_SHA,relation_summary_sha256=RELATION_SHA,
                runner_sha256=sha256(Path(__file__)),output_sha256={target.name:sha256(target),fold_table.name:sha256(fold_table)},
                method='Incremental weighted residual R-squared for four image-quality variables '
                       'predicting member concentration after existing model covariates. '
                       'One 5-degree connected-block layout, 5000 bootstrap draws and leave-one-block-out '
                       'held-out squared prediction error for reduced and expanded models.',
                limitation='Exploratory spatial holdout and in-sample association, not causal confounding. '
                           'Adjacent blocks may remain correlated and only one grid layout is tested.')
    (output/'summary.json').write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    run(parser.parse_args().output.resolve())
