#!/usr/bin/env python3
"""Paired spatial bootstrap of ICL–satellite image-quality adjustment."""
import argparse
import json
from pathlib import Path

import numpy as np

from audit_s23b_icl_satellite_image_quality import CONTROL, CONTROL_SHA
from audit_s23b_icl_satellite_field_omission import (
    DATA, FLOORS, N_BOOTSTRAP, ORIGIN, RELATION, RELATION_SHA,
    connected_pair_blocks, model_matrix, read_tsv, sha256, write_tsv,
)

PREVIOUS=DATA/'s23b_icl_satellite_image_quality_2026-11-21_01'
PREVIOUS_SHA=''


def fit(design,y,weights):
    mask=weights>0
    matrix=design[mask]*np.sqrt(weights[mask,None])
    response=y[mask]*np.sqrt(weights[mask])
    beta,_,rank,singular=np.linalg.lstsq(matrix,response,rcond=None)
    if rank!=design.shape[1] or singular[-1]<=singular[0]*1e-5:
        return None
    return float(beta[1])


def run(output):
    if output.exists():
        raise FileExistsError(output)
    relation_path=RELATION/'summary.json'
    if sha256(relation_path)!=RELATION_SHA or sha256(CONTROL)!=CONTROL_SHA:
        raise ValueError('Changed pinned input')
    relation=json.loads(relation_path.read_text())
    table=RELATION/'paired_icl_satellite_profiles.tsv'
    if sha256(table)!=relation['output_sha256'][table.name]:
        raise ValueError('Changed relation table')
    controls={c['id']:c for c in json.loads(CONTROL.read_text())['clusters']}
    rows=[r for r in read_tsv(table) if r['icl_relative_contrast_njy_arcsec2']]
    if len(rows)!=988:
        raise ValueError('Unexpected supported sample')
    for row in rows:
        c=controls[row['control_id']]
        if c['target_cluster_id']!=row['target_id']:
            raise ValueError('Target mismatch')
        row['quality']=[float(c[k]) for k in ('target_depth_mag','delta_depth_mag','target_seeing_arcsec','seeing_ratio')]
    prev=json.loads((PREVIOUS/'summary.json').read_text())
    prior_file=PREVIOUS/'quality_covariate_sensitivity.tsv'
    if sha256(prior_file)!=prev['output_sha256'][prior_file.name]:
        raise ValueError('Changed prior quality table')
    prior={(r['filter'],float(r['support_floor']),r['model']):r for r in read_tsv(prior_file)}
    results=[]
    for bi,band in enumerate(('i','i2')):
        for fi,floor in enumerate(FLOORS):
            sample=[r for r in rows if r['filter']==band and min(float(r['inner_common_support_fraction']),float(r['outer_common_support_fraction']))>=floor]
            X,y,_=model_matrix(sample)
            quality=np.asarray([r['quality'] for r in sample])
            if not np.all(np.isfinite(quality)) or np.any(quality.std(axis=0)==0):
                raise ValueError('Invalid quality columns')
            quality=(quality-quality.mean(axis=0))/quality.std(axis=0)
            extended=np.column_stack((X,quality))
            reference=fit(X,y,np.ones(len(y)))
            updated=fit(extended,y,np.ones(len(y)))
            if reference is None or updated is None:
                raise ValueError('Full model rank deficient')
            for key,value in (('baseline',reference),('image_quality',updated)):
                saved=prior[band,floor,key]
                if saved['status']=='complete' and not np.isclose(value,float(saved['slope']),rtol=0,atol=1e-9):
                    raise ValueError('Prior point estimate differs')
            control={r['target_id']:dict(ra_deg=r['control_ra_deg'],dec_deg=r['control_dec_deg']) for r in sample}
            for wi,width in enumerate((2.5,5.,10.)):
                for pi,phase in enumerate((0.,0.5)):
                    labels,blocks=connected_pair_blocks(sample,control,width,phase,ORIGIN)
                    unique,row_group=np.unique(labels,return_inverse=True)
                    rec=dict(filter=band,support_floor=floor,width_deg=width,phase=phase,n_pairs=len(sample),
                             n_blocks=len(blocks),baseline=reference,adjusted=updated,
                             delta=updated-reference,low=None,high=None,valid_draws=0,status='insufficient_blocks')
                    if len(unique)>=extended.shape[1]+2:
                        rng=np.random.default_rng(2026112200+bi*100+fi*10+wi*2+pi)
                        draws=rng.multinomial(len(unique),np.full(len(unique),1/len(unique)),size=N_BOOTSTRAP)
                        changes=[]
                        for draw in draws:
                            weights=draw[row_group]
                            a=fit(X,y,weights)
                            b=fit(extended,y,weights)
                            if a is not None and b is not None:
                                changes.append(b-a)
                        rec['valid_draws']=len(changes)
                        rec['status']='insufficient_valid_draws'
                        if len(changes)>=0.9*N_BOOTSTRAP:
                            rec['low'],rec['high']=map(float,np.percentile(changes,[2.5,97.5]))
                            rec['status']='complete'
                    results.append(rec)
    output.mkdir(parents=True)
    dest=output/'quality_paired_bootstrap.tsv'
    write_tsv(dest,results)
    summary=dict(science_measurement=False,replicas=N_BOOTSTRAP,n_results=len(results),
                 relation_sha256=RELATION_SHA,control_manifest_sha256=CONTROL_SHA,
                 previous_summary_sha256=sha256(PREVIOUS/'summary.json'),runner_sha256=sha256(Path(__file__)),
                 output_sha256={dest.name:sha256(dest)},
                 limitation='Exploratory paired connected-block percentiles on fixed samples. '
                            'Selection, PSF wings, pixel covariance and causal interpretation are unresolved.')
    (output/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    run(parser.parse_args().output.resolve())
