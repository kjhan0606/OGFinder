#!/usr/bin/env python3
"""Conditional within-bin bootstrap for the two COSMOS tract comparison."""
import argparse
import csv
import json
from pathlib import Path

import numpy as np

from audit_s16a_cosmos_lens_response import INTAKE, sha, write_table
from ogfkit.wlcluster import Cosmology

BASE=INTAKE.parent
SUPPORT=BASE/'s16a_cosmos_bin_support_2026-11-21_01'
PRIOR=BASE/'s16a_cosmos_tract_common_support_2026-11-21_01'
SEED=2026112102
REPLICAS=2000


def run(output):
    if output.exists():
        raise FileExistsError(output)
    reports={}
    for name,folder in (('intake',INTAKE),('support',SUPPORT),('prior',PRIOR)):
        report=json.loads((folder/'summary.json').read_text())
        for filename,digest in report['output_sha256'].items():
            if sha(folder/filename)!=digest:
                raise ValueError('Changed input artifact')
        reports[name]=report
    if reports['support']['intake_sha256']!=sha(INTAKE/'summary.json') or reports['prior']['input_support_sha256']!=sha(SUPPORT/'summary.json'):
        raise ValueError('Input lineage differs')
    with (INTAKE/'cosmos_join.tsv').open() as stream:
        data=[r for r in csv.DictReader(stream,delimiter='\t') if r['eligible_for_selection_review']=='True']
    if len(data)!=reports['intake']['eligible_rows']:
        raise ValueError('Eligibility mismatch')
    photo,ref,shape,som=np.array([[float(r[k]) for k in ('photo_z','reference_z','shape_weight','som_weight')] for r in data]).T
    tract=np.array([r['tract'] for r in data])
    with (SUPPORT/'photoz_bin_support.tsv').open() as stream:
        bins=list(csv.DictReader(stream,delimiter='\t'))
    with (SUPPORT/'lens_photoz_reweight.tsv').open() as stream:
        lenses=list(csv.DictReader(stream,delimiter='\t'))
    with (PRIOR/'tract_common_support.tsv').open() as stream:
        prior={r['lens_id']:r for r in csv.DictReader(stream,delimiter='\t')}
    edges=np.asarray(reports['support']['edges'])
    cosmo=Cosmology()
    rng=np.random.default_rng(SEED)
    results=[]
    for lens in lenses:
        zl=float(lens['lens_z'])
        target=np.asarray([float(r['observed_weight_fraction']) for r in bins if r['lens_id']==lens['lens_id']])
        if len(target)!=len(edges)-1 or not np.isclose(target.sum(),1):
            raise ValueError('Invalid target histogram')
        prepared={}
        for region in ('9812','9813'):
            mask=(tract==region)&(photo>zl+0.2)
            p,z,w,sw=photo[mask],ref[mask],shape[mask],som[mask]
            chi=float(cosmo.comoving(zl))
            est=1-chi/cosmo.comoving(p)
            true=np.zeros_like(z)
            background=z>zl
            true[background]=1-chi/cosmo.comoving(z[background])
            weights=w*sw*est**2
            prepared[region]=(np.searchsorted(edges,p,side='right')-1,weights,true/est)
        common=(target>0)
        for index,weights,_ in prepared.values():
            common &= np.bincount(index,weights=weights,minlength=len(target))>0
        retained=float(target[common].sum())
        if retained<=0 or not np.isclose(retained,float(prior[lens['lens_id']]['retained_weight_fraction']),atol=1e-12):
            raise ValueError('Common support mismatch')
        normalized=target/retained
        samples={}
        for region,(index,weights,ratio) in prepared.items():
            full=0.
            draws=np.zeros(REPLICAS)
            for k in np.flatnonzero(common):
                chosen=(index==k)&(weights>0)
                w,q=weights[chosen],ratio[chosen]
                if len(w)<1:
                    raise ValueError('Empty common bin')
                full+=normalized[k]*float(w@q/w.sum())
                # Fixed number of objects per original bin avoids silently
                # dropping support; weights are re-fitted in every draw.
                picks=rng.integers(0,len(w),size=(REPLICAS,len(w)))
                sampled_w=w[picks]
                draws+=normalized[k]*np.sum(sampled_w*q[picks],axis=1)/sampled_w.sum(axis=1)
            multiplier=1/full
            if not np.isclose(multiplier,float(prior[lens['lens_id']]['multiplier_'+region]),rtol=1e-12):
                raise ValueError('Prior point estimate differs')
            samples[region]=1/draws
        difference=samples['9812']/samples['9813']-1
        results.append(dict(lens_id=lens['lens_id'],retained_weight_fraction=retained,
                            point_difference=float(prior[lens['lens_id']]['relative_difference']),
                            conditional_sd=float(np.std(difference,ddof=1)),
                            interval_low=float(np.percentile(difference,2.5)),
                            interval_high=float(np.percentile(difference,97.5))))
    output.mkdir(parents=True)
    table=output/'tract_bootstrap.tsv'
    write_table(table,results)
    report=dict(science_measurement=False,calibration_applied=False,n_lenses=len(results),
                seed=SEED,replicas=REPLICAS,input_intake_sha256=sha(INTAKE/'summary.json'),
                input_support_sha256=sha(SUPPORT/'summary.json'),input_prior_sha256=sha(PRIOR/'summary.json'),
                runner_sha256=sha(Path(__file__)),output_sha256={table.name:sha(table)},
                limitation='Conditional independent-object resampling within fixed tract/redshift bins; '
                           'does not include spatial clustering, cosmic variance, reference-z uncertainty, '
                           'observed histogram error, selection mismatch or excluded bins.')
    (output/'summary.json').write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True)
    run(p.parse_args().output.resolve())
