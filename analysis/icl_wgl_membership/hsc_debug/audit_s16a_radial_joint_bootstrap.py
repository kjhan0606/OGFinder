#!/usr/bin/env python3
"""Joint Poisson bootstrap of observed sources and shared COSMOS references."""
import argparse
from collections import defaultdict
import csv
import json
from pathlib import Path

import numpy as np

from audit_s16a_radial_photoz_response import INTAKE, SUPPORT, CALIB, GEOMETRY, GEOMETRY_SHA, PLAN, OBJECT_QA, finite_point_redshifts, load_inputs, read, sha, write_table
from ogfkit.wlcluster import Cosmology

PREVIOUS=INTAKE.parent/'s16a_radial_photoz_response_2026-11-21_01'
REPLICAS=2000
SEED=2026112104


def run(output):
    if output.exists():
        raise FileExistsError(output)
    reports={}
    for folder in (INTAKE,SUPPORT,CALIB,PREVIOUS):
        report=json.loads((folder/'summary.json').read_text())
        for name,digest in report['output_sha256'].items():
            if sha(folder/name)!=digest:
                raise ValueError('Changed source product')
        reports[folder]=report
    if reports[CALIB]['intake_sha256']!=sha(INTAKE/'summary.json') or reports[PREVIOUS]['support_sha256']!=sha(SUPPORT/'summary.json') or reports[PREVIOUS]['calibration_sha256']!=sha(CALIB/'summary.json'):
        raise ValueError('Different source lineage')
    centers,pairs,selected,available,_,provenance=load_inputs(PLAN,OBJECT_QA)
    for key,value in provenance.items():
        if reports[SUPPORT][key]!=value or reports[CALIB][key]!=value:
            raise ValueError('Different pilot lineage')
    point=finite_point_redshifts(available)
    if sha(GEOMETRY/'summary.json')!=GEOMETRY_SHA:
        raise ValueError('Changed geometry summary')
    geom=json.loads((GEOMETRY/'summary.json').read_text())
    geom_table=GEOMETRY/'pair_geometry.tsv'
    if sha(geom_table)!=geom['output_sha256'][geom_table.name]:
        raise ValueError('Changed geometry rows')
    groups=defaultdict(list)
    with geom_table.open() as stream:
        for pair,geo in zip(pairs,csv.DictReader(stream,delimiter='\t'),strict=True):
            if (pair['center_id'],pair['object_id'])!=(geo['center_id'],geo['object_id']):
                raise ValueError('Pair geometry mismatch')
            obj=pair['object_id']
            if obj in selected and obj in point:
                groups[pair['center_id'],int(geo['radial_bin'])].append((obj,point[obj],float(pair['weight'])))
    data=[r for r in read(INTAKE/'cosmos_join.tsv') if r['eligible_for_selection_review']=='True']
    if len(data)!=reports[INTAKE]['eligible_rows']:
        raise ValueError('COSMOS eligibility mismatch')
    photo,ref,shape,som=np.asarray([[float(r[k]) for k in ('photo_z','reference_z','shape_weight','som_weight')] for r in data]).T
    edges=np.asarray(reports[CALIB]['edges'])
    common=reports[SUPPORT]['per_annulus_common_lens_ids']['0.2']
    saved_weights={(r['lens_id'],int(r['annulus'])):float(r['pair_weight_sum'])
                   for r in read(SUPPORT/'center_annulus_pointz_support.tsv')
                   if r['role']=='lens' and float(r['gap'])==0.2}
    prior={(r['lens_id'],int(r['annulus'])):r for r in read(PREVIOUS/'lens_radial_photoz_response.tsv')}
    cosmo=Cosmology()
    prepared=[]
    for center in centers:
        if center['role']!='lens':
            continue
        lid=center['matched_lens_id']; zl=float(center['z_lens']); chi=float(cosmo.comoving(zl))
        active=photo>zl+0.2
        p,r=photo[active],ref[active]
        est=1-chi/cosmo.comoving(p)
        true=np.zeros_like(r); back=r>zl
        true[back]=1-chi/cosmo.comoving(r[back])
        cal_weight=shape[active]*som[active]*est**2
        cal_bin=np.searchsorted(edges,p,side='right')-1
        mass=np.bincount(cal_bin,weights=cal_weight,minlength=len(edges)-1)
        num=np.bincount(cal_bin,weights=cal_weight*true/est,minlength=len(edges)-1)
        means=np.divide(num,mass,out=np.zeros_like(mass),where=mass>0)
        for annulus,members in enumerate(common):
            if lid not in members:
                continue
            sample=[(obj,z,w) for obj,z,w in groups[center['center_id'],annulus] if z>zl+0.2]
            objs=[x[0] for x in sample]
            z=np.array([x[1] for x in sample])
            w=np.array([x[2] for x in sample])*(1-chi/cosmo.comoving(z))**2
            obin=np.searchsorted(edges,z,side='right')-1
            if len(set(objs))!=len(objs) or not np.all((obin>=0)&(obin<len(mass))) or np.any(mass[obin]<=0):
                raise ValueError('Invalid observed support')
            baseline=float(1/(w@means[obin]/w.sum()))
            if not np.isclose(baseline,float(prior[lid,annulus]['radial_multiplier']),rtol=1e-12):
                raise ValueError('Prior radial multiplier mismatch')
            pair_scale=saved_weights[lid,annulus]/w.sum()
            prepared.append(dict(lid=lid,annulus=annulus,objs=objs,w=w,obin=obin,pair_scale=pair_scale,
                                 active=active,cal_bin=cal_bin,cal_weight=cal_weight,
                                 cal_num=cal_weight*true/est,baseline=baseline))
    if len(prepared)!=24:
        raise ValueError('Unexpected supported sample')
    unique=sorted({obj for case in prepared for obj in case['objs']})
    obj_index={obj:i for i,obj in enumerate(unique)}
    for case in prepared:
        case['obs_index']=np.asarray([obj_index[obj] for obj in case['objs']])
    rng=np.random.default_rng(SEED)
    draws=np.full((REPLICAS,len(prepared)),np.nan)
    pair_weight=np.full((REPLICAS,len(prepared)),np.nan)
    reasons=np.zeros((REPLICAS,len(prepared)),dtype=np.uint8)
    for k in range(REPLICAS):
        reference_mult=rng.poisson(1,len(data))
        observed_mult=rng.poisson(1,len(unique))
        for j,case in enumerate(prepared):
            cm=reference_mult[case['active']]
            mass=np.bincount(case['cal_bin'],weights=case['cal_weight']*cm,minlength=len(edges)-1)
            num=np.bincount(case['cal_bin'],weights=case['cal_num']*cm,minlength=len(edges)-1)
            obs_w=case['w']*observed_mult[case['obs_index']]
            if obs_w.sum()<=0:
                reasons[k,j]=1
                continue
            if np.any((obs_w>0)&(mass[case['obin']]<=0)):
                reasons[k,j]=2
                continue
            means=np.divide(num,mass,out=np.zeros_like(mass),where=mass>0)
            response=obs_w@means[case['obin']]/obs_w.sum()
            if response>0:
                draws[k,j]=1/response
                pair_weight[k,j]=obs_w.sum()*case['pair_scale']
            else:
                reasons[k,j]=3
    valid=np.isfinite(draws).all(axis=1)
    if valid.sum()<REPLICAS*0.8:
        raise ValueError(f'Too few jointly supported draws: {valid.sum()}')
    x=draws[valid]
    cov=np.cov(x,rowvar=False,ddof=1)
    sampled_weights=pair_weight[valid]
    fixed_weights=np.asarray([c['w'].sum()*c['pair_scale'] for c in prepared])
    stack_draws=np.empty((valid.sum(),8))
    fixed_draws=np.empty((valid.sum(),8))
    stack_baseline=[]
    for annulus in range(8):
        positions=[j for j,c in enumerate(prepared) if c['annulus']==annulus]
        base_multipliers=np.asarray([prepared[j]['baseline'] for j in positions])
        fw=fixed_weights[positions]
        sw=sampled_weights[:,positions]
        stack_draws[:,annulus]=sw.sum(axis=1)/np.sum(sw/x[:,positions],axis=1)
        fixed_draws[:,annulus]=fw.sum()/np.sum(fw/x[:,positions],axis=1)
        stack_baseline.append(float(fw.sum()/np.sum(fw/base_multipliers)))
    stack_covariance=np.cov(stack_draws,rowvar=False,ddof=1)
    output.mkdir(parents=True)
    table=output/'joint_bootstrap_lens_annulus.tsv'
    matrix=output/'joint_bootstrap_covariance.tsv'
    stack_table=output/'joint_bootstrap_stack.tsv'
    stack_matrix=output/'joint_bootstrap_stack_covariance.tsv'
    rejection_table=output/'bootstrap_support_sensitivity.tsv'
    write_table(table,[dict(lens_id=c['lid'],annulus=c['annulus'],n_pairs=len(c['objs']),baseline=c['baseline'],
                            sd=float(np.std(x[:,i],ddof=1)),low=float(np.percentile(x[:,i],2.5)),
                            high=float(np.percentile(x[:,i],97.5))) for i,c in enumerate(prepared)])
    write_table(matrix,[dict(lens_a=a['lid'],annulus_a=a['annulus'],lens_b=b['lid'],annulus_b=b['annulus'],covariance=float(cov[i,j]))
                        for i,a in enumerate(prepared) for j,b in enumerate(prepared)])
    write_table(stack_table,[dict(annulus=i,baseline=stack_baseline[i],
                                  varying_weight_sd=float(np.std(stack_draws[:,i],ddof=1)),
                                  fixed_weight_sd=float(np.std(fixed_draws[:,i],ddof=1)),
                                  low=float(np.percentile(stack_draws[:,i],2.5)),
                                  high=float(np.percentile(stack_draws[:,i],97.5)))
                             for i in range(8)])
    write_table(stack_matrix,[dict(annulus_a=i,annulus_b=j,covariance=float(stack_covariance[i,j]))
                              for i in range(8) for j in range(8)])
    rejection=[]
    for annulus in range(8):
        positions=[j for j,c in enumerate(prepared) if c['annulus']==annulus]
        marginal=np.isfinite(draws[:,positions]).all(axis=1)
        mw=pair_weight[marginal][:,positions]
        mx=draws[marginal][:,positions]
        marginal_stack=mw.sum(axis=1)/np.sum(mw/mx,axis=1)
        rejection.append(dict(annulus=annulus,valid_joint=int(valid.sum()),
                              valid_annulus=int(marginal.sum()),
                              annulus_only_sd=float(np.std(marginal_stack,ddof=1)),
                              joint_valid_sd=float(np.std(stack_draws[:,annulus],ddof=1)),
                              annulus_only_low=float(np.percentile(marginal_stack,2.5)),
                              annulus_only_high=float(np.percentile(marginal_stack,97.5)),
                              zero_observed_draws=int(np.any(reasons[:,positions]==1,axis=1).sum()),
                              missing_reference_draws=int(np.any(reasons[:,positions]==2,axis=1).sum()),
                              zero_response_draws=int(np.any(reasons[:,positions]==3,axis=1).sum())))
    write_table(rejection_table,rejection)
    report=dict(science_measurement=False,calibration_applied=False,replicas=REPLICAS,seed=SEED,
                valid_joint=int(valid.sum()),n_cases=len(prepared),unique_observed_sources=len(unique),
                inputs={str(f/'summary.json'):sha(f/'summary.json') for f in reports},
                runner_sha256=sha(Path(__file__)),output_sha256={p.name:sha(p) for p in (table,matrix,stack_table,stack_matrix,rejection_table)},
                limitation='Pair weights are normalized to saved annular totals; within-annulus weight variation '
                           'uses the interpolated distance grid. Independent Poisson object weights, shared by all appearances of each observed source '
                           'and COSMOS object. Jointly invalid draws discarded. No spatial clustering, reference-z '
                           'errors, selection uncertainty or cosmic variance. Sparse cases may be unstable.')
    (output/'summary.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(report,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    run(parser.parse_args().output.resolve())
