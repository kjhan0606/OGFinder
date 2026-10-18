#!/usr/bin/env python3
"""Independent scaled SVD replay and template observability audit."""
import argparse
import csv
import json
from pathlib import Path
import platform
import numpy as np
from inject_neighbor_recovery import profile,design
import prepare_detection_psf as prep


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--work',type=Path,required=True);a=p.parse_args()
    if not platform.node().lower().startswith('lageunha'):p.error('Run on lageunha CPU')
    s=json.loads((a.work/'summary.json').read_text())
    with (a.work/'recovery.tsv').open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
    lookup={(r['cluster_id'],r['scenario'],r['repeat'],r['true_diffuse_flux'],r['true_neighbor_halo_flux'],r['method']):r for r in rows}
    audit=[];observability=[];repeat=s['repeats']-1
    for case in s['cases']:
        cid=case['id'];d=np.load(a.work/'inputs'/(cid+'.npz'));var=d['variance'];base=d['base'];center=d['center'];members=d['members'];shape=base.shape;kpp=case['kpc_per_pixel'];n=len(members)
        core,diffuse,cores,halos,sky=design(shape,center,members,kpp)
        valid=~base;weight=1/np.sqrt(var[valid]);noise=np.random.default_rng(2026101300+int(cid)*10+repeat).normal(size=shape)*np.sqrt(var)
        for scenario in s['scenarios']:
            wing=scenario=='psf_wing_mismatch'
            injected=4000*profile(shape,center,5/kpp,q=.85,theta=.3,wing=wing)+4500*profile(shape,center,35/kpp,wing=wing)
            injected+=sum((2000+1000*i)*profile(shape,xy,2/kpp,q=.75,theta=.5*i,wing=wing) for i,xy in enumerate(members))
            injected+=sum(1000*profile(shape,xy,(12 if scenario=='extended_neighbor_halo' else 8)/kpp,q=.8,theta=.5*i,wing=wing) for i,xy in enumerate(members))
            bg=.03*sky[0]+.02*sky[1]-.01*sky[2]
            if scenario=='curved_sky':
                curve=sky[1]**2+sky[2]**2;bg+=.05*(curve-curve.mean())
            image=np.asarray(injected+bg+noise,dtype='f4')
            for method,cols in [('neighbor_core',[core,diffuse]+cores+sky),('neighbor_core_halo',[core,diffuse]+cores+halos+sky)]:
                matrix=np.column_stack([c[valid]*weight for c in cols]);norms=np.linalg.norm(matrix,axis=0)
                # Scaling makes rank tests invariant to flux-amplitude vs per-pixel sky units.
                scaled=matrix/norms;solution,_,rank,_=np.linalg.lstsq(scaled,image[valid]*weight,rcond=None);coef=solution/norms
                assert rank==len(cols)
                row=lookup[(cid,scenario,str(repeat),'4500.0',str(float(1000*n)),method)]
                target=float(row['diffuse_estimate']);assert np.isclose(coef[1],target,atol=1e-5,rtol=1e-8),(cid,scenario,method)
                if method=='neighbor_core_halo':assert np.isclose(coef[2+n:2+2*n].sum(),float(row['neighbor_halo_estimate']),atol=1e-5,rtol=1e-8)
                audit.append(dict(id=cid,scenario=scenario,method=method,diffuse_absolute_difference=float(abs(coef[1]-target))))
                if scenario=='matched':
                    _,singular,vt=np.linalg.svd(scaled,full_matrices=False)
                    covariance=((vt.T/singular**2)@vt)/norms[:,None]/norms[None,:]
                    visible=[float(x[valid].sum()) for x in cols[:2+n+(n if method=='neighbor_core_halo' else 0)]]
                    record=dict(id=cid,method=method,visible_template_flux_fractions=visible,
                        diffuse_nominal_sigma=float(np.sqrt(covariance[1,1])),
                        any_component_below_one_percent_visible=any(x<.01 for x in visible),
                        unit_scaled_full_rank=True,science_validated=False)
                    if method=='neighbor_core_halo':
                        hcov=covariance[2+n:2+2*n,2+n:2+2*n];record['summed_halo_nominal_sigma']=float(np.sqrt(hcov.sum()))
                    observability.append(record)
    prep.write_json(a.work/'numerical_replay_audit.json',dict(passed=True,n_independent_svd_replays=len(audit),
        max_diffuse_absolute_difference=max(x['diffuse_absolute_difference'] for x in audit),records=audit,
        observability=observability,science_measurement=False,
        note='Column-normalized independent SVD agrees with production normal-equation solve. Raw-unit rank truncation is inappropriate when hidden compact templates have tiny column norms. Full rank does not establish information or astrophysical identifiability.'))
    print(f'Passed {len(audit)} independent SVD replays; max diffuse difference {max(x["diffuse_absolute_difference"] for x in audit):.3g}')
    for r in observability:
        if r['method']=='neighbor_core_halo':print(r['id'],'diffuse sigma',r['diffuse_nominal_sigma'],'halo sigma',r['summed_halo_nominal_sigma'],'min visible fraction',min(r['visible_template_flux_fractions']))


if __name__=='__main__':main()
