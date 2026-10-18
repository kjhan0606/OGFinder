#!/usr/bin/env python3
"""Audit a bounded synthetic experiment and summarize transfer, not science validity."""
import argparse
import collections
import csv
import hashlib
import json
from pathlib import Path
import platform
import numpy as np
import prepare_detection_psf as prep


def transfer(rows):
    """Complete paired amplitude grid only; never estimate from detection survivors."""
    if len(rows)!=6 or any(r['status']!='ok' for r in rows):return None
    truth=np.array([[1,float(r['true_diffuse_flux']),float(r['true_neighbor_halo_flux'])] for r in rows])
    if np.linalg.matrix_rank(truth)!=3:return None
    estimated=np.array([float(r['diffuse_estimate']) for r in rows])
    coef=np.linalg.lstsq(truth,estimated,rcond=None)[0]
    result=dict(diffuse_offset=float(coef[0]),diffuse_gain=float(coef[1]),halo_to_diffuse_leakage=float(coef[2]))
    if all(r['neighbor_halo_estimate']!='' for r in rows):
        hcoef=np.linalg.lstsq(truth,np.array([float(r['neighbor_halo_estimate']) for r in rows]),rcond=None)[0]
        result.update(halo_offset=float(hcoef[0]),diffuse_to_halo_leakage=float(hcoef[1]),halo_gain=float(hcoef[2]))
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--work',type=Path,required=True);a=p.parse_args()
    if not platform.node().lower().startswith('lageunha'):p.error('Run on lageunha CPU')
    summary=json.loads((a.work/'summary.json').read_text());assert summary['status']=='complete'
    with (a.work/'recovery.tsv').open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
    assert len(rows)==summary['n_method_rows']==summary['n_scenes']*5
    keys=[tuple(r[k] for k in ('cluster_id','scenario','repeat','true_diffuse_flux','true_neighbor_halo_flux','method')) for r in rows]
    assert len(set(keys))==len(keys)
    cases={c['id']:c for c in summary['cases']};by_scene=collections.defaultdict(list)
    for cid,c in cases.items():assert hashlib.sha256((a.work/'inputs'/(cid+'.npz')).read_bytes()).hexdigest()==c['saved_input_sha256']
    for r,key in zip(rows,keys):
        by_scene[key[:-1]].append(r)
        assert r['scenario'] in summary['scenarios']
        assert float(r['true_diffuse_flux']) in [1500.,4500.]
        assert float(r['true_neighbor_halo_flux'])/int(r['n_neighbors']) in [0,500,1000]
        if r['status']=='ok':
            for field in ('diffuse_estimate','diffuse_error','nominal_diffuse_sigma','condition_number'):assert np.isfinite(float(r[field]))
            assert abs(float(r['diffuse_estimate'])-float(r['true_diffuse_flux'])-float(r['diffuse_error']))<1e-6
            assert 0<=float(r['coverage_10_80kpc'])<=1
        else:assert r['diffuse_estimate']==''
    for rr in by_scene.values():
        assert len(rr)==5 and len({r['seed'] for r in rr})==1
        masks=[next(r for r in rr if r['method']==name) for name in ['mask_k2.5','mask_k3.5','mask_k5.0']]
        if all(r['status']=='ok' for r in masks):assert all(float(x['coverage_10_80kpc'])>=float(y['coverage_10_80kpc']) for x,y in zip(masks[:-1],masks[1:]))
    aggregate=[];transfers=[]
    for scenario in summary['scenarios']:
        for method in sorted({r['method'] for r in rows}):
            selected=[r for r in rows if r['scenario']==scenario and r['method']==method];good=[r for r in selected if r['status']=='ok']
            paired=collections.defaultdict(list)
            for r in selected:paired[(r['cluster_id'],r['repeat'])].append(r)
            trs=[]
            for (cid,repeat),rr in paired.items():
                tr=transfer(rr)
                if tr is not None:
                    tr=dict(cluster_id=cid,repeat=int(repeat),scenario=scenario,method=method,**tr);trs.append(tr);transfers.append(tr)
            biases=[float(r['diffuse_bias_fraction']) for r in good]
            aggregate.append(dict(scenario=scenario,method=method,n_expected=len(selected),n_estimated=len(good),
                blocked=dict(collections.Counter(r['status'] for r in selected if r['status']!='ok')),
                median_bias_fraction=float(np.median(biases)) if biases else None,
                median_absolute_bias_fraction=float(np.median(np.abs(biases))) if biases else None,
                error_p16_p84=[float(x) for x in np.percentile(biases,[16,84])] if biases else None,
                n_complete_paired_grids=len(trs),median_diffuse_gain=float(np.median([x['diffuse_gain'] for x in trs])) if trs else None,
                median_halo_to_diffuse_leakage=float(np.median([x['halo_to_diffuse_leakage'] for x in trs])) if trs else None,
                median_diffuse_to_halo_leakage=float(np.median([x['diffuse_to_halo_leakage'] for x in trs if 'diffuse_to_halo_leakage' in x])) if method=='neighbor_core_halo' and trs else None))
    correlations=[]
    for scenario in summary['scenarios']:
        groups=collections.defaultdict(list)
        for r in rows:
            if r['scenario']==scenario and r['method']=='neighbor_core_halo' and r['status']=='ok':
                groups[(r['cluster_id'],r['true_diffuse_flux'],r['true_neighbor_halo_flux'])].append([float(r['diffuse_error']),float(r['neighbor_halo_error'])])
        centered=[np.array(v)-np.mean(v,axis=0) for v in groups.values() if len(v)>1]
        rho=float(np.corrcoef(np.concatenate(centered).T)[0,1]) if centered else None
        correlations.append(dict(scenario=scenario,within_cell_error_correlation=rho,
            note='Descriptive only: amplitude cells share seeds; no independent-sample p-value.'))
    prep.write_json(a.work/'artifact_audit.json',dict(passed=True,n_scenes=len(by_scene),n_method_rows=len(rows),
        input_hashes_verified=True,grid_uniqueness=True,paired_noise_seeds=True,mask_growth_checked=True,science_measurement=False))
    prep.write_json(a.work/'transfer_report.json',dict(aggregate=aggregate,transfers=transfers,error_correlations=correlations,
        science_measurement=False,observational_validation_passed=False,
        interpretation='Finite-stamp transfer in a conditional toy model. Matched template family is an optimistic control; no observed ICL/satellite relation measured.'))
    from matplotlib import pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(12,4.5));methods=sorted({r['method'] for r in rows});x=np.arange(len(methods))
    for j,scenario in enumerate(summary['scenarios']):
        rr=[next(r for r in aggregate if r['scenario']==scenario and r['method']==method) for method in methods]
        axes[0].plot(x+j*.025,[np.nan if r['median_bias_fraction'] is None else 100*r['median_bias_fraction'] for r in rr],'o-',label=scenario)
        axes[1].plot(x+j*.025,[np.nan if r['median_halo_to_diffuse_leakage'] is None else r['median_halo_to_diffuse_leakage'] for r in rr],'o-')
    axes[0].set_ylabel('Median diffuse flux error (%)');axes[1].set_ylabel('Recovered diffuse / injected halo slope')
    for ax in axes:ax.axhline(0,color='gray',ls=':');ax.set_xticks(x,methods,rotation=25,ha='right');ax.grid(alpha=.2)
    axes[0].legend(fontsize=7);fig.suptitle('Conditional synthetic recovery; known geometry and assumed PSF');fig.tight_layout();fig.savefig(a.work/'transfer_comparison.png',dpi=140);plt.close(fig)
    print(f'Audited {len(by_scene)} scenes / {len(rows)} method records')
    for r in aggregate:
        print(r['scenario'],r['method'],'n=',r['n_estimated'],
              'median_bias=',r['median_bias_fraction'],'halo_leakage=',r['median_halo_to_diffuse_leakage'])
    print('Error correlations:',json.dumps(correlations))


if __name__=='__main__':main()
