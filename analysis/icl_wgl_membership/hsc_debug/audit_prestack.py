#!/usr/bin/env python3
"""Independent pixel, refit, bootstrap and covariance audit of pre-stack validation."""
import argparse
import csv
import json
from pathlib import Path
import platform
import sys

import numpy as np
import astropy
import sep
from astropy.io import fits

import prepare_detection_psf as prep
from download_pilot import file_sha256 as sha


def read(path):
    return json.loads(path.read_text())


def unpack(z, key, shape):
    return np.unpackbits(z[key], axis=-1, count=int(shape[1]), bitorder='little').astype(bool)


def independent_profile(image, variance, mask, center, kpp):
    y,x=np.ogrid[:image.shape[0],:image.shape[1]]
    radius=np.hypot(x-center[0],y-center[1])*kpp
    valid=~mask & np.isfinite(image) & np.isfinite(variance) & (variance>0)
    raw=[];design=[];coverage=[];edges=[50,100,150,200,300,400,500]
    for lo,hi in zip(edges[:-1],edges[1:]):
        keep=valid & (radius>=lo) & (radius<hi); yy,xx=np.nonzero(keep);n=len(xx)
        raw.append(float(image[keep].mean()) if n else np.nan)
        design.append([1.,(xx.mean()-center[0])*kpp/1000 if n else 0,
                       (yy.mean()-center[1])*kpp/1000 if n else 0])
        coverage.append(min(1.,n/(np.pi*(hi**2-lo**2)/kpp**2)))
    tiles=[]
    for iy in range(0,image.shape[0],64):
        for ix in range(0,image.shape[1],64):
            sl=np.s_[iy:iy+64,ix:ix+64];region=(radius[sl]>=600)&(radius[sl]<900)
            keep=valid[sl]&region;n=int(keep.sum());area=int(region.sum())
            if area<128 or n<128 or n<.5*area:continue
            yy,xx=np.nonzero(keep)
            tiles.append([1.,(ix+xx.mean()-center[0])*kpp/1000,(iy+yy.mean()-center[1])*kpp/1000,
                          image[sl][keep].astype(float).mean(),variance[sl][keep].astype(float).sum()/n**2])
    corrected=np.full(6,np.nan)
    if len(tiles)>=10:
        tile=np.array(tiles);X=tile[:,:3];w=1/tile[:,4]
        normal=X.T@(w[:,None]*X)
        if np.linalg.matrix_rank(normal)==3:
            coef=np.linalg.solve(normal,X.T@(w*tile[:,3]))
            corrected=np.array(raw)-np.array(design)@coef
    return dict(raw=np.array(raw),corrected=corrected,coverage=np.array(coverage))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--extraction',type=Path,required=True);p.add_argument('--analysis',type=Path,required=True)
    p.add_argument('--work',type=Path,required=True);a=p.parse_args()
    if not platform.node().lower().startswith('lageunha'):p.error('Run on lageunha CPU')
    a.work.mkdir(parents=True,exist_ok=False)
    doc=read(a.extraction/'extraction_summary.json');summary=read(a.analysis/'summary.json')
    assert doc['status']==summary['status']=='complete'
    assert not summary['science_measurement'] and not summary['sky_correction_approved']
    code=dict(doc['script_sha256'],**summary['script_sha256'])
    if summary.get('paired_cluster_input'):
        code.update(read(Path(summary['paired_cluster_input'])/'summary.json')['script_sha256'])
    code[str(Path(__file__))]=sha(Path(__file__))
    snapshots=a.work/'source_snapshots';snapshots.mkdir()
    for path,digest in code.items():
        assert sha(path)==digest,path
        with (snapshots/(digest+'_'+Path(path).name)).open('xb') as stream:stream.write(Path(path).read_bytes())
    for path,digest in dict(code,**doc['source_sha256'],**summary['source_sha256']).items():
        assert sha(path)==digest,path
    records=doc['records'];hash_count=0
    for record in records:
        assert sha(record['source_image'])==record['source_sha256']
        baseline=a.extraction/record['control_id']/'baseline.npz'
        assert sha(baseline)==record['baseline_sha256']
        for path,digest in record['injection_sha256'].items():assert sha(path)==digest
        hash_count+=2+len(record['injection_sha256'])
    old_mask_checks=0
    for record in records:
        if not record['n_injection_models']:continue
        baseline=a.extraction/record['control_id']/'baseline.npz'
        previous=Path(record['source_image']).parents[2]/'control_outer_2026-10-14_01'/record['control_id']/'masks.fits.gz'
        with np.load(baseline) as z,fits.open(previous,memmap=False) as h:
            masks=unpack(z,'mask_bits',z['image_shape'])
            for j,factor in enumerate(('2P5','3P5','5P0')):
                np.testing.assert_array_equal(masks[j],h['MASK_K'+factor].data.astype(bool));old_mask_checks+=1
    chosen=[r for r in records if r['target_cluster_id'] in ('5724','619','5814','5104','3505','6940')]
    if not chosen:chosen=[r for r in records if r['n_injection_models']][:1]
    checks=0;max_difference=0.;common_checks=0;fixed_checks=0
    for record in chosen:
        folder=a.extraction/record['control_id']
        with np.load(folder/'baseline.npz') as z:
            shape=z['image_shape'];center=z['center'];kpp=float(z['kpp']);factor=float(z['native_to_sb'])
            masks=unpack(z,'mask_bits',shape);baseline={k:z[k].copy() for k in ('raw','corrected','coverage')}
        with fits.open(record['source_image'],memmap=False) as h:
            image=h[1].data.astype(float);variance=h[3].data.astype(float)
        for j,mask in enumerate(masks):
            result=independent_profile(image,variance,mask,center,kpp)
            for mode in ('raw','corrected'):
                delta=np.nanmax(abs(result[mode]*factor-baseline[mode][j]));max_difference=max(max_difference,float(delta))
                np.testing.assert_allclose(result[mode]*factor,baseline[mode][j],rtol=2e-9,atol=1e-9,equal_nan=True)
                checks+=6
        for mi in range(record['n_injection_models']):
            with np.load(folder/f'injection_{mi:02d}.npz') as z:
                hscale,q,amplitude=z['parameters'];after_masks=unpack(z,'after_mask_bits',shape)
                y,x=np.indices(image.shape,dtype=float);x=(x-center[0])*kpp;y=(y-center[1])*kpp
                theta=37*np.pi/180;rho=np.sqrt((x*np.cos(theta)+y*np.sin(theta))**2+((-x*np.sin(theta)+y*np.cos(theta))/q)**2)
                model=amplitude*np.exp(-(rho-100)/hscale)/factor
                for j,mask in enumerate(after_masks):
                    result=independent_profile(image+model,variance,mask,center,kpp)
                    fixed=independent_profile(model,variance,masks[j],center,kpp)
                    union=mask|masks[j]
                    common_before=independent_profile(image,variance,union,center,kpp)
                    common_after=independent_profile(image+model,variance,union,center,kpp)
                    for mode in ('raw','corrected'):
                        expected=z[f'after_{mode}_{j}'];actual=result[mode]*factor
                        np.testing.assert_allclose(actual,expected,rtol=2e-9,atol=1e-9,equal_nan=True)
                        np.testing.assert_allclose(fixed[mode]*factor,z[f'fixed_signal_{mode}_{j}'],rtol=2e-9,atol=1e-9,equal_nan=True)
                        difference=(common_after[mode]-common_before[mode])*factor
                        np.testing.assert_allclose(difference,z[f'common_signal_{mode}_{j}'],rtol=2e-8,atol=1e-9,equal_nan=True)
                        max_difference=max(max_difference,float(np.nanmax(abs(actual-expected))))
                        checks+=6;common_checks+=6;fixed_checks+=6
        print('Pixel replay',record['target_cluster_id'],'complete',flush=True)
    profiles=[]
    for record in records:
        with np.load(a.extraction/record['control_id']/'baseline.npz') as z:
            profiles.append({key:z[key].copy() for key in ('raw','corrected','coverage','nominal_raw','nominal_corrected')})
    report=read(a.analysis/'covariance_report.json');matrices=np.load(a.analysis/'covariances.npz')
    cluster=None
    if summary.get('paired_cluster_input'):
        cp=Path(summary['paired_cluster_input']);cd=read(cp/'summary.json')
        for path,digest in cd['source_sha256'].items():assert sha(path)==digest
        cluster=np.load(cp/'profiles.npz')
    cov_checks=0
    for item in report:
        if item['status']!='pilot_estimate':continue
        name=item['name'];selected=matrices[name+'__selected_indices'];fi=(2.5,3.5,5.).index(item['kron_factor']);bins=item['bin_indices']
        if item.get('sample_kind')=='paired_cluster_minus_control':
            linked=matrices['paired_control_indices']
            values=np.array([cluster[item['extraction']][i,fi][bins]-profiles[linked[i]][item['extraction']][fi][bins] for i in selected])
        else:
            values=np.array([profiles[i][item['extraction']][fi][bins] for i in selected])
        n=len(values)
        centered=values-values.mean(axis=0);expected=centered.T@centered/(n-1)
        np.testing.assert_allclose(expected,matrices[name+'__sample_covariance'],rtol=1e-11,atol=1e-11)
        np.testing.assert_allclose(expected/n,matrices[name+'__independent_stack_covariance'],rtol=1e-11,atol=1e-11)
        for key in ('sample_covariance','iid_bootstrap_covariance','spatial_bootstrap_covariance','delete_cell_covariance'):
            matrix=matrices[name+'__'+key]
            if np.isfinite(matrix).all():
                np.testing.assert_allclose(matrix,matrix.T,atol=1e-12)
                assert np.linalg.eigvalsh(matrix).min()>=-1e-9*max(1.,np.linalg.norm(matrix))
        cov_checks+=1
    result=dict(passed=True,n_source_and_output_hashes_verified=hash_count,n_pixel_profile_values_checked=checks,
        n_common_support_increment_checks=common_checks,n_fixed_support_signal_checks=fixed_checks,
        n_covariance_configs_replayed=cov_checks,n_previous_masks_reproduced=old_mask_checks,
        maximum_profile_absolute_difference=max_difference,script_sha256=sha(Path(__file__)),
        audit_environment=dict(python=sys.version,numpy=np.__version__,astropy=astropy.__version__,sep=sep.__version__),
        extraction=str(a.extraction),analysis=str(a.analysis),science_measurement=False,sky_correction_approved=False)
    prep.write_json(a.work/'audit.json',result);print(json.dumps(result))


if __name__=='__main__':main()
