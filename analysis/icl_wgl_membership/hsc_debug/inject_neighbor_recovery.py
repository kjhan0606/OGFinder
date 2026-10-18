#!/usr/bin/env python3
"""Controlled neighbour/BCG+diffuse transfer test with real HSC noise/mask geometry.

Synthetic pixels only. Fixed catalogue positions, fixed model scales, linear signed
amplitudes, and an explicitly assumed PSF: this is not an observational measurement
or an end-to-end validation of the HSC reduction or survey background transfer.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import platform
import time

import numpy as np
from scipy.ndimage import gaussian_filter
from astropy.io import fits
from astropy.wcs import WCS
import prepare_detection_psf as prep
from review_bcg_association import peak_candidate

CASE_IDS=('5724','619','5814','5104','3505','6940')
SCENARIOS=('matched','extended_neighbor_halo','psf_wing_mismatch','curved_sky')
SIZE=321


def profile(shape, center, length, q=1., theta=0., wing=False):
    """Unit flux in the finite stamp, exponential radius, assumed Gaussian PSF."""
    yy,xx=np.indices(shape,dtype=float);dx=xx-center[0];dy=yy-center[1]
    x=dx*np.cos(theta)+dy*np.sin(theta);y=-dx*np.sin(theta)+dy*np.cos(theta)
    source=np.exp(-np.hypot(x,y/q)/length)
    blurred=gaussian_filter(source,2.,mode='constant')
    if wing:blurred=.95*blurred+.05*gaussian_filter(source,10.,mode='constant')
    return blurred/blurred.sum()


def design(shape,center,members,kpc_per_pixel):
    core=profile(shape,center,5/kpc_per_pixel,q=.85,theta=.3)
    diffuse=profile(shape,center,35/kpc_per_pixel)
    neighbors=[profile(shape,xy,2/kpc_per_pixel,q=.75,theta=.5*i) for i,xy in enumerate(members)]
    halos=[profile(shape,xy,8/kpc_per_pixel,q=.8,theta=.5*i) for i,xy in enumerate(members)]
    y,x=np.indices(shape,dtype=float);x=(x-center[0])/(SIZE/2);y=(y-center[1])/(SIZE/2)
    sky=[np.ones(shape),x,y]
    return core,diffuse,neighbors,halos,sky


def fit_operator(columns,variance,mask):
    """Weighted linear fit; normalize columns and reject rank/conditioning failures."""
    valid=(~mask)&np.isfinite(variance)&(variance>0)
    if valid.sum()<10*len(columns):return None,'too_few_pixels'
    weight=1/np.sqrt(variance[valid]);matrix=np.column_stack([c[valid]*weight for c in columns])
    norms=np.linalg.norm(matrix,axis=0)
    if np.any(norms<=0):return None,'zero_column'
    scaled=matrix/norms
    singular=np.linalg.svd(scaled,compute_uv=False)
    condition=float(singular[0]/singular[-1]) if singular[-1]>0 else float('inf')
    if condition>1e6:return None,'ill_conditioned'
    gram=scaled.T@scaled
    covariance=np.linalg.inv(gram)/norms[:,None]/norms[None,:]
    # The projector is conditional on fixed templates, positions and diagonal noise.
    projection=np.linalg.solve(gram,scaled.T)/norms[:,None]*weight[None,:]
    return dict(valid=valid,projection=projection,covariance=covariance,condition=condition),'ok'


def estimate(operator,image):
    return operator['projection']@image[operator['valid']]


def load_case(frame):
    info=prep.validate_fits(Path(frame['image']),frame)
    with fits.open(frame['image'],memmap=False) as h:
        wcs=WCS(h[info['science_hdu']].header).celestial
        cx,cy=info['center_pixel_xy'];ix,iy=round(cx),round(cy);half=SIZE//2
        xs=slice(ix-half,ix+half+1);ys=slice(iy-half,iy+half+1)
        variance=np.array(h[info['variance_hdu']].data[ys,xs],dtype=float)
        bits=np.array(h[info['mask_hdu']].data[ys,xs],dtype=np.int32)
    if variance.shape!=(SIZE,SIZE):raise ValueError('Crop does not fit source image')
    bad=((bits&sum(1<<b for b in prep.BAD_BIT_POSITIONS))!=0)|~np.isfinite(variance)|(variance<=0)
    base=bad|((bits&prep.BRIGHT_OBJECT_BIT)!=0)
    sigma=float(np.sqrt(np.median(variance[~base])))
    variance=np.where(bad,1.,variance/sigma**2)
    center=(cx-(ix-half),cy-(iy-half))
    members=wcs.all_world2pix(np.asarray(frame['members_radec_deg']).reshape(-1,2),0)-[ix-half,iy-half]
    return dict(variance=variance,bad=bad,base=base,center=center,members=members,
                sigma_original=sigma,pixel_scale=info['pixel_scale_arcsec'],source_image=frame['image'],
                crop_origin_xy=[ix-half,iy-half])


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--manifest',type=Path,required=True);p.add_argument('--work',type=Path,required=True)
    p.add_argument('--repeats',type=int,default=3);p.add_argument('--case-limit',type=int,default=6)
    args=p.parse_args()
    if not platform.node().lower().startswith('lageunha'):p.error('Run on lageunha CPU')
    if not 1<=args.repeats<=10 or not 1<=args.case_limit<=6:p.error('Bounded engineering matrix only')
    doc=json.loads(args.manifest.read_text());lookup={str(x['id']):x for x in doc['clusters']}
    args.work.mkdir(parents=True,exist_ok=False);(args.work/'inputs').mkdir();started=time.time()
    cosmo=prep.Cosmology(H0=doc['h0'],Om=doc['omega_m']);rows=[];cases=[]
    for cid in CASE_IDS[:args.case_limit]:
        f=lookup[cid];data=load_case(f);shape=data['base'].shape;center=data['center']
        kpp=float(cosmo.angular_diameter(f['z']))*prep.ARCSEC*1000*data['pixel_scale']
        distances=np.hypot(data['members'][:,0]-center[0],data['members'][:,1]-center[1])*kpp
        eligible=np.flatnonzero(np.isfinite(distances)&(distances>3)&(distances<60))
        ids=sorted(eligible,key=lambda i:distances[i])[:3];members=data['members'][ids]
        if len(members)==0:raise ValueError('No eligible catalogue neighbours')
        core,diffuse,cores,halos,sky=design(shape,center,members,kpp)
        count=len(members);noisevar=data['variance'];base=data['base']
        # Same fixed recovery family across all scenarios. Truth changes independently.
        columns=dict(neighbor_core=[core,diffuse]+cores+sky,
                     neighbor_core_halo=[core,diffuse]+cores+halos+sky)
        operators={key:fit_operator(value,noisevar,base) for key,value in columns.items()}
        yy,xx=np.indices(shape);radius=np.hypot(xx-center[0],yy-center[1])*kpp;ring=(radius>=10)&(radius<80)
        np.savez_compressed(args.work/'inputs'/(cid+'.npz'),variance=noisevar,bad=data['bad'],base=base,center=center,members=members)
        case=dict(id=cid,cohort=f['preparation_cohort'],kpc_per_pixel=kpp,source_image=f['image'],crop_origin_xy=data['crop_origin_xy'],
            normalization_sigma_original=data['sigma_original'],n_neighbors=count,neighbor_radius_kpc=[float(distances[i]) for i in ids],
            central_10_80kpc_base_coverage=float(np.mean(~base[ring])),pixel_scale_arcsec=data['pixel_scale'],
            saved_input_sha256=hashlib.sha256((args.work/'inputs'/(cid+'.npz')).read_bytes()).hexdigest())
        cases.append(case)
        for scenario in SCENARIOS:
            wing=scenario=='psf_wing_mismatch'
            tcore=profile(shape,center,5/kpp,q=.85,theta=.3,wing=wing)
            tdiffuse=profile(shape,center,35/kpp,wing=wing)
            tcores=[profile(shape,xy,2/kpp,q=.75,theta=.5*i,wing=wing) for i,xy in enumerate(members)]
            thalos=[profile(shape,xy,(12 if scenario=='extended_neighbor_halo' else 8)/kpp,q=.8,theta=.5*i,wing=wing) for i,xy in enumerate(members)]
            true_sky=.03*sky[0]+.02*sky[1]-.01*sky[2]
            if scenario=='curved_sky':
                curvature=sky[1]**2+sky[2]**2;true_sky=true_sky+.05*(curvature-curvature.mean())
            for repeat in range(args.repeats):
                seed=2026101300+int(cid)*10+repeat
                noise=np.random.default_rng(seed).normal(size=shape)*np.sqrt(noisevar)
                for diffuse_flux in (1500.,4500.):
                    for halo_factor in (0.,.5,1.):
                        halo_flux=1000*halo_factor
                        injected=4000*tcore+diffuse_flux*tdiffuse+sum((2000+1000*i)*v for i,v in enumerate(tcores))+sum(halo_flux*v for v in thalos)
                        image=np.asarray(injected+true_sky+noise,dtype='f4')
                        detected=prep.detect(image,np.asarray(noisevar,dtype='f4'),data['bad'],center,members,data['pixel_scale'])
                        obj,seg,kron,kflags,_,_,_=detected;bcg,bcgreason=peak_candidate(obj,seg,data['bad'],center,data['pixel_scale'])
                        masks,_=prep.source_masks(obj,seg,kron,kflags,base,bcg,members,data['pixel_scale'])
                        methods=list(operators.items())
                        for factor,mask in zip(prep.FACTORS,masks):
                            op,status=fit_operator([core,diffuse]+sky,noisevar,mask) if bcg>=0 else (None,'bcg_detection_blocked')
                            methods.append(('mask_k'+str(factor),(op,status)))
                        for method,(operator,status) in methods:
                            row=dict(cluster_id=cid,scenario=scenario,repeat=repeat,seed=seed,method=method,
                                true_diffuse_flux=diffuse_flux,true_neighbor_halo_flux=halo_flux*count,
                                status=status,bcg_detection_status=bcgreason,n_detected=len(obj),n_neighbors=count,
                                diffuse_estimate=None,diffuse_error=None,diffuse_bias_fraction=None,
                                neighbor_halo_estimate=None,neighbor_halo_error=None,nominal_diffuse_sigma=None,
                                condition_number=None,retained_fraction=None,coverage_10_80kpc=None)
                            if operator is not None:
                                coef=estimate(operator,image);err=float(coef[1]-diffuse_flux)
                                row.update(diffuse_estimate=float(coef[1]),diffuse_error=err,diffuse_bias_fraction=err/diffuse_flux,
                                    nominal_diffuse_sigma=float(np.sqrt(operator['covariance'][1,1])),condition_number=operator['condition'],
                                    retained_fraction=float(operator['valid'].mean()),coverage_10_80kpc=float(np.mean(operator['valid'][ring])))
                                if method=='neighbor_core_halo':
                                    total=float(coef[2+count:2+2*count].sum());row.update(neighbor_halo_estimate=total,neighbor_halo_error=total-halo_flux*count)
                            rows.append(row)
                        if cid in ('5724','619') and repeat==0 and diffuse_flux==4500 and halo_factor==1:
                            from matplotlib import pyplot as plt
                            fig,axes=plt.subplots(1,3,figsize=(11,3.5))
                            axes[0].imshow(np.arcsinh(image),origin='lower',cmap='gray',vmin=-1,vmax=4);axes[0].set_title('Synthetic scene / noise sigma')
                            axes[1].imshow(masks[0],origin='lower',vmin=0,vmax=1,cmap='gray');axes[1].set_title('Detected Kron 2.5 mask')
                            op,status=operators['neighbor_core_halo']
                            residual=image if op is None else image-sum(c*t for c,t in zip(estimate(op,image),columns['neighbor_core_halo']))
                            axes[2].imshow(np.where(base,np.nan,residual),origin='lower',cmap='coolwarm',vmin=-2,vmax=2);axes[2].set_title('Joint model residual')
                            for ax in axes:ax.plot(center[0],center[1],'r+')
                            fig.suptitle(cid+' / '+scenario+' / assumed PSF; no science approval');fig.tight_layout();fig.savefig(args.work/(cid+'_'+scenario+'.png'),dpi=110);plt.close(fig)
            print(f'{cid}: {scenario} done',flush=True)
    prep.table(args.work/'recovery.tsv',rows,list(rows[0]))
    prep.write_json(args.work/'summary.json',dict(status='complete',n_scenes=len(rows)//5,n_method_rows=len(rows),cases=cases,
        elapsed_seconds=time.time()-started,science_measurement=False,observational_validation_passed=False,
        scenarios=SCENARIOS,repeats=args.repeats,manifest_sha256=hashlib.sha256(args.manifest.read_bytes()).hexdigest(),
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),dependency_sha256=hashlib.sha256(Path(prep.__file__).read_bytes()).hexdigest(),
        association_sha256=hashlib.sha256((Path(__file__).parent/'review_bcg_association.py').read_bytes()).hexdigest(),
        software=dict(numpy=np.__version__,sep=prep.sep.__version__),
        settings=dict(stamp_size=SIZE,neighbor_count_max=3,psf_sigma_pix=2,wing_fraction=.05,wing_sigma_pix=10,
            model_diffuse_scale_kpc=35,model_neighbor_halo_scale_kpc=8,wrong_true_halo_scale_kpc=12,
            diffuse_flux_grid=[1500,4500],per_neighbor_halo_flux_grid=[0,500,1000],curved_sky_amplitude_sigma=.05),
        limitations=['All science image pixels are synthetic; real mask, variance, redshift and member positions only.',
        'Fluxes are finite-stamp integrated amplitudes in units of median input pixel-noise sigma, not AB fluxes.',
        'Known source positions and fixed profile scales are conditional assumptions, not blind fitting.',
        'Independent Gaussian noise; coadd correlation, unresolved population and true PSF unvalidated.',
        'Input variance is held fixed; additional shot noise from the injected sources is not simulated.',
        'Only up to three known neighbours, 321-pixel crop, no full 0-500 kpc validation.',
        'Matched scenario is an optimistic control; mismatch scenarios test transfer failure.',
        'Repeated amplitude/scenario cells share noise seeds; they are not independent replicates.',
        'Linear signed coefficients and nominal conditional covariance do not include model/systematic errors.']))


if __name__=='__main__':main()
