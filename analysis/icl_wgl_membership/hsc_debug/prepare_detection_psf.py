#!/usr/bin/env python3
"""Offline engineering inputs: detections, BCG-aware masks, provisional PSF stars.

Uses the existing research SEP dependency. No science profiles, PSF fallback,
extinction acceptance, random-field substitution or edits to original images.
"""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import platform
import sys
import time

import numpy as np
import sep
from astropy.io import fits
from astropy.wcs import WCS
from scipy.spatial import cKDTree

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[2]))
from ogfkit.wlcluster import ARCSEC, Cosmology
from download_pilot import validate_fits, BAD_BIT_POSITIONS, BRIGHT_OBJECT_BIT

FACTORS = (2.5, 3.5, 5.0)


def plain(value):
    if isinstance(value, dict): return {k:plain(v) for k,v in value.items()}
    if isinstance(value, (list,tuple,np.ndarray)): return [plain(v) for v in value]
    if isinstance(value, (np.bool_,bool)): return bool(value)
    if isinstance(value, (np.integer,int)): return int(value)
    if isinstance(value, (np.floating,float)): return float(value) if np.isfinite(value) else None
    return value


def write_json(path, value):
    with path.open('x') as f: json.dump(plain(value), f, indent=2, allow_nan=False);f.write('\n')


def table(path, rows, fields):
    with path.open('x', newline='') as f:
        w=csv.DictWriter(f, fieldnames=fields, delimiter='\t', extrasaction='ignore')
        w.writeheader();w.writerows(plain(rows))


def circle(mask, x, y, radius):
    ny,nx=mask.shape
    x0,x1=max(0,int(math.floor(x-radius))),min(nx,int(math.ceil(x+radius))+1)
    y0,y1=max(0,int(math.floor(y-radius))),min(ny,int(math.ceil(y+radius))+1)
    if x0>=x1 or y0>=y1: return
    yy,xx=np.ogrid[y0:y1,x0:x1]
    mask[y0:y1,x0:x1] |= (xx-x)**2+(yy-y)**2 <= radius**2


def source_masks(objects, segmentation, kron, kron_flags, base, bcg, members, pixel_scale):
    """Preserve detector defects even on BCG; neighbouring ellipses may overlap it."""
    pos=np.column_stack([objects['x'],objects['y']])
    match=max(1/pixel_scale,1.5)
    nearest=cKDTree(pos).query(members)[0] if len(pos) and len(members) else np.full(len(members),np.inf)
    fallback=[]
    for xy,distance in zip(members,nearest):
        if distance<=match: continue
        fallback.append(xy)
    masks=[]
    for factor in FACTORS:
        mask=base.copy() | ((segmentation>0)&(segmentation!=bcg+1) if bcg>=0 else segmentation>0)
        keep=np.arange(len(objects))!=bcg
        if keep.any():
            scale=factor*np.where((kron_flags==0)&np.isfinite(kron)&(kron>0),kron,1.2)
            sep.mask_ellipse(mask,objects['x'][keep],objects['y'][keep],
                             np.maximum(objects['a'][keep],.5),np.maximum(objects['b'][keep],.5),
                             objects['theta'][keep],r=scale[keep])
        for x,y in fallback: circle(mask,float(x),float(y),max(2/pixel_scale,1.5))
        masks.append(mask)
    return masks,len(fallback)


def detect(image, variance, bad, center, members, pixel_scale, deblend_cont=.005):
    valid=(~bad)&np.isfinite(image)&np.isfinite(variance)&(variance>0)
    if np.count_nonzero(valid)<100: raise ValueError('Too few valid pixels')
    data=np.ascontiguousarray(np.where(valid,image,np.median(image[valid])),dtype=np.float32)
    invalid=np.ascontiguousarray(~valid)
    background=sep.Background(data,mask=invalid,bw=64,bh=64)
    residual=np.ascontiguousarray(data-background.back())
    error=np.ascontiguousarray(np.sqrt(np.where(valid,variance,np.median(variance[valid]))),dtype=np.float32)
    objects,seg=sep.extract(residual,1.5,err=error,mask=invalid,minarea=5,
                            deblend_nthresh=32,deblend_cont=deblend_cont,segmentation_map=True)
    if len(objects):
        kron,flags=sep.kron_radius(residual,objects['x'],objects['y'],objects['a'],objects['b'],objects['theta'],6,mask=invalid)
    else: kron,flags=np.array([]),np.array([],dtype=int)
    distance=np.hypot(objects['x']-center[0],objects['y']-center[1])
    candidates=np.flatnonzero(distance<=max(1/pixel_scale,1.5))
    # Ambiguous and absent BCG detections remain an explicit inspection gate.
    bcg=int(candidates[0]) if len(candidates)==1 else -1
    noise=np.sqrt(np.bincount(seg.ravel(),weights=np.where(valid,variance,0).ravel(),minlength=len(objects)+1)[1:])
    snr=np.divide(objects['flux'],noise,out=np.zeros(len(objects)),where=noise>0)
    return objects,seg,kron,flags,bcg,snr,dict(bcg_matches=len(candidates),
        nearest_bcg_detection_arcsec=float(distance.min()*pixel_scale) if len(distance) else None,
        bcg_position_segment=int(seg[round(center[1]),round(center[0])]),
        background_global=float(background.globalback),background_rms=float(background.globalrms),
        detection_threshold_sigma=1.5,variance_snr_model='nominal independent-pixel; coadd correlation not calibrated')


def psf_candidates(objects,seg,image,variance,base,snr,bcg,members,pixel_scale):
    n=len(objects);positions=np.column_stack([objects['x'],objects['y']])
    fwhm=2.35482*np.sqrt((objects['a']**2+objects['b']**2)/2)
    nearest=cKDTree(positions).query(positions,k=2)[0][:,1] if n>1 else np.full(n,np.inf)
    memberdist=cKDTree(members).query(positions)[0] if len(members) and n else np.full(n,np.inf)
    usable=(snr>=50)&(objects['flag']==0)&(objects['b']/np.maximum(objects['a'],1e-6)>=.8)
    usable &= (fwhm*pixel_scale>=.35)&(fwhm*pixel_scale<=1.5)&(memberdist>1/pixel_scale)
    if bcg>=0: usable[bcg]=False
    usable &= nearest>=np.maximum(12,4*fwhm)
    if np.any(usable):
        edges=np.arange(0,12.25,.25);counts,_=np.histogram(fwhm[usable],edges)
        mode=float((edges[np.argmax(counts)]+edges[np.argmax(counts)+1])/2)
        usable &= abs(fwhm-mode)<=max(.5,.2*mode)
    else: mode=None
    accepted=[];stamps=[];profiles=[];ny,nx=image.shape
    yy,xx=np.indices((101,101));rr=np.hypot(xx-50,yy-50)
    for i in sorted(np.flatnonzero(usable),key=lambda i:-snr[i]):
        x,y=float(objects['x'][i]),float(objects['y'][i]);ix,iy=round(x),round(y)
        if ix<50 or iy<50 or ix+50>=nx or iy+50>=ny: continue
        cut=image[iy-50:iy+51,ix-50:ix+51].copy()
        labels=seg[iy-50:iy+51,ix-50:ix+51]
        masked=base[iy-50:iy+51,ix-50:ix+51] | ((labels>0)&(labels!=i+1))
        core=rr<=12
        if np.mean(masked[core])>.01: continue
        cut[masked]=np.nan
        accepted.append(dict(NUMBER=int(i+1),X_IMAGE=x+1,Y_IMAGE=y+1,SNR_NOMINAL=float(snr[i]),
                             FWHM_MOMENT_PIX=float(fwhm[i]),SUBPIXEL_DX=x-ix,SUBPIXEL_DY=y-iy,
                             STAMP_INDEX=len(stamps),WING_VALID_FRACTION=float(np.isfinite(cut).mean()),
                             CLASSIFICATION='provisional_compact_candidate'))
        stamps.append(cut)
        for lo,hi in zip((0,2,4,6,8,12,16,24,32,40),(2,4,6,8,12,16,24,32,40,50)):
            ring=(rr>=lo)&(rr<hi);valid=ring&np.isfinite(cut)
            profiles.append(dict(NUMBER=int(i+1),R_MIN_PIX=lo,R_MAX_PIX=hi,N_VALID=int(valid.sum()),
                                 AREA_FRACTION=float(valid.sum()/ring.sum()),
                                 RAW_MEAN=float(np.mean(cut[valid])) if valid.any() else None))
        if len(accepted)>=25: break
    return accepted,stamps,profiles,mode


def prepare_frame(frame, cosmology, out):
    info=validate_fits(Path(frame['image']),frame)
    with fits.open(frame['image'],memmap=False) as h:
        image=np.array(h[info['science_hdu']].data,dtype=np.float32)
        variance=np.array(h[info['variance_hdu']].data,dtype=np.float32)
        bits=np.array(h[info['mask_hdu']].data,dtype=np.int32)
        header=h[info['science_hdu']].header.copy()
        header.extend(h[0].header,unique=True)
    wcs=WCS(header).celestial;scale=info['pixel_scale_arcsec'];center=np.array(info['center_pixel_xy'])
    sky=np.asarray(frame['members_radec_deg'],float).reshape(-1,2)
    members=np.column_stack(wcs.all_world2pix(sky[:,0],sky[:,1],0)) if len(sky) else np.empty((0,2))
    inside=np.isfinite(members).all(axis=1)&(members[:,0]>=0)&(members[:,0]<image.shape[1])&(members[:,1]>=0)&(members[:,1]<image.shape[0])
    members_in=members[inside]
    bad=(bits&sum(1<<bit for bit in BAD_BIT_POSITIONS))!=0
    bad |= ~np.isfinite(image)|~np.isfinite(variance)|(variance<=0)
    bright=(bits&BRIGHT_OBJECT_BIT)!=0;base=bad|bright
    objects,seg,kron,kflags,bcg,snr,extra=detect(image,variance,bad,center,members_in,scale)
    masks,fallback=source_masks(objects,seg,kron,kflags,base,bcg,members_in,scale)
    assert all(np.all(mask[base]) for mask in masks)
    assert np.all(~masks[0]|masks[1]) and np.all(~masks[1]|masks[2])
    stars,stamps,profiles,mode=psf_candidates(objects,seg,image,variance,base,snr,bcg,members_in,scale)
    out.mkdir()
    pix=np.column_stack([objects['x'],objects['y']])
    world=np.column_stack(wcs.all_pix2world(pix[:,0],pix[:,1],0)) if len(pix) else np.empty((0,2))
    mdist=cKDTree(members_in).query(pix)[0] if len(members_in) and len(pix) else np.full(len(pix),np.inf)
    rows=[dict(NUMBER=i+1,X_IMAGE=float(o['x'])+1,Y_IMAGE=float(o['y'])+1,RA=world[i,0],DEC=world[i,1],
        A_IMAGE=o['a'],B_IMAGE=o['b'],THETA_RAD=o['theta'],FLUX_ISO=o['flux'],SNR_NOMINAL=snr[i],
        FLAGS=int(o['flag']),KRON_RADIUS=kron[i],KRON_FLAGS=int(kflags[i]),BCG_CANDIDATE=i==bcg,
        CAMIRA_POSITION_MATCH=bool(mdist[i]<=1/scale)) for i,o in enumerate(objects)]
    table(out/'detections.tsv',rows,['NUMBER','X_IMAGE','Y_IMAGE','RA','DEC','A_IMAGE','B_IMAGE','THETA_RAD','FLUX_ISO','SNR_NOMINAL','FLAGS','KRON_RADIUS','KRON_FLAGS','BCG_CANDIDATE','CAMIRA_POSITION_MATCH'])
    positions=cKDTree(pix).query(members)[0] if len(pix) and len(members) else np.full(len(members),np.inf)
    table(out/'members.tsv',[dict(INDEX=i,RA=xy[0],DEC=xy[1],X_IMAGE=members[i,0]+1,Y_IMAGE=members[i,1]+1,
        IN_FRAME=inside[i],DETECTION_SEPARATION_ARCSEC=positions[i]*scale) for i,xy in enumerate(sky)],
        ['INDEX','RA','DEC','X_IMAGE','Y_IMAGE','IN_FRAME','DETECTION_SEPARATION_ARCSEC'])
    table(out/'psf_candidates.tsv',stars,['NUMBER','X_IMAGE','Y_IMAGE','SNR_NOMINAL','FWHM_MOMENT_PIX','SUBPIXEL_DX','SUBPIXEL_DY','STAMP_INDEX','WING_VALID_FRACTION','CLASSIFICATION'])
    table(out/'psf_raw_profiles.tsv',profiles,['NUMBER','R_MIN_PIX','R_MAX_PIX','N_VALID','AREA_FRACTION','RAW_MEAN'])
    wh=wcs.to_header();wh['SCIVAL']=False;wh['BCGOK']=bool(bcg>=0)
    hdus=[fits.PrimaryHDU(),fits.ImageHDU(seg.astype(np.int32),header=wh,name='SEGMENTATION'),
          fits.ImageHDU(base.astype(np.uint8),header=wh,name='BASE_MASK')]
    for factor,mask in zip(FACTORS,masks): hdus.append(fits.ImageHDU(mask.astype(np.uint8),header=wh,name=f'MASK_K{str(factor).replace(".","P")}'))
    fits.HDUList(hdus).writeto(out/'masks.fits.gz')
    if stamps:
        ph=fits.Header();ph['SCIVAL']=False;ph['PSFFIT']=False;ph['BGSUB']=False
        ph['COMMENT']='Raw masked 101x101 candidate stamps; no recentering or normalization.'
        fits.writeto(out/'psf_candidate_stamps.fits',np.asarray(stamps,np.float32),ph)
    kpp=float(cosmology.angular_diameter(frame['z']))*ARCSEC*1000*scale
    yy,xx=np.indices(image.shape,dtype=np.float32);bin_index=np.floor(np.hypot(xx-center[0],yy-center[1])*kpp/25).astype(np.int16)
    covers=[]
    for factor,mask in zip(FACTORS,masks):
        for b in range(20):
            used=int(np.count_nonzero((bin_index==b)&~mask));expected=math.pi*((b+1)**2-b*b)*25**2/kpp**2
            covers.append(dict(KRON_FACTOR=factor,R_MIN_KPC=b*25,R_MAX_KPC=(b+1)*25,N_VALID=used,COVERAGE=min(1,used/expected)))
    table(out/'mask_radial_coverage.tsv',covers,['KRON_FACTOR','R_MIN_KPC','R_MAX_KPC','N_VALID','COVERAGE'])
    from matplotlib import pyplot as plt
    step=max(1,math.ceil(max(image.shape)/650));small=image[::step,::step];finite=small[np.isfinite(small)]
    med=np.median(finite);std=max(1e-6,1.4826*np.median(abs(finite-med)))
    fig,axes=plt.subplots(2,2,figsize=(10,10))
    axes[0,0].imshow(np.arcsinh((small-med)/std),origin='lower',cmap='gray',vmin=-1,vmax=6)
    for star in stars: axes[0,0].plot((star['X_IMAGE']-1)/step,(star['Y_IMAGE']-1)/step,'go',ms=5,fillstyle='none')
    axes[0,0].plot(center[0]/step,center[1]/step,'r+',ms=12);axes[0,0].set_title(f"Raw display; PSF candidates {len(stars)}")
    for ax,mask,title in zip([axes[0,1],axes[1,0],axes[1,1]],[base,masks[0],masks[-1]],['Original quality + bright mask','Kron 2.5 exclusion','Kron 5.0 exclusion']):
        ax.imshow(mask[::step,::step],origin='lower',vmin=0,vmax=1,cmap='gray');ax.plot(center[0]/step,center[1]/step,'r+',ms=12);ax.set_title(title)
    fig.suptitle(f"Cluster {frame['id']} / {frame['preparation_cohort']} / engineering review")
    fig.tight_layout();fig.savefig(out/'review.png',dpi=110);plt.close(fig)
    ix,iy=round(center[0]),round(center[1]);x0,x1=max(0,ix-100),min(image.shape[1],ix+101);y0,y1=max(0,iy-100),min(image.shape[0],iy+101)
    fig,axes=plt.subplots(1,2,figsize=(10,5))
    axes[0].imshow(np.arcsinh((image[y0:y1,x0:x1]-med)/std),origin='lower',cmap='gray',vmin=-1,vmax=7)
    axes[1].imshow(seg[y0:y1,x0:x1],origin='lower',cmap='nipy_spectral',interpolation='nearest')
    for ax in axes:
        ax.plot(center[0]-x0,center[1]-y0,'r+',ms=12)
        close=(objects['x']>=x0)&(objects['x']<x1)&(objects['y']>=y0)&(objects['y']<y1)
        ax.plot(objects['x'][close]-x0,objects['y'][close]-y0,'co',ms=5,fillstyle='none')
    axes[0].set_title('BCG catalogue position + detections');axes[1].set_title('Detection segmentation; manual review')
    fig.tight_layout();fig.savefig(out/'bcg_review.png',dpi=110);plt.close(fig)
    record=dict(id=str(frame['id']),cohort=frame['preparation_cohort'],n_detections=len(objects),n_psf_candidates=len(stars),
        psf_status='provisional_inputs' if len(stars)>=4 else 'insufficient_candidates',psf_fwhm_locus_pix=mode,
        bcg_status='unique_positional_candidate' if bcg>=0 else 'missing_or_ambiguous',
        bcg_base_masked=bool(base[round(center[1]),round(center[0])]),n_members=len(members),n_members_in_frame=int(inside.sum()),
        member_fallback_count=fallback,minimum_coverage_by_kron={str(f):min(c['COVERAGE'] for c in covers if c['KRON_FACTOR']==f) for f in FACTORS},
        source_image=frame['image'],science_measurement=False,detection_background_subtracted_from_science=False,
        mask_ready_for_science=False,bcg_preservation_gate='candidate_requires_review' if bcg>=0 else 'blocked_missing_or_ambiguous_bcg',
        bright_mask_used_in_detection=False,bright_mask_retained_in_outputs=True,psf_wings_validated=False,
        psf_stamps_recentered=False,psf_stamps_normalized=False,psf_stamps_background_subtracted=False,**extra)
    write_json(out/'preparation.json',record)
    return record


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--manifest',type=Path,required=True)
    p.add_argument('--work',type=Path,required=True);p.add_argument('--limit',type=int)
    a=p.parse_args()
    if not platform.node().lower().startswith('lageunha'): p.error('Run on lageunha CPU')
    doc=json.loads(a.manifest.read_text());frames=doc['clusters'][:a.limit] if a.limit else doc['clusters']
    a.work.mkdir(parents=True,exist_ok=False)
    records=[];errors=[];started=time.time();cosmo=Cosmology(H0=doc['h0'],Om=doc['omega_m'])
    for i,frame in enumerate(frames,1):
        print(f'[{i}/{len(frames)}] preparing cluster {frame["id"]}',flush=True)
        try:
            record=prepare_frame(frame,cosmo,a.work/('cluster_'+str(frame['id'])));records.append(record)
            print(f'  detections={record["n_detections"]}, PSF candidates={record["n_psf_candidates"]}, BCG={record["bcg_status"]}',flush=True)
        except Exception as e:
            errors.append(dict(id=str(frame['id']),error_type=type(e).__name__,error=str(e)))
            print(f'  ERROR {type(e).__name__}: {e}',flush=True)
    write_json(a.work/'summary.json',dict(status='complete' if not errors else 'partial',n_requested=len(frames),n_prepared=len(records),
        elapsed_seconds=time.time()-started,errors=errors,records=records,science_measurement=False,
        manifest_sha256=hashlib.sha256(a.manifest.read_bytes()).hexdigest(),
        preparation_script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        software=dict(sep=sep.__version__,numpy=np.__version__),
        notes=['PSF stars are morphology-only candidates, not externally identified stars.',
               'No source-mask completeness, sky-transfer, PSF-wing or ICL inference certified.',
               'Kron factors 2.5/3.5/5.0 and star cuts are engineering sensitivity choices.',
               'Bright-star bits are excluded from detection only to inspect covered sources; never cleared from output masks.']))
    return 1 if errors else 0


if __name__=='__main__':sys.exit(main())
