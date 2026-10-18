#!/usr/bin/env python3
"""Native-pixel Gaia-supported star diagnostics, not a certified PSF model.

All signed profiles use a measured finite core aperture, never total PSF flux.
Background/normalization/correlated-noise uncertainties are not calibrated.
"""
import argparse
from collections import Counter, defaultdict
import csv
import json
from pathlib import Path
import platform
import time

import numpy as np
from astropy.io import fits

from fetch_ancillary_catalogues import sha
from download_pilot import BAD_BIT_POSITIONS, BRIGHT_OBJECT_BIT
import prepare_detection_psf as prep

EDGES = (0,2,4,6,8,12,16,24,32,48,64,96)
ANNULI = ((120,150),(160,190))
FACTORS = (2.5,3.5)


def read_table(path):
    with path.open() as f: return list(csv.DictReader(f,delimiter='\t'))


def coordinates(shape, x, y):
    yy,xx=np.indices(shape)
    return np.hypot(xx-x,yy-y),np.minimum(7,((np.arctan2(yy-y,xx-x)+np.pi)*4/np.pi).astype(int))


def measure(image, variance, mask, x, y, annulus):
    radius,sector=coordinates(image.shape,x,y)
    good=~mask & np.isfinite(image) & np.isfinite(variance) & (variance>0)
    sky=(radius>=annulus[0])&(radius<annulus[1]); sky_good=sky&good
    means=[]; sector_counts=[]; sector_fractions=[]
    for k in range(8):
        region=sky&(sector==k); valid=region&good; n=int(valid.sum())
        sector_counts.append(n); sector_fractions.append(n/max(1,int(region.sum())))
        if n: means.append(float(image[valid].mean()))
    sky_fraction=float(sky_good.sum()/max(1,sky.sum()))
    core=radius<12; core_good=core&good; core_fraction=float(core_good.sum()/core.sum())
    result=dict(status='diagnostic_only',background_r_min_pix=annulus[0],background_r_max_pix=annulus[1],
        sky_n_valid=int(sky_good.sum()),sky_area_fraction=sky_fraction,
        sky_min_sector_n=min(sector_counts),sky_min_sector_fraction=min(sector_fractions),
        core_n_valid=int(core_good.sum()),core_area_fraction=core_fraction,
        background=None,background_sector_rms=None,core_flux=None,core_snr_nominal=None,
        psf_approved=False,science_measurement=False)
    if sky_good.sum()<300 or sky_fraction<.3 or min(sector_counts)<20 or min(sector_fractions)<.25:
        result['status']='blocked_background_coverage';return result,[]
    # Equal azimuth weights, no clipping: background choice remains a diagnostic assumption.
    background=float(np.mean(means));scatter=float(np.std(means,ddof=1))
    flux=float(np.sum(image[core_good]-background));noise=float(np.sqrt(variance[core_good].sum()))
    result.update(background=background,background_sector_rms=scatter,core_flux=flux,core_snr_nominal=flux/noise if noise else None)
    if core_fraction<.99:
        result['status']='blocked_core_mask';return result,[]
    if not np.isfinite(flux) or flux<=0 or noise<=0 or flux/noise<30:
        result['status']='blocked_core_flux';return result,[]
    profiles=[]
    for lo,hi in zip(EDGES[:-1],EDGES[1:]):
        ring=(radius>=lo)&(radius<hi);valid=ring&good;n=int(valid.sum())
        fraction=float(n/max(1,ring.sum()))
        sectors=sum(np.count_nonzero(valid&(sector==k))>=max(1,.25*np.count_nonzero(ring&(sector==k))) for k in range(8))
        raw=float(image[valid].mean()) if n else None
        residual=raw-background if n else None
        profiles.append(dict(r_min_pix=lo,r_max_pix=hi,n_valid=n,area_fraction=fraction,n_covered_sectors=sectors,
            radial_coverage_pass=fraction>=.5 and sectors==8,
            raw_mean=raw,sky_subtracted_mean=residual,
            mean_per_core_flux=residual/flux if n else None,
            observed_ring_flux_per_core_flux=residual*n/flux if n else None,
            nominal_raw_mean_sigma=float(np.sqrt(variance[valid].sum())/n) if n else None))
    return result,profiles


def objects_from_rows(rows, x0, y0):
    if [int(r['NUMBER']) for r in rows]!=list(range(1,len(rows)+1)):
        raise ValueError('Detection numbering changed')
    objects=np.zeros(len(rows),dtype=[(s,float) for s in ('x','y','a','b','theta')])
    for key,column,shift in (('x','X_IMAGE',1+x0),('y','Y_IMAGE',1+y0),('a','A_IMAGE',0),('b','B_IMAGE',0),('theta','THETA_RAD',0)):
        objects[key]=[float(r[column])-shift for r in rows]
    return objects,np.array([float(r['KRON_RADIUS']) for r in rows]),np.array([int(r['KRON_FLAGS']) for r in rows])


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('sample','work'):p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args()
    if not platform.node().lower().startswith('lageunha'):p.error('Run on lageunha CPU')
    a.work.mkdir(parents=True,exist_ok=False); started=time.time();hashes={}
    def record(path):
        path=Path(path);hashes[str(path)]=sha(path);return path
    evidence_dir=a.sample/'psf_cds_2026-10-13_01/evidence'
    evidence=json.loads(record(evidence_dir/'evidence_summary.json').read_text())
    requests=json.loads(record(evidence_dir/'candidates.json').read_text())
    gs=json.loads(record(evidence_dir/'status.json').read_text())
    gp=record(evidence_dir/'gaia_matches.csv')
    if hashes[str(gp)]!=gs['response_sha256']:raise ValueError('Gaia response changed')
    with gp.open() as f: gaia={int(r['candidate_id']):r for r in csv.DictReader(f)}
    selected=[r for r in evidence['records'] if r['status']=='stellar_evidence_requires_epoch_and_psf_review']
    manifest=json.loads(record(a.sample/'validation_sample_2026-10-13_01/validation_manifest.json').read_text())
    frames={str(r['id']):r for r in manifest['clusters']}
    ds=json.loads(record(a.sample/'download_01/download_status.json').read_text())
    download={str(r['id']):r for r in ds['completed']}
    groups=defaultdict(list)
    for star in selected:groups[str(star['cluster_id'])].append(star)
    measures=[];profiles=[];star_records=[]
    for cid,stars in groups.items():
        frame=frames[cid];raw=record(frame['image']);review=download[cid]['review']
        if hashes[str(raw)]!=download[cid]['sha256']:raise ValueError('Original image hash changed')
        source=a.sample/'detection_psf_2026-10-13_01'/('cluster_'+cid)
        detections=read_table(record(source/'detections.tsv'));mp=record(source/'masks.fits.gz')
        with fits.open(raw,memmap=False) as h,fits.open(mp,memmap=False) as mh:
            image=np.asarray(h[review['science_hdu']].data,dtype=float)
            variance=np.asarray(h[review['variance_hdu']].data,dtype=float)
            bits=np.asarray(h[review['mask_hdu']].data)
            seg=mh['SEGMENTATION'].data;original_base=mh['BASE_MASK'].data.astype(bool)
            base=(bits&(sum(1<<b for b in BAD_BIT_POSITIONS)|BRIGHT_OBJECT_BIT))!=0
            base|=~np.isfinite(image)|~np.isfinite(variance)|(variance<=0)
            if not np.array_equal(base,original_base):raise ValueError('Base quality mask mismatch')
            for star in stars:
                sid=star['candidate_id'];number=star['number'];d=detections[number-1];g=gaia[sid];req=requests[sid]
                if str(req['cluster_id'])!=cid or req['number']!=number or g['source_id']!=star['source_ids'][0]:
                    raise ValueError('Gaia/HSC candidate identity mismatch')
                np.testing.assert_allclose([float(d['RA']),float(d['DEC'])],[req['ra_deg'],req['dec_deg']],rtol=0,atol=1e-9)
                x,y=float(d['X_IMAGE'])-1,float(d['Y_IMAGE'])-1;ix,iy=round(x),round(y)
                x0,y0=ix-200,iy-200;scale=review['pixel_scale_arcsec']
                common=dict(candidate_id=sid,cluster_id=cid,number=number,source_id=g['source_id'],
                    gaia_variability=g['phot_variable_flag'],epoch_review_required=True,pixel_scale_arcsec=scale)
                sr=dict(**common,x0=x0,y0=y0,native_x=x,native_y=y,stamp=None,status='diagnostic_only',psf_approved=False)
                block=None
                if x0<0 or y0<0 or ix+200>=image.shape[1] or iy+200>=image.shape[0]:block='blocked_stamp_edge'
                elif str(g['phot_variable_flag']).upper()=='VARIABLE':block='blocked_known_variable'
                if block:
                    sr['status']=block
                    for factor in FACTORS:
                        for ann in ANNULI:measures.append(dict(**common,kron_factor=factor,background_r_min_pix=ann[0],background_r_max_pix=ann[1],status=block,psf_approved=False,science_measurement=False))
                    star_records.append(sr);continue
                sl=np.s_[y0:y0+401,x0:x0+401];cut=image[sl].copy();var=variance[sl].copy();bas=base[sl].copy();labels=seg[sl].copy()
                objects,kron,kflags=objects_from_rows(detections,x0,y0)
                masks,_=prep.source_masks(objects,labels,kron,kflags,bas,number-1,np.empty((0,2)),scale)
                assert all(np.all(mask[bas]) for mask in masks) and np.all(~masks[0]|masks[1])
                path=a.work/f'star_{sid}.fits';head=fits.Header()
                head['SCIVAL']=False;head['PSFOK']=False;head['X0']=x0;head['Y0']=y0;head['XCENTER']=x-x0;head['YCENTER']=y-y0
                head['PIXSCALE']=scale;head['CANDID']=sid;head['SOURCEID']=g['source_id']
                fits.HDUList([fits.PrimaryHDU(header=head),fits.ImageHDU(cut,name='SCIENCE'),fits.ImageHDU(var,name='VARIANCE'),
                    fits.ImageHDU(bas.astype('u1'),name='BASE_MASK'),fits.ImageHDU(labels,name='SEGMENTATION')]+
                    [fits.ImageHDU(mask.astype('u1'),name='MASK_K'+str(f).replace('.','P')) for f,mask in zip(FACTORS,masks)]).writeto(path)
                sr.update(stamp=str(path),stamp_sha256=sha(path));star_records.append(sr)
                for factor,mask in zip(FACTORS,masks):
                    for ann in ANNULI:
                        stats,radial=measure(cut,var,mask,x-x0,y-y0,ann)
                        identity=dict(**common,kron_factor=factor,background_r_min_pix=ann[0],background_r_max_pix=ann[1])
                        measures.append(dict(identity,**stats))
                        profiles.extend(dict(**identity,**r) for r in radial)
        print('cluster',cid,'stars',len(stars),flush=True)
    for name,rows in (('star_ledger.tsv',star_records),('measurements.tsv',measures),('radial_profiles.tsv',profiles)):
        if rows:prep.table(a.work/name,rows,list(dict.fromkeys(k for r in rows for k in r)))
    good={r['candidate_id'] for r in measures if r['status']=='diagnostic_only'}
    all_four=[sid for sid in good if sum(r['candidate_id']==sid and r['status']=='diagnostic_only' for r in measures)==4]
    by_key={(r['candidate_id'],r['kron_factor'],r['background_r_min_pix']):r for r in measures}
    sensitivity=[]
    for sid in all_four:
        for factor in FACTORS:
            first=by_key[sid,factor,120];second=by_key[sid,factor,160]
            sensitivity.append(dict(candidate_id=sid,kron_factor=factor,
                background_change=second['background']-first['background'],
                fractional_core_change=second['core_flux']/first['core_flux']-1))
    if sensitivity:prep.table(a.work/'background_sensitivity.tsv',sensitivity,list(sensitivity[0]))
    outer=[r for r in profiles if r['r_min_pix']==64 and r['radial_coverage_pass']]
    prep.write_json(a.work/'summary.json',dict(status='complete',n_candidates=len(selected),n_frames=len(groups),
        n_stamps=sum(r['stamp'] is not None for r in star_records),n_method_rows=len(measures),n_profile_rows=len(profiles),
        method_status_counts=dict(Counter(r['status'] for r in measures)),n_stars_any_valid=len(good),n_stars_all_four_valid=len(all_four),
        n_outer_bins_coverage_pass=len(outer),n_outer_bins_nonpositive=sum(r['sky_subtracted_mean']<=0 for r in outer),
        median_absolute_core_change=float(np.median([abs(r['fractional_core_change']) for r in sensitivity])) if sensitivity else None,
        max_absolute_core_change=max((abs(r['fractional_core_change']) for r in sensitivity),default=None),
        science_measurement=False,psf_approved=False,script_sha256=sha(Path(__file__)),source_sha256=hashes,
        stamp_side_pixels=401,core_radius_pixels=12,radial_edges_pixels=EDGES,background_annuli_pixels=ANNULI,
        elapsed_seconds=time.time()-started,
        limitations=['Gaia 2016 coordinates have not been propagated to the HSC effective coadd epoch.',
            'Unknown Gaia variability remains unknown; only known VARIABLE entries are blocked.',
            'Measured native HSC center; no interpolation, extrapolated wings or unit-total-flux normalization.',
            'Background annuli may contain stellar wings, undetected sources and cluster light.',
            'Signed ring means retain negatives; coverage tests are geometry checks, not detections.',
            'Variance errors condition on the measured background and core; no coadd covariance or sky/core uncertainty calibration.',
            'Changing masks changes the retained pixels. Between-mask differences are not pure contaminant removal.',
            'No universal PSF model is built; these stars do not establish PSF coverage of the 800-cluster sample.']))
    plot(profiles,a.work/'radial_review.png')
    print(json.dumps(dict(candidates=len(selected),stamps=sum(r['stamp'] is not None for r in star_records),
        methods=dict(Counter(r['status'] for r in measures)),stars_all_four=len(all_four))))


def plot(rows,path):
    from matplotlib import pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(12,4.5))
    for factor,ax in zip(FACTORS,axes):
        for sid in sorted({r['candidate_id'] for r in rows}):
            for ann,style in ((120,'-'),(160,'--')):
                rr=[r for r in rows if r['candidate_id']==sid and r['kron_factor']==factor and r['background_r_min_pix']==ann and r['radial_coverage_pass']]
                if rr:ax.plot([(r['r_min_pix']+r['r_max_pix'])/2*r['pixel_scale_arcsec'] for r in rr],
                    [r['mean_per_core_flux'] for r in rr],style,alpha=.35,lw=.8)
        ax.set_yscale('symlog',linthresh=1e-6);ax.axhline(0,c='gray',lw=.6)
        ax.set_xlabel('Native radius (arcsec)');ax.set_ylabel('Signed mean / measured r<12 px core flux')
        ax.set_title(f'Kron {factor}; solid sky 120-150 px / dashed 160-190 px')
    fig.suptitle('Empirical star diagnostics; background and PSF approval pending',fontsize=11)
    fig.tight_layout();fig.savefig(path,dpi=140);plt.close(fig)


if __name__=='__main__':main()
