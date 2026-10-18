#!/usr/bin/env python3
"""Forward reference PSFs through existing masks/backgrounds; diagnostic only."""
import argparse
from collections import Counter
import json
from pathlib import Path
import platform

import numpy as np
from astropy.io import fits

from fetch_ancillary_catalogues import sha
from diagnose_empirical_psf import read_table,measure,coordinates
import prepare_detection_psf as prep


def radial_template(rows,shape,x,y):
    rows=[r for r in rows if str(r['full_annulus'])=='True']
    radius,_=coordinates(shape,x,y)
    knots=np.array([float(r['mean_radius_pix']) for r in rows]);values=np.array([float(r['mean']) for r in rows])
    if not len(knots) or knots[0]!=0 or np.any(np.diff(knots)<=0) or radius.max()>knots[-1]:
        raise ValueError('Reference radial support is incomplete')
    # Documented derivative: annular averaging then native subpixel radial interpolation.
    # No clipping of signed reference pixels and no crop-to-unit-flux renormalization.
    return np.interp(radius,knots,values)


def forward_case(image,variance,mask,template,x,y,annulus):
    measured,radial=measure(image,variance,mask,x,y,annulus)
    if measured['status']!='diagnostic_only':return measured,[],[]
    model,model_radial=measure(template,np.full(template.shape,1e-30),mask,x,y,annulus)
    if model['status']!='diagnostic_only':raise ValueError('Reference model cannot pass the same geometric operator')
    radius,_=coordinates(image.shape,x,y)
    good=~mask&np.isfinite(image)&np.isfinite(variance)&(variance>0)
    uncorrected_core=float(template[good&(radius<12)].sum())
    if uncorrected_core<=0 or model['core_flux']<=0:raise ValueError('Reference core is not positive')
    scale=measured['core_flux']/model['core_flux'];profile_rows=[];injections=[]
    for observed,ref in zip(radial,model_radial):
        if observed['n_valid']!=ref['n_valid']:raise ValueError('Reference/data support mismatch')
        expected=scale*ref['sky_subtracted_mean'] if ref['n_valid'] else None
        profile_rows.append(dict(r_min_pix=ref['r_min_pix'],r_max_pix=ref['r_max_pix'],n_valid=ref['n_valid'],
            radial_coverage_pass=observed['radial_coverage_pass'],area_fraction=observed['area_fraction'],
            observed_per_core_flux=observed['mean_per_core_flux'],reference_per_core_flux=ref['mean_per_core_flux'],
            reference_raw_ring_mean=ref['raw_mean'],reference_subtracted_background=model['background'],
            reference_sky_removal_fraction=model['background']/ref['raw_mean'] if ref['raw_mean'] else None,
            observed_minus_reference=(observed['mean_per_core_flux']-ref['mean_per_core_flux']) if ref['n_valid'] else None,
            observed_raw_mean_sigma=observed['nominal_raw_mean_sigma'],expected_mean_in_native_flux=expected,
            expected_over_nominal_raw_sigma=expected/observed['nominal_raw_mean_sigma'] if ref['n_valid'] else None))
    # Paired deterministic additions retain real data, variance and masks. This tests
    # the extraction operator, not HSC upstream sky or re-detection completeness.
    for level in (.25,1.,4.):
        added,after=measure(image+scale*level*template,variance,mask,x,y,annulus)
        if added['status']!='diagnostic_only':raise ValueError('Positive injection unexpectedly fails')
        expected_core=level*measured['core_flux']
        difference=added['core_flux']-measured['core_flux']
        np.testing.assert_allclose(difference,expected_core,rtol=1e-10,atol=1e-10)
        for before,post,ref in zip(radial,after,model_radial):
            if not before['n_valid']:continue
            expected=scale*level*ref['sky_subtracted_mean']
            observed=post['sky_subtracted_mean']-before['sky_subtracted_mean']
            np.testing.assert_allclose(observed,expected,rtol=1e-8,atol=1e-10)
            injections.append(dict(relative_core_injection=level,r_min_pix=ref['r_min_pix'],r_max_pix=ref['r_max_pix'],
                predicted_increment=expected,measured_increment=observed,absolute_difference=abs(observed-expected),
                core_increment=difference,expected_core_increment=expected_core))
    measured.update(reference_core_transfer=model['core_flux']/uncorrected_core,
        reference_background=model['background'],reference_amplitude=scale,
        reference_approved=False,noise_added=False,masks_redetected=False)
    return measured,profile_rows,injections


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('empirical','reference','work'):p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args()
    if not platform.node().lower().startswith('lageunha'):p.error('Run on lageunha CPU')
    a.work.mkdir(parents=True,exist_ok=False)
    empirical=json.loads((a.empirical/'summary.json').read_text());reference=json.loads((a.reference/'summary.json').read_text())
    if empirical['script_sha256']!=sha(Path(__file__).with_name('diagnose_empirical_psf.py')):raise ValueError('Empirical extraction changed')
    assert not empirical['psf_approved'] and not reference['psf_approved']
    hashes={str(d/'summary.json'):sha(d/'summary.json') for d in (a.empirical,a.reference)}
    ledger=read_table(a.empirical/'star_ledger.tsv');previous=read_table(a.empirical/'measurements.tsv')
    for name in ('star_ledger.tsv','measurements.tsv'):hashes[str(a.empirical/name)]=sha(a.empirical/name)
    methods=[];profiles=[];injections=[];templates=[]
    for model in reference['records']:
        if not model['name'].endswith('.fits'):continue
        path=Path(model['path']);assert sha(path)==model['sha256'];hashes[str(path)]=model['sha256']
        radial_path=a.reference/(model['name']+'.radial.tsv');hashes[str(radial_path)]=sha(radial_path);curve=read_table(radial_path)
        for star in ledger:
            group=[r for r in previous if r['candidate_id']==star['candidate_id']]
            if not star['stamp']:
                for row in group:methods.append(dict(candidate_id=star['candidate_id'],cluster_id=star['cluster_id'],model=model['name'],
                    kron_factor=row['kron_factor'],background_r_min_pix=row['background_r_min_pix'],status=row['status'],psf_approved=False))
                continue
            stamp=Path(star['stamp']);assert sha(stamp)==star['stamp_sha256'];hashes[str(stamp)]=star['stamp_sha256']
            with fits.open(stamp,memmap=False) as h:
                if not np.isclose(float(star['pixel_scale_arcsec']),model['pixel_scale_arcsec'],atol=1e-9,rtol=0):raise ValueError('Pixel scale mismatch')
                image=h['SCIENCE'].data.astype(float);variance=h['VARIANCE'].data.astype(float)
                x,y=h[0].header['XCENTER'],h[0].header['YCENTER'];template=radial_template(curve,image.shape,x,y)
                # Small derived template saved separately with its construction metadata.
                tp=a.work/(model['name'].removesuffix('.fits')+'_star_'+star['candidate_id']+'.fits')
                hdr=fits.Header();hdr['SCIVAL']=False;hdr['PSFOK']=False;hdr['RADAVG']=True;hdr['PIXSCALE']=model['pixel_scale_arcsec']
                fits.writeto(tp,template,hdr);templates.append(dict(path=str(tp),sha256=sha(tp),candidate_id=star['candidate_id'],model=model['name']))
                for row in group:
                    factor=row['kron_factor'];ann=(int(row['background_r_min_pix']),int(row['background_r_max_pix']))
                    mask=h['MASK_K'+factor.replace('.','P')].data.astype(bool)
                    measured,radial,added=forward_case(image,variance,mask,template,x,y,ann)
                    assert measured['status']==row['status']
                    identity=dict(candidate_id=star['candidate_id'],cluster_id=star['cluster_id'],model=model['name'],
                        kron_factor=factor,background_r_min_pix=ann[0],background_r_max_pix=ann[1])
                    methods.append(dict(identity,**measured));profiles.extend(dict(**identity,**r) for r in radial)
                    injections.extend(dict(**identity,**r) for r in added)
        print(model['name'],'forward comparison complete',flush=True)
    for name,rows in (('method_status.tsv',methods),('profile_comparison.tsv',profiles),('differential_injections.tsv',injections)):
        prep.table(a.work/name,rows,list(dict.fromkeys(k for r in rows for k in r)))
    outer=[r for r in profiles if r['r_min_pix']==64 and r['radial_coverage_pass']]
    summaries=[]
    for name in sorted({r['model'] for r in methods}):
        rows=[r for r in outer if r['model']==name]
        summaries.append(dict(model=name,n_outer_coverage_rows=len(rows),
            n_outer_observed_nonpositive=sum(r['observed_per_core_flux']<=0 for r in rows),
            n_outer_expected_positive=sum(r['reference_per_core_flux']>0 for r in rows),
            sky_removal_fraction_percentiles=np.percentile([r['reference_sky_removal_fraction'] for r in rows],[0,50,100]).tolist(),
            expected_over_nominal_sigma_percentiles=np.percentile([r['expected_over_nominal_raw_sigma'] for r in rows],[0,50,100]).tolist()))
    prep.write_json(a.work/'summary.json',dict(status='complete',n_candidates=len(ledger),n_method_rows=len(methods),
        method_status_counts=dict(Counter(r['status'] for r in methods)),n_profile_rows=len(profiles),n_injection_rows=len(injections),
        max_differential_absolute_error=max(r['absolute_difference'] for r in injections),outer_summaries=summaries,
        templates=templates,source_sha256=hashes,script_sha256=sha(Path(__file__)),
        science_measurement=False,psf_approved=False,reference_approved=False,
        limitations=['Annular averaging removes the released 2D anisotropy; native radial interpolation is a documented model transformation.',
            'Symmetric/asymmetric products share input stars and are not independent reference measurements.',
            'Only the existing good observed methods contribute profiles; all original failures remain in the ledger.',
            'Core normalization uses the same background and mask operator on observation and reference.',
            'No position/epoch/seeing match or fitted spatial PSF model; comparison does not certify reference applicability.',
            'Deterministic injection holds variance and masks fixed and adds no photon noise; it validates the downstream operator only.',
            'Nominal sigma omits sky, core normalization and pixel covariance; ratios are not detection significances.',
            'Reference encircled flux is relative to its finite published grid, not an infinite PSF integral.']))
    plot(profiles,a.work/'reference_comparison.png')
    print(json.dumps(dict(methods=dict(Counter(r['status'] for r in methods)),outer=summaries,
        n_injection_rows=len(injections),max_error=max(r['absolute_difference'] for r in injections))))


def plot(rows,path):
    from matplotlib import pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(12,4.5));ids=sorted({r['candidate_id'] for r in rows});colors=plt.cm.tab20(np.linspace(0,1,len(ids)))
    for factor,ax in zip(('2.5','3.5'),axes):
        reference=[]
        for sid,color in zip(ids,colors):
            for ann,style in ((120,'-'),(160,'--')):
                group=[r for r in rows if r['model']=='HSC_PSF_i.fits' and r['candidate_id']==sid and r['kron_factor']==factor and r['background_r_min_pix']==ann and r['radial_coverage_pass']]
                if not group:continue
                x=[(r['r_min_pix']+r['r_max_pix'])/2*.168 for r in group]
                ax.plot(x,[r['observed_per_core_flux'] for r in group],style,color=color,alpha=.45,lw=.8)
                reference.extend(group)
        for edge in sorted({(r['r_min_pix'],r['r_max_pix']) for r in reference}):
            ys=[r['reference_per_core_flux'] for r in reference if (r['r_min_pix'],r['r_max_pix'])==edge]
            ax.plot(np.mean(edge)*.168,np.median(ys),'ko',ms=3,label='Forward reference median' if edge==(0,2) else None)
        ax.axhline(0,c='gray',lw=.5);ax.set_yscale('symlog',linthresh=1e-6);ax.set_xlabel('Radius (arcsec)')
        ax.set_ylabel('Signed mean / measured core flux');ax.set_title('Kron '+factor);ax.legend(fontsize=8)
    fig.suptitle('Observed stars and published radial reference through the same masks/sky; diagnostic only',fontsize=10)
    fig.tight_layout();fig.savefig(path,dpi=140);plt.close(fig)


if __name__=='__main__':main()
