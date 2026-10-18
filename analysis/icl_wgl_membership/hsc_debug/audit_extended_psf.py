#!/usr/bin/env python3
"""Re-read published FITS and replay extraction without the measurement function."""
import argparse
import hashlib
import json
from pathlib import Path
import platform

import numpy as np
from astropy.io import fits

from fetch_ancillary_catalogues import sha
from diagnose_empirical_psf import read_table
import prepare_detection_psf as prep


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('reference','comparison','empirical','support','work'):p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args()
    if not platform.node().lower().startswith('lageunha'):p.error('Run on lageunha CPU')
    a.work.mkdir(parents=True,exist_ok=False)
    ref=json.loads((a.reference/'summary.json').read_text());comp=json.loads((a.comparison/'summary.json').read_text())
    support=json.loads((a.support/'summary.json').read_text())
    assert not ref['psf_approved'] and not comp['psf_approved'] and not support['psf_approved']
    hashes=dict(comp['source_sha256'],**support['source_sha256'])
    for path,digest in hashes.items():assert sha(Path(path))==digest,path
    for name,report in (('prepare_extended_psf.py',ref),('compare_extended_psf.py',comp),('inventory_psf_support.py',support)):
        assert sha(Path(__file__).with_name(name))==report['script_sha256']
    assert sha(a.reference/'zenodo_record.json')==ref['metadata_sha256']
    templates={(r['model'],r['candidate_id']):r for r in comp['templates']}
    methods=read_table(a.comparison/'method_status.tsv');profiles=read_table(a.comparison/'profile_comparison.tsv')
    additions=read_table(a.comparison/'differential_injections.tsv');stars={r['candidate_id']:r for r in read_table(a.empirical/'star_ledger.tsv')}
    assert len(methods)==comp['n_method_rows'] and len(profiles)==comp['n_profile_rows'] and len(additions)==comp['n_injection_rows']
    replay=0;injection_replay=0;model_checks=[];max_template=0.;max_increment=0.
    for model in ref['records']:
        path=Path(model['path']);assert sha(path)==model['sha256']
        assert hashlib.md5(path.read_bytes()).hexdigest()==model['published_md5']
        if not model['name'].endswith('.fits'):continue
        with fits.open(path,memmap=False) as h:
            data=np.asarray(h[model['hdu']].data,dtype=float)
            assert data.shape==(4001,4001) and np.isfinite(data).all()
            np.testing.assert_allclose(data.sum(),model['finite_grid_sum'],rtol=1e-12)
            assert np.count_nonzero(data<0)==model['negative_pixels']
            # Reconstruct the local radial curve by direct ring selection, without bincount.
            cut=data[1700:2301,1700:2301];yy,xx=np.indices(cut.shape);rr=np.sqrt((xx-300)**2+(yy-300)**2)
            radii=[];values=[]
            for i in range(286):
                chosen=(rr>=i)&(rr<i+1);radii.append(float(rr[chosen].mean()));values.append(float(cut[chosen].mean()))
            full_y,full_x=np.ogrid[-2000:2001,-2000:2001];r2=full_x**2+full_y**2
            for radius,fraction in model['encircled_grid_fractions'].items():
                actual=data[r2<(float(radius)/model['pixel_scale_arcsec'])**2].sum()/data.sum()
                np.testing.assert_allclose(actual,fraction,rtol=1e-11)
        model_checks.append(dict(model=model['name'],published_md5_verified=True,negative_pixels=model['negative_pixels']))
        relevant=[r for r in methods if r['model']==model['name'] and r['status']=='diagnostic_only']
        for row in relevant:
            sid=row['candidate_id'];star=stars[sid];tp=templates[model['name'],sid]
            assert sha(Path(tp['path']))==tp['sha256']
            with fits.open(star['stamp'],memmap=False) as h,fits.open(tp['path'],memmap=False) as mh:
                image=h['SCIENCE'].data.astype(float);variance=h['VARIANCE'].data.astype(float)
                mask=h['MASK_K'+row['kron_factor'].replace('.','P')].data.astype(bool)
                yy,xx=np.indices(image.shape);dx=xx-h[0].header['XCENTER'];dy=yy-h[0].header['YCENTER']
                radius=np.sqrt(dx*dx+dy*dy);sector=np.minimum(7,np.floor((np.angle(dx+1j*dy)+np.pi)*4/np.pi).astype(int))
                template=np.interp(radius,radii,values);delta=float(np.max(abs(template-mh[0].data)));max_template=max(max_template,delta)
                assert delta<1e-12
                good=~mask&np.isfinite(image)&np.isfinite(variance)&(variance>0)
                sky=good&(radius>=int(row['background_r_min_pix']))&(radius<int(row['background_r_max_pix']))
                bg=np.mean([image[sky&(sector==k)].mean() for k in range(8)])
                mbg=np.mean([template[sky&(sector==k)].mean() for k in range(8)])
                core=good&(radius<12);flux=np.sum(image[core]-bg);mf=np.sum(template[core]-mbg);amplitude=flux/mf
                np.testing.assert_allclose([flux,mbg,amplitude],[float(row['core_flux']),float(row['reference_background']),float(row['reference_amplitude'])],rtol=1e-10,atol=1e-12)
                subset=[r for r in profiles if all(r[k]==row[k] for k in ('model','candidate_id','kron_factor','background_r_min_pix'))]
                assert len(subset)==11
                for profile in subset:
                    ring=good&(radius>=int(profile['r_min_pix']))&(radius<int(profile['r_max_pix']))
                    assert ring.sum()==int(profile['n_valid'])
                    raw=image[ring].mean();mraw=template[ring].mean()
                    np.testing.assert_allclose([(raw-bg)/flux,(mraw-mbg)/mf,mbg/mraw],
                        [float(profile['observed_per_core_flux']),float(profile['reference_per_core_flux']),float(profile['reference_sky_removal_fraction'])],rtol=1e-9,atol=1e-13)
                    replay+=1
                    injected=[r for r in additions if all(r[k]==profile[k] for k in ('model','candidate_id','kron_factor','background_r_min_pix','r_min_pix'))]
                    assert len(injected)==3
                    for item in injected:
                        level=float(item['relative_core_injection']);addition=amplitude*level
                        actual=(image[ring]+addition*template[ring]).mean()-raw-addition*mbg
                        predicted=addition*(mraw-mbg)
                        max_increment=max(max_increment,abs(actual-predicted))
                        np.testing.assert_allclose([actual,predicted],[float(item['measured_increment']),float(item['predicted_increment'])],rtol=1e-8,atol=1e-10)
                        injection_replay+=1
    assert replay==len(profiles) and injection_replay==len(additions)
    result=dict(passed=True,input_hashes_verified=len(hashes),model_checks=model_checks,
        radial_comparisons_replayed=replay,differential_injections_replayed=injection_replay,
        max_template_pixel_difference=max_template,max_independent_injection_error=max_increment,
        script_sha256=sha(Path(__file__)),science_measurement=False,psf_approved=False,
        report_sha256={str(d/'summary.json'):sha(d/'summary.json') for d in (a.reference,a.comparison,a.support)})
    prep.write_json(a.work/'audit.json',result);print(json.dumps(result))


if __name__=='__main__':main()
