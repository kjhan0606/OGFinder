#!/usr/bin/env python3
"""Profile the existing 50-cluster masks for paired pre-stack covariance checks."""
import argparse
import json
from pathlib import Path
import platform

import numpy as np
from astropy.cosmology import FlatLambdaCDM
from astropy.io import fits

import prepare_detection_psf as prep
from download_pilot import file_sha256 as sha, validate_fits
from validate_control_images import nJy_factor
from validate_prestack import surface_result, vectors
from prestack_core import geometry, measure, FACTORS


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--sample',type=Path,required=True);p.add_argument('--work',type=Path,required=True)
    a=p.parse_args()
    if not platform.node().lower().startswith('lageunha'):p.error('Run on lageunha CPU')
    a.work.mkdir(parents=True,exist_ok=False)
    vp=a.sample/'validation_sample_2026-10-13_01/validation_manifest.json'
    doc=json.loads(vp.read_text());cosmo=FlatLambdaCDM(H0=doc['h0'],Om0=doc['omega_m'],Tcmb0=0)
    old_path=a.sample/'outer_sky_diagnostics_2026-10-14_02/summary.json'
    old=json.loads(old_path.read_text());hashes={str(vp):sha(vp),str(old_path):sha(old_path)}
    arrays=[];records=[]
    for frame in doc['clusters']:
        path=Path(frame['image']);mp=a.sample/'bcg_masks_2026-10-13_01'/('cluster_'+str(frame['id']))/'proposed_masks.fits.gz'
        for artifact in (path,mp):
            digest=sha(artifact)
            if old['source_sha256'][str(artifact)]!=digest:raise ValueError('Changed cluster input')
            hashes[str(artifact)]=digest
        info=validate_fits(path,frame);factor=nJy_factor(info);center=info['center_pixel_xy']
        kpp=cosmo.kpc_proper_per_arcmin(frame['z']).value*info['pixel_scale_arcsec']/60
        with fits.open(path,memmap=False) as h,fits.open(mp,memmap=False) as mh:
            image=h[info['science_hdu']].data.astype('f4');var=h[info['variance_hdu']].data.astype('f4')
            geom=geometry(image.shape,center,kpp)
            rows=[surface_result(measure(image,var,mh['MASK_K'+str(f).replace('.','P')].data.astype(bool),geom,center,kpp),factor) for f in FACTORS]
        arrays.append(vectors(rows));records.append(dict(id=str(frame['id']),ra_deg=frame['ra_deg'],dec_deg=frame['dec_deg'],cohort=frame['preparation_cohort']))
    payload={key:np.array([r[key] for r in arrays]) for key in arrays[0]}
    out=a.work/'profiles.npz';np.savez_compressed(out,**payload)
    prep.write_json(a.work/'summary.json',dict(status='complete',n_clusters=len(records),records=records,source_sha256=hashes,
        output_sha256=sha(out),script_sha256={str(Path(__file__)):sha(Path(__file__)),str(Path(__file__).with_name('prestack_core.py')):sha(Path(__file__).with_name('prestack_core.py'))},
        science_measurement=False,sky_correction_approved=False))
    print('Existing cluster profiles complete:',len(records),flush=True)


if __name__=='__main__':main()
