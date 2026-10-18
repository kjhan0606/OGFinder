#!/usr/bin/env python3
"""Check known star support against the existing 800 cutouts using native WCS."""
import argparse
import hashlib
import json
from pathlib import Path
import platform

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS
import astropy.units as u

from diagnose_empirical_psf import read_table
from fetch_ancillary_catalogues import sha
import prepare_detection_psf as prep


def full_stamp(x,y,shape,half=200):
    if not np.isfinite([x,y]).all():return False
    ix,iy=round(float(x)),round(float(y))
    return ix-half>=0 and iy-half>=0 and ix+half<shape[1] and iy+half<shape[0]


def source_path(manifest_path,image):
    path=Path(image)
    if path.is_absolute():return path
    if len(path.parts)!=1:raise ValueError('Unexpected relative downloader image path')
    return manifest_path.parent/'images'/path


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('sample','work'):p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args()
    if not platform.node().lower().startswith('lageunha'):p.error('Run on lageunha CPU')
    a.work.mkdir(parents=True,exist_ok=False)
    manifest_path=a.sample/'download_01/download_manifest.json';status_path=a.sample/'download_01/download_status.json'
    candidate_path=a.sample/'psf_cds_2026-10-13_01/evidence/candidates.json'
    ledger_path=a.sample/'empirical_psf_2026-10-14_01/star_ledger.tsv'
    frames=json.loads(manifest_path.read_text())['clusters'];status={str(f['id']):f for f in json.loads(status_path.read_text())['completed']}
    requests=json.loads(candidate_path.read_text());stars=read_table(ledger_path)
    world=SkyCoord([requests[int(s['candidate_id'])]['ra_deg'] for s in stars]*u.deg,
        [requests[int(s['candidate_id'])]['dec_deg'] for s in stars]*u.deg)
    pairs=[];headers=[]
    for frame in frames:
        cid=str(frame['id']);review=status[cid]['review']
        center=SkyCoord(frame['ra_deg']*u.deg,frame['dec_deg']*u.deg)
        # Generous candidate circle; exact support is decided only by WCS and dimensions.
        bound=np.sqrt(2)*frame['requested_half_width_arcsec']+60
        nearby=np.flatnonzero(center.separation(world).arcsec<=bound)
        if not len(nearby):continue
        path=source_path(manifest_path,frame['image'])
        with fits.open(path,memmap=True,lazy_load_hdus=True) as h:
            header=h[review['science_hdu']].header;wcs=WCS(header).celestial
            shape=(header['NAXIS2'],header['NAXIS1']);assert list(shape)==review['shape']
            scale=3600*np.sqrt(abs(np.linalg.det(wcs.pixel_scale_matrix)))
            np.testing.assert_allclose(scale,review['pixel_scale_arcsec'],atol=1e-9,rtol=0)
            x,y=wcs.all_world2pix(world[nearby].ra.deg,world[nearby].dec.deg,0)
            dates={str(i)+':'+k:v.header[k] for i,v in enumerate(h) for k in v.header if 'DATE' in k.upper() or 'MJD' in k.upper()}
            if len(h)!=4:raise ValueError('Unexpected additional HDUs: inspect PSF/epoch payload before proceeding')
            header_record=dict(cluster_id=cid,path=str(path),shape=list(shape),date_keys=dates,
                psf_archive_identifier=h[0].header.get('ARCHIVE_ID_PSF'),hdu_count=len(h),
                image_header_sha256=hashlib.sha256(header.tostring().encode()).hexdigest(),
                actual_psf_model_present=False,effective_epoch_verified=False)
            headers.append(header_record)
            for index,px,py in zip(nearby,x,y):
                star=stars[index]
                pairs.append(dict(candidate_id=star['candidate_id'],original_cluster_id=star['cluster_id'],candidate_frame_id=cid,
                    is_original=cid==star['cluster_id'],x_pixel=float(px),y_pixel=float(py),
                    center_inside=bool(0<=px<shape[1] and 0<=py<shape[0]),
                    full_401_stamp=full_stamp(px,py,shape),source_image=str(path),pixel_quality_verified=False))
    recovery=[]
    for star in stars:
        own=next(r for r in pairs if r['candidate_id']==star['candidate_id'] and r['is_original'])
        assert own['full_401_stamp']!=(star['status']=='blocked_stamp_edge')
        alternatives=[r['candidate_frame_id'] for r in pairs if r['candidate_id']==star['candidate_id'] and not r['is_original'] and r['full_401_stamp']]
        recovery.append(dict(candidate_id=star['candidate_id'],original_cluster_id=star['cluster_id'],
            originally_blocked_by_edge=star['status']=='blocked_stamp_edge',alternative_frame_ids=','.join(alternatives),n_alternatives=len(alternatives)))
    prep.table(a.work/'footprint_pairs.tsv',pairs,list(pairs[0]));prep.table(a.work/'recovery_candidates.tsv',recovery,list(recovery[0]))
    prep.write_json(a.work/'header_evidence.json',headers)
    prep.write_json(a.work/'summary.json',dict(status='complete',n_cutouts_screened=len(frames),n_known_stars=len(stars),
        n_nearby_cutout_headers=len(headers),n_wcs_pairs=len(pairs),n_original_full_stamps=sum(bool(s['stamp']) for s in stars),
        n_edge_stars_with_alternative=sum(r['originally_blocked_by_edge'] and r['n_alternatives']>0 for r in recovery),
        n_stars_with_alternative=sum(r['n_alternatives']>0 for r in recovery),
        n_headers_with_date_keys=sum(bool(h['date_keys']) for h in headers),
        n_headers_with_psf_archive_id=sum(h['psf_archive_identifier'] is not None for h in headers),
        source_sha256={str(p):sha(p) for p in (manifest_path,status_path,candidate_path,ledger_path)},
        script_sha256=sha(Path(__file__)),science_measurement=False,psf_approved=False,
        limitations=['Existing manifest/status supply the angular screening and dimensions; full image hashes are not re-read in this header inventory.',
                     'Alternative support is geometry only; pixels and local detection masks still require validation.',
                     'A PSF archive integer is not a PSF image or serialized model; no model is reconstructed from it.',
                     'Date keywords, when present, are evidence inventory only; an effective coadd epoch requires contributing visits/weights.']))
    print(json.dumps(dict(cutouts=len(frames),stars=len(stars),nearby_headers=len(headers),
        edge_recoverable=sum(r['originally_blocked_by_edge'] and r['n_alternatives']>0 for r in recovery),
        headers_with_dates=sum(bool(h['date_keys']) for h in headers))))


if __name__=='__main__':main()
