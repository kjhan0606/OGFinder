#!/usr/bin/env python3
"""Rebuild reviewed masks without treating the uniquely matched BCG as a satellite."""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import platform
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
import prepare_detection_psf as prep


def satellite_members(members, center, scale, stable):
    """Exclude only one uniquely position-matched member after BCG association."""
    near=np.flatnonzero(np.hypot(members[:,0]-center[0],members[:,1]-center[1])*scale<=1.)
    if stable and len(near)==1:
        return np.delete(members,near[0],axis=0),1
    return members,0


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('manifest','prepared','reviewed','work'):p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args()
    if not platform.node().lower().startswith('lageunha'):p.error('Run on lageunha CPU')
    review=json.loads((a.reviewed/'summary.json').read_text())
    assert review['status']=='complete'
    assert hashlib.sha256(a.manifest.read_bytes()).hexdigest()==review['manifest_sha256']
    frames={str(f['id']):f for f in json.loads(a.manifest.read_text())['clusters']}
    a.work.mkdir(parents=True,exist_ok=False);records=[]
    for r in review['records']:
        cid=r['id'];d=a.prepared/('cluster_'+cid);out=a.work/('cluster_'+cid);out.mkdir()
        with (d/'detections.tsv').open() as f:rows=list(csv.DictReader(f,delimiter='\t'))
        objects=np.zeros(len(rows),dtype=[(k,'f8') for k in ('x','y','a','b','theta')])
        for k,col in [('x','X_IMAGE'),('y','Y_IMAGE'),('a','A_IMAGE'),('b','B_IMAGE'),('theta','THETA_RAD')]:
            objects[k]=[float(row[col])-(1 if k in ('x','y') else 0) for row in rows]
        kron=np.array([float(row['KRON_RADIUS']) for row in rows]);flags=np.array([int(row['KRON_FLAGS']) for row in rows])
        with fits.open(d/'masks.fits.gz') as h:
            seg=h['SEGMENTATION'].data.copy();base=h['BASE_MASK'].data.astype(bool);header=h['BASE_MASK'].header.copy()
        wcs=WCS(header).celestial;f=frames[cid];center=wcs.all_world2pix(f['ra_deg'],f['dec_deg'],0)
        scale=3600*np.sqrt(abs(np.linalg.det(wcs.pixel_scale_matrix)))
        with (d/'members.tsv').open() as stream:members=np.array([[float(x['X_IMAGE'])-1,float(x['Y_IMAGE'])-1] for x in csv.DictReader(stream,delimiter='\t') if x['IN_FRAME']=='True']).reshape(-1,2)
        satellites,excluded=satellite_members(members,center,scale,r['stable_peak_candidate'])
        bcg=r['candidate_labels'][0]-1 if r['stable_peak_candidate'] else -1
        masks,nfallback=prep.source_masks(objects,seg,kron,flags,base,bcg,satellites,scale)
        assert all(np.all(mask[base]) for mask in masks)
        assert np.all(~masks[0]|masks[1]) and np.all(~masks[1]|masks[2])
        header['SCIVAL']=False;header['BCGOK']=r['stable_peak_candidate']
        hdus=[fits.PrimaryHDU(),fits.ImageHDU(base.astype('u1'),header=header,name='BASE_MASK')]
        for factor,mask in zip(prep.FACTORS,masks):hdus.append(fits.ImageHDU(mask.astype('u1'),header=header,name='MASK_K'+str(factor).replace('.','P')))
        fits.HDUList(hdus).writeto(out/'proposed_masks.fits.gz')
        ix,iy=[round(float(c)) for c in center]
        record=dict(r,bcg_member_excluded_from_fallback=excluded,remaining_member_fallbacks=nfallback,
            center_masked_before=r['center_masked_by_kron'],center_masked_by_kron=[bool(mask[iy,ix]) for mask in masks])
        prep.write_json(out/'review.json',record);records.append(record)
    prep.write_json(a.work/'summary.json',dict(status='complete',records=records,science_measurement=False,mask_ready_for_science=False,
        manifest_sha256=review['manifest_sha256'],association_summary_sha256=hashlib.sha256((a.reviewed/'summary.json').read_bytes()).hexdigest(),
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),dependency_sha256=hashlib.sha256(Path(prep.__file__).read_bytes()).hexdigest()))
    print(f'Refined {len(records)} masks; science approval remains false')


if __name__=='__main__':main()
