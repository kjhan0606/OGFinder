#!/usr/bin/env python3
"""Retrieve attributed public i-band reference PSF data, not science calibration."""
import argparse
import hashlib
import json
from pathlib import Path
import platform
import urllib.request

import numpy as np
from astropy.io import fits

from fetch_ancillary_catalogues import sha
import prepare_detection_psf as prep

RECORD='https://zenodo.org/api/records/14915578'
EXPECTED={'HSC_PSF_i.fits':'5bf3317bcc41a83c81c25acab2faab09',
          'HSC_PSF_i_asymmetric.fits':'665e7b459cda203d6610133c3626194d',
          'README.md':'feb630894b1183bcb77f2932d71fb9ad'}
SCALE=.168  # Publication Section 2; not inferred from a missing FITS WCS.


def radial_summary(data):
    if data.ndim!=2 or data.shape[0]!=data.shape[1] or data.shape[0]%2!=1:
        raise ValueError('Reference PSF must have an odd square native grid')
    if not np.isfinite(data).all():raise ValueError('Nonfinite reference PSF')
    center=(data.shape[0]-1)//2
    if np.unravel_index(np.argmax(data),data.shape)!=(center,center):
        raise ValueError('Reference PSF peak is not at the documented grid center')
    total=float(data.sum(dtype=float))
    if not np.isfinite(total) or abs(total-1)>1e-3:raise ValueError('Reference finite-grid normalization differs from one')
    yy,xx=np.ogrid[:data.shape[0],:data.shape[1]];r=np.hypot(xx-center,yy-center)
    bins=np.floor(r).astype('i4');n=np.bincount(bins.ravel());s=np.bincount(bins.ravel(),weights=data.ravel())
    radii=np.bincount(bins.ravel(),weights=r.ravel())/n
    rows=[dict(radius_bin_pix=i,mean_radius_pix=float(radii[i]),mean=float(s[i]/n[i]),n_pixels=int(n[i]),
               sum=float(s[i]),cumulative_sum=float(s[:i+1].sum()),full_annulus=i<center)
          for i in range(len(n)) if n[i]]
    info=dict(shape=list(data.shape),finite_grid_sum=total,negative_pixels=int(np.count_nonzero(data<0)),
        center_pixel_xy=[center,center],pixel_scale_arcsec=SCALE,scale_source='Garate-Nunez et al. 2024 Section 2',
        inscribed_radius_arcsec=center*SCALE,encircled_grid_fractions={str(rad):float(data[r<rad/SCALE].sum()/total)
           for rad in (2.016,6.,16.128,30.,60.,120.,300.)},
        normalization_scope='finite published square grid; flux outside it is not measured',psf_approved=False)
    return rows,info


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--work',type=Path,required=True);a=p.parse_args()
    if not platform.node().lower().startswith('lageunha'):p.error('Run on lageunha CPU')
    a.work.mkdir(parents=True,exist_ok=False)
    with urllib.request.urlopen(RECORD,timeout=60) as response:raw=response.read(4*1024**2)
    metadata=json.loads(raw);(a.work/'zenodo_record.json').write_bytes(raw)
    if metadata['id']!=14915578 or metadata['metadata']['license']['id']!='cc-by-4.0':raise ValueError('Dataset metadata changed')
    files={f['key']:f for f in metadata['files']};records=[]
    for name,digest in EXPECTED.items():
        item=files[name];url=item['links']['self']
        if item['checksum']!='md5:'+digest or item['size']>64*1024**2:raise ValueError('Pinned dataset checksum/size changed')
        if not url.startswith(RECORD+'/files/'):raise ValueError('Unexpected dataset URL')
        path=a.work/name;md5=hashlib.md5();count=0
        with urllib.request.urlopen(url,timeout=90) as response,path.open('xb') as f:
            for block in iter(lambda:response.read(1024**2),b''):
                count+=len(block)
                if count>item['size']:raise ValueError('Download exceeds published size')
                f.write(block);md5.update(block)
        if count!=item['size'] or md5.hexdigest()!=digest:raise ValueError('Downloaded dataset checksum mismatch')
        entry=dict(name=name,path=str(path),url=url,bytes=count,published_md5=digest,sha256=sha(path))
        if name.endswith('.fits'):
            with fits.open(path,memmap=False) as h:
                planes=[(i,v) for i,v in enumerate(h) if v.header.get('NAXIS')==2]
                if len(planes)!=1:raise ValueError('Expected exactly one reference image plane')
                index,plane=planes[0];data=np.asarray(plane.data,dtype=float)
                if data.shape!=(4001,4001):raise ValueError('Published image dimensions changed')
                rows,info=radial_summary(data);entry.update(hdu=index,**info)
                (a.work/(name+'.header.txt')).write_text(plane.header.tostring(sep='\n',endcard=True,padding=False)+'\n')
                prep.table(a.work/(name+'.radial.tsv'),rows,list(rows[0]))
        records.append(entry);print(name,count,'bytes verified',flush=True)
    prep.write_json(a.work/'summary.json',dict(status='complete',dataset_doi=metadata['doi'],record_url='https://zenodo.org/records/14915578',
        paper='https://arxiv.org/abs/2309.16244',authors='Lucia P. Garate-Nunez et al. (2024)',license='CC-BY-4.0',
        attribution='Garate-Nunez et al. 2024, The Hyper Suprime-Cam extended Point Spread Functions and applications; doi:10.5281/zenodo.14915578',
        records=records,metadata_sha256=sha(a.work/'zenodo_record.json'),script_sha256=sha(Path(__file__)),
        science_measurement=False,psf_approved=False,
        limitations=['Published survey reference; no per-position/epoch/seeing match to the cluster frames.',
                     'Raw i-band models retain their finite-grid normalization and signed pixels.',
                     'Dataset stored outside the application repository; external R scripts are not imported or executed.']))


if __name__=='__main__':main()
