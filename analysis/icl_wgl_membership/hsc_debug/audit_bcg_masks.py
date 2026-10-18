#!/usr/bin/env python3
"""Independent saved-artifact checks for experimental BCG masks."""
import argparse
import csv
import json
from pathlib import Path
import platform
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('prepared','reviewed','refined'):p.add_argument('--'+key,type=Path,required=True)
    a=p.parse_args()
    if not platform.node().lower().startswith('lageunha'):p.error('Run on lageunha CPU')
    summary=json.loads((a.refined/'summary.json').read_text());assert summary['status']=='complete'
    original=json.loads((a.prepared/'summary.json').read_text())
    assert {r['id'] for r in original['records']}=={r['id'] for r in summary['records']}
    results=[]
    for r in summary['records']:
        cid=r['id'];sub='cluster_'+cid
        with fits.open(a.prepared/sub/'masks.fits.gz') as h:
            original_base=h['BASE_MASK'].data.copy();original_header=h['BASE_MASK'].header.copy()
        with fits.open(a.reviewed/sub/'proposed_masks.fits.gz') as h:
            before=[x.data.astype(bool) for x in h[2:]]
        with fits.open(a.refined/sub/'proposed_masks.fits.gz') as h:
            assert [x.name for x in h][1:]==['BASE_MASK','MASK_K2P5','MASK_K3P5','MASK_K5P0']
            np.testing.assert_array_equal(h['BASE_MASK'].data,original_base)
            after=[x.data.astype(bool) for x in h[2:]]
            assert all(x.shape==original_base.shape for x in after)
            base=original_base.astype(bool)
            assert all(np.all(x[base]) for x in after)
            assert np.all(~after[0]|after[1]) and np.all(~after[1]|after[2])
            assert all(np.all(~new|old) for new,old in zip(after,before))
            if not r['bcg_member_excluded_from_fallback']:
                for old,new in zip(before,after):np.testing.assert_array_equal(old,new)
            assert not h[1].header['SCIVAL'] and h[1].header['BCGOK']==r['stable_peak_candidate']
            ny,nx=base.shape;xy=np.array([[0.,0.],[nx-1.,ny-1.],[nx/2.,ny/2.]])
            sky=WCS(original_header).celestial.all_pix2world(xy,0)
            assert np.max(abs(WCS(h[1].header).celestial.all_world2pix(sky,0)-xy))<1e-5
        results.append(dict(id=cid,passed=True,removed_fallback_pixels=[int(np.count_nonzero(old&~new)) for old,new in zip(before,after)]))
    requests=json.loads((a.reviewed/'psf_external_match_requests.json').read_text())['candidates']
    assert len(requests)==sum(r['n_psf_candidates'] for r in original['records'])
    assert len({(r['cluster_id'],r['number']) for r in requests})==len(requests)
    lookup={}
    for r in original['records']:
        with (a.prepared/('cluster_'+r['id'])/'detections.tsv').open() as f:
            lookup[r['id']]={int(x['NUMBER']):x for x in csv.DictReader(f,delimiter='\t')}
    for req in requests:
        src=lookup[req['cluster_id']][req['number']]
        assert abs(float(src['RA'])-req['ra_deg'])<1e-8 and abs(float(src['DEC'])-req['dec_deg'])<1e-8
        assert req['classification']=='unconfirmed'
    report=dict(all_passed=True,n_masks=len(results),n_psf_coordinates=len(requests),records=results,
        science_measurement=False,checks=['unchanged original quality mask','HDU/shape/WCS','monotonic mask growth',
        'no newly masked pixels in member refinement','unchanged masks without unique BCG member','BCG status',
        'PSF request uniqueness and original detection sky coordinates'])
    with (a.refined/'artifact_audit.json').open('x') as f:json.dump(report,f,indent=2)
    print(f'Audited {len(results)} masks and {len(requests)} PSF coordinates')
    print('Stable BCG candidates',sum(r['stable_peak_candidate'] for r in summary['records']))
    print('BCG member excluded',sum(r['bcg_member_excluded_from_fallback'] for r in summary['records']))
    print('Center masked by Kron',[sum(r['center_masked_by_kron'][i] for r in summary['records']) for i in range(3)])


if __name__=='__main__':main()
