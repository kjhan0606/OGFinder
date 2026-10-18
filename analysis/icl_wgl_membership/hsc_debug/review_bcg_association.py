#!/usr/bin/env python3
"""Offline peak/segment BCG association sensitivity review; no science approval."""
import argparse
import hashlib
import json
from pathlib import Path
import platform

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
import prepare_detection_psf as prep


def peak_candidate(objects, segmentation, bad, center, scale):
    """Require a unique raw AND filtered peak and agreement with central label."""
    cx, cy = center
    ix, iy = round(float(cx)), round(float(cy))
    if not (0 <= ix < bad.shape[1] and 0 <= iy < bad.shape[0]):
        return -1, 'center_outside'
    if bad[iy, ix]:
        return -1, 'center_invalid'
    raw = np.hypot(objects['xpeak']-cx, objects['ypeak']-cy)*scale
    filtered = np.hypot(objects['xcpeak']-cx, objects['ycpeak']-cy)*scale
    candidates = np.flatnonzero((raw <= 1.) & (filtered <= 1.))
    if len(candidates) != 1:
        return -1, 'no_peak' if not len(candidates) else 'multiple_peaks'
    i = int(candidates[0])
    if segmentation[iy, ix] != i+1:
        return -1, 'central_segment_disagrees'
    return i, 'peak_segment_candidate'


def stable_candidate(first, second, pixel_scale, max_shift_arcsec=.336):
    # Two HSC pixels is an engineering stability tolerance, not a centering error.
    return bool(np.hypot(first['xpeak']-second['xpeak'],
                         first['ypeak']-second['ypeak'])*pixel_scale <= max_shift_arcsec)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--manifest', type=Path, required=True)
    parser.add_argument('--prepared', type=Path, required=True)
    parser.add_argument('--work', type=Path, required=True)
    args = parser.parse_args()
    if not platform.node().lower().startswith('lageunha'):
        parser.error('Run on lageunha CPU')
    frames = json.loads(args.manifest.read_text())['clusters']
    previous = json.loads((args.prepared/'summary.json').read_text())
    if hashlib.sha256(args.manifest.read_bytes()).hexdigest() != previous['manifest_sha256']:
        raise ValueError('Manifest differs from previous preparation')
    old = {r['id']: r for r in previous['records']}
    args.work.mkdir(parents=True, exist_ok=False)
    records = []; errors = []; star_requests = []
    import csv
    for n, frame in enumerate(frames, 1):
        cid = str(frame['id']); out = args.work/('cluster_'+cid)
        try:
            info = prep.validate_fits(Path(frame['image']), frame)
            with fits.open(frame['image'], memmap=False) as hdus:
                image = np.array(hdus[info['science_hdu']].data, dtype='f4')
                var = np.array(hdus[info['variance_hdu']].data, dtype='f4')
                bits = np.array(hdus[info['mask_hdu']].data, dtype='i4')
                header = hdus[info['science_hdu']].header.copy()
            wcs = WCS(header).celestial; center = info['center_pixel_xy']; scale = info['pixel_scale_arcsec']
            bad = ((bits & sum(1<<b for b in prep.BAD_BIT_POSITIONS)) != 0)
            bad |= ~np.isfinite(image) | ~np.isfinite(var) | (var <= 0)
            base = bad | ((bits & prep.BRIGHT_OBJECT_BIT) != 0)
            sky = np.asarray(frame['members_radec_deg']).reshape(-1,2)
            members = wcs.all_world2pix(sky,0)
            inside = np.isfinite(members).all(axis=1) & (members[:,0]>=0) & (members[:,1]>=0) & (members[:,0]<image.shape[1]) & (members[:,1]<image.shape[0])
            members = members[inside]
            runs = []; candidates = []; reasons = []
            for contrast in (.005, .001):
                run = prep.detect(image,var,bad,center,members,scale,deblend_cont=contrast)
                candidate, reason = peak_candidate(run[0],run[1],bad,center,scale)
                runs.append(run); candidates.append(candidate); reasons.append(reason)
            stable = all(i>=0 for i in candidates) and stable_candidate(runs[0][0][candidates[0]],runs[1][0][candidates[1]],scale)
            bcg = candidates[0] if stable else -1
            objects,seg,kron,kflags,_,_,extra = runs[0]
            # Reject any accidental difference from the recorded reference run.
            assert len(objects)==old[cid]['n_detections']
            assert extra['bcg_matches']==old[cid]['bcg_matches']
            masks,_ = prep.source_masks(objects,seg,kron,kflags,base,bcg,members,scale)
            assert all(np.all(mask[base]) for mask in masks)
            assert np.all(~masks[0]|masks[1]) and np.all(~masks[1]|masks[2])
            out.mkdir()
            wh = wcs.to_header(); wh['SCIVAL']=False; wh['BCGOK']=bool(stable)
            wh['COMMENT']='Experimental peak/segment association; requires review.'
            hdus=[fits.PrimaryHDU(),fits.ImageHDU(base.astype('u1'),header=wh,name='BASE_MASK')]
            for factor,mask in zip(prep.FACTORS,masks):
                hdus.append(fits.ImageHDU(mask.astype('u1'),header=wh,name='MASK_K'+str(factor).replace('.','P')))
            fits.HDUList(hdus).writeto(out/'proposed_masks.fits.gz')
            ix,iy=map(lambda x:round(float(x)),center)
            r=dict(id=cid,cohort=frame['preparation_cohort'],old_centroid_matches=old[cid]['bcg_matches'],
                   peak_segment_reasons=reasons,stable_peak_candidate=bool(stable),
                   candidate_labels=[i+1 if i>=0 else None for i in candidates],
                   bcg_base_masked=bool(base[iy,ix]),
                   center_masked_by_kron=[bool(mask[iy,ix]) for mask in masks],
                   centroid_offset_arcsec=float(np.hypot(objects[bcg]['x']-center[0],objects[bcg]['y']-center[1])*scale) if stable else None,
                   mask_ready_for_science=False,science_measurement=False)
            # A small comparison makes an association inspectable without unmasking defects.
            from matplotlib import pyplot as plt
            x0,x1=max(0,ix-100),min(image.shape[1],ix+101); y0,y1=max(0,iy-100),min(image.shape[0],iy+101)
            cut=image[y0:y1,x0:x1]; med=np.nanmedian(cut); spread=max(1e-6,1.4826*np.nanmedian(abs(cut-med)))
            fig,axes=plt.subplots(1,3,figsize=(12,4))
            axes[0].imshow(np.arcsinh((cut-med)/spread),origin='lower',cmap='gray',vmin=-1,vmax=6)
            for run,i,color in zip(runs,candidates,['cyan','lime']):
                if i>=0:
                    obj=run[0][i]; axes[0].plot(obj['xpeak']-x0,obj['ypeak']-y0,'x',color=color)
                    axes[0].plot(obj['x']-x0,obj['y']-y0,'o',fillstyle='none',color=color)
            for ax,mask,label in zip(axes[1:],[base,masks[0]],['Preserved quality mask','Proposed Kron 2.5']):
                ax.imshow(mask[y0:y1,x0:x1],origin='lower',vmin=0,vmax=1,cmap='gray');ax.set_title(label)
            axes[0].set_title('Peaks (x), centroids (o)')
            for ax in axes:ax.plot(center[0]-x0,center[1]-y0,'r+',ms=12)
            fig.suptitle(f'{cid}: stable peak candidate={stable}; engineering only')
            fig.tight_layout();fig.savefig(out/'association_review.png',dpi=110);plt.close(fig)
            prep.write_json(out/'review.json',r);records.append(r)
            with (args.prepared/('cluster_'+cid)/'psf_candidates.tsv').open() as f:
                stars=list(csv.DictReader(f,delimiter='\t'))
            for star in stars:
                ra,dec=wcs.all_pix2world(float(star['X_IMAGE'])-1,float(star['Y_IMAGE'])-1,0)
                star_requests.append(dict(cluster_id=cid,number=int(star['NUMBER']),ra_deg=float(ra),dec_deg=float(dec),
                    classification='unconfirmed',required_evidence=['external point-source classification','quality flags','positional uniqueness','isolation and wide-wing coverage']))
            print(f'[{n}/{len(frames)}] {cid}: stable={stable}, base_masked={r["bcg_base_masked"]}',flush=True)
        except Exception as exc:
            errors.append(dict(id=cid,error=str(exc))); print(f'{cid}: ERROR {exc}',flush=True)
    prep.write_json(args.work/'summary.json',dict(status='complete' if not errors else 'partial',records=records,errors=errors,
        science_measurement=False,mask_ready_for_science=False,deblend_contrasts=[.005,.001],peak_match_arcsec=1.,stability_arcsec=.336,
        manifest_sha256=hashlib.sha256(args.manifest.read_bytes()).hexdigest(),
        script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        dependency_sha256=hashlib.sha256(Path(prep.__file__).read_bytes()).hexdigest()))
    prep.write_json(args.work/'psf_external_match_requests.json',dict(candidates=star_requests,external_matches_performed=False,
        note='Existing extinction CSV has no stellar classification; random CSV is not a star catalogue. No credentials or network used.'))
    return int(bool(errors))


if __name__=='__main__':raise SystemExit(main())
