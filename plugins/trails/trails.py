#!/usr/bin/env python3
"""Satellite / aircraft trail removal: detect linear trails, build a mask, fill or exclude, flag catalogue objects.

    trails.py IMAGE --work DIR [--mask MASKFILE] [--catalog TSV] [--threshold 8 --min-length 0 --max-trails 8 --smooth 1]
              [--margin 2 --edge-sigma 2 --end-extend 3] [--fill mask|interpolate|stack] [--fill-noise] [--touch-k 2.5]
              [--catalog-check] [--no-flag-catalog] [--pixel-scale S --mag-zeropoint Z] [--meta-out FILE]
    trails.py IMAGE --task stack --frames F1 F2 ... --work DIR [--stack-method sigclip|median|mean] [--frame-masks M1 M2 ...]

What it writes into DIR (prefix trails_): trails.json (all trails: pixel end points, angle, length, width, amplitude, scores), mask.fits (0/1
uint8, 1 = trail), flags.tsv (NUMBER TRAIL_FLAG TRAIL_ID TRAIL_DIST, with --catalog), <base>_trailfree.fits (--fill interpolate),
<base>_trails_mef.fits (--fill stack: SCI + TRAILMASK extension, for stackers that read mask extensions), trails_overlay.reg (ds9 regions).
With --mask the trail bit (32) of the shared mask manager's flag mask (ds9/library/ds9_mask.py --mode trails) is set, so the mask can be undone
with the mask manager's undo.  stdout: '#TRAILS key=value ...' summary and, with --catalog, the NUMBER + TRAIL_* columns (add_columns contract).
The image file itself is never modified.
"""
import argparse
import json
import math
import os
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
import warnings  # noqa: E402
warnings.filterwarnings('ignore')
from ogfkit import trails as T, trails_ext as TX, tsvio, imageio, meta as ometa  # noqa: E402


FLAGCOLS = ['TRAIL_FLAG', 'TRAIL_ID', 'TRAIL_DIST', 'TRAIL_FLUX', 'TRAIL_FRAC']


def base_of(path):
    b = os.path.basename(path)
    for suf in ('.gz', '.fits', '.fit', '.fts'):
        if b.endswith(suf):
            b = b[:-len(suf)]
    return b


def write_fits(path, arr, header=None, **kw):
    from astropy.io import fits
    h = fits.PrimaryHDU(arr, header=header)
    for k, v in kw.items():
        h.header[k] = v
    h.writeto(path, overwrite=True)


def read_cat(path):
    if not path or not os.path.isfile(path):
        return None
    cols, rows = tsvio.read_catalog(path)
    if not rows or 'X_IMAGE' not in cols:
        return None
    cat = {'NUMBER': [r.get('NUMBER', str(i + 1)) for i, r in enumerate(rows)]}
    for c in ('X_IMAGE', 'Y_IMAGE', 'A_IMAGE', 'B_IMAGE', 'THETA_IMAGE', 'FLUX_AUTO', 'FLUX_ISO'):
        if c in cols:
            cat[c] = np.array([tsvio.fnum(r.get(c)) for r in rows])
    return cat


def sky_header(hdr):
    """WCS keys only (the output images keep the geometry of the input)"""
    from astropy.io import fits
    keep = fits.Header()
    for k in ('CRPIX1', 'CRPIX2', 'CRVAL1', 'CRVAL2', 'CD1_1', 'CD1_2', 'CD2_1', 'CD2_2', 'PC1_1', 'PC1_2', 'PC2_1', 'PC2_2', 'CDELT1', 'CDELT2',
              'CTYPE1', 'CTYPE2', 'CUNIT1', 'CUNIT2', 'EQUINOX', 'RADESYS', 'LONPOLE', 'LATPOLE'):
        if k in hdr:
            keep[k] = hdr[k]
    return keep


def detect(a, img, hdr, cat):
    return T.detect_trails(img, threshold=a.threshold, min_length=a.min_length or None, max_trails=a.max_trails, smooth=a.smooth,
                           margin=a.margin, edge_sigma=a.edge_sigma, min_aspect=a.min_aspect, reject_bleeds=a.reject_bleeds, spike_fac=a.spike_fac, max_fwhm=a.max_fwhm, catalog=cat if a.catalog_check else None, catalog_rescue=0.7 if a.catalog_check else None)


def run_detect(a):
    os.makedirs(a.work, exist_ok=True)
    img, hdr = imageio.load_image(a.image)
    base = base_of(a.image)
    cat = read_cat(a.catalog)
    res = detect(a, img, hdr, cat)
    trails = res['trails']
    if a.curved:                     # curved / segmented trails: tile-wise detection + chain linking (segments share a `group` id)
        for t in TX.detect_curved(img, trails, margin=a.margin, edge_sigma=a.edge_sigma, min_aspect=a.min_aspect, max_fwhm=a.max_fwhm):
            t['id'] = len(trails) + 1
            trails.append(t)
    for t in trails:                 # flicker: duty cycle along the track
        try:
            ap_ = TX.along_profile(img, t)
            t['duty'] = float(ap_['duty']); t['n_gaps'] = int(ap_['n_gaps']); t['flicker'] = bool(ap_['duty'] < 0.8 and ap_['n_gaps'] >= 2)
        except Exception:
            pass
    mask = T.trail_mask(img.shape, trails, end_extend=a.end_extend) if trails else np.zeros(img.shape, bool)
    out = lambda n: os.path.join(a.work, 'trails_' + n)
    kw = sky_header(hdr)
    write_fits(out('mask.fits'), mask.astype(np.uint8), kw, OGFTRAIL=(1, 'trail mask 1 = trail'))
    summary = dict(image=os.path.abspath(a.image), n_trails=len(trails), masked_pixels=int(mask.sum()), masked_fraction=float(mask.mean()),
                   sigma_pix=res['sigma_pix'], threshold=a.threshold, best_zscore=res.get('best_zscore'), fill=a.fill, binning=res['binning'],
                   trails=trails)
    for t in summary['trails']:
        for k in list(t):
            if isinstance(t[k], (np.floating, np.integer)):
                t[k] = float(t[k])
    if a.pixel_scale:
        for t in summary['trails']:
            t['width_arcsec'] = t['fwhm'] * a.pixel_scale
            if a.mag_zeropoint is not None and t['amp'] > 0:
                t['peak_mu_mag_arcsec2'] = a.mag_zeropoint - 2.5 * math.log10(t['amp'] / a.pixel_scale ** 2)
    msg = ''
    # ds9 regions for the overlay
    with open(out('overlay.reg'), 'w') as f:
        f.write('# Region file format: DS9 version 4.1\nglobal color=cyan width=2\nimage\n')
        for t in trails:
            f.write('line(%.2f,%.2f,%.2f,%.2f) # line=0 0 text={trail %d%s}\n' % (t['x1'], t['y1'], t['x2'], t['y2'], t['id'], (' g%d' % t['group']) if t.get('group') else ''))
    # fills
    if trails and a.fill == 'interpolate':
        fixed = T.fill_interpolate(img, trails, mask=mask, noise=a.fill_noise)
        p = os.path.join(a.work, '%s_trailfree.fits' % base)
        write_fits(p, fixed.astype(np.float32), hdr)
        summary['filled_image'] = p
    if a.fill == 'stack' or a.write_mef:
        from astropy.io import fits
        p = os.path.join(a.work, '%s_trails_mef.fits' % base)
        h0 = fits.PrimaryHDU(header=hdr)
        sci = fits.ImageHDU(img.astype(np.float32), header=hdr, name='SCI')
        mk = fits.ImageHDU(mask.astype(np.uint8), header=sky_header(hdr), name='TRAILMASK')
        mk.header['COMMENT'] = '1 = satellite/aircraft trail (exclude from stacks)'
        fits.HDUList([h0, sci, mk]).writeto(p, overwrite=True)
        summary['mef'] = p
    # catalogue flags
    flags_rows = None
    if cat is not None and a.flag_catalog:
        fl, tid, dist = T.flag_catalog(cat, trails, touch_k=a.touch_k) if trails else (np.zeros(len(cat['X_IMAGE']), int), np.zeros(len(cat['X_IMAGE']), int), np.full(len(cat['X_IMAGE']), -1.0))
        tfx, tfr = TX.trail_flux(cat, trails) if trails else (np.zeros(len(fl)), np.full(len(fl), -1.0))
        flags_rows = [(n, dict(TRAIL_FLAG=int(f), TRAIL_ID=int(i), TRAIL_DIST=float(d), TRAIL_FLUX=float(x), TRAIL_FRAC=float(q)))
                      for n, f, i, d, x, q in zip(cat['NUMBER'], fl, tid, dist, tfx, tfr)]
        with open(out('flags.tsv'), 'w') as f:
            tsvio.write_columns(FLAGCOLS, flags_rows, f)
        summary['n_flux_contaminated'] = int((tfx > 0).sum())
        summary['n_flagged'] = int((fl > 0).sum())
        summary['n_flagged_centre_inside'] = int(((fl & 2) > 0).sum())
        summary['n_trail_fragments'] = int(((fl & 4) > 0).sum())
    json.dump(summary, open(os.path.join(a.work, 'trails.json'), 'w'), indent=1, default=float)
    # shared mask manager (undoable)
    if a.mask:
        if trails:
            script = os.path.join(ROOT, 'ds9', 'library', 'ds9_mask.py')
            r = subprocess.run([sys.executable, script, a.image, '--mode', 'trails', '--mask', a.mask, '--file', out('mask.fits'), '--fresh'],
                               capture_output=True, text=True)
            if r.returncode:
                sys.stderr.write(r.stderr)
                sys.exit('ds9_mask.py failed')
            msg = [l for l in r.stdout.split('\n') if l.startswith('#MASK_STATS')]
            summary['mask_stats'] = msg[-1] if msg else ''
        else:
            summary['mask_stats'] = 'no trail: shared mask left unchanged'
    print('#TRAILS n=%d masked_fraction=%.5f n_flagged=%s fill=%s' % (len(trails), mask.mean(), summary.get('n_flagged', '-'), a.fill))
    for t in trails:
        print('#TRAIL id=%d x1=%.1f y1=%.1f x2=%.1f y2=%.1f angle=%.2f length=%.0f fwhm=%.1f halfwidth=%.1f z=%.1f' %
              (t['id'], t['x1'], t['y1'], t['x2'], t['y2'], t['theta_deg'], t['length'], t['fwhm'], t['halfwidth'], t['zscore']))
    if flags_rows is not None:
        tsvio.write_columns(FLAGCOLS, flags_rows)
    if a.meta_out:
        ometa.update(a.meta_out, 'trails', dict(n_trails=len(trails), masked_fraction=float(mask.mean()), fill=a.fill))
    return 0


def _wcs(hdr):
    try:
        from astropy.wcs import WCS
        w = WCS(hdr).celestial
        return w if w.has_celestial else None
    except Exception:
        return None


def _differs(w0, w1, shape):
    """True when the two WCS map the image centre and corners to pixels that differ by more than 0.05 px"""
    ny, nx = shape
    pts = np.array([[nx / 2, ny / 2], [0, 0], [nx - 1, ny - 1]], float)
    ra, de = w1.all_pix2world(pts[:, 0], pts[:, 1], 0)
    x, y = w0.all_world2pix(ra, de, 0)
    return float(np.max(np.hypot(x - pts[:, 0], y - pts[:, 1]))) > 0.05


def run_stack(a):
    """combine frames excluding each frame's trails; with --register (default auto) frames with a different WCS are resampled (bilinear) onto the first
    frame's grid together with their trail masks (detection runs on the native frame, before the resampling)"""
    os.makedirs(a.work, exist_ok=True)
    frames, masks, ntr = [], [], []
    hdr0, w0, shape0, nreg = None, None, None, 0
    for i, f in enumerate(a.frames):
        img, hdr = imageio.load_image(f)
        if hdr0 is None:
            hdr0, shape0 = hdr, img.shape
            w0 = _wcs(hdr) if a.register != 'none' else None
        if a.frame_masks and i < len(a.frame_masks) and a.frame_masks[i] != '-':
            m = imageio.load_mask(a.frame_masks[i], img.shape)
        else:
            r = detect(a, img, hdr, None)
            m = T.trail_mask(img.shape, r['trails'], end_extend=a.end_extend) if r['trails'] else np.zeros(img.shape, bool)
            ntr.append(len(r['trails']))
        bad = m | (img == 0)
        wi = _wcs(hdr) if (w0 is not None and i > 0) else None
        if wi is not None and (a.register == 'wcs' or img.shape != shape0 or _differs(w0, wi, img.shape)):
            im2 = TX.wcs_resample(np.where(bad & (img == 0), np.nan, img), wi, w0, shape0, order=1, cval=np.nan)
            mk, nodata = TX.register_masks(m, wi, w0, shape0, grow=1)
            img, bad = im2, mk | nodata | ~np.isfinite(im2)
            nreg += 1
        elif a.register == 'wcs' and i > 0 and wi is None:
            sys.stderr.write('frame %s has no usable WCS: not resampled\n' % f)
        frames.append(img); masks.append(bad)
    out, n = T.stack_frames(frames, masks, method=a.stack_method)
    p = os.path.join(a.work, 'trails_stack.fits')
    write_fits(p, np.nan_to_num(out).astype(np.float32), hdr0)
    write_fits(os.path.join(a.work, 'trails_stack_n.fits'), n.astype(np.int16), sky_header(hdr0))
    print('#TRAILS_STACK frames=%d method=%s trails_per_frame=%s registered=%d out=%s' % (len(frames), a.stack_method, ','.join(map(str, ntr)), nreg, p))
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('image')
    ap.add_argument('--task', choices=['detect', 'stack'], default='detect')
    ap.add_argument('--work', required=True)
    ap.add_argument('--mask', help='flag mask of the shared mask manager (~/.ds9/mask_<base>.fits): the trail bit is set there (undoable)')
    ap.add_argument('--catalog')
    ap.add_argument('--threshold', type=float, default=8.0)
    ap.add_argument('--min-length', type=float, default=0.0, help='minimum trail length in pixels (0 = 12 percent of the diagonal)')
    ap.add_argument('--max-trails', type=int, default=8)
    ap.add_argument('--smooth', type=float, default=1.0)
    ap.add_argument('--min-aspect', type=float, default=12.0, help='reject candidates shorter than this many times their FWHM (galaxy chains)')
    ap.add_argument('--max-fwhm', type=float, default=40.0)
    ap.add_argument('--curved', action='store_true', help='also look for curved / broken (segmented) trails: tile-wise detection and chain linking (slower)')
    ap.add_argument('--spike-fac', type=float, default=30.0, help='diffraction-spike rule: a run through a core brighter than this x the run amplitude whose brightness falls off along the line is rejected (300 = the pre-SDSS value: fewer spikes rejected, ~10 %% more faint HUDF trails kept)')
    ap.add_argument('--keep-bleeds', dest='reject_bleeds', action='store_false', help='do not reject column/row-aligned runs with a saturated core (CCD bleeds)')
    ap.add_argument('--margin', type=float, default=2.0)
    ap.add_argument('--edge-sigma', type=float, default=2.0)
    ap.add_argument('--end-extend', type=float, default=3.0)
    ap.add_argument('--fill', choices=['mask', 'interpolate', 'stack'], default='mask')
    ap.add_argument('--fill-noise', action='store_true')
    ap.add_argument('--write-mef', action='store_true')
    ap.add_argument('--no-flag-catalog', dest='flag_catalog', action='store_false', help='do not write the TRAIL_FLAG / TRAIL_ID / TRAIL_DIST columns')
    ap.add_argument('--touch-k', type=float, default=2.5)
    ap.add_argument('--catalog-check', action='store_true', help='use elongated aligned catalogue objects as extra evidence (lowers the threshold to 70 percent for supported candidates)')
    ap.add_argument('--pixel-scale', type=float, default=0.0)
    ap.add_argument('--mag-zeropoint', type=float, default=None)
    ap.add_argument('--meta-out')
    ap.add_argument('--frames', nargs='*')
    ap.add_argument('--frame-masks', nargs='*')
    ap.add_argument('--register', choices=['none', 'wcs', 'auto'], default='auto', help='stack: resample frames and their trail masks onto the first frame through the WCS (auto = when the WCS differ)')
    ap.add_argument('--stack-method', default='sigclip', choices=['sigclip', 'median', 'mean'])
    a = ap.parse_args(argv)
    return run_stack(a) if a.task == 'stack' else run_detect(a)


if __name__ == '__main__':
    sys.exit(main())
