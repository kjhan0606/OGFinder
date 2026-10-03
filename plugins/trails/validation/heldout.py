#!/usr/bin/env python3
"""Held-out check on real images that were not used to tune the detector (the detector code is frozen at 58f69ec69 for the straight finder):
  sdss   SDSS Stripe 82 run 94 camcol 3 frames, 5 bands x N fields (ground based, 1489x2048).  Every detection is measured in its own frame and looked for
         along the same sky line in the other four bands (WCS); present in one band only = transient candidate (satellite / aircraft / meteor),
         present in several = static feature (spike, bleed, scattered light, sky structure).  PNG cutouts of all bands are rendered.
  hudf   unseen 700x700 windows of the HUDF F160W mosaic (not overlapping the 40 tuning windows) and of the F105W / F125W mosaics (other filters, other noise).
usage: heldout.py --exp sdss,hudf --out DIR [--workers 7] [--curved]
Nothing here looks at images by eye: the numbers (length, FWHM, amplitude, band consistency, along-track uniformity) are the inspection."""
import argparse, glob, json, math, os, sys, time, warnings
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.abspath(os.path.join(HERE, '..', '..', '..')))
sys.path.insert(0, HERE)
warnings.simplefilter('ignore')
from ogfkit import trails as T, trails_ext as X
import trails_validate as V

SDSS_DIR = os.environ.get('TRAILS_SDSS', '/workspace/work/item4/sdss')
BANDS = 'ugriz'


def load_sdss(band, field):
    from astropy.io import fits
    p = os.path.join(SDSS_DIR, 'frame-%s-000094-3-%04d.fits' % (band, field))
    if not os.path.exists(p):
        return None
    with fits.open(p) as h:
        return np.asarray(h[0].data, np.float32), h[0].header


def to_line(t, wA, wB, shapeB):
    """trail end points (1-based pixels of frame A) -> theta, rho, s0, s1 in frame B (None when outside)"""
    ra, de = wA.all_pix2world([t['x1'] - 1, t['x2'] - 1], [t['y1'] - 1, t['y2'] - 1], 0)
    xb, yb = wB.all_world2pix(ra, de, 0)
    ny, nx = shapeB
    th = math.degrees(math.atan2(yb[1] - yb[0], xb[1] - xb[0]))
    if th < 0 or th >= 180:
        xb, yb = xb[::-1], yb[::-1]; th = math.degrees(math.atan2(yb[1] - yb[0], xb[1] - xb[0]))
    thr = math.radians(th); cx, cy = (nx - 1) / 2, (ny - 1) / 2
    rho = -(xb[0] - cx) * math.sin(thr) + (yb[0] - cy) * math.cos(thr)
    s0 = (xb[0] - cx) * math.cos(thr) + (yb[0] - cy) * math.sin(thr); s1 = (xb[1] - cx) * math.cos(thr) + (yb[1] - cy) * math.sin(thr)
    return th, rho, s0, s1


def clip_run(th, rho, s0, s1, shape):
    ny, nx = shape
    ch = T.chord(th, rho, nx, ny)
    if ch is None:
        return None
    a, b = max(s0, ch[0] + 1), min(s1, ch[1] - 1)
    return (a, b) if b - a > 30 else None


def measure_other(img, valid, sig, th, rho, s0, s1):
    r = clip_run(th, rho, s0, s1, img.shape)
    if r is None:
        return None
    p = T.measure_profile(img, valid, th, rho, r[0], r[1], sig)
    if p is None:
        return dict(snr=0.0, covered=r[1] - r[0])
    return dict(snr=float(p['profile_snr']), amp=float(p['amp']), fwhm=float(p['fwhm']), t0=float(p['t0']), covered=float(r[1] - r[0]))


def sdss_frame(args):
    band, field, curved, png_dir = args
    from astropy.wcs import WCS
    d = load_sdss(band, field)
    if d is None:
        return None
    img, hdr = d
    t0 = time.time()
    res = T.detect_trails(img.copy())
    trails = list(res['trails'])
    cur = X.detect_curved(img, trails) if curved else []
    out = dict(band=band, field=field, sec=time.time() - t0, best_z=res.get('best_zscore'), n_straight=len(trails), n_curved=len(cur), det=[],
               rejected=[(r['reason'], round(r['zscore'], 1)) for r in res.get('rejected', [])])
    wA = WCS(hdr).celestial
    for t in trails + cur:
        row = dict(id=t['id'], curved=bool(t.get('group')), theta=t['theta_deg'], length=t['length'], fwhm=t['fwhm'], amp=t['amp'], amp_snr=t['amp_snr'], z=t['zscore'],
                   x1=t['x1'], y1=t['y1'], x2=t['x2'], y2=t['y2'])
        ap = X.along_profile(img, t)
        row.update(duty=ap['duty'], cv=ap['cv'])
        others = {}
        for b in BANDS:
            if b == band:
                continue
            e = load_sdss(b, field)
            if e is None:
                continue
            imB, hB = e
            wB = WCS(hB).celestial
            th, rho, s0, s1 = to_line(t, wA, wB, imB.shape)
            sigB = T.standardise(imB, np.isfinite(imB) & (imB != 0))[1]
            m = measure_other(imB, np.isfinite(imB) & (imB != 0), sigB, th, rho, s0, s1)
            others[b] = m
        row['others'] = others
        cov = {b: m for b, m in others.items() if m is not None and m.get('covered', 0) > 0.3 * t['length']}
        row['n_cov'] = len(cov); row['n_present'] = sum(1 for m in cov.values() if m['snr'] >= 6.0)
        row['class'] = 'no-coverage' if not cov else ('transient-candidate (one band)' if row['n_present'] == 0 else 'static (%d other bands)' % row['n_present'])
        row['tai'] = hdr.get('TAI')
        out['det'].append(row)
        if png_dir:
            try:
                render(png_dir, band, field, t, row, img, hdr)
            except Exception as ex:
                row['png_error'] = str(ex)[:80]
    return out


def render(png_dir, band, field, t, row, img, hdr):
    import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
    from astropy.wcs import WCS
    wA = WCS(hdr).celestial
    mx, my = (t['x1'] + t['x2']) / 2 - 1, (t['y1'] + t['y2']) / 2 - 1
    fig, axs = plt.subplots(1, 5, figsize=(20, 4.2))
    for ax, b in zip(axs, BANDS):
        e = load_sdss(b, field)
        if e is None:
            ax.axis('off'); continue
        im, h = e
        wB = WCS(h).celestial
        ra, de = wA.all_pix2world([mx, t['x1'] - 1, t['x2'] - 1], [my, t['y1'] - 1, t['y2'] - 1], 0)
        xb, yb = wB.all_world2pix(ra, de, 0)
        r = 150
        x0, y0 = int(xb[0]) - r, int(yb[0]) - r
        sub = im[max(y0, 0):y0 + 2 * r, max(x0, 0):x0 + 2 * r]
        med = np.median(sub); s = 1.4826 * np.median(abs(sub - med)) + 1e-6
        ax.imshow(sub, origin='lower', cmap='gray', vmin=med - 2 * s, vmax=med + 6 * s, extent=(max(x0, 0), max(x0, 0) + sub.shape[1], max(y0, 0), max(y0, 0) + sub.shape[0]))
        ax.plot(xb[1:], yb[1:], 'r-', lw=.6, alpha=.4)
        ax.set_xlim(x0, x0 + 2 * r); ax.set_ylim(y0, y0 + 2 * r)
        ax.set_title('%s%s %s' % (b, ' (det)' if b == band else '', ('snr %.1f' % row['others'][b]['snr']) if b != band and row['others'].get(b) else ''), fontsize=9)
    fig.suptitle('%s field %d trail %d: theta %.1f len %.0f fwhm %.1f z %.1f -> %s' % (band, field, t['id'], t['theta_deg'], t['length'], t['fwhm'], t['zscore'], row['class']), fontsize=9)
    os.makedirs(png_dir, exist_ok=True)
    fig.savefig(os.path.join(png_dir, '%s_%04d_%d.png' % (band, field, t['id'])), dpi=45)
    plt.close(fig)


def hudf_task(args):
    path, win, curved, tag = args
    from astropy.io import fits
    d = V.fits_data(path)
    y, x = win
    img = d[y:y + 700, x:x + 700].copy()
    t0 = time.time()
    res = T.detect_trails(img)
    cur = X.detect_curved(img, res['trails']) if curved else []
    return dict(tag=tag, win=win, sec=time.time() - t0, best_z=res.get('best_zscore'), n_straight=len(res['trails']), n_curved=len(cur),
                det=[dict(theta=t['theta_deg'], length=t['length'], fwhm=t['fwhm'], z=t['zscore'], curved=bool(t.get('group')), x1=t['x1'], y1=t['y1'], x2=t['x2'], y2=t['y2']) for t in res['trails'] + cur],
                rejected=[(r['reason'], round(r['zscore'], 1)) for r in res.get('rejected', [])])


def unseen_windows(path, n, seed, exclude_tuning=False):
    """Non-overlapping fully covered 700x700 windows of `path`.  The 40 F160W tuning windows cover essentially all of the covered F160W
    area, so no F160W window is unseen; F105W / F125W frames were never used for tuning (same sky, independent pixels and noise)."""
    tuning = V.hudf_windows(700, 40) if exclude_tuning else []
    d = V.fits_data(path)
    v = (d != 0).astype(np.float32)
    ii = np.cumsum(np.cumsum(np.pad(v, ((1, 0), (1, 0))), 0), 1)
    pos = []
    for y in range(0, d.shape[0] - 700, 100):
        for x in range(0, d.shape[1] - 700, 100):
            if (ii[y + 700, x + 700] - ii[y, x + 700] - ii[y + 700, x] + ii[y, x]) > 0.9995 * 700 * 700 and \
                    all(abs(y - ty) >= 700 or abs(x - tx) >= 700 for ty, tx in tuning):
                pos.append((y, x))
    rng = np.random.default_rng(seed)
    out = []
    for i in rng.permutation(len(pos)):
        p = pos[i]
        if all(abs(p[0] - q[0]) >= 700 or abs(p[1] - q[1]) >= 700 for q in out):
            out.append(p)
        if len(out) >= n:
            break
    return out, len(pos)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--exp', default='sdss,hudf')
    ap.add_argument('--out', required=True)
    ap.add_argument('--workers', type=int, default=7)
    ap.add_argument('--curved', action='store_true')
    ap.add_argument('--fields', type=int, default=0, help='limit the number of SDSS fields (0 = all)')
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    from multiprocessing import Pool
    md, allres = ['# held-out false-positive / real-trail check (heldout.py%s)' % (' --curved' if a.curved else '')], {}
    if 'sdss' in a.exp:
        fields = sorted(set(int(os.path.basename(f).split('-')[-1][:4]) for f in glob.glob(os.path.join(SDSS_DIR, 'frame-r-*.fits'))))
        if a.fields:
            fields = fields[:a.fields]
        tasks = [(b, f, a.curved, os.path.join(a.out, 'png')) for f in fields for b in BANDS]
        t0 = time.time()
        with Pool(a.workers) as p:
            R = [r for r in p.map(sdss_frame, tasks, chunksize=1) if r]
        allres['sdss'] = R
        dets = [dict(d, band=r['band'], field=r['field']) for r in R for d in r['det']]
        md.append('\n## SDSS run 94 camcol 3: %d frames (%d fields x 5 bands), %.0f s\n' % (len(R), len(fields), time.time() - t0))
        md.append('frames with >= 1 detection: %d; detections: %d (straight %d, curved-chain segments %d); best z of detection-free frames: median %.1f, 95%% %.1f, max %.1f' % (
            sum(1 for r in R if r['det']), len(dets), sum(1 for d in dets if not d['curved']), sum(1 for d in dets if d['curved']),
            np.median([r['best_z'] for r in R if not r['det'] and r['best_z']]), np.percentile([r['best_z'] for r in R if not r['det'] and r['best_z']], 95), max([r['best_z'] for r in R if not r['det'] and r['best_z']] or [0])))
        cl = {}
        for d in dets:
            cl.setdefault(d['class'].split(' (')[0], []).append(d)
        for k, v in cl.items():
            md.append('* %s: %d' % (k, len(v)))
        md.append('\n| band | field | id | theta | length | FWHM | amp [sigma_pix] | z | duty | other bands present / covered (snr) | class |')
        md.append('|---|---|---|---|---|---|---|---|---|---|---|')
        for d in sorted(dets, key=lambda d: (d['field'], d['band'])):
            md.append('| %s | %d | %d%s | %.1f | %.0f | %.1f | %.1f | %.1f | %.2f | %d/%d (%s) | %s |' % (d['band'], d['field'], d['id'], 'c' if d['curved'] else '', d['theta'], d['length'], d['fwhm'], d['amp_snr'], d['z'], d['duty'],
                      d['n_present'], d['n_cov'], ' '.join('%s%.0f' % (b, m['snr']) for b, m in d['others'].items() if m), d['class']))
    if 'hudf' in a.exp:
        tasks, npools = [], {}
        for f in ('f105w', 'f125w'):
            path = '/workspace/fits/hudf_%s.fits' % f
            wins, npools[f] = unseen_windows(path, 20, 5)
            tasks += [(path, w, a.curved, f.upper()) for w in wins]
        t0 = time.time()
        with Pool(a.workers) as p:
            R = p.map(hudf_task, tasks, chunksize=1)
        allres['hudf'] = R
        md.append('\n## HUDF held-out windows (non-overlapping 700x700 windows of the F105W / F125W frames, never used for tuning; candidate positions %s), %.0f s\n' % (npools, time.time() - t0))
        md.append('| set | windows | with a detection | detections (straight / curved segments) | best z median / max |')
        md.append('|---|---|---|---|---|')
        for tag in ('F105W', 'F125W'):
            g = [r for r in R if r['tag'] == tag]
            bz = [r['best_z'] for r in g if r['best_z']]
            md.append('| %s | %d | %d | %d / %d | %.1f / %.1f |' % (tag, len(g), sum(1 for r in g if r['det']), sum(r['n_straight'] for r in g), sum(r['n_curved'] for r in g), np.median(bz or [0]), max(bz or [0])))
        for r in R:
            for d in r['det']:
                md.append('* %s window %s: theta %.1f length %.0f FWHM %.1f z %.1f%s' % (r['tag'], r['win'], d['theta'], d['length'], d['fwhm'], d['z'], ' (curved chain)' if d['curved'] else ''))
    json.dump(allres, open(os.path.join(a.out, 'heldout.json'), 'w'), default=float)
    open(os.path.join(a.out, 'heldout.md'), 'w').write('\n'.join(md) + '\n')
    print('\n'.join(md))


if __name__ == '__main__':
    main()
