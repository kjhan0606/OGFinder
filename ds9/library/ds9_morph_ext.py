#!/usr/bin/env python3
"""Extended non-parametric morphology: Petrosian radii (3 eta), Petrosian flux / R20..R90 / concentration, Kron radius and magnitude, smoothness, Petrosian-segmented Gini and M20,
growth curves.  Complements ds9_morphometry.py (unchanged).

    ds9_morph_ext.py IMAGE --catalog TSV [--mask FITS] [--work DIR] [options]       stdout: NUMBER + MX_* columns (add_columns contract), logs on stderr

Catalog columns used: X_IMAGE Y_IMAGE (1-based), A_IMAGE B_IMAGE THETA_IMAGE (aperture shape), KRON_RADIUS, FLUX_RADIUS (optional).
Other catalog objects are masked inside the aperture of the target (--mask-neighbours); the image is background subtracted (sep.Background).
Files in WORK: morph_ext_growth.tsv (growth curves, NUMBER r flux eta), morph_ext_curves.png (growth curve + eta plot of the brightest objects).
"""
import argparse
import math
import multiprocessing as mp
import os
import sys

import numpy as np

_script_dir = os.path.dirname(os.path.abspath(__file__))
_project_root = os.path.abspath(os.path.join(_script_dir, '..', '..'))
if _project_root not in sys.path:
    sys.path.insert(0, _project_root)
import warnings  # noqa: E402
warnings.filterwarnings('ignore')

from ogfkit import tsvio, imageio  # noqa: E402
from morphometry import extended as X  # noqa: E402

COLUMNS = ['MX_RP', 'MX_RP_LO', 'MX_RP_HI', 'MX_FLUX_P', 'MX_MAG_P', 'MX_R20', 'MX_R50', 'MX_R80', 'MX_R90', 'MX_CONC', 'MX_KRON_R', 'MX_KRON_MAG', 'MX_SMOOTH', 'MX_GINI_P', 'MX_M20_P', 'MX_FLAG']
G = {}


def num(r, k, default=float('nan')):
    return tsvio.fnum(r.get(k), default)


def shape_of(r, circular):
    A, B = num(r, 'A_IMAGE', 2.0), num(r, 'B_IMAGE', 2.0)
    if not (A > 0):
        A = 2.0
    if not (B > 0):
        B = A
    th = num(r, 'THETA_IMAGE', 0.0)
    q = 1.0 if circular else min(max(B / A, 0.15), 1.0)
    return A, q, math.radians(th if np.isfinite(th) else 0.0)


def ell_mask(shape, x, y, a, q, theta, box=None):
    ny, nx = shape
    yy, xx = np.mgrid[:ny, :nx]
    dx, dy = xx - x, yy - y
    u = dx * math.cos(theta) + dy * math.sin(theta)
    v = -dx * math.sin(theta) + dy * math.cos(theta)
    return (u / max(a, 0.5)) ** 2 + (v / max(a * q, 0.5)) ** 2 <= 1.0


def one(i):
    a = G['args']
    r = G['rows'][i]
    data = G['data']
    zp = a.mag_zeropoint
    x, y = num(r, 'X_IMAGE') - 1, num(r, 'Y_IMAGE') - 1
    out = dict(NUMBER=r['NUMBER'], MX_FLAG=0)
    if not (np.isfinite(x) and np.isfinite(y)):
        out['MX_FLAG'] = 1024
        return out
    A, q, th = shape_of(r, a.circular)
    kr = num(r, 'KRON_RADIUS', 3.5)
    if not (kr > 0):
        kr = 3.5
    rmax = min(max(a.rmax_scale * A * kr / 3.5, 12.0), a.rmax)
    # cutout around the object (square, +-1.3 rmax), mask of the neighbours inside it
    h = int(math.ceil(1.3 * rmax)) + 2
    ny, nx = data.shape
    x0, x1, y0, y1 = max(0, int(x) - h), min(nx, int(x) + h + 1), max(0, int(y) - h), min(ny, int(y) + h + 1)
    cut = data[y0:y1, x0:x1]
    m = G['mask'][y0:y1, x0:x1].copy() if G['mask'] is not None else np.zeros(cut.shape, bool)
    m |= ~np.isfinite(cut)
    nn = 0
    if a.mask_neighbours and G['tree'] is not None:
        for j in G['tree'].query_ball_point([x, y], h * 1.5 + 60):
            if j == i:
                continue
            t = G['geo'][j]
            if not (x0 - t['A'] * a.neighbour_radius <= t['x'] <= x1 + t['A'] * a.neighbour_radius and y0 - t['A'] * a.neighbour_radius <= t['y'] <= y1 + t['A'] * a.neighbour_radius):
                continue
            # a neighbour fainter than 3 % of the target is ignored (its light is part of the noise); the target itself is never masked
            if t['flux'] < a.neighbour_min_frac * G['geo'][i]['flux']:
                continue
            mm = ell_mask(cut.shape, t['x'] - x0, t['y'] - y0, a.neighbour_radius * t['A'], t['q'], t['th'])
            m |= mm
            nn += 1
        tgt = ell_mask(cut.shape, x - x0, y - y0, 1.5, 1.0, 0.0)
        m &= ~tgt
    res = X.measure_object(cut, x - x0, y - y0, q, th, G['rms'], rmax=rmax, mask=m, etas=(a.eta_lo, a.eta, a.eta_hi), main_eta=a.eta, kron_scale=a.kron_scale,
                           kron_min=a.kron_min, smooth_frac=a.smooth_frac, seed=int(r['NUMBER']) if str(r['NUMBER']).isdigit() else i, curve_fracs=(0.2, 0.5, 0.8, 0.9), want_curve=a.want_curve)
    fp = res['flux_p']
    out.update(MX_RP=res['rp'][a.eta], MX_RP_LO=res['rp'][a.eta_lo], MX_RP_HI=res['rp'][a.eta_hi], MX_FLUX_P=fp, MX_MAG_P=(zp - 2.5 * math.log10(fp)) if fp and fp > 0 else None,
               MX_R20=res['r_frac'][0.2], MX_R50=res['r_frac'][0.5], MX_R80=res['r_frac'][0.8], MX_R90=res['r_frac'][0.9], MX_CONC=res['conc'], MX_KRON_R=res['kron_r'],
               MX_KRON_MAG=(zp - 2.5 * math.log10(res['kron_flux'])) if res['kron_flux'] and res['kron_flux'] > 0 else None, MX_SMOOTH=res['smooth'], MX_GINI_P=res['gini_p'],
               MX_M20_P=res['m20_p'], MX_FLAG=res['flag'] | (8 if nn else 0))
    if a.want_curve and 'curve' in res:
        out['curve'] = res['curve']
    return out


def plot_curves(results, path, n):
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except Exception:
        return
    res = [r for r in results if 'curve' in r][:n]
    if not res:
        return
    fig, ax = plt.subplots(len(res), 2, figsize=(7.5, 2.3 * len(res)), squeeze=False)
    for k, r in enumerate(res):
        rg, F, (re_, eta) = r['curve']
        rg, F = np.array(rg), np.array(F)
        ax[k][0].semilogx(rg, F / max(np.nanmax(F), 1e-30), 'k-')
        ax[k][1].semilogx(re_, eta, 'b-')
        ax[k][1].axhline(0.2, color='r', lw=0.6)
        for kk in (0, 1):
            if r.get('MX_RP') and np.isfinite(r['MX_RP']):
                ax[k][kk].axvline(r['MX_RP'], color='r', lw=0.6)
        ax[k][0].set_ylabel('#%s  F(<r)/Fmax' % r['NUMBER'], fontsize=7)
        ax[k][1].set_ylabel('eta(r)', fontsize=7)
        ax[k][1].set_ylim(0, 1.2)
        ax[k][0].tick_params(labelsize=6); ax[k][1].tick_params(labelsize=6)
    ax[-1][0].set_xlabel('semi-major axis (px)', fontsize=7); ax[-1][1].set_xlabel('semi-major axis (px)', fontsize=7)
    fig.tight_layout()
    fig.savefig(path, dpi=80)
    plt.close(fig)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('image')
    ap.add_argument('--catalog', required=True)
    ap.add_argument('--mask', default='')
    ap.add_argument('--work', default='')
    ap.add_argument('--eta', type=float, default=0.2)
    ap.add_argument('--eta-lo', type=float, default=0.1)
    ap.add_argument('--eta-hi', type=float, default=0.3)
    ap.add_argument('--kron-scale', type=float, default=2.5)
    ap.add_argument('--kron-min', type=float, default=3.5)
    ap.add_argument('--smooth-frac', type=float, default=0.25)
    ap.add_argument('--rmax', type=float, default=150.0)
    ap.add_argument('--rmax-scale', type=float, default=10.0, help='search radius = scale x A_IMAGE x KRON_RADIUS / 3.5 (px, capped by --rmax)')
    ap.add_argument('--circular', action='store_true')
    ap.add_argument('--mask-neighbours', action='store_true')
    ap.add_argument('--neighbour-radius', type=float, default=3.0)
    ap.add_argument('--neighbour-min-frac', type=float, default=0.03)
    ap.add_argument('--max-sources', type=int, default=500)
    ap.add_argument('--min-flux-radius', type=float, default=0.0, help='skip objects with FLUX_RADIUS below this (px)')
    ap.add_argument('--curves', type=int, default=4, help='growth-curve plot for this many brightest objects (0 = none)')
    ap.add_argument('--mag-zeropoint', type=float, default=25.0)
    ap.add_argument('--n-workers', type=int, default=0)
    a = ap.parse_args(argv)
    data0, hdr = imageio.load_image(a.image)
    mask0 = imageio.load_mask(a.mask, data0.shape) if a.mask and os.path.isfile(a.mask) else None
    sep = X._sep()
    d = np.array(data0, np.float64)
    bad = ~np.isfinite(d) | (d == 0)
    if mask0 is not None:
        bad |= mask0
    d[~np.isfinite(d)] = 0.0
    bkg = sep.Background(d, mask=bad.astype(np.uint8), bw=64, bh=64, fw=3, fh=3)
    data = d - bkg.back()
    rms = float(np.median(bkg.rms()[~bad])) if (~bad).any() else imageio.robust_sigma(data)
    cols, rows = tsvio.read_catalog(a.catalog)
    geo = []
    for r in rows:
        A, q, th = shape_of(r, a.circular)
        fl = num(r, 'FLUX_AUTO', 1.0)
        geo.append(dict(x=num(r, 'X_IMAGE') - 1, y=num(r, 'Y_IMAGE') - 1, A=A, q=q, th=th, flux=fl if np.isfinite(fl) else 1.0))
    from scipy.spatial import cKDTree
    pts = np.array([[g['x'] if np.isfinite(g['x']) else -1e9, g['y'] if np.isfinite(g['y']) else -1e9] for g in geo]) if geo else None
    idx = [i for i in range(len(rows)) if np.isfinite(geo[i]['x']) and (a.min_flux_radius <= 0 or num(rows[i], 'FLUX_RADIUS', 99) >= a.min_flux_radius)]
    idx.sort(key=lambda i: -geo[i]['flux'])
    idx = idx[:a.max_sources]
    a.want_curve = a.curves > 0 and bool(a.work)
    G.update(args=a, rows=rows, data=data, mask=(bad if mask0 is not None else None), rms=rms, geo=geo, tree=cKDTree(pts) if pts is not None and len(pts) else None)
    sys.stderr.write('morph_ext: %d of %d objects, rms %.4g, eta %.2f (%.2f, %.2f), Kron %.1f x r_k (min %.1f px), smoothness box %.2f R_P\n' % (len(idx), len(rows), rms, a.eta, a.eta_lo, a.eta_hi,
                     a.kron_scale, a.kron_min, a.smooth_frac))
    nw = a.n_workers if a.n_workers > 0 else min(8, os.cpu_count() or 1)
    if nw > 1 and len(idx) > 4:
        with mp.get_context('fork').Pool(nw) as pool:
            results = pool.map(one, idx, chunksize=4)
    else:
        results = [one(i) for i in idx]
    by = {rr['NUMBER']: rr for rr in results}
    if a.work:
        os.makedirs(a.work, exist_ok=True)
        recs = []
        for rr in results:
            if 'curve' in rr:
                r_, F_, (re_, eta_) = rr['curve']
                eta_i = np.interp(r_, re_, eta_, left=np.nan, right=np.nan)
                for k in range(0, len(r_), 4):
                    recs.append(dict(NUMBER=rr['NUMBER'], R=r_[k], FLUX=F_[k], ETA=eta_i[k]))
        tsvio.write_table(os.path.join(a.work, 'morph_ext_growth.tsv'), ['NUMBER', 'R', 'FLUX', 'ETA'], recs)
        if a.curves > 0:
            plot_curves(results, os.path.join(a.work, 'morph_ext_curves.png'), a.curves)
    out = [(r['NUMBER'], {c: (by[r['NUMBER']].get(c) if r['NUMBER'] in by else None) for c in COLUMNS}) for r in rows]
    nok = sum(1 for rr in results if rr.get('MX_RP') is not None and np.isfinite(rr.get('MX_RP')))
    sys.stderr.write('morph_ext: %d with a Petrosian radius\n' % nok)
    tsvio.write_columns(COLUMNS, out)
    return 0


if __name__ == '__main__':
    sys.exit(main())
