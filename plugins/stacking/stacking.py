#!/usr/bin/env python3
"""Stacking and mean radial profiles: positions -> aligned cutouts -> median/mean with bootstrap errors -> profile.

    stacking.py IMAGE --work DIR (--catalog TSV [--numbers 1,5-9] [--sel-col MAG_AUTO --sel-min 26 --sel-max 28]
                                  | --positions FILE) [--mask FITS] [--bkg-file FITS] [--half 30] [--align sub|int]
        [--bkg annulus|none|image] [--norm none|flux|peak|column] [--method mean|median|wmean|clipmean] [--n-boot 200]
        [--mask-detect 3] [--mask-catalog 1.5] [--null 300 [--subtract-null]] [--scale-col FLUX_RADIUS --scale-unit 5]
        [--angle-col THETA_IMAGE] [--q 1 --pa 0] [--mag-zeropoint 25 --pixel-scale 0.06] [--save-cube]

With --catalog the add_columns contract is followed (stdout = NUMBER + ST_* columns of the selected rows: ST_USED,
ST_BKG, ST_NORM, ST_APFLUX, ST_MASKFRAC).  Files in DIR/stacking_*: stack.fits, err.fits (bootstrap per-pixel error),
cutouts.fits (--save-cube), profile.tsv, summary.json, plot.png.
Positions file: columns X_IMAGE/Y_IMAGE (1-based), x/y, or RA/DEC (ALPHA_J2000/DELTA_J2000) in degrees; otherwise the first
two numeric columns are 1-based pixel coordinates.
"""
import argparse
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)
import warnings  # noqa: E402
warnings.filterwarnings('ignore')

from ogfkit import tsvio, imageio, stacking as st, meta as ometa  # noqa: E402

COLUMNS = ['ST_USED', 'ST_BKG', 'ST_NORM', 'ST_APFLUX', 'ST_MASKFRAC']


def read_positions(path, header=None):
    """-> (x0, y0 0-based, extra dict name->array)."""
    with open(path, encoding='utf-8') as f:
        lines = [l.strip() for l in f if l.strip() and not l.startswith('#')]
    sep = '\t' if '\t' in lines[0] else (',' if ',' in lines[0] else None)
    rows = [l.split(sep) for l in lines]
    first = rows[0]
    try:
        [float(v) for v in first[:2]]
        names = None
    except ValueError:
        names = [c.strip() for c in first]
        rows = rows[1:]
    arr = {}
    if names is None:
        for j in range(len(rows[0])):
            arr['COL%d' % j] = np.array([tsvio.fnum(r[j]) for r in rows])
        x, y = arr['COL0'] - 1, arr['COL1'] - 1
        return x, y, arr
    for j, n in enumerate(names):
        arr[n] = np.array([tsvio.fnum(r[j]) if j < len(r) else np.nan for r in rows])
    up = {n.upper(): n for n in names}
    for xn, yn in (('X_IMAGE', 'Y_IMAGE'), ('X', 'Y')):
        if xn in up and yn in up:
            off = 1 if xn == 'X_IMAGE' else 1
            return arr[up[xn]] - off, arr[up[yn]] - off, arr
    for an, dn in (('RA', 'DEC'), ('ALPHA_J2000', 'DELTA_J2000')):
        if an in up and dn in up:
            from astropy.wcs import WCS
            if header is None:
                raise ValueError('RA/DEC positions need an image with a WCS')
            w = WCS(header)
            x, y = w.all_world2pix(arr[up[an]], arr[up[dn]], 0)
            return np.asarray(x), np.asarray(y), arr
    raise ValueError('positions file: no X_IMAGE/Y_IMAGE, X/Y or RA/DEC columns')


def neighbours_from_catalog(rows, mask_scale, mask_min):
    x = np.array([tsvio.fnum(r.get('X_IMAGE')) - 1 for r in rows])
    y = np.array([tsvio.fnum(r.get('Y_IMAGE')) - 1 for r in rows])
    A = np.array([tsvio.fnum(r.get('A_IMAGE'), 2.0) for r in rows])
    K = np.array([tsvio.fnum(r.get('KRON_RADIUS'), 3.0) for r in rows])
    ext = np.where(np.isfinite(A), A, 2.0) * np.where(np.isfinite(K), K, 3.0)
    return np.column_stack([x, y, np.maximum(mask_min, mask_scale * ext)])


def write_profile(path, res, cfg, nullres, subtract_null, zp, pixscale, scale_unit):
    prof, perr = res['profile'].copy(), res['profile_err'].copy()
    recs = []
    corr = None
    if nullres is not None:
        nprof, nerr = nullres['profile'], nullres['profile_err']
        corr, cerr = prof - nprof, np.hypot(perr, nerr)
    use, usee = (corr, cerr) if (corr is not None and subtract_null) else (prof, perr)
    if zp is not None and pixscale:
        mu, mue, lim = st.surface_brightness(use, usee, zp, pixscale)
    for b in range(len(prof)):
        r = dict(R=res['r'][b], R_LO=res['edges'][b], R_HI=res['edges'][b + 1], N_PIX=res['npix'][b], I=prof[b], I_ERR=perr[b])
        if scale_unit:
            r['R_UNIT'] = res['r'][b] / scale_unit
        if nullres is not None:
            r.update(NULL=nprof[b], NULL_ERR=nerr[b], I_CORR=corr[b], I_CORR_ERR=cerr[b])
        r['SNR'] = use[b] / usee[b] if usee[b] > 0 else float('nan')
        if zp is not None and pixscale:
            r.update(MU=mu[b], MU_ERR=mue[b], MU_LIM3=lim[b])
        recs.append(r)
    cols = ['R'] + (['R_UNIT'] if scale_unit else []) + ['R_LO', 'R_HI', 'N_PIX', 'I', 'I_ERR']
    if nullres is not None:
        cols += ['NULL', 'NULL_ERR', 'I_CORR', 'I_CORR_ERR']
    cols += ['SNR'] + (['MU', 'MU_ERR', 'MU_LIM3'] if (zp is not None and pixscale) else [])
    tsvio.write_table(path, cols, recs)
    return use, usee


def make_plot(path, res, nullres, use, usee, cfg, title, zp, pixscale):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 2, figsize=(10, 4.2), gridspec_kw=dict(width_ratios=[1, 1.5]))
    s = res['stack']
    v = np.nanpercentile(s, [1, 99.5])
    ax[0].imshow(np.arcsinh(s / max(abs(v[1]) * 0.05, 1e-30)), origin='lower', cmap='gray')
    ax[0].add_patch(plt.Circle((s.shape[1] // 2, s.shape[0] // 2), cfg.ap_r, fill=False, color='c', lw=0.8))
    ax[0].set_title('%s of %d (aperture r=%g px)' % (cfg.method, res['n'], cfg.ap_r), fontsize=9)
    r = res['r']
    ax[1].errorbar(r, np.where(use > 0, use, np.nan), usee, fmt='o-', ms=3, color='k', label='stack')
    neg = use <= 0
    if neg.any():
        ax[1].plot(r[neg], usee[neg] * 3, 'v', color='r', ms=4, label='3 sigma (I<=0)')
    if nullres is not None:
        ax[1].errorbar(r, np.abs(nullres['profile']), nullres['profile_err'], fmt='s', ms=3, color='gray', label='|null|')
    ax[1].set_xscale('log')
    ax[1].set_yscale('log')
    ax[1].set_xlabel('r (px)')
    ax[1].set_ylabel('I (per pixel)')
    ax[1].legend(fontsize=8)
    ax[1].set_title(title, fontsize=9)
    fig.tight_layout()
    fig.savefig(path, dpi=100)
    plt.close(fig)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('image')
    ap.add_argument('--work', required=True)
    ap.add_argument('--catalog')
    ap.add_argument('--positions')
    ap.add_argument('--numbers', default='')
    ap.add_argument('--sel-col', default='')
    ap.add_argument('--sel-min', type=float, default=float('-inf'))
    ap.add_argument('--sel-max', type=float, default=float('inf'))
    ap.add_argument('--isolate', type=float, default=0.0, help='reject sources with another catalog source closer than this (px)')
    ap.add_argument('--mask')
    ap.add_argument('--bkg-file')
    ap.add_argument('--half', type=int, default=30)
    ap.add_argument('--align', default='sub', choices=['sub', 'int'])
    ap.add_argument('--bkg', default='annulus', choices=['none', 'annulus', 'image'])
    ap.add_argument('--bkg-in', type=float, default=0.7)
    ap.add_argument('--norm', default='none', choices=['none', 'flux', 'peak', 'column'])
    ap.add_argument('--norm-r', type=float, default=5.0)
    ap.add_argument('--norm-col', default='')
    ap.add_argument('--method', default='mean', choices=list(st.METHODS))
    ap.add_argument('--clip', type=float, default=3.0)
    ap.add_argument('--weights', default='none', choices=['none', 'invvar', 'column'])
    ap.add_argument('--weight-col', default='')
    ap.add_argument('--scale-col', default='')
    ap.add_argument('--scale-unit', type=float, default=5.0, help='output pixels per scale-col size')
    ap.add_argument('--angle-col', default='')
    ap.add_argument('--mask-detect', type=float, default=0.0)
    ap.add_argument('--protect-r', type=float, default=3.0)
    ap.add_argument('--mask-catalog', type=float, default=0.0)
    ap.add_argument('--mask-min', type=float, default=3.0)
    ap.add_argument('--min-valid', type=float, default=0.6)
    ap.add_argument('--n-boot', type=int, default=200)
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--null', type=int, default=0)
    ap.add_argument('--subtract-null', action='store_true')
    ap.add_argument('--q', type=float, default=1.0)
    ap.add_argument('--pa', type=float, default=0.0)
    ap.add_argument('--ap-r', type=float, default=5.0)
    ap.add_argument('--estimator', default='mean', choices=['mean', 'median'])
    ap.add_argument('--rlin', type=float, default=6.0)
    ap.add_argument('--rfac', type=float, default=1.25)
    ap.add_argument('--mag-zeropoint', type=float, default=float('nan'))
    ap.add_argument('--pixel-scale', type=float, default=float('nan'))
    ap.add_argument('--save-cube', action='store_true')
    ap.add_argument('--meta-out')
    a = ap.parse_args(argv)

    os.makedirs(a.work, exist_ok=True)
    W = lambda n: os.path.join(a.work, 'stacking_' + n)
    data, hdr = imageio.load_image(a.image)
    bad = imageio.load_mask(a.mask, data.shape) if a.mask and os.path.isfile(a.mask) else None
    bkg_map = None
    if a.bkg_file and os.path.isfile(a.bkg_file):
        bkg_map = imageio.load_image(a.bkg_file)[0]
        if bkg_map.shape != data.shape:
            raise SystemExit('bkg-file shape differs from the image')
    elif a.bkg == 'image':
        raise SystemExit('--bkg image needs --bkg-file')

    rows, extra, numbers = None, {}, None
    cat_x = cat_y = None
    if a.catalog and os.path.isfile(a.catalog) and not a.positions:
        cols, rows = tsvio.read_catalog(a.catalog)
        cat_x = np.array([tsvio.fnum(r.get('X_IMAGE')) - 1 for r in rows])
        cat_y = np.array([tsvio.fnum(r.get('Y_IMAGE')) - 1 for r in rows])
        sel = np.isfinite(cat_x) & np.isfinite(cat_y)
        if a.numbers:
            want = set(tsvio.parse_numbers(a.numbers))
            sel &= np.array([int(tsvio.fnum(r.get('NUMBER'), -1)) in want for r in rows])
        if a.sel_col and (a.sel_min > -1e29 or a.sel_max < 1e29):
            v = np.array([tsvio.fnum(r.get(a.sel_col)) for r in rows])
            sel &= np.isfinite(v) & (v >= a.sel_min) & (v <= a.sel_max)
        ids = np.where(sel)[0]
        xs, ys = cat_x[ids], cat_y[ids]
        get = lambda c: np.array([tsvio.fnum(rows[i].get(c)) for i in ids])
    elif a.positions:
        xs, ys, extra = read_positions(a.positions, hdr)
        ids = np.arange(len(xs))
        get = lambda c: extra[c]
        if a.catalog and os.path.isfile(a.catalog):
            cols, rows = tsvio.read_catalog(a.catalog)
            cat_x = np.array([tsvio.fnum(r.get('X_IMAGE')) - 1 for r in rows])
            cat_y = np.array([tsvio.fnum(r.get('Y_IMAGE')) - 1 for r in rows])
    else:
        raise SystemExit('need --catalog or --positions')
    if len(xs) == 0:
        raise SystemExit('no positions selected')

    keep = np.ones(len(xs), bool)
    if a.isolate > 0 and cat_x is not None:
        from scipy.spatial import cKDTree
        g = np.isfinite(cat_x) & np.isfinite(cat_y)
        t = cKDTree(np.column_stack([cat_x[g], cat_y[g]]))
        d, _ = t.query(np.column_stack([xs, ys]), k=min(3, int(g.sum())))
        d = np.atleast_2d(d.reshape(len(xs), -1))
        d = np.where(d < 1.5, np.inf, d)               # the source itself
        keep = d.min(axis=1) >= a.isolate
    cfg_kw = dict(half=a.half, align=a.align, bkg=a.bkg, bkg_in=a.bkg_in, norm=a.norm, norm_r=a.norm_r, method=a.method, clip=a.clip,
                  weights=a.weights, n_boot=a.n_boot, seed=a.seed, q=a.q, pa=a.pa, ap_r=a.ap_r, estimator=a.estimator,
                  min_valid=a.min_valid, mask_detect=a.mask_detect, protect_r=a.protect_r, mask_cat_scale=a.mask_catalog, mask_cat_min=a.mask_min)
    cfg = st.StackConfig(**cfg_kw)
    scales = angles = None
    if a.scale_col:
        R = get(a.scale_col)
        scales = np.where(np.isfinite(R) & (R > 0), R / a.scale_unit, np.nan)
    if a.angle_col:
        angles = get(a.angle_col)         # THETA_IMAGE (deg from +x towards +y): the grid is rotated so the major axis becomes +x
    kw = dict(scales=scales, angles=angles, bkg_map=bkg_map)
    if a.norm == 'column':
        kw['norm_values'] = get(a.norm_col)
    if a.weights == 'column':
        kw['weight_values'] = get(a.weight_col)
    neigh = None
    if a.mask_catalog > 0 and rows is not None:
        neigh = neighbours_from_catalog(rows, a.mask_catalog, a.mask_min)
        kw['neigh'] = neigh
        if not a.positions:
            kw['self_index'] = ids
    sx, sy = xs[keep], ys[keep]
    for k in ('scales', 'angles', 'norm_values', 'weight_values'):
        if kw.get(k) is not None:
            kw[k] = np.asarray(kw[k])[keep]
    if kw.get('self_index') is not None:
        kw['self_index'] = kw['self_index'][keep]
    edges = st.radial_bins(a.half, a.rlin, a.rfac)
    cut = st.make_cutouts(data, bad, sx, sy, cfg, **kw)
    prof = st.Profiler(2 * a.half + 1, edges, cfg.q, cfg.pa)
    if not cut['ok'].any():
        why = {}
        for r_ in cut['reason']:
            k_ = r_.split(' ')[0] + ' ' + (r_.split(' ')[1] if len(r_.split(' ')) > 1 else '')
            why[k_] = why.get(k_, 0) + 1
        raise SystemExit('stacking: no valid cutouts among %d positions (%s)' % (len(sx), ', '.join('%s: %d' % kv for kv in sorted(why.items()))))
    res = st.stack_cutouts(cut['cube'], cut['ok'], cfg, cut['weight'], profiler=prof)
    nullres = None
    if a.null > 0:
        avoid = np.column_stack([cat_x, cat_y]) if cat_x is not None else None
        nx_, ny_ = st.null_positions(data, bad, a.null, cfg, avoid=avoid)
        nkw = {k: v for k, v in kw.items() if k in ('bkg_map', 'neigh')}
        nc = st.make_cutouts(data, bad, nx_, ny_, cfg, **nkw)
        if nc['ok'].sum() >= 5:
            nullres = st.stack_cutouts(nc['cube'], nc['ok'], cfg, nc['weight'], profiler=prof)
    zp = a.mag_zeropoint if np.isfinite(a.mag_zeropoint) else None
    ps = a.pixel_scale if np.isfinite(a.pixel_scale) else None
    res['null'] = nullres
    use, usee = write_profile(W('profile.tsv'), res, cfg, nullres, a.subtract_null, zp, ps, a.scale_unit if a.scale_col else None)
    imageio.save_fits(W('stack.fits'), np.nan_to_num(res['stack']).astype(np.float32), extra=dict(NSTACK=res['n'], STMETHOD=a.method))
    imageio.save_fits(W('err.fits'), np.nan_to_num(res['err_map']).astype(np.float32))
    if a.save_cube:
        imageio.save_fits(W('cutouts.fits'), cut['cube'][cut['ok']])
    summ = st.summarize(res, cfg, zp)
    if a.null > 0 and nullres is None:
        summ['null_warning'] = 'fewer than 5 usable blank positions: no null stack'
    summ.update(n_selected=int(len(xs)), n_isolated=int(keep.sum()), n_rejected=int(keep.sum() - res['n']),
                rejection_reasons={k: int(v) for k, v in zip(*np.unique([r for r, o in zip(cut['reason'], cut['ok']) if not o], return_counts=True))},
                half=a.half, align=a.align, norm=a.norm, subtract_null=bool(a.subtract_null and nullres is not None),
                bkg_median=float(np.nanmedian(cut['bkg'][cut['ok']])), rms_median=float(np.nanmedian(cut['rms'][cut['ok']])))
    if nullres is not None and a.subtract_null:
        apc = res['aper'] - nullres['aper']
        summ['aperture_sum_corrected'] = float(apc)
        summ['aperture_err_corrected'] = float(math.hypot(res['aper_err'], nullres['aper_err']))
    with open(W('summary.json'), 'w') as f:
        json.dump(summ, f, indent=1)
    try:
        make_plot(W('plot.png'), res, nullres, use, usee, cfg, os.path.basename(a.image), zp, ps)
    except Exception as e:                                           # plotting must never fail the step
        print('plot failed: %s' % e, file=sys.stderr)
    if a.meta_out:
        ometa.update(a.meta_out, 'stacking', dict(n=res['n'], aperture_sum=summ['aperture_sum'], aperture_err=summ['aperture_err'], method=a.method))
    if rows is not None and not a.positions:
        out = []
        used = {int(ids[keep][j]): j for j in range(len(sx))}
        ap_f = None
        for j, i in enumerate(ids):
            d = dict()
            num = rows[i].get('NUMBER', str(i + 1))
            k = used.get(int(i))
            if k is None:
                d['ST_USED'] = 0
            else:
                ok = bool(cut['ok'][k])
                d['ST_USED'] = int(ok)
                d['ST_MASKFRAC'] = float(1 - cut['valid_frac'][k])
                if ok:
                    d['ST_BKG'] = float(cut['bkg'][k])
                    d['ST_NORM'] = float(cut['norm'][k])
                    c = cut['cube'][k] / cut['norm'][k]
                    d['ST_APFLUX'] = st.aperture(np.asarray(c, float), a.ap_r)
            out.append((num, d))
        tsvio.write_columns(COLUMNS, out)
    else:
        print('stacked %d of %d positions; aperture sum %.5g +- %.3g' % (res['n'], len(xs), summ['aperture_sum'], summ['aperture_err']), file=sys.stderr)
    return 0


if __name__ == '__main__':
    sys.exit(main())
