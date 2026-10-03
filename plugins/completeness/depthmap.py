#!/usr/bin/env python3
"""Depth and completeness maps of an image, built on the completeness injection / recovery engine.

    depthmap.py IMAGE --work DIR --mode depth|compmap [--catalog TSV] [--mask FITS] [--mag-zeropoint 25] [--pixel-scale 0.06]
        [--aper-radius 3] [--nsigma 5] [--apcorr 0] [--box-arcsec 10] [--box-nsigma 3] [--tile 128] [--bw 64]
        [--rms-file F | --weight-file F]
      compmap only: [--regions grid|rms|file] [--grid 3x3] [--rms-classes 4] [--region-file F] [--maps-at 26,27,28]
        [--kind star|galaxy] [--detector sep] [--mag-min --mag-max --n-bins --per-bin --per-image --psf-fwhm --match-radius ...]

depth:    local pixel rms -> point-source limiting magnitude map (nsigma, aperture radius) and surface-brightness limit map (box x box arcsec^2,
          box-nsigma) through the blank-aperture / blank-box noise laws; per-tile empirical check; area-vs-depth table.
compmap:  depth maps + injection / recovery over the whole image; recovered fractions are binned by region (grid / depth classes / label
          image) -> lim50 / lim90 per region, completeness curve as a function of (m - local depth) (the "universal" curve),
          per-pixel lim50 / lim90 = depth + offset, completeness maps at chosen magnitudes.
Files in DIR/depth_* and DIR/compmap_*; catalog columns (add_columns): DEPTH_RMS, DEPTH_LIM, DEPTH_SBLIM (mode depth) and COMPL_REGION, COMPL_LOC,
COMPL_LIM50_LOC, COMPL_LIM90_LOC (compmap).
"""
import argparse
import json
import math
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..', '..'))
for p_ in (ROOT, HERE):
    if p_ not in sys.path:
        sys.path.insert(0, p_)
import warnings  # noqa: E402
warnings.filterwarnings('ignore')

from ogfkit import tsvio, imageio, depth as dp, meta as ometa, noise as nz  # noqa: E402
import completeness as cp  # noqa: E402

DEPTH_COLS = ['DEPTH_RMS', 'DEPTH_LIM', 'DEPTH_SBLIM']
COMP_COLS = ['COMPL_REGION', 'COMPL_LOC', 'COMPL_LIM50_LOC', 'COMPL_LIM90_LOC']


def clean(o):
    return cp.json_clean(o)


def save_map(path, arr, hdr, **cards):
    imageio.save_fits(path, np.asarray(arr, np.float32), header=hdr, extra={k: v for k, v in cards.items() if v is not None})


def do_depth(data, hdr, bad, a, W):
    P = dp.prepare(data, bad, bw=a.bw)
    valid = ~P['bad']
    sigma_pix = nz.robust_std(P['sub'][~P['mask']][:2000000])
    rms = P['rms']
    if a.rms_file and os.path.isfile(a.rms_file):
        rms = dp.rms_from_file(imageio.load_image(a.rms_file)[0], valid, P['rms'], 'rms')
    elif a.weight_file and os.path.isfile(a.weight_file):
        rms = dp.rms_from_file(imageio.load_image(a.weight_file)[0], valid, P['rms'], 'weight')
    valid &= np.isfinite(rms) & (rms > 0)
    valid &= rms > 0.05 * np.median(rms[valid])           # empty / zero-weight margins
    r = a.aper_radius
    f_ap, n_ap = dp.circle_factor(P['sub'], P['mask'], rms, r, a.n_aper, a.seed)
    if not np.isfinite(f_ap):
        raise SystemExit('depth: too few blank apertures (image mostly masked?)')
    radii = sorted({max(1.0, round(r * f, 2)) for f in (0.5, 0.75, 1.0, 1.5, 2.0, 3.0)})
    law = dp.aperture_law(P['sub'], P['mask'], rms, radii, max(300, a.n_aper // 2), a.seed)
    mag = dp.mag_limit_map(rms, f_ap, a.mag_zeropoint, a.nsigma, a.apcorr)
    mag = np.where(valid, mag, np.nan)
    out = dict(sigma_pix=float(sigma_pix), factor=dict(aperture=f_ap, n_apertures=n_ap),
               law=dict(alpha=law['alpha'], beta=law['beta'], ok=bool(law.get('ok'))), aper_radius_px=r, nsigma=a.nsigma,
               apcorr=a.apcorr, zeropoint=a.mag_zeropoint, valid_fraction=float(valid.mean()), source_masked_fraction=float(P['mask'].mean()))
    v = mag[valid]
    out['mag_limit'] = dict(median=float(np.median(v)), p5=float(np.percentile(v, 5)), p95=float(np.percentile(v, 95)), min=float(v.min()), max=float(v.max()))
    L = a.box_arcsec / a.pixel_scale if a.pixel_scale > 0 else 0.0
    sb = None
    f_box = float('nan')
    if L >= 4:
        f_box, n_box = dp.box_factor(P['sub'], P['mask'], rms, int(round(L)))
        sizes = [s_ for s_ in sorted({int(round(L * f)) for f in (0.2, 0.35, 0.5, 0.7, 1.0)}) if s_ >= 3]
        blaw = dp.box_law(P['sub'], P['mask'], rms, sizes)
        f_law = blaw['gamma'] * float(int(round(L)) ** 2) ** blaw['delta'] if blaw['ok'] else float('nan')
        n_indep = n_box / 4.0                                   # boxes on a half-box stride overlap
        f_direct = f_box
        if np.isfinite(f_law) and (not np.isfinite(f_box) or n_indep < 300):
            f_box = f_law                                        # few independent boxes: use the power law through the smaller boxes
        if np.isfinite(f_box):
            sb = dp.sb_limit_map(rms, f_box, int(round(L)), a.pixel_scale, a.mag_zeropoint, a.box_nsigma)
            sb = np.where(valid, sb, np.nan)
            out['sb_limit'] = dict(box_arcsec=a.box_arcsec, nsigma=a.box_nsigma, factor=f_box, factor_direct=f_direct, factor_law=f_law, n_boxes=n_box, gamma=blaw['gamma'], delta=blaw['delta'], law_ok=bool(blaw['ok']),
                                   median=float(np.nanmedian(sb)), p5=float(np.nanpercentile(sb, 5)), p95=float(np.nanpercentile(sb, 95)))
        else:
            out['sb_limit_warning'] = 'fewer than 12 blank boxes of %.0f arcsec: no surface-brightness limit' % a.box_arcsec
    tiles = dp.tile_table(P['sub'], P['mask'], rms, f_ap, f_box, a.tile, r, int(round(L)) if sb is not None else 0, n_aper=a.tile_aper, seed=a.seed)
    rat = np.array([t['ratio'] for t in tiles])
    rat = rat[np.isfinite(rat)]
    out['tile_check'] = dict(n_tiles=len(tiles), n_valid=int(len(rat)), ratio_median=float(np.median(rat)) if len(rat) else None,
                             ratio_std=float(rat.std()) if len(rat) else None,
                             ratio_min=float(rat.min()) if len(rat) else None, ratio_max=float(rat.max()) if len(rat) else None)
    area = dp.area_depth(mag, valid)
    out['area_vs_depth'] = [dict(mag=m_, fraction_deeper=f_) for m_, f_ in area]
    save_map(W('depth_rms_map.fits'), np.where(valid, rms, np.nan), hdr)
    save_map(W('depth_mag_map.fits'), mag, hdr, DPNSIG=a.nsigma, DPRAD=r, DPZP=a.mag_zeropoint)
    if sb is not None:
        save_map(W('depth_sb_map.fits'), sb, hdr, DPBOX=a.box_arcsec, DPBNSIG=a.box_nsigma)
    cols = ['tile', 'xc', 'yc', 'valid_frac', 'rms_pix', 'n_aper', 'sigma_emp', 'sigma_emp_err', 'sigma_model', 'ratio', 'n_box', 'box_sigma_emp', 'box_sigma_model', 'box_ratio']
    for t in tiles:
        t['mag_model'] = float(a.mag_zeropoint - 2.5 * math.log10(a.nsigma * t['sigma_model']) - a.apcorr) if t['sigma_model'] > 0 else float('nan')
        t['mag_emp'] = float(a.mag_zeropoint - 2.5 * math.log10(a.nsigma * t['sigma_emp']) - a.apcorr) if t['sigma_emp'] > 0 else float('nan')
    tsvio.write_table(W('depth_tiles.tsv'), cols + ['mag_model', 'mag_emp'], tiles)
    tsvio.write_table(W('depth_area.tsv'), ['mag', 'fraction_deeper'], out['area_vs_depth'])
    return out, dict(P=P, rms=rms, mag=mag, sb=sb, valid=valid, law=law, tiles=tiles)


def sample(m, x, y):
    return nz.sample_map(m, x, y)


def do_compmap(data, hdr, bad, a, W, ctx, out):
    mask = ~ctx['valid'] if bad is None else (bad | ~ctx['valid'])
    det = cp.load_detector(a.detector, zp=a.mag_zeropoint, thresh=a.detect_thresh, minarea=a.detect_minarea, local_rms=not a.global_threshold)
    psf = None
    res = cp.run_completeness(data, det, mask=mask, kind=a.kind, mag_min=a.mag_min, mag_max=a.mag_max, n_bins=a.n_bins, per_bin=a.per_bin,
                              per_image=a.per_image, zp=a.mag_zeropoint, psf=psf, psf_fwhm=a.psf_fwhm, match_radius=a.match_radius,
                              re_edges=[float(v) for v in a.re_bins.split(',')], seed=a.seed, pixel_scale=(a.pixel_scale or None),
                              avoid_detected=not a.allow_blends, n_workers=a.n_workers, return_records=True)
    recs = res.pop('records')
    mag, valid = ctx['mag'], ctx['valid']
    gx, gy = [int(v) for v in a.grid.lower().split('x')]
    lab_file = imageio.load_image(a.region_file)[0] if (a.regions == 'file' and a.region_file) else None
    lab = dp.region_labels(a.regions, data.shape, depth=mag, valid=valid, grid=(gx, gy), classes=a.rms_classes, label_file=lab_file)
    x = np.array([r['x'] for r in recs]); y = np.array([r['y'] for r in recs])
    m = np.array([r['mag'] for r in recs])
    rec = np.array([r['recovered'] for r in recs])
    reg = sample(lab, x, y)
    dep = sample(mag, x, y)
    edges = np.linspace(a.mag_min, a.mag_max, a.n_bins + 1)
    area_pix = int(valid.sum())
    regions = []
    for k in sorted(set(int(v) for v in np.unique(lab) if v > 0)):
        sel = [i for i in range(len(recs)) if reg[i] == k]
        area_k = int((lab == k).sum())
        rr = dict(region=k, area_pix=area_k, depth_median=float(np.nanmedian(mag[lab == k])), n_inj=len(sel))
        if len(sel) >= a.min_region:
            sm = cp.summarize([recs[i] for i in sel], edges, 0, area_k, 0, 1, None)
            rr.update(lim50=sm['lim50'], lim90=sm['lim90'], n_rec=sm['n_recovered'], fit_ok=bool(sm['fit'].get('ok')))
        else:
            rr.update(lim50=float('nan'), lim90=float('nan'), n_rec=int(rec[[i for i in sel]].sum()) if sel else 0, fit_ok=False)
        regions.append(rr)
    # universal curve in x = m - depth
    ok = np.isfinite(dep)
    xs = m[ok] - dep[ok]
    xe = np.linspace(np.percentile(xs, 1), np.percentile(xs, 99), a.n_bins + 1)
    kk, nn, cc = [], [], []
    for i in range(a.n_bins):
        s = (xs >= xe[i]) & (xs < xe[i + 1] if i < a.n_bins - 1 else xs <= xe[i + 1])
        nn.append(int(s.sum())); kk.append(int(rec[ok][s].sum())); cc.append(0.5 * (xe[i] + xe[i + 1]))
    fit = cp.fit_curve(np.array(cc), np.array(kk), np.array(nn))
    fd = np.array(kk) / np.maximum(np.array(nn), 1)
    d50 = fit['m50'] if np.isfinite(fit['m50']) else cp.direct_limit(np.array(cc), fd, 0.5)
    d90 = fit['m90'] if np.isfinite(fit['m90']) else cp.direct_limit(np.array(cc), fd, 0.9)
    # collapse check: per-region lim50 - depth_median should be constant = d50
    off = np.array([r_['lim50'] - r_['depth_median'] for r_ in regions if np.isfinite(r_['lim50'])])
    coll = dict(n_regions=int(len(off)), offset_median=float(np.median(off)) if len(off) else None, offset_std=float(off.std()) if len(off) > 1 else None,
                universal_delta50=float(d50), universal_delta90=float(d90), universal_fit=clean(fit))
    # regression lim50 = c + s * depth over regions
    dd = np.array([r_['depth_median'] for r_ in regions if np.isfinite(r_['lim50'])])
    l50 = np.array([r_['lim50'] for r_ in regions if np.isfinite(r_['lim50'])])
    if len(dd) >= 3 and dd.std() > 0:
        s_, c_ = np.polyfit(dd, l50, 1)
        coll['regression_slope'] = float(s_); coll['regression_intercept'] = float(c_)
    lim50_pix = mag + d50
    lim90_pix = mag + d90
    save_map(W('compmap_lim50_map.fits'), lim50_pix, hdr)
    save_map(W('compmap_lim90_map.fits'), lim90_pix, hdr)
    save_map(W('compmap_regions.fits'), lab.astype(np.float32), hdr)
    if fit.get('ok') and np.isfinite(fit.get('mh', np.nan)):
        F = lambda xx: np.clip(cp.logistic(xx, fit['mh'], fit['w'], fit['fmax']), 0, 1)
    else:
        F = lambda xx: np.interp(xx, cc, fd)
    for mm in [float(v) for v in a.maps_at.split(',') if v.strip()]:
        fm = np.where(valid, F(mm - mag), np.nan)
        save_map(W('compmap_frac_m%g.fits' % mm), fm, hdr, COMPMAG=mm)
    coll['fraction_area_complete50'] = {}
    recs_out = dict(config=res['config'], global_lim50=res['lim50'], global_lim90=res['lim90'], n_injected=res['n_injected'], n_recovered=res['n_recovered'],
                    regions=regions, collapse=coll, universal_curve=dict(x=[float(v) for v in cc], n=nn, k=kk, frac=[float(v) for v in fd]),
                    false_positive=res['false_positive'])
    tsvio.write_table(W('compmap_regions.tsv'), ['region', 'area_pix', 'depth_median', 'n_inj', 'n_rec', 'lim50', 'lim90', 'fit_ok'], regions)
    tsvio.write_table(W('compmap_universal.tsv'), ['x', 'n', 'k', 'frac'], [dict(x=cc[i], n=nn[i], k=kk[i], frac=fd[i]) for i in range(len(cc))])
    ctx.update(lab=lab, lim50_pix=lim50_pix, lim90_pix=lim90_pix, F=F, d50=d50, d90=d90)
    return recs_out


def plot(path, ctx, out, comp):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    n = 3 if comp else 1
    fig, ax = plt.subplots(1, n, figsize=(4.6 * n, 4))
    ax = np.atleast_1d(ax)
    im = ax[0].imshow(ctx['mag'], origin='lower', cmap='viridis')
    plt.colorbar(im, ax=ax[0], fraction=0.046)
    ax[0].set_title('%g-sigma limiting mag (r=%g px)' % (out['nsigma'], out['aper_radius_px']), fontsize=9)
    if comp:
        c = out['compmap']
        rg = [r for r in c['regions'] if np.isfinite(r['lim50'])]
        ax[1].plot([r['depth_median'] for r in rg], [r['lim50'] for r in rg], 'ko', label='lim50')
        ax[1].plot([r['depth_median'] for r in rg], [r['lim90'] for r in rg], 'rs', label='lim90')
        xx = np.array(ax[1].get_xlim())
        ax[1].plot(xx, xx + c['collapse']['universal_delta50'], 'k--', lw=0.7)
        ax[1].set_xlabel('depth map (mag)'); ax[1].set_ylabel('injection limit (mag)'); ax[1].legend(fontsize=8)
        u = c['universal_curve']
        ax[2].plot(u['x'], u['frac'], 'ko')
        ax[2].set_xlabel('m - local depth'); ax[2].set_ylabel('completeness')
        ax[2].axhline(0.5, color='gray', lw=0.5)
    fig.tight_layout()
    fig.savefig(path, dpi=100)
    plt.close(fig)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('image')
    ap.add_argument('--work', required=True)
    ap.add_argument('--mode', default='depth', choices=['depth', 'compmap'])
    ap.add_argument('--catalog', default='')
    ap.add_argument('--mask', default='')
    ap.add_argument('--rms-file', default='')
    ap.add_argument('--weight-file', default='')
    ap.add_argument('--mag-zeropoint', type=float, default=25.0)
    ap.add_argument('--pixel-scale', type=float, default=0.0)
    ap.add_argument('--aper-radius', type=float, default=3.0)
    ap.add_argument('--nsigma', type=float, default=5.0)
    ap.add_argument('--apcorr', type=float, default=0.0)
    ap.add_argument('--box-arcsec', type=float, default=10.0)
    ap.add_argument('--box-nsigma', type=float, default=3.0)
    ap.add_argument('--tile', type=int, default=128)
    ap.add_argument('--tile-aper', type=int, default=200)
    ap.add_argument('--n-aper', type=int, default=3000)
    ap.add_argument('--bw', type=int, default=64)
    ap.add_argument('--regions', default='grid', choices=['grid', 'rms', 'depth', 'file'])
    ap.add_argument('--grid', default='3x3')
    ap.add_argument('--rms-classes', type=int, default=4)
    ap.add_argument('--region-file', default='')
    ap.add_argument('--min-region', type=int, default=60)
    ap.add_argument('--maps-at', default='')
    ap.add_argument('--detector', default='sep')
    ap.add_argument('--kind', default='star', choices=['star', 'galaxy'])
    ap.add_argument('--mag-min', type=float, default=22.0)
    ap.add_argument('--mag-max', type=float, default=29.0)
    ap.add_argument('--n-bins', type=int, default=14)
    ap.add_argument('--per-bin', type=int, default=200)
    ap.add_argument('--per-image', type=int, default=25)
    ap.add_argument('--psf-fwhm', type=float, default=3.0)
    ap.add_argument('--match-radius', type=float, default=3.0)
    ap.add_argument('--re-bins', default='2,4,8,16')
    ap.add_argument('--detect-thresh', type=float, default=1.5)
    ap.add_argument('--detect-minarea', type=int, default=5)
    ap.add_argument('--global-threshold', action='store_true', help='sep detector: threshold relative to the global rms (default: local rms map)')
    ap.add_argument('--allow-blends', action='store_true')
    ap.add_argument('--n-workers', type=int, default=0)
    ap.add_argument('--seed', type=int, default=1)
    ap.add_argument('--meta-out', default='')
    a = ap.parse_args(argv)
    os.makedirs(a.work, exist_ok=True)
    W = lambda n: os.path.join(a.work, n)
    data, hdr = imageio.load_image(a.image)
    bad = imageio.load_mask(a.mask, data.shape) if a.mask and os.path.isfile(a.mask) else None
    out, ctx = do_depth(data, hdr, bad, a, W)
    if a.mode == 'compmap':
        out['compmap'] = do_compmap(data, hdr, bad, a, W, ctx, out)
    out['image'] = os.path.basename(a.image)
    with open(W('depth_summary.json' if a.mode == 'depth' else 'compmap_summary.json'), 'w') as f:
        json.dump(clean(out), f, indent=1, allow_nan=False)
    try:
        plot(W('depth_plot.png' if a.mode == 'depth' else 'compmap_plot.png'), ctx, out, a.mode == 'compmap')
    except Exception as e:                                                  # never fail the step for the plot
        sys.stderr.write('plot failed: %s\n' % e)
    ml = out['mag_limit']
    sys.stderr.write('depth: %g-sigma limiting mag (r=%g px) median %.2f [%.2f..%.2f]%s\n' % (a.nsigma, a.aper_radius, ml['median'], ml['p5'], ml['p95'],
                     (', SB limit %.2f' % out['sb_limit']['median']) if 'sb_limit' in out else ''))
    if a.catalog:
        cols, rows = tsvio.read_catalog(a.catalog)
        x = np.array([tsvio.fnum(r.get('X_IMAGE')) - 1 for r in rows]); y = np.array([tsvio.fnum(r.get('Y_IMAGE')) - 1 for r in rows])
        good = np.isfinite(x) & np.isfinite(y)
        xi = np.where(good, x, 0); yi = np.where(good, y, 0)
        d_rms = sample(ctx['rms'], xi, yi); d_lim = sample(ctx['mag'], xi, yi)
        d_sb = sample(ctx['sb'], xi, yi) if ctx['sb'] is not None else np.full(len(rows), np.nan)
        if a.mode == 'compmap':
            lab = sample(ctx['lab'], xi, yi)
            mags = np.array([tsvio.fnum(r.get('MAG_AUTO')) for r in rows])
            cl = ctx['F'](np.where(np.isfinite(mags), mags, 99.0) - d_lim)
            l50 = sample(ctx['lim50_pix'], xi, yi); l90 = sample(ctx['lim90_pix'], xi, yi)
        res_rows = []
        for i, r in enumerate(rows):
            d = {}
            if good[i]:
                if a.mode == 'depth':
                    d.update(DEPTH_RMS=float(d_rms[i]), DEPTH_LIM=float(d_lim[i]), DEPTH_SBLIM=float(d_sb[i]))
                else:
                    d.update(COMPL_REGION=int(lab[i]), COMPL_LOC=(float(cl[i]) if np.isfinite(mags[i]) else None), COMPL_LIM50_LOC=float(l50[i]), COMPL_LIM90_LOC=float(l90[i]))
            res_rows.append((r['NUMBER'], d))
        tsvio.write_columns(COMP_COLS if a.mode == 'compmap' else DEPTH_COLS, res_rows)
        if a.meta_out:
            ometa.update(a.meta_out, 'depth', clean(dict(mag_limit_median=ml['median'], nsigma=a.nsigma, aper_radius=a.aper_radius,
                         sb_limit_median=(out['sb_limit']['median'] if 'sb_limit' in out else None), mag_zeropoint=a.mag_zeropoint,
                         lim50_offset=(out['compmap']['collapse']['universal_delta50'] if 'compmap' in out else None))), nrows=len(rows))
    else:
        print(json.dumps(clean(dict(mag_limit=ml, sb_limit=out.get('sb_limit'), tile_check=out['tile_check']))))
    return 0


if __name__ == '__main__':
    sys.exit(main())
