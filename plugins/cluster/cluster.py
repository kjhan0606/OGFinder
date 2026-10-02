#!/usr/bin/env python3
"""Cluster and lensing tools for the current catalog.

    cluster.py --task members --catalog TSV --work DIR [--image FITS] [--mag-column C --blue-column C --red-column C ...]
    cluster.py --task arcs    --catalog TSV --image FITS --work DIR [...]
    cluster.py --task density --catalog TSV --work DIR [--image FITS] [...]

Every task follows the add_columns contract (stdout = NUMBER + the task's columns) and records what it found in --meta-out (key "cluster").
Files in DIR: cluster_summary.json (all tasks merge into it), cluster_rs.png (colour-magnitude diagram), cluster_density.fits (significance map),
cluster_peaks.tsv.  Positions: catalog X_IMAGE / Y_IMAGE (1-based); --center-x/-y are 1-based too.
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

from ogfkit import tsvio, imageio, cluster as cl, meta as ometa  # noqa: E402

COLUMNS = {
    'members': ['CL_RS_RES', 'CL_RS_PULL', 'CL_RS', 'CL_R', 'CL_PMEM', 'CL_MEMBER'],
    'arcs': ['ARC_LEN', 'ARC_WID', 'ARC_LW', 'ARC_RCURV', 'ARC_SPAN', 'ARC_ALIGN', 'ARC_CC', 'ARC_SNR', 'ARC_SCORE', 'ARC_FLAG'],
    'density': ['CL_SIGMA'],
}


def col(rows, name, default=float('nan')):
    return np.array([tsvio.fnum(r.get(name), default) for r in rows]) if name else np.full(len(rows), default)


def positions(rows):
    return col(rows, 'X_IMAGE') - 1.0, col(rows, 'Y_IMAGE') - 1.0


def extent(x, y, shape=None):
    if shape is not None:
        return shape
    ok = np.isfinite(x) & np.isfinite(y)
    if not ok.any():
        return (1, 1)
    return (int(math.ceil(np.nanmax(y[ok]))) + 1, int(math.ceil(np.nanmax(x[ok]))) + 1)


def image_shape_and_valid(path, want_valid):
    if not path or not os.path.isfile(path):
        return None, None
    data, _ = imageio.load_image(path)
    valid = (np.isfinite(data) & (data != 0)) if want_valid else None
    return data.shape, valid


def load_meta(path):
    return ometa.load(path).get('cluster', {}) if path else {}


def resolve_centre(a, x, y, shape, mags, mag_lo, mag_hi, valid=None, meta=None):
    """-> (cx, cy, how) in 0-based pixels."""
    if a.center_mode == 'manual' and a.center_x > 0 and a.center_y > 0:
        return a.center_x - 1.0, a.center_y - 1.0, 'manual'
    if meta and meta.get('centre') and a.task != 'members':
        return float(meta['centre'][0]), float(meta['centre'][1]), 'previous step'
    ok = np.isfinite(x) & np.isfinite(y) & np.isfinite(mags) & (mags >= mag_lo) & (mags <= mag_hi)
    sig = a.smooth_px if a.smooth_px > 0 else 0.03 * min(shape)
    if ok.sum() < 20:
        return 0.5 * shape[1], 0.5 * shape[0], 'image centre (too few objects)'
    dm = cl.density_map(x[ok], y[ok], shape, sig, binsize=max(2, int(round(sig / 4))), valid=valid)
    pk = cl.find_peaks(dm, 0.0, max_peaks=1)
    if not pk:
        return 0.5 * shape[1], 0.5 * shape[0], 'image centre (no density peak)'
    return pk[0]['x'], pk[0]['y'], 'density peak (%.1f sigma)' % pk[0]['sig']


def mag_window(a, mags, mask=None):
    m = mags if mask is None else mags[mask]
    m = m[np.isfinite(m)]
    lo, hi = a.mag_min, a.mag_max
    if lo == 0 and hi == 0 and m.size:
        lo, hi = float(np.percentile(m, 3)), float(np.percentile(m, 85))
    return lo, hi


def plot_cmd(path, mag, color, rs, flag, member, centre_info):
    try:
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
    except Exception:
        return
    fig, ax = plt.subplots(1, 1, figsize=(5.2, 4.0))
    ok = np.isfinite(mag) & np.isfinite(color)
    ax.plot(mag[ok], color[ok], '.', ms=2, color='0.6', label='all')
    k = ok & flag
    ax.plot(mag[k], color[k], 'o', ms=3, color='tab:red', label='red sequence')
    k = ok & member
    ax.plot(mag[k], color[k], 's', ms=5, mfc='none', mec='k', label='members')
    if rs.get('ok'):
        mm = np.linspace(np.nanmin(mag[ok]), np.nanmax(mag[ok]), 20)
        ln = rs['a'] + rs['b'] * (mm - rs['m0'])
        ax.plot(mm, ln, 'b-')
        ax.fill_between(mm, ln - 2 * rs['scatter'], ln + 2 * rs['scatter'], color='b', alpha=0.15)
    ax.set_xlabel('magnitude'); ax.set_ylabel('colour (blue - red)'); ax.legend(fontsize=7)
    ax.set_ylim(np.nanpercentile(color[ok], 1) - 0.2, np.nanpercentile(color[ok], 99) + 0.2)
    ax.set_title('colour = %.3f %+.3f (m - %.2f), scatter %.3f' % (rs['a'], rs['b'], rs['m0'], rs['scatter']) if rs.get('ok') else 'no red sequence', fontsize=8)
    fig.tight_layout(); fig.savefig(path, dpi=80); plt.close(fig)


def task_members(a, cols, rows, W):
    x, y = positions(rows)
    ng = len(rows)
    if a.mag_column not in cols or a.blue_column not in cols or a.red_column not in cols:
        raise SystemExit('cluster: missing column(s): need %s' % ', '.join(c for c in (a.mag_column, a.blue_column, a.red_column) if c not in cols))
    shape, valid = image_shape_and_valid(a.image, True)
    shape = extent(x, y, shape)
    mag = col(rows, a.mag_column)
    blue, red = col(rows, a.blue_column), col(rows, a.red_column)
    color = blue - red
    eb = col(rows, a.blue_err_column if a.blue_err_column in cols else '', 0.0)
    er = col(rows, a.red_err_column if a.red_err_column in cols else '', 0.0)
    err = np.sqrt(np.nan_to_num(eb) ** 2 + np.nan_to_num(er) ** 2 + a.min_color_err ** 2)
    lo, hi = mag_window(a, mag)
    cx, cy, how = resolve_centre(a, x, y, shape, mag, lo, hi)
    r = np.hypot(x - cx, y - cy)
    rfit = a.fit_radius if a.fit_radius > 0 else 0.12 * min(shape)
    sel = (r < rfit) & np.isfinite(color)
    rs = cl.red_sequence(color[sel], mag[sel], err[sel], mag_range=(lo, hi), fit_slope=not a.fixed_slope, slope_range=(a.slope_min, a.slope_max),
                         m0=a.mag_ref if a.mag_ref != 0 else None, window=a.rs_window)
    out = [dict() for _ in rows]
    summary = dict(centre=[cx, cy], centre_how=how, fit_radius=rfit, mag_range=[lo, hi], n_fit=int(sel.sum()), red_sequence=rs, columns=dict(mag=a.mag_column, blue=a.blue_column, red=a.red_column))
    flag = np.zeros(ng, bool); member = np.zeros(ng, bool)
    if rs['ok']:
        res, pull = cl.rs_pull(color, mag, rs, err)
        flag = np.isfinite(pull) & (np.abs(pull) < a.rs_nsig)
        pd, pw = cl.pixel_distances(shape, cx, cy, valid)
        inrange = (mag >= lo) & (mag <= hi)
        pmem, info = cl.radial_pmem(r, flag & inrange, pd, pw)
        member = flag & inrange & (pmem >= max(a.pmem_min, 1e-9))
        zinfo = None
        if a.z_column and a.z_column in cols:
            z = col(rows, a.z_column)
            zerr = None
            zc = a.z_center
            if zc <= 0:
                zz = z[flag & inrange & (r < rfit) & np.isfinite(z)]
                zc = float(np.median(zz)) if zz.size >= 3 else 0.0
            if zc > 0:
                inw, has = cl.redshift_window(z, zc, a.z_width, zerr)
                member = member & (~has | inw)
                zinfo = dict(zc=zc, width=a.z_width, n_with_z=int(has.sum()), n_rejected=int(np.sum(flag & inrange & has & ~inw)))
        if a.max_radius > 0:
            member &= r <= a.max_radius
        summary.update(n_rs=int(np.sum(flag & inrange)), n_members=int(member.sum()), radial_profile=info, z=zinfo,
                       expected_contamination=float(np.sum(1.0 - pmem[member])) if member.any() else 0.0)
        for i in range(ng):
            out[i] = dict(CL_RS_RES=res[i], CL_RS_PULL=pull[i], CL_RS=int(flag[i] and inrange[i]) if np.isfinite(pull[i]) else None,
                          CL_R=r[i], CL_PMEM=pmem[i], CL_MEMBER=int(member[i]))
    else:
        for i in range(ng):
            out[i] = dict(CL_R=r[i])
    plot_cmd(W('rs.png'), mag, color, rs, flag, member, how)
    return out, summary


def task_arcs(a, cols, rows, W):
    if not a.image or not os.path.isfile(a.image):
        raise SystemExit('cluster: arc detection needs the image')
    data, _ = imageio.load_image(a.image)
    data = np.asarray(data, np.float64)
    x, y = positions(rows)
    A = col(rows, 'A_IMAGE', 3.0); B = col(rows, 'B_IMAGE', 3.0)
    elong = A / np.maximum(B, 1e-3)
    from scipy import ndimage as ndi
    sm = ndi.gaussian_filter(np.nan_to_num(data), a.arc_smooth)
    sig_sm = imageio.robust_sigma(sm)
    meta = load_meta(a.meta_out)
    centre = None
    how = 'none'
    if a.center_mode == 'manual' and a.center_x > 0 and a.center_y > 0:
        centre = (a.center_x - 1.0, a.center_y - 1.0); how = 'manual'
    elif meta.get('centre'):
        centre = (float(meta['centre'][0]), float(meta['centre'][1])); how = 'previous cluster step'
    cand = np.isfinite(x) & np.isfinite(y) & (elong >= a.arc_elong_min) & (A * 2.0 >= a.arc_min_size)
    idx = np.nonzero(cand)[0]
    half = np.clip(4 * np.nan_to_num(A[idx], nan=3.0) + 10, 20, a.arc_max_cutout)
    search = np.maximum(4.0, 1.5 * np.nan_to_num(A[idx], nan=3.0))
    geo = cl.measure_arcs(data, x[idx], y[idx], sig_sm, centre=centre, half=half, search=search, smooth=a.arc_smooth, nsig=a.arc_nsig)
    out = [dict() for _ in rows]
    nflag = 0
    for i, g in zip(idx, geo):
        if not g.get('ok'):
            continue
        flag, score = cl.arc_score(g, a.arc_lw_min, a.arc_align_max, a.arc_cc_max, a.arc_snr_min, not a.arc_allow_straight, have_centre=centre is not None)
        nflag += flag
        out[i] = dict(ARC_LEN=g['length'], ARC_WID=g['width'], ARC_LW=g['lw'], ARC_RCURV=g['rcurv'] if g['curved'] else None, ARC_SPAN=g['span_deg'] if g['curved'] else None,
                      ARC_ALIGN=g.get('align_deg'), ARC_CC=g.get('cc_ratio') if g['curved'] else None, ARC_SNR=g['snr'], ARC_SCORE=score, ARC_FLAG=flag)
    summary = dict(arcs=dict(n_measured=int(sum(1 for g in geo if g.get('ok'))), n_candidates=int(nflag), prefilter_elongation=a.arc_elong_min, centre=list(centre) if centre else None,
                             centre_how=how, sigma_smoothed=float(sig_sm), lw_min=a.arc_lw_min, align_max=a.arc_align_max))
    return out, summary


def task_density(a, cols, rows, W):
    x, y = positions(rows)
    shape, valid = image_shape_and_valid(a.image, True)
    shape = extent(x, y, shape)
    mag = col(rows, a.mag_column) if a.mag_column in cols else np.full(len(rows), np.nan)
    lo, hi = mag_window(a, mag)
    use = np.isfinite(x) & np.isfinite(y)
    if np.isfinite(mag).any() and (a.mag_min != 0 or a.mag_max != 0):
        use &= (mag >= lo) & (mag <= hi)
    which = a.dens_sample
    if which in ('rs', 'members'):
        key = 'CL_RS' if which == 'rs' else 'CL_MEMBER'
        if key not in cols:
            raise SystemExit('cluster: sample "%s" needs the column %s (run the red-sequence step first)' % (which, key))
        use &= col(rows, key, 0.0) > 0
    sig = a.dens_sigma if a.dens_sigma > 0 else 0.03 * min(shape)
    b = max(1, int(a.dens_bin))
    if use.sum() < 10:
        raise SystemExit('cluster: only %d objects selected for the density map' % int(use.sum()))
    dm = cl.density_map(x[use], y[use], shape, sig, binsize=b, valid=valid)
    pk = cl.find_peaks(dm, a.peak_threshold)
    sigs = cl.sample_map(dm, x, y)
    out = [dict(CL_SIGMA=s) for s in sigs]
    from astropy.io import fits
    hdr = fits.Header()
    hdr['BUNIT'] = 'sigma'; hdr['CDELT1'] = float(b); hdr['CDELT2'] = float(b); hdr['CRPIX1'] = 0.5; hdr['CRPIX2'] = 0.5
    hdr['CRVAL1'] = 0.5 * b + 0.5; hdr['CRVAL2'] = 0.5 * b + 0.5
    hdr['CTYPE1'] = 'LINEAR'; hdr['CTYPE2'] = 'LINEAR'
    hdr['COMMENT'] = 'cluster_density: Poisson significance; map cell (i,j) covers pixels (i*b+1 .. (i+1)*b) of the catalog image; background %.4g per cell' % dm['lam']
    imageio.save_fits(W('density.fits'), dm['sig'].astype(np.float32), header=hdr)
    tsvio.write_table(W('peaks.tsv'), ['X_IMAGE', 'Y_IMAGE', 'SIGMA', 'EXCESS_CELLS'], [dict(X_IMAGE=p['x'] + 1, Y_IMAGE=p['y'] + 1, SIGMA=p['sig'], EXCESS_CELLS=p['excess_cells']) for p in pk])
    summary = dict(density=dict(n_used=int(use.sum()), sample=which, sigma_px=sig, bin=b, lam_per_cell=dm['lam'], peaks=pk[:10], threshold=a.peak_threshold, max_sigma=float(dm['sig'].max())))
    if pk:
        summary['density']['centre'] = [pk[0]['x'], pk[0]['y']]
    return out, summary


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--task', required=True, choices=['members', 'arcs', 'density'])
    ap.add_argument('--catalog', required=True)
    ap.add_argument('--image', default='')
    ap.add_argument('--work', default='.')
    ap.add_argument('--meta-out', default='')
    ap.add_argument('--mag-column', default='MAG_AUTO')
    ap.add_argument('--blue-column', default='MAG_AUTO')
    ap.add_argument('--red-column', default='MAG_ISOCOR')
    ap.add_argument('--blue-err-column', default='')
    ap.add_argument('--red-err-column', default='')
    ap.add_argument('--min-color-err', type=float, default=0.02)
    ap.add_argument('--mag-min', type=float, default=0.0)
    ap.add_argument('--mag-max', type=float, default=0.0)
    ap.add_argument('--mag-ref', type=float, default=0.0)
    ap.add_argument('--center-mode', default='auto', choices=['auto', 'manual'])
    ap.add_argument('--center-x', type=float, default=0.0)
    ap.add_argument('--center-y', type=float, default=0.0)
    ap.add_argument('--smooth-px', type=float, default=0.0)
    ap.add_argument('--fit-radius', type=float, default=0.0)
    ap.add_argument('--fixed-slope', action='store_true')
    ap.add_argument('--slope-min', type=float, default=-0.2)
    ap.add_argument('--slope-max', type=float, default=0.1)
    ap.add_argument('--rs-window', type=float, default=0.6)
    ap.add_argument('--rs-nsig', type=float, default=2.0)
    ap.add_argument('--z-column', default='')
    ap.add_argument('--z-center', type=float, default=0.0)
    ap.add_argument('--z-width', type=float, default=0.05)
    ap.add_argument('--max-radius', type=float, default=0.0)
    ap.add_argument('--pmem-min', type=float, default=0.5)
    ap.add_argument('--arc-nsig', type=float, default=2.0)
    ap.add_argument('--arc-smooth', type=float, default=1.0)
    ap.add_argument('--arc-elong-min', type=float, default=2.0)
    ap.add_argument('--arc-min-size', type=float, default=6.0)
    ap.add_argument('--arc-lw-min', type=float, default=4.0)
    ap.add_argument('--arc-align-max', type=float, default=30.0)
    ap.add_argument('--arc-cc-max', type=float, default=0.5)
    ap.add_argument('--arc-snr-min', type=float, default=5.0)
    ap.add_argument('--arc-allow-straight', action='store_true')
    ap.add_argument('--arc-max-cutout', type=float, default=300.0)
    ap.add_argument('--dens-sample', default='all', choices=['all', 'rs', 'members'])
    ap.add_argument('--dens-sigma', type=float, default=0.0)
    ap.add_argument('--dens-bin', type=int, default=4)
    ap.add_argument('--peak-threshold', type=float, default=3.0)
    a = ap.parse_args(argv)
    os.makedirs(a.work, exist_ok=True)
    W = lambda n: os.path.join(a.work, 'cluster_' + n)
    cols, rows = tsvio.read_catalog(a.catalog)
    fn = dict(members=task_members, arcs=task_arcs, density=task_density)[a.task]
    out, summary = fn(a, cols, rows, W)
    sp = W('summary.json')
    try:
        with open(sp) as fh:
            old = json.load(fh)
    except (OSError, ValueError):
        old = {}
    old.update(summary)
    with open(sp, 'w') as fh:
        json.dump(old, fh, indent=1, default=lambda o: None)
    if a.meta_out:
        m = dict(load_meta(a.meta_out))
        if a.task == 'members' and summary.get('red_sequence', {}).get('ok'):
            rs = summary['red_sequence']
            m.update(centre=summary['centre'], rs_a=rs['a'], rs_b=rs['b'], rs_m0=rs['m0'], rs_scatter=rs['scatter'], n_members=summary.get('n_members'))
        elif a.task == 'members':
            m.update(centre=summary['centre'])
        if a.task == 'density' and summary['density'].get('centre'):
            m['density_peak_sigma'] = summary['density']['peaks'][0]['sig']
            m.setdefault('centre', summary['density']['centre'])
        if a.task == 'arcs':
            m['n_arc_candidates'] = summary['arcs']['n_candidates']
        ometa.update(a.meta_out, 'cluster', m, nrows=len(rows))
    cs = COLUMNS[a.task]
    tsvio.write_columns(cs, [(r['NUMBER'], d) for r, d in zip(rows, out)])
    sys.stderr.write('cluster %s: %s\n' % (a.task, json.dumps({k: v for k, v in summary.items() if k != 'red_sequence'}, default=lambda o: None)[:400]))
    if a.task == 'members' and summary.get('red_sequence'):
        rs = summary['red_sequence']
        sys.stderr.write('red sequence: ' + ('colour = %.3f %+.4f (m - %.2f), intrinsic scatter %.3f, %d in sequence, %d members\n' % (rs['a'], rs['b'], rs['m0'], rs['scatter'], summary.get('n_rs', 0), summary.get('n_members', 0)) if rs['ok'] else rs['note'] + '\n'))
    return 0


if __name__ == '__main__':
    sys.exit(main())
