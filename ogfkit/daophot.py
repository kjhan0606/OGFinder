"""DAOPHOT / ALLSTAR / DAOMASTER-style crowded-field photometry - pure functions, JSON-serialisable results.

    find        FIND      matched-filter detection with sharpness / roundness cuts            (photutils DAOStarFinder)
    phot        PHOT      multi-aperture photometry with a sky annulus (median / mode / clipped mean / global map)
    pick_psf    PICKPSF   automatic PSF-star selection with neighbour / quality rejection + manual add / remove
    (psf)       PSF       ogfkit.psfmodel.build_psf_model: analytic (gaussian, moffat, lorentz, penny) + lookup table, order 0-3,
                          iterative neighbour subtraction
    allstar     NSTAR / ALLSTAR  simultaneous group fits, iterative re-find on the residual, sky re-fit, chi / sharp / error / flags
    substar / addstar   SUBSTAR / ADDSTAR
    match_lists / master_list   DAOMATCH / DAOMASTER: triangle matching + polynomial coordinate transforms, master list
    aperture_correction   growth curve of the PSF stars (neighbours subtracted) -> PSF magnitude -> aperture magnitude
    DaophotDetector   detection callback for the completeness (artificial-star) tests

All positions are 0-based pixel centres inside the functions; the CLI prints 1-based X/Y like the catalog.
Quality flags (bit field) - FLAG_*: see the constants below.
"""
import math
import warnings

import numpy as np
from scipy import ndimage as ndi
from scipy.spatial import cKDTree

from . import psfmodel as pm

FLAG_NOCONV = 1      # fit not converged
FLAG_CHI = 2         # chi above the limit
FLAG_SHARP = 4       # sharp outside the range
FLAG_SNR = 8         # flux / error below the limit
FLAG_EDGE = 16       # fit radius reaches the image edge
FLAG_GROUP = 32      # group larger than max_group (fitted in chunks)
FLAG_NEG = 64        # flux <= 0 or no valid error
FLAG_NEW = 128       # found on the residual image (not in the first list)
FLAG_ROUND = 256     # round outside the range
BAD_MASK = FLAG_NOCONV | FLAG_CHI | FLAG_SHARP | FLAG_SNR | FLAG_NEG | FLAG_ROUND


# ------------------------------------------------------------------------------------------------------------------ FIND
def find(data, bkg, rms, fwhm=3.0, thresh=4.0, sharp=(0.2, 1.0), round_=(-1.0, 1.0), mask=None, max_stars=None):
    """DAOFIND: Gaussian-kernel matched filter, local maxima above thresh * sigma, DAOFIND sharpness / roundness cuts.

    Returns dict of arrays x, y (0-based), peak, flux, sharp, round1, round2, npix."""
    from photutils.detection import DAOStarFinder
    warnings.filterwarnings('ignore', module='photutils')
    img = np.nan_to_num(np.asarray(data, float) - bkg, nan=0.0)
    sig = float(np.median(np.asarray(rms)[np.isfinite(rms)]))
    try:
        f = DAOStarFinder(threshold=thresh * sig, fwhm=fwhm, sharpness_range=tuple(sharp), roundness_range=tuple(round_), n_brightest=max_stars)
    except TypeError:                         # photutils < 3
        f = DAOStarFinder(threshold=thresh * sig, fwhm=fwhm, sharplo=sharp[0], sharphi=sharp[1], roundlo=round_[0], roundhi=round_[1], brightest=max_stars)
    bad = ~np.isfinite(data)
    if mask is not None:
        bad = bad | mask
    t = f(img, mask=bad)
    if t is None or len(t) == 0:
        z = np.zeros(0)
        return dict(x=z, y=z, peak=z, flux=z, sharp=z, round1=z, round2=z, npix=z)
    cn = lambda a, b: a if a in t.colnames else b
    return dict(x=np.asarray(t[cn('x_centroid', 'xcentroid')], float), y=np.asarray(t[cn('y_centroid', 'ycentroid')], float), peak=np.asarray(t['peak'], float),
                flux=np.asarray(t['flux'], float), sharp=np.asarray(t['sharpness'], float), round1=np.asarray(t['roundness1'], float),
                round2=np.asarray(t['roundness2'], float), npix=np.asarray(t[cn('n_pixels', 'npix')], float))


# ------------------------------------------------------------------------------------------------------------------ PHOT
SKY_MODES = ('median', 'mode', 'mean', 'global')


def _sky_stats(vals, mode):
    from astropy.stats import sigma_clipped_stats
    if vals.size < 8:
        return float('nan'), float('nan'), int(vals.size)
    mean, med, sd = sigma_clipped_stats(vals, sigma=3.0, maxiters=5)
    if mode == 'mean':
        v = mean
    elif mode == 'mode':
        v = 3.0 * med - 2.0 * mean              # DAOPHOT / SExtractor 'mode' estimate for skewed (crowded) annuli
    else:
        v = med
    return float(v), float(sd), int(vals.size)


def phot(data, xy, radii=(3.0, 5.0, 8.0), sky_inner=12.0, sky_outer=18.0, sky_mode='median', bkg=None, rms=None, zp=25.0, gain=None):
    """PHOT: circular-aperture photometry (exact pixel overlap) with a local sky annulus.

    Returns a list of dicts: sky, skysig, nsky, flux[], fluxerr[], mag[], magerr[] (one entry per aperture)."""
    from photutils.aperture import CircularAperture, CircularAnnulus
    d = np.asarray(data, float)
    xy = np.asarray(xy, float).reshape(-1, 2)
    out = []
    ann = CircularAnnulus(xy, sky_inner, sky_outer) if len(xy) else None
    ann_masks = ann.to_mask(method='center') if ann is not None else []
    ann_masks = ann_masks if isinstance(ann_masks, list) else [ann_masks]
    aps = [CircularAperture(xy, r).to_mask(method='exact') for r in radii] if len(xy) else []
    aps = [a if isinstance(a, list) else [a] for a in aps]
    for i in range(len(xy)):
        sky = sd = float('nan')
        nsky = 0
        if sky_mode == 'global' and bkg is not None:
            ix, iy = int(round(xy[i][0])), int(round(xy[i][1]))
            if 0 <= iy < d.shape[0] and 0 <= ix < d.shape[1]:
                sky = float(bkg[iy, ix])
                sd = float(rms[iy, ix]) if rms is not None else float('nan')
        else:
            cut = ann_masks[i].cutout(d, fill_value=np.nan)
            w = ann_masks[i].data
            if cut is not None:
                v = cut[(w > 0) & np.isfinite(cut)]
                sky, sd, nsky = _sky_stats(v, sky_mode)
        fl, fe, mg, me = [], [], [], []
        for r, mkl in zip(radii, aps):
            mk = mkl[i]
            cut = mk.cutout(d, fill_value=np.nan)
            if cut is None or not np.isfinite(sky):
                fl.append(float('nan')); fe.append(float('nan')); mg.append(float('nan')); me.append(float('nan'))
                continue
            w = mk.data
            ok = np.isfinite(cut)
            area = float((w * ok).sum())
            f = float(((cut - sky) * w)[ok].sum())
            var = area * (sd ** 2) * (1.0 + (area / max(nsky, 1)) if nsky else 1.0)
            if gain:
                var += max(f, 0.0) / gain
            e = math.sqrt(max(var, 0.0))
            fl.append(f); fe.append(e)
            mg.append(-2.5 * math.log10(f) + zp if f > 0 else float('nan'))
            me.append(1.0857 * e / f if f > 0 else float('nan'))
        out.append(dict(sky=sky, skysig=sd, nsky=nsky, flux=fl, fluxerr=fe, mag=mg, magerr=me))
    return out


# --------------------------------------------------------------------------------------------------------------- PICKPSF
def pick_psf(stars, shape, n_max=60, half=10, snr_min=15.0, saturation=None, sharp=(0.3, 0.95), round_=(-0.5, 0.5), neighbour_radius=None,
             neighbour_ratio=0.1, add=None, remove=None, add_tol=3.0):
    """PICKPSF: the brightest stars that are unsaturated, not at the edge, pass the sharp / round cuts and have no neighbour
    brighter than `neighbour_ratio` x their flux inside `neighbour_radius` (default 1.5 x half).

    stars: dict of arrays x, y, flux, peak, sharp, round1 (output of find()).  add / remove: lists of (x, y) 0-based positions
    (matched to the nearest star within add_tol pixels; an `add` position without a star is used as it is).  Returns (positions
    array (N,2), info dict with the per-star reason of every rejection)."""
    x, y, fl = np.asarray(stars['x']), np.asarray(stars['y']), np.asarray(stars['flux'])
    n = len(x)
    rad = neighbour_radius or 1.5 * half
    why = np.array([''] * n, dtype=object)
    ok = np.ones(n, bool)
    def rej(mask, text):
        nonlocal ok
        new = ok & mask
        why[new] = text
        ok = ok & ~mask
    snr_ok = np.asarray(stars.get('snr', np.full(n, np.inf))) >= snr_min
    rej(~snr_ok, 'snr')
    rej((x < half) | (y < half) | (x > shape[1] - 1 - half) | (y > shape[0] - 1 - half), 'edge')
    if saturation is not None:
        rej(np.asarray(stars['peak']) >= saturation, 'saturated')
    rej((np.asarray(stars['sharp']) < sharp[0]) | (np.asarray(stars['sharp']) > sharp[1]), 'sharp')
    rej((np.asarray(stars['round1']) < round_[0]) | (np.asarray(stars['round1']) > round_[1]), 'round')
    if n > 1:
        tr = cKDTree(np.c_[x, y])
        for i in np.where(ok)[0]:
            for j in tr.query_ball_point([x[i], y[i]], rad):
                if j != i and fl[j] > neighbour_ratio * fl[i]:
                    ok[i] = False
                    why[i] = 'neighbour'
                    break
    order = np.where(ok)[0]
    order = order[np.argsort(-fl[order])][:n_max]
    sel = [i for i in order]
    def nearest(p):
        if n == 0:
            return None, 1e9
        d = np.hypot(x - p[0], y - p[1])
        j = int(np.argmin(d))
        return j, float(d[j])
    removed = []
    for p in (remove or []):
        j, d = nearest(p)
        if j is not None and d <= add_tol and j in sel:
            sel.remove(j)
            removed.append(int(j))
    pos = [(float(x[i]), float(y[i])) for i in sel]
    added = []
    for p in (add or []):
        j, d = nearest(p)
        q = (float(x[j]), float(y[j])) if (j is not None and d <= add_tol) else (float(p[0]), float(p[1]))
        if all(math.hypot(q[0] - a, q[1] - b) > 1.0 for a, b in pos):
            pos.append(q)
            added.append(q)
    return np.array(pos, float).reshape(-1, 2), dict(n_input=int(n), n_auto=int(len(order)), n_added=len(added), n_removed=len(removed),
                                                    reasons={k: int((why == k).sum()) for k in ('snr', 'edge', 'saturated', 'sharp', 'round', 'neighbour')})


# ------------------------------------------------------------------------------------------------------------------ fitting
_SH = {}


def render(shape, model, xs, ys, fs):
    """Image of the PSF-model stars (float32)."""
    img = np.zeros(shape, np.float32)
    for x, y, f in zip(xs, ys, fs):
        if f > 0 and np.isfinite(f):
            st, x0, y0 = model.stamp_at(x, y)
            pm._paste(img, (st * f).astype(np.float32), x0, y0)
    return img


def _region(shape, x, y, pad):
    return (max(0, int(math.floor(min(x) - pad))), max(0, int(math.floor(min(y) - pad))),
            min(shape[1], int(math.ceil(max(x) + pad)) + 1), min(shape[0], int(math.ceil(max(y) + pad)) + 1))


def _stamp_region(model, x, y, bb):
    x0, y0, x1, y1 = bb
    out = np.zeros((y1 - y0, x1 - x0))
    st, sx, sy = model.stamp_at(x, y)
    h, w = st.shape
    ya, xa = max(y0, sy), max(x0, sx)
    yb, xb = min(y1, sy + h), min(x1, sx + w)
    if yb > ya and xb > xa:
        out[ya - y0:yb - y0, xa - x0:xb - x0] = st[ya - sy:yb - sy, xa - sx:xb - sx]
    return out


def fit_group(resid, rms, model, group, R, maxshift=2.0, fit_sky=True, gain=None, maxit=20, mask=None):
    """NSTAR: simultaneous (flux, x, y [, common sky]) fit of the stars of one group (list of dicts x, y, flux).

    `resid` is the image with ALL current stars subtracted; the group's own light is added back inside.  Levenberg-Marquardt
    on the pixels within R of any group member, weights 1/rms (plus Poisson with `gain`).  Returns the list of updated dicts
    (x, y, flux, fluxerr, xerr, yerr, chi, sharp, sky, npix, nit, conv, edge)."""
    n = len(group)
    xs = np.array([g['x'] for g in group], float)
    ys = np.array([g['y'] for g in group], float)
    fs = np.array([max(g['flux'], 1e-6) for g in group], float)
    pad = R + model.size // 2 + 1
    bb = _region(resid.shape, xs, ys, pad)
    x0, y0, x1, y1 = bb
    yy, xx = np.mgrid[y0:y1, x0:x1]
    sel = np.zeros(yy.shape, bool)
    for x, y in zip(xs, ys):
        sel |= np.hypot(xx - x, yy - y) <= R
    sig = np.asarray(rms[y0:y1, x0:x1], float)
    ok = sel & np.isfinite(resid[y0:y1, x0:x1]) & np.isfinite(sig) & (sig > 0)
    if mask is not None:
        ok &= ~mask[y0:y1, x0:x1]
    npix = int(ok.sum())
    edge = any(x - R < 0 or y - R < 0 or x + R > resid.shape[1] - 1 or y + R > resid.shape[0] - 1 for x, y in zip(xs, ys))
    npar = 3 * n + (1 if fit_sky else 0)
    if npix < npar + 2:
        return [dict(x=float(xs[i]), y=float(ys[i]), flux=float(fs[i]), fluxerr=float('nan'), xerr=float('nan'), yerr=float('nan'), chi=float('nan'),
                     sharp=float('nan'), sky=0.0, npix=npix, nit=0, conv=False, edge=edge) for i in range(n)]
    x_init, y_init = xs.copy(), ys.copy()
    base = resid[y0:y1, x0:x1].astype(float).copy()
    for i in range(n):                      # add the group's own light back
        base += fs[i] * _stamp_region(model, xs[i], ys[i], bb)
    d = base[ok]
    sky = 0.0
    lam = 1e-3
    h = 0.05
    def build(xs, ys):
        P = np.array([_stamp_region(model, xs[i], ys[i], bb)[ok] for i in range(n)])
        dPx = np.array([(_stamp_region(model, xs[i] + h, ys[i], bb)[ok] - _stamp_region(model, xs[i] - h, ys[i], bb)[ok]) / (2 * h) for i in range(n)])
        dPy = np.array([(_stamp_region(model, xs[i], ys[i] + h, bb)[ok] - _stamp_region(model, xs[i], ys[i] - h, bb)[ok]) / (2 * h) for i in range(n)])
        return P, dPx, dPy
    def sigma_of(mod):
        s2 = sig[ok] ** 2
        if gain:
            s2 = s2 + np.clip(mod, 0, None) / gain
        return np.sqrt(s2)
    P, dPx, dPy = build(xs, ys)
    mod = fs @ P + sky
    s = sigma_of(mod)
    chi2 = float((((d - mod) / s) ** 2).sum())
    conv = False
    nit = 0
    for nit in range(1, maxit + 1):
        J = np.zeros((d.size, npar))
        for i in range(n):
            J[:, 3 * i] = P[i] / s
            J[:, 3 * i + 1] = fs[i] * dPx[i] / s
            J[:, 3 * i + 2] = fs[i] * dPy[i] / s
        if fit_sky:
            J[:, -1] = 1.0 / s
        r = (d - mod) / s
        A = J.T @ J
        g = J.T @ r
        improved = False
        for _try in range(8):
            M = A + lam * np.diag(np.diag(A) + 1e-12)
            try:
                dp = np.linalg.solve(M, g)
            except np.linalg.LinAlgError:
                lam *= 10
                continue
            nxs, nys, nfs = xs.copy(), ys.copy(), fs.copy()
            for i in range(n):
                nfs[i] = max(fs[i] + dp[3 * i], 1e-6)
                nxs[i] = np.clip(xs[i] + np.clip(dp[3 * i + 1], -0.7, 0.7), x_init[i] - maxshift, x_init[i] + maxshift)
                nys[i] = np.clip(ys[i] + np.clip(dp[3 * i + 2], -0.7, 0.7), y_init[i] - maxshift, y_init[i] + maxshift)
            nsky = sky + (dp[-1] if fit_sky else 0.0)
            Pn, dPxn, dPyn = build(nxs, nys)
            modn = nfs @ Pn + nsky
            sn = sigma_of(modn)
            chi2n = float((((d - modn) / sn) ** 2).sum())
            if chi2n <= chi2 or lam > 1e4:
                improved = chi2n <= chi2
                step_small = (np.abs(nxs - xs).max() < 0.005 and np.abs(nys - ys).max() < 0.005 and np.all(np.abs(nfs - fs) < 2e-3 * np.maximum(fs, 1e-3) + 1e-9))
                xs, ys, fs, sky, P, dPx, dPy, mod, s, chi2 = nxs, nys, nfs, nsky, Pn, dPxn, dPyn, modn, sn, chi2n
                lam = max(lam / 5.0, 1e-7)
                if step_small or not improved:
                    conv = True
                break
            lam *= 8.0
        else:
            conv = False
            break
        if conv:
            break
    # covariance and diagnostics
    J = np.zeros((d.size, npar))
    for i in range(n):
        J[:, 3 * i] = P[i] / s
        J[:, 3 * i + 1] = fs[i] * dPx[i] / s
        J[:, 3 * i + 2] = fs[i] * dPy[i] / s
    if fit_sky:
        J[:, -1] = 1.0 / s
    dof = max(npix - npar, 1)
    chi_red = chi2 / dof
    try:
        cov = np.linalg.pinv(J.T @ J) * max(chi_red, 1.0)
        err = np.sqrt(np.clip(np.diag(cov), 0, None))
    except np.linalg.LinAlgError:
        err = np.full(npar, float('nan'))
    res_img = np.zeros(yy.shape)
    res_img[ok] = d - mod
    out = []
    for i in range(n):
        Pi = _stamp_region(model, xs[i], ys[i], bb)
        core = (np.hypot(xx - xs[i], yy - ys[i]) <= 0.7 * max(model.fwhm_estimate(), 2.0)) & ok
        pk = Pi.max()
        sh = float((res_img[core] * Pi[core]).sum() / max((fs[i] * (Pi[core] ** 2)).sum(), 1e-12)) if core.any() else float('nan')
        out.append(dict(x=float(xs[i]), y=float(ys[i]), flux=float(fs[i]), fluxerr=float(err[3 * i]), xerr=float(err[3 * i + 1]), yerr=float(err[3 * i + 2]),
                        chi=float(math.sqrt(chi_red)), sharp=sh, sky=float(sky), npix=npix, nit=nit, conv=bool(conv), edge=edge))
    return out


def _init_worker(shared):
    # spawn start method (Windows/macOS): the globals of the parent are not inherited, so they are passed through the initializer
    _SH.update(shared)


def _worker(chunk):
    S = _SH
    res = []
    for grp in chunk:
        res.append(fit_group(S['resid'], S['rms'], S['model'], grp, S['R'], S['maxshift'], S['fit_sky'], S['gain'], mask=S['mask']))
    return res


def make_groups(stars, crit, max_group=20):
    """Connected components of the pair graph (distance < crit); components larger than max_group are cut in BFS order
    (returned as separate groups; their members get FLAG_GROUP).  Returns list of index lists, set of chunked indices."""
    n = len(stars)
    if n == 0:
        return [], set()
    pts = np.array([[s['x'], s['y']] for s in stars])
    tr = cKDTree(pts)
    nb = [set() for _ in range(n)]
    for i, j in tr.query_pairs(crit):
        nb[i].add(j)
        nb[j].add(i)
    seen = np.zeros(n, bool)
    groups, chunked = [], set()
    for s0 in range(n):
        if seen[s0]:
            continue
        comp, q = [], [s0]
        seen[s0] = True
        while q:
            k = q.pop(0)
            comp.append(k)
            for m in sorted(nb[k]):
                if not seen[m]:
                    seen[m] = True
                    q.append(m)
        if len(comp) <= max_group:
            groups.append(comp)
        else:
            chunked.update(comp)
            for a in range(0, len(comp), max_group):
                groups.append(comp[a:a + max_group])
    return groups, chunked


def _fit_all(stars, resid, rms, model, R, maxshift, fit_sky, gain, max_group, n_workers, mask):
    groups, chunked = make_groups(stars, 2.0 * R, max_group)
    _SH.update(resid=resid, rms=rms, model=model, R=R, maxshift=maxshift, fit_sky=fit_sky, gain=gain, mask=mask)
    work = [[[stars[i] for i in g]] for g in groups]
    nw = max(1, int(n_workers))
    if nw > 1 and len(groups) >= 8:
        import multiprocessing as mp
        chunks = [[[stars[i] for i in g] for g in groups[k::nw]] for k in range(nw)]
        ctx = mp.get_context('fork' if 'fork' in mp.get_all_start_methods() else 'spawn')
        with ctx.Pool(nw, initializer=_init_worker, initargs=(dict(_SH),)) as pool:
            parts = pool.map(_worker, chunks)
        results = {}
        for k in range(nw):
            for g, r in zip(groups[k::nw], parts[k]):
                results[tuple(g)] = r
        fitted = [results[tuple(g)] for g in groups]
    else:
        fitted = [fit_group(resid, rms, model, [stars[i] for i in g], R, maxshift, fit_sky, gain, mask=mask) for g in groups]
    for g, r in zip(groups, fitted):
        for i, rr in zip(g, r):
            stars[i].update(rr)
            stars[i]['gsize'] = len(g)
            stars[i]['chunked'] = i in chunked
    return stars


def _flag(star, zp, cuts):
    f = star.get('flag0', 0)
    if not star.get('conv', True):
        f |= FLAG_NOCONV
    fl, fe = star['flux'], star.get('fluxerr', float('nan'))
    if not (np.isfinite(fl) and fl > 1e-5 and np.isfinite(fe) and fe > 0):
        f |= FLAG_NEG
    else:
        if fl / fe < cuts['snr_min']:
            f |= FLAG_SNR
    if np.isfinite(star.get('chi', float('nan'))) and star['chi'] > cuts['chi_max']:
        f |= FLAG_CHI
    sh = star.get('sharp', float('nan'))
    if np.isfinite(sh) and not (cuts['sharp'][0] <= sh <= cuts['sharp'][1]):
        f |= FLAG_SHARP
    rd = star.get('round', float('nan'))
    if np.isfinite(rd) and not (cuts['round'][0] <= rd <= cuts['round'][1]):
        f |= FLAG_ROUND
    if star.get('edge'):
        f |= FLAG_EDGE
    if star.get('chunked'):
        f |= FLAG_GROUP
    if star.get('found_iter', 0) > 0:
        f |= FLAG_NEW
    return f


def allstar(data, model, bkg=None, rms=None, mask=None, init_xy=None, init_stars=None, fwhm=None, thresh=4.0, n_iter=3, fit_radius=None, fit_sky=True,
            gain=None, zp=25.0, max_group=20, snr_min=3.0, chi_max=3.0, sharp_cut=(-1.0, 1.0), round_cut=(-1.0, 1.0), find_sharp=(0.2, 1.0),
            find_round=(-1.0, 1.0), min_sep=None, maxshift=2.0, n_workers=1, refind=True):
    """ALLSTAR: iterate find (on the residual) -> group fit -> subtract.  Returns (result dict, residual image float32).

    result['stars'] = list of dicts: id, x, y (0-based), flux, fluxerr, mag, magerr, chi, sharp, round, snr, sky, gsize, nit, found_iter, flag, good."""
    data = np.asarray(data, np.float32)
    if bkg is None or rms is None:
        bkg, rms = pm.background(data, mask)
    sub = (np.nan_to_num(data, nan=0.0) - bkg).astype(np.float32)
    fw = float(fwhm or model.fwhm_estimate())
    R = float(fit_radius or max(2.5, 1.5 * fw))
    min_sep = float(min_sep or 0.9 * fw)
    cuts = dict(snr_min=snr_min, chi_max=chi_max, sharp=sharp_cut, round=round_cut)
    stars = []
    def add(x, y, fl, it, rd=float('nan')):
        stars.append(dict(x=float(x), y=float(y), flux=float(max(fl, 1e-3)), found_iter=it, **{'round': rd}, flux0=float(fl)))
    if init_stars is not None:
        for x, y, fl, rd in zip(init_stars['x'], init_stars['y'], init_stars['flux'], init_stars.get('round1', np.full(len(init_stars['x']), np.nan))):
            add(x, y, fl, 0, rd)
    elif init_xy is not None:
        pk = model.stamp(0.5 * sum(model.xr), 0.5 * sum(model.yr)).max()
        for x, y in np.asarray(init_xy, float).reshape(-1, 2):
            ix, iy = int(round(x)), int(round(y))
            v = sub[iy, ix] if (0 <= iy < sub.shape[0] and 0 <= ix < sub.shape[1]) else 1.0
            add(x, y, max(v, 1.0) / pk, 0)
    resid = sub.copy()
    nfound = []
    for it in range(max(1, n_iter)):
        if refind or it == 0 and init_xy is None and init_stars is None:
            f = find(resid + 0.0, np.zeros_like(resid), rms, fwhm=fw, thresh=thresh, sharp=find_sharp, round_=find_round, mask=mask)
            pk = model.stamp(0.5 * sum(model.xr), 0.5 * sum(model.yr)).max()
            if len(stars):
                tr = cKDTree(np.array([[s['x'], s['y']] for s in stars]))
                d, _ = tr.query(np.c_[f['x'], f['y']]) if len(f['x']) else (np.zeros(0), None)
            else:
                d = np.full(len(f['x']), 1e9)
            nn = 0
            for k in range(len(f['x'])):
                if d[k] > min_sep:
                    add(f['x'][k], f['y'][k], max(f['peak'][k], 0.1) / pk, it, f['round1'][k])
                    nn += 1
            nfound.append(nn)
        if not stars:
            break
        _fit_all(stars, resid, rms, model, R, maxshift, fit_sky, gain, max_group, n_workers, mask)
        # drop stars that went negative / invisible so that they do not subtract light
        stars = [s for s in stars if s['flux'] > 1e-5 and np.isfinite(s.get('fluxerr', np.nan)) and s['flux'] / max(s['fluxerr'], 1e-12) >= 0.5 * snr_min]
        resid = sub - render(sub.shape, model, [s['x'] for s in stars], [s['y'] for s in stars], [s['flux'] for s in stars])
    if stars:                                                           # final settling pass (no new stars)
        _fit_all(stars, resid, rms, model, R, maxshift, fit_sky, gain, max_group, n_workers, mask)
        resid = sub - render(sub.shape, model, [s['x'] for s in stars], [s['y'] for s in stars], [s['flux'] for s in stars])
    out = []
    for k, s in enumerate(stars):
        fl, fe = s['flux'], s.get('fluxerr', float('nan'))
        s['flag'] = _flag(s, zp, cuts)
        out.append(dict(id=k + 1, x=s['x'], y=s['y'], flux=fl, fluxerr=fe, mag=-2.5 * math.log10(fl) + zp if fl > 0 else float('nan'),
                        magerr=1.0857 * fe / fl if (fl > 0 and np.isfinite(fe)) else float('nan'), chi=s.get('chi'), sharp=s.get('sharp'),
                        round=s.get('round'), snr=fl / fe if (np.isfinite(fe) and fe > 0) else float('nan'), sky=s.get('sky', 0.0) , gsize=s.get('gsize', 1),
                        nit=s.get('nit', 0), found_iter=s['found_iter'], flag=int(s['flag']), good=bool((s['flag'] & BAD_MASK) == 0),
                        xerr=s.get('xerr'), yerr=s.get('yerr')))
    return dict(stars=out, n_found_per_iter=nfound, fit_radius=R, fwhm=fw, n=len(out), n_good=int(sum(o['good'] for o in out))), resid.astype(np.float32)


# ------------------------------------------------------------------------------------------------------- SUBSTAR / ADDSTAR
def substar(data, model, stars, keep=None, bkg=None):
    """SUBSTAR: subtract the fitted stars (dicts with x, y, flux) from `data`; `keep` = ids that stay in the image."""
    keep = set(keep or [])
    sel = [s for s in stars if s.get('id') not in keep]
    return (np.asarray(data, np.float32) - render(data.shape, model, [s['x'] for s in sel], [s['y'] for s in sel], [s['flux'] for s in sel])).astype(np.float32)


def addstar(data, model, xy, mags, zp=25.0, gain=None, seed=1):
    """ADDSTAR: add artificial PSF stars (0-based positions, magnitudes); Poisson noise of the added light when `gain` is given.
    Returns (image float32, truth list)."""
    rng = np.random.default_rng(seed)
    fl = 10.0 ** (-0.4 * (np.asarray(mags, float) - zp))
    add = render(data.shape, model, [p[0] for p in xy], [p[1] for p in xy], fl).astype(float)
    if gain:
        add = rng.poisson(np.clip(add, 0, None) * gain) / gain
    truth = [dict(x=float(p[0]), y=float(p[1]), mag=float(m), flux=float(f)) for p, m, f in zip(xy, mags, fl)]
    return (np.asarray(data, np.float32) + add.astype(np.float32)), truth


def random_positions(shape, n, margin, seed=1, avoid=None, min_dist=0.0):
    rng = np.random.default_rng(seed)
    pts = []
    tr = cKDTree(np.asarray(avoid, float)) if avoid is not None and len(avoid) else None
    guard = 0
    while len(pts) < n and guard < 100 * n:
        guard += 1
        p = (rng.uniform(margin, shape[1] - 1 - margin), rng.uniform(margin, shape[0] - 1 - margin))
        if tr is not None and tr.query(p)[0] < min_dist:
            continue
        pts.append(p)
    return pts


# ------------------------------------------------------------------------------------------------------ aperture correction
def aperture_correction(sub, model, stars, psf_xy, radii=(3, 4, 5, 6, 8, 10, 12), sky_mode='median', zp=25.0, tol=2.0, ref_radius=None):
    """Growth curve from the PSF stars after subtraction of all other fitted stars, and the correction PSF-mag -> aperture-mag.

    corr(R) = median( m_ap(R) - m_psf ) over the PSF stars (neighbours subtracted, sky 0 because `sub` is sky-subtracted);
    adding corr(R_ref) to a PSF magnitude gives the magnitude inside R_ref.  Returns dict(radii, corr, scatter, n, model_growth)."""
    from photutils.aperture import CircularAperture, aperture_photometry
    xs = np.array([s['x'] for s in stars]); ys = np.array([s['y'] for s in stars])
    if len(xs) == 0 or len(psf_xy) == 0:
        return dict(radii=list(radii), corr=[float('nan')] * len(radii), scatter=[float('nan')] * len(radii), n=0)
    tr = cKDTree(np.c_[xs, ys])
    allmod = render(sub.shape, model, xs, ys, [s['flux'] for s in stars])
    per = []
    for p in psf_xy:
        d, j = tr.query(p)
        if d > tol:
            continue
        s = stars[j]
        if not (s['flux'] > 0):
            continue
        own = np.zeros_like(allmod)
        st, x0, y0 = model.stamp_at(s['x'], s['y'])
        pm._paste(own, (st * s['flux']).astype(np.float32), x0, y0)
        clean = sub - (allmod - own)                       # only this star left
        ap = [CircularAperture([(s['x'], s['y'])], r) for r in radii]
        t = [float(aperture_photometry(clean, a, method='exact')['aperture_sum'][0]) for a in ap]
        per.append([(-2.5 * math.log10(f / s['flux'])) if f > 0 else float('nan') for f in t])
    A = np.array(per).reshape(-1, len(radii))
    corr, scat = [], []
    for k in range(len(radii)):
        v = A[:, k][np.isfinite(A[:, k])]
        if len(v) >= 3:
            q = v[np.abs(v - np.median(v)) < 3 * 1.4826 * np.median(np.abs(v - np.median(v))) + 1e-6]
            corr.append(float(np.median(q))); scat.append(float(np.std(q) / math.sqrt(max(len(q), 1))))
        else:
            corr.append(float(np.median(v)) if len(v) else float('nan')); scat.append(float('nan'))
    rr, gm = pm.growth(model, radii=radii)
    return dict(radii=[float(r) for r in radii], corr=corr, err=scat, n=int(len(A)), model_growth=gm, model_corr=[float(-2.5 * math.log10(g)) if g > 0 else float('nan') for g in gm])


# ------------------------------------------------------------------------------------------------------ DAOMATCH / DAOMASTER
def _tri_features(pts, idx):
    out, tri = [], []
    n = len(idx)
    for a in range(n):
        for b in range(a + 1, n):
            for c in range(b + 1, n):
                P = pts[[idx[a], idx[b], idx[c]]]
                s = sorted([np.hypot(*(P[0] - P[1])), np.hypot(*(P[1] - P[2])), np.hypot(*(P[0] - P[2]))])
                if s[2] < 1e-6 or s[0] < 0.1 * s[2]:
                    continue
                out.append((s[1] / s[2], s[0] / s[2]))
                tri.append((idx[a], idx[b], idx[c]))
    return np.array(out), tri


def _similarity(P, Q):
    """Least-squares similarity (rotation + scale + shift) mapping P -> Q, as 2x3 matrix."""
    mp, mq = P.mean(0), Q.mean(0)
    p, q = P - mp, Q - mq
    a = (p * q).sum()
    b = (p[:, 0] * q[:, 1] - p[:, 1] * q[:, 0]).sum()
    n = (p ** 2).sum()
    c, s = a / n, b / n
    Rm = np.array([[c, -s], [s, c]])
    t = mq - Rm @ mp
    return np.c_[Rm, t]


def _poly_terms(order):
    return [(i, d - i) for d in range(order + 1) for i in range(d, -1, -1)]


def fit_transform(P, Q, order=1, model='affine'):
    """Fit Q = T(P).  model 'shift' | 'similarity' | 'affine' (order 1 general) | 'poly2' | 'poly3'.  Returns coefficient dict (JSON-able)."""
    if model == 'shift':
        t = (Q - P).mean(0)
        return dict(kind='shift', M=[[1, 0, float(t[0])], [0, 1, float(t[1])]])
    if model == 'similarity':
        return dict(kind='similarity', M=_similarity(P, Q).tolist())
    deg = {'affine': 1, 'poly2': 2, 'poly3': 3}[model]
    c = P.mean(0); sc = max(P.std(0).max(), 1.0)
    u = (P - c) / sc
    A = np.array([u[:, 0] ** i * u[:, 1] ** j for i, j in _poly_terms(deg)]).T
    cx = np.linalg.lstsq(A, Q[:, 0], rcond=None)[0]
    cy = np.linalg.lstsq(A, Q[:, 1], rcond=None)[0]
    return dict(kind='poly', degree=deg, centre=c.tolist(), scale=float(sc), cx=cx.tolist(), cy=cy.tolist())


def apply_transform(T, P):
    P = np.asarray(P, float).reshape(-1, 2)
    if T['kind'] in ('shift', 'similarity'):
        M = np.array(T['M'])
        return P @ M[:, :2].T + M[:, 2]
    u = (P - np.array(T['centre'])) / T['scale']
    A = np.array([u[:, 0] ** i * u[:, 1] ** j for i, j in _poly_terms(T['degree'])]).T
    return np.c_[A @ np.array(T['cx']), A @ np.array(T['cy'])]


def match_lists(A, B, model='affine', radius=2.0, n_bright=25, guess=None, max_iter=10, clip=3.0):
    """DAOMATCH: match star lists A and B (dict/arrays x, y [, mag]; B plays the reference) and fit the transform A -> B.

    1. initial similarity from triangle matching of the n_bright brightest stars (or `guess`, a transform dict);
    2. iterative nearest-neighbour matching with growing transform order and sigma clipping.
    Returns dict(ok, transform, n_match, rms, pairs (ia, ib, sep), dmag_median, dmag_scatter)."""
    PA = np.c_[A['x'], A['y']].astype(float)
    PB = np.c_[B['x'], B['y']].astype(float)
    ma = np.asarray(A.get('mag', np.arange(len(PA))), float)
    mb = np.asarray(B.get('mag', np.arange(len(PB))), float)
    if len(PA) < 3 or len(PB) < 3:
        return dict(ok=False, reason='fewer than 3 stars', n_match=0)
    ia = np.argsort(ma)[:n_bright]
    ib = np.argsort(mb)[:n_bright]
    trB = cKDTree(PB)
    def score(T):
        q = apply_transform(T, PA)
        d, _ = trB.query(q)
        return int((d < radius).sum())
    if guess is not None:
        T0 = guess
    else:
        fa, ta = _tri_features(PA, ia[:min(len(ia), 20)])
        fb, tb = _tri_features(PB, ib[:min(len(ib), 20)])
        best, T0 = -1, None
        if len(fa) and len(fb):
            kd = cKDTree(fb)
            cand = kd.query_ball_point(fa, 0.01)
            tried = 0
            for i, lst in enumerate(cand):
                for j in lst:
                    # vertices in sorted-side order are ambiguous: use the similarity from the 3 points in both orders
                    Pa = PA[list(ta[i])]
                    Pb = PB[list(tb[j])]
                    for perm in ((0, 1, 2), (0, 2, 1), (1, 0, 2), (1, 2, 0), (2, 0, 1), (2, 1, 0)):
                        M = _similarity(Pa, Pb[list(perm)])
                        if abs(np.linalg.det(M[:, :2])) < 1e-6:
                            continue
                        T = dict(kind='similarity', M=M.tolist())
                        sc_ = score(T)
                        if sc_ > best:
                            best, T0 = sc_, T
                    tried += 1
                    if tried > 3000:
                        break
                if tried > 3000:
                    break
        if T0 is None or best < 4:
            # fall back: no rotation / scale, shift from the brightest-star offset histogram
            dxy = (PB[ib][None, :, :] - PA[ia][:, None, :]).reshape(-1, 2)
            H, xe, ye = np.histogram2d(dxy[:, 0], dxy[:, 1], bins=int(max(20, 2 * np.ptp(dxy[:, 0]) / radius)))
            k = np.unravel_index(np.argmax(H), H.shape)
            T0 = dict(kind='shift', M=[[1, 0, float(0.5 * (xe[k[0]] + xe[k[0] + 1]))], [0, 1, float(0.5 * (ye[k[1]] + ye[k[1] + 1]))]])
    T = T0
    models = ['similarity'] if model in ('shift', 'similarity') else ['similarity', 'affine'] + (['poly2'] if model in ('poly2', 'poly3') else []) + (['poly3'] if model == 'poly3' else [])
    if model == 'shift':
        models = ['shift']
    pairs = []
    for mname in models:
        for _ in range(max_iter):
            q = apply_transform(T, PA)
            d, j = trB.query(q)
            # mutual nearest neighbour
            trQ = cKDTree(q)
            _, back = trQ.query(PB)
            ok = (d < radius) & (back[j] == np.arange(len(PA)))
            if ok.sum() < 3:
                break
            ii, jj = np.where(ok)[0], j[ok]
            for _c in range(3):                  # sigma clip
                Tn = fit_transform(PA[ii], PB[jj], model=mname)
                res = np.hypot(*(apply_transform(Tn, PA[ii]) - PB[jj]).T)
                s = max(1.4826 * np.median(res), 0.05)
                keep = res < max(clip * s, 0.2)
                if keep.all() or keep.sum() < 3:
                    break
                ii, jj = ii[keep], jj[keep]
            T = fit_transform(PA[ii], PB[jj], model=mname)
            if len(pairs) and len(ii) == len(pairs) and set(ii) == set(p[0] for p in pairs):
                pairs = [(int(a), int(b)) for a, b in zip(ii, jj)]
                break
            pairs = [(int(a), int(b)) for a, b in zip(ii, jj)]
    if not pairs:
        return dict(ok=False, reason='no stable match', n_match=0)
    ia_, ib_ = np.array([p[0] for p in pairs]), np.array([p[1] for p in pairs])
    sep = np.hypot(*(apply_transform(T, PA[ia_]) - PB[ib_]).T)
    dm = (mb[ib_] - ma[ia_]) if 'mag' in A and 'mag' in B else np.array([])
    out = dict(ok=True, transform=T, n_match=int(len(pairs)), rms=float(np.sqrt(np.mean(sep ** 2))), pairs=[(int(a), int(b), float(s_)) for a, b, s_ in zip(ia_, ib_, sep)])
    if dm.size:
        med = float(np.median(dm))
        out.update(dmag_median=med, dmag_scatter=float(1.4826 * np.median(np.abs(dm - med))))
    return out


def master_list(frames, ref=0, model='affine', radius=2.0):
    """DAOMASTER: match every frame (dicts x, y, mag) to frame `ref`, fit the per-frame magnitude zero-point offset, and return the master
    list (mean position in the reference system, mean magnitude, scatter, number of frames, per-frame magnitudes)."""
    R = frames[ref]
    nref = len(R['x'])
    tab = {i: dict(x=[R['x'][i]], y=[R['y'][i]], m=[R['mag'][i]], n=1) for i in range(nref)}
    mags = [[float('nan')] * nref for _ in frames]
    for i in range(nref):
        mags[ref][i] = R['mag'][i]
    info = []
    for k, F in enumerate(frames):
        if k == ref:
            info.append(dict(frame=k, transform=None, n_match=nref, dmag=0.0))
            continue
        m = match_lists(F, R, model=model, radius=radius)
        info.append(dict(frame=k, ok=m['ok'], n_match=m.get('n_match', 0), rms=m.get('rms'), dmag=m.get('dmag_median'), transform=m.get('transform')))
        if not m['ok']:
            continue
        off = m.get('dmag_median', 0.0)
        for a, b, _s in m['pairs']:
            q = apply_transform(m['transform'], [[F['x'][a], F['y'][a]]])[0]
            mg = F['mag'][a] + off
            mags[k][b] = mg
            t = tab[b]; t['x'].append(q[0]); t['y'].append(q[1]); t['m'].append(mg); t['n'] += 1
    rows = []
    for i in range(nref):
        t = tab[i]
        rows.append(dict(id=i + 1, x=float(np.mean(t['x'])), y=float(np.mean(t['y'])), mag=float(np.mean(t['m'])),
                         sigma=float(np.std(t['m'], ddof=1)) if t['n'] > 1 else float('nan'), nframes=t['n'], mags=[mags[k][i] for k in range(len(frames))]))
    return dict(master=rows, frames=info)


# ---------------------------------------------------------------------------------------------------------- completeness hook
class DaophotDetector:
    """Detection callback for ogfkit/plugins completeness: runs find + allstar with a frozen PSF model (dict) on the image."""
    def __init__(self, model_dict, thresh=4.0, n_iter=2, zp=25.0, fwhm=None, snr_min=3.0, chi_max=3.0, aperture_mag=False):
        self.model_dict, self.thresh, self.n_iter, self.zp, self.fwhm, self.snr_min, self.chi_max = model_dict, thresh, n_iter, zp, fwhm, snr_min, chi_max
        self._model = None

    def __call__(self, img, mask=None):
        if self._model is None:
            self._model = pm.PSFModel.from_dict(self.model_dict)
        res, _ = allstar(img, self._model, mask=mask, fwhm=self.fwhm, thresh=self.thresh, n_iter=self.n_iter, zp=self.zp, snr_min=self.snr_min,
                         chi_max=self.chi_max, n_workers=1)
        g = [s for s in res['stars'] if s['good'] and np.isfinite(s['mag'])]
        return dict(x=np.array([s['x'] for s in g]), y=np.array([s['y'] for s in g]), mag=np.array([s['mag'] for s in g]))


# ------------------------------------------------------------------------------------------------------------------ pipeline
def _parse_xy_list(text):
    """'x,y;x,y' (1-based, as shown in the catalog) -> list of 0-based (x, y)."""
    out = []
    for part in str(text or '').replace('\n', ';').split(';'):
        part = part.strip()
        if part:
            a = part.replace(' ', ',').split(',')
            a = [v for v in a if v != '']
            out.append((float(a[0]) - 1.0, float(a[1]) - 1.0))
    return out


def run_pipeline(data, mask=None, zp=25.0, fwhm=3.0, find_thresh=4.0, find_sharp=(0.2, 1.0), find_round=(-1.0, 1.0), radii=(3.0, 5.0, 8.0),
                 sky_inner=12.0, sky_outer=18.0, sky_mode='median', psf_kind='empirical', psf_order=2, psf_lookup=True, psf_size=0,
                 psf_nstars=60, psf_snr=20.0, saturation=None, psf_add=None, psf_remove=None, psf_model=None, neighbour_iter=2,
                 fit_radius=0.0, fit_sky=True, gain=None, n_iter=3, chi_max=3.0, snr_min=3.0, sharp_cut=(-1.0, 1.0), round_cut=(-1.0, 1.0),
                 max_group=20, apcorr_radius=10.0, n_workers=1, log=None):
    """FIND -> PHOT -> PICKPSF -> PSF -> ALLSTAR -> aperture correction.  Returns (result dict, residual image, PSFModel).

    `psf_model` (PSFModel) skips the PSF stages.  psf_add / psf_remove: lists of 0-based (x, y)."""
    log = log or (lambda *a: None)
    data = np.asarray(data, np.float32)
    bkg, rms = pm.background(data, mask)
    sub = np.nan_to_num(data - bkg, nan=0.0)
    f0 = find(data, bkg, rms, fwhm=fwhm, thresh=find_thresh, sharp=find_sharp, round_=find_round, mask=mask)
    ix = np.clip(np.round(f0['y']).astype(int), 0, data.shape[0] - 1)
    jx = np.clip(np.round(f0['x']).astype(int), 0, data.shape[1] - 1)
    f0['snr'] = f0['peak'] / rms[ix, jx] if len(ix) else np.zeros(0)
    log('FIND: %d candidates (thresh %.1f sigma, sharp %s, round %s)' % (len(f0['x']), find_thresh, find_sharp, find_round))
    res = dict(n_find=int(len(f0['x'])))
    psf_info = {}
    pick_info = {}
    psf_xy = np.zeros((0, 2))
    model = psf_model
    if model is None:
        half = int(max(8, math.ceil(2.6 * fwhm)))
        psf_xy, pick_info = pick_psf(f0, data.shape, n_max=psf_nstars, half=half, snr_min=psf_snr, saturation=saturation, add=psf_add, remove=psf_remove)
        log('PICKPSF: %d stars selected (%s)' % (len(psf_xy), pick_info['reasons']))

        def nb_fit(mdl, _sub):
            r, _ = allstar(data, mdl, bkg=bkg, rms=rms, mask=mask, init_stars=f0, fwhm=fwhm, thresh=find_thresh, n_iter=1, refind=False,
                           fit_sky=False, zp=zp, n_workers=n_workers, snr_min=snr_min)
            s = [q for q in r['stars'] if q['flux'] > 0]
            return [q['x'] for q in s], [q['y'] for q in s], [q['flux'] for q in s]

        model, psf_info = pm.build_psf_model(data, bkg=bkg, rms=rms, mask=mask, xy=psf_xy if len(psf_xy) else None, fwhm_prior=fwhm,
                                             size=int(psf_size) or None, degree=psf_order, saturation=saturation, snr_min=psf_snr,
                                             neighbour_fit=nb_fit if (neighbour_iter > 0 and len(psf_xy) >= 4) else None, n_iter=max(1, neighbour_iter + 1),
                                             kind=psf_kind, lookup=psf_lookup)
        log('PSF: mode %s, %s stars used, order %s, FWHM(centre) %.2f px' % (psf_info.get('mode'), psf_info.get('n_used', psf_info.get('n_stars')), psf_info.get('degree'),
                                                                          model.fwhm_estimate()))
    res['psf'] = {k: v for k, v in psf_info.items() if k != 'stars'}
    res['psf_stars'] = psf_info.get('stars', [])
    res['pick'] = pick_info
    fres, resid = allstar(data, model, bkg=bkg, rms=rms, mask=mask, init_stars=f0, fwhm=fwhm, thresh=find_thresh, n_iter=n_iter, fit_radius=fit_radius or None,
                          fit_sky=fit_sky, gain=gain, zp=zp, max_group=max_group, snr_min=snr_min, chi_max=chi_max, sharp_cut=sharp_cut, round_cut=round_cut,
                          find_sharp=find_sharp, find_round=find_round, n_workers=n_workers)
    stars = fres['stars']
    log('ALLSTAR: %d stars (%d good), found per iteration %s' % (fres['n'], fres['n_good'], fres['n_found_per_iter']))
    # PHOT at the fitted positions on the raw image
    ph = phot(data, [(s['x'], s['y']) for s in stars], radii=radii, sky_inner=sky_inner, sky_outer=sky_outer, sky_mode=sky_mode, bkg=bkg, rms=rms, zp=zp, gain=gain)
    # aperture correction
    rr = sorted(set([float(r) for r in (3, 4, 5, 6, 8, 10, 12, 15)] + [float(apcorr_radius)]))
    rr = [r for r in rr if r <= model.size // 2 + 4]
    ac = aperture_correction(sub, model, stars, psf_xy if len(psf_xy) else np.array([(s['x'], s['y']) for s in stars if s['good'] and s['snr'] > 30][:40]), radii=rr, zp=zp)
    cref = float('nan')
    if ac['n']:
        k = int(np.argmin(np.abs(np.array(ac['radii']) - apcorr_radius)))
        cref = ac['corr'][k]
    for s, p in zip(stars, ph):
        s['sky_ap'], s['skysig'] = p['sky'], p['skysig']
        s['ap_mag'], s['ap_magerr'], s['ap_flux'] = p['mag'], p['magerr'], p['flux']
        s['mag_apc'] = s['mag'] + cref if np.isfinite(cref) else float('nan')
    res.update(stars=stars, n_stars=fres['n'], n_good=fres['n_good'], n_found_per_iter=fres['n_found_per_iter'], fit_radius=fres['fit_radius'], fwhm=fwhm,
               apcorr=ac, apcorr_radius=float(apcorr_radius), apcorr_value=cref, zp=zp,
               psf_diag=pm.diagnostics(model, 5), bkg_median=float(np.median(bkg)), rms_median=float(np.median(rms)))
    return res, resid, model
