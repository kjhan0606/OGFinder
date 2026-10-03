"""PSF model extensions: reduced-rank (PCA) spatial variation, bright-star wings / halo, convolution error budget, and stacked-star PSFs.

* reduce_rank(model, rank)      the polynomial coefficient images of a PSFModel (mean + K-1 polynomial variation terms) are replaced by the mean plus the `rank` leading singular
                                images of the variation: PSF(x, y) = mean + sum_j a_j(x, y) E_j with polynomial a_j (PCA in position).  Returns a new PSFModel (so psfex / daophot / multifit read it).
* measure_wings / fit_wings     azimuthally averaged profile of bright, unsaturated, isolated stars out to rmax, normalised to the flux inside r_core, and a power-law fit of the outer part.
* extended_psf(model, wings)    core from the model, wings from the fitted power law (circular), joined smoothly at r_join; unit sum including the tail beyond the stamp.
* error_budget(...)             what a wrong / truncated / mis-centred PSF does to galaxy photometry: max |dModel|/peak, flux error and the bias of (mag, R_e, n) of a Sersic fit.
* stack_star_psf(data, xy, ...)  median stack of recentred stars with azimuthally averaged wings (the pixel noise of faint-star wings would otherwise make the PSF sum negative).
"""
import math

import numpy as np

from . import models as M, multifit as MF, psfmodel as PM


# ------------------------------------------------------------------------------------------------------------------ reduced rank
def reduce_rank(model, rank):
    """New PSFModel with the coefficient images of the polynomial terms (without the constant) restricted to `rank` singular components."""
    K, nf, _ = model.coeffs.shape
    if K <= 1 or rank >= K - 1:
        return PM.PSFModel(model.size, model.oversample, model.degree, model.coeffs.copy(), model.xr, model.yr, dict(model.meta, rank=int(min(rank, K - 1))))
    C = model.coeffs[1:].reshape(K - 1, nf * nf)
    U, S, Vt = np.linalg.svd(C, full_matrices=False)
    r = int(max(rank, 0))
    Cr = (U[:, :r] * S[:r]) @ Vt[:r] if r > 0 else np.zeros_like(C)
    new = np.concatenate([model.coeffs[:1], Cr.reshape(K - 1, nf, nf)], axis=0)
    meta = dict(model.meta, rank=r, rank_variance_kept=float((S[:r] ** 2).sum() / max((S ** 2).sum(), 1e-300)))
    return PM.PSFModel(model.size, model.oversample, model.degree, new, model.xr, model.yr, meta)


def singular_values(model):
    K, nf, _ = model.coeffs.shape
    if K <= 1:
        return np.zeros(0)
    return np.linalg.svd(model.coeffs[1:].reshape(K - 1, nf * nf), compute_uv=False)


# ------------------------------------------------------------------------------------------------------------------ wings
def measure_wings(data, xy, bkg=None, rcore=6.0, rmax=40.0, mask=None, saturation=None, min_dist=None, nmax=40):
    """Median normalised radial profile of isolated bright stars: returns (r_centres, profile [per pixel, per unit flux inside rcore], n_stars).
    The sky is `bkg` if given (preferred: pass the background MAP of the image) else the median of the star's outermost annulus (0.9-1.0 rmax; this absorbs part of a halo), neighbours are rejected by a 3 sigma clip in each annulus."""
    data = np.asarray(data, float)
    if bkg is not None:                                              # a supplied sky (scalar or map) is used as is: the outer-ring sky below would absorb a halo
        data = data - bkg
    ny, nx = data.shape
    xy = np.asarray(xy, float).reshape(-1, 2)
    h = int(math.ceil(rmax)) + 2
    edges = np.concatenate([np.arange(0, 3, 0.5), np.arange(3, rmax + 1e-6, 1.0)]) if rmax > 3 else np.arange(0, rmax + 1e-6, 0.5)
    rc = 0.5 * (edges[1:] + edges[:-1])
    md = min_dist if min_dist is not None else 2.0 * rmax / 3.0
    profs = []
    for i, (x, y) in enumerate(xy):
        ix, iy = int(round(x)), int(round(y))
        if ix - h < 0 or iy - h < 0 or ix + h >= nx or iy + h >= ny:
            continue
        if md > 0 and len(xy) > 1:
            d = np.hypot(xy[:, 0] - x, xy[:, 1] - y)
            d[i] = 1e9
            if d.min() < md:
                continue
        st = data[iy - h:iy + h + 1, ix - h:ix + h + 1]
        if saturation is not None and np.nanmax(st) > saturation:
            continue
        yy, xx = np.mgrid[-h:h + 1, -h:h + 1]
        r = np.hypot(xx - (x - ix), yy - (y - iy))
        m = np.isfinite(st) if mask is None else (np.isfinite(st) & ~mask[iy - h:iy + h + 1, ix - h:ix + h + 1])
        out = (r > 0.9 * rmax) & (r < rmax + 1) & m
        if out.sum() < 20:
            continue
        sky = 0.0 if bkg is not None else np.median(st[out])
        s = st - sky
        fc = s[(r < rcore) & m].sum()
        if fc <= 0:
            continue
        pr = []
        for a, b in zip(edges[:-1], edges[1:]):
            sel = (r >= a) & (r < b) & m
            v = s[sel]
            if len(v) > 4:
                med = np.median(v)
                sg = 1.4826 * np.median(np.abs(v - med)) + 1e-30
                v = v[np.abs(v - med) < 3 * max(sg, 1e-30)] if b > rcore else v
            pr.append(v.mean() / fc if len(v) else np.nan)
        profs.append(pr)
        if len(profs) >= nmax:
            break
    if not profs:
        return rc, np.full(len(rc), np.nan), 0
    P = np.array(profs)
    return rc, np.nanmedian(P, axis=0), len(profs)


def powerlaw(r, A, gamma, r0):
    return A * (np.asarray(r, float) / r0) ** (-gamma)


def fit_wings(r, prof, rmin=8.0, rmax=None):
    """Power-law fit I(r) = A (r / r0)^-gamma (r0 = rmin) of the normalised profile for r in [rmin, rmax] (weighted log-space least squares; a Moffat tail is a power law of index 2 beta far out).
    -> dict(A, gamma, r0, r_fit, rms_dex)  or None.  A is per pixel per unit flux inside r_core.  `valid` = index > 2.1 (a flatter fit has an unbounded flux: the wing is not measurable)."""
    sel = np.isfinite(prof) & (prof > 0) & (r >= rmin) & ((r <= rmax) if rmax else True)
    if sel.sum() < 4:
        return None
    x, y = np.log(r[sel] / rmin), np.log(prof[sel])
    b, a = np.polyfit(x, y, 1)
    res = y - (a + b * x)
    return dict(A=float(math.exp(a)), gamma=float(-b), valid=bool(-b > 2.1), r0=float(rmin), rms_dex=float(np.sqrt(np.mean(res ** 2)) / math.log(10)), r_fit=[float(r[sel][0]), float(r[sel][-1])])


def wing_flux(w, r1, r2=np.inf):
    """Flux (per unit flux inside r_core) of the fitted power-law wing between r1 and r2: 2 pi A r0^2 / (gamma - 2) [(r1/r0)^(2-gamma) - (r2/r0)^(2-gamma)]."""
    g, r0 = w['gamma'], w['r0']
    f = lambda r: (r / r0) ** (2.0 - g) if np.isfinite(r) else 0.0
    return 2.0 * math.pi * w['A'] * r0 ** 2 / (g - 2.0) * (f(r1) - f(r2))


def extended_psf(model, wings, size=121, x=None, y=None, r_join=None):
    """Extended unit-sum PSF stamp (size x size, odd): the model core inside r_join, the circular power-law wing outside (scaled to the azimuthal mean of the core at r_join, blended over 2 px);
    the analytic tail beyond the stamp is included in the normalisation.  -> (stamp, info dict with tail_fraction, join_radius)"""
    x = 0.5 * sum(model.xr) if x is None else x
    y = 0.5 * sum(model.yr) if y is None else y
    core = model.stamp(x, y, 0.0, 0.0)
    n = core.shape[0]
    cy = n // 2
    rj = r_join if r_join is not None else max(0.5 * cy, 3.0)
    h = size // 2
    yy, xx = np.mgrid[-h:h + 1, -h:h + 1]
    rr = np.hypot(xx, yy)
    cyy, cxx = np.mgrid[-cy:cy + 1, -cy:cy + 1]
    cr = np.hypot(cxx, cyy)
    target = core[(cr >= rj - 1) & (cr < rj + 1)].mean()
    g = wings['gamma']
    wing = target * (np.maximum(rr, 0.5) / rj) ** (-g)
    full = np.zeros((size, size))
    sl = slice(h - cy, h + cy + 1)
    full[sl, sl] = core
    inside = np.zeros((size, size), bool)
    inside[sl, sl] = True
    wgt = np.clip((rr - rj) / 2.0, 0.0, 1.0)
    ext = np.where(inside & (rr <= cy), (1 - wgt) * full + wgt * wing, wing)
    tail = 2.0 * math.pi * target * rj ** 2 / (g - 2.0) * (h / rj) ** (2.0 - g)
    tot = ext.sum() + tail
    return ext / tot, dict(tail_fraction=float(tail / tot), join_radius=float(rj), gamma=float(g))


# ------------------------------------------------------------------------------------------------------------------ convolution error budget
def _fit_bias(psf_true, psf_used, comp, shape=(101, 101)):
    t = MF._normalise(dict(comp), 25.0)
    img = MF.render_model([t], shape, psf=psf_true)
    st = dict(comp)
    st.update(flux=comp['flux'] * 1.05, re=comp['re'] * 1.1)
    if 'n' in comp:
        st['n'] = comp['n'] * 1.1
    res = MF.fit(img, [st], psf=psf_used, rms=max(img.max() * 1e-3, 1e-6), sky='fixed', sky_value=0.0, max_nfev=150, zp=25.0)
    c = res['components'][0]
    out = dict(dmag=c['mag'] - (25.0 - 2.5 * math.log10(comp['flux'])), dre_rel=c['re'] / comp['re'] - 1)
    if 'n' in comp:
        out['dn_rel'] = c['n'] / comp['n'] - 1
    return out


def error_budget(psf_true, psf_used, profiles=((1.0, 3.0), (1.0, 8.0), (4.0, 3.0), (4.0, 8.0)), fit=True):
    """Effect of using `psf_used` where the image was formed with `psf_true` (2-D unit-sum arrays): for Sersic (n, R_e) profiles the maximum |dModel|/peak, the relative flux change of the
    convolved model inside the 101 x 101 box and (fit=True) the bias of a free fit (mag, R_e, n) of the true-PSF image with the used PSF.  -> list of dicts"""
    out = []
    for n, re in profiles:
        comp = dict(kind='sersic', x=50.0, y=50.0, flux=1e5, re=re, n=n, q=0.8, pa=30.0)
        t = MF._normalise(dict(comp), 25.0)
        a = MF.render_model([t], (101, 101), psf=psf_true)
        b = MF.render_model([t], (101, 101), psf=psf_used)
        d = dict(n=n, re=re, max_abs_over_peak=float(np.abs(a - b).max() / a.max()), flux_rel=float(b.sum() / a.sum() - 1))
        if fit:
            d.update(_fit_bias(psf_true, psf_used, comp))
        out.append(d)
    return out


def truncation_loss(wings, size, core_fraction=1.0):
    """Flux fraction (of the unit-core-flux PSF) in the fitted wing beyond a (size x size) stamp (circular approximation, r = size/2)."""
    return float(wing_flux(wings, size / 2.0) / (core_fraction + wing_flux(wings, wings['r0'])))


# ------------------------------------------------------------------------------------------------------------------ stacked-star PSF
def denoise_wings(p, core=4.0, rmax=12):
    """Keep the measured pixels inside `core` px, replace the wings by the monotonic azimuthal median profile, taper to zero at rmax, clip negatives, unit sum."""
    n = p.shape[0]
    yy, xx = np.mgrid[:n, :n] - n // 2
    rr = np.hypot(xx, yy)
    prof = np.array([np.median(p[(rr >= k - 0.5) & (rr < k + 0.5)]) if ((rr >= k - 0.5) & (rr < k + 0.5)).any() else 0.0 for k in range(n)])
    for k in range(2, len(prof)):
        prof[k] = min(prof[k], prof[k - 1])
    wing = np.interp(rr, np.arange(n), np.clip(prof, 0, None))
    w = np.clip((rr - core) / 2.0, 0.0, 1.0)
    q = (1 - w) * p + w * wing
    q = np.clip(q * np.clip((rmax - rr) / 3.0, 0.0, 1.0), 0, None)
    return q / q.sum()


def stack_star_psf(data, xy, size=31, bkg=None, core=4.0, flux_radius=6.0, nmax=60):
    """Median stack of the given stars (0-based x, y), each recentred by a windowed centroid + cubic shift, sky = median of the outer ring, normalised inside flux_radius; wings from the
    azimuthal median (see denoise_wings).  -> (psf [size x size, unit sum], n_used)"""
    import sep
    from scipy.ndimage import shift as ndshift
    d = np.ascontiguousarray(np.nan_to_num(np.asarray(data, float)))
    if bkg is not None:
        d = d - bkg
    ny, nx = d.shape
    h = size // 2 + 4
    cut = []
    for x, y in np.asarray(xy, float).reshape(-1, 2):
        ix, iy = int(round(x)), int(round(y))
        if ix - h < 1 or iy - h < 1 or ix + h >= nx - 1 or iy + h >= ny - 1:
            continue
        xw, yw, fl = sep.winpos(d, [x], [y], [1.5], subpix=5)
        x, y = float(xw[0]), float(yw[0])
        ix, iy = int(round(x)), int(round(y))
        st = d[iy - h:iy + h + 1, ix - h:ix + h + 1].copy()
        st = ndshift(st, (-(y - iy), -(x - ix)), order=3, mode='nearest')
        yy, xx = np.mgrid[:2 * h + 1, :2 * h + 1] - h
        rr = np.hypot(xx, yy)
        st = st - np.median(st[rr > h - 1])
        s = st[rr < flux_radius].sum()
        if s > 0:
            cut.append(st / s)
        if len(cut) >= nmax:
            break
    if len(cut) < 3:
        raise ValueError('too few usable stars (%d)' % len(cut))
    p = np.median(np.array(cut), axis=0)[4:-4, 4:-4]
    return denoise_wings(p, core=core, rmax=min(size // 2, 12)), len(cut)
