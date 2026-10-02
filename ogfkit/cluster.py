"""Cluster / lensing tools: red sequence from a colour-magnitude diagram, member selection, arc candidates, overdensity maps.

Pure functions on numpy arrays / plain dicts (JSON-serialisable results), no GUI or file access, so a service can call them directly.
Conventions: positions in pixels (0-based unless said otherwise), magnitudes AB, colour = blue - red band magnitude.
"""
import math

import numpy as np

try:  # scipy is a hard requirement of the plugin (declared in plugin.json)
    from scipy import ndimage as ndi
    from scipy.optimize import least_squares
except ImportError:  # pragma: no cover
    ndi = None

SQRT12 = math.sqrt(12.0)


# ------------------------------------------------------------------ red sequence

def _trunc_var_factor(k):
    """Variance of a unit Gaussian truncated at +-k (<1): corrects the scatter of a clipped sample."""
    phi = math.exp(-0.5 * k * k) / math.sqrt(2 * math.pi)
    cdf = 0.5 * (1 + math.erf(k / math.sqrt(2)))
    return 1.0 - 2.0 * k * phi / (2.0 * cdf - 1.0)


def _wls(m, c, w):
    """Weighted straight-line fit c = a + b m; returns (a, b, cov)."""
    A = np.vstack([np.ones_like(m), m]).T
    W = np.sqrt(w)[:, None]
    sol, *_ = np.linalg.lstsq(A * W, c * W[:, 0], rcond=None)
    cov = np.linalg.inv((A * w[:, None]).T @ A)
    return sol[0], sol[1], cov


def red_sequence(color, mag, color_err=None, mag_range=None, color_range=(-1.0, 3.0), slope_range=(-0.20, 0.10), fit_slope=True,
                 nsig=2.5, niter=20, m0=None, bin_width=0.04, min_n=8, window=0.6):
    """Fit the red sequence  colour = a + b (mag - m0)  of a colour-magnitude diagram.

    Start: a grid of slopes; for every slope the colours are de-trended and the highest narrow-bin peak of the histogram wins (a tilted
    sequence is narrowest at the true slope).  Then iterated clipped least squares gives the starting line and a maximum-likelihood fit of a Gaussian sequence
    (intrinsic scatter + measurement error per object) on a uniform background inside +-window of the line gives colour, slope, scatter and their errors.
    Returns a dict: ok, a (colour at m0), b (slope), m0, scatter (intrinsic), rms (total), n, a_err, b_err, nsig, note."""
    c = np.asarray(color, float)
    m = np.asarray(mag, float)
    e = np.zeros_like(c) if color_err is None else np.nan_to_num(np.asarray(color_err, float), nan=0.0)
    ok = np.isfinite(c) & np.isfinite(m) & (c >= color_range[0]) & (c <= color_range[1])
    if mag_range is not None:
        ok &= (m >= mag_range[0]) & (m <= mag_range[1])
    out = dict(ok=False, a=float('nan'), b=float('nan'), m0=float('nan'), scatter=float('nan'), rms=float('nan'), n=0, a_err=float('nan'), b_err=float('nan'), nsig=nsig, note='')
    if ok.sum() < max(min_n, 5):
        out['note'] = 'too few objects in the magnitude / colour range (%d)' % int(ok.sum())
        return out
    c, m, e = c[ok], m[ok], e[ok]
    m0 = float(np.median(m)) if m0 is None else float(m0)
    dm = m - m0
    slopes = np.linspace(slope_range[0], slope_range[1], 61) if fit_slope else np.array([0.0])
    edges = np.arange(color_range[0], color_range[1] + bin_width, bin_width)
    best = (-1, 0.0, 0.0)
    for s in slopes:
        h, _ = np.histogram(c - s * dm, bins=edges)
        hs = np.convolve(h, np.ones(3), mode='same')
        i = int(np.argmax(hs))
        if hs[i] > best[0]:
            best = (hs[i], s, 0.5 * (edges[i] + edges[i + 1]))
    s0, a0 = best[1], best[2]
    near = np.abs(c - s0 * dm - a0) < 0.1
    a0 = float(np.median((c - s0 * dm)[near])) if near.any() else a0
    a, b, sig = a0, (s0 if fit_slope else 0.0), 0.08
    sel_prev = None
    tf = _trunc_var_factor(nsig)
    cov = np.eye(2)
    for _ in range(niter):
        res = c - (a + b * dm)
        tot = np.sqrt(sig ** 2 + e ** 2)
        sel = np.abs(res) < nsig * tot
        if sel.sum() < max(min_n, 4):
            out['note'] = 'red sequence lost during clipping'
            return out
        w = 1.0 / tot[sel] ** 2
        if fit_slope:
            a_, b_, cov = _wls(dm[sel], c[sel], w)
        else:
            a_ = float(np.sum(w * c[sel]) / np.sum(w)); b_ = 0.0
            cov = np.array([[1.0 / np.sum(w), 0], [0, 0]])
        rms = math.sqrt(float(np.sum(w * (c[sel] - (a_ + b_ * dm[sel])) ** 2) / np.sum(w)) / tf)
        sig_new = math.sqrt(max(rms ** 2 - float(np.mean(e[sel] ** 2)), 1e-4))
        conv = sel_prev is not None and np.array_equal(sel, sel_prev) and abs(sig_new - sig) < 1e-4
        a, b, sig, sel_prev = a_, b_, sig_new, sel
        if conv:
            break
    # final stage: maximum likelihood of  f N(res; 0, sigma^2 + err^2) + (1 - f) / (2 W)  inside the window |res| < W (the clipped fit above is only the
    # starting point; the mixture does not let the scatter collapse when the window shrinks)
    W = float(window)
    inwin = np.abs(c - (a + b * dm)) < W
    ci, di, ei = c[inwin], dm[inwin], e[inwin]
    f0 = min(max(sel.sum() / max(inwin.sum(), 1), 0.05), 0.95)
    def nll(p):
        a_, b_, ls_, lf_ = p
        if not fit_slope:
            b_ = 0.0
        s2 = (0.005 + math.exp(ls_)) ** 2 + ei ** 2
        f = 1.0 / (1.0 + math.exp(-lf_))
        res = ci - (a_ + b_ * di)
        # the sample is fixed (objects inside the window of the starting line), so the likelihood is smooth in the parameters
        L = f * np.exp(-0.5 * res ** 2 / s2) / np.sqrt(2 * np.pi * s2) + (1 - f) / (2 * W)
        return -float(np.sum(np.log(np.maximum(L, 1e-300))))
    from scipy.optimize import minimize
    r = None
    for sg0 in (0.03, 0.06, 0.12):                    # several starts: the likelihood has a narrow local mode on a few tight objects
        p0 = np.array([a, b, math.log(sg0), math.log(f0 / (1 - f0))])
        rr = minimize(nll, p0, method='Nelder-Mead', options=dict(xatol=1e-5, fatol=1e-7, maxiter=4000))
        if r is None or rr.fun < r.fun:
            r = rr
    pbest = r.x
    H = np.zeros((4, 4)); h = np.array([1e-3, 1e-3, 1e-2, 1e-2])
    f00 = nll(pbest)
    for i_ in range(4):
        for j_ in range(i_, 4):
            ei_, ej_ = np.eye(4)[i_] * h[i_], np.eye(4)[j_] * h[j_]
            H[i_, j_] = H[j_, i_] = (nll(pbest + ei_ + ej_) - nll(pbest + ei_ - ej_) - nll(pbest - ei_ + ej_) + nll(pbest - ei_ - ej_)) / (4 * h[i_] * h[j_])
    idx = [0, 1, 2, 3] if fit_slope else [0, 2, 3]
    try:
        cv = np.linalg.inv(H[np.ix_(idx, idx)])
        a_err = math.sqrt(max(cv[0, 0], 0)); b_err = math.sqrt(max(cv[1, 1], 0)) if fit_slope else 0.0
    except np.linalg.LinAlgError:
        a_err = b_err = float('nan')
    a, b = float(pbest[0]), float(pbest[1]) if fit_slope else 0.0
    sig = 0.005 + math.exp(pbest[2]); frac = 1.0 / (1.0 + math.exp(-pbest[3]))
    res_all = c - (a + b * dm)
    inside = np.abs(res_all) < nsig * np.sqrt(sig ** 2 + e ** 2)
    out.update(ok=True, a=float(a), b=float(b), m0=m0, scatter=float(sig), rms=float(np.sqrt(np.mean(res_all[inside] ** 2))), n=int(inside.sum()), n_model=float(frac * inwin.sum()),
               a_err=a_err, b_err=b_err, converged=bool(r.success))
    return out


def rs_pull(color, mag, rs, color_err=None):
    """(colour - red sequence) and its pull in units of sqrt(scatter^2 + err^2).  NaN where colour/mag are missing."""
    c = np.asarray(color, float); m = np.asarray(mag, float)
    e = np.zeros_like(c) if color_err is None else np.nan_to_num(np.asarray(color_err, float), nan=0.0)
    res = c - (rs['a'] + rs['b'] * (m - rs['m0']))
    return res, res / np.sqrt(rs['scatter'] ** 2 + e ** 2)


def radial_pmem(r, flag, pix_dist, pix_weight=1.0, min_n=25):
    """Statistical membership probability of colour-selected objects from their radial distribution around the cluster centre.

    The surface density of the selected objects is fitted by a cluster + uniform field model  S(r) = S_bg + S_0 (1 + (r/a)^2)^-beta  (unbinned
    Poisson maximum likelihood, the integral over the area that really lies inside the image: pix_dist = sorted distances of a subsample of
    valid pixels, each worth pix_weight px^2).  p(r) = 1 - S_bg / S(r) is unbiased for the expected member count (no clipping of negative
    bin excesses).  Objects that are not selected get 0.  Returns (p, info dict with S_bg per px^2, S_0, a, beta, n_fit, converged)."""
    from scipy.optimize import minimize
    r = np.asarray(r, float)
    flag = np.asarray(flag, bool) & np.isfinite(r)
    p = np.zeros(len(r))
    rr = r[flag]
    info = dict(sigma_bg=float('nan'), s0=float('nan'), a=float('nan'), beta=float('nan'), n_fit=int(flag.sum()), converged=False)
    if len(rr) < min_n:
        info['note'] = 'too few selected objects for a radial fit'
        return p, info
    A_tot = pix_weight * len(pix_dist)
    pd = pix_dist[::max(1, len(pix_dist) // 20000)]
    pdw = pix_weight * len(pix_dist) / len(pd)

    def model(rv, th):
        sbg, s0, a, beta = np.exp(th[0]), np.exp(th[1]), np.exp(th[2]), th[3]
        return sbg + s0 * (1.0 + (rv / a) ** 2) ** (-beta)

    def nll(th):
        if not (0.6 < th[3] < 3.0):
            return 1e30
        return -float(np.sum(np.log(model(rr, th) + 1e-300))) + float(pdw * np.sum(model(pd, th)))

    sbg0 = max(float(np.sum(rr > np.quantile(pd, 0.6))) / max(pix_weight * np.sum(pix_dist > np.quantile(pd, 0.6)), 1.0), 1e-12)
    a0 = max(float(np.median(rr)) * 0.5, 10.0)
    s00 = max(len(rr) / (math.pi * a0 ** 2), 1e-9)
    best = None
    for beta0 in (1.0, 2.0):
        th0 = np.array([math.log(sbg0), math.log(s00), math.log(a0), beta0])
        res = minimize(nll, th0, method='Nelder-Mead', options=dict(xatol=1e-4, fatol=1e-6, maxiter=3000))
        if best is None or res.fun < best.fun:
            best = res
    th = best.x
    sb = float(np.exp(th[0]))
    pm = 1.0 - sb / np.maximum(model(r, th), 1e-300)
    p = np.where(flag, np.clip(pm, 0.0, 1.0), 0.0)
    info.update(sigma_bg=sb, s0=float(np.exp(th[1])), a=float(np.exp(th[2])), beta=float(th[3]), converged=bool(best.success), n_members_model=float(np.sum(p)))
    return p, info


def pixel_distances(shape, cx, cy, valid=None, step=4):
    """Sorted distances from (cx, cy) of every step-th valid pixel (rows and columns) -> (array, area per sample in px^2)."""
    yy, xx = np.mgrid[0:shape[0]:step, 0:shape[1]:step]
    d = np.hypot(xx - cx, yy - cy)
    if valid is not None:
        d = d[valid[0:shape[0]:step, 0:shape[1]:step]]
    return np.sort(d.ravel()), float(step * step)


def area_beyond(pix_dist, pix_weight, r_bg):
    return float(pix_weight * (len(pix_dist) - np.searchsorted(pix_dist, r_bg)))


def distance_quantile(pix_dist, q):
    return float(pix_dist[min(int(q * len(pix_dist)), len(pix_dist) - 1)])


def redshift_window(z, zc, dz, zerr=None):
    """True where |z - zc| < dz (1 + zc), also accepting objects whose own 1 sigma interval overlaps; NaN z -> None (undecided)."""
    z = np.asarray(z, float)
    half = dz * (1 + zc)
    pad = 0.0 if zerr is None else np.nan_to_num(np.asarray(zerr, float), nan=0.0)
    ok = np.abs(z - zc) < half + pad
    return np.where(np.isfinite(z), ok, False), np.isfinite(z)


# ------------------------------------------------------------------ overdensity

def _gauss1d(sigma):
    h = int(math.ceil(4 * sigma))
    x = np.arange(-h, h + 1, dtype=float)
    k = np.exp(-0.5 * (x / sigma) ** 2)
    return k / k.sum()


def _sepfilter(a, k):
    return ndi.correlate1d(ndi.correlate1d(a, k, axis=0, mode='constant'), k, axis=1, mode='constant')


def density_map(x, y, shape, sigma_px, binsize=4, valid=None, weights=None, clip_sig=3.0, niter=6):
    """Poisson significance map of a point sample on a regular grid.

    Counts per cell (cell = binsize px) are smoothed with a Gaussian of sigma_px.  The background rate lambda per cell is the clipped mean of
    the cells outside excesses > clip_sig; the variance of the smoothed value under a uniform Poisson process is lambda * sum(k^2) (exact for
    the discrete kernel, with the validity mask included), so significance = (smoothed - lambda) / sqrt(that).  `weights` (e.g. 0/1) are not
    used for the variance.  Returns dict(sig, dens (smoothed counts per cell), lam, binsize, ny, nx, kernel_sigma_cells)."""
    ny, nx = int(math.ceil(shape[0] / binsize)), int(math.ceil(shape[1] / binsize))
    x = np.asarray(x, float); y = np.asarray(y, float)
    ok = np.isfinite(x) & np.isfinite(y)
    w = np.ones(ok.sum()) if weights is None else np.asarray(weights, float)[ok]
    H, _, _ = np.histogram2d(y[ok], x[ok], bins=[ny, nx], range=[[0, ny * binsize], [0, nx * binsize]], weights=w)
    if valid is None:
        V = np.ones((ny, nx))
    else:
        b = int(round(binsize))
        pv = np.zeros((ny * b, nx * b)); pv[:valid.shape[0], :valid.shape[1]] = valid
        V = (pv.reshape(ny, b, nx, b).mean(axis=(1, 3)) > 0.5).astype(float)
    ks = max(sigma_px / binsize, 0.5)
    k = _gauss1d(ks)
    norm = _sepfilter(V, k)
    s = _sepfilter(H * V, k) / np.maximum(norm, 1e-9)
    k2 = k ** 2
    sumk2 = _sepfilter(V, k2) / np.maximum(norm, 1e-9) ** 2
    good = (norm > 0.5) & (V > 0)
    use = good.copy()
    lam = float(np.mean(H[use])) if use.any() else 0.0
    for _ in range(niter):
        sig = (s - lam) / np.sqrt(np.maximum(lam * sumk2, 1e-12))
        use = good & (sig < clip_sig)
        if not use.any():
            break
        lam_new = float(np.sum((H * V)[use]) / np.sum(V[use]))
        if abs(lam_new - lam) < 1e-6 * max(lam, 1e-6):
            lam = lam_new
            break
        lam = lam_new
    sig = np.where(good, (s - lam) / np.sqrt(np.maximum(lam * sumk2, 1e-12)), 0.0)
    return dict(sig=sig, dens=s, lam=lam, binsize=float(binsize), ny=ny, nx=nx, kernel_sigma_cells=float(ks), good=good)


def find_peaks(dm, threshold=3.0, min_sep_px=None, max_peaks=20):
    """Local maxima of the significance map above threshold; centroid refined over the 3x3 cells (significance-weighted).  Returns a list of dicts
    (x, y in pixels, 0-based pixel centres), sorted by significance."""
    sig = dm['sig']; b = dm['binsize']
    sep = int(max(2, round((min_sep_px or 4 * dm['kernel_sigma_cells'] * b) / b)))
    mx = ndi.maximum_filter(sig, size=2 * sep + 1, mode='constant', cval=-1e9)
    ii = np.argwhere((sig == mx) & (sig >= threshold))
    pk = []
    for iy, ix in ii:
        y0, y1, x0, x1 = max(iy - 1, 0), min(iy + 2, sig.shape[0]), max(ix - 1, 0), min(ix + 2, sig.shape[1])
        wgt = np.clip(sig[y0:y1, x0:x1], 0, None)
        yy, xx = np.mgrid[y0:y1, x0:x1]
        cx = float(np.sum(wgt * xx) / np.sum(wgt)); cy = float(np.sum(wgt * yy) / np.sum(wgt))
        pk.append(dict(x=(cx + 0.5) * b - 0.5, y=(cy + 0.5) * b - 0.5, sig=float(sig[iy, ix]), excess_cells=float(dm['dens'][iy, ix] - dm['lam'])))
    pk.sort(key=lambda d: -d['sig'])
    return pk[:max_peaks]


def sample_map(dm, x, y, key='sig'):
    """Map value at catalog positions (nearest cell); NaN outside."""
    m = dm[key]; b = dm['binsize']
    ix = np.floor((np.asarray(x, float) + 0.5) / b).astype(int); iy = np.floor((np.asarray(y, float) + 0.5) / b).astype(int)
    ok = (ix >= 0) & (ix < m.shape[1]) & (iy >= 0) & (iy < m.shape[0])
    out = np.full(len(ix), np.nan)
    out[ok] = m[iy[ok], ix[ok]]
    return out


# ------------------------------------------------------------------ arcs

def _fit_circle(x, y, w):
    """Weighted algebraic (Kasa) circle fit followed by a few geometric Gauss-Newton steps -> (cx, cy, R)."""
    xm, ym = np.average(x, weights=w), np.average(y, weights=w)
    u, v = x - xm, y - ym
    A = np.vstack([u, v, np.ones_like(u)]).T
    sw = np.sqrt(w)
    sol, *_ = np.linalg.lstsq(A * sw[:, None], -(u * u + v * v) * sw, rcond=None)
    cx, cy = -sol[0] / 2, -sol[1] / 2
    R2 = cx * cx + cy * cy - sol[2]
    R = math.sqrt(R2) if R2 > 0 else float('inf')
    if np.isfinite(R) and R < 1e4:
        f = lambda p: sw * (np.hypot(u - p[0], v - p[1]) - p[2])
        try:
            r = least_squares(f, [cx, cy, R], method='lm', max_nfev=60)
            cx, cy, R = r.x
        except Exception:
            pass
    return cx + xm, cy + ym, abs(R)


def arc_geometry(cut, sigma_sm, x0, y0, smooth=1.0, nsig=2.0, minpix=12, centre=None, curved_gain=0.85, max_curved_span=300.0, search=4.0):
    """Shape of the object at (x0, y0) (0-based, in cutout pixels) of a background-subtracted cutout.

    The object is the connected region of the Gaussian-smoothed cutout (sigma `smooth`) above nsig * sigma_sm that contains / is nearest to the
    position.  Measured: flux-weighted principal axes (a, b, theta) and elongation a/b; a circle fit to the flux-weighted pixels
    (centre, radius R, radial rms sigma_r).  The object is called *curved* when the circle shrinks the width: sigma_r < curved_gain * b.
    Lengths: L = sqrt(12) * rms of the along-ridge coordinate (arc length for curved objects) = the length of a uniform ridge; width
    W = 2.355 * (radial rms for curved, minor-axis rms for straight) = FWHM including the PSF; LW = L / W.
    angle = position angle of the major axis (deg, from +x towards +y); align = angle between the major axis and the tangential direction at
    `centre` (deg, 0 = tangential) if centre is given.  `search` = how far (px) the object may be from the given position.  Returns dict (ok False + reason when nothing is found)."""
    out = dict(ok=False, reason='')
    sm = ndi.gaussian_filter(cut, smooth)
    lab, n = ndi.label(sm > nsig * sigma_sm)
    if n:
        sizes = np.bincount(lab.ravel())
        lab = np.where(np.isin(lab, np.nonzero(sizes >= minpix)[0]) & (lab > 0), lab, 0)
        n = int(lab.max()) if (lab > 0).any() else 0
    if n == 0:
        out['reason'] = 'nothing above threshold'
        return out
    iy, ix = int(round(y0)), int(round(x0))
    l = 0
    if 0 <= iy < lab.shape[0] and 0 <= ix < lab.shape[1]:
        l = lab[iy, ix]
    if l == 0:      # the centroid of a curved arc lies off the arc: take the nearest region within `search` px
        yy, xx = np.nonzero(lab)
        d = np.hypot(xx - x0, yy - y0)
        j = int(np.argmin(d))
        if d[j] > search:
            out['reason'] = 'no region at the position'
            return out
        l = lab[yy[j], xx[j]]
    reg = lab == l
    npix = int(reg.sum())
    if npix < minpix:
        out['reason'] = 'region too small (%d px)' % npix
        return out
    ys, xs = np.nonzero(reg)
    w = np.clip(sm[reg], 1e-9, None)
    xc, yc = np.average(xs, weights=w), np.average(ys, weights=w)
    cxx = np.average((xs - xc) ** 2, weights=w); cyy = np.average((ys - yc) ** 2, weights=w); cxy = np.average((xs - xc) * (ys - yc), weights=w)
    tr, det = cxx + cyy, cxx * cyy - cxy ** 2
    disc = math.sqrt(max(tr * tr / 4 - det, 0.0))
    la, lb = tr / 2 + disc, max(tr / 2 - disc, 1e-6)
    theta = 0.5 * math.atan2(2 * cxy, cxx - cyy)                      # major axis direction (rad), +x towards +y
    ux, uy = math.cos(theta), math.sin(theta)
    s_ax = (xs - xc) * ux + (ys - yc) * uy
    sa = math.sqrt(la); sb = math.sqrt(lb)
    flux = float(np.sum(cut[reg]))
    out.update(ok=True, npix=npix, x=float(xc), y=float(yc), a=sa, b=sb, elong=sa / sb, theta_deg=math.degrees(theta) % 180.0, flux=flux)
    cx, cy, R = _fit_circle(xs.astype(float), ys.astype(float), w)
    dist = np.hypot(xs - cx, ys - cy)
    mu = np.average(dist, weights=w)
    sr = math.sqrt(np.average((dist - mu) ** 2, weights=w))
    ang = np.arctan2(ys - cy, xs - cx)
    am = math.atan2(np.average(np.sin(ang), weights=w), np.average(np.cos(ang), weights=w))
    dphi = np.angle(np.exp(1j * (ang - am)))
    span = float(np.degrees((np.percentile(dphi, 98) - np.percentile(dphi, 2)) / 0.96))
    curved = bool(np.isfinite(R) and sr < curved_gain * sb and R < 1e3 and span < max_curved_span)     # a closed blob is not an arc
    if curved:
        s = mu * dphi
        ls = SQRT12 * math.sqrt(np.average((s - np.average(s, weights=w)) ** 2, weights=w))
        out.update(curved=True, rcurv=float(mu), cx=float(cx), cy=float(cy), span_deg=span, width=2.3548 * sr, length=ls)
    else:
        ls = SQRT12 * math.sqrt(np.average((s_ax - np.average(s_ax, weights=w)) ** 2, weights=w))
        out.update(curved=False, rcurv=float('nan'), cx=float('nan'), cy=float('nan'), span_deg=0.0, width=2.3548 * sb, length=ls)
    out['lw'] = out['length'] / max(out['width'], 1e-6)
    out['snr'] = flux / (sigma_sm * 2 * math.sqrt(math.pi) * smooth * math.sqrt(npix)) if sigma_sm > 0 else float('nan')
    if centre is not None:
        rx, ry = xc - centre[0], yc - centre[1]
        rad = math.atan2(ry, rx)
        d = (theta - (rad + math.pi / 2)) % math.pi                    # major axis vs tangential direction
        out['align_deg'] = math.degrees(min(d, math.pi - d))
        out['centre_dist'] = float(math.hypot(rx, ry))
        out['cc_ratio'] = float(math.hypot(out['cx'] - centre[0], out['cy'] - centre[1]) / out['rcurv']) if curved else float('nan')
    return out


def arc_score(g, lw_min=4.0, align_max=30.0, cc_max=0.5, snr_min=5.0, require_curved=True, max_span=270.0, have_centre=True):
    """Candidate flag and a ranking score in [0, 1] (heuristic, not a probability) from an arc_geometry dict."""
    if not g.get('ok'):
        return 0, 0.0
    lw = g['lw']
    s = min(lw / 10.0, 1.0) * min(max(g['snr'], 0.0) / 20.0, 1.0)
    s *= 1.0 if g['curved'] else 0.5
    flag = lw >= lw_min and g['snr'] >= snr_min
    if require_curved:
        flag = flag and g['curved'] and g['span_deg'] < max_span
    if have_centre and 'align_deg' in g:
        s *= max(0.0, 1.0 - g['align_deg'] / 90.0)
        flag = flag and g['align_deg'] <= align_max
        if g['curved'] and np.isfinite(g.get('cc_ratio', float('nan'))):
            flag = flag and g['cc_ratio'] <= cc_max
    return int(bool(flag)), float(s)


def measure_arcs(img, x, y, sigma_sm, centre=None, half=30, **kw):
    """arc_geometry for catalog positions (x, y 0-based) on a full image: local median background is removed in each cutout.
    half may be an array (per object).  Returns a list of dicts in cutout-independent (image) coordinates."""
    ny, nx = img.shape
    halfs = np.broadcast_to(np.asarray(half, float), (len(x),))
    searches = np.broadcast_to(np.asarray(kw.pop('search', 4.0), float), (len(x),))
    res = []
    for xi, yi, h, sr_ in zip(x, y, halfs, searches):
        if not (np.isfinite(xi) and np.isfinite(yi)):
            res.append(dict(ok=False, reason='no position')); continue
        h = int(max(h, 8))
        ix, iy = int(round(xi)), int(round(yi))
        x0, x1, y0, y1 = max(ix - h, 0), min(ix + h + 1, nx), max(iy - h, 0), min(iy + h + 1, ny)
        cut = np.array(img[y0:y1, x0:x1], float)
        cut = np.nan_to_num(cut, nan=0.0)
        bkg = float(np.median(cut))
        cut -= bkg
        c = None if centre is None else (centre[0] - x0, centre[1] - y0)
        g = arc_geometry(cut, sigma_sm, xi - x0, yi - y0, centre=c, search=float(sr_), **kw)
        if g.get('ok'):
            g['x'] += x0; g['y'] += y0
            if g.get('curved'):
                g['cx'] += x0; g['cy'] += y0
        res.append(g)
    return res
