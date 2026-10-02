"""Emission-line kinematics from 2D long-slit spectra and IFU cubes (extension of ogfkit/spectra.py, round 2 item 8).

* ``fit_slit_line``       per-row (adaptively binned) Gaussian fits of one emission line along the spatial axis of a 2D spectrum -> v(y), sigma(y), flux(y)
* ``fit_rotation_curve``  arctan rotation curve v(y) = v_sys + (2/pi) v_c arctan((y - y0)/r_t) with errors, optional inclination / slit-offset deprojection
* ``fit_cube_line``       per-spaxel (accretion-binned) Gaussian fits of a line in an IFU cube -> flux, velocity, dispersion maps with errors and S/N
* ``fit_velocity_field``  2D disc model (arctan rotation curve, position angle, inclination, centre, v_sys) fitted to a velocity map

Conventions: velocities in km/s relative to the systemic wavelength lam0 (1+z) (v = c (lam/lam_sys - 1)); positions in pixels (0-based); the position angle PA
is the direction (degrees counter-clockwise from +x of the array) towards the *receding* (redshifted) side of the major axis; inclination 0 = face-on, 90 = edge-on;
v_c is the deprojected circular-velocity amplitude (the observed asymptote is v_c sin i).  Dispersions have the instrumental width removed in quadrature when it is given.
"""
import math

import numpy as np

from .spectra import fit_line, noise_level

C_KMS = 299792.458
SQRT2PI = math.sqrt(2 * math.pi)


# ------------------------------------------------------------------------------------------------------------- line fitting
def _side_continuum(wave, spec, lam0, half, gap, width):
    """mean of the two side bands (|dl| from half+gap to half+gap+width) -> per-spectrum continuum value(s); works on the last axis"""
    d = np.abs(wave - lam0)
    side = (d > half + gap) & (d <= half + gap + width)
    if side.sum() < 4:
        return np.zeros(spec.shape[:-1]), np.full(spec.shape[:-1], np.nan), side
    cont = np.nanmedian(spec[..., side], axis=-1)
    sd = 1.4826 * np.nanmedian(np.abs(spec[..., side] - cont[..., None]), axis=-1)
    return cont, sd, side


def _fit_gauss(x, y, e, mu0, s0):
    """Gaussian + constant on a continuum-subtracted window.  -> dict(ok, flux, flux_err, mu, mu_err, sigma, sigma_err, snr, chi2r)"""
    from scipy.optimize import least_squares
    out = dict(ok=False)
    m = np.isfinite(y) & np.isfinite(e) & (e > 0)
    if m.sum() < 8:
        return out
    x, y, e = x[m], y[m], e[m]
    dl = float(np.median(np.diff(x)))
    A0 = float(np.max(y) - np.median(y))
    if A0 <= 0:
        return out
    p0 = [0.0, A0, x[np.argmax(y)] if abs(x[np.argmax(y)] - mu0) < 4 * s0 else mu0, s0]
    lo = [-np.inf, 0.0, mu0 - 6 * s0 - 2 * dl, 0.35 * dl]
    hi = [np.inf, np.inf, mu0 + 6 * s0 + 2 * dl, max(10 * s0, 12 * dl)]
    p0[2] = float(np.clip(p0[2], lo[2] + 1e-9, hi[2] - 1e-9))
    try:
        r = least_squares(lambda p: (p[0] + p[1] * np.exp(-0.5 * ((x - p[2]) / p[3]) ** 2) - y) / e, p0, bounds=(lo, hi), x_scale=[max(abs(A0), 1e-12) * 0.1, max(abs(A0), 1e-12), s0, s0])
    except Exception:
        return out
    p = r.x
    dof = max(len(x) - 4, 1)
    chi2r = float(np.sum(r.fun ** 2) / dof)
    try:
        cov = np.linalg.inv(r.jac.T @ r.jac) * max(chi2r, 1.0)
    except np.linalg.LinAlgError:
        return out
    A, mu, s = p[1], p[2], p[3]
    gA, gS = s * SQRT2PI, A * SQRT2PI
    var_f = gA ** 2 * cov[1, 1] + gS ** 2 * cov[3, 3] + 2 * gA * gS * cov[1, 3]
    flux = A * s * SQRT2PI
    ferr = math.sqrt(max(var_f, 0.0))
    out.update(ok=True, flux=float(flux), flux_err=ferr, mu=float(mu), mu_err=float(math.sqrt(max(cov[2, 2], 0))), sigma=float(s), sigma_err=float(math.sqrt(max(cov[3, 3], 0))),
               snr=float(flux / ferr) if ferr > 0 else float('nan'), chi2r=chi2r, amp=float(A))
    return out


def _to_kms(fit, lam_sys, inst_sigma_A):
    """wavelength-space Gaussian -> (v, v_err, sigma_int, sigma_err) in km/s; the instrumental width (Angstrom sigma, at lam_sys) is removed in quadrature"""
    v = (fit['mu'] / lam_sys - 1.0) * C_KMS
    ve = fit['mu_err'] / lam_sys * C_KMS
    s2 = fit['sigma'] ** 2 - inst_sigma_A ** 2
    s_int = math.sqrt(s2) if s2 > 0 else 0.0
    # d(s_int) = (s / s_int) d(s)
    se = fit['sigma_err'] * fit['sigma'] / max(s_int, 0.3 * fit['sigma']) / lam_sys * C_KMS
    return v, ve, s_int / lam_sys * C_KMS, se


def inst_sigma_from_fwhm(fwhm_A):
    return float(fwhm_A) / 2.3548200450309493


# ------------------------------------------------------------------------------------------------------------- 2D slit
def adaptive_bins_1d(signal, noise, target):
    """greedy contiguous binning along one axis: rows are accumulated (starting at the brightest unbinned row, growing to the brighter neighbour) until
    sum(signal)/sqrt(sum(noise^2)) >= target.  -> list of index arrays; rows that never reach the target are left out"""
    n = len(signal)
    used = np.zeros(n, bool)
    bins = []
    sig = np.where(np.isfinite(signal), signal, 0.0)
    nz = np.where(np.isfinite(noise) & (noise > 0), noise, np.inf)
    order = np.argsort(-sig / nz)
    for i0 in order:
        if used[i0] or sig[i0] <= 0:
            continue
        lo = hi = i0
        s, v = sig[i0], nz[i0] ** 2
        while s / math.sqrt(v) < target:
            cl = lo - 1 if lo - 1 >= 0 and not used[lo - 1] else None
            ch = hi + 1 if hi + 1 < n and not used[hi + 1] else None
            if cl is None and ch is None:
                break
            pick = ch if cl is None else (cl if ch is None else (cl if sig[cl] / nz[cl] >= sig[ch] / nz[ch] else ch))
            s += sig[pick]
            v += nz[pick] ** 2
            lo, hi = min(lo, pick), max(hi, pick)
        if s / math.sqrt(v) >= target:
            idx = np.arange(lo, hi + 1)
            used[idx] = True
            bins.append(idx)
    bins.sort(key=lambda b: b[0])
    return bins


def fit_slit_line(wave, data, err, lam_sys, window_kms=900.0, bin_snr=5.0, snr_row=3.0, inst_fwhm_A=0.0, sigma_guess_kms=100.0, gap_kms=1200.0, side_kms=2500.0):
    """Gaussian fits of the emission line at lam_sys (observed wavelength) along the slit.

    data[ny, nw] (dispersion along axis 1); err same shape or None (then the noise is the MAD in the side bands, per row).  Continuum per row = median of the side bands.
    Rows are binned along the slit until the line S/N of the bin reaches ``bin_snr`` (0 = fit every row separately; rows below ``snr_row`` are then skipped).
    -> dict(y, v, v_err, sigma, sigma_err, flux, flux_err, snr, nrows, rows) arrays over the bins that gave a fit"""
    w = np.asarray(wave, float)
    D = np.asarray(data, float)
    ny = D.shape[0]
    half = window_kms / C_KMS * lam_sys
    gap, side = gap_kms / C_KMS * lam_sys, side_kms / C_KMS * lam_sys
    cont, sd, sidemask = _side_continuum(w, D, lam_sys, half, gap, side)
    R = D - cont[:, None]
    if err is not None:
        E = np.asarray(err, float)
    else:
        E = np.broadcast_to(sd[:, None], D.shape)
    inl = np.abs(w - lam_sys) <= half
    dl = np.gradient(w)
    s0 = max(sigma_guess_kms / C_KMS * lam_sys, 1.2 * float(np.median(np.diff(w))))
    inst = inst_sigma_from_fwhm(inst_fwhm_A)
    # row-wise line flux and its noise (boxcar over the window, weights = the Gaussian guess at the brightest wavelength of the collapsed profile)
    coll = np.nansum(R[:, inl], axis=0)
    ic = np.argmax(coll)
    lc = w[inl][ic] if coll[ic] > 0 else lam_sys
    wt = np.exp(-0.5 * ((w[inl] - lc) / max(s0, 1e-9)) ** 2)
    wt /= wt.sum()
    sig_row = np.nansum(R[:, inl] * wt[None], axis=1) * (1.0 / np.sum(wt ** 2))      # amplitude estimate x ... (matched filter)
    nz_row = np.sqrt(np.nansum((E[:, inl] ** 2) * (wt[None] ** 2), axis=1)) * (1.0 / np.sum(wt ** 2))
    snr_rows = sig_row / np.where(nz_row > 0, nz_row, np.inf)
    sel = np.where(np.isfinite(snr_rows) & (snr_rows >= snr_row))[0] if bin_snr <= 0 else np.arange(ny)
    if bin_snr > 0:
        bins = adaptive_bins_1d(np.where(snr_rows > 0.5, sig_row, 0.0), nz_row, bin_snr)
    else:
        bins = [np.array([i]) for i in sel]
    res = dict(y=[], v=[], v_err=[], sigma=[], sigma_err=[], flux=[], flux_err=[], snr=[], nrows=[], rows=[])
    xw = w[inl]
    for b in bins:
        spec = np.nansum(R[b][:, inl], axis=0)
        ee = np.sqrt(np.nansum(E[b][:, inl] ** 2, axis=0))
        # centre the window on the previous guess: the line may be shifted by the rotation; use the bin's own peak
        f = _fit_gauss(xw, spec, ee, lc, s0)
        if not f['ok'] or not np.isfinite(f['snr']) or f['snr'] < max(2.0, 0.5 * (bin_snr if bin_snr > 0 else snr_row)):
            continue
        v, ve, s, se = _to_kms(f, lam_sys, inst)
        wts = np.clip(np.nan_to_num(sig_row[b]), 1e-12, None)
        res['y'].append(float(np.sum(b * wts) / np.sum(wts)))
        res['v'].append(v); res['v_err'].append(ve); res['sigma'].append(s); res['sigma_err'].append(se)
        res['flux'].append(f['flux']); res['flux_err'].append(f['flux_err']); res['snr'].append(f['snr']); res['nrows'].append(int(len(b)))
        res['rows'].append([int(b[0]), int(b[-1])])
    for k in ('y', 'v', 'v_err', 'sigma', 'sigma_err', 'flux', 'flux_err', 'snr', 'nrows'):
        res[k] = np.array(res[k], float)
    res['lam_sys'] = float(lam_sys)
    res['inst_sigma_A'] = inst
    return res


# ------------------------------------------------------------------------------------------------------------- rotation curve
def arctan_rc(r, vc, rt):
    return (2.0 / math.pi) * vc * np.arctan(np.asarray(r, float) / max(rt, 1e-9))


def fit_rotation_curve(y, v, ev, inc_deg=None, slit_offset_deg=0.0, y0=None, free_centre=True, vsys_free=True, n_starts=12, seed=0):
    """v(y) = v_sys + (2/pi) v_c arctan((y - y0)/r_t), weighted least squares with bounds (v_c may be negative: the sign says which side recedes).
    ``inc_deg``: when given, v_c is deprojected: v_c,true = v_c,obs / (sin i cos psi) with psi = slit_offset_deg (slit-to-major-axis angle).
    -> dict(ok, vsys, vc_obs, rt, y0, errors, chi2r, dof, vc (deprojected or None), vflat_obs = |vc_obs|)"""
    from scipy.optimize import least_squares
    y = np.asarray(y, float); v = np.asarray(v, float); ev = np.asarray(ev, float)
    m = np.isfinite(y) & np.isfinite(v) & np.isfinite(ev) & (ev > 0)
    y, v, ev = y[m], v[m], ev[m]
    out = dict(ok=False, n=int(len(y)))
    npar = 2 + (1 if vsys_free else 0) + (1 if free_centre else 0) + 1
    if len(y) < max(4, npar):
        out['reason'] = 'too few points (%d)' % len(y)
        return out
    span = float(np.ptp(y)) or 1.0
    yc0 = float(y0) if y0 is not None else float(np.median(y))
    rng = np.random.default_rng(seed)
    names = (['vsys'] if vsys_free else []) + ['vc', 'rt'] + (['y0'] if free_centre else [])
    fixed = dict(vsys=0.0, y0=yc0)

    def unpack(p):
        d = dict(fixed)
        for k, val in zip(names, p):
            d[k] = val
        return d

    def resid(p):
        d = unpack(p)
        return (d['vsys'] + arctan_rc(y - d['y0'], d['vc'], d['rt']) - v) / ev

    best = None
    vr = float(np.ptp(v)) or 50.0
    lb = np.array([{'vsys': -np.inf, 'vc': -5e3, 'rt': 0.02 * span, 'y0': y.min()}[k] for k in names])
    ub = np.array([{'vsys': np.inf, 'vc': 5e3, 'rt': 3.0 * span, 'y0': y.max()}[k] for k in names])
    for i in range(n_starts):
        p0 = []
        for k in names:
            if k == 'vsys':
                p0.append(float(np.median(v)))
            elif k == 'vc':
                p0.append(vr * (0.5 if i % 2 == 0 else -0.5) * rng.uniform(0.6, 1.6))
            elif k == 'rt':
                p0.append(span * rng.uniform(0.05, 0.5))
            else:
                p0.append(yc0 + rng.normal(0, 0.05 * span))
        p0 = np.clip(p0, lb + 1e-9, ub - 1e-9)
        try:
            r = least_squares(resid, p0, bounds=(lb, ub), x_scale=np.maximum(np.abs(p0), 1.0))
        except Exception:
            continue
        if best is None or r.cost < best.cost:
            best = r
    if best is None:
        out['reason'] = 'fit failed'
        return out
    d = unpack(best.x)
    dof = max(len(y) - len(names), 1)
    chi2r = float(2 * best.cost / dof)
    try:
        cov = np.linalg.inv(best.jac.T @ best.jac) * max(chi2r, 1.0)
        errs = {k: float(math.sqrt(max(cov[i, i], 0.0))) for i, k in enumerate(names)}
    except np.linalg.LinAlgError:
        errs = {}
    out.update(ok=True, vsys=float(d['vsys']), vc_obs=float(d['vc']), rt=float(d['rt']), y0=float(d['y0']), errors=errs, chi2r=chi2r, dof=int(dof), vflat_obs=float(abs(d['vc'])))
    if inc_deg:
        f = math.sin(math.radians(inc_deg)) * math.cos(math.radians(slit_offset_deg))
        out['vc'] = float(d['vc'] / f) if abs(f) > 1e-3 else None
        out['vc_err'] = float(errs.get('vc', float('nan')) / abs(f)) if abs(f) > 1e-3 and errs else None
    return out


# ------------------------------------------------------------------------------------------------------------- IFU cube
def accretion_bins_2d(signal, noise, target, max_radius=6.0):
    """simple accretion binning (Cappellari & Copin 2003 step 1, without regularisation): start at the highest-S/N free spaxel, add the nearest free neighbour of the
    bin centroid until S/N = sum(s)/sqrt(sum(n^2)) >= target or the bin exceeds ``max_radius`` pixels.  -> bin-id map (-1 = not binned), list of index lists"""
    ny, nx = signal.shape
    sig = np.where(np.isfinite(signal), signal, 0.0)
    nz = np.where(np.isfinite(noise) & (noise > 0), noise, np.inf)
    ids = -np.ones((ny, nx), int)
    bins = []
    order = np.dstack(np.unravel_index(np.argsort(-(sig / nz).ravel()), (ny, nx)))[0]
    for (j0, i0) in order:
        if ids[j0, i0] >= 0 or sig[j0, i0] <= 0:
            continue
        mem = [(j0, i0)]
        s, v = sig[j0, i0], nz[j0, i0] ** 2
        ids[j0, i0] = len(bins)
        while s / math.sqrt(v) < target:
            cj = np.mean([m[0] for m in mem]); ci = np.mean([m[1] for m in mem])
            cand = None
            bd = 1e9
            for (j, i) in mem:
                for dj in (-1, 0, 1):
                    for di in (-1, 0, 1):
                        jj, ii = j + dj, i + di
                        if 0 <= jj < ny and 0 <= ii < nx and ids[jj, ii] < 0 and sig[jj, ii] > 0:
                            d = math.hypot(jj - cj, ii - ci)
                            if d < bd:
                                bd, cand = d, (jj, ii)
            if cand is None or bd > max_radius:
                break
            mem.append(cand)
            ids[cand] = len(bins)
            s += sig[cand]
            v += nz[cand] ** 2
        if s / math.sqrt(v) >= target:
            bins.append(mem)
        else:
            for m in mem:
                ids[m] = -1
    return ids, bins


def fit_cube_line(cube, wave, lam_sys, err=None, window_kms=900.0, snr_pix=3.0, bin_snr=0.0, inst_fwhm_A=0.0, sigma_guess_kms=80.0, gap_kms=1200.0, side_kms=2500.0, max_bin_radius=6.0):
    """Gaussian line fits per spaxel (or per accretion bin when bin_snr > 0) of an IFU cube cube[nw, ny, nx].
    Continuum per spaxel = median of the side bands; a spaxel is fitted when its matched-filter line S/N >= snr_pix (bin_snr = 0), or the bin reaches bin_snr.
    -> dict of maps [ny, nx]: flux, flux_err, vel, vel_err, sigma, sigma_err, snr, bin (id or -1) + counts"""
    w = np.asarray(wave, float)
    C = np.asarray(cube, float)
    nw, ny, nx = C.shape
    half = window_kms / C_KMS * lam_sys
    gap, side = gap_kms / C_KMS * lam_sys, side_kms / C_KMS * lam_sys
    flat = C.reshape(nw, -1).T                                           # [npix, nw]
    cont, sd, sidemask = _side_continuum(w, flat, lam_sys, half, gap, side)
    R = flat - cont[:, None]
    if err is not None:
        E = np.asarray(err, float).reshape(nw, -1).T
    else:
        E = np.broadcast_to(sd[:, None], flat.shape)
    inl = np.abs(w - lam_sys) <= half
    xw = w[inl]
    dlam = float(np.median(np.diff(w)))
    s0 = max(sigma_guess_kms / C_KMS * lam_sys, 1.2 * dlam)
    inst = inst_sigma_from_fwhm(inst_fwhm_A)
    # matched filter on the integrated profile (spatially summed line profile defines the template position)
    coll = np.nansum(R[:, inl], axis=0)
    lc = float(xw[np.argmax(coll)]) if np.nanmax(coll) > 0 else lam_sys
    # per-spaxel matched filter at the centroid of its own brightest pixel would bias the S/N; use a boxcar of +-2.5 s0 around the spaxel's moment centroid
    Rin = R[:, inl]
    pos = np.clip(Rin, 0, None)
    tot = pos.sum(axis=1)
    cen = np.where(tot > 0, (pos * xw[None]).sum(axis=1) / np.where(tot > 0, tot, 1), lc)
    box = np.abs(xw[None] - cen[:, None]) <= 2.5 * s0
    nbox = np.maximum(box.sum(axis=1), 1)
    sig_pix = np.nansum(np.where(box, Rin, 0.0), axis=1)
    nz_pix = np.sqrt(np.nansum(np.where(box, E[:, inl] ** 2, 0.0), axis=1))
    snr_map = np.where(nz_pix > 0, sig_pix / np.where(nz_pix > 0, nz_pix, np.inf), 0.0).reshape(ny, nx)
    maps = {k: np.full((ny, nx), np.nan) for k in ('flux', 'flux_err', 'vel', 'vel_err', 'sigma', 'sigma_err', 'snr')}
    binid = -np.ones((ny, nx), int)
    if bin_snr > 0:
        ids, bins = accretion_bins_2d(sig_pix.reshape(ny, nx), nz_pix.reshape(ny, nx), bin_snr, max_bin_radius)
        groups = [[j * nx + i for (j, i) in b] for b in bins]
        binid = ids
    else:
        good = np.where(snr_map.ravel() >= snr_pix)[0]
        groups = [[int(k)] for k in good]
    nfit = 0
    for gi, g in enumerate(groups):
        spec = np.nansum(Rin[g], axis=0)
        ee = np.sqrt(np.nansum(E[g][:, inl] ** 2, axis=0))
        # per-group window centred on its moment centroid keeps the line inside the window for large velocities
        pc = np.clip(spec, 0, None)
        c0 = float((pc * xw).sum() / pc.sum()) if pc.sum() > 0 else lc
        f = _fit_gauss(xw, spec, ee, c0, s0)
        if not f['ok'] or not np.isfinite(f['snr']) or f['snr'] < (0.6 * bin_snr if bin_snr > 0 else 0.6 * snr_pix) or f['flux'] <= 0:
            continue
        v, ve, s, se = _to_kms(f, lam_sys, inst)
        for k in g:
            j, i = divmod(k, nx)
            maps['flux'][j, i] = f['flux'] / len(g); maps['flux_err'][j, i] = f['flux_err'] / len(g)
            maps['vel'][j, i] = v; maps['vel_err'][j, i] = ve; maps['sigma'][j, i] = s; maps['sigma_err'][j, i] = se; maps['snr'][j, i] = f['snr']
        nfit += 1
    maps['bin'] = binid
    maps['snr_matched'] = snr_map
    maps['n_fits'] = nfit
    maps['n_spaxels'] = int(np.isfinite(maps['vel']).sum())
    maps['lam_sys'] = float(lam_sys)
    return maps


# ------------------------------------------------------------------------------------------------------------- disc model
def disc_velocity(x, y, x0, y0, vsys, vc, rt, pa_deg, inc_deg):
    """line-of-sight velocity of a thin disc with an arctan rotation curve at pixel (x, y)"""
    ph, inc = math.radians(pa_deg), math.radians(inc_deg)
    dx, dy = np.asarray(x, float) - x0, np.asarray(y, float) - y0
    xr = dx * math.cos(ph) + dy * math.sin(ph)                 # along the major axis (towards the receding side for vc > 0)
    yr = -dx * math.sin(ph) + dy * math.cos(ph)
    r = np.sqrt(xr ** 2 + (yr / max(math.cos(inc), 1e-3)) ** 2)
    cth = np.where(r > 0, xr / np.where(r > 0, r, 1.0), 0.0)
    return vsys + arctan_rc(r, vc, rt) * math.sin(inc) * cth


def fit_velocity_field(vel, vel_err, x0=None, y0=None, inc_fixed=None, pa0=None, free_centre=True, n_starts=18, seed=0, weights=None):
    """Fit the disc model to a velocity map (NaN = no data).  Free: v_sys, v_c, r_t, PA, inclination (unless inc_fixed), centre (unless free_centre False).
    The sign convention puts v_c >= 0 (PA = receding side).  -> dict(ok, vsys, vc, rt, pa, inc, x0, y0, errors, chi2r, dof, n, rms, corr_vc_inc, model [ny, nx])"""
    from scipy.optimize import least_squares
    V = np.asarray(vel, float)
    ny, nx = V.shape
    E = np.asarray(vel_err, float) if vel_err is not None else np.full_like(V, 10.0)
    ok = np.isfinite(V) & np.isfinite(E) & (E > 0)
    if weights is not None:
        ok &= np.isfinite(weights)
    out = dict(ok=False, n=int(ok.sum()))
    if ok.sum() < 20:
        out['reason'] = 'too few velocity pixels (%d)' % ok.sum()
        return out
    jj, ii = np.nonzero(ok)
    xs, ys, v, e = ii.astype(float), jj.astype(float), V[ok], np.maximum(E[ok], 1.0)
    cx0 = float(x0) if x0 is not None else float(np.mean(xs))
    cy0 = float(y0) if y0 is not None else float(np.mean(ys))
    rmax = float(np.hypot(xs - cx0, ys - cy0).max())
    # initial position angle from the linear gradient of the velocity field
    A = np.column_stack([np.ones_like(xs), xs - cx0, ys - cy0])
    coef = np.linalg.lstsq(A / e[:, None], v / e, rcond=None)[0]
    pa_init = math.degrees(math.atan2(coef[2], coef[1])) % 360 if pa0 is None else float(pa0)
    vspan = 0.5 * float(np.ptp(v))
    names = ['vsys', 'vc', 'rt', 'pa'] + ([] if inc_fixed else ['inc']) + (['x0', 'y0'] if free_centre else [])
    fixed = dict(x0=cx0, y0=cy0, inc=float(inc_fixed) if inc_fixed else 60.0)

    def unpack(p):
        d = dict(fixed)
        for k, val in zip(names, p):
            d[k] = val
        return d

    def resid(p):
        d = unpack(p)
        m = disc_velocity(xs, ys, d['x0'], d['y0'], d['vsys'], d['vc'], d['rt'], d['pa'], d['inc'])
        return (m - v) / e

    lb = np.array([{'vsys': np.median(v) - 3 * vspan - 50, 'vc': 0.0, 'rt': 0.1, 'pa': -720.0, 'inc': 5.0, 'x0': cx0 - 0.3 * rmax - 2, 'y0': cy0 - 0.3 * rmax - 2}[k] for k in names])
    ub = np.array([{'vsys': np.median(v) + 3 * vspan + 50, 'vc': 5e3, 'rt': 2.0 * rmax, 'pa': 720.0, 'inc': 85.0, 'x0': cx0 + 0.3 * rmax + 2, 'y0': cy0 + 0.3 * rmax + 2}[k] for k in names])
    rng = np.random.default_rng(seed)
    best = None
    for s in range(n_starts):
        p0 = []
        for k in names:
            if k == 'vsys':
                p0.append(float(coef[0]))
            elif k == 'vc':
                p0.append(vspan * rng.uniform(0.8, 2.5))
            elif k == 'rt':
                p0.append(rmax * rng.uniform(0.1, 0.6))
            elif k == 'pa':
                p0.append(pa_init + (0.0 if s == 0 else rng.uniform(-60, 60)))
            elif k == 'inc':
                p0.append((30.0, 50.0, 70.0)[s % 3])
            elif k == 'x0':
                p0.append(cx0)
            else:
                p0.append(cy0)
        p0 = np.clip(p0, lb + 1e-6, ub - 1e-6)
        try:
            r = least_squares(resid, p0, bounds=(lb, ub), x_scale=np.maximum(np.abs(p0), 1.0), max_nfev=400)
        except Exception:
            continue
        if best is None or r.cost < best.cost:
            best = r
    if best is None:
        out['reason'] = 'fit failed'
        return out
    d = unpack(best.x)
    d['pa'] = float(d['pa'] % 360.0)
    dof = max(len(v) - len(names), 1)
    chi2r = float(2 * best.cost / dof)
    errs, corr = {}, None
    try:
        cov = np.linalg.inv(best.jac.T @ best.jac) * max(chi2r, 1.0)
        errs = {k: float(math.sqrt(max(cov[i, i], 0.0))) for i, k in enumerate(names)}
        if 'inc' in names:
            a, b = names.index('vc'), names.index('inc')
            corr = float(cov[a, b] / math.sqrt(cov[a, a] * cov[b, b]))
    except np.linalg.LinAlgError:
        pass
    yy, xx = np.mgrid[0:ny, 0:nx]
    model = disc_velocity(xx, yy, d['x0'], d['y0'], d['vsys'], d['vc'], d['rt'], d['pa'], d['inc'])
    rms = float(np.sqrt(np.mean((model[ok] - V[ok]) ** 2)))
    out.update(ok=True, vsys=float(d['vsys']), vc=float(d['vc']), rt=float(d['rt']), pa=float(d['pa']), inc=float(d['inc']), x0=float(d['x0']), y0=float(d['y0']), errors=errs,
               chi2r=chi2r, dof=int(dof), rms=rms, corr_vc_inc=corr, model=model, inc_fixed=bool(inc_fixed), vsini=float(d['vc'] * math.sin(math.radians(d['inc']))))
    return out


def angle_diff(a, b, period=360.0):
    return (a - b + period / 2) % period - period / 2


def dispersion_summary(sigma, sigma_err, flux=None, snr_min=3.0):
    """flux-weighted mean / median intrinsic dispersion (km/s) over the pixels with sigma/err > snr_min"""
    s, e = np.asarray(sigma, float), np.asarray(sigma_err, float)
    ok = np.isfinite(s) & np.isfinite(e) & (e > 0) & (s / np.where(e > 0, e, 1) > snr_min)
    if not ok.any():
        return dict(n=0, mean=float('nan'), median=float('nan'))
    w = np.nan_to_num(np.asarray(flux, float))[ok] if flux is not None else np.ones(ok.sum())
    w = np.clip(w, 1e-12, None)
    return dict(n=int(ok.sum()), mean=float(np.sum(w * s[ok]) / np.sum(w)), median=float(np.median(s[ok])))
