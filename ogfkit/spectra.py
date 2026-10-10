"""Spectra: 1D / 2D / IFU-cube readers and extraction, continuum, line detection, Gaussian line fits, redshift from emission lines.

Pure numpy/scipy/astropy functions, JSON-serialisable results.  Wavelengths in Angstrom (vacuum) internally; flux in whatever unit the file has.
"""
import math

import numpy as np

SQRT2PI = math.sqrt(2 * math.pi)
FWHM = 2.3548200450309493

# rest wavelengths (vacuum, Angstrom) and a rough strength weight used when matching detected lines to a redshift
LINES = [
    ('Lya', 1215.67, 1.0), ('CIV', 1549.06, 0.5), ('HeII', 1640.42, 0.2), ('CIII]', 1908.73, 0.5), ('MgII', 2798.0, 0.6),
    ('[OII]', 3727.4, 0.9), ('[NeIII]', 3869.86, 0.3), ('Hd', 4102.89, 0.2), ('Hg', 4341.69, 0.3), ('Hb', 4862.68, 0.6),
    ('[OIII]4960', 4960.30, 0.5), ('[OIII]5008', 5008.24, 1.0), ('HeI', 5877.25, 0.2), ('[NII]6550', 6549.86, 0.2), ('Ha', 6564.61, 1.0),
    ('[NII]6585', 6585.27, 0.4), ('[SII]6718', 6718.29, 0.4), ('[SII]6733', 6732.67, 0.4), ('[SIII]9069', 9071.1, 0.3), ('[SIII]9532', 9533.2, 0.3),
    ('Pa-b', 12821.6, 0.3), ('Pa-a', 18756.4, 0.3),
]
LINE_DICT = {n: w for n, w, _ in LINES}

_UNIT_TO_A = {'angstrom': 1.0, 'a': 1.0, 'ang': 1.0, 'nm': 10.0, 'um': 1e4, 'micron': 1e4, 'microns': 1e4, 'mum': 1e4, 'm': 1e10}


def unit_factor(unit):
    u = (unit or 'angstrom').strip().lower().replace(' ', '')
    for k in ('angstroms', 'angstrom'):
        u = u.replace(k, 'angstrom')
    return _UNIT_TO_A.get(u, 1.0)


# ------------------------------------------------------------------ readers

def _wcs_axis(hdr, axis, n):
    """Wavelength (Angstrom) of the n pixels of FITS axis `axis` (1-based) from CRVAL/CDELT/CD/CRPIX, with CUNIT conversion; log axes (CTYPE ...-LOG) supported."""
    cr = float(hdr.get('CRVAL%d' % axis, 1.0)); cp = float(hdr.get('CRPIX%d' % axis, 1.0))
    cd = hdr.get('CD%d_%d' % (axis, axis), hdr.get('CDELT%d' % axis, 1.0))
    pix = np.arange(n) + 1.0
    ctype = str(hdr.get('CTYPE%d' % axis, ''))
    if ctype.upper().endswith('LOG'):
        w = cr * np.exp(float(cd) / cr * (pix - cp))
    else:
        w = cr + float(cd) * (pix - cp)
    unit = hdr.get('CUNIT%d' % axis, None)
    if unit is None and w.max() < 100 and w.min() > 0.05:
        unit = 'um'                                          # microns without CUNIT
    return w * unit_factor(unit)


def read_spectrum(path, wave_unit=None):
    """-> dict(wave, flux, err (or None)).  Accepts FITS tables (wave/flux/err columns by common names), 1D FITS images (WCS or 'WAVELENGTH' column),
    ASCII/CSV with 2 or 3 numeric columns (wave flux [err])."""
    p = str(path)
    if p.lower().endswith(('.fits', '.fit', '.fits.gz', '.fz')):
        from astropy.io import fits
        with fits.open(p) as hl:
            for h in hl:
                if hasattr(h, 'columns') and h.data is not None and len(h.columns) >= 2:
                    names = {c.name.lower(): c.name for c in h.columns}
                    wk = next((names[k] for k in ('wave', 'wavelength', 'lambda', 'loglam', 'lam', 'wav') if k in names), None)
                    fk = next((names[k] for k in ('flux', 'f_lambda', 'fnu', 'flam', 'spec') if k in names), None)
                    ek = next((names[k] for k in ('err', 'error', 'flux_err', 'sigma', 'noise', 'fluxerr', 'ferr') if k in names), None)
                    ik = next((names[k] for k in ('ivar', 'invvar') if k in names), None)
                    if wk and fk:
                        w = np.asarray(h.data[wk], float).ravel(); f = np.asarray(h.data[fk], float).ravel()
                        if wk.lower() == 'loglam':
                            w = 10 ** w
                        e = np.asarray(h.data[ek], float).ravel() if ek else None
                        if e is None and ik:
                            iv = np.asarray(h.data[ik], float).ravel()
                            e = np.where(iv > 0, 1.0 / np.sqrt(np.where(iv > 0, iv, 1.0)), np.nan)
                        unit = None
                        try:
                            unit = h.columns[wk].unit
                        except Exception:
                            pass
                        fac = unit_factor(wave_unit or (str(unit) if unit else None))
                        return dict(wave=w * fac, flux=f, err=e)
            for h in hl:
                if h.data is not None and getattr(h.data, 'ndim', 0) == 1:
                    f = np.asarray(h.data, float)
                    w = _wcs_axis(h.header, 1, len(f))
                    if wave_unit:
                        w = w / unit_factor(h.header.get('CUNIT1')) * unit_factor(wave_unit)
                    e = None
                    for g in hl:
                        if g is not h and g.data is not None and g.data.ndim == 1 and len(g.data) == len(f) and str(g.header.get('EXTNAME', '')).upper() in ('ERR', 'ERROR', 'SIGMA', 'NOISE'):
                            e = np.asarray(g.data, float)
                    return dict(wave=w, flux=f, err=e)
        raise ValueError('no spectrum found in %s' % p)
    rows = []
    with open(p) as fh:
        for line in fh:
            t = line.replace(',', ' ').replace('\t', ' ').split()
            if not t or t[0].startswith('#'):
                continue
            try:
                rows.append([float(v) for v in t[:3]])
            except ValueError:
                continue
    if not rows:
        raise ValueError('no numeric rows in %s' % p)
    n = min(len(r) for r in rows)
    a = np.array([r[:n] for r in rows])
    return dict(wave=a[:, 0] * unit_factor(wave_unit), flux=a[:, 1], err=a[:, 2] if n > 2 else None)


def read_2d(path, wave_unit=None, ext=None):
    """2D spectrum: image with dispersion along axis 1 (FITS NAXIS1) and space along axis 2.  If NAXIS2 is the dispersion axis (CTYPE2 wavelength-like) it is transposed.
    Optional 'ERR'/'VAR' extension of the same shape.  -> dict(wave, data[ny, nw], err or None)."""
    from astropy.io import fits
    with fits.open(path) as hl:
        h0 = next(h for h in hl if h.data is not None and h.data.ndim == 2)
        data = np.asarray(h0.data, float)
        hdr = h0.header
        err = None
        for h in hl:
            if h is not h0 and h.data is not None and h.data.shape == data.shape:
                nm = str(h.header.get('EXTNAME', '')).upper()
                if nm in ('ERR', 'ERROR', 'SIGMA', 'NOISE'):
                    err = np.asarray(h.data, float)
                elif nm in ('VAR', 'VARIANCE'):
                    err = np.sqrt(np.maximum(np.asarray(h.data, float), 0))
        ct1, ct2 = str(hdr.get('CTYPE1', '')).upper(), str(hdr.get('CTYPE2', '')).upper()
        if ('WAVE' in ct2 or 'AWAV' in ct2) and 'WAVE' not in ct1 and 'AWAV' not in ct1:
            data = data.T; err = None if err is None else err.T
            w = _wcs_axis(hdr, 2, data.shape[1])
        else:
            w = _wcs_axis(hdr, 1, data.shape[1])
    if wave_unit:
        w = w / 1.0 * unit_factor(wave_unit)
    return dict(wave=w, data=data, err=err)


def read_cube(path):
    """IFU cube (NAXIS3 = wavelength, then y, x): -> dict(wave, data[nw, ny, nx], err or None)."""
    from astropy.io import fits
    with fits.open(path) as hl:
        h0 = next(h for h in hl if h.data is not None and h.data.ndim == 3)
        data = np.asarray(h0.data, float)
        w = _wcs_axis(h0.header, 3, data.shape[0])
        err = None
        for h in hl:
            if h is not h0 and h.data is not None and h.data.shape == data.shape and str(h.header.get('EXTNAME', '')).upper() in ('ERR', 'ERROR', 'SIGMA', 'NOISE'):
                err = np.asarray(h.data, float)
    return dict(wave=w, data=data, err=err)


# ------------------------------------------------------------------ 2D / cube extraction

def extract_2d(data, center, half_width=3.0, err=None, sky=None, mode='boxcar', profile_sigma=None):
    """Collapse a 2D spectrum (rows = space) around row `center` (0-based).  sky = (inner, outer) offsets in rows: median of rows with inner<|dy|<=outer is subtracted
    per column.  mode boxcar = sum over |dy| <= half_width (fractional edge weights); optimal = Horne weights with a Gaussian profile of sigma `profile_sigma`
    (default measured from the column-collapsed profile).  Returns (flux, err, info)."""
    ny, nw = data.shape
    rows = np.arange(ny, dtype=float)
    d = np.array(data, float)
    var = np.full_like(d, np.nan) if err is None else np.asarray(err, float) ** 2
    skyv = np.zeros(nw)
    if sky:
        m = (np.abs(rows - center) > sky[0]) & (np.abs(rows - center) <= sky[1])
        if m.sum() >= 2:
            skyv = np.nanmedian(d[m], axis=0)
            d = d - skyv[None, :]
    if mode == 'optimal':
        prof = np.nansum(np.where(np.isfinite(d), d, 0), axis=1)
        if profile_sigma is None:
            from scipy.optimize import curve_fit
            m = np.abs(rows - center) <= max(3 * half_width, 6)
            yy = rows[m]
            try:
                popt, _ = curve_fit(lambda y, A, s_, b: A * np.exp(-0.5 * ((y - center) / s_) ** 2) + b, yy, prof[m], p0=[np.max(prof[m]), max(half_width / 2, 1.0), 0.0],
                                    bounds=([0, 0.4, -np.inf], [np.inf, 3 * half_width, np.inf]))
                profile_sigma = float(popt[1])
            except Exception:
                profile_sigma = float(max(half_width / 2, 1.0))
        P = np.exp(-0.5 * ((rows - center) / profile_sigma) ** 2)
        P[np.abs(rows - center) > 3 * profile_sigma + 1] = 0
        P /= P.sum()
        if err is None:                                      # variance from the scatter of the sky rows
            v = np.nanvar(d[(np.abs(rows - center) > 3 * profile_sigma + 1)], axis=0) if (np.abs(rows - center) > 3 * profile_sigma + 1).sum() > 2 else np.ones(nw)
            var = np.broadcast_to(v[None, :], d.shape).copy()
        var = np.where(np.isfinite(var) & (var > 0), var, np.nanmedian(var[np.isfinite(var)]) if np.isfinite(var).any() else 1.0)
        num = np.nansum(P[:, None] * d / var, axis=0)
        den = np.nansum(P[:, None] ** 2 / var, axis=0)
        flux = num / np.maximum(den, 1e-30)
        e = 1.0 / np.sqrt(np.maximum(den, 1e-30))
        return flux, e, dict(profile_sigma=float(profile_sigma), sky=skyv)
    wgt = np.clip(half_width + 0.5 - np.abs(rows - center), 0.0, 1.0)
    flux = np.nansum(wgt[:, None] * d, axis=0)
    if err is None:
        v = np.nanvar(d[(np.abs(rows - center) > half_width + 2)], axis=0) if (np.abs(rows - center) > half_width + 2).sum() > 2 else np.nanvar(d, axis=0)
        e = np.sqrt(np.sum(wgt ** 2) * v)
    else:
        e = np.sqrt(np.nansum((wgt[:, None] ** 2) * var, axis=0))
    return flux, e, dict(sky=skyv)


def extract_cube(cube, x, y, radius, err=None):
    """Circular-aperture spectrum of an IFU cube at (x, y) (0-based spaxel coordinates), fractional edge weights by 5x5 sub-pixel sampling.  Returns (flux, err)."""
    nw, ny, nx = cube.shape
    sub = 5
    gy, gx = np.mgrid[0:ny * sub, 0:nx * sub]
    w = (np.hypot((gx + 0.5) / sub - 0.5 - x, (gy + 0.5) / sub - 0.5 - y) <= radius).reshape(ny, sub, nx, sub).mean(axis=(1, 3))
    spec = np.nansum(cube * w[None], axis=(1, 2))
    if err is None:
        return spec, None
    return spec, np.sqrt(np.nansum((err ** 2) * (w[None] ** 2), axis=(1, 2)))


def white_light(cube, wave=None, wrange=None):
    m = np.ones(cube.shape[0], bool) if (wave is None or wrange is None) else (wave >= wrange[0]) & (wave <= wrange[1])
    return np.nansum(cube[m], axis=0)


def line_moments(cube, wave, lam0, half=15.0, cont_gap=15.0, cont_width=40.0, snr_min=3.0):
    """Flux, velocity (km/s) and dispersion (km/s) maps of an emission line in a cube from continuum-subtracted moments inside +-half Angstrom around lam0 (the observed
    wavelength).  Continuum = mean of two side bands (cont_gap..cont_gap+cont_width away).  Spaxels whose line flux is below snr_min x its rms are NaN."""
    inl = np.abs(wave - lam0) <= half
    side = ((np.abs(wave - lam0) > half + cont_gap) & (np.abs(wave - lam0) <= half + cont_gap + cont_width))
    cont = np.nanmean(cube[side], axis=0) if side.any() else np.zeros(cube.shape[1:])
    sd = np.nanstd(cube[side], axis=0) if side.sum() > 3 else np.full(cube.shape[1:], np.nan)
    sub = cube[inl] - cont[None]
    dl = np.gradient(wave)[inl]
    flux = np.nansum(sub * dl[:, None, None], axis=0)
    ferr = sd * np.sqrt(np.sum(dl ** 2)) if np.isfinite(sd).any() else np.full_like(flux, np.nan)
    ok = flux > snr_min * np.where(np.isfinite(ferr), ferr, np.inf)
    cen = np.nansum(sub * (wave[inl] * dl)[:, None, None], axis=0) / np.where(ok, flux, np.nan)
    var = np.nansum(sub * ((wave[inl][:, None, None] - cen[None]) ** 2 * dl[:, None, None]), axis=0) / np.where(ok, flux, np.nan)
    c = 299792.458
    vel = (cen / lam0 - 1) * c
    sig = np.sqrt(np.clip(var, 0, None)) / lam0 * c
    return dict(flux=np.where(ok, flux, np.nan), velocity=np.where(ok, vel, np.nan), sigma=np.where(ok, sig, np.nan))


# ------------------------------------------------------------------ continuum and noise

def continuum(wave, flux, width=101, niter=4, clip=2.5):
    """Running-median continuum with iterative masking of positive (emission) and strong negative outliers; interpolated through masked pixels."""
    f = np.asarray(flux, float)
    ok = np.isfinite(f)
    mask = ~ok
    from scipy.ndimage import median_filter
    width = int(max(5, width)) | 1
    cont = np.zeros_like(f)
    for _ in range(niter):
        g = np.where(mask, np.nan, f)
        # nan-aware running median via interpolation across masked pixels
        idx = np.arange(len(f))
        gi = np.interp(idx, idx[~mask], f[~mask]) if (~mask).sum() > 3 else np.zeros_like(f)
        cont = median_filter(gi, size=width, mode='nearest')
        res = f - cont
        s = 1.4826 * np.nanmedian(np.abs(res[~mask] - np.nanmedian(res[~mask]))) if (~mask).sum() > 3 else 1.0
        newmask = ~ok | (res > clip * s) | (res < -4 * s)
        # grow the mask by 2 px around emission
        newmask = np.convolve(newmask.astype(float), np.ones(5), mode='same') > 0
        if np.array_equal(newmask, mask):
            break
        mask = newmask
    return cont, mask


def noise_level(resid, mask=None):
    r = np.asarray(resid, float)
    ok = np.isfinite(r) if mask is None else np.isfinite(r) & ~mask
    if ok.sum() < 5:
        return float('nan')
    return float(1.4826 * np.median(np.abs(r[ok] - np.median(r[ok]))))


# ------------------------------------------------------------------ lines

def detect_lines(wave, flux, err=None, cont=None, kernel_sigma_px=2.0, snr_min=4.0, min_sep_px=4, sign=1.0):
    """Matched-filter (Gaussian) detection of emission lines.  sign=-1 detects absorption lines instead.
    Returns a list of dicts (wave, snr, index), strongest first; snr is reported positive for both."""
    from scipy.ndimage import maximum_filter1d
    f = np.asarray(flux, float)
    if cont is None:
        cont, m = continuum(wave, f)
    else:
        m = None
    res = np.nan_to_num(f - cont)
    if sign < 0:
        res = -res
    sig = noise_level(res, m) if err is None else float(np.nanmedian(err))
    h = int(math.ceil(4 * kernel_sigma_px))
    x = np.arange(-h, h + 1)
    g = np.exp(-0.5 * (x / kernel_sigma_px) ** 2)
    num = np.convolve(res, g, mode='same')
    if err is None:
        den = sig * math.sqrt(np.sum(g ** 2)) if sig > 0 else np.inf
        snr = num / den
    else:
        e = np.where(np.isfinite(err) & (err > 0), err, np.nanmedian(err))
        snr = np.convolve(res / e ** 2, g, mode='same') / np.sqrt(np.convolve(1.0 / e ** 2, g ** 2, mode='same'))
    mx = maximum_filter1d(snr, size=2 * min_sep_px + 1)
    idx = np.nonzero((snr == mx) & (snr >= snr_min))[0]
    out = [dict(index=int(i), wave=float(wave[i]), snr=float(snr[i])) for i in idx]
    out.sort(key=lambda d: -d['snr'])
    return out


def _gauss_lin(p, x, x0):
    c0, c1, A, mu, s = p
    return c0 + c1 * (x - x0) + A * np.exp(-0.5 * ((x - mu) / s) ** 2)


def fit_line(wave, flux, err, mu0, sigma0=None, window=None, fit_cont_slope=True, min_pts=8):
    """Gaussian + local linear continuum fit around mu0 (weighted least squares, errors from the covariance scaled by sqrt(chi2_red) if > 1).
    Returns dict(ok, mu, mu_err, sigma, sigma_err, flux, flux_err, amp, cont, fwhm, ew, chi2r, npts) (wave units as input)."""
    from scipy.optimize import least_squares
    w = np.asarray(wave, float); f = np.asarray(flux, float)
    dl = float(np.median(np.diff(w)))
    s0 = sigma0 if sigma0 else 2.5 * dl
    win = window if window else 6 * s0 + 4 * dl
    m = (np.abs(w - mu0) <= win) & np.isfinite(f)
    out = dict(ok=False, npts=int(m.sum()))
    if m.sum() < min_pts:
        out['reason'] = 'too few points'
        return out
    x, y = w[m], f[m]
    e = np.asarray(err, float)[m] if err is not None else np.full(m.sum(), max(noise_level(y - np.median(y)), 1e-30))
    e = np.where(np.isfinite(e) & (e > 0), e, np.nanmedian(e[e > 0]) if (e > 0).any() else 1.0)
    c0 = float(np.median(np.concatenate([y[:3], y[-3:]])))
    A0 = float(y[np.argmin(np.abs(x - mu0))] - c0)
    p0 = [c0, 0.0, A0 if A0 != 0 else 1e-3 * abs(c0 + 1), mu0, s0]
    lo = [-np.inf, -np.inf, -np.inf, mu0 - 3 * s0 - dl, 0.3 * dl]
    hi = [np.inf, np.inf, np.inf, mu0 + 3 * s0 + dl, win / 2]
    if not fit_cont_slope:
        lo[1], hi[1], p0[1] = -1e-30, 1e-30, 0.0
    try:
        r = least_squares(lambda p: (_gauss_lin(p, x, mu0) - y) / e, p0, bounds=(lo, hi))
    except Exception as ex:                                   # pragma: no cover
        out['reason'] = str(ex)
        return out
    p = r.x
    dof = max(len(x) - 5, 1)
    chi2r = float(np.sum(r.fun ** 2) / dof)
    try:
        cov = np.linalg.inv(r.jac.T @ r.jac) * max(chi2r, 1.0)
    except np.linalg.LinAlgError:
        cov = np.full((5, 5), np.nan)
    sd = np.sqrt(np.clip(np.diag(cov), 0, None))
    A, s = p[2], p[4]
    flux_ = A * s * SQRT2PI
    # flux = A s sqrt(2pi): propagate with the A-s covariance
    gA, gS = s * SQRT2PI, A * SQRT2PI
    var_f = gA ** 2 * cov[2, 2] + gS ** 2 * cov[4, 4] + 2 * gA * gS * cov[2, 4]
    cont_at = p[0] + p[1] * (p[3] - mu0)
    out.update(ok=True, mu=float(p[3]), mu_err=float(sd[3]), sigma=float(s), sigma_err=float(sd[4]), amp=float(A), cont=float(cont_at), flux=float(flux_), flux_err=float(math.sqrt(max(var_f, 0.0))),
               fwhm=float(FWHM * s), ew=float(flux_ / cont_at) if cont_at else float('nan'), chi2r=chi2r)
    return out


def estimate_redshift(det, wave_range, tol_px=3.0, dl=1.0, zmin=0.0, zmax=7.0, dz=0.0005, lines=LINES, snr_single=8.0, min_confidence=1.2):
    """Blind redshift from the detected emission lines.  For every trial z each rest line whose observed wavelength falls in the spectral range is matched to the
    nearest detection within tol_px*dl; score = sum over matched lines of weight x snr(capped at 30) - 0.5 x sum of the weights of strong (>=0.9) rest lines in range
    that were not detected.  Returns dict(z, score, n_match, matched [(name, rest, obs, snr)], quality 0/1/2/3 (0 none or ambiguous, 1 single line,
    2 two lines, 3 three or more), confidence = score / score of the best alternative (< min_confidence -> quality 0), runner_up (z, score))."""
    if not det:
        return dict(z=float('nan'), score=0.0, n_match=0, matched=[], quality=0)
    obs = np.array([d['wave'] for d in det]); snr = np.array([d['snr'] for d in det])
    zs = np.arange(zmin, zmax, dz)
    tol = tol_px * dl
    scores = np.zeros(len(zs))
    for iz, z in enumerate(zs):
        sc = 0.0
        for name, lam, wt in lines:
            lo = lam * (1 + z)
            if lo < wave_range[0] + tol or lo > wave_range[1] - tol:
                continue
            j = int(np.argmin(np.abs(obs - lo)))
            if abs(obs[j] - lo) <= tol * (1 + 0 * z):
                sc += wt * min(snr[j], 30.0)
            elif wt >= 0.9:
                sc -= 0.5 * wt * min(float(np.max(snr)), 30.0) * 0.2
        scores[iz] = sc
    iz = int(np.argmax(scores))
    z = float(zs[iz])
    # peak refinement + matched list
    matched = []
    for name, lam, wt in lines:
        lo = lam * (1 + z)
        if lo < wave_range[0] or lo > wave_range[1]:
            continue
        j = int(np.argmin(np.abs(obs - lo)))
        if abs(obs[j] - lo) <= tol:
            matched.append((name, lam, float(obs[j]), float(snr[j])))
    # runner-up: best score farther than 0.02(1+z) from the winner
    far = np.abs(zs - z) > 0.02 * (1 + z)
    ru = int(np.argmax(np.where(far, scores, -1e9)))
    n = len(matched)
    q = 0 if n == 0 else (1 if n == 1 else (2 if n == 2 else 3))
    if n == 1 and matched[0][3] < snr_single:
        q = 0
    conf = float(scores[iz] / max(scores[ru], 1e-9)) if scores[ru] > 0 else float('inf')
    if conf < min_confidence:                                # the best redshift is hardly better than an alternative: ambiguous
        q = 0
    return dict(z=z, score=float(scores[iz]), n_match=n, matched=matched, quality=q, confidence=conf, runner_up=dict(z=float(zs[ru]), score=float(scores[ru])))


def refine_redshift(matched_fits):
    """Weighted mean redshift from fitted lines: list of (rest, mu, mu_err) -> (z, z_err, chi2r)."""
    rest = np.array([m[0] for m in matched_fits]); mu = np.array([m[1] for m in matched_fits]); me = np.array([max(m[2], 1e-6) for m in matched_fits])
    zi = mu / rest - 1.0
    ze = me / rest
    w = 1.0 / ze ** 2
    z = float(np.sum(w * zi) / np.sum(w))
    err = float(1.0 / math.sqrt(np.sum(w)))
    chi2r = float(np.sum(w * (zi - z) ** 2) / max(len(zi) - 1, 1)) if len(zi) > 1 else float('nan')
    if len(zi) > 1 and chi2r > 1:
        err *= math.sqrt(chi2r)
    return z, err, chi2r


def analyse_spectrum(wave, flux, err=None, snr_min=4.0, kernel_sigma_px=2.0, cont_width=101, z_known=None, zmax=7.0, line_list=None):
    """Complete 1D analysis: continuum, line detection, Gaussian fits, redshift.  Returns a JSON-serialisable dict."""
    wave = np.asarray(wave, float); flux = np.asarray(flux, float)
    cont, mask = continuum(wave, flux, cont_width)
    res = flux - cont
    sig = noise_level(res, mask)
    e = err if err is not None else np.full_like(flux, sig)
    det = detect_lines(wave, flux, err, cont, kernel_sigma_px, snr_min)
    dl = float(np.median(np.diff(wave)))
    cont_snr = float(np.nanmedian(cont[~mask] / sig)) if (~mask).any() and sig > 0 else float('nan')
    out = dict(n_pix=int(len(wave)), dlambda=dl, noise=sig, cont_snr_per_pix=cont_snr, detections=det, cont_median=float(np.nanmedian(cont)))
    ll = LINES if line_list is None else line_list
    if z_known is None:
        zr = estimate_redshift(det, (float(wave[0]), float(wave[-1])), tol_px=max(3.0, kernel_sigma_px * 1.5), dl=dl, zmax=zmax, lines=ll)
    else:
        zr = dict(z=float(z_known), quality=4, n_match=0, matched=[], score=0.0)
        zr['matched'] = [(n, l, l * (1 + z_known), float('nan')) for n, l, _ in ll if wave[0] < l * (1 + z_known) < wave[-1]]
        zr['n_match'] = len(zr['matched'])
    fits = []
    for name, lam, mu_obs, s in zr['matched']:
        ft = fit_line(wave, flux - cont + np.nanmedian(cont), e, mu_obs, sigma0=kernel_sigma_px * dl)
        ft['name'] = name; ft['rest'] = lam
        fits.append(ft)
    good = [(f['rest'], f['mu'], f['mu_err']) for f in fits if f.get('ok') and f['mu_err'] > 0 and f['flux'] / max(f['flux_err'], 1e-30) > 2.5 and abs(f['mu'] - f['rest'] * (1 + zr['z'])) < 6 * dl]
    if good:
        z, ze, c2 = refine_redshift(good)
        zr.update(z=z, z_err=ze, z_chi2r=c2, n_lines_fit=len(good))
    else:
        zr.update(z_err=float('nan'), n_lines_fit=0)
    out.update(redshift=zr, lines=fits)
    return out
