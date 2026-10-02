"""Background and noise modelling: source-free sky, large-scale background (mesh + polynomial, scattered light / gradients), flat-field residual diagnostics, rms maps,
pixel correlation (autocorrelation, correlation length), empirical aperture-noise law sigma_N = sigma_pix * alpha * N^beta (blank apertures; Casertano et al. 2000, Whitaker et al. 2019)
and the resulting photometric errors.  Pure functions (numpy / scipy / sep), JSON-serialisable results.

Conventions: `data` is the science image (any sky level); `mask` True = pixel excluded (sources, bad pixels); apertures are circular with radius r (px), N = pi r^2.
"""
import math

import numpy as np
from scipy import ndimage as ndi


def _sep():
    try:
        import sep
    except ImportError:
        import sep_pjw as sep
    return sep


def robust_std(a):
    a = np.asarray(a, float)
    a = a[np.isfinite(a)]
    if a.size < 5:
        return float('nan')
    m = np.median(a)
    return float(1.4826 * np.median(np.abs(a - m)))


def clipped_std(a, nsig=3.0, iters=5):
    a = np.asarray(a, float)
    a = a[np.isfinite(a)]
    for _ in range(iters):
        if a.size < 5:
            break
        m, s = np.median(a), np.std(a)
        keep = np.abs(a - m) < nsig * s
        if keep.all():
            break
        a = a[keep]
    return float(np.std(a)) if a.size > 4 else float('nan')


# ---------------------------------------------------------------------------------------------------------------- source masking / background
def source_mask(data, bad=None, thresh=2.0, minarea=8, grow=3.0, bw=64):
    """Boolean mask of detected sources (grown by `grow` px) plus bad pixels: two passes (rough background, then re-detect)."""
    sep = _sep()
    d = np.ascontiguousarray(np.nan_to_num(data, nan=0.0), np.float64)
    bad = np.zeros(d.shape, bool) if bad is None else np.asarray(bad, bool) | ~np.isfinite(data)
    m = bad.copy()
    for _ in range(2):
        bkg = sep.Background(d, mask=m.astype(np.uint8), bw=bw, bh=bw, fw=3, fh=3)
        sub = d - bkg.back()
        objs, seg = sep.extract(sub, thresh, err=bkg.globalrms, minarea=minarea, segmentation_map=True, mask=bad.astype(np.uint8))
        det = seg > 0
        if grow > 0:
            det = ndi.binary_dilation(det, structure=ndi.generate_binary_structure(2, 1), iterations=int(round(grow)))
        m = bad | det
    return m


def mesh_background(data, mask, bw=64, fw=3):
    """SEP mesh background and rms (maps), masked pixels ignored."""
    sep = _sep()
    d = np.ascontiguousarray(np.nan_to_num(data, nan=0.0), np.float64)
    b = sep.Background(d, mask=np.asarray(mask, np.uint8), bw=bw, bh=bw, fw=fw, fh=fw)
    return b.back(), b.rms()


def poly_surface(shape, coeffs, order):
    ny, nx = shape
    y, x = np.mgrid[:ny, :nx]
    u = 2.0 * x / max(nx - 1, 1) - 1.0
    v = 2.0 * y / max(ny - 1, 1) - 1.0
    out = np.zeros(shape)
    k = 0
    for d in range(order + 1):
        for i in range(d, -1, -1):
            out += coeffs[k] * u ** i * v ** (d - i)
            k += 1
    return out


def fit_poly_background(data, mask, order=2, step=8, nsig=3.0, iters=4):
    """Polynomial surface (total degree `order`) fitted to the unmasked, block-averaged sky with sigma clipping.  Returns (surface, coeffs, rms of the block means)."""
    ny, nx = data.shape
    by, bx = ny // step, nx // step
    d = np.where(mask, np.nan, data)[:by * step, :bx * step].reshape(by, step, bx, step)
    cnt = np.sum(np.isfinite(d), axis=(1, 3))
    with np.errstate(invalid='ignore'):
        mean = np.nanmean(d, axis=(1, 3))
    ok = (cnt >= 0.5 * step * step) & np.isfinite(mean)
    yy, xx = np.mgrid[:by, :bx]
    u = 2.0 * (xx * step + step / 2.0) / max(nx - 1, 1) - 1.0
    v = 2.0 * (yy * step + step / 2.0) / max(ny - 1, 1) - 1.0
    terms = [(i, dd - i) for dd in range(order + 1) for i in range(dd, -1, -1)]
    A = np.stack([u[ok] ** i * v[ok] ** j for i, j in terms], axis=1)
    y = mean[ok]
    sel = np.ones(y.size, bool)
    c = np.zeros(len(terms))
    for _ in range(iters):
        if sel.sum() < len(terms) + 3:
            break
        c, *_ = np.linalg.lstsq(A[sel], y[sel], rcond=None)
        r = y - A @ c
        s = robust_std(r[sel])
        new = np.abs(r) < nsig * max(s, 1e-12)
        if (new == sel).all():
            break
        sel = new
    rms = float(np.std((y - A @ c)[sel])) if sel.sum() > 3 else float('nan')
    return poly_surface(data.shape, c, order), c, rms


def flatness(data, mask, grid=8):
    """Block medians of the (sky-subtracted) unmasked pixels on a grid x grid layout: peak-to-peak, rms of the block means vs the expectation from the pixel noise,
    and the plane fitted to them.  Returns a dict + the grid of medians."""
    ny, nx = data.shape
    ys = np.linspace(0, ny, grid + 1).astype(int)
    xs = np.linspace(0, nx, grid + 1).astype(int)
    med = np.full((grid, grid), np.nan)
    err = np.full((grid, grid), np.nan)
    for j in range(grid):
        for i in range(grid):
            blk = data[ys[j]:ys[j + 1], xs[i]:xs[i + 1]]
            mk = mask[ys[j]:ys[j + 1], xs[i]:xs[i + 1]]
            v = blk[~mk]
            if v.size > 100:
                med[j, i] = np.median(v)
                err[j, i] = 1.2533 * robust_std(v) / math.sqrt(v.size)
    ok = np.isfinite(med)
    yy, xx = np.mgrid[:grid, :grid]
    A = np.stack([np.ones(ok.sum()), xx[ok] / (grid - 1.0) - 0.5, yy[ok] / (grid - 1.0) - 0.5], axis=1)
    c, *_ = np.linalg.lstsq(A, med[ok], rcond=None)
    res = med[ok] - A @ c
    chi2 = float(np.sum((med[ok] - np.median(med[ok])) ** 2 / err[ok] ** 2)) / max(ok.sum() - 1, 1)
    return dict(grid=grid, level=float(np.median(med[ok])), peak_to_peak=float(np.nanmax(med) - np.nanmin(med)), rms_blocks=float(np.nanstd(med)), mean_block_err=float(np.nanmean(err)),
                chi2_flat=chi2, plane_dx=float(c[1]), plane_dy=float(c[2]), plane_amp=float(math.hypot(c[1], c[2])), resid_rms=float(np.std(res)), medians=med.tolist())


# ---------------------------------------------------------------------------------------------------------------- pixel correlation
def autocorrelation(sub, mask, maxlag=12):
    """Normalised 2-D autocorrelation of the unmasked, mean-subtracted sky (masked pixels set to 0 and the pair counts corrected).  Returns the (2 maxlag+1)^2 array."""
    d = np.where(mask | ~np.isfinite(sub), 0.0, sub)
    w = (~mask & np.isfinite(sub)).astype(float)
    mu = d.sum() / max(w.sum(), 1.0)
    d = (d - mu) * w
    ny, nx = d.shape
    sh = (2 * ny, 2 * nx)
    F = np.fft.rfft2(d, sh)
    W = np.fft.rfft2(w, sh)
    cc = np.fft.irfft2(np.abs(F) ** 2, sh)
    cw = np.fft.irfft2(np.abs(W) ** 2, sh)
    acf = np.fft.fftshift(np.where(cw > 0.5, cc / np.maximum(cw, 1.0), 0.0))
    cy, cx = ny, nx
    c = acf[cy - maxlag:cy + maxlag + 1, cx - maxlag:cx + maxlag + 1]
    return c / c[maxlag, maxlag]


def correlation_summary(acf):
    """Correlation length (the lag where the radial mean ACF falls below 1/e... reported as sqrt(2 ln2)-equivalent FWHM of a Gaussian ACF fit) and nearest-neighbour correlations."""
    m = acf.shape[0] // 2
    yy, xx = np.indices(acf.shape)
    rr = np.hypot(xx - m, yy - m)
    prof = np.array([acf[(rr >= k - 0.5) & (rr < k + 0.5)].mean() for k in range(m + 1)])
    # Gaussian-ACF width from r=1 value: rho(1) = exp(-1/(2 s^2)) -> s; FWHM of the implied noise-kernel
    rho1 = float(np.clip(prof[1], 1e-6, 0.999))
    s_acf = math.sqrt(-1.0 / (2.0 * math.log(rho1)))
    corr_fwhm = 2.355 * s_acf / math.sqrt(2.0)       # an ACF of std s_acf comes from a kernel of std s_acf / sqrt(2)
    return dict(rho1=rho1, rho1_x=float(acf[m, m + 1]), rho1_y=float(acf[m + 1, m]), rho_diag=float(acf[m + 1, m + 1]), kernel_fwhm_px=float(corr_fwhm), profile=prof.tolist())


# ---------------------------------------------------------------------------------------------------------------- blank apertures
def blank_positions(mask, r, n, rng, margin=None, max_tries=60):
    """Random circle centres whose aperture (radius r, plus margin) contains no masked pixel and lies fully inside the image."""
    sep = _sep()
    ny, nx = mask.shape
    margin = r + 1.0 if margin is None else margin
    out = []
    tries = 0
    m8 = np.ascontiguousarray(mask, np.float64)
    while len(out) < n and tries < max_tries:
        k = max(4 * (n - len(out)), 200)
        x = rng.uniform(margin, nx - 1 - margin, k)
        y = rng.uniform(margin, ny - 1 - margin, k)
        bad, _, _ = sep.sum_circle(m8, x, y, r, subpix=1)
        good = bad < 0.5
        out.extend(zip(x[good], y[good]))
        tries += 1
    return np.array(out[:n]).reshape(-1, 2)


def aperture_noise(sub, mask, radii, n_aper=600, seed=1, local_rms=None):
    """Std of blank-aperture sums vs radius.  Returns dict(radii, N, sigma, sigma_sqrtN (sigma_pix sqrt(N)), n_used, sigma_pix)."""
    sep = _sep()
    rng = np.random.default_rng(seed)
    d = np.ascontiguousarray(np.nan_to_num(sub, nan=0.0), np.float64)
    sig_pix = robust_std(sub[~mask][:2000000] if (~mask).sum() > 2000000 else sub[~mask])
    out = dict(radii=[], N=[], sigma=[], n_used=[], sigma_pix=sig_pix)
    for r in radii:
        pos = blank_positions(mask, r, n_aper, rng)
        if len(pos) < 30:
            continue
        s, _, _ = sep.sum_circle(d, pos[:, 0], pos[:, 1], r, subpix=5)
        out['radii'].append(float(r)); out['N'].append(math.pi * r * r); out['sigma'].append(clipped_std(s, 4.0)); out['n_used'].append(int(len(pos)))
    return out


def fit_noise_law(N, sigma, sigma_pix):
    """sigma_N = sigma_pix * alpha * N^beta (least squares in log space on N >= 4).  Also the uncorrelated expectation beta = 0.5, alpha = 1."""
    N = np.asarray(N, float); s = np.asarray(sigma, float)
    ok = (N >= 4) & np.isfinite(s) & (s > 0)
    if ok.sum() < 3:
        return dict(alpha=1.0, beta=0.5, ok=False)
    A = np.stack([np.ones(ok.sum()), np.log(N[ok])], axis=1)
    c, *_ = np.linalg.lstsq(A, np.log(s[ok] / sigma_pix), rcond=None)
    return dict(alpha=float(math.exp(c[0])), beta=float(c[1]), ok=True)


def noise_law_sigma(N, law, sigma_pix):
    return sigma_pix * law['alpha'] * np.asarray(N, float) ** law['beta']


# ---------------------------------------------------------------------------------------------------------------- photometric errors
def flux_error(flux, N, sigma_pix, law=None, gain=None, correlated=True):
    """Aperture flux error: sky term (correlated law, or sigma_pix sqrt(N) when `correlated` is False) plus Poisson term of the source (flux / gain)."""
    N = np.asarray(N, float)
    if correlated and law is not None:
        sky = noise_law_sigma(N, law, sigma_pix)
    else:
        sky = sigma_pix * np.sqrt(N)
    var = sky ** 2
    if gain:
        var = var + np.clip(flux, 0, None) / gain
    return np.sqrt(var)


def sample_map(m, x, y):
    ny, nx = m.shape
    ix = np.clip(np.round(x).astype(int), 0, nx - 1)
    iy = np.clip(np.round(y).astype(int), 0, ny - 1)
    return m[iy, ix]
