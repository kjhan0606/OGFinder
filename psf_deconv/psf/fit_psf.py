"""PSF fitting via analytical models (Gaussian, Moffat)."""

import numpy as np
from scipy.optimize import least_squares
from .psf_utils import extract_star_cutout, centroid_com, normalize_cutout, shift_to_center


def _gaussian_2d(params, x, y):
    """2D Gaussian model: amplitude, x0, y0, sigma_x, sigma_y, theta, bg."""
    amp, x0, y0, sx, sy, theta, bg = params
    cos_t = np.cos(theta)
    sin_t = np.sin(theta)
    dx = x - x0
    dy = y - y0
    a = cos_t**2 / (2 * sx**2) + sin_t**2 / (2 * sy**2)
    b = -np.sin(2 * theta) / (4 * sx**2) + np.sin(2 * theta) / (4 * sy**2)
    c = sin_t**2 / (2 * sx**2) + cos_t**2 / (2 * sy**2)
    return amp * np.exp(-(a * dx**2 + 2 * b * dx * dy + c * dy**2)) + bg


def _moffat_2d(params, x, y):
    """2D Moffat model: amplitude, x0, y0, alpha, beta, bg."""
    amp, x0, y0, alpha, beta, bg = params
    r2 = (x - x0)**2 + (y - y0)**2
    return amp * (1.0 + r2 / alpha**2)**(-beta) + bg


def _fit_single_star(cutout, model='gaussian'):
    """Fit a 2D model to a single star cutout.

    Returns fitted parameters.
    """
    ny, nx = cutout.shape
    yy, xx = np.mgrid[0:ny, 0:nx]
    x_flat = xx.ravel().astype(float)
    y_flat = yy.ravel().astype(float)
    z_flat = cutout.ravel().astype(float)

    cx, cy = nx / 2.0, ny / 2.0
    bg_est = np.median(cutout)
    amp_est = cutout.max() - bg_est

    if model == 'gaussian':
        p0 = [amp_est, cx, cy, 2.0, 2.0, 0.0, bg_est]
        bounds_lo = [0, cx - 5, cy - 5, 0.5, 0.5, -np.pi, -np.inf]
        bounds_hi = [np.inf, cx + 5, cy + 5, nx / 2, ny / 2, np.pi, np.inf]

        def residuals(p):
            return (_gaussian_2d(p, x_flat, y_flat) - z_flat)

    elif model == 'moffat':
        p0 = [amp_est, cx, cy, 3.0, 2.5, bg_est]
        bounds_lo = [0, cx - 5, cy - 5, 0.5, 1.0, -np.inf]
        bounds_hi = [np.inf, cx + 5, cy + 5, nx / 2, 10.0, np.inf]

        def residuals(p):
            return (_moffat_2d(p, x_flat, y_flat) - z_flat)

    else:
        raise ValueError(f"Unknown model: {model}")

    result = least_squares(residuals, p0, bounds=(bounds_lo, bounds_hi),
                           method='trf', max_nfev=1000)
    return result.x


def fit_gaussian_psf(data, positions, size=51):
    """Build PSF by fitting 2D Gaussians to stars and averaging parameters.

    Parameters
    ----------
    data : 2D ndarray
        Image data.
    positions : list of (x, y)
        Star positions.
    size : int
        Cutout/PSF size.

    Returns
    -------
    psf : 2D ndarray
        Gaussian PSF model, normalized to unit sum.
    n_used : int
        Number of stars successfully fitted.
    params : dict
        Averaged fit parameters.
    """
    fit_params = []
    for x, y in positions:
        c = extract_star_cutout(data, x, y, size)
        if c is None:
            continue
        try:
            p = _fit_single_star(c, model='gaussian')
            fit_params.append(p)
        except Exception:
            continue

    if len(fit_params) == 0:
        raise ValueError("No stars successfully fitted")

    params = np.median(np.array(fit_params), axis=0)
    # Generate model PSF with centered position
    half = size // 2
    yy, xx = np.mgrid[0:size, 0:size]
    # Use averaged sigma_x, sigma_y, theta but center at half, half
    model_params = [1.0, half, half, params[3], params[4], params[5], 0.0]
    psf = _gaussian_2d(model_params, xx.astype(float), yy.astype(float))
    psf = np.maximum(psf, 0.0)
    total = psf.sum()
    if total > 0:
        psf /= total

    param_dict = {
        'sigma_x': float(params[3]),
        'sigma_y': float(params[4]),
        'theta': float(params[5]),
        'fwhm_x': float(params[3] * 2.3548),
        'fwhm_y': float(params[4] * 2.3548),
    }
    return psf, len(fit_params), param_dict


def fit_moffat_psf(data, positions, size=51):
    """Build PSF by fitting 2D Moffat profiles to stars.

    Parameters
    ----------
    data : 2D ndarray
        Image data.
    positions : list of (x, y)
        Star positions.
    size : int
        Cutout/PSF size.

    Returns
    -------
    psf : 2D ndarray
        Moffat PSF model, normalized to unit sum.
    n_used : int
        Number of stars successfully fitted.
    params : dict
        Averaged fit parameters (alpha, beta, FWHM).
    """
    fit_params = []
    for x, y in positions:
        c = extract_star_cutout(data, x, y, size)
        if c is None:
            continue
        try:
            p = _fit_single_star(c, model='moffat')
            fit_params.append(p)
        except Exception:
            continue

    if len(fit_params) == 0:
        raise ValueError("No stars successfully fitted")

    params = np.median(np.array(fit_params), axis=0)
    half = size // 2
    yy, xx = np.mgrid[0:size, 0:size]
    model_params = [1.0, half, half, params[3], params[4], 0.0]
    psf = _moffat_2d(model_params, xx.astype(float), yy.astype(float))
    psf = np.maximum(psf, 0.0)
    total = psf.sum()
    if total > 0:
        psf /= total

    alpha = float(params[3])
    beta = float(params[4])
    fwhm = 2.0 * alpha * np.sqrt(2.0**(1.0 / beta) - 1.0)

    param_dict = {
        'alpha': alpha,
        'beta': beta,
        'fwhm': fwhm,
    }
    return psf, len(fit_params), param_dict
