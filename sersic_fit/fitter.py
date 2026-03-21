"""2D Sérsic profile fitting via scipy.optimize.least_squares."""

import numpy as np
from scipy.optimize import least_squares
from .models import sersic_2d
from .config import SersicConfig


def fit_sersic(cutout, x0, y0, a, b, theta, flux, cfg=None):
    """Fit a 2D Sérsic profile to a source cutout.

    Parameters
    ----------
    cutout : 2D array
    x0, y0 : float, center relative to cutout
    a, b : float, semi-major/minor axes (initial guess for re)
    theta : float, position angle (radians)
    flux : float, total flux (initial guess)
    cfg : SersicConfig

    Returns
    -------
    dict with keys: n, re, Ie, ellip, theta, chi2 (or None on failure)
    """
    if cfg is None:
        cfg = SersicConfig()

    ny, nx = cutout.shape
    yy, xx = np.mgrid[0:ny, 0:nx]

    # Initial guesses
    ellip0 = max(0.0, min(0.9, 1.0 - b / max(a, 0.5)))
    re0 = max(cfg.re_min, np.sqrt(a * b))
    bg0 = np.median(cutout[:3, :].ravel())
    Ie0 = max(1.0, flux / (2.0 * np.pi * re0 * re0 * (1.0 - ellip0)))
    n0 = 1.5  # Start between disk (1) and elliptical (4)

    p0 = [Ie0, re0, n0, x0, y0, ellip0, theta, bg0]

    bounds_lower = [0, cfg.re_min, cfg.n_min, x0 - 5, y0 - 5, 0.0, -np.pi, -np.inf]
    bounds_upper = [np.inf, max(a * 5, 50), cfg.n_max, x0 + 5, y0 + 5, 0.95, np.pi, np.inf]

    # Check for C fast path: full LM in C
    try:
        from .fast import is_available, fit_sersic_2d_c
        use_fast = is_available()
    except ImportError:
        use_fast = False

    if use_fast:
        xx_flat = np.ascontiguousarray(xx.ravel(), dtype=np.float64)
        yy_flat = np.ascontiguousarray(yy.ravel(), dtype=np.float64)
        data_flat = np.ascontiguousarray(cutout.ravel(), dtype=np.float64)

        try:
            params_out, chi2 = fit_sersic_2d_c(
                data_flat, xx_flat, yy_flat, p0,
                bounds_lower, bounds_upper, max_iter=500
            )
        except Exception:
            return None

        Ie, re, n, xc, yc, ellip, theta_fit, bg = params_out
    else:
        def residuals(params):
            model = sersic_2d(params, xx, yy)
            return (cutout - model).ravel()

        try:
            result = least_squares(residuals, p0,
                                    bounds=(bounds_lower, bounds_upper),
                                    method='trf', max_nfev=500)
        except Exception:
            return None

        Ie, re, n, xc, yc, ellip, theta_fit, bg = result.x
        chi2 = np.sum(result.fun**2) / max(1, len(result.fun) - 8)

    return {
        'n': n,
        're': re,
        'Ie': Ie,
        'ellip': ellip,
        'theta': theta_fit * 180.0 / np.pi,
        'chi2': chi2,
    }
