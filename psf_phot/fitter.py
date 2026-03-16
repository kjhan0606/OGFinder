"""Single source PSF fitting."""

import numpy as np
from scipy.optimize import least_squares
from scipy.ndimage import shift


def fit_psf_single(cutout, psf, x0=0.0, y0=0.0, max_shift=3.0):
    """Fit a PSF model to a single source cutout.

    Model: amplitude * shifted_PSF + background

    Parameters
    ----------
    cutout : 2D array, source cutout
    psf : 2D array, normalized PSF image
    x0, y0 : float, initial sub-pixel offset
    max_shift : float, max centroid shift

    Returns
    -------
    dict with keys: amplitude, dx, dy, bg, flux, fluxerr, chi2
    """
    ny, nx = cutout.shape
    psf_sum = np.sum(psf)
    if psf_sum <= 0:
        return None

    # Normalize PSF
    psf_norm = psf / psf_sum

    # Initial guesses
    bg0 = np.median(cutout[:3, :].ravel())
    amp0 = np.sum(cutout - bg0)

    p0 = [amp0, x0, y0, bg0]

    def residuals(params):
        amp, dx, dy, bg = params
        model_psf = shift(psf_norm, [dy, dx], order=1, mode='constant', cval=0)
        # Trim to cutout size if needed
        py, px = model_psf.shape
        cy, cx = cutout.shape
        sy = min(py, cy)
        sx = min(px, cx)
        model = amp * model_psf[:sy, :sx] + bg
        return (cutout[:sy, :sx] - model).ravel()

    bounds_lower = [0, -max_shift, -max_shift, -np.inf]
    bounds_upper = [np.inf, max_shift, max_shift, np.inf]

    try:
        result = least_squares(residuals, p0,
                                bounds=(bounds_lower, bounds_upper),
                                method='trf', max_nfev=200)
    except Exception:
        return None

    amp, dx, dy, bg = result.x
    flux = amp * psf_sum
    chi2 = np.sum(result.fun**2) / max(1, len(result.fun) - 4)

    # Flux error from Jacobian
    try:
        J = result.jac
        cov = np.linalg.inv(J.T @ J) * chi2
        amp_err = np.sqrt(cov[0, 0])
        fluxerr = amp_err * psf_sum
    except (np.linalg.LinAlgError, ValueError):
        fluxerr = np.nan

    return {
        'amplitude': amp,
        'dx': dx,
        'dy': dy,
        'bg': bg,
        'flux': flux,
        'fluxerr': fluxerr,
        'chi2': chi2,
    }
