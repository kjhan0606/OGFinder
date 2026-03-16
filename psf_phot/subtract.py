"""PSF subtraction for iterative cleaning."""

import numpy as np
from scipy.ndimage import shift


def subtract_psf(data, psf, results):
    """Subtract fitted PSFs from image.

    Parameters
    ----------
    data : 2D array, image data (modified in-place)
    psf : 2D array, PSF image
    results : list of fit result dicts

    Returns
    -------
    residual : 2D array
    """
    residual = data.copy()
    psf_sum = np.sum(psf)
    if psf_sum <= 0:
        return residual

    psf_norm = psf / psf_sum
    psf_hy, psf_hx = psf.shape[0] // 2, psf.shape[1] // 2
    ny, nx = data.shape

    for r in results:
        if np.isnan(r.get('FLUX_PSF', np.nan)):
            continue

        x = r['X_PSF'] - 1.0  # back to 0-indexed
        y = r['Y_PSF'] - 1.0
        amp = r['FLUX_PSF'] / psf_sum

        # Stamp region
        x0 = int(x) - psf_hx
        y0 = int(y) - psf_hy
        x1 = x0 + psf.shape[1]
        y1 = y0 + psf.shape[0]

        if x0 < 0 or y0 < 0 or x1 > nx or y1 > ny:
            continue

        dx = x - int(x)
        dy = y - int(y)
        shifted_psf = shift(psf_norm, [dy, dx], order=1, mode='constant', cval=0)
        residual[y0:y1, x0:x1] -= amp * shifted_psf

    return residual
