"""PSF subtraction for crowded field iteration."""

import numpy as np
from scipy.ndimage import shift


def subtract_sources(data, psf, results):
    """Subtract fitted PSFs from image data.

    Parameters
    ----------
    data : 2D array (modified copy returned)
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
        flux = r.get('FLUX_CROWD', np.nan)
        if np.isnan(flux) or flux <= 0:
            continue

        x = r['X_CROWD'] - 1.0
        y = r['Y_CROWD'] - 1.0
        amp = flux / psf_sum

        px0 = int(x) - psf_hx
        py0 = int(y) - psf_hy
        px1 = px0 + psf.shape[1]
        py1 = py0 + psf.shape[0]

        if px0 < 0 or py0 < 0 or px1 > nx or py1 > ny:
            continue

        frac_dx = x - int(x)
        frac_dy = y - int(y)
        shifted_psf = shift(psf_norm, [frac_dy, frac_dx], order=1,
                             mode='constant', cval=0)
        residual[py0:py1, px0:px1] -= amp * shifted_psf

    return residual
