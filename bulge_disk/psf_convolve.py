"""PSF convolution utilities for Bulge+Disk decomposition."""

import sys
import numpy as np
from scipy.signal import fftconvolve


def convolve_model(model, psf):
    """Convolve a 2D model image with a PSF kernel.

    Parameters
    ----------
    model : 2D array
        Model image to convolve.
    psf : 2D array or None
        PSF kernel (should be normalized to sum=1).
        If None, model is returned unchanged.

    Returns
    -------
    convolved : 2D array, same shape as model
    """
    if psf is None:
        return model

    return fftconvolve(model, psf, mode='same')


def load_psf(path):
    """Load PSF from FITS file and normalize to unit sum.

    Parameters
    ----------
    path : str or None
        Path to PSF FITS file. If None, returns None with a warning.

    Returns
    -------
    psf : 2D array normalized to sum=1, or None
    """
    if path is None:
        print("WARNING: No PSF provided; fitting without PSF convolution",
              file=sys.stderr)
        return None

    from astropy.io import fits as pyfits
    with pyfits.open(path) as hdul:
        for hdu in hdul:
            if hdu.data is not None and hdu.data.ndim == 2:
                psf = hdu.data.astype(np.float64)
                total = np.sum(psf)
                if total > 0:
                    psf /= total
                return psf

    print("WARNING: Could not load PSF from {}; fitting without convolution"
          .format(path), file=sys.stderr)
    return None
