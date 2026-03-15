"""Deconvolution utility functions."""

import numpy as np
from scipy.signal import fftconvolve


def pad_psf(psf, image_shape):
    """Pad PSF to match image size, centered.

    Parameters
    ----------
    psf : 2D ndarray
        PSF array.
    image_shape : tuple
        Target (ny, nx) shape.

    Returns
    -------
    padded : 2D ndarray
        Zero-padded PSF with same shape as image.
    """
    ny, nx = image_shape
    py, px = psf.shape
    padded = np.zeros((ny, nx), dtype=psf.dtype)

    # Place PSF centered
    sy = ny // 2 - py // 2
    sx = nx // 2 - px // 2
    padded[sy:sy + py, sx:sx + px] = psf
    return padded


def psf_fft(psf, image_shape):
    """Compute FFT of PSF padded to image size.

    Parameters
    ----------
    psf : 2D ndarray
        PSF array.
    image_shape : tuple
        Target shape.

    Returns
    -------
    H : 2D complex ndarray
        FFT of the padded, centered PSF.
    """
    padded = pad_psf(psf, image_shape)
    # Roll to put center at (0,0) for correct FFT convolution
    padded = np.roll(padded, -(padded.shape[0] // 2), axis=0)
    padded = np.roll(padded, -(padded.shape[1] // 2), axis=1)
    return np.fft.rfft2(padded)


def convergence_metric(prev, current):
    """Compute relative change between iterations.

    Parameters
    ----------
    prev, current : 2D ndarray
        Successive estimates.

    Returns
    -------
    metric : float
        RMS relative change.
    """
    diff = current - prev
    norm = np.sqrt(np.mean(prev**2))
    if norm < 1e-12:
        return float('inf')
    return float(np.sqrt(np.mean(diff**2)) / norm)
