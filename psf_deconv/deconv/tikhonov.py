"""Tikhonov regularized deconvolution."""

import numpy as np
from .utils import psf_fft


def tikhonov_deconvolve(image, psf, lam=0.001):
    """Tikhonov regularized deconvolution.

    Computes: F^{-1}[ H* / (|H|^2 + lambda) * G ]

    Similar to Wiener but with a constant regularization parameter
    instead of noise-to-signal ratio.

    Parameters
    ----------
    image : 2D ndarray
        Observed image.
    psf : 2D ndarray
        Point spread function.
    lam : float
        Regularization parameter. Higher = smoother result.

    Returns
    -------
    result : 2D ndarray
        Deconvolved image.
    """
    image = image.astype(np.float64)

    G = np.fft.rfft2(image)
    H = psf_fft(psf, image.shape)

    H_conj = np.conj(H)
    H_abs2 = np.abs(H)**2
    T = H_conj / (H_abs2 + lam)

    result = np.fft.irfft2(T * G, s=image.shape)
    return result.real
