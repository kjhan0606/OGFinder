"""Wiener deconvolution in Fourier domain."""

import numpy as np
from .utils import psf_fft


def wiener_deconvolve(image, psf, nsr=0.01):
    """Wiener deconvolution.

    Computes: F^{-1}[ H* / (|H|^2 + NSR) * G ]

    where H is the PSF OTF, G is the image spectrum, and NSR is the
    noise-to-signal ratio.

    Parameters
    ----------
    image : 2D ndarray
        Observed image.
    psf : 2D ndarray
        Point spread function.
    nsr : float
        Noise-to-signal power ratio. Higher values = more smoothing.

    Returns
    -------
    result : 2D ndarray
        Deconvolved image.
    """
    image = image.astype(np.float64)

    G = np.fft.rfft2(image)
    H = psf_fft(psf, image.shape)

    # Wiener filter: H* / (|H|^2 + NSR)
    H_conj = np.conj(H)
    H_abs2 = np.abs(H)**2
    W = H_conj / (H_abs2 + nsr)

    result = np.fft.irfft2(W * G, s=image.shape)
    return result.real
