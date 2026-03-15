"""Richardson-Lucy deconvolution: standard, accelerated, and TV-regularized."""

import sys
import numpy as np
from scipy.signal import fftconvolve


def richardson_lucy(image, psf, iterations=30, callback=None):
    """Standard Richardson-Lucy deconvolution.

    Iterative ML estimate for Poisson noise:
        x_{k+1} = x_k * conv(PSF^T, y / conv(PSF, x_k))

    Parameters
    ----------
    image : 2D ndarray
        Observed image.
    psf : 2D ndarray
        Point spread function (normalized to unit sum).
    iterations : int
        Number of iterations.
    callback : callable, optional
        Called as callback(iteration, estimate) each step.

    Returns
    -------
    estimate : 2D ndarray
        Deconvolved image.
    """
    image = np.maximum(image.astype(np.float64), 0.0)
    psf = psf.astype(np.float64)
    psf_flip = psf[::-1, ::-1]

    estimate = np.full_like(image, image.mean())
    eps = 1e-12

    for i in range(iterations):
        blurred = fftconvolve(estimate, psf, mode='same')
        ratio = image / np.maximum(blurred, eps)
        correction = fftconvolve(ratio, psf_flip, mode='same')
        estimate = np.maximum(estimate * correction, 0.0)

        if callback:
            callback(i, estimate)
        print(f"  RL iteration {i + 1}/{iterations}", file=sys.stderr)

    return estimate


def richardson_lucy_accelerated(image, psf, iterations=30, callback=None):
    """Accelerated Richardson-Lucy using Biggs-Andrews vector extrapolation.

    Converges 3-10x faster than standard RL by predicting the next
    estimate using momentum from previous iterations.

    Parameters
    ----------
    image : 2D ndarray
        Observed image.
    psf : 2D ndarray
        Point spread function.
    iterations : int
        Number of iterations.
    callback : callable, optional
        Called each iteration.

    Returns
    -------
    estimate : 2D ndarray
        Deconvolved image.
    """
    image = np.maximum(image.astype(np.float64), 0.0)
    psf = psf.astype(np.float64)
    psf_flip = psf[::-1, ::-1]

    estimate = np.full_like(image, image.mean())
    estimate_prev = estimate.copy()
    eps = 1e-12

    for i in range(iterations):
        # Standard RL step
        blurred = fftconvolve(estimate, psf, mode='same')
        ratio = image / np.maximum(blurred, eps)
        correction = fftconvolve(ratio, psf_flip, mode='same')
        estimate_new = np.maximum(estimate * correction, 0.0)

        # Biggs-Andrews acceleration (from iteration 2 onward)
        if i >= 1:
            # Vector extrapolation
            diff_new = estimate_new - estimate
            diff_old = estimate - estimate_prev

            # Compute acceleration factor
            dot_num = np.sum(diff_old * diff_new)
            dot_den = np.sum(diff_old * diff_old)
            if dot_den > eps:
                alpha = min(max(dot_num / dot_den, 0.0), 1.0)
            else:
                alpha = 0.0

            estimate_accel = estimate_new + alpha * (estimate_new - estimate)
            estimate_accel = np.maximum(estimate_accel, 0.0)

            estimate_prev = estimate.copy()
            estimate = estimate_accel
        else:
            estimate_prev = estimate.copy()
            estimate = estimate_new

        if callback:
            callback(i, estimate)
        print(f"  Accelerated RL iteration {i + 1}/{iterations}", file=sys.stderr)

    return estimate


def richardson_lucy_tv(image, psf, iterations=30, tv_lambda=0.001, callback=None):
    """Richardson-Lucy with Total Variation regularization.

    Adds TV regularization to suppress noise amplification:
        x_{k+1} = x_k / (1 + lambda * div(grad(x_k)/|grad(x_k)|))
                  * conv(PSF^T, y / conv(PSF, x_k))

    Parameters
    ----------
    image : 2D ndarray
        Observed image.
    psf : 2D ndarray
        Point spread function.
    iterations : int
        Number of iterations.
    tv_lambda : float
        TV regularization strength.
    callback : callable, optional
        Called each iteration.

    Returns
    -------
    estimate : 2D ndarray
        Deconvolved image.
    """
    image = np.maximum(image.astype(np.float64), 0.0)
    psf = psf.astype(np.float64)
    psf_flip = psf[::-1, ::-1]

    estimate = np.full_like(image, image.mean())
    eps = 1e-12

    for i in range(iterations):
        # Standard RL correction
        blurred = fftconvolve(estimate, psf, mode='same')
        ratio = image / np.maximum(blurred, eps)
        correction = fftconvolve(ratio, psf_flip, mode='same')

        # TV regularization term: divergence of normalized gradient
        tv_term = _tv_regularization(estimate, eps)
        denominator = 1.0 + tv_lambda * tv_term
        denominator = np.maximum(denominator, eps)

        estimate = np.maximum(estimate * correction / denominator, 0.0)

        if callback:
            callback(i, estimate)
        print(f"  TV-RL iteration {i + 1}/{iterations}", file=sys.stderr)

    return estimate


def _tv_regularization(image, eps=1e-12):
    """Compute TV regularization term: -div(grad(u)/|grad(u)|).

    Parameters
    ----------
    image : 2D ndarray
        Current estimate.
    eps : float
        Small constant for numerical stability.

    Returns
    -------
    tv : 2D ndarray
        TV regularization term.
    """
    # Gradients (forward differences)
    gy = np.zeros_like(image)
    gx = np.zeros_like(image)
    gy[:-1, :] = image[1:, :] - image[:-1, :]
    gx[:, :-1] = image[:, 1:] - image[:, :-1]

    # Gradient magnitude
    grad_mag = np.sqrt(gx**2 + gy**2 + eps**2)

    # Normalized gradients
    nx = gx / grad_mag
    ny = gy / grad_mag

    # Divergence (backward differences)
    div = np.zeros_like(image)
    div[:, 1:] += nx[:, 1:] - nx[:, :-1]
    div[:, 0] += nx[:, 0]
    div[1:, :] += ny[1:, :] - ny[:-1, :]
    div[0, :] += ny[0, :]

    return -div
