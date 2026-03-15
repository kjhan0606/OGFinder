"""Maximum Entropy Method (MEM) deconvolution."""

import sys
import numpy as np
from scipy.signal import fftconvolve


def mem_deconvolve(image, psf, lam=0.1, niter=100, callback=None):
    """Maximum Entropy deconvolution.

    Maximizes the entropy S = -sum(x * log(x/m)) subject to the
    chi-squared constraint, using a multiplicative RL-like update
    with entropy regularization.

    Parameters
    ----------
    image : 2D ndarray
        Observed image.
    psf : 2D ndarray
        Point spread function.
    lam : float
        Entropy regularization strength.
        Smaller = less smoothing (closer to RL).
        Larger = more smoothing (stronger entropy prior).
    niter : int
        Maximum number of iterations.
    callback : callable, optional
        Called as callback(iteration, estimate) periodically.

    Returns
    -------
    estimate : 2D ndarray
        Deconvolved image.
    """
    image = np.maximum(image.astype(np.float64), 0.0)
    psf = psf.astype(np.float64)
    psf_flip = psf[::-1, ::-1]
    eps = 1e-12

    # Default model: flat image at mean level
    img_mean = max(image.mean(), eps)
    model = np.full_like(image, img_mean)

    # Initialize with RL-like start (better than flat)
    estimate = image.copy()
    estimate = np.maximum(estimate, eps)

    report_interval = max(1, niter // 20)

    for i in range(niter):
        # RL-style data fidelity update
        predicted = fftconvolve(estimate, psf, mode='same')
        ratio = image / np.maximum(predicted, eps)
        rl_correction = fftconvolve(ratio, psf_flip, mode='same')

        # Entropy regularization: multiplicative damping toward model
        # entropy_factor = (model/estimate)^lam → pulls estimate toward model
        entropy_factor = np.power(
            np.maximum(model, eps) / np.maximum(estimate, eps), lam)

        # Combined multiplicative update
        estimate = estimate * rl_correction * entropy_factor
        estimate = np.maximum(estimate, eps)

        if (i + 1) % report_interval == 0:
            residual = image - predicted
            chi2 = np.sum(residual**2)
            log_ratio = np.log(np.maximum(estimate, eps) /
                               np.maximum(model, eps))
            entropy = -np.sum(estimate * log_ratio)
            print(f"  MEM iteration {i + 1}/{niter}, chi2={chi2:.2e}, "
                  f"entropy={entropy:.2e}", file=sys.stderr)
            if callback:
                callback(i, estimate)

    return estimate
