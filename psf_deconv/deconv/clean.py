"""CLEAN deconvolution algorithm (radio astronomy standard)."""

import sys
import numpy as np
from scipy.signal import fftconvolve


def clean_deconvolve(image, psf, gain=0.1, niter=1000, threshold=0.0,
                     callback=None):
    """CLEAN deconvolution algorithm.

    Iteratively finds and subtracts the PSF at the brightest residual pixel.
    Produces a clean component model convolved with a clean beam (Gaussian).

    Parameters
    ----------
    image : 2D ndarray
        Observed (dirty) image.
    psf : 2D ndarray
        Point spread function (dirty beam).
    gain : float
        Loop gain (fraction of peak subtracted each iteration).
    niter : int
        Maximum number of iterations.
    threshold : float
        Stop when peak residual drops below this value.
        If 0, uses 3x the image noise level.
    callback : callable, optional
        Called as callback(iteration, residual) periodically.

    Returns
    -------
    restored : 2D ndarray
        Restored (clean) image.
    model : 2D ndarray
        Clean component model (delta functions).
    residual : 2D ndarray
        Final residual image.
    """
    image = image.astype(np.float64)
    psf = psf.astype(np.float64)

    ny, nx = image.shape
    py, px = psf.shape

    # Center of PSF
    pcy, pcx = py // 2, px // 2

    # Normalize PSF peak to 1
    psf_peak = psf[pcy, pcx]
    if psf_peak > 0:
        psf = psf / psf_peak

    residual = image.copy()
    model = np.zeros_like(image)

    # Auto-threshold: 3 * noise estimate
    if threshold <= 0:
        # Estimate noise from edges of the image
        border = max(10, min(ny // 10, nx // 10))
        edge_pixels = np.concatenate([
            residual[:border, :].ravel(),
            residual[-border:, :].ravel(),
            residual[:, :border].ravel(),
            residual[:, -border:].ravel(),
        ])
        threshold = 3.0 * np.std(edge_pixels)

    report_interval = max(1, niter // 20)

    for i in range(niter):
        # Find peak in residual
        peak_idx = np.argmax(np.abs(residual))
        peak_y, peak_x = np.unravel_index(peak_idx, residual.shape)
        peak_val = residual[peak_y, peak_x]

        if np.abs(peak_val) < threshold:
            print(f"  CLEAN converged at iteration {i + 1}, "
                  f"peak={peak_val:.4f} < threshold={threshold:.4f}",
                  file=sys.stderr)
            break

        # Add clean component
        component = gain * peak_val
        model[peak_y, peak_x] += component

        # Subtract shifted PSF from residual
        y_lo = max(0, peak_y - pcy)
        y_hi = min(ny, peak_y + py - pcy)
        x_lo = max(0, peak_x - pcx)
        x_hi = min(nx, peak_x + px - pcx)

        py_lo = pcy - (peak_y - y_lo)
        py_hi = pcy + (y_hi - peak_y)
        px_lo = pcx - (peak_x - x_lo)
        px_hi = pcx + (x_hi - peak_x)

        residual[y_lo:y_hi, x_lo:x_hi] -= component * psf[py_lo:py_hi, px_lo:px_hi]

        if (i + 1) % report_interval == 0:
            print(f"  CLEAN iteration {i + 1}/{niter}, "
                  f"peak={peak_val:.4f}", file=sys.stderr)
            if callback:
                callback(i, residual)

    # Restore: convolve model with clean beam + add residual
    clean_beam = _make_clean_beam(psf)
    restored = fftconvolve(model, clean_beam, mode='same') + residual

    return restored, model, residual


def _make_clean_beam(psf):
    """Create a Gaussian clean beam by fitting the PSF core.

    Parameters
    ----------
    psf : 2D ndarray
        Dirty beam (PSF).

    Returns
    -------
    beam : 2D ndarray
        Gaussian clean beam (same size as PSF).
    """
    py, px = psf.shape
    cy, cx = py // 2, px // 2

    # Estimate FWHM from PSF by finding half-max radius
    peak = psf[cy, cx]
    half_max = peak / 2.0

    # Find half-max along x and y axes
    def find_hwhm(profile, center):
        for r in range(1, len(profile) - center):
            if center + r < len(profile) and profile[center + r] < half_max:
                # Linear interpolation
                if r > 0:
                    f = (profile[center + r - 1] - half_max) / \
                        (profile[center + r - 1] - profile[center + r] + 1e-12)
                    return r - 1 + f
                return float(r)
        return float(min(len(profile) - center - 1, 5))

    hwhm_x = find_hwhm(psf[cy, :], cx)
    hwhm_y = find_hwhm(psf[:, cx], cy)

    sigma_x = hwhm_x / 1.1774  # HWHM to sigma
    sigma_y = hwhm_y / 1.1774

    sigma_x = max(sigma_x, 1.0)
    sigma_y = max(sigma_y, 1.0)

    yy, xx = np.mgrid[0:py, 0:px]
    beam = np.exp(-((xx - cx)**2 / (2 * sigma_x**2) +
                    (yy - cy)**2 / (2 * sigma_y**2)))
    beam /= beam.sum()
    return beam
