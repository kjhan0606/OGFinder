"""PSF building via stacking methods."""

import numpy as np
from .psf_utils import extract_star_cutout, centroid_com, normalize_cutout, shift_to_center


def _collect_cutouts(data, positions, size=51):
    """Extract, center, and normalize star cutouts.

    Parameters
    ----------
    data : 2D ndarray
        Image data.
    positions : list of (x, y)
        Star positions (0-indexed).
    size : int
        Cutout size.

    Returns
    -------
    cutouts : list of 2D ndarray
        Centered, normalized cutouts.
    """
    cutouts = []
    for x, y in positions:
        c = extract_star_cutout(data, x, y, size)
        if c is None:
            continue
        dx, dy = centroid_com(c)
        c = shift_to_center(c, dx, dy)
        c = normalize_cutout(c)
        cutouts.append(c)
    return cutouts


def build_psf_median(data, positions, size=51):
    """Build PSF using pixel-by-pixel median of star cutouts.

    Parameters
    ----------
    data : 2D ndarray
        Image data.
    positions : list of (x, y)
        Star positions.
    size : int
        Cutout size (odd).

    Returns
    -------
    psf : 2D ndarray
        Median-stacked PSF, normalized to unit sum.
    n_used : int
        Number of stars used.
    """
    cutouts = _collect_cutouts(data, positions, size)
    if len(cutouts) == 0:
        raise ValueError("No valid star cutouts extracted")

    stack = np.array(cutouts)
    psf = np.median(stack, axis=0)
    psf = np.maximum(psf, 0.0)
    total = psf.sum()
    if total > 0:
        psf /= total
    return psf, len(cutouts)


def build_psf_mean(data, positions, size=51, sigma_clip=3.0):
    """Build PSF using sigma-clipped mean of star cutouts.

    Parameters
    ----------
    data : 2D ndarray
        Image data.
    positions : list of (x, y)
        Star positions.
    size : int
        Cutout size.
    sigma_clip : float
        Sigma for clipping outliers.

    Returns
    -------
    psf : 2D ndarray
        Mean-stacked PSF, normalized to unit sum.
    n_used : int
        Number of stars used.
    """
    cutouts = _collect_cutouts(data, positions, size)
    if len(cutouts) == 0:
        raise ValueError("No valid star cutouts extracted")

    stack = np.array(cutouts)
    med = np.median(stack, axis=0)
    mad = np.median(np.abs(stack - med[np.newaxis, :, :]), axis=0)
    mad = np.where(mad > 0, mad, 1e-10)
    sigma_scale = 1.4826  # MAD to sigma conversion

    mask = np.abs(stack - med[np.newaxis, :, :]) < sigma_clip * sigma_scale * mad[np.newaxis, :, :]
    counts = mask.sum(axis=0)
    counts = np.where(counts > 0, counts, 1)
    psf = (stack * mask).sum(axis=0) / counts

    psf = np.maximum(psf, 0.0)
    total = psf.sum()
    if total > 0:
        psf /= total
    return psf, len(cutouts)


def build_psf_epsf(data, positions, size=51, n_iters=10, oversample=2):
    """Build effective PSF (ePSF) via iterative sub-pixel refinement.

    Implements a simplified Anderson & King (2000) approach:
    1. Create initial PSF from median stack
    2. For each iteration:
       a. Re-measure centroids against current PSF
       b. Shift cutouts to sub-pixel grid
       c. Average on oversampled grid
       d. Downsample back

    Parameters
    ----------
    data : 2D ndarray
        Image data.
    positions : list of (x, y)
        Star positions.
    size : int
        Output PSF size.
    n_iters : int
        Number of refinement iterations.
    oversample : int
        Oversampling factor.

    Returns
    -------
    psf : 2D ndarray
        ePSF, normalized to unit sum.
    n_used : int
        Number of stars used.
    """
    from scipy.ndimage import zoom

    # Initial median PSF
    cutouts_raw = []
    for x, y in positions:
        c = extract_star_cutout(data, x, y, size)
        if c is not None:
            cutouts_raw.append((c, x, y))

    if len(cutouts_raw) == 0:
        raise ValueError("No valid star cutouts extracted")

    # Start with median stack
    normed = []
    for c, x, y in cutouts_raw:
        dx, dy = centroid_com(c)
        shifted = shift_to_center(c, dx, dy)
        normed.append(normalize_cutout(shifted))
    psf = np.median(np.array(normed), axis=0)

    osize = size * oversample

    for iteration in range(n_iters):
        # Oversample current PSF for comparison
        psf_over = zoom(psf, oversample, order=3)

        # Re-center each cutout against current PSF
        refined = []
        for c_orig, x, y in cutouts_raw:
            c = c_orig.copy()
            # Cross-correlate with PSF to find precise offset
            dx, dy = _cross_correlate_offset(normalize_cutout(c), psf)
            shifted = shift_to_center(c, dx, dy)
            refined.append(normalize_cutout(shifted))

        # Stack on oversampled grid
        over_stack = np.zeros((len(refined), osize, osize))
        for i, r in enumerate(refined):
            over_stack[i] = zoom(r, oversample, order=3)

        psf_over = np.median(over_stack, axis=0)
        psf_over = np.maximum(psf_over, 0.0)

        # Downsample back
        psf = zoom(psf_over, 1.0 / oversample, order=3)
        # Ensure correct size
        if psf.shape[0] != size or psf.shape[1] != size:
            psf = psf[:size, :size]

    psf = np.maximum(psf, 0.0)
    total = psf.sum()
    if total > 0:
        psf /= total
    return psf, len(cutouts_raw)


def _cross_correlate_offset(cutout, reference):
    """Find sub-pixel offset between cutout and reference via cross-correlation.

    Parameters
    ----------
    cutout, reference : 2D ndarray
        Same-sized arrays.

    Returns
    -------
    dx, dy : float
        Offset (column, row) to shift cutout to match reference.
    """
    from scipy.signal import fftconvolve

    # Cross-correlation
    ref_flip = reference[::-1, ::-1]
    cc = fftconvolve(cutout, ref_flip, mode='same')

    # Find peak
    cy, cx = np.unravel_index(np.argmax(cc), cc.shape)
    center_y, center_x = cc.shape[0] // 2, cc.shape[1] // 2

    # Sub-pixel refinement via quadratic fit around peak
    if 1 <= cy < cc.shape[0] - 1 and 1 <= cx < cc.shape[1] - 1:
        dy_sub = 0.5 * (cc[cy - 1, cx] - cc[cy + 1, cx]) / \
                 (cc[cy - 1, cx] - 2 * cc[cy, cx] + cc[cy + 1, cx] + 1e-12)
        dx_sub = 0.5 * (cc[cy, cx - 1] - cc[cy, cx + 1]) / \
                 (cc[cy, cx - 1] - 2 * cc[cy, cx] + cc[cy, cx + 1] + 1e-12)
    else:
        dy_sub, dx_sub = 0.0, 0.0

    dx = (cx - center_x) + dx_sub
    dy = (cy - center_y) + dy_sub
    return dx, dy
