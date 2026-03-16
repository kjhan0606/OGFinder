"""Asymmetry index: A = min(sum|I - I_180|) / (2 * sum|I|)."""

import numpy as np
from scipy.ndimage import rotate


def measure_asymmetry(cutout, max_shift=3):
    """Measure rotational asymmetry of a cutout.

    Parameters
    ----------
    cutout : 2D array, background-subtracted cutout centered on source
    max_shift : int, maximum pixel shift for center optimization

    Returns
    -------
    A : float, asymmetry index (0 = symmetric, 1 = fully asymmetric)
    """
    if cutout is None or cutout.size == 0:
        return np.nan

    ny, nx = cutout.shape
    cy, cx = ny // 2, nx // 2

    sum_abs = np.nansum(np.abs(cutout))
    if sum_abs == 0:
        return np.nan

    best_a = np.inf

    for dy in range(-max_shift, max_shift + 1):
        for dx in range(-max_shift, max_shift + 1):
            # Shift center
            shifted = np.roll(np.roll(cutout, dy, axis=0), dx, axis=1)
            # Rotate 180 degrees
            rotated = rotate(shifted, 180.0, reshape=False, order=1)
            # Compute asymmetry
            residual = np.nansum(np.abs(shifted - rotated))
            a_val = residual / (2.0 * sum_abs)
            if a_val < best_a:
                best_a = a_val

    return best_a if np.isfinite(best_a) else np.nan
