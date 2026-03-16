"""Gini coefficient for galaxy morphology."""

import numpy as np


def measure_gini(cutout, segmap_cutout=None):
    """Compute Gini coefficient.

    Gini = (2 * sum(i * x_i)) / (mean * n * (n-1)) - (n+1)/(n-1)

    Parameters
    ----------
    cutout : 2D array
    segmap_cutout : 2D bool array, optional mask (True = source pixel)

    Returns
    -------
    G : float, Gini coefficient (0 = uniform, 1 = all flux in one pixel)
    """
    if cutout is None or cutout.size == 0:
        return np.nan

    if segmap_cutout is not None:
        pixels = cutout[segmap_cutout].ravel()
    else:
        # Use pixels above background (> 0)
        pixels = cutout[cutout > 0].ravel()

    n = len(pixels)
    if n < 2:
        return np.nan

    pixels = np.sort(np.abs(pixels))
    mean_val = np.mean(pixels)
    if mean_val == 0:
        return np.nan

    indices = np.arange(1, n + 1)
    gini = (2.0 * np.sum(indices * pixels)) / (mean_val * n * (n - 1)) - (n + 1.0) / (n - 1.0)

    return np.clip(gini, 0.0, 1.0)
