"""PSF utility functions: cutout extraction, centroiding, normalization."""

import numpy as np
from scipy.ndimage import shift as ndimage_shift


def extract_star_cutout(data, x, y, size=51):
    """Extract a square cutout centered on (x, y).

    Parameters
    ----------
    data : 2D ndarray
        Image data.
    x, y : float
        Center position (0-indexed, x=column, y=row).
    size : int
        Cutout size (odd number preferred).

    Returns
    -------
    cutout : 2D ndarray or None
        Cutout array, or None if too close to edge.
    """
    half = size // 2
    ny, nx = data.shape
    ix, iy = int(round(x)), int(round(y))

    if ix - half < 0 or ix + half >= nx or iy - half < 0 or iy + half >= ny:
        return None

    cutout = data[iy - half:iy + half + 1, ix - half:ix + half + 1].copy()
    return cutout.astype(np.float64)


def centroid_com(cutout):
    """Compute center-of-mass centroid of a cutout.

    Parameters
    ----------
    cutout : 2D ndarray
        Image cutout.

    Returns
    -------
    dx, dy : float
        Offset from geometric center (column, row).
    """
    cy, cx = cutout.shape[0] / 2.0, cutout.shape[1] / 2.0
    d = cutout - np.median(cutout)
    d = np.maximum(d, 0.0)
    total = d.sum()
    if total <= 0:
        return 0.0, 0.0

    yy, xx = np.mgrid[0:cutout.shape[0], 0:cutout.shape[1]]
    com_x = (xx * d).sum() / total
    com_y = (yy * d).sum() / total
    return com_x - cx, com_y - cy


def normalize_cutout(cutout):
    """Normalize cutout to unit total flux.

    Parameters
    ----------
    cutout : 2D ndarray
        Image cutout.

    Returns
    -------
    normalized : 2D ndarray
        Cutout with sum = 1.0.
    """
    bg = np.median(cutout)
    out = cutout - bg
    total = out.sum()
    if total > 0:
        out /= total
    return out


def shift_to_center(cutout, dx, dy):
    """Sub-pixel shift cutout so the star is centered.

    Parameters
    ----------
    cutout : 2D ndarray
        Image cutout.
    dx, dy : float
        Offset to shift (column, row).

    Returns
    -------
    shifted : 2D ndarray
        Shifted cutout.
    """
    return ndimage_shift(cutout, [-dy, -dx], order=3, mode='constant', cval=0.0)
