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


def extract_star_cutout_masked(data, x, y, size=201, sat_limit=60000.0):
    """Extract cutout with saturated pixels masked to NaN.

    Parameters
    ----------
    data : 2D ndarray
        Image data.
    x, y : float
        Center position (0-indexed).
    size : int
        Cutout size.
    sat_limit : float
        Saturation threshold in ADU.

    Returns
    -------
    cutout : 2D ndarray or None
        Cutout with saturated pixels set to NaN, or None if too close to edge.
    """
    cutout = extract_star_cutout(data, x, y, size)
    if cutout is None:
        return None
    cutout = cutout.copy()
    cutout[cutout > sat_limit] = np.nan
    return cutout


def radial_blend_weight(size, r_inner, r_outer):
    """2D cosine blend weight map for core/wing PSF blending.

    Returns weight=1 inside r_inner (core region),
    cosine taper from 1→0 between r_inner and r_outer,
    weight=0 outside r_outer (wing region).

    Parameters
    ----------
    size : int
        Output array size (size x size).
    r_inner : float
        Inner blend radius (pixels).
    r_outer : float
        Outer blend radius (pixels).

    Returns
    -------
    weight : 2D ndarray
        Weight map (size x size), values in [0, 1].
    """
    center = size / 2.0
    yy, xx = np.mgrid[0:size, 0:size]
    r = np.sqrt((xx - center) ** 2 + (yy - center) ** 2)

    weight = np.ones((size, size), dtype=np.float64)
    # Transition zone: cosine taper
    mask = (r >= r_inner) & (r <= r_outer)
    weight[mask] = 0.5 * (1.0 + np.cos(np.pi * (r[mask] - r_inner) / (r_outer - r_inner)))
    # Outside outer radius: wing only
    weight[r > r_outer] = 0.0

    return weight
