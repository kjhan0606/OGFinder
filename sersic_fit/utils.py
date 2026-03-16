"""Utility functions for Sérsic fitting."""

import numpy as np


def sersic_bn(n):
    """Approximate bn for Sérsic profile.

    bn satisfies: gamma(2n, bn) = 0.5 * Gamma(2n)
    Approximation from Ciotti & Bertin (1999).
    """
    return 1.9992 * n - 0.3271 if n > 0.36 else 0.01


def ellipse_radius(x, y, xc, yc, ellip, theta):
    """Compute elliptical radius.

    Parameters
    ----------
    x, y : arrays, pixel coordinates
    xc, yc : float, center
    ellip : float, ellipticity (1 - b/a)
    theta : float, position angle (radians, from x-axis)

    Returns
    -------
    r : array, elliptical radius
    """
    dx = x - xc
    dy = y - yc
    cos_t = np.cos(theta)
    sin_t = np.sin(theta)

    x_rot = dx * cos_t + dy * sin_t
    y_rot = -dx * sin_t + dy * cos_t

    q = 1.0 - ellip  # axis ratio
    if q < 0.05:
        q = 0.05

    return np.sqrt(x_rot**2 + (y_rot / q)**2)
