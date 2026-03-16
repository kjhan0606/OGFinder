"""Sérsic profile model: I(r) = Ie * exp(-bn * ((r/re)^(1/n) - 1))."""

import numpy as np
from .utils import sersic_bn, ellipse_radius


def sersic_2d(params, x, y):
    """Evaluate 2D Sérsic model.

    Parameters
    ----------
    params : [Ie, re, n, xc, yc, ellip, theta, bg]
    x, y : 2D arrays of pixel coordinates

    Returns
    -------
    model : 2D array
    """
    Ie, re, n, xc, yc, ellip, theta, bg = params

    r = ellipse_radius(x, y, xc, yc, ellip, theta)
    bn = sersic_bn(n)

    # Avoid overflow
    exponent = bn * ((r / max(re, 0.1))**(1.0 / max(n, 0.1)) - 1.0)
    exponent = np.clip(exponent, -50, 50)

    return Ie * np.exp(-exponent) + bg
