"""Concentration index: C = 5 * log10(r80 / r20)."""

import numpy as np


def measure_concentration(data, x, y, a, b, theta, rmax=50.0):
    """Measure concentration using circular flux radii.

    Parameters
    ----------
    data : 2D array, background-subtracted
    x, y : float, center (0-indexed)
    a, b : float, semi-major/minor axes
    theta : float, position angle in radians
    rmax : float, max radius for flux measurement

    Returns
    -------
    C : float, concentration index
    """
    try:
        import sep
    except ImportError:
        return np.nan

    if data.dtype != np.float64:
        data = data.astype(np.float64)

    # Get total flux first
    flux_tot, _, _ = sep.sum_circle(data, [x], [y], rmax, subpix=5)
    flux_tot = flux_tot[0]
    if flux_tot <= 0:
        return np.nan

    # Get r20 and r80
    fracs = np.array([0.2, 0.8])
    radii, _ = sep.flux_radius(data, [x], [y], [rmax], fracs, subpix=5)
    r20 = radii[0][0]
    r80 = radii[0][1]

    if r20 <= 0 or r80 <= 0 or not np.isfinite(r20) or not np.isfinite(r80):
        return np.nan

    return 5.0 * np.log10(r80 / r20)
