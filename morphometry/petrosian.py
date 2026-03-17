"""Petrosian radius: where eta(r) = SB(r) / <SB(<r)> = 0.2."""

import numpy as np


def measure_petrosian(data, x, y, eta_thresh=0.2, rmax=100.0, nsteps=50):
    """Measure Petrosian radius.

    Parameters
    ----------
    data : 2D array, background-subtracted
    x, y : float, center (0-indexed)
    eta_thresh : float, Petrosian ratio threshold (default 0.2)
    rmax : float, maximum search radius (pixels)
    nsteps : int, number of radial bins

    Returns
    -------
    r_petro : float, Petrosian radius in pixels
    """
    try:
        import sep
    except ImportError:
        try:
            import sep_pjw as sep
        except ImportError:
            return np.nan

    if data.dtype != np.float64:
        data = data.astype(np.float64)

    radii = np.linspace(2.0, rmax, nsteps)

    prev_flux = 0.0
    for i, r in enumerate(radii):
        flux_r, _, _ = sep.sum_circle(data, [x], [y], r, subpix=5)
        flux_r = flux_r[0]

        if i == 0:
            prev_flux = flux_r
            continue

        # Mean SB within r
        area_r = np.pi * r * r
        mean_sb = flux_r / area_r if area_r > 0 else 0

        # Annular SB
        r_prev = radii[i - 1]
        annular_area = np.pi * (r * r - r_prev * r_prev)
        annular_flux = flux_r - prev_flux
        annular_sb = annular_flux / annular_area if annular_area > 0 else 0

        prev_flux = flux_r

        if mean_sb <= 0:
            continue

        eta = annular_sb / mean_sb

        if eta <= eta_thresh:
            # Interpolate
            return r

    return np.nan
