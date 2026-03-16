"""M20: second-order moment of brightest 20% of flux."""

import numpy as np


def measure_m20(cutout, segmap_cutout=None):
    """Compute M20 statistic.

    M20 = log10(sum(fi * ri^2) for brightest 20% / sum(fi * ri^2) total)

    Parameters
    ----------
    cutout : 2D array
    segmap_cutout : 2D bool array, optional mask (True = source pixel)

    Returns
    -------
    M20 : float (typically -2.5 to 0)
    """
    if cutout is None or cutout.size == 0:
        return np.nan

    ny, nx = cutout.shape

    if segmap_cutout is not None:
        yy, xx = np.where(segmap_cutout)
        fluxes = cutout[segmap_cutout]
    else:
        yy, xx = np.where(cutout > 0)
        fluxes = cutout[cutout > 0]

    if len(fluxes) < 2:
        return np.nan

    total_flux = np.sum(fluxes)
    if total_flux <= 0:
        return np.nan

    # Flux-weighted center
    xc = np.sum(xx * fluxes) / total_flux
    yc = np.sum(yy * fluxes) / total_flux

    # Second-order moments
    ri2 = (xx - xc)**2 + (yy - yc)**2
    m_tot = np.sum(fluxes * ri2)

    if m_tot <= 0:
        return np.nan

    # Sort by flux (descending) and find brightest 20%
    sort_idx = np.argsort(fluxes)[::-1]
    cumflux = np.cumsum(fluxes[sort_idx])
    thresh = 0.2 * total_flux
    n_bright = np.searchsorted(cumflux, thresh) + 1
    n_bright = min(n_bright, len(fluxes))

    bright_idx = sort_idx[:n_bright]
    m_bright = np.sum(fluxes[sort_idx[:n_bright]] * ri2[sort_idx[:n_bright]])

    if m_bright <= 0:
        return np.nan

    return np.log10(m_bright / m_tot)
