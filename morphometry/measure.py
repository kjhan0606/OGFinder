"""Integrated morphometry measurement pipeline."""

import numpy as np
from .config import MorphometryConfig
from .concentration import measure_concentration
from .asymmetry import measure_asymmetry
from .gini import measure_gini
from .m20 import measure_m20
from .petrosian import measure_petrosian


def extract_cutout(data, x, y, halfsize):
    """Extract a square cutout centered on (x, y)."""
    ny, nx = data.shape
    x0 = max(0, int(x - halfsize))
    x1 = min(nx, int(x + halfsize + 1))
    y0 = max(0, int(y - halfsize))
    y1 = min(ny, int(y + halfsize + 1))

    if x1 <= x0 or y1 <= y0:
        return None

    return data[y0:y1, x0:x1].copy()


def measure_source(data, src, cfg=None):
    """Measure all non-parametric morphology indices for a source.

    Parameters
    ----------
    data : 2D array, background-subtracted
    src : dict with keys: x, y, a, b, theta, kron_radius, npix
    cfg : MorphometryConfig

    Returns
    -------
    dict with keys: CONC, ASYM, GINI, M20, R_PETRO
    """
    if cfg is None:
        cfg = MorphometryConfig()

    x = src['x']
    y = src['y']
    a = src['a']
    theta = src.get('theta', 0.0)

    result = {
        'CONC': np.nan,
        'ASYM': np.nan,
        'GINI': np.nan,
        'M20': np.nan,
        'R_PETRO': np.nan,
    }

    npix = src.get('npix', 0)
    if npix < cfg.min_npix:
        return result

    # Cutout size
    kron_r = src.get('kron_radius', 3.5)
    halfsize = int(cfg.cutout_scale * kron_r * max(a, 3.0))
    halfsize = max(halfsize, 20)
    halfsize = min(halfsize, 200)

    cutout = extract_cutout(data, x, y, halfsize)
    if cutout is None or cutout.size == 0:
        return result

    # Concentration
    rmax = min(3.0 * kron_r * a, 100.0)
    if rmax < 5.0:
        rmax = 5.0
    result['CONC'] = measure_concentration(data, x, y, a, src.get('b', a),
                                            theta, rmax=rmax)

    # Asymmetry
    result['ASYM'] = measure_asymmetry(cutout, max_shift=cfg.asym_center_maxshift)

    # Gini and M20
    result['GINI'] = measure_gini(cutout)
    result['M20'] = measure_m20(cutout)

    # Petrosian radius
    result['R_PETRO'] = measure_petrosian(data, x, y,
                                           eta_thresh=cfg.petro_eta,
                                           rmax=cfg.petro_rmax,
                                           nsteps=cfg.petro_nsteps)

    return result


def measure_catalog(data, catalog, cfg=None):
    """Measure morphometry for all sources in catalog.

    Parameters
    ----------
    data : 2D array, background-subtracted
    catalog : list of dicts with source properties
    cfg : MorphometryConfig

    Returns
    -------
    list of dicts with morphometry results
    """
    if cfg is None:
        cfg = MorphometryConfig()

    results = []
    for src in catalog:
        result = measure_source(data, src, cfg)
        result['NUMBER'] = src.get('number', 0)
        results.append(result)

    return results
