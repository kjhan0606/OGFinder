"""LSBG candidate filtering and confidence scoring.

Applies surface brightness, size, shape, and SNR criteria to
select genuine LSBG candidates from the detection catalog.
"""

import numpy as np


def filter_lsbg_candidates(catalog, config):
    """Filter detected sources by LSBG selection criteria.

    Parameters
    ----------
    catalog : list of dict
        Photometry results from measure_photometry.
    config : LSBGConfig
        Pipeline configuration with filtering thresholds.

    Returns
    -------
    filtered : list of dict
        Sources passing all LSBG criteria (with LSBG_CONF added).
    rejected : list of dict
        Sources failing criteria.
    """
    filtered = []
    rejected = []

    for src in catalog:
        mu_eff = src.get('mu_eff', np.nan)
        r_eff = src.get('r_eff_arcsec', np.nan)
        ellip = src.get('ellipticity', 1.0)
        snr = src.get('snr_total', 0.0)

        passes = True

        # Surface brightness cut
        if not np.isfinite(mu_eff):
            passes = False
        elif mu_eff < config.mu_eff_min or mu_eff > config.mu_eff_max:
            passes = False

        # Size cut
        if not np.isfinite(r_eff):
            passes = False
        elif r_eff < config.r_eff_min or r_eff > config.r_eff_max:
            passes = False

        # Ellipticity cut
        if ellip > config.ellipticity_max:
            passes = False

        # SNR cut
        if snr < config.min_snr:
            passes = False

        if passes:
            src['lsbg_conf'] = classify_confidence(src, config)
            filtered.append(src)
        else:
            src['lsbg_conf'] = 0.0
            rejected.append(src)

    return filtered, rejected


def classify_confidence(source, config):
    """Compute LSBG confidence score from source properties.

    Combines surface brightness, size, shape, and SNR into a
    0-1 confidence score. Higher = more likely a genuine LSBG.

    Parameters
    ----------
    source : dict
        Source photometry measurements.
    config : LSBGConfig
        Pipeline configuration.

    Returns
    -------
    confidence : float
        Score between 0 and 1.
    """
    scores = []

    # Surface brightness score: peaks at center of mu_eff range
    mu = source.get('mu_eff', np.nan)
    if np.isfinite(mu):
        mu_center = (config.mu_eff_min + config.mu_eff_max) / 2.0
        mu_range = (config.mu_eff_max - config.mu_eff_min) / 2.0
        mu_score = 1.0 - abs(mu - mu_center) / mu_range
        scores.append(max(0.0, min(1.0, mu_score)))
    else:
        scores.append(0.0)

    # Size score: larger is better (within range)
    r = source.get('r_eff_arcsec', np.nan)
    if np.isfinite(r) and config.r_eff_max > config.r_eff_min:
        r_score = (r - config.r_eff_min) / (config.r_eff_max - config.r_eff_min)
        scores.append(max(0.0, min(1.0, r_score)))
    else:
        scores.append(0.0)

    # Roundness score: rounder is better
    ellip = source.get('ellipticity', 1.0)
    e_score = 1.0 - ellip / config.ellipticity_max
    scores.append(max(0.0, min(1.0, e_score)))

    # SNR score: higher is better (saturates at 10x threshold)
    snr = source.get('snr_total', 0.0)
    if config.min_snr > 0:
        snr_score = min(1.0, snr / (10.0 * config.min_snr))
    else:
        snr_score = min(1.0, snr / 20.0)
    scores.append(max(0.0, snr_score))

    # Weighted combination
    weights = [0.35, 0.25, 0.20, 0.20]
    confidence = sum(s * w for s, w in zip(scores, weights))

    return round(confidence, 3)
