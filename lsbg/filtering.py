"""LSBG candidate filtering and confidence scoring.

Applies surface brightness, size, shape, SNR, and Sérsic profile
criteria to select genuine LSBG candidates from the detection catalog.
Assigns A/B/C/D confidence grades.
"""

import numpy as np


def filter_lsbg_candidates(catalog, config, sersic_results=None):
    """Filter detected sources by LSBG selection criteria.

    Parameters
    ----------
    catalog : list of dict
        Photometry results from measure_photometry.
    config : LSBGConfig
        Pipeline configuration with filtering thresholds.
    sersic_results : list of dict or None
        Sérsic fit results (parallel to catalog). None entries = fit failed.

    Returns
    -------
    filtered : list of dict
        Sources passing all LSBG criteria (with LSBG_CONF/GRADE added).
    rejected : list of dict
        Sources failing criteria.
    """
    filtered = []
    rejected = []

    for i, src in enumerate(catalog):
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

        # Sérsic-based filtering (if available)
        sersic = None
        if sersic_results is not None and i < len(sersic_results):
            sersic = sersic_results[i]

        if sersic is not None and passes:
            n = sersic.get('n', np.nan)
            chi2 = sersic.get('chi2', np.nan)

            # Sérsic index cut: reject artifacts (n too low) or
            # compact sources (n too high for LSBG)
            if np.isfinite(n):
                if n < config.sersic_n_filter_min:
                    passes = False
                elif n > config.sersic_n_filter_max:
                    passes = False

            # Fit quality cut
            if np.isfinite(chi2) and chi2 > config.sersic_chi2_max:
                passes = False

        # Add Sérsic results to source dict
        if sersic is not None:
            src['sersic_n'] = sersic.get('n', np.nan)
            src['sersic_re'] = sersic.get('re', np.nan)
            src['sersic_chi2'] = sersic.get('chi2', np.nan)
            # Sérsic total magnitude
            flux_total = sersic.get('flux_total', np.nan)
            if np.isfinite(flux_total) and flux_total > 0:
                zp = config.mag_zeropoint
                src['mag_sersic'] = -2.5 * np.log10(flux_total) + zp
                # Sérsic effective surface brightness
                re_arcsec = sersic.get('re', 0) * config.pixel_scale
                ellip_s = sersic.get('ellip', 0)
                area = np.pi * re_arcsec**2 * (1.0 - ellip_s)
                if area > 0:
                    src['mu_eff_sersic'] = src['mag_sersic'] + 2.5 * np.log10(2.0 * area)
                else:
                    src['mu_eff_sersic'] = np.nan
            else:
                src['mag_sersic'] = np.nan
                src['mu_eff_sersic'] = np.nan
        else:
            src['sersic_n'] = np.nan
            src['sersic_re'] = np.nan
            src['sersic_chi2'] = np.nan
            src['mag_sersic'] = np.nan
            src['mu_eff_sersic'] = np.nan

        if passes:
            src['lsbg_conf'] = classify_confidence(src, config)
            src['lsbg_grade'] = assign_grade(src['lsbg_conf'])
            filtered.append(src)
        else:
            src['lsbg_conf'] = 0.0
            src['lsbg_grade'] = 'D'
            rejected.append(src)

    return filtered, rejected


def classify_confidence(source, config):
    """Compute LSBG confidence score from source properties.

    Combines surface brightness, size, shape, SNR, and Sérsic profile
    quality into a 0-1 confidence score.

    Parameters
    ----------
    source : dict
        Source photometry + Sérsic measurements.
    config : LSBGConfig
        Pipeline configuration.

    Returns
    -------
    confidence : float
        Score between 0 and 1.
    """
    scores = []
    weights = []

    # Surface brightness score: peaks at center of mu_eff range
    mu = source.get('mu_eff', np.nan)
    if np.isfinite(mu):
        mu_center = (config.mu_eff_min + config.mu_eff_max) / 2.0
        mu_range = (config.mu_eff_max - config.mu_eff_min) / 2.0
        mu_score = 1.0 - abs(mu - mu_center) / mu_range
        scores.append(max(0.0, min(1.0, mu_score)))
    else:
        scores.append(0.0)
    weights.append(0.25)

    # Size score: larger is better (within range)
    r = source.get('r_eff_arcsec', np.nan)
    if np.isfinite(r) and config.r_eff_max > config.r_eff_min:
        r_score = (r - config.r_eff_min) / (config.r_eff_max - config.r_eff_min)
        scores.append(max(0.0, min(1.0, r_score)))
    else:
        scores.append(0.0)
    weights.append(0.20)

    # Roundness score: rounder is better
    ellip = source.get('ellipticity', 1.0)
    e_score = 1.0 - ellip / config.ellipticity_max
    scores.append(max(0.0, min(1.0, e_score)))
    weights.append(0.15)

    # SNR score: higher is better (saturates at 10x threshold)
    snr = source.get('snr_total', 0.0)
    if config.min_snr > 0:
        snr_score = min(1.0, snr / (10.0 * config.min_snr))
    else:
        snr_score = min(1.0, snr / 20.0)
    scores.append(max(0.0, snr_score))
    weights.append(0.15)

    # Sérsic profile quality score
    sersic_n = source.get('sersic_n', np.nan)
    sersic_chi2 = source.get('sersic_chi2', np.nan)
    if np.isfinite(sersic_n) and np.isfinite(sersic_chi2):
        # Sérsic n score: favor typical LSBG range (0.5-2.0)
        if 0.5 <= sersic_n <= 2.0:
            n_score = 1.0
        elif sersic_n < 0.5:
            n_score = max(0.0, sersic_n / 0.5)
        else:
            n_score = max(0.0, 1.0 - (sersic_n - 2.0) / 4.0)
        # Chi2 score: lower is better (saturates at 1.0)
        chi2_score = max(0.0, 1.0 - sersic_chi2 / config.sersic_chi2_max)
        sersic_score = 0.6 * n_score + 0.4 * chi2_score
        scores.append(max(0.0, min(1.0, sersic_score)))
        weights.append(0.25)
    else:
        # No Sérsic fit — redistribute weight
        scores.append(0.5)  # neutral score
        weights.append(0.10)

    # Normalize weights
    total_w = sum(weights)
    confidence = sum(s * w for s, w in zip(scores, weights)) / total_w

    return round(confidence, 3)


def assign_grade(confidence):
    """Assign A/B/C/D letter grade based on confidence score.

    A (>=0.7): High-confidence LSBG candidate
    B (>=0.5): Probable LSBG candidate
    C (>=0.3): Marginal candidate, needs visual inspection
    D (<0.3): Likely artifact or non-LSBG
    """
    if confidence >= 0.7:
        return 'A'
    elif confidence >= 0.5:
        return 'B'
    elif confidence >= 0.3:
        return 'C'
    else:
        return 'D'
