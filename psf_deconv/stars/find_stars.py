"""Star finding algorithms for PSF extraction."""

import numpy as np


def find_stars_class_star(catalog, thresh=0.8):
    """Select stars using CLASS_STAR threshold.

    Parameters
    ----------
    catalog : dict
        Catalog with 'class_star' array.
    thresh : float
        Minimum CLASS_STAR value.

    Returns
    -------
    mask : ndarray of bool
        True for sources classified as stars.
    """
    cs = np.asarray(catalog['class_star'], dtype=float)
    return cs >= thresh


def find_stars_fwhm(catalog, sigma=2.0):
    """Select stars using FWHM clustering (histogram peak detection).

    Stars cluster at a characteristic FWHM. We find the histogram peak
    and select sources within sigma * MAD of that peak.

    Parameters
    ----------
    catalog : dict
        Catalog with 'a_image' and 'b_image' arrays.
    sigma : float
        Number of MADs around the peak FWHM to accept.

    Returns
    -------
    mask : ndarray of bool
        True for sources near the FWHM peak (likely stars).
    """
    a = np.asarray(catalog['a_image'], dtype=float)
    b = np.asarray(catalog['b_image'], dtype=float)
    fwhm = np.sqrt(a * b) * 2.3548  # FWHM from geometric mean of axes

    # Histogram-based peak detection
    valid = np.isfinite(fwhm) & (fwhm > 0)
    if valid.sum() < 5:
        return np.zeros(len(a), dtype=bool)

    fwhm_valid = fwhm[valid]
    n_bins = max(10, min(50, int(np.sqrt(len(fwhm_valid)))))
    counts, edges = np.histogram(fwhm_valid, bins=n_bins)
    peak_idx = np.argmax(counts)
    peak_fwhm = 0.5 * (edges[peak_idx] + edges[peak_idx + 1])

    # MAD-based selection around peak
    mad = np.median(np.abs(fwhm_valid - peak_fwhm))
    if mad < 0.1:
        mad = 0.1  # minimum spread

    mask = np.abs(fwhm - peak_fwhm) < sigma * mad
    return mask & valid


def find_stars_combined(catalog, config=None):
    """Combined star selection using CLASS_STAR + FWHM + ellipticity + flux S/N.

    Parameters
    ----------
    catalog : dict
        Catalog with keys: class_star, a_image, b_image, ellipticity,
        flux_auto, fluxerr_auto (or flux_snr).
    config : PSFDeconvConfig, optional
        Configuration object. Uses defaults if None.

    Returns
    -------
    mask : ndarray of bool
        True for sources identified as stars.
    info : dict
        Diagnostic info: n_candidates per criterion, median FWHM.
    """
    if config is None:
        from ..config import PSFDeconvConfig
        config = PSFDeconvConfig()

    n = len(catalog['a_image'])

    # Criterion 1: CLASS_STAR
    cs = np.asarray(catalog['class_star'], dtype=float)
    mask_cs = cs >= config.class_star_thresh

    # Criterion 2: Ellipticity
    ell = np.asarray(catalog['ellipticity'], dtype=float)
    mask_ell = ell < config.max_ellipticity

    # Criterion 3: FWHM clustering
    mask_fwhm = find_stars_fwhm(catalog, sigma=config.fwhm_sigma)

    # Criterion 4: Flux S/N
    if 'flux_snr' in catalog:
        snr = np.asarray(catalog['flux_snr'], dtype=float)
    else:
        flux = np.asarray(catalog['flux_auto'], dtype=float)
        if 'fluxerr_auto' in catalog:
            fluxerr = np.asarray(catalog['fluxerr_auto'], dtype=float)
            fluxerr = np.where(fluxerr > 0, fluxerr, 1.0)
            snr = flux / fluxerr
        else:
            # Fallback: use flux rank as proxy
            snr = flux / max(np.median(flux[flux > 0]), 1.0) * 20.0
    mask_snr = snr >= config.min_flux_snr

    # Criterion 5: Not saturated (FLAGS == 0 or < 4)
    if 'flags' in catalog:
        flags = np.asarray(catalog['flags'], dtype=int)
        mask_flags = flags < 4
    else:
        mask_flags = np.ones(n, dtype=bool)

    # Combined: require at least 3 of 4 soft criteria + flags
    score = (mask_cs.astype(int) + mask_ell.astype(int) +
             mask_fwhm.astype(int) + mask_snr.astype(int))
    mask = (score >= 3) & mask_flags

    # If too few stars, relax to score >= 2
    if mask.sum() < 3:
        mask = (score >= 2) & mask_flags

    # Compute diagnostic FWHM
    a = np.asarray(catalog['a_image'], dtype=float)
    b = np.asarray(catalog['b_image'], dtype=float)
    fwhm = np.sqrt(a * b) * 2.3548
    median_fwhm = float(np.median(fwhm[mask])) if mask.any() else 0.0

    info = {
        'n_class_star': int(mask_cs.sum()),
        'n_ellipticity': int(mask_ell.sum()),
        'n_fwhm': int(mask_fwhm.sum()),
        'n_snr': int(mask_snr.sum()),
        'n_combined': int(mask.sum()),
        'median_fwhm': median_fwhm,
    }

    return mask, info
