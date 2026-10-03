"""Feature extraction for photometric redshift estimation."""

import numpy as np


# Sentinel value for missing data
MISSING_VALUE = -99.0


def _is_missing(val):
    """Check if a value is missing/invalid."""
    if val is None:
        return True
    try:
        v = float(val)
        if np.isnan(v) or np.isinf(v):
            return True
        if v >= 99.0 or v <= -99.0:
            return True
        return False
    except (ValueError, TypeError):
        return True


def extract_features(mag_dict, mag_err_dict, band_names,
                     extra_features=None):
    """
    Extract feature vectors for photometric redshift estimation.

    For each source, features consist of:
    - Per-band: magnitude, magnitude error, presence flag (1=valid, 0=missing)
    - Color indices: adjacent band differences (b1-b2, b2-b3, ...)
    - Extra features: e.g. R_EFF, ELLIPTICITY if available

    Missing values (NaN, >=99, None) are replaced with MISSING_VALUE (-99)
    and the corresponding presence flag is set to 0.

    Parameters
    ----------
    mag_dict : dict
        {band_name: array_of_magnitudes} for each band.
        Each array has shape (n_sources,).
    mag_err_dict : dict
        {band_name: array_of_magnitude_errors} for each band.
        Each array has shape (n_sources,).
    band_names : list of str
        Ordered list of band names (e.g. ['u', 'g', 'r', 'i', 'z']).
        Order determines color index computation.
    extra_features : dict or None
        {feature_name: array} of additional features per source.
        Each array has shape (n_sources,).

    Returns
    -------
    features : numpy.ndarray
        Shape (n_sources, n_features). Feature array ready for model input.
    feature_names : list of str
        Names of each feature column.
    """
    # Determine number of sources
    first_band = band_names[0]
    n_sources = len(mag_dict[first_band])

    feature_list = []
    feature_names = []

    # Per-band features: magnitude, error, flag
    for band in band_names:
        mags = np.array(mag_dict[band], dtype=np.float64)
        errs = np.array(mag_err_dict.get(band, np.full(n_sources, np.nan)),
                        dtype=np.float64)

        clean_mags = np.full(n_sources, MISSING_VALUE, dtype=np.float64)
        clean_errs = np.full(n_sources, MISSING_VALUE, dtype=np.float64)
        flags = np.zeros(n_sources, dtype=np.float64)

        for i in range(n_sources):
            if not _is_missing(mags[i]):
                clean_mags[i] = mags[i]
                flags[i] = 1.0
                if not _is_missing(errs[i]):
                    clean_errs[i] = errs[i]
                # If mag is valid but err is missing, err stays MISSING_VALUE

        feature_list.append(clean_mags)
        feature_names.append(f"MAG_{band}")

        feature_list.append(clean_errs)
        feature_names.append(f"MAGERR_{band}")

        feature_list.append(flags)
        feature_names.append(f"FLAG_{band}")

    # Color indices: adjacent band differences
    for i in range(len(band_names) - 1):
        b1 = band_names[i]
        b2 = band_names[i + 1]

        mags1 = np.array(mag_dict[b1], dtype=np.float64)
        mags2 = np.array(mag_dict[b2], dtype=np.float64)

        colors = np.full(n_sources, MISSING_VALUE, dtype=np.float64)
        for j in range(n_sources):
            if not _is_missing(mags1[j]) and not _is_missing(mags2[j]):
                colors[j] = float(mags1[j]) - float(mags2[j])

        feature_list.append(colors)
        feature_names.append(f"COLOR_{b1}-{b2}")

    # Extra features
    if extra_features is not None:
        for name, values in extra_features.items():
            vals = np.array(values, dtype=np.float64)
            clean_vals = np.full(n_sources, MISSING_VALUE, dtype=np.float64)
            for i in range(n_sources):
                if not _is_missing(vals[i]):
                    clean_vals[i] = vals[i]
            feature_list.append(clean_vals)
            feature_names.append(name)

    features = np.column_stack(feature_list)
    return features, feature_names


def extract_features_block(mags, mag_errs):
    """Feature matrix in the column order of the training table (photo_z/data/sdss_specphoto.h5):
    [mag_1..n, magerr_1..n, flag_1..n, colour_12 .. colour_(n-1)n]  (block layout; `extract_features` is interleaved per band).
    Further training columns (e.g. petroR50_r) are padded by `predict_photoz` with the training mean.

    mags, mag_errs : arrays (n_sources, n_bands); missing = NaN or >= 99."""
    mags = np.asarray(mags, dtype=np.float64)
    errs = np.asarray(mag_errs, dtype=np.float64)
    bad = ~np.isfinite(mags) | (mags >= 99.0) | (mags <= -99.0)
    badr = bad | ~np.isfinite(errs) | (errs >= 99.0) | (errs <= -99.0)
    m = np.where(bad, MISSING_VALUE, mags)
    e = np.where(badr, MISSING_VALUE, errs)
    flags = (~bad).astype(np.float64)
    cols = [m, e, flags]
    if mags.shape[1] > 1:
        col = mags[:, :-1] - mags[:, 1:]
        cbad = bad[:, :-1] | bad[:, 1:]
        cols.append(np.where(cbad, MISSING_VALUE, col))
    return np.hstack(cols)
