"""SVM-based LSBG candidate classification.

Extracts photometric + structural features from LSBG candidates and
classifies them using a trained SVM pipeline (SimpleImputer -> StandardScaler -> SVC).
"""

import os
import numpy as np


FEATURE_NAMES = [
    'MU_EFF', 'R_EFF_ARCSEC', 'ELLIPTICITY', 'MAG_AUTO',
    'SERSIC_N', 'SERSIC_CHI2', 'SNR_TOTAL', 'KRON_RADIUS',
    'MU_EFF_SERSIC', 'CONC_5_20', 'CONC_10_40',
]

# Default checkpoint path (relative to this file's package)
_DEFAULT_CHECKPOINT = os.path.join(
    os.path.dirname(__file__), 'data', 'checkpoints', 'svm_lsbg_best.joblib')


def _safe_float(val):
    """Convert value to float, returning NaN on failure."""
    try:
        v = float(val)
        return v if np.isfinite(v) else np.nan
    except (ValueError, TypeError):
        return np.nan


def _safe_ratio(a, b):
    """Compute a/b, returning NaN if b is zero or either is NaN."""
    a = _safe_float(a)
    b = _safe_float(b)
    if not np.isfinite(a) or not np.isfinite(b) or b == 0:
        return np.nan
    return a / b


def extract_features(catalog_rows):
    """Extract feature matrix from catalog rows.

    Parameters
    ----------
    catalog_rows : list of dict
        Each dict has column names as keys (case-insensitive lookup).

    Returns
    -------
    features : np.ndarray, shape (N, 11)
        Feature matrix. Missing values are NaN (handled by pipeline imputer).
    """
    n = len(catalog_rows)
    X = np.full((n, len(FEATURE_NAMES)), np.nan)

    for i, row in enumerate(catalog_rows):
        # Case-insensitive key lookup
        r = {k.upper(): v for k, v in row.items()}

        X[i, 0] = _safe_float(r.get('MU_EFF'))
        X[i, 1] = _safe_float(r.get('R_EFF_ARCSEC'))
        X[i, 2] = _safe_float(r.get('ELLIPTICITY'))
        X[i, 3] = _safe_float(r.get('MAG_AUTO'))
        X[i, 4] = _safe_float(r.get('SERSIC_N'))

        # log1p transforms
        chi2 = _safe_float(r.get('SERSIC_CHI2'))
        X[i, 5] = np.log1p(chi2) if np.isfinite(chi2) else np.nan

        snr = _safe_float(r.get('SNR_TOTAL'))
        X[i, 6] = np.log1p(snr) if np.isfinite(snr) else np.nan

        X[i, 7] = _safe_float(r.get('KRON_RADIUS'))
        X[i, 8] = _safe_float(r.get('MU_EFF_SERSIC'))

        # Concentration indices
        X[i, 9] = _safe_ratio(r.get('FLUX_APER_5'), r.get('FLUX_APER_20'))
        X[i, 10] = _safe_ratio(r.get('FLUX_APER_10'), r.get('FLUX_APER_40'))

    return X


def load_model(checkpoint_path=None):
    """Load trained SVM pipeline from joblib checkpoint.

    Parameters
    ----------
    checkpoint_path : str or None
        Path to .joblib file. None uses default path.

    Returns
    -------
    pipeline : sklearn.pipeline.Pipeline
        Fitted pipeline (SimpleImputer -> StandardScaler -> SVC).

    Raises
    ------
    FileNotFoundError
        If checkpoint file does not exist.
    """
    import joblib

    path = checkpoint_path or _DEFAULT_CHECKPOINT
    if not os.path.exists(path):
        raise FileNotFoundError(f"SVM checkpoint not found: {path}")
    return joblib.load(path)


def classify(catalog_rows, checkpoint_path=None, threshold=0.5):
    """Classify LSBG candidates using trained SVM.

    Parameters
    ----------
    catalog_rows : list of dict
        Catalog entries with photometric columns.
    checkpoint_path : str or None
        Path to model checkpoint.
    threshold : float
        Probability threshold for LSBG classification (default 0.5).

    Returns
    -------
    results : list of dict
        Each dict: {number, svm_lsbg (0/1), svm_conf (probability)}.
    """
    if len(catalog_rows) == 0:
        return []

    model = load_model(checkpoint_path)
    X = extract_features(catalog_rows)
    proba = model.predict_proba(X)[:, 1]  # probability of class 1 (LSBG)

    results = []
    for i, row in enumerate(catalog_rows):
        r = {k.upper(): v for k, v in row.items()}
        number = r.get('NUMBER', i + 1)
        conf = float(proba[i])
        results.append({
            'number': int(float(number)),
            'svm_lsbg': 1 if conf >= threshold else 0,
            'svm_conf': round(conf, 4),
        })

    return results
