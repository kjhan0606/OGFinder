"""Feature processing for SED fitting models."""

import numpy as np
import torch

from .config import PHYS_PARAM_NAMES


def prepare_emulator_input(phys_params_dict):
    """
    Convert a dict of physical parameters to a tensor for the emulator.

    Parameters
    ----------
    phys_params_dict : dict
        Keys are parameter names (z, log_mass, log_age, log_Z, Av, log_tau).
        Values are scalars or 1-D arrays.

    Returns
    -------
    tensor : torch.Tensor
        Shape (N, 6) or (1, 6) for scalar inputs.
    """
    arrays = []
    for name in PHYS_PARAM_NAMES:
        val = phys_params_dict[name]
        arr = np.atleast_1d(np.asarray(val, dtype=np.float32))
        arrays.append(arr)

    # Stack into (N, 6)
    data = np.column_stack(arrays)
    return torch.from_numpy(data)


def prepare_inverse_input(mag_dict, photo_z, band_names):
    """
    Convert observed magnitudes + photo-z to a tensor for the inverse model.

    Parameters
    ----------
    mag_dict : dict
        Keys are band names, values are magnitudes (scalars or 1-D arrays).
        Missing bands should be NaN.
    photo_z : float or array
        Photometric redshift(s).
    band_names : list of str
        Ordered list of band names matching the model's training order.

    Returns
    -------
    tensor : torch.Tensor
        Shape (N, n_bands + 1). Normalized magnitudes followed by photo_z.
    """
    n = None
    mag_arrays = []
    for band in band_names:
        val = mag_dict.get(band, np.nan)
        arr = np.atleast_1d(np.asarray(val, dtype=np.float32))
        if n is None:
            n = len(arr)
        mag_arrays.append(arr)

    mags = np.column_stack(mag_arrays)  # (N, n_bands)

    # Normalize magnitudes (subtract per-source median, handle NaN)
    mags_norm = normalize_magnitudes(mags)

    # Append photo_z
    z_arr = np.atleast_1d(np.asarray(photo_z, dtype=np.float32))
    if len(z_arr) == 1 and n > 1:
        z_arr = np.full(n, z_arr[0], dtype=np.float32)

    features = np.column_stack([mags_norm, z_arr])  # (N, n_bands + 1)
    return torch.from_numpy(features)


def normalize_magnitudes(mags, mag_errors=None):
    """
    Normalize magnitudes for network input.

    Subtracts the per-source median magnitude (color information preserved)
    and replaces NaN values with 0 (neutral after median subtraction).
    Optionally uses mag_errors to down-weight uncertain measurements.

    Parameters
    ----------
    mags : np.ndarray
        Shape (N, n_bands). Raw magnitudes, may contain NaN.
    mag_errors : np.ndarray, optional
        Shape (N, n_bands). Magnitude errors. If provided, NaN bands
        are set to the median and their effective weight is zero.

    Returns
    -------
    normalized : np.ndarray
        Shape (N, n_bands). Median-subtracted, NaN-replaced magnitudes.
    """
    mags = np.asarray(mags, dtype=np.float32).copy()

    # Per-source median (ignoring NaN)
    median = np.nanmedian(mags, axis=1, keepdims=True)

    # Handle all-NaN rows
    all_nan = np.all(np.isnan(mags), axis=1)
    median[all_nan] = 0.0

    # Subtract median (preserves colors)
    normalized = mags - median

    # Replace remaining NaN with 0
    nan_mask = np.isnan(normalized)
    normalized[nan_mask] = 0.0

    return normalized
