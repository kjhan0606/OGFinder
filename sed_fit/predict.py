"""Inference for SED fitting: observed photometry -> physical parameters."""

import os
import sys
import numpy as np
import torch

from .config import SEDConfig, TrainConfig, OUTPUT_PARAM_NAMES
from .models.inverse import InverseSEDMLP
from .features import normalize_magnitudes


def predict_sed(features, photo_z, checkpoint_inverse, n_bands,
                config=None):
    """
    Run inverse SED model to estimate physical parameters from photometry.

    Parameters
    ----------
    features : np.ndarray
        Shape (N, n_bands). Observed magnitudes (raw, will be normalized).
    photo_z : np.ndarray
        Shape (N,). Photometric redshifts.
    checkpoint_inverse : str
        Path to inverse model checkpoint (.pt).
    n_bands : int
        Number of photometric bands.
    config : SEDConfig, optional
        Model configuration. Uses defaults if None.

    Returns
    -------
    results : dict
        Keys: LOG_MASS, LOG_MASS_ERR, LOG_AGE, LOG_AGE_ERR, LOG_Z,
        LOG_Z_ERR, AV, AV_ERR, LOG_TAU, LOG_TAU_ERR, SFR, SFR_ERR.
        Values: np.ndarray of shape (N,).
    """
    if config is None:
        config = SEDConfig()

    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

    # Load checkpoint
    if not os.path.exists(checkpoint_inverse):
        print(f"ERROR: checkpoint not found: {checkpoint_inverse}",
              file=sys.stderr)
        return None

    ckpt = torch.load(checkpoint_inverse, map_location=device,
                      weights_only=True)

    # Reconstruct model from checkpoint config
    ckpt_cfg = ckpt.get('config', {})
    n_input = ckpt_cfg.get('n_input', n_bands + 1)
    n_output = ckpt_cfg.get('n_output', config.n_output_params)
    n_hidden = ckpt_cfg.get('n_hidden', config.n_hidden_inverse)
    n_components = ckpt_cfg.get('n_components', config.n_components)
    dropout = ckpt_cfg.get('dropout', config.dropout)

    model = InverseSEDMLP(
        n_input=n_input,
        n_output=n_output,
        n_hidden=n_hidden,
        n_components=n_components,
        dropout=dropout,
    ).to(device)

    model.load_state_dict(ckpt['model_state'])
    model.eval()

    # Prepare input
    features = np.asarray(features, dtype=np.float32)
    photo_z = np.atleast_1d(np.asarray(photo_z, dtype=np.float32))
    N = features.shape[0]

    # Normalize magnitudes
    mags_norm = normalize_magnitudes(features)

    # Apply training normalization if stored in checkpoint
    if 'norm_mean' in ckpt and 'norm_std' in ckpt:
        norm_mean = np.asarray(ckpt['norm_mean'], dtype=np.float32)
        norm_std = np.asarray(ckpt['norm_std'], dtype=np.float32)
        # norm_mean/std cover [n_bands + 1] features
        if len(norm_mean) == n_bands + 1:
            full_input = np.column_stack([mags_norm, photo_z])
            full_input = (full_input - norm_mean) / np.maximum(norm_std, 1e-8)
            mags_norm = full_input[:, :n_bands]
            photo_z = full_input[:, n_bands]

    # Concatenate: [normalized_mags, photo_z]
    input_data = np.column_stack([mags_norm, photo_z])
    input_tensor = torch.from_numpy(input_data).to(device)

    # Run inference in batches
    batch_size = config.batch_size
    all_results = {name: {'value': [], 'error': []} for name in OUTPUT_PARAM_NAMES}

    with torch.no_grad():
        for i in range(0, N, batch_size):
            batch = input_tensor[i:i + batch_size]
            pred = model.predict(batch)
            for name in OUTPUT_PARAM_NAMES:
                all_results[name]['value'].append(pred[name]['value'])
                all_results[name]['error'].append(pred[name]['error'])

    # Concatenate batches
    results = {}
    for name in OUTPUT_PARAM_NAMES:
        val = np.concatenate(all_results[name]['value'])
        err = np.concatenate(all_results[name]['error'])
        key = name.upper()
        results[key] = val
        results[key + '_ERR'] = err

    # Estimate SFR from mass and age: SFR ~ M* / age
    # log_mass is log10(M*/Msun), log_age is log10(age/yr)
    log_mass = results['LOG_MASS']
    log_age = results['LOG_AGE']
    log_mass_err = results['LOG_MASS_ERR']
    log_age_err = results['LOG_AGE_ERR']

    # SFR = 10^log_mass / 10^log_age  [Msun/yr]
    # log10(SFR) = log_mass - log_age
    log_sfr = log_mass - log_age
    results['SFR'] = 10.0 ** log_sfr

    # Error propagation (in log space, then exponentiated)
    log_sfr_err = np.sqrt(log_mass_err ** 2 + log_age_err ** 2)
    # Fractional error: sigma_SFR / SFR = ln(10) * sigma_log_SFR
    results['SFR_ERR'] = results['SFR'] * np.log(10.0) * log_sfr_err

    return results
