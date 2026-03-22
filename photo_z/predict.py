"""Inference for photometric redshift estimation."""

import numpy as np
import torch

from .config import PhotoZConfig
from .models.mlp_mdn import PhotoZMDN


def predict_photoz(features, checkpoint_path, config=None):
    """
    Run photometric redshift inference on extracted features.

    Parameters
    ----------
    features : numpy.ndarray
        Shape (n_sources, n_features). Feature array from extract_features().
    checkpoint_path : str
        Path to model checkpoint (.pt file).
    config : PhotoZConfig or None
        Model configuration. If None, loaded from checkpoint metadata
        or defaults are used.

    Returns
    -------
    results : dict
        Dictionary with keys:
        - 'z_point': numpy.ndarray, shape (n_sources,), point estimates
        - 'z_err': numpy.ndarray, shape (n_sources,), uncertainty (weighted std)
        - 'z_q68': numpy.ndarray, shape (n_sources,), 68% credible interval width
        - 'p_reliable': numpy.ndarray or None, shape (n_sources,),
          P(reliable) from Quality Head v2. None if quality_head not enabled.
    """
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

    # Load checkpoint
    ckpt = torch.load(checkpoint_path, map_location=device, weights_only=False)

    # Determine model config from checkpoint or argument
    if config is None:
        config = PhotoZConfig()
    ckpt_config = ckpt.get('config', {})
    n_features = ckpt_config.get('n_features', features.shape[1])
    n_hidden = ckpt_config.get('n_hidden', config.n_hidden)
    n_components = ckpt_config.get('n_components', config.n_components)
    dropout = ckpt_config.get('dropout', config.dropout)
    quality_head = ckpt_config.get('quality_head', False)

    # Build model and load weights
    model = PhotoZMDN(
        n_features=n_features,
        n_hidden=n_hidden,
        n_components=n_components,
        dropout=dropout,
        quality_head=quality_head,
    ).to(device)
    model.load_state_dict(ckpt['model_state'])
    model.eval()

    # Load normalization parameters if saved
    feat_mean = ckpt.get('feat_mean', None)
    feat_std = ckpt.get('feat_std', None)

    # Prepare features
    features = np.array(features, dtype=np.float32)
    if feat_mean is not None and feat_std is not None:
        feat_mean = np.array(feat_mean, dtype=np.float32)
        feat_std = np.array(feat_std, dtype=np.float32)
        feat_std[feat_std < 1e-8] = 1.0
        features = (features - feat_mean) / feat_std

    n_sources = features.shape[0]
    batch_size = config.batch_size

    z_point_all = []
    z_err_all = []
    z_q68_all = []
    p_reliable_all = []

    for start in range(0, n_sources, batch_size):
        end = min(start + batch_size, n_sources)
        batch = torch.from_numpy(features[start:end]).to(device)

        z_point, z_err, z_q68, p_reliable = model.predict(batch)

        z_point_all.append(z_point.cpu().numpy())
        z_err_all.append(z_err.cpu().numpy())
        z_q68_all.append(z_q68.cpu().numpy())
        if p_reliable is not None:
            p_reliable_all.append(p_reliable.cpu().numpy())

    results = {
        'z_point': np.concatenate(z_point_all),
        'z_err': np.concatenate(z_err_all),
        'z_q68': np.concatenate(z_q68_all),
        'p_reliable': np.concatenate(p_reliable_all) if p_reliable_all else None,
    }

    return results
