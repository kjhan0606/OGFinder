"""Model loading and inference for star/galaxy classification."""

import os
import numpy as np
import torch

from ..config import StarFinderConfig
from ..models.cnn import StarFinderCNN


def load_star_model(checkpoint_path, device='cpu'):
    """Load trained star/galaxy model from checkpoint.

    Args:
        checkpoint_path: path to .pt checkpoint file
        device: 'cpu' or 'cuda'

    Returns:
        model: StarFinderCNN in eval mode
    """
    ckpt = torch.load(checkpoint_path, map_location=device, weights_only=False)
    cfg = ckpt['config']

    model = StarFinderCNN(in_channels=cfg.get('in_channels', 1))
    model.load_state_dict(ckpt['model_state'])
    model.eval()
    model.to(device)

    return model


def predict_batch(model, cutouts, device='cpu', batch_size=256):
    """Run batch inference on cutouts.

    Args:
        model: StarFinderCNN in eval mode
        cutouts: (N, 1, H, W) float32 numpy array
        device: torch device
        batch_size: inference batch size

    Returns:
        proba: (N,) float32 numpy array of P(star)
    """
    tensor = torch.FloatTensor(cutouts)

    model.eval()
    all_proba = []

    with torch.no_grad():
        for start in range(0, len(tensor), batch_size):
            batch = tensor[start:start + batch_size].to(device)
            outputs = model(batch)  # (B, 1)
            all_proba.append(outputs.cpu().numpy().flatten())

    return np.concatenate(all_proba, axis=0)


def classify_sources(data, catalog, checkpoint_path=None,
                     config=None, device='cpu'):
    """Full pipeline: cutout extraction -> prediction -> results.

    Args:
        data: 2D FITS image array (background-subtracted)
        catalog: dict with 'number', 'x', 'y' arrays
        checkpoint_path: path to model checkpoint
        config: StarFinderConfig
        device: 'cpu' or 'cuda'

    Returns:
        results: list of (number, ai_star_prob, ai_star_label)
    """
    from ..cutout.extract import batch_extract_cutouts

    if config is None:
        config = StarFinderConfig()

    if checkpoint_path is None:
        checkpoint_path = os.path.join(
            os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
            "star_finder", "data", "checkpoints", "star_finder_best.pt"
        )

    # Extract cutouts
    cutouts, valid_mask = batch_extract_cutouts(
        data, catalog, size=config.cutout_size)

    if cutouts.shape[0] == 0:
        return []

    # Load model and predict
    model = load_star_model(checkpoint_path, device)
    proba = predict_batch(model, cutouts, device)

    # Build results
    results = []
    valid_count = 0
    for i in range(len(catalog['number'])):
        if not valid_mask[i]:
            continue
        p = float(proba[valid_count])
        valid_count += 1
        label = "star" if p >= config.star_threshold else "galaxy"
        results.append((int(catalog['number'][i]), p, label))

    return results
