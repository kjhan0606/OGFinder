"""Model loading and inference for galaxy morphology classification."""

import os
import numpy as np
import torch

from ..config import GALAXY10_CLASSES, MorphConfig
from ..models.cnn import GalaxyMorphCNN


def load_morph_model(checkpoint_path, device='cpu'):
    """Load trained morphology model from checkpoint.

    Args:
        checkpoint_path: path to .pt checkpoint file
        device: 'cpu' or 'cuda'

    Returns:
        model: GalaxyMorphCNN in eval mode
        norm_stats: dict with 'mean' and 'std' (numpy arrays, per-channel)
        class_names: list of class name strings
    """
    ckpt = torch.load(checkpoint_path, map_location=device, weights_only=False)
    cfg = ckpt['config']

    in_channels = cfg.get('in_channels', 3)
    model = GalaxyMorphCNN(
        n_classes=cfg['n_classes'],
        in_channels=in_channels,
        dropout=cfg['dropout'],
        dropout2=cfg['dropout2'],
        hidden_dim=cfg['hidden_dim'],
    )
    model.load_state_dict(ckpt['model_state'])
    model.eval()
    model.to(device)

    mean = ckpt['norm_mean']
    std = ckpt['norm_std']
    # Convert tensors to numpy if needed
    if hasattr(mean, 'cpu'):
        mean = mean.cpu().numpy()
    if hasattr(std, 'cpu'):
        std = std.cpu().numpy()
    mean = np.asarray(mean, dtype=np.float32)
    std = np.asarray(std, dtype=np.float32)
    norm_stats = {'mean': mean, 'std': std, 'in_channels': in_channels}

    class_names = ckpt.get('class_names', GALAXY10_CLASSES)

    return model, norm_stats, class_names


def predict_batch(model, cutouts, norm_stats, device='cpu', batch_size=256):
    """Run batch inference on cutouts.

    Args:
        model: GalaxyMorphCNN in eval mode
        cutouts: (N, 1, H, W) float32 numpy array (single-channel FITS cutouts)
        norm_stats: dict with 'mean', 'std', 'in_channels'
        device: torch device
        batch_size: inference batch size

    Returns:
        proba: (N, n_classes) float32 numpy array of softmax probabilities
    """
    in_channels = norm_stats.get('in_channels', 3)
    mean = norm_stats['mean']
    std = norm_stats['std']

    # For FITS single-channel input with RGB-trained model:
    # duplicate channel 3x, then apply per-channel normalization
    if cutouts.shape[1] == 1 and in_channels == 3:
        cutouts = np.repeat(cutouts, 3, axis=1)  # (N, 3, H, W)

    # Per-channel standardization
    cutouts_norm = cutouts.copy()
    if mean.ndim == 0:
        # Scalar mean/std (old checkpoint)
        cutouts_norm = (cutouts_norm - float(mean)) / (float(std) + 1e-8)
    else:
        # Per-channel mean/std
        for c in range(cutouts_norm.shape[1]):
            c_idx = min(c, len(mean) - 1)
            cutouts_norm[:, c] = (cutouts_norm[:, c] - mean[c_idx]) / (std[c_idx] + 1e-8)

    tensor = torch.FloatTensor(cutouts_norm)

    model.eval()
    all_proba = []

    with torch.no_grad():
        for start in range(0, len(tensor), batch_size):
            batch = tensor[start:start + batch_size].to(device)
            logits = model(batch)
            proba = torch.softmax(logits, dim=1)
            all_proba.append(proba.cpu().numpy())

    return np.concatenate(all_proba, axis=0)


def classify_sources(cutouts, checkpoint_path, device='cpu', batch_size=256):
    """End-to-end classification: load model + predict.

    Args:
        cutouts: (N, 1, H, W) float32 numpy array (already asinh-normalized)
        checkpoint_path: path to model checkpoint
        device: 'cpu' or 'cuda'
        batch_size: inference batch size

    Returns:
        results: list of dicts, one per source, with keys:
            'class_id', 'class_name', 'confidence',
            'top3_classes', 'top3_probs', 'color'
    """
    from ..config import MORPH_COLORS

    model, norm_stats, class_names = load_morph_model(checkpoint_path, device)
    proba = predict_batch(model, cutouts, norm_stats, device, batch_size)

    results = []
    for i in range(len(proba)):
        p = proba[i]
        top3_idx = np.argsort(p)[::-1][:3]
        top1_id = top3_idx[0]

        results.append({
            'class_id': int(top1_id),
            'class_name': class_names[top1_id],
            'confidence': float(p[top1_id]),
            'top3_classes': [class_names[j] for j in top3_idx],
            'top3_probs': [float(p[j]) for j in top3_idx],
            'color': MORPH_COLORS[top1_id],
            'proba': p.tolist(),
        })

    return results
