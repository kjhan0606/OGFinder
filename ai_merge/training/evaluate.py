"""Evaluate trained MLP on real FITS data (e.g., HUDF)."""

import os
import sys
import argparse
import numpy as np
import torch
from astropy.io import fits

from ..config import CandidateConfig, TrainConfig
from ..models.mlp import MergeMLP
from ..utils.sep_runner import run_sep_on_fits
from ..catalog.find_candidates import find_candidate_pairs, build_merge_groups
from ..catalog.extract_features import extract_pair_features, FEATURE_NAMES


def load_model(checkpoint_path, device='cpu'):
    """Load trained model and normalization stats from checkpoint."""
    ckpt = torch.load(checkpoint_path, map_location=device, weights_only=True)
    cfg = ckpt['config']
    model = MergeMLP(
        n_features=cfg['n_features'],
        n_classes=cfg['n_classes'],
        hidden_dims=cfg['hidden_dims'],
        dropout_rates=cfg['dropout_rates'],
    )
    model.load_state_dict(ckpt['model_state'])
    model.eval()
    model.to(device)

    mean = ckpt['norm_mean']
    std = ckpt['norm_std']
    # Handle both tensor and numpy formats
    if hasattr(mean, 'cpu'):
        mean = mean.cpu().numpy()
    if hasattr(std, 'cpu'):
        std = std.cpu().numpy()
    norm_stats = {'mean': mean, 'std': std}
    return model, norm_stats


def evaluate_fits(fits_path, checkpoint_path=None, output_dir=None,
                  candidate_cfg=None, verbose=True):
    """
    Run merge/separate evaluation on a real FITS image.

    Parameters
    ----------
    fits_path : str
        Path to FITS image.
    checkpoint_path : str, optional
        Path to model checkpoint.
    output_dir : str, optional
        Directory for output files.
    candidate_cfg : CandidateConfig, optional
    verbose : bool

    Returns
    -------
    results : list of dict
        Per-pair predictions.
    """
    if checkpoint_path is None:
        checkpoint_path = TrainConfig().checkpoint_path
    if output_dir is None:
        output_dir = "SAOImageDS9/ai_merge/data/evaluation"
    if candidate_cfg is None:
        candidate_cfg = CandidateConfig()

    os.makedirs(output_dir, exist_ok=True)

    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')

    # Load model
    if verbose:
        print(f"Loading model from {checkpoint_path}...")
    model, norm_stats = load_model(checkpoint_path, device)

    # Load FITS
    if verbose:
        print(f"Loading FITS: {fits_path}...")
    hdulist = fits.open(fits_path)
    data = None
    for hdu in hdulist:
        if hdu.data is not None and hdu.data.ndim == 2:
            data = hdu.data.astype(np.float64)
            break
    if data is None:
        data = hdulist[0].data
        if data is not None and data.ndim > 2:
            while data.ndim > 2:
                data = data[0]
            data = data.astype(np.float64)
    hdulist.close()

    if data is None:
        print("ERROR: No 2D image data found")
        return []

    data = np.ascontiguousarray(data, dtype=np.float64)
    if verbose:
        print(f"Image shape: {data.shape}")

    # Run SEP
    if verbose:
        print("Running SEP extraction...")
    sep_result = run_sep_on_fits(data)
    if sep_result is None:
        print("No sources detected")
        return []

    n = sep_result['n']
    if verbose:
        print(f"Detected {n} sources")

    # Find candidate pairs
    if verbose:
        print("Finding candidate pairs...")
    pairs = find_candidate_pairs(
        sep_result['objects'],
        sep_result['kron_radius'],
        sep_result['flux_auto'],
        alpha=candidate_cfg.alpha,
        radius_scale=candidate_cfg.radius_scale,
        min_flux=candidate_cfg.min_flux,
    )

    if verbose:
        print(f"Found {len(pairs)} candidate pairs")

    if not pairs:
        print("No candidate pairs found")
        return []

    # Extract features for all pairs first
    valid_pairs = []
    all_feats = []
    n_failed = 0
    norm_mean = norm_stats['mean']
    norm_std = np.where(norm_stats['std'] < 1e-8, 1.0, norm_stats['std'])

    if verbose:
        print("Extracting features...")
    n_pairs = len(pairs)
    for p_idx, (idx_i, idx_j) in enumerate(pairs):
        feats, _ = extract_pair_features(idx_i, idx_j, sep_result)
        if feats is None or not np.all(np.isfinite(feats)):
            n_failed += 1
            continue
        valid_pairs.append((idx_i, idx_j, feats))
        feats_norm = (feats - norm_mean) / norm_std
        all_feats.append(feats_norm)
        if verbose and (p_idx + 1) % 500 == 0:
            print(f"  Features: {p_idx+1}/{n_pairs} pairs processed", flush=True)

    # Batch CUDA inference
    if verbose:
        print(f"Running batch inference on {device} ({len(all_feats)} pairs)...")
    results = []

    if all_feats:
        feat_tensor = torch.FloatTensor(np.array(all_feats)).to(device)

        model.eval()
        with torch.no_grad():
            # Process in batches for memory efficiency
            batch_size = 1024
            all_proba = []
            for start in range(0, len(feat_tensor), batch_size):
                batch = feat_tensor[start:start + batch_size]
                if device.type == 'cuda':
                    with torch.amp.autocast('cuda'):
                        proba = torch.softmax(model(batch), dim=1)
                else:
                    proba = torch.softmax(model(batch), dim=1)
                all_proba.append(proba.cpu().numpy())
            all_proba = np.concatenate(all_proba, axis=0)

        obj = sep_result['objects']
        for k, (idx_i, idx_j, feats) in enumerate(valid_pairs):
            proba = all_proba[k]
            pred = int(proba.argmax())
            conf = float(proba.max())

            result = {
                'idx_i': idx_i,
                'idx_j': idx_j,
                'x_i': float(obj['x'][idx_i]),
                'y_i': float(obj['y'][idx_i]),
                'x_j': float(obj['x'][idx_j]),
                'y_j': float(obj['y'][idx_j]),
                'prediction': 'merge' if pred == 0 else 'separate',
                'pred_label': pred,
                'confidence': conf,
                'prob_merge': float(proba[0]),
                'prob_separate': float(proba[1]),
                'flux_i': float(sep_result['flux_auto'][idx_i]),
                'flux_j': float(sep_result['flux_auto'][idx_j]),
            }
            for m, name in enumerate(FEATURE_NAMES):
                result[f'feat_{name}'] = float(feats[m])
            results.append(result)

    if verbose:
        print(f"Predictions: {len(results)} pairs ({n_failed} failed)")

    # Summary
    n_merge = sum(1 for r in results if r['pred_label'] == 0)
    n_sep = sum(1 for r in results if r['pred_label'] == 1)
    if verbose:
        print(f"  Merge: {n_merge}, Separate: {n_sep}")

    # Build merge groups
    merge_pairs = [(r['idx_i'], r['idx_j']) for r in results if r['pred_label'] == 0]
    groups = build_merge_groups(merge_pairs, n)
    if verbose and groups:
        print(f"  Merge groups: {len(groups)}")
        for g_idx, g in enumerate(groups[:10]):
            print(f"    Group {g_idx+1}: sources {g}")

    # Save CSV
    csv_path = os.path.join(output_dir, "predictions.csv")
    _save_csv(results, csv_path, verbose)

    # Save summary
    summary_path = os.path.join(output_dir, "summary.txt")
    with open(summary_path, 'w') as f:
        f.write(f"FITS: {fits_path}\n")
        f.write(f"Sources detected: {n}\n")
        f.write(f"Candidate pairs: {len(pairs)}\n")
        f.write(f"Predictions: {len(results)}\n")
        f.write(f"  Merge: {n_merge}\n")
        f.write(f"  Separate: {n_sep}\n")
        f.write(f"Merge groups: {len(groups)}\n")
        for g_idx, g in enumerate(groups):
            f.write(f"  Group {g_idx+1}: sources {g}\n")

    if verbose:
        print(f"Results saved to {output_dir}/")

    return results


def _save_csv(results, path, verbose=True):
    """Save results to CSV."""
    if not results:
        return

    keys = list(results[0].keys())
    with open(path, 'w') as f:
        f.write(','.join(keys) + '\n')
        for r in results:
            vals = []
            for k in keys:
                v = r[k]
                if isinstance(v, float):
                    vals.append(f'{v:.6f}')
                else:
                    vals.append(str(v))
            f.write(','.join(vals) + '\n')

    if verbose:
        print(f"CSV saved: {path} ({len(results)} rows)")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Evaluate merge classifier on FITS')
    parser.add_argument('--fits', type=str, required=True,
                        help='Path to FITS image')
    parser.add_argument('--checkpoint', type=str, default=None,
                        help='Path to model checkpoint')
    parser.add_argument('--output-dir', type=str, default=None)
    parser.add_argument('--alpha', type=float, default=1.5)
    args = parser.parse_args()

    cfg = CandidateConfig(alpha=args.alpha)
    evaluate_fits(args.fits, checkpoint_path=args.checkpoint,
                  output_dir=args.output_dir, candidate_cfg=cfg)
