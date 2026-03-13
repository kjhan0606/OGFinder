#!/usr/bin/env python3
"""
DS9 AI Merge CLI — runs AI merge prediction on a FITS file.

Outputs TSV to stdout for Tcl parsing in layout.tcl.

Usage:
    ds9_ai_merge.py FITSFILE [--checkpoint PATH] [--threshold 0.7]
                    [--detect-thresh 1.5] [--detect-minarea 5]
                    [--deblend-nthresh 32] [--deblend-mincont 0.005]
                    [--mag-zeropoint 25.0] [--back-size 64]
                    [--back-filtersize 3] [--alpha 1.5]
                    [--radius-scale 2.0] [--min-flux 0.0]
"""

import sys
import os
import argparse
import numpy as np

# Add ai_merge parent to path
_script_dir = os.path.dirname(os.path.abspath(__file__))
_project_root = os.path.abspath(os.path.join(_script_dir, '..', '..'))
if _project_root not in sys.path:
    sys.path.insert(0, _project_root)

from ai_merge.utils.sep_runner import run_sep_on_fits
from ai_merge.catalog.find_candidates import find_candidate_pairs, build_merge_groups
from ai_merge.catalog.extract_features import extract_pair_features


def main():
    parser = argparse.ArgumentParser(description='DS9 AI Merge prediction')
    parser.add_argument('fits', help='Path to FITS file')
    parser.add_argument('--checkpoint', default=None,
                        help='Path to model checkpoint (.pt)')
    parser.add_argument('--threshold', type=float, default=0.7,
                        help='Confidence threshold for merge prediction')
    # SEP parameters
    parser.add_argument('--detect-thresh', type=float, default=1.5)
    parser.add_argument('--detect-minarea', type=int, default=5)
    parser.add_argument('--deblend-nthresh', type=int, default=32)
    parser.add_argument('--deblend-mincont', type=float, default=0.005)
    parser.add_argument('--mag-zeropoint', type=float, default=25.0)
    parser.add_argument('--back-size', type=int, default=64)
    parser.add_argument('--back-filtersize', type=int, default=3)
    # Candidate search parameters
    parser.add_argument('--alpha', type=float, default=1.5)
    parser.add_argument('--radius-scale', type=float, default=2.0)
    parser.add_argument('--min-flux', type=float, default=0.0)

    args = parser.parse_args()

    # --- Load FITS ---
    try:
        from astropy.io import fits
    except ImportError:
        print("ERROR: astropy not installed", file=sys.stderr)
        sys.exit(1)

    if not os.path.exists(args.fits):
        print(f"ERROR: File not found: {args.fits}", file=sys.stderr)
        sys.exit(1)

    print("Loading FITS...", file=sys.stderr)
    hdulist = fits.open(args.fits)
    data = None
    for hdu in hdulist:
        if hdu.data is not None and hdu.data.ndim == 2:
            data = hdu.data.astype(np.float64)
            break
    if data is None:
        d = hdulist[0].data
        if d is not None and d.ndim > 2:
            while d.ndim > 2:
                d = d[0]
            data = d.astype(np.float64)
    hdulist.close()

    if data is None:
        print("ERROR: No 2D image data found", file=sys.stderr)
        sys.exit(1)

    data = np.ascontiguousarray(data, dtype=np.float64)
    print(f"Image: {data.shape[1]}x{data.shape[0]}", file=sys.stderr)

    # --- Run SEP ---
    print("Running SEP extraction...", file=sys.stderr)
    sep_result = run_sep_on_fits(
        data,
        back_size=args.back_size,
        back_filtersize=args.back_filtersize,
        detect_thresh=args.detect_thresh,
        detect_minarea=args.detect_minarea,
        deblend_nthresh=args.deblend_nthresh,
        deblend_mincont=args.deblend_mincont,
        mag_zeropoint=args.mag_zeropoint,
    )
    if sep_result is None:
        print("ERROR: No sources detected", file=sys.stderr)
        sys.exit(1)

    n = sep_result['n']
    print(f"Detected {n} sources", file=sys.stderr)

    # --- Find candidate pairs ---
    print("Finding candidate pairs...", file=sys.stderr)
    pairs = find_candidate_pairs(
        sep_result['objects'],
        sep_result['kron_radius'],
        sep_result['flux_auto'],
        alpha=args.alpha,
        radius_scale=args.radius_scale,
        min_flux=args.min_flux,
    )
    print(f"Found {len(pairs)} candidate pairs", file=sys.stderr)

    if not pairs:
        # Output empty header
        print(f"#AI_MERGE\tN_GROUPS=0\tTHRESHOLD={args.threshold:.2f}")
        sys.exit(0)

    # --- Extract features ---
    print("Extracting features...", file=sys.stderr)
    valid_pairs = []
    all_feats = []
    for idx_i, idx_j in pairs:
        feats, _ = extract_pair_features(idx_i, idx_j, sep_result)
        if feats is not None and np.all(np.isfinite(feats)):
            valid_pairs.append((idx_i, idx_j))
            all_feats.append(feats)

    if not valid_pairs:
        print(f"#AI_MERGE\tN_GROUPS=0\tTHRESHOLD={args.threshold:.2f}")
        sys.exit(0)

    print(f"Valid pairs: {len(valid_pairs)}", file=sys.stderr)

    # --- Load model and predict ---
    try:
        import torch
    except ImportError:
        print("ERROR: torch not installed", file=sys.stderr)
        sys.exit(1)

    from ai_merge.training.evaluate import load_model

    checkpoint = args.checkpoint
    if checkpoint is None:
        # Default: look relative to project root
        checkpoint = os.path.join(_project_root, 'ai_merge', 'data',
                                  'checkpoints', 'mlp_best.pt')
    if not os.path.exists(checkpoint):
        print(f"ERROR: Checkpoint not found: {checkpoint}", file=sys.stderr)
        sys.exit(1)

    print("Loading model...", file=sys.stderr)
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    model, norm_stats = load_model(checkpoint, device)

    norm_mean = norm_stats['mean']
    norm_std = np.where(norm_stats['std'] < 1e-8, 1.0, norm_stats['std'])

    # Normalize and batch predict
    feat_array = np.array(all_feats)
    feat_norm = (feat_array - norm_mean) / norm_std
    feat_tensor = torch.FloatTensor(feat_norm).to(device)

    print(f"Running inference on {device}...", file=sys.stderr)
    model.eval()
    with torch.no_grad():
        batch_size = 1024
        all_proba = []
        for start in range(0, len(feat_tensor), batch_size):
            batch = feat_tensor[start:start + batch_size]
            proba = torch.softmax(model(batch), dim=1)
            all_proba.append(proba.cpu().numpy())
        all_proba = np.concatenate(all_proba, axis=0)

    # --- Filter merge predictions above threshold ---
    merge_pairs = []
    merge_confs = []
    obj = sep_result['objects']
    threshold = args.threshold

    for k, (idx_i, idx_j) in enumerate(valid_pairs):
        prob_merge = float(all_proba[k, 0])
        if prob_merge >= threshold:
            merge_pairs.append((idx_i, idx_j))
            merge_confs.append(prob_merge)

    print(f"Merge predictions above threshold: {len(merge_pairs)}", file=sys.stderr)

    if not merge_pairs:
        print(f"#AI_MERGE\tN_GROUPS=0\tTHRESHOLD={threshold:.2f}")
        sys.exit(0)

    # --- Build merge groups ---
    groups = build_merge_groups(merge_pairs, n)
    print(f"Merge groups: {len(groups)}", file=sys.stderr)

    # Compute per-group confidence (average of pairwise confidences for members)
    pair_conf_map = {}
    for k, (idx_i, idx_j) in enumerate(merge_pairs):
        pair_conf_map[(min(idx_i, idx_j), max(idx_i, idx_j))] = merge_confs[k]

    # --- Output TSV ---
    print(f"#AI_MERGE\tN_GROUPS={len(groups)}\tTHRESHOLD={threshold:.2f}")
    print("GROUP\tN_MEMBERS\tCONFIDENCE\tMEMBERS_X\tMEMBERS_Y")

    for g_idx, group in enumerate(groups):
        n_members = len(group)
        # Gather coordinates
        xs = [f"{obj['x'][m]:.2f}" for m in group]
        ys = [f"{obj['y'][m]:.2f}" for m in group]

        # Average confidence of all pairs within this group
        confs = []
        for a in range(len(group)):
            for b in range(a + 1, len(group)):
                key = (min(group[a], group[b]), max(group[a], group[b]))
                if key in pair_conf_map:
                    confs.append(pair_conf_map[key])
        avg_conf = np.mean(confs) if confs else 0.0

        print(f"{g_idx+1}\t{n_members}\t{avg_conf:.4f}\t{','.join(xs)}\t{','.join(ys)}")

    print("Done.", file=sys.stderr)


if __name__ == '__main__':
    main()
