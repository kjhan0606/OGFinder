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

from ai_merge.utils.sep_runner import run_sep_on_fits, run_sep_light
from ai_merge.catalog.find_candidates import find_candidate_pairs, build_merge_groups
from ai_merge.catalog.extract_features import extract_pair_features


def parse_catalog_tsv(path):
    """Parse a TSV catalog file exported from ds9_sextract panel.

    Returns dict of numpy arrays with keys:
        number, x, y, a, b, theta, flux_auto, kron_radius,
        flux_radius, ellipticity, class_star
    Coordinates are converted from 1-based (ds9) to 0-based (SEP).
    Theta is converted from degrees to radians.
    """
    required = ['NUMBER', 'X_IMAGE', 'Y_IMAGE', 'A_IMAGE', 'B_IMAGE',
                'THETA_IMAGE', 'FLUX_AUTO', 'KRON_RADIUS', 'FLUX_RADIUS',
                'ELLIPTICITY', 'CLASS_STAR']

    with open(path, 'r') as f:
        lines = f.readlines()

    if not lines:
        return None

    # Find header
    header = lines[0].strip().split('\t')
    col_idx = {}
    for col in required:
        for i, h in enumerate(header):
            if h.strip() == col:
                col_idx[col] = i
                break

    missing = [c for c in required if c not in col_idx]
    if missing:
        print(f"WARNING: Missing catalog columns: {missing}", file=sys.stderr)
        # NUMBER is critical
        if 'NUMBER' not in col_idx or 'X_IMAGE' not in col_idx or 'Y_IMAGE' not in col_idx:
            return None

    # Parse data rows
    rows = {col: [] for col in required if col in col_idx}
    for line in lines[1:]:
        line = line.strip()
        if not line:
            continue
        fields = line.split('\t')
        try:
            for col in rows:
                val = float(fields[col_idx[col]])
                rows[col].append(val)
        except (IndexError, ValueError):
            continue

    if not rows.get('NUMBER'):
        return None

    result = {
        'number': np.array(rows['NUMBER'], dtype=np.int64),
        'x': np.array(rows['X_IMAGE']) - 1.0,  # 1-based → 0-based
        'y': np.array(rows['Y_IMAGE']) - 1.0,
        'a': np.array(rows.get('A_IMAGE', [0.5] * len(rows['NUMBER']))),
        'b': np.array(rows.get('B_IMAGE', [0.5] * len(rows['NUMBER']))),
        'theta': np.deg2rad(np.array(rows.get('THETA_IMAGE',
                                               [0.0] * len(rows['NUMBER'])))),
        'flux_auto': np.array(rows.get('FLUX_AUTO',
                                       [1.0] * len(rows['NUMBER']))),
        'kron_radius': np.array(rows.get('KRON_RADIUS',
                                         [3.5] * len(rows['NUMBER']))),
        'flux_radius': np.array(rows.get('FLUX_RADIUS',
                                         [2.0] * len(rows['NUMBER']))),
        'ellipticity': np.array(rows.get('ELLIPTICITY',
                                         [0.0] * len(rows['NUMBER']))),
        'class_star': np.array(rows.get('CLASS_STAR',
                                        [0.5] * len(rows['NUMBER']))),
    }
    # Filter out sources with NaN/Inf in critical fields
    valid = (np.isfinite(result['x']) & np.isfinite(result['y']) &
             np.isfinite(result['a']) & np.isfinite(result['b']))
    if not np.all(valid):
        n_bad = int(np.sum(~valid))
        print(f"WARNING: Dropping {n_bad} sources with NaN/Inf coordinates",
              file=sys.stderr)
        for key in result:
            result[key] = result[key][valid]

    # Sanitize
    result['a'] = np.maximum(result['a'], 0.5)
    result['b'] = np.clip(result['b'], 0.5, result['a'])
    return result


def map_catalog_to_segmap(catalog, segmap):
    """Map catalog sources to SEP segmap labels. Returns remapped segmap.

    For each catalog source, looks up the SEP segmap label at its (x,y)
    position, then relabels that region as k+1 (catalog index + 1).
    """
    h, w = segmap.shape
    n = len(catalog['x'])
    new_segmap = np.zeros_like(segmap)

    for k in range(n):
        if not (np.isfinite(catalog['x'][k]) and np.isfinite(catalog['y'][k])):
            continue
        ix = int(round(catalog['x'][k]))
        iy = int(round(catalog['y'][k]))
        if 0 <= ix < w and 0 <= iy < h:
            old_label = segmap[iy, ix]
            if old_label > 0:
                new_segmap[segmap == old_label] = k + 1

    return new_segmap


def build_hybrid_sep_result(catalog, sep_light):
    """Build a sep_result dict combining catalog photometry with SEP segmap.

    Parameters
    ----------
    catalog : dict from parse_catalog_tsv()
    sep_light : dict from run_sep_light()

    Returns
    -------
    dict compatible with extract_pair_features() and find_candidate_pairs()
    """
    n = len(catalog['x'])
    data_sub = sep_light['data_sub']
    h, w = data_sub.shape

    # Build structured array matching SEP objects interface
    dtype = np.dtype([('x', np.float64), ('y', np.float64),
                      ('a', np.float64), ('b', np.float64),
                      ('theta', np.float64), ('peak', np.float64)])
    objects = np.zeros(n, dtype=dtype)
    objects['x'] = catalog['x']
    objects['y'] = catalog['y']
    objects['a'] = catalog['a']
    objects['b'] = catalog['b']
    objects['theta'] = catalog['theta']

    # Peak: 3x3 max around each source position in data_sub
    for k in range(n):
        if not (np.isfinite(catalog['x'][k]) and np.isfinite(catalog['y'][k])):
            objects['peak'][k] = 0.0
            continue
        ix = int(round(catalog['x'][k]))
        iy = int(round(catalog['y'][k]))
        y0 = max(0, iy - 1)
        y1 = min(h, iy + 2)
        x0 = max(0, ix - 1)
        x1 = min(w, ix + 2)
        if y1 > y0 and x1 > x0:
            objects['peak'][k] = data_sub[y0:y1, x0:x1].max()
        else:
            objects['peak'][k] = 0.0

    # Remap segmap to catalog indices
    new_segmap = map_catalog_to_segmap(catalog, sep_light['segmap'])

    return {
        'objects': objects,
        'segmap': new_segmap,
        'n': n,
        'flux_auto': catalog['flux_auto'],
        'kron_radius': catalog['kron_radius'],
        'flux_radius': catalog['flux_radius'],
        'ellipticity': catalog['ellipticity'],
        'class_star': catalog['class_star'],
        'data_sub': data_sub,
        'bkg_rms': sep_light['bkg_rms'],
    }


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
    # Catalog mode
    parser.add_argument('--catalog', default=None,
                        help='Path to catalog TSV from ds9_sextract panel')
    # Candidate search parameters
    parser.add_argument('--alpha', type=float, default=1.1)
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

    # --- Run SEP or use catalog ---
    catalog = None
    if args.catalog:
        print(f"Loading catalog from {args.catalog}...", file=sys.stderr)
        catalog = parse_catalog_tsv(args.catalog)
        if catalog is None:
            print("ERROR: Failed to parse catalog TSV", file=sys.stderr)
            sys.exit(1)
        n_cat = len(catalog['x'])
        print(f"Catalog: {n_cat} sources", file=sys.stderr)

        print("Running lightweight SEP (segmap only)...", file=sys.stderr)
        sep_light = run_sep_light(
            data,
            back_size=args.back_size,
            back_filtersize=args.back_filtersize,
            detect_thresh=args.detect_thresh,
            detect_minarea=args.detect_minarea,
            deblend_nthresh=args.deblend_nthresh,
            deblend_mincont=args.deblend_mincont,
        )

        print("Building hybrid SEP result...", file=sys.stderr)
        sep_result = build_hybrid_sep_result(catalog, sep_light)
        n = sep_result['n']
        print(f"Hybrid result: {n} sources", file=sys.stderr)
    else:
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
    print("GROUP\tN_MEMBERS\tCONFIDENCE\tMEMBERS_X\tMEMBERS_Y\tMEMBERS_NUM")

    for g_idx, group in enumerate(groups):
        n_members = len(group)
        # Gather coordinates (1-based for ds9 display)
        xs = [f"{obj['x'][m] + 1.0:.2f}" for m in group]
        ys = [f"{obj['y'][m] + 1.0:.2f}" for m in group]

        # Member NUMBERs from catalog
        if catalog is not None:
            nums = [str(catalog['number'][m]) for m in group]
        else:
            nums = [str(m + 1) for m in group]

        # Average confidence of all pairs within this group
        confs = []
        for a in range(len(group)):
            for b in range(a + 1, len(group)):
                key = (min(group[a], group[b]), max(group[a], group[b]))
                if key in pair_conf_map:
                    confs.append(pair_conf_map[key])
        avg_conf = np.mean(confs) if confs else 0.0

        print(f"{g_idx+1}\t{n_members}\t{avg_conf:.4f}\t{','.join(xs)}\t{','.join(ys)}\t{','.join(nums)}")

    print("Done.", file=sys.stderr)


if __name__ == '__main__':
    main()
