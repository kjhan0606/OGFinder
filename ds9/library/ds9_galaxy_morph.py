#!/usr/bin/env python3
"""
DS9 Galaxy Morphology Classification CLI.

Runs CNN-based morphology classification on extracted sources.
Outputs TSV to stdout for Tcl parsing in layout.tcl.

Usage:
    ds9_galaxy_morph.py FITSFILE --catalog catalog.tsv [--checkpoint path]
                        [--min-confidence 0.3] [--max-class-star 0.8]
                        [--min-a-image 2.0] [--cutout-size 64]
"""

import sys
import os
import argparse
import numpy as np

# Add project root to path
_script_dir = os.path.dirname(os.path.abspath(__file__))
_project_root = os.path.abspath(os.path.join(_script_dir, '..', '..'))
if _project_root not in sys.path:
    sys.path.insert(0, _project_root)

from galaxy_morph.config import GALAXY10_CLASSES, MORPH_COLORS, MorphConfig
from galaxy_morph.cutout.extract import batch_extract_cutouts
from galaxy_morph.training.evaluate import load_morph_model, predict_batch
from galaxy_morph.describe import build_descriptor


def parse_catalog_tsv(path):
    """Parse a TSV catalog file exported from ds9_sextract panel.

    Returns dict of numpy arrays with keys:
        number, x, y, a_image, class_star
    Coordinates: 1-based (ds9) -> 0-based.
    """
    required = ['NUMBER', 'X_IMAGE', 'Y_IMAGE']
    optional = ['A_IMAGE', 'B_IMAGE', 'CLASS_STAR', 'FLUX_AUTO']

    with open(path, 'r') as f:
        lines = f.readlines()

    if not lines:
        return None

    # Find header
    header = lines[0].strip().split('\t')
    col_idx = {}
    for col in required + optional:
        for i, h in enumerate(header):
            if h.strip() == col:
                col_idx[col] = i
                break

    missing_req = [c for c in required if c not in col_idx]
    if missing_req:
        print(f"ERROR: Missing required columns: {missing_req}", file=sys.stderr)
        return None

    # Parse data rows
    all_cols = [c for c in required + optional if c in col_idx]
    rows = {col: [] for col in all_cols}
    for line in lines[1:]:
        line = line.strip()
        if not line:
            continue
        fields = line.split('\t')
        try:
            for col in all_cols:
                val = float(fields[col_idx[col]])
                rows[col].append(val)
        except (IndexError, ValueError):
            continue

    if not rows.get('NUMBER'):
        return None

    n = len(rows['NUMBER'])
    result = {
        'number': np.array(rows['NUMBER'], dtype=np.int64),
        'x': np.array(rows['X_IMAGE']) - 1.0,  # 1-based -> 0-based
        'y': np.array(rows['Y_IMAGE']) - 1.0,
        'a_image': np.array(rows.get('A_IMAGE', [5.0] * n)),
        'class_star': np.array(rows.get('CLASS_STAR', [0.0] * n)),
        'flux_auto': np.array(rows.get('FLUX_AUTO', [1.0] * n)),
    }
    return result


def load_fits(fits_path):
    """Load FITS data, handling multi-dimensional cubes."""
    from astropy.io import fits as pyfits

    with pyfits.open(fits_path) as hdul:
        # Find image HDU
        data = None
        for hdu in hdul:
            if hdu.data is not None and hdu.data.ndim >= 2:
                data = hdu.data
                break

        if data is None:
            return None

        # Handle multi-dimensional data (take first 2D slice)
        while data.ndim > 2:
            data = data[0]

        return np.ascontiguousarray(data, dtype=np.float64)


def subtract_background(data):
    """Simple background subtraction using SEP."""
    try:
        import sep
    except ImportError:
        try:
            import sep_pjw as sep
        except ImportError:
            # Fallback: simple median subtraction
            print("WARNING: SEP not available, using median background", file=sys.stderr)
            return data - np.median(data)
    data = np.ascontiguousarray(data, dtype=np.float64)
    bkg = sep.Background(data, bw=64, bh=64, fw=3, fh=3)
    return data - bkg.back()


def main():
    parser = argparse.ArgumentParser(
        description='DS9 Galaxy Morphology Classification')
    parser.add_argument('fits', help='Path to FITS file')
    parser.add_argument('--catalog', required=True,
                        help='Path to catalog TSV from ds9_sextract')
    parser.add_argument('--checkpoint', default=None,
                        help='Path to CNN model checkpoint')
    parser.add_argument('--min-confidence', type=float, default=0.3,
                        help='Minimum classification confidence')
    parser.add_argument('--max-class-star', type=float, default=0.8,
                        help='Max CLASS_STAR (filter stars)')
    parser.add_argument('--min-a-image', type=float, default=2.0,
                        help='Min A_IMAGE (filter small sources)')
    parser.add_argument('--cutout-size', type=int, default=64,
                        help='Cutout size in pixels')
    args = parser.parse_args()

    # Find checkpoint
    checkpoint = args.checkpoint
    if checkpoint is None:
        # Default location
        checkpoint = os.path.join(_project_root, 'galaxy_morph', 'data',
                                   'checkpoints', 'cnn_morph_best.pt')
    if not os.path.exists(checkpoint):
        print(f"ERROR: Checkpoint not found: {checkpoint}", file=sys.stderr)
        print("Train the model first with: python3 -m galaxy_morph.training.train",
              file=sys.stderr)
        sys.exit(1)

    # Load FITS
    print(f"Loading FITS: {args.fits}", file=sys.stderr)
    data = load_fits(args.fits)
    if data is None:
        print("ERROR: Cannot load FITS data", file=sys.stderr)
        sys.exit(1)

    # Background subtraction
    print("Subtracting background ...", file=sys.stderr)
    data_sub = subtract_background(data)

    # Parse catalog
    print(f"Parsing catalog: {args.catalog}", file=sys.stderr)
    catalog = parse_catalog_tsv(args.catalog)
    if catalog is None:
        print("ERROR: Cannot parse catalog TSV", file=sys.stderr)
        sys.exit(1)

    n_sources = len(catalog['number'])
    print(f"  Total sources: {n_sources}", file=sys.stderr)

    # Filter: exclude stars and small sources
    galaxy_mask = np.ones(n_sources, dtype=bool)
    galaxy_mask &= catalog['class_star'] <= args.max_class_star
    galaxy_mask &= catalog['a_image'] >= args.min_a_image
    n_galaxies = np.sum(galaxy_mask)
    print(f"  Galaxy candidates (CLASS_STAR<={args.max_class_star}, "
          f"A_IMAGE>={args.min_a_image}): {n_galaxies}", file=sys.stderr)

    if n_galaxies == 0:
        # Output header with zero results
        print(f"#GALAXY_MORPH\tN_CLASSIFIED=0\tN_SOURCES={n_sources}")
        print("NUMBER\tMORPH_TYPE\tMORPH_DESC\tMORPH_CONF\tTOP1_CLASS\tTOP1_PROB\t"
              "TOP2_CLASS\tTOP2_PROB\tTOP3_CLASS\tTOP3_PROB\tMORPH_COLOR")
        sys.exit(0)

    # Build filtered catalog for cutout extraction
    galaxy_indices = np.where(galaxy_mask)[0]
    galaxy_catalog = {
        'x': catalog['x'][galaxy_indices],
        'y': catalog['y'][galaxy_indices],
        'a_image': catalog['a_image'][galaxy_indices],
    }

    # Extract cutouts
    print(f"Extracting {n_galaxies} cutouts ({args.cutout_size}x{args.cutout_size}) ...",
          file=sys.stderr)
    cutouts, valid_mask = batch_extract_cutouts(
        data_sub, galaxy_catalog, size=args.cutout_size)
    n_valid = int(np.sum(valid_mask))
    print(f"  Valid cutouts: {n_valid}", file=sys.stderr)

    if n_valid == 0:
        print(f"#GALAXY_MORPH\tN_CLASSIFIED=0\tN_SOURCES={n_sources}")
        print("NUMBER\tMORPH_TYPE\tMORPH_DESC\tMORPH_CONF\tTOP1_CLASS\tTOP1_PROB\t"
              "TOP2_CLASS\tTOP2_PROB\tTOP3_CLASS\tTOP3_PROB\tMORPH_COLOR")
        sys.exit(0)

    # Load model and run inference
    import torch
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    print(f"Loading model: {checkpoint} (device={device})", file=sys.stderr)
    model, norm_stats, class_names = load_morph_model(checkpoint, device)

    print(f"Running inference on {n_valid} cutouts ...", file=sys.stderr)
    proba = predict_batch(model, cutouts, norm_stats, device, batch_size=256)

    # Build results
    n_classified = 0
    results = []

    valid_count = 0
    for gi, galaxy_idx in enumerate(galaxy_indices):
        if not valid_mask[gi]:
            continue

        p = proba[valid_count]
        cutout = cutouts[valid_count]  # (1, H, W)
        valid_count += 1

        top3_idx = np.argsort(p)[::-1][:3]
        top1_id = top3_idx[0]
        top1_conf = float(p[top1_id])

        if top1_conf < args.min_confidence:
            continue

        n_classified += 1
        src_num = int(catalog['number'][galaxy_idx])
        color = MORPH_COLORS[top1_id]

        # Build compact morphological descriptor
        descriptor = build_descriptor(top1_id, p, cutout)

        results.append({
            'number': src_num,
            'morph_type': class_names[top1_id],
            'morph_conf': top1_conf,
            'morph_desc': descriptor,
            'top1_class': class_names[top3_idx[0]],
            'top1_prob': float(p[top3_idx[0]]),
            'top2_class': class_names[top3_idx[1]],
            'top2_prob': float(p[top3_idx[1]]),
            'top3_class': class_names[top3_idx[2]],
            'top3_prob': float(p[top3_idx[2]]),
            'color': color,
        })

    print(f"  Classified: {n_classified}/{n_valid} "
          f"(min_confidence={args.min_confidence})", file=sys.stderr)

    # Output TSV to stdout
    print(f"#GALAXY_MORPH\tN_CLASSIFIED={n_classified}\tN_SOURCES={n_sources}")
    print("NUMBER\tMORPH_TYPE\tMORPH_DESC\tMORPH_CONF\tTOP1_CLASS\tTOP1_PROB\t"
          "TOP2_CLASS\tTOP2_PROB\tTOP3_CLASS\tTOP3_PROB\tMORPH_COLOR")

    for r in results:
        print(f"{r['number']}\t{r['morph_type']}\t{r['morph_desc']}\t"
              f"{r['morph_conf']:.4f}\t"
              f"{r['top1_class']}\t{r['top1_prob']:.4f}\t"
              f"{r['top2_class']}\t{r['top2_prob']:.4f}\t"
              f"{r['top3_class']}\t{r['top3_prob']:.4f}\t"
              f"{r['color']}")


if __name__ == '__main__':
    main()
