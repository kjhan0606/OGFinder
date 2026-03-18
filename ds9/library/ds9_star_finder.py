#!/usr/bin/env python3
"""
DS9 AI Star/Galaxy Classification CLI.

Runs CNN-based binary star/galaxy classification on extracted sources.
Outputs TSV to stdout for Tcl parsing in layout.tcl.

Usage:
    ds9_star_finder.py FITSFILE --catalog catalog.tsv [--checkpoint path]
                       [--threshold 0.5] [--cutout-size 64] [--device cpu]
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

from star_finder.config import StarFinderConfig
from star_finder.cutout.extract import batch_extract_cutouts
from star_finder.training.evaluate import load_star_model, predict_batch


def parse_catalog_tsv(path):
    """Parse TSV catalog (NUMBER, X_IMAGE, Y_IMAGE required)."""
    required = ['NUMBER', 'X_IMAGE', 'Y_IMAGE']

    with open(path, 'r') as f:
        lines = f.readlines()

    if not lines:
        return None

    header = lines[0].strip().split('\t')
    col_idx = {}
    for col in required:
        for i, h in enumerate(header):
            if h.strip() == col:
                col_idx[col] = i
                break

    missing = [c for c in required if c not in col_idx]
    if missing:
        print(f"ERROR: Missing required columns: {missing}", file=sys.stderr)
        return None

    numbers = []
    xs = []
    ys = []
    for line in lines[1:]:
        line = line.strip()
        if not line:
            continue
        fields = line.split('\t')
        try:
            numbers.append(float(fields[col_idx['NUMBER']]))
            xs.append(float(fields[col_idx['X_IMAGE']]))
            ys.append(float(fields[col_idx['Y_IMAGE']]))
        except (IndexError, ValueError):
            continue

    if not numbers:
        return None

    return {
        'number': np.array(numbers, dtype=np.int64),
        'x': np.array(xs) - 1.0,  # 1-based -> 0-based
        'y': np.array(ys) - 1.0,
    }


def load_fits(fits_path):
    """Load FITS data, handling multi-dimensional cubes."""
    from astropy.io import fits as pyfits

    with pyfits.open(fits_path) as hdul:
        data = None
        for hdu in hdul:
            if hdu.data is not None and hdu.data.ndim >= 2:
                data = hdu.data
                break
        if data is None:
            return None
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
            print("WARNING: SEP not available, using median background",
                  file=sys.stderr)
            return data - np.median(data)
    data = np.ascontiguousarray(data, dtype=np.float64)
    bkg = sep.Background(data, bw=64, bh=64, fw=3, fh=3)
    return data - bkg.back()


def main():
    parser = argparse.ArgumentParser(
        description='DS9 AI Star/Galaxy Classification')
    parser.add_argument('fits', help='Path to FITS file')
    parser.add_argument('--catalog', required=True,
                        help='Path to catalog TSV from ds9_sextract')
    parser.add_argument('--checkpoint', default=None,
                        help='Path to CNN model checkpoint')
    parser.add_argument('--threshold', type=float, default=0.5,
                        help='Star classification threshold (default: 0.5)')
    parser.add_argument('--cutout-size', type=int, default=64,
                        help='Cutout size in pixels')
    parser.add_argument('--device', default='cpu',
                        help='Device: cpu or cuda')
    args = parser.parse_args()

    # Find checkpoint
    checkpoint = args.checkpoint
    if checkpoint is None:
        checkpoint = os.path.join(_project_root, 'star_finder', 'data',
                                   'checkpoints', 'star_finder_best.pt')
    if not os.path.exists(checkpoint):
        print(f"ERROR: Checkpoint not found: {checkpoint}", file=sys.stderr)
        print("Train the model first with: python3 -m star_finder.training.train",
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

    # Extract cutouts
    print(f"Extracting {n_sources} cutouts ({args.cutout_size}x{args.cutout_size}) ...",
          file=sys.stderr)
    cutouts, valid_mask = batch_extract_cutouts(
        data_sub, catalog, size=args.cutout_size)
    n_valid = int(np.sum(valid_mask))
    print(f"  Valid cutouts: {n_valid}", file=sys.stderr)

    if n_valid == 0:
        print(f"#STAR_FINDER\tN_CLASSIFIED=0\tN_SOURCES={n_sources}")
        print("NUMBER\tAI_STAR\tAI_STAR_CONF\tAI_STAR_COLOR")
        sys.exit(0)

    # Load model and run inference
    import torch
    device = torch.device(args.device if args.device != 'cpu' and
                          torch.cuda.is_available() else 'cpu')
    print(f"Loading model: {checkpoint} (device={device})", file=sys.stderr)
    model = load_star_model(checkpoint, device)

    print(f"Running inference on {n_valid} cutouts ...", file=sys.stderr)
    proba = predict_batch(model, cutouts, device, batch_size=256)

    # Build results
    n_stars = 0
    n_galaxies = 0
    results = []

    valid_count = 0
    for i in range(n_sources):
        if not valid_mask[i]:
            continue
        p = float(proba[valid_count])
        valid_count += 1

        src_num = int(catalog['number'][i])
        if p >= args.threshold:
            label = "star"
            color = "cyan"
            n_stars += 1
        else:
            label = "galaxy"
            color = "yellow"
            n_galaxies += 1

        results.append({
            'number': src_num,
            'label': label,
            'conf': p,
            'color': color,
        })

    n_classified = len(results)
    print(f"  Classified: {n_classified} (star={n_stars}, galaxy={n_galaxies})",
          file=sys.stderr)

    # Output TSV to stdout
    print(f"#STAR_FINDER\tN_CLASSIFIED={n_classified}\tN_SOURCES={n_sources}")
    print("NUMBER\tAI_STAR\tAI_STAR_CONF\tAI_STAR_COLOR")

    for r in results:
        print(f"{r['number']}\t{r['label']}\t{r['conf']:.4f}\t{r['color']}")


if __name__ == '__main__':
    main()
