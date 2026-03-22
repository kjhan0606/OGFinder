#!/usr/bin/env python3
"""
DS9 Photo-z (Photometric Redshift) CLI.

Uses MLP+MDN model for photo-z estimation from multi-band photometry.
Falls back to color-based empirical estimation if no checkpoint available.

Usage:
    ds9_photo_z.py FITSFILE --catalog catalog.tsv
                   --bands g,r,i,z --mag-columns MAG_g,MAG_r,MAG_i,MAG_z
                   [--checkpoint photo_z/data/checkpoints/photo_z_mdn_best.pt]
                   [--n-workers 4]
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


def parse_catalog_tsv(path):
    """Parse TSV catalog."""
    with open(path, 'r') as f:
        lines = f.readlines()
    if not lines:
        return None, None

    header = lines[0].strip().split('\t')
    rows = []
    for line in lines[1:]:
        line = line.strip()
        if not line:
            continue
        rows.append(line.split('\t'))

    return header, rows


def get_column_index(header, name):
    """Find column index by name."""
    for i, h in enumerate(header):
        if h.strip() == name:
            return i
    return -1


def estimate_photoz_empirical(mags, mag_errs, band_names):
    """Empirical photo-z from color indices (no model needed).

    Uses simple color-redshift relations calibrated on SDSS data.
    This is a rough fallback when no trained model is available.
    """
    n_sources = mags.shape[0]
    n_bands = mags.shape[1]

    z_point = np.full(n_sources, np.nan)
    z_err = np.full(n_sources, np.nan)
    z_q68 = np.full(n_sources, np.nan)
    z_outlier = np.zeros(n_sources, dtype=int)

    for i in range(n_sources):
        m = mags[i]
        valid = np.isfinite(m) & (m < 90) & (m > 0)

        if np.sum(valid) < 2:
            continue

        # Compute all available colors
        colors = []
        for j in range(n_bands - 1):
            if valid[j] and valid[j + 1]:
                colors.append(m[j] - m[j + 1])

        if not colors:
            continue

        colors = np.array(colors)

        # Simple empirical relation: z ~ 0.3 * mean_color + 0.1
        # Calibrated roughly on SDSS ugriz
        mean_color = np.mean(colors)
        z_est = max(0.0, 0.3 * mean_color + 0.1)

        # Scatter from color variations
        if len(colors) > 1:
            color_scatter = np.std(colors)
        else:
            color_scatter = 0.3

        # Magnitude-based refinement: fainter -> likely higher z
        valid_mags = m[valid]
        mean_mag = np.mean(valid_mags)
        if mean_mag > 22:
            z_est += 0.1 * (mean_mag - 22)

        z_est = max(0.001, min(z_est, 6.0))
        err = max(0.02, 0.15 * (1 + z_est) + 0.1 * color_scatter)

        z_point[i] = round(z_est, 4)
        z_err[i] = round(err, 4)
        z_q68[i] = round(err * 0.85, 4)

        if err > 0.5 * (1 + z_est):
            z_outlier[i] = 1

    return z_point, z_err, z_q68, z_outlier


def estimate_photoz_mdn(mags, mag_errs, band_names, checkpoint):
    """Photo-z with trained MDN model."""
    try:
        from photo_z.features import extract_features
        from photo_z.predict import predict_photoz
        from photo_z.config import PhotoZConfig

        cfg = PhotoZConfig()
        features = extract_features(mags, mag_errs, band_names)
        result = predict_photoz(features, checkpoint, cfg)

        # Quality Head v2: flag outliers using P(reliable)
        p_reliable = result.get('p_reliable', None)
        if p_reliable is not None:
            z_outlier = np.where(p_reliable < 0.5, 1, 0)
        else:
            z_outlier = result.get('z_outlier', None)

        return (result['z_point'], result['z_err'], result['z_q68'],
                z_outlier, p_reliable)
    except Exception as e:
        print(f"WARNING: MDN model failed ({e}), falling back to empirical",
              file=sys.stderr)
        return (*estimate_photoz_empirical(mags, mag_errs, band_names), None)


def main():
    parser = argparse.ArgumentParser(description='DS9 Photo-z Estimation')
    parser.add_argument('fits', help='Path to FITS file')
    parser.add_argument('--catalog', required=True, help='Path to catalog TSV')
    parser.add_argument('--bands', default='g,r,i,z',
                        help='Band names (comma-separated)')
    parser.add_argument('--mag-columns', default='',
                        help='Magnitude column names (comma-separated)')
    parser.add_argument('--checkpoint', default=None,
                        help='Path to MDN checkpoint')
    parser.add_argument('--n-workers', type=int, default=1,
                        help='Number of workers (unused for photo-z)')
    args = parser.parse_args()

    # Parse catalog
    print(f"Loading catalog: {args.catalog}", file=sys.stderr)
    header, rows = parse_catalog_tsv(args.catalog)
    if header is None or not rows:
        print("ERROR: Cannot parse catalog", file=sys.stderr)
        sys.exit(1)

    num_idx = get_column_index(header, 'NUMBER')
    if num_idx < 0:
        print("ERROR: NUMBER column not found", file=sys.stderr)
        sys.exit(1)

    # Determine magnitude columns
    bands = [b.strip() for b in args.bands.split(',')]
    if args.mag_columns:
        mag_col_names = [c.strip() for c in args.mag_columns.split(',')]
    else:
        # Auto-detect: look for MAG_AUTO, MAG_APER, or MAG_{band}
        mag_col_names = []
        for b in bands:
            for pattern in [f'MAG_{b}', f'MAG_AUTO_{b}', f'mag_{b}']:
                idx = get_column_index(header, pattern)
                if idx >= 0:
                    mag_col_names.append(pattern)
                    break
        if not mag_col_names:
            # Single-band: use MAG_AUTO
            mag_auto = get_column_index(header, 'MAG_AUTO')
            if mag_auto >= 0:
                mag_col_names = ['MAG_AUTO']
                bands = ['auto']

    if not mag_col_names:
        print("ERROR: No magnitude columns found", file=sys.stderr)
        sys.exit(1)

    # Get column indices
    mag_indices = []
    for mc in mag_col_names:
        idx = get_column_index(header, mc)
        if idx < 0:
            print(f"WARNING: Column {mc} not found, skipping", file=sys.stderr)
        mag_indices.append(idx)

    # Build magnitude arrays
    n_sources = len(rows)
    n_bands = len(mag_col_names)
    mags = np.full((n_sources, n_bands), np.nan)
    mag_errs = np.full((n_sources, n_bands), np.nan)
    numbers = []

    # Check for error columns
    err_indices = []
    for mc in mag_col_names:
        err_name = mc.replace('MAG_', 'MAGERR_')
        err_indices.append(get_column_index(header, err_name))

    for i, row in enumerate(rows):
        numbers.append(row[num_idx])
        for j, mi in enumerate(mag_indices):
            if mi >= 0 and mi < len(row):
                try:
                    val = float(row[mi])
                    if val < 90 and val > 0:
                        mags[i, j] = val
                except (ValueError, IndexError):
                    pass
            if err_indices[j] >= 0 and err_indices[j] < len(row):
                try:
                    mag_errs[i, j] = float(row[err_indices[j]])
                except (ValueError, IndexError):
                    pass

    print(f"  Sources: {n_sources}, Bands: {n_bands} ({','.join(mag_col_names)})",
          file=sys.stderr)

    # Estimate photo-z
    p_reliable = None
    if args.checkpoint and os.path.exists(args.checkpoint):
        print(f"Using MDN model: {args.checkpoint}", file=sys.stderr)
        result = estimate_photoz_mdn(
            mags, mag_errs, bands, args.checkpoint)
        z_point, z_err, z_q68, z_outlier, p_reliable = result
    else:
        print("Using empirical color-z relation (no checkpoint)", file=sys.stderr)
        z_point, z_err, z_q68, z_outlier = estimate_photoz_empirical(
            mags, mag_errs, bands)

    # Count valid
    n_valid = np.sum(np.isfinite(z_point))
    print(f"  Valid photo-z: {n_valid}/{n_sources}", file=sys.stderr)
    if p_reliable is not None:
        n_reliable = np.sum(p_reliable > 0.5)
        print(f"  Reliable (P>0.5): {n_reliable}/{n_sources}", file=sys.stderr)

    # Output TSV
    has_quality = p_reliable is not None
    print(f"#PHOTO_Z\tN_ESTIMATED={n_valid}\tN_SOURCES={n_sources}")
    if has_quality:
        print("NUMBER\tPHOTO_Z\tPHOTO_Z_ERR\tPHOTO_Z_Q68\t"
              "PHOTO_Z_OUTLIER\tPHOTO_Z_RELIABLE")
    else:
        print("NUMBER\tPHOTO_Z\tPHOTO_Z_ERR\tPHOTO_Z_Q68\tPHOTO_Z_OUTLIER")

    for i in range(n_sources):
        num = numbers[i]
        zp = f"{z_point[i]:.4f}" if np.isfinite(z_point[i]) else "-99"
        ze = f"{z_err[i]:.4f}" if np.isfinite(z_err[i]) else "-99"
        zq = f"{z_q68[i]:.4f}" if np.isfinite(z_q68[i]) else "-99"
        zo = str(z_outlier[i]) if z_outlier is not None else "0"
        if has_quality:
            pr = f"{p_reliable[i]:.4f}"
            print(f"{num}\t{zp}\t{ze}\t{zq}\t{zo}\t{pr}")
        else:
            print(f"{num}\t{zp}\t{ze}\t{zq}\t{zo}")


if __name__ == '__main__':
    main()
