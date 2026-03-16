#!/usr/bin/env python3
"""Measure non-parametric morphology (C, A, Gini, M20, Petrosian) for sources."""

import sys
import os
import argparse
import numpy as np

# Add project root to path
_script_dir = os.path.dirname(os.path.abspath(__file__))
_project_root = os.path.abspath(os.path.join(_script_dir, '..', '..'))
if _project_root not in sys.path:
    sys.path.insert(0, _project_root)


def parse_catalog(catalog_file):
    """Parse TSV catalog and return list of source dicts."""
    sources = []
    with open(catalog_file, 'r') as f:
        lines = [l.strip() for l in f if l.strip()]

    if len(lines) < 2:
        return sources

    headers = [h.strip() for h in lines[0].split('\t')]

    # Find column indices
    col_map = {}
    for i, h in enumerate(headers):
        col_map[h] = i

    required = ['NUMBER', 'X_IMAGE', 'Y_IMAGE', 'A_IMAGE', 'B_IMAGE']
    for r in required:
        if r not in col_map:
            print(f"ERROR: Missing column {r} in catalog", file=sys.stderr)
            return sources

    for line in lines[1:]:
        fields = line.split('\t')
        try:
            src = {
                'number': int(fields[col_map['NUMBER']].strip()),
                'x': float(fields[col_map['X_IMAGE']].strip()) - 1.0,  # to 0-indexed
                'y': float(fields[col_map['Y_IMAGE']].strip()) - 1.0,
                'a': float(fields[col_map['A_IMAGE']].strip()),
                'b': float(fields[col_map['B_IMAGE']].strip()),
            }
            if 'THETA_IMAGE' in col_map:
                src['theta'] = float(fields[col_map['THETA_IMAGE']].strip()) * np.pi / 180.0
            if 'KRON_RADIUS' in col_map:
                src['kron_radius'] = float(fields[col_map['KRON_RADIUS']].strip())
            if 'NPIX_ISO' in col_map:
                src['npix'] = int(fields[col_map['NPIX_ISO']].strip())
            else:
                src['npix'] = 200  # default
            sources.append(src)
        except (ValueError, IndexError):
            continue

    return sources


def main():
    parser = argparse.ArgumentParser(description='Non-parametric morphometry (CAS/Gini/M20)')
    parser.add_argument('fitsfile', help='Input FITS image')
    parser.add_argument('--catalog', required=True, help='Input TSV catalog')
    parser.add_argument('--min-npix', type=int, default=100,
                        help='Minimum pixels for morphometry (default: 100)')
    args = parser.parse_args()

    try:
        import sep
    except ImportError:
        print("ERROR: sep is required", file=sys.stderr)
        sys.exit(1)

    try:
        from astropy.io import fits
    except ImportError:
        print("ERROR: astropy is required", file=sys.stderr)
        sys.exit(1)

    from morphometry.config import MorphometryConfig
    from morphometry.measure import measure_catalog

    # Load FITS
    with fits.open(args.fitsfile) as hdul:
        data = None
        for hdu in hdul:
            if hdu.data is not None and hdu.data.ndim >= 2:
                data = hdu.data.astype(np.float64)
                break
        if data is None:
            print("ERROR: No 2D image data found", file=sys.stderr)
            sys.exit(1)

    # Background subtraction
    mask = ~np.isfinite(data) | (data == 0)
    data[~np.isfinite(data)] = 0.0
    bkg = sep.Background(data, mask=mask.astype(np.uint8), bw=64, bh=64, fw=3, fh=3)
    data_sub = data - bkg

    # Parse catalog
    sources = parse_catalog(args.catalog)
    if not sources:
        print("ERROR: No sources in catalog", file=sys.stderr)
        sys.exit(1)

    print(f"Measuring morphometry for {len(sources)} sources...", file=sys.stderr)

    # Measure
    cfg = MorphometryConfig(min_npix=args.min_npix)
    results = measure_catalog(data_sub, sources, cfg)

    # Output TSV
    print("NUMBER\tCONC\tASYM\tGINI\tM20\tR_PETRO")
    for r in results:
        def fmt(v):
            return f"{v:.4f}" if np.isfinite(v) else "nan"
        print(f"{r['NUMBER']}\t{fmt(r['CONC'])}\t{fmt(r['ASYM'])}\t"
              f"{fmt(r['GINI'])}\t{fmt(r['M20'])}\t{fmt(r['R_PETRO'])}")

    print(f"Done: {len(results)} sources measured", file=sys.stderr)


if __name__ == '__main__':
    main()
