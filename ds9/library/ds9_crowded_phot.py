#!/usr/bin/env python3
"""Crowded field photometry: DAOPHOT-style iterative multi-PSF fitting."""

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
    col_map = {h: i for i, h in enumerate(headers)}

    for line in lines[1:]:
        fields = line.split('\t')
        try:
            src = {
                'number': int(fields[col_map['NUMBER']].strip()),
                'x': float(fields[col_map['X_IMAGE']].strip()) - 1.0,
                'y': float(fields[col_map['Y_IMAGE']].strip()) - 1.0,
                'a': float(fields[col_map.get('A_IMAGE', -1)].strip()) if 'A_IMAGE' in col_map else 2.0,
                'kron_radius': float(fields[col_map.get('KRON_RADIUS', -1)].strip()) if 'KRON_RADIUS' in col_map else 3.5,
            }
            sources.append(src)
        except (ValueError, IndexError, KeyError):
            continue

    return sources


def main():
    parser = argparse.ArgumentParser(description='Crowded field photometry')
    parser.add_argument('fitsfile', help='Input FITS image')
    parser.add_argument('--catalog', required=True, help='Input TSV catalog')
    parser.add_argument('--psf', required=True, help='PSF FITS image')
    parser.add_argument('--max-iter', type=int, default=3,
                        help='Maximum iterations (default: 3)')
    parser.add_argument('--mag-zeropoint', type=float, default=25.0)
    args = parser.parse_args()

    try:
        from astropy.io import fits
    except ImportError:
        print("ERROR: astropy is required", file=sys.stderr)
        sys.exit(1)

    from crowded_phot.config import CrowdedPhotConfig
    from crowded_phot.iterate import crowded_photometry

    # Load FITS image
    with fits.open(args.fitsfile) as hdul:
        data = None
        for hdu in hdul:
            if hdu.data is not None and hdu.data.ndim >= 2:
                data = hdu.data.astype(np.float64)
                break
    if data is None:
        print("ERROR: No 2D image data found", file=sys.stderr)
        sys.exit(1)

    # Load PSF
    with fits.open(args.psf) as hdul:
        psf = None
        for hdu in hdul:
            if hdu.data is not None and hdu.data.ndim >= 2:
                psf = hdu.data.astype(np.float64)
                break
    if psf is None:
        print("ERROR: No PSF data found", file=sys.stderr)
        sys.exit(1)

    # Parse catalog
    sources = parse_catalog(args.catalog)
    if not sources:
        print("ERROR: No sources in catalog", file=sys.stderr)
        sys.exit(1)

    print(f"Crowded field photometry for {len(sources)} sources ({args.max_iter} iterations)...",
          file=sys.stderr)

    cfg = CrowdedPhotConfig(max_iterations=args.max_iter,
                             mag_zeropoint=args.mag_zeropoint)
    results = crowded_photometry(data, psf, sources, cfg)

    # Output TSV
    print("NUMBER\tFLUX_CROWD\tFLUXERR_CROWD\tMAG_CROWD\tX_CROWD\tY_CROWD\tN_NEIGHBORS")
    for r in sorted(results, key=lambda x: x.get('MAG_CROWD', 99.0)):
        def fmt(v):
            return f"{v:.4f}" if np.isfinite(v) else "nan"
        print(f"{r['NUMBER']}\t{fmt(r['FLUX_CROWD'])}\t{fmt(r['FLUXERR_CROWD'])}\t"
              f"{fmt(r['MAG_CROWD'])}\t{fmt(r['X_CROWD'])}\t{fmt(r['Y_CROWD'])}\t"
              f"{r['N_NEIGHBORS']}")

    print(f"Done: {len(results)} sources measured", file=sys.stderr)


if __name__ == '__main__':
    main()
