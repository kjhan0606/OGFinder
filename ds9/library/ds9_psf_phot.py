#!/usr/bin/env python3
"""PSF photometry: fit PSF model to each source."""

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
            }
            sources.append(src)
        except (ValueError, IndexError, KeyError):
            continue

    return sources


def main():
    parser = argparse.ArgumentParser(description='PSF photometry')
    parser.add_argument('fitsfile', nargs='?', help='Input FITS image')
    parser.add_argument('--catalog', required=True, help='Input TSV catalog')
    parser.add_argument('--psf', required=True, help='PSF FITS image')
    parser.add_argument('--fit-radius', type=int, default=10,
                        help='Fitting radius in pixels (default: 10)')
    parser.add_argument('--mag-zeropoint', type=float, default=25.0)

    from parallel import add_batch_args
    add_batch_args(parser)

    args = parser.parse_args()

    try:
        from astropy.io import fits
    except ImportError:
        print("ERROR: astropy is required", file=sys.stderr)
        sys.exit(1)

    from psf_phot.config import PSFPhotConfig
    from psf_phot.photometry import do_psf_photometry

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

    print(f"PSF photometry for {len(sources)} sources...", file=sys.stderr)

    cfg = PSFPhotConfig(fit_radius=args.fit_radius, mag_zeropoint=args.mag_zeropoint)
    results = do_psf_photometry(data, psf, sources, cfg,
                                 n_workers=args.n_workers)

    # Output TSV
    print("NUMBER\tFLUX_PSF\tFLUXERR_PSF\tMAG_PSF\tMAGERR_PSF\tCHI2_PSF\tX_PSF\tY_PSF")
    for r in results:
        def fmt(v):
            return f"{v:.4f}" if np.isfinite(v) else "nan"
        print(f"{r['NUMBER']}\t{fmt(r['FLUX_PSF'])}\t{fmt(r['FLUXERR_PSF'])}\t"
              f"{fmt(r['MAG_PSF'])}\t{fmt(r['MAGERR_PSF'])}\t{fmt(r['CHI2_PSF'])}\t"
              f"{fmt(r['X_PSF'])}\t{fmt(r['Y_PSF'])}")

    print(f"Done: {len(results)} sources measured", file=sys.stderr)


if __name__ == '__main__':
    main()
