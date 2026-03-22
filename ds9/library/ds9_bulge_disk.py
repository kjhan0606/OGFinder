#!/usr/bin/env python3
"""DS9 2D Bulge+Disk Decomposition CLI.

Fits each source with a 2-component model: Bulge (Sérsic n~4) + Disk (Sérsic n~1).
Outputs B/T ratio, component magnitudes, and effective radii.

Usage:
    ds9_bulge_disk.py FITSFILE --catalog cat.tsv [--psf psf.fits]
                      [--mag-zeropoint 25.0] [--pixel-scale 0.263]
                      [--max-sources 100] [--free-bulge-n] [--n-workers 4]
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


def _parse_sources(catalog_file, max_sources):
    """Parse catalog TSV and return list of source info dicts."""
    with open(catalog_file, 'r') as f:
        lines = [l.strip() for l in f if l.strip()]
    if not lines:
        return []

    headers = [h.strip() for h in lines[0].split('\t')]
    col_map = {h: i for i, h in enumerate(headers)}

    required = ['NUMBER', 'X_IMAGE', 'Y_IMAGE']
    for col in required:
        if col not in col_map:
            print(f"ERROR: Missing required column {col}", file=sys.stderr)
            return []

    sources = []
    for line in lines[1:]:
        if len(sources) >= max_sources:
            break
        fields = line.split('\t')
        try:
            src = {
                'num': int(fields[col_map['NUMBER']].strip()),
                'x': float(fields[col_map['X_IMAGE']].strip()) - 1.0,
                'y': float(fields[col_map['Y_IMAGE']].strip()) - 1.0,
                'a': float(fields[col_map['A_IMAGE']].strip()) if 'A_IMAGE' in col_map else 3.0,
                'b': float(fields[col_map['B_IMAGE']].strip()) if 'B_IMAGE' in col_map else 2.0,
                'theta': (float(fields[col_map['THETA_IMAGE']].strip()) * np.pi / 180.0
                          if 'THETA_IMAGE' in col_map else 0.0),
                'kron_r': (float(fields[col_map['KRON_RADIUS']].strip())
                           if 'KRON_RADIUS' in col_map else 3.5),
                'flux': (float(fields[col_map['FLUX_AUTO']].strip())
                         if 'FLUX_AUTO' in col_map else 1000.0),
            }
            sources.append(src)
        except (ValueError, IndexError, KeyError):
            continue

    return sources


def load_fits(fits_path):
    """Load FITS data as float64."""
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
    """Background subtraction using SEP."""
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
    mask = ~np.isfinite(data) | (data == 0)
    data[~np.isfinite(data)] = 0.0
    bkg = sep.Background(data, mask=mask.astype(np.uint8), bw=64, bh=64)
    return data - bkg.back()


def process_file(fitsfile, catalog_file=None, psf_path=None,
                 mag_zeropoint=25.0, pixel_scale=0.263, max_sources=100,
                 free_bulge_n=False, n_workers=0):
    """Process a single FITS file and return TSV string."""
    from bulge_disk.config import BulgeDiskConfig
    from bulge_disk.psf_convolve import load_psf
    from bulge_disk.decompose import decompose_catalog
    from parallel import resolve_n_workers

    # Load FITS
    print(f"Loading FITS: {fitsfile}", file=sys.stderr)
    data = load_fits(fitsfile)
    if data is None:
        print("ERROR: Cannot load FITS data", file=sys.stderr)
        return ""

    # Background subtraction
    print("Subtracting background ...", file=sys.stderr)
    data_sub = subtract_background(data)

    # Parse catalog
    print(f"Parsing catalog: {catalog_file}", file=sys.stderr)
    sources = _parse_sources(catalog_file, max_sources)
    if not sources:
        print("ERROR: No sources in catalog", file=sys.stderr)
        return ""

    n_total = len(sources)
    print(f"  Sources to decompose: {n_total}", file=sys.stderr)

    # Load PSF
    psf = load_psf(psf_path)
    if psf is not None:
        print(f"  PSF shape: {psf.shape}", file=sys.stderr)

    # Config
    cfg = BulgeDiskConfig(
        max_sources=max_sources,
        mag_zeropoint=mag_zeropoint,
        pixel_scale=pixel_scale,
        free_bulge_n=free_bulge_n,
    )

    # Decompose
    print(f"Bulge+Disk decomposition ({n_total} sources) ...", file=sys.stderr)
    results = decompose_catalog(data_sub, sources, psf=psf, cfg=cfg,
                                 n_workers=n_workers)

    # Count successes
    n_fitted = sum(1 for r in results if r['flag'] == 0)

    # Build TSV output
    lines = []
    lines.append(f"#BULGE_DISK\tN_FITTED={n_fitted}\tN_SOURCES={n_total}")
    lines.append("NUMBER\tBT_RATIO\tBULGE_RE\tBULGE_MAG\tDISK_RS\tDISK_MAG\tBD_CHI2\tBD_FLAG")

    for r in results:
        def fmt(v):
            return f"{v:.4f}" if np.isfinite(v) else "nan"
        def fmt3(v):
            return f"{v:.3f}" if np.isfinite(v) else "nan"
        def fmt2(v):
            return f"{v:.2f}" if np.isfinite(v) else "nan"

        lines.append(
            f"{r['num']}\t{fmt(r['bt_ratio'])}\t{fmt3(r['bulge_re'])}\t"
            f"{fmt2(r['bulge_mag'])}\t{fmt3(r['disk_rs'])}\t"
            f"{fmt2(r['disk_mag'])}\t{fmt(r['chi2'])}\t{r['flag']}"
        )

    print(f"Done: {n_fitted}/{n_total} sources fitted successfully",
          file=sys.stderr)
    return '\n'.join(lines)


def main():
    parser = argparse.ArgumentParser(
        description='DS9 2D Bulge+Disk Decomposition')
    parser.add_argument('fitsfile', nargs='?', help='Input FITS image')
    parser.add_argument('--catalog', required=True, help='Input TSV catalog')
    parser.add_argument('--psf', default=None, help='Path to PSF FITS file')
    parser.add_argument('--mag-zeropoint', type=float, default=25.0,
                        help='Magnitude zeropoint (default: 25.0)')
    parser.add_argument('--pixel-scale', type=float, default=0.263,
                        help='Pixel scale in arcsec/pixel (default: 0.263)')
    parser.add_argument('--max-sources', type=int, default=100,
                        help='Maximum sources to fit (default: 100)')
    parser.add_argument('--free-bulge-n', action='store_true',
                        help='Allow bulge Sérsic index to vary')

    from parallel import add_batch_args
    add_batch_args(parser)

    args = parser.parse_args()

    # Batch mode
    if args.batch is not None:
        from parallel import batch_run
        if not args.batch:
            print("ERROR: --batch requires FITS file list", file=sys.stderr)
            sys.exit(1)
        batch_run(args.batch,
                  lambda f, n_workers=0, **kw: process_file(
                      f, catalog_file=args.catalog,
                      psf_path=args.psf,
                      mag_zeropoint=args.mag_zeropoint,
                      pixel_scale=args.pixel_scale,
                      max_sources=args.max_sources,
                      free_bulge_n=args.free_bulge_n,
                      n_workers=n_workers),
                  n_workers=args.n_workers,
                  output_dir=args.output_dir)
        return

    if not args.fitsfile:
        print("ERROR: fitsfile required (or use --batch)", file=sys.stderr)
        sys.exit(1)

    output = process_file(args.fitsfile, catalog_file=args.catalog,
                          psf_path=args.psf,
                          mag_zeropoint=args.mag_zeropoint,
                          pixel_scale=args.pixel_scale,
                          max_sources=args.max_sources,
                          free_bulge_n=args.free_bulge_n,
                          n_workers=args.n_workers)
    if output:
        print(output)
    else:
        sys.exit(1)


if __name__ == '__main__':
    main()
