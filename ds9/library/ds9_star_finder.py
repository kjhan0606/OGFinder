#!/usr/bin/env python3
"""
DS9 AI Star/Galaxy Classification CLI.

Uses PSF cross-correlation: builds empirical PSF from the image itself,
then classifies each source by how well it matches the PSF.
No pre-trained model needed — image-adaptive.

Usage:
    ds9_star_finder.py FITSFILE --catalog catalog.tsv
                       [--threshold 0.6] [--cutout-size 64]
                       [--psf psf.fits]
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

from star_finder.psf_match import classify_sources, build_empirical_psf


def parse_catalog_tsv(path):
    """Parse TSV catalog. Required: NUMBER, X_IMAGE, Y_IMAGE."""
    columns_map = {
        'NUMBER': 'number', 'X_IMAGE': 'x', 'Y_IMAGE': 'y',
        'FLUX_AUTO': 'flux', 'A_IMAGE': 'a_image', 'B_IMAGE': 'b_image',
        'CLASS_STAR': 'class_star', 'FWHM_IMAGE': 'fwhm_image',
    }

    with open(path, 'r') as f:
        lines = f.readlines()
    if not lines:
        return None

    header = lines[0].strip().split('\t')
    col_idx = {}
    for i, h in enumerate(header):
        h = h.strip()
        if h in columns_map:
            col_idx[columns_map[h]] = i

    if 'number' not in col_idx or 'x' not in col_idx or 'y' not in col_idx:
        print("ERROR: Missing NUMBER/X_IMAGE/Y_IMAGE columns", file=sys.stderr)
        return None

    result = {k: [] for k in col_idx}
    for line in lines[1:]:
        line = line.strip()
        if not line:
            continue
        fields = line.split('\t')
        try:
            for key, ci in col_idx.items():
                result[key].append(float(fields[ci]))
        except (IndexError, ValueError):
            continue

    if not result['number']:
        return None

    # Convert to numpy
    for key in result:
        result[key] = np.array(result[key])

    result['number'] = result['number'].astype(np.int64)
    result['x'] -= 1.0  # 1-based -> 0-based
    result['y'] -= 1.0

    return result


def load_fits(fits_path):
    """Load FITS data."""
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


def load_psf_fits(psf_path):
    """Load PSF from FITS file."""
    from astropy.io import fits as pyfits
    with pyfits.open(psf_path) as hdul:
        for hdu in hdul:
            if hdu.data is not None and hdu.data.ndim == 2:
                return hdu.data.astype(np.float64)
    return None


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
    bkg = sep.Background(data, bw=64, bh=64, fw=3, fh=3)
    return data - bkg.back()


def main():
    parser = argparse.ArgumentParser(
        description='DS9 Star/Galaxy Classification (PSF matching)')
    parser.add_argument('fits', help='Path to FITS file')
    parser.add_argument('--catalog', required=True,
                        help='Path to catalog TSV')
    parser.add_argument('--threshold', type=float, default=0.6,
                        help='PSF correlation threshold (default: 0.6)')
    parser.add_argument('--cutout-size', type=int, default=64,
                        help='Cutout size in pixels')
    parser.add_argument('--psf', default=None,
                        help='Path to PSF FITS file (auto-detect if omitted)')
    args = parser.parse_args()

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

    # Load PSF if provided
    psf_template = None
    if args.psf and os.path.exists(args.psf):
        print(f"Loading PSF: {args.psf}", file=sys.stderr)
        psf_template = load_psf_fits(args.psf)
        if psf_template is not None:
            # Resize PSF to match cutout size if needed
            if (psf_template.shape[0] != args.cutout_size or
                    psf_template.shape[1] != args.cutout_size):
                from star_finder.cutout.extract import _resize_2d
                psf_template = _resize_2d(psf_template, args.cutout_size,
                                           args.cutout_size)
            print(f"  PSF shape: {psf_template.shape}", file=sys.stderr)

    # Classify
    print("Classifying sources by PSF correlation ...", file=sys.stderr)
    results = classify_sources(data_sub, catalog, size=args.cutout_size,
                                psf_template=psf_template,
                                star_threshold=args.threshold)

    # Count
    n_stars = sum(1 for r in results if r['label'] == 'star')
    n_galaxies = sum(1 for r in results if r['label'] == 'galaxy')
    n_classified = len(results)
    print(f"  Classified: {n_classified} (star={n_stars}, galaxy={n_galaxies})",
          file=sys.stderr)

    # Output TSV
    print(f"#STAR_FINDER\tN_CLASSIFIED={n_classified}\tN_SOURCES={n_sources}")
    print("NUMBER\tAI_STAR\tAI_STAR_CONF\tAI_STAR_COLOR")

    for r in results:
        print(f"{r['number']}\t{r['label']}\t{r['correlation']:.4f}\t"
              f"{r['color']}")


if __name__ == '__main__':
    main()
