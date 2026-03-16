#!/usr/bin/env python3
"""
DS9 LSBG (Low Surface Brightness Galaxy) Detection CLI.

Modes:
    ds9_lsbg.py FITS --mode mask [--catalog cat.tsv] [--mask-mag-threshold 22.0] ...
    ds9_lsbg.py FITS --mode clean --mask mask.fits [--bkg-method sep_large] ...
    ds9_lsbg.py FITS --mode detect --cleaned cleaned.fits [--detect-thresh 0.8] ...
    ds9_lsbg.py FITS --mode photometry --cleaned cleaned.fits --segmap seg.fits [--catalog detect.tsv] ...
    ds9_lsbg.py FITS --mode run [--n-workers 4] [--mu-eff-min 24.0] ...
"""

import sys
import os
import argparse
import numpy as np

_script_dir = os.path.dirname(os.path.abspath(__file__))
_project_root = os.path.abspath(os.path.join(_script_dir, '..', '..'))
if _project_root not in sys.path:
    sys.path.insert(0, _project_root)

from lsbg.config import LSBGConfig
from lsbg.masking import create_lsbg_mask
from lsbg.background import iterative_background
from lsbg.detection import detect_lsbg_candidates
from lsbg.photometry import measure_photometry
from lsbg.filtering import filter_lsbg_candidates


def load_fits_data(path):
    """Load FITS image data from primary or first image extension."""
    from astropy.io import fits
    with fits.open(path) as hdul:
        if hdul[0].data is not None:
            data = hdul[0].data.astype(np.float64)
            header = hdul[0].header
        else:
            for ext in hdul[1:]:
                if ext.data is not None and ext.data.ndim >= 2:
                    data = ext.data.astype(np.float64)
                    header = ext.header
                    break
            else:
                print("ERROR: No image data found in FITS file", file=sys.stderr)
                sys.exit(1)
        while data.ndim > 2:
            data = data[0]
        return data, header


def save_fits(data, header, path):
    """Save data as FITS with WCS header."""
    from astropy.io import fits
    hdu = fits.PrimaryHDU(data.astype(np.float32))
    wcs_keys = ['CRPIX1', 'CRPIX2', 'CRVAL1', 'CRVAL2',
                'CD1_1', 'CD1_2', 'CD2_1', 'CD2_2',
                'CDELT1', 'CDELT2', 'CTYPE1', 'CTYPE2',
                'CUNIT1', 'CUNIT2', 'EQUINOX', 'RADESYS',
                'NAXIS1', 'NAXIS2']
    for key in wcs_keys:
        if key in header:
            hdu.header[key] = header[key]
    hdu.writeto(path, overwrite=True)
    print(f"Saved: {path}", file=sys.stderr)


def build_config(args):
    """Build LSBGConfig from parsed arguments."""
    return LSBGConfig(
        mask_detect_thresh=args.mask_detect_thresh,
        mask_detect_minarea=args.mask_detect_minarea,
        mask_expand_factor=args.mask_expand_factor,
        bright_star_mag_limit=args.bright_star_mag_limit,
        bright_star_radius_scale=args.bright_star_radius_scale,
        mask_mag_threshold=args.mask_mag_threshold,
        interp_method=args.interp_method,
        bkg_method=args.bkg_method,
        bkg_mesh_size=args.bkg_mesh_size,
        bkg_poly_order=args.bkg_poly_order,
        bkg_sigma_clip=args.bkg_sigma_clip,
        bkg_n_iterations=args.bkg_n_iterations,
        bkg_refine_thresh=args.bkg_refine_thresh,
        detect_thresh=args.detect_thresh,
        detect_minarea=args.detect_minarea,
        detect_filter_kernel=args.detect_filter_kernel,
        deblend_nthresh=args.deblend_nthresh,
        deblend_mincont=args.deblend_mincont,
        phot_apertures=args.phot_apertures,
        mag_zeropoint=args.mag_zeropoint,
        pixel_scale=args.pixel_scale,
        mu_eff_min=args.mu_eff_min,
        mu_eff_max=args.mu_eff_max,
        r_eff_min=args.r_eff_min,
        r_eff_max=args.r_eff_max,
        ellipticity_max=args.ellipticity_max,
        min_snr=args.min_snr,
        n_workers=args.n_workers,
    )


# Column names for output TSV
COLUMNS = [
    'NUMBER', 'X_IMAGE', 'Y_IMAGE', 'A_IMAGE', 'B_IMAGE',
    'THETA_IMAGE', 'ELLIPTICITY',
    'FLUX_AUTO', 'FLUXERR_AUTO', 'MAG_AUTO', 'MAGERR_AUTO',
    'KRON_RADIUS', 'R_EFF_PIX', 'R_EFF_ARCSEC', 'MU_EFF',
    'SNR_TOTAL',
    'FLUX_APER_5', 'FLUX_APER_10', 'FLUX_APER_20', 'FLUX_APER_40',
    'LSBG_CONF', 'FLAGS',
]


def format_catalog_tsv(catalog):
    """Format photometry catalog as TSV string."""
    lines = ["\t".join(COLUMNS)]
    for i, src in enumerate(catalog):
        vals = [
            str(i + 1),
            f"{src.get('x', 0) + 1:.3f}",
            f"{src.get('y', 0) + 1:.3f}",
            f"{src.get('a', 0):.3f}",
            f"{src.get('b', 0):.3f}",
            f"{np.degrees(src.get('theta', 0)):.3f}",
            f"{src.get('ellipticity', 0):.4f}",
            f"{src.get('flux_auto', 0):.6e}",
            f"{src.get('fluxerr_auto', 0):.6e}",
            f"{src.get('mag_auto', 99):.4f}",
            f"{src.get('magerr_auto', 99):.4f}",
            f"{src.get('kron_radius', 0):.3f}",
            f"{src.get('r_eff_pix', 0):.3f}" if np.isfinite(src.get('r_eff_pix', np.nan)) else "NaN",
            f"{src.get('r_eff_arcsec', 0):.4f}" if np.isfinite(src.get('r_eff_arcsec', np.nan)) else "NaN",
            f"{src.get('mu_eff', 0):.4f}" if np.isfinite(src.get('mu_eff', np.nan)) else "NaN",
            f"{src.get('snr_total', 0):.3f}",
            f"{src.get('flux_aper_5', 0):.6e}",
            f"{src.get('flux_aper_10', 0):.6e}",
            f"{src.get('flux_aper_20', 0):.6e}",
            f"{src.get('flux_aper_40', 0):.6e}",
            f"{src.get('lsbg_conf', 0):.3f}",
            str(int(src.get('flags', 0))),
        ]
        lines.append("\t".join(vals))
    return "\n".join(lines)


# ---- Mode: mask ----

def mode_mask(args):
    """Create magnitude-based source mask."""
    data, header = load_fits_data(args.fits)
    ny, nx = data.shape
    print(f"Image: {nx}x{ny}", file=sys.stderr)

    config = build_config(args)
    mask, masked_data = create_lsbg_mask(data, config)

    n_masked = int(mask.sum())

    mask_path = args.mask_output
    save_fits(mask.astype(np.int32), header, mask_path)

    masked_path = args.masked_output
    save_fits(masked_data, header, masked_path)

    print(f"#LSBG_MASK\tN_MASKED={n_masked}\tMASK={mask_path}\tMASKED={masked_path}")


# ---- Mode: clean ----

def mode_clean(args):
    """Iterative background refinement."""
    data, header = load_fits_data(args.fits)
    ny, nx = data.shape

    # Load mask
    if args.mask:
        from astropy.io import fits as pyfits
        with pyfits.open(args.mask) as hdul:
            mask = hdul[0].data.astype(bool)
        print(f"Loaded mask: {int(mask.sum())} masked pixels", file=sys.stderr)
    else:
        print("ERROR: --mask required for clean mode", file=sys.stderr)
        sys.exit(1)

    config = build_config(args)

    background, final_mask, cleaned = iterative_background(data, mask, config)

    bkg_path = args.bkg_output
    save_fits(background, header, bkg_path)

    cleaned_path = args.cleaned_output
    save_fits(cleaned, header, cleaned_path)

    method = config.bkg_method
    n_iter = config.bkg_n_iterations
    print(f"#LSBG_CLEAN\tMETHOD={method}\tITER={n_iter}\t"
          f"BKG={bkg_path}\tCLEANED={cleaned_path}")


# ---- Mode: detect ----

def mode_detect(args):
    """Low-threshold source detection on cleaned image."""
    if args.cleaned:
        cleaned, header = load_fits_data(args.cleaned)
    else:
        print("ERROR: --cleaned required for detect mode", file=sys.stderr)
        sys.exit(1)

    config = build_config(args)

    objects, segmap, rms = detect_lsbg_candidates(
        cleaned,
        detect_thresh=config.detect_thresh,
        detect_minarea=config.detect_minarea,
        filter_kernel=config.detect_filter_kernel,
        deblend_nthresh=config.deblend_nthresh,
        deblend_mincont=config.deblend_mincont,
    )

    # Save segmap
    segmap_path = args.segmap_output
    save_fits(segmap.astype(np.int32), header, segmap_path)

    # Output basic detection catalog as TSV
    cols = ['NUMBER', 'X_IMAGE', 'Y_IMAGE', 'A_IMAGE', 'B_IMAGE',
            'THETA_IMAGE', 'ELLIPTICITY', 'FLUX_AUTO', 'FLAGS']
    print("\t".join(cols))
    for i, obj in enumerate(objects):
        a = float(obj['a'])
        b = float(obj['b'])
        ellip = 1.0 - b / a if a > 0 else 0.0
        vals = [
            str(i + 1),
            f"{obj['x'] + 1:.3f}",
            f"{obj['y'] + 1:.3f}",
            f"{obj['a']:.3f}",
            f"{obj['b']:.3f}",
            f"{np.degrees(obj['theta']):.3f}",
            f"{ellip:.4f}",
            f"{obj['flux']:.6e}",
            str(int(obj['flag'])),
        ]
        print("\t".join(vals))


# ---- Mode: photometry ----

def mode_photometry(args):
    """Measure photometry and filter LSBG candidates."""
    # Load original data for photometry
    data, header = load_fits_data(args.fits)

    if args.cleaned:
        cleaned, _ = load_fits_data(args.cleaned)
    else:
        print("ERROR: --cleaned required for photometry mode", file=sys.stderr)
        sys.exit(1)

    if args.segmap:
        from astropy.io import fits as pyfits
        with pyfits.open(args.segmap) as hdul:
            segmap = hdul[0].data.astype(np.int32)
    else:
        print("ERROR: --segmap required for photometry mode", file=sys.stderr)
        sys.exit(1)

    config = build_config(args)

    # Re-detect on cleaned image to get objects array
    import sep
    cleaned_c = np.ascontiguousarray(cleaned, dtype=np.float64)
    try:
        bkg = sep.Background(cleaned_c)
        rms = bkg.rms()
    except Exception:
        rms = np.full_like(cleaned_c, np.std(cleaned_c))

    objects, _, _ = detect_lsbg_candidates(
        cleaned,
        bkg_rms=rms,
        detect_thresh=config.detect_thresh,
        detect_minarea=config.detect_minarea,
        filter_kernel=config.detect_filter_kernel,
        deblend_nthresh=config.deblend_nthresh,
        deblend_mincont=config.deblend_mincont,
    )

    if len(objects) == 0:
        print("No sources detected for photometry", file=sys.stderr)
        print("\t".join(COLUMNS))
        return

    # Measure photometry on cleaned image
    catalog = measure_photometry(cleaned, rms, objects, segmap, config)
    print(f"Measured {len(catalog)} sources", file=sys.stderr)

    # Filter LSBG candidates
    filtered, rejected = filter_lsbg_candidates(catalog, config)
    print(f"LSBG candidates: {len(filtered)} passed, {len(rejected)} rejected",
          file=sys.stderr)

    # Output filtered catalog
    tsv = format_catalog_tsv(filtered)
    print(tsv)


# ---- Mode: run ----

def mode_run(args):
    """Full LSBG detection pipeline."""
    data, header = load_fits_data(args.fits)
    ny, nx = data.shape
    print(f"Image: {nx}x{ny}", file=sys.stderr)

    config = build_config(args)

    # Step 1: Mask bright sources
    print("=== Step 1: Masking bright sources ===", file=sys.stderr)
    mask, masked_data = create_lsbg_mask(data, config)

    mask_path = args.mask_output
    save_fits(mask.astype(np.int32), header, mask_path)
    masked_path = args.masked_output
    save_fits(masked_data, header, masked_path)

    # Step 2: Iterative background refinement
    print("=== Step 2: Iterative background ===", file=sys.stderr)
    background, final_mask, cleaned = iterative_background(data, mask, config)

    bkg_path = args.bkg_output
    save_fits(background, header, bkg_path)
    cleaned_path = args.cleaned_output
    save_fits(cleaned, header, cleaned_path)

    # Step 3: Low-threshold detection
    print("=== Step 3: Low-threshold detection ===", file=sys.stderr)
    objects, segmap, rms = detect_lsbg_candidates(
        cleaned,
        detect_thresh=config.detect_thresh,
        detect_minarea=config.detect_minarea,
        filter_kernel=config.detect_filter_kernel,
        deblend_nthresh=config.deblend_nthresh,
        deblend_mincont=config.deblend_mincont,
    )

    segmap_path = args.segmap_output
    save_fits(segmap.astype(np.int32), header, segmap_path)

    if len(objects) == 0:
        print("No sources detected", file=sys.stderr)
        print("\t".join(COLUMNS))
        return

    # Step 4: Photometry + filtering
    print("=== Step 4: Photometry + LSBG filtering ===", file=sys.stderr)
    catalog = measure_photometry(cleaned, rms, objects, segmap, config)
    filtered, rejected = filter_lsbg_candidates(catalog, config)

    print(f"LSBG candidates: {len(filtered)} passed, {len(rejected)} rejected",
          file=sys.stderr)

    # Save full catalog
    catalog_path = args.catalog_output
    if catalog_path:
        tsv = format_catalog_tsv(filtered)
        with open(catalog_path, 'w') as f:
            f.write(tsv + "\n")
        print(f"Catalog saved: {catalog_path}", file=sys.stderr)

    # Output to stdout
    tsv = format_catalog_tsv(filtered)
    print(tsv)


def main():
    parser = argparse.ArgumentParser(
        description='DS9 LSBG (Low Surface Brightness Galaxy) Detection Pipeline')
    parser.add_argument('fits', help='Input FITS file')
    parser.add_argument('--mode', required=True,
                        choices=['mask', 'clean', 'detect', 'photometry', 'run'],
                        help='Operation mode')

    # Masking args
    parser.add_argument('--catalog', help='TSV catalog file')
    parser.add_argument('--mask-detect-thresh', type=float, default=1.5)
    parser.add_argument('--mask-detect-minarea', type=int, default=5)
    parser.add_argument('--mask-expand-factor', type=float, default=3.0)
    parser.add_argument('--bright-star-mag-limit', type=float, default=18.0)
    parser.add_argument('--bright-star-radius-scale', type=float, default=12.0)
    parser.add_argument('--mask-mag-threshold', type=float, default=22.0)
    parser.add_argument('--interp-method', default='linear',
                        choices=['linear', 'cubic', 'nearest'])
    parser.add_argument('--mask-output',
                        default=os.path.expanduser('~/.ds9/lsbg_mask.fits'))
    parser.add_argument('--masked-output',
                        default=os.path.expanduser('~/.ds9/lsbg_masked.fits'))

    # Background / clean args
    parser.add_argument('--mask', help='Mask FITS file for background fitting')
    parser.add_argument('--bkg-method', default='sep_large',
                        choices=['sep_large', 'polynomial', 'chebyshev'])
    parser.add_argument('--bkg-mesh-size', type=int, default=256)
    parser.add_argument('--bkg-poly-order', type=int, default=3)
    parser.add_argument('--bkg-sigma-clip', type=float, default=3.0)
    parser.add_argument('--bkg-n-iterations', type=int, default=3)
    parser.add_argument('--bkg-refine-thresh', type=float, default=2.0)
    parser.add_argument('--bkg-output',
                        default=os.path.expanduser('~/.ds9/lsbg_background.fits'))
    parser.add_argument('--cleaned-output',
                        default=os.path.expanduser('~/.ds9/lsbg_cleaned.fits'))

    # Detection args
    parser.add_argument('--cleaned', help='Cleaned (background-subtracted) FITS')
    parser.add_argument('--detect-thresh', type=float, default=0.8)
    parser.add_argument('--detect-minarea', type=int, default=50)
    parser.add_argument('--detect-filter-kernel', default='gauss5x5',
                        choices=['none', 'gauss3x3', 'gauss5x5', 'gauss7x7',
                                 'gauss9x9', 'tophat5', 'tophat7', 'mexhat'])
    parser.add_argument('--deblend-nthresh', type=int, default=32)
    parser.add_argument('--deblend-mincont', type=float, default=0.005)
    parser.add_argument('--segmap', help='Segmentation map FITS')
    parser.add_argument('--segmap-output',
                        default=os.path.expanduser('~/.ds9/lsbg_segmap.fits'))

    # Photometry args
    parser.add_argument('--phot-apertures', default='5,10,20,40')
    parser.add_argument('--mag-zeropoint', type=float, default=25.0)
    parser.add_argument('--pixel-scale', type=float, default=0.06)

    # Filtering args
    parser.add_argument('--mu-eff-min', type=float, default=24.0)
    parser.add_argument('--mu-eff-max', type=float, default=30.0)
    parser.add_argument('--r-eff-min', type=float, default=2.5)
    parser.add_argument('--r-eff-max', type=float, default=60.0)
    parser.add_argument('--ellipticity-max', type=float, default=0.7)
    parser.add_argument('--min-snr', type=float, default=2.0)

    # Output
    parser.add_argument('--catalog-output',
                        default=os.path.expanduser('~/.ds9/lsbg_catalog.tsv'))

    # Parallel
    parser.add_argument('--n-workers', type=int, default=0)

    args = parser.parse_args()

    if args.mode == 'mask':
        mode_mask(args)
    elif args.mode == 'clean':
        mode_clean(args)
    elif args.mode == 'detect':
        mode_detect(args)
    elif args.mode == 'photometry':
        mode_photometry(args)
    elif args.mode == 'run':
        mode_run(args)


if __name__ == '__main__':
    main()
