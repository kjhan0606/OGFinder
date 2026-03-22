#!/usr/bin/env python3
"""
DS9 LSBG (Low Surface Brightness Galaxy) Detection CLI.

Modes:
    ds9_lsbg.py FITS --mode mask [--catalog cat.tsv] [--mask-mag-threshold 22.0] ...
    ds9_lsbg.py FITS --mode clean --mask mask.fits [--bkg-method sep_large] ...
    ds9_lsbg.py FITS --mode detect --cleaned cleaned.fits [--detect-thresh 0.8] ...
    ds9_lsbg.py FITS --mode photometry --cleaned cleaned.fits --segmap seg.fits [--catalog detect.tsv] ...
    ds9_lsbg.py FITS --mode run [--n-workers 4] [--mu-eff-min 24.0] ...

Enhanced pipeline features:
    - LSB structure protection (--lsb-protect)
    - Multi-scale detection (--multiscale, --multiscale-factors)
    - Sérsic profile fitting (--sersic-fit)
    - Lower-quantile RMS background (--bkg-rms-quantile)
    - RMS convergence check (--bkg-convergence-tol)
    - Sérsic-based filtering + A/B/C/D grading
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
from lsbg.detection import detect_lsbg_candidates, detect_multiscale
from lsbg.photometry import measure_photometry
from lsbg.filtering import filter_lsbg_candidates, svm_classify_candidates


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
        max_dilate_radius=args.max_dilate_radius,
        bright_star_mag_limit=args.bright_star_mag_limit,
        bright_star_radius_scale=args.bright_star_radius_scale,
        mask_mag_threshold=args.mask_mag_threshold,
        interp_method=args.interp_method,
        lsb_protect=args.lsb_protect,
        lsb_mu_threshold=args.lsb_mu_threshold,
        bkg_method=args.bkg_method,
        bkg_mesh_size=args.bkg_mesh_size,
        bkg_poly_order=args.bkg_poly_order,
        bkg_sigma_clip=args.bkg_sigma_clip,
        bkg_n_iterations=args.bkg_n_iterations,
        bkg_refine_thresh=args.bkg_refine_thresh,
        bkg_rms_quantile=args.bkg_rms_quantile,
        bkg_convergence_tol=args.bkg_convergence_tol,
        detect_thresh=args.detect_thresh,
        detect_minarea=args.detect_minarea,
        detect_filter_kernel=args.detect_filter_kernel,
        deblend_nthresh=args.deblend_nthresh,
        deblend_mincont=args.deblend_mincont,
        multiscale=args.multiscale,
        multiscale_factors=args.multiscale_factors,
        sersic_fit=args.sersic_fit,
        sersic_n_min=args.sersic_n_min,
        sersic_n_max=args.sersic_n_max,
        sersic_re_min=args.sersic_re_min,
        sersic_cutout_scale=args.sersic_cutout_scale,
        sersic_max_nfev=args.sersic_max_nfev,
        phot_apertures=args.phot_apertures,
        mag_zeropoint=args.mag_zeropoint,
        pixel_scale=args.pixel_scale,
        mu_eff_min=args.mu_eff_min,
        mu_eff_max=args.mu_eff_max,
        r_eff_min=args.r_eff_min,
        r_eff_max=args.r_eff_max,
        ellipticity_max=args.ellipticity_max,
        min_snr=args.min_snr,
        sersic_n_filter_min=args.sersic_n_filter_min,
        sersic_n_filter_max=args.sersic_n_filter_max,
        sersic_chi2_max=args.sersic_chi2_max,
        svm_classify=args.svm_classify,
        svm_threshold=args.svm_threshold,
        svm_checkpoint=args.svm_checkpoint,
        n_workers=args.n_workers,
    )


# Column names for output TSV
COLUMNS = [
    'NUMBER', 'X_IMAGE', 'Y_IMAGE', 'RA', 'DEC',
    'A_IMAGE', 'B_IMAGE',
    'THETA_IMAGE', 'ELLIPTICITY',
    'FLUX_AUTO', 'FLUXERR_AUTO', 'MAG_AUTO', 'MAGERR_AUTO',
    'KRON_RADIUS', 'R_EFF_PIX', 'R_EFF_ARCSEC', 'MU_EFF',
    'SNR_TOTAL',
    'FLUX_APER_5', 'FLUX_APER_10', 'FLUX_APER_20', 'FLUX_APER_40',
    'SERSIC_N', 'SERSIC_RE', 'SERSIC_CHI2',
    'MAG_SERSIC', 'MU_EFF_SERSIC',
    'LSBG_CONF', 'LSBG_GRADE',
    'SVM_LSBG', 'SVM_CONF',
    'FLAGS',
]


def add_wcs_coords(catalog, header):
    """Add RA/DEC to catalog sources using FITS WCS header."""
    try:
        from astropy.wcs import WCS
        w = WCS(header)
    except Exception:
        return
    for src in catalog:
        x = src.get('x', 0)
        y = src.get('y', 0)
        try:
            coord = w.pixel_to_world(x, y)
            src['ra'] = float(coord.ra.deg)
            src['dec'] = float(coord.dec.deg)
        except Exception:
            src['ra'] = np.nan
            src['dec'] = np.nan


def format_catalog_tsv(catalog):
    """Format photometry catalog as TSV string."""
    lines = ["\t".join(COLUMNS)]
    for i, src in enumerate(catalog):
        sn = src.get('sersic_n', np.nan)
        sre = src.get('sersic_re', np.nan)
        schi2 = src.get('sersic_chi2', np.nan)
        mag_s = src.get('mag_sersic', np.nan)
        mu_s = src.get('mu_eff_sersic', np.nan)
        vals = [
            str(i + 1),
            f"{src.get('x', 0) + 1:.3f}",
            f"{src.get('y', 0) + 1:.3f}",
            f"{src.get('ra', 0):.7f}" if np.isfinite(src.get('ra', np.nan)) else "NaN",
            f"{src.get('dec', 0):.7f}" if np.isfinite(src.get('dec', np.nan)) else "NaN",
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
            f"{sn:.3f}" if np.isfinite(sn) else "NaN",
            f"{sre:.3f}" if np.isfinite(sre) else "NaN",
            f"{schi2:.3f}" if np.isfinite(schi2) else "NaN",
            f"{mag_s:.4f}" if np.isfinite(mag_s) else "NaN",
            f"{mu_s:.4f}" if np.isfinite(mu_s) else "NaN",
            f"{src.get('lsbg_conf', 0):.3f}",
            src.get('lsbg_grade', 'D'),
            str(int(src.get('svm_lsbg', -1))),
            f"{src.get('svm_conf', 0):.4f}" if np.isfinite(src.get('svm_conf', np.nan)) else "NaN",
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
    mask, masked_data = create_lsbg_mask(data, config,
                                         n_workers=args.n_workers)

    n_masked = int(mask.sum())

    mask_path = args.mask_output
    save_fits(mask.astype(np.int32), header, mask_path)

    masked_path = args.masked_output
    save_fits(masked_data, header, masked_path)

    print(f"#LSBG_MASK\tN_MASKED={n_masked}\tMASK={mask_path}\tMASKED={masked_path}")


# ---- Mode: import-mask ----

def mode_import_mask(args):
    """Import external mask FITS with WCS reprojection."""
    from astropy.io import fits as pyfits
    from icl.masking import reproject_mask, interpolate_masked

    data, header = load_fits_data(args.fits)
    ny, nx = data.shape
    print(f"Target image: {nx}x{ny}", file=sys.stderr)

    if not args.import_mask_file:
        print("ERROR: --import-mask-file required for import-mask mode",
              file=sys.stderr)
        sys.exit(1)

    with pyfits.open(args.import_mask_file) as hdul:
        mask_data = hdul[0].data
        mask_header = hdul[0].header
    print(f"External mask: {mask_data.shape[1]}x{mask_data.shape[0]}",
          file=sys.stderr)

    mask = reproject_mask(mask_data, mask_header, header)

    n_masked = int(mask.sum())
    print(f"Masked {n_masked} pixels ({100*n_masked/(ny*nx):.1f}%)",
          file=sys.stderr)

    mask_path = args.mask_output
    save_fits(mask.astype(np.int32), header, mask_path)

    print(f"Interpolating masked regions ({args.interp_method})...",
          file=sys.stderr)
    masked_data = interpolate_masked(data, mask, method=args.interp_method)

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

    if config.multiscale:
        objects, segmap, rms = detect_multiscale(
            cleaned,
            detect_thresh=config.detect_thresh,
            detect_minarea=config.detect_minarea,
            filter_kernel=config.detect_filter_kernel,
            deblend_nthresh=config.deblend_nthresh,
            deblend_mincont=config.deblend_mincont,
            scale_factors=config.multiscale_factors,
        )
    else:
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
    if segmap is not None:
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

def _detect_on_cleaned(cleaned, config):
    """Re-detect sources on cleaned image and return objects + rms."""
    try:
        import sep
    except ImportError:
        import sep_pjw as sep
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
    return objects, rms


def mode_photometry(args):
    """Measure photometry on detected LSBG candidates (Step 4)."""
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
    objects, rms = _detect_on_cleaned(cleaned, config)

    if len(objects) == 0:
        print("No sources detected for photometry", file=sys.stderr)
        print("\t".join(COLUMNS))
        return

    # Measure photometry only (no Sérsic, no filtering)
    catalog = measure_photometry(cleaned, rms, objects, segmap, config)
    print(f"Measured {len(catalog)} sources", file=sys.stderr)

    # Output unfiltered catalog with NaN Sérsic + neutral grade
    for src in catalog:
        src['sersic_n'] = np.nan
        src['sersic_re'] = np.nan
        src['sersic_chi2'] = np.nan
        src['mag_sersic'] = np.nan
        src['mu_eff_sersic'] = np.nan
        src['lsbg_conf'] = 0.0
        src['lsbg_grade'] = '-'
    tsv = format_catalog_tsv(catalog)
    print(tsv)


# ---- Mode: sersic ----

def mode_sersic(args):
    """Sérsic profile fitting on detected sources (Step 5)."""
    if args.cleaned:
        cleaned, _ = load_fits_data(args.cleaned)
    else:
        print("ERROR: --cleaned required for sersic mode", file=sys.stderr)
        sys.exit(1)

    config = build_config(args)
    objects, rms = _detect_on_cleaned(cleaned, config)

    if len(objects) == 0:
        print("No sources for Sérsic fitting", file=sys.stderr)
        print("\t".join(COLUMNS))
        return

    # Sérsic fitting
    from lsbg.sersic import fit_sources_sersic
    sersic_results = fit_sources_sersic(cleaned, objects, config)

    # Read photometry catalog from stdin or file
    if args.segmap:
        from astropy.io import fits as pyfits
        with pyfits.open(args.segmap) as hdul:
            segmap = hdul[0].data.astype(np.int32)
    else:
        print("ERROR: --segmap required for sersic mode", file=sys.stderr)
        sys.exit(1)

    data, _ = load_fits_data(args.fits)
    catalog = measure_photometry(cleaned, rms, objects, segmap, config)

    # Add Sérsic results to catalog
    zp = config.mag_zeropoint
    ps = config.pixel_scale
    for i, src in enumerate(catalog):
        sersic = None
        if sersic_results is not None and i < len(sersic_results):
            sersic = sersic_results[i]
        if sersic is not None:
            src['sersic_n'] = sersic.get('n', np.nan)
            src['sersic_re'] = sersic.get('re', np.nan)
            src['sersic_chi2'] = sersic.get('chi2', np.nan)
            flux_total = sersic.get('flux_total', np.nan)
            if np.isfinite(flux_total) and flux_total > 0:
                src['mag_sersic'] = -2.5 * np.log10(flux_total) + zp
                re_arcsec = sersic.get('re', 0) * ps
                ellip_s = sersic.get('ellip', 0)
                area = np.pi * re_arcsec**2 * (1.0 - ellip_s)
                if area > 0:
                    src['mu_eff_sersic'] = src['mag_sersic'] + 2.5 * np.log10(2.0 * area)
                else:
                    src['mu_eff_sersic'] = np.nan
            else:
                src['mag_sersic'] = np.nan
                src['mu_eff_sersic'] = np.nan
        else:
            src['sersic_n'] = np.nan
            src['sersic_re'] = np.nan
            src['sersic_chi2'] = np.nan
            src['mag_sersic'] = np.nan
            src['mu_eff_sersic'] = np.nan
        src['lsbg_conf'] = 0.0
        src['lsbg_grade'] = '-'

    tsv = format_catalog_tsv(catalog)
    print(tsv)


# ---- Mode: filter ----

def mode_filter(args):
    """Filter + grade LSBG candidates (Step 6)."""
    if args.cleaned:
        cleaned, _ = load_fits_data(args.cleaned)
    else:
        print("ERROR: --cleaned required for filter mode", file=sys.stderr)
        sys.exit(1)

    if args.segmap:
        from astropy.io import fits as pyfits
        with pyfits.open(args.segmap) as hdul:
            segmap = hdul[0].data.astype(np.int32)
    else:
        print("ERROR: --segmap required for filter mode", file=sys.stderr)
        sys.exit(1)

    config = build_config(args)
    objects, rms = _detect_on_cleaned(cleaned, config)

    if len(objects) == 0:
        print("No sources for filtering", file=sys.stderr)
        print("\t".join(COLUMNS))
        return

    # Photometry
    catalog = measure_photometry(cleaned, rms, objects, segmap, config)

    # Sérsic fitting (if enabled)
    sersic_results = None
    if config.sersic_fit:
        from lsbg.sersic import fit_sources_sersic
        sersic_results = fit_sources_sersic(cleaned, objects, config)

    # Filter + grade
    filtered, rejected = filter_lsbg_candidates(catalog, config,
                                                sersic_results=sersic_results)

    grades = {}
    for src in filtered:
        g = src.get('lsbg_grade', 'D')
        grades[g] = grades.get(g, 0) + 1
    grade_str = ', '.join(f"{g}:{n}" for g, n in sorted(grades.items()))

    print(f"LSBG candidates: {len(filtered)} passed, {len(rejected)} rejected",
          file=sys.stderr)
    if grade_str:
        print(f"  Grades: {grade_str}", file=sys.stderr)

    tsv = format_catalog_tsv(filtered)
    print(tsv)


# ---- Mode: svm-classify ----

def mode_svm_classify(args):
    """Apply SVM classification to an existing LSBG catalog."""
    if not args.catalog:
        print("ERROR: --catalog required for svm-classify mode", file=sys.stderr)
        sys.exit(1)

    # Load catalog TSV
    with open(args.catalog, 'r') as f:
        lines = f.read().strip().split('\n')
    if len(lines) < 2:
        print("ERROR: Catalog has no data rows", file=sys.stderr)
        sys.exit(1)

    headers = [h.strip() for h in lines[0].split('\t')]
    rows = []
    for line in lines[1:]:
        fields = line.split('\t')
        if len(fields) < len(headers):
            continue
        row = {}
        for h, v in zip(headers, fields):
            row[h] = v.strip()
        rows.append(row)

    if len(rows) == 0:
        print("NUMBER\tSVM_LSBG\tSVM_CONF")
        return

    print(f"SVM classifying {len(rows)} candidates...", file=sys.stderr)

    from lsbg.classify import classify
    checkpoint = args.svm_checkpoint if args.svm_checkpoint else None
    results = classify(rows, checkpoint_path=checkpoint,
                       threshold=args.svm_threshold)

    n_lsbg = sum(1 for r in results if r['svm_lsbg'] == 1)
    print(f"SVM: {n_lsbg} LSBG, {len(results) - n_lsbg} rejected",
          file=sys.stderr)

    # Output: NUMBER, SVM_LSBG, SVM_CONF (for AddColumnsFromTSV)
    print("NUMBER\tSVM_LSBG\tSVM_CONF")
    for res in results:
        print(f"{res['number']}\t{res['svm_lsbg']}\t{res['svm_conf']:.4f}")


# ---- Mode: run ----

def mode_run(args):
    """Full LSBG detection pipeline."""
    data, header = load_fits_data(args.fits)
    ny, nx = data.shape
    print(f"Image: {nx}x{ny}", file=sys.stderr)

    config = build_config(args)

    # Step 1: Mask bright sources (with LSB protection)
    print("=== Step 1: Masking bright sources ===", file=sys.stderr)
    mask, masked_data = create_lsbg_mask(data, config,
                                         n_workers=args.n_workers)

    mask_path = args.mask_output
    save_fits(mask.astype(np.int32), header, mask_path)
    masked_path = args.masked_output
    save_fits(masked_data, header, masked_path)

    # Step 2: Iterative background refinement (with convergence check)
    print("=== Step 2: Iterative background ===", file=sys.stderr)
    background, final_mask, cleaned = iterative_background(data, mask, config)

    bkg_path = args.bkg_output
    save_fits(background, header, bkg_path)
    cleaned_path = args.cleaned_output
    save_fits(cleaned, header, cleaned_path)

    # Step 3: Detection (multi-scale or single-scale)
    if config.multiscale:
        print("=== Step 3: Multi-scale detection ===", file=sys.stderr)
        # Use lower-quantile RMS for more sensitive detection
        from lsbg.background import compute_quantile_rms
        rms_q = compute_quantile_rms(cleaned, final_mask,
                                     quantile=config.bkg_rms_quantile)
        objects, segmap, rms = detect_multiscale(
            cleaned, bkg_rms=rms_q,
            detect_thresh=config.detect_thresh,
            detect_minarea=config.detect_minarea,
            filter_kernel=config.detect_filter_kernel,
            deblend_nthresh=config.deblend_nthresh,
            deblend_mincont=config.deblend_mincont,
            scale_factors=config.multiscale_factors,
        )
    else:
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
    if segmap is not None:
        save_fits(segmap.astype(np.int32), header, segmap_path)

    if len(objects) == 0:
        print("No sources detected", file=sys.stderr)
        print("\t".join(COLUMNS))
        return

    # Step 4: Photometry
    print("=== Step 4: Photometry ===", file=sys.stderr)
    catalog = measure_photometry(cleaned, rms, objects, segmap, config)

    # Add WCS coordinates (RA/DEC) for cross-matching
    add_wcs_coords(catalog, header)

    # Step 5: Sérsic profile fitting
    sersic_results = None
    if config.sersic_fit:
        print("=== Step 5: Sérsic profile fitting ===", file=sys.stderr)
        from lsbg.sersic import fit_sources_sersic
        sersic_results = fit_sources_sersic(cleaned, objects, config)

    # Step 6: Filtering + confidence grading
    print("=== Step 6: LSBG filtering + grading ===", file=sys.stderr)
    filtered, rejected = filter_lsbg_candidates(catalog, config,
                                                sersic_results=sersic_results)

    # Count grades
    grades = {}
    for src in filtered:
        g = src.get('lsbg_grade', 'D')
        grades[g] = grades.get(g, 0) + 1
    grade_str = ', '.join(f"{g}:{n}" for g, n in sorted(grades.items()))

    print(f"LSBG candidates: {len(filtered)} passed, {len(rejected)} rejected",
          file=sys.stderr)
    if grade_str:
        print(f"  Grades: {grade_str}", file=sys.stderr)

    # Step 7: SVM classification (optional)
    if config.svm_classify:
        print("=== Step 7: SVM classification ===", file=sys.stderr)
        classified, svm_ok, svm_rej = svm_classify_candidates(filtered, config)
        print(f"SVM: {len(svm_ok)} LSBG, {len(svm_rej)} rejected",
              file=sys.stderr)
        filtered = svm_ok

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


# ---- Mode: forced ----

def mode_forced(args):
    """Forced photometry at LSBG positions on another band image."""
    data, header = load_fits_data(args.fits)
    ny, nx = data.shape
    config = build_config(args)
    band = args.band_name if args.band_name else 'BAND'

    # Load source catalog
    if not args.catalog:
        print("ERROR: --catalog required for forced mode", file=sys.stderr)
        sys.exit(1)

    with open(args.catalog, 'r') as f:
        cat_lines = f.read().strip().split('\n')
    if len(cat_lines) < 2:
        print("ERROR: Catalog has no data rows", file=sys.stderr)
        sys.exit(1)

    # Parse catalog header and data
    hdr_fields = cat_lines[0].split('\t')
    col_idx = {name: i for i, name in enumerate(hdr_fields)}
    required = ['X_IMAGE', 'Y_IMAGE', 'A_IMAGE', 'B_IMAGE', 'THETA_IMAGE']
    for req in required:
        if req not in col_idx:
            print(f"ERROR: Catalog missing column {req}", file=sys.stderr)
            sys.exit(1)

    sources = []
    for line in cat_lines[1:]:
        fields = line.split('\t')
        if len(fields) < len(hdr_fields):
            continue
        src = {}
        for name, idx in col_idx.items():
            try:
                src[name] = float(fields[idx])
            except ValueError:
                src[name] = fields[idx]
        sources.append(src)

    if len(sources) == 0:
        print("No sources in catalog", file=sys.stderr)
        print(f"NUMBER\tFLUX_{band}\tFLUXERR_{band}\tMAG_{band}\tMAGERR_{band}")
        return

    print(f"Forced photometry on {len(sources)} sources, band={band}",
          file=sys.stderr)

    # Background estimate
    try:
        import sep
    except ImportError:
        import sep_pjw as sep
    data_c = np.ascontiguousarray(data, dtype=np.float64)
    try:
        bkg = sep.Background(data_c, bw=64, bh=64)
        data_sub = data_c - bkg.back()
        rms = bkg.rms()
    except Exception:
        med = np.median(data_c)
        data_sub = data_c - med
        rms = np.full_like(data_c, np.std(data_c))

    zp = config.mag_zeropoint

    # Output header
    out_cols = ['NUMBER',
                f'FLUX_{band}', f'FLUXERR_{band}',
                f'MAG_{band}', f'MAGERR_{band}']
    print("\t".join(out_cols))

    for i, src in enumerate(sources):
        # Convert from 1-indexed DS9 to 0-indexed
        x = src['X_IMAGE'] - 1.0
        y = src['Y_IMAGE'] - 1.0
        a = src['A_IMAGE']
        b = src['B_IMAGE']
        theta = np.radians(src['THETA_IMAGE'])

        # Kron photometry at fixed position
        kronrad = src.get('KRON_RADIUS', 2.5)
        if not np.isfinite(kronrad) or kronrad < 1.0:
            kronrad = 2.5

        try:
            flux, fluxerr, _ = sep.sum_ellipse(
                data_sub, [x], [y], [a], [b], [theta],
                kronrad * 2.5, err=rms, subpix=5
            )
            flux = float(flux[0])
            fluxerr = float(fluxerr[0])
        except Exception:
            flux = np.nan
            fluxerr = np.nan

        if np.isfinite(flux) and flux > 0:
            mag = -2.5 * np.log10(flux) + zp
            if np.isfinite(fluxerr) and fluxerr > 0:
                magerr = 2.5 / np.log(10.0) * fluxerr / flux
            else:
                magerr = 99.0
        else:
            mag = 99.0
            magerr = 99.0

        vals = [
            str(i + 1),
            f"{flux:.6e}",
            f"{fluxerr:.6e}",
            f"{mag:.4f}",
            f"{magerr:.4f}",
        ]
        print("\t".join(vals))


def main():
    parser = argparse.ArgumentParser(
        description='DS9 LSBG (Low Surface Brightness Galaxy) Detection Pipeline')
    parser.add_argument('fits', help='Input FITS file')
    parser.add_argument('--mode', required=True,
                        choices=['mask', 'import-mask', 'clean', 'detect',
                                 'photometry', 'sersic', 'filter',
                                 'svm-classify', 'run', 'forced'],
                        help='Operation mode')

    # Masking args
    parser.add_argument('--catalog', help='TSV catalog file')
    parser.add_argument('--mask-detect-thresh', type=float, default=1.5)
    parser.add_argument('--mask-detect-minarea', type=int, default=5)
    parser.add_argument('--mask-expand-factor', type=float, default=1.5)
    parser.add_argument('--max-dilate-radius', type=int, default=30)
    parser.add_argument('--bright-star-mag-limit', type=float, default=18.0)
    parser.add_argument('--bright-star-radius-scale', type=float, default=12.0)
    parser.add_argument('--mask-mag-threshold', type=float, default=22.0)
    parser.add_argument('--interp-method', default='linear',
                        choices=['linear', 'cubic', 'nearest'])
    parser.add_argument('--lsb-protect', action='store_true', default=True,
                        help='Protect low-SB sources from masking')
    parser.add_argument('--no-lsb-protect', dest='lsb_protect',
                        action='store_false',
                        help='Disable LSB structure protection')
    parser.add_argument('--lsb-mu-threshold', type=float, default=24.0,
                        help='SB threshold for LSB protection (mag/arcsec^2)')
    parser.add_argument('--mask-output',
                        default=os.path.expanduser('~/.ds9/lsbg_mask.fits'))
    parser.add_argument('--masked-output',
                        default=os.path.expanduser('~/.ds9/lsbg_masked.fits'))
    parser.add_argument('--import-mask-file',
                        help='External mask FITS to import')

    # Background / clean args
    parser.add_argument('--mask', help='Mask FITS file for background fitting')
    parser.add_argument('--bkg-method', default='sep_large',
                        choices=['sep_large', 'polynomial', 'chebyshev'])
    parser.add_argument('--bkg-mesh-size', type=int, default=256)
    parser.add_argument('--bkg-poly-order', type=int, default=3)
    parser.add_argument('--bkg-sigma-clip', type=float, default=3.0)
    parser.add_argument('--bkg-n-iterations', type=int, default=3)
    parser.add_argument('--bkg-refine-thresh', type=float, default=2.0)
    parser.add_argument('--bkg-rms-quantile', type=float, default=0.25,
                        help='Lower-quantile RMS for sensitive detection')
    parser.add_argument('--bkg-convergence-tol', type=float, default=0.01,
                        help='Fractional RMS change threshold for convergence')
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
                                 'gauss9x9', 'gauss15x15', 'gauss21x21',
                                 'gauss33x33', 'tophat5', 'tophat7', 'mexhat'])
    parser.add_argument('--deblend-nthresh', type=int, default=32)
    parser.add_argument('--deblend-mincont', type=float, default=0.005)
    parser.add_argument('--multiscale', action='store_true', default=True,
                        help='Enable multi-scale detection')
    parser.add_argument('--no-multiscale', dest='multiscale',
                        action='store_false',
                        help='Disable multi-scale detection')
    parser.add_argument('--multiscale-factors', default='1,2,4',
                        help='Comma-separated rebinning factors')
    parser.add_argument('--segmap', help='Segmentation map FITS')
    parser.add_argument('--segmap-output',
                        default=os.path.expanduser('~/.ds9/lsbg_segmap.fits'))

    # Sérsic fitting args
    parser.add_argument('--sersic-fit', action='store_true', default=True,
                        help='Enable Sérsic profile fitting')
    parser.add_argument('--no-sersic-fit', dest='sersic_fit',
                        action='store_false',
                        help='Disable Sérsic profile fitting')
    parser.add_argument('--sersic-n-min', type=float, default=0.2)
    parser.add_argument('--sersic-n-max', type=float, default=10.0)
    parser.add_argument('--sersic-re-min', type=float, default=0.5)
    parser.add_argument('--sersic-cutout-scale', type=float, default=5.0)
    parser.add_argument('--sersic-max-nfev', type=int, default=500)

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
    parser.add_argument('--sersic-n-filter-min', type=float, default=0.3,
                        help='Min Sérsic n for LSBG (reject artifacts)')
    parser.add_argument('--sersic-n-filter-max', type=float, default=6.0,
                        help='Max Sérsic n for LSBG (reject compact sources)')
    parser.add_argument('--sersic-chi2-max', type=float, default=10.0,
                        help='Max reduced chi2 for Sérsic fit quality')

    # SVM classification
    parser.add_argument('--svm-classify', action='store_true', default=False,
                        help='Apply SVM classification after filtering')
    parser.add_argument('--no-svm-classify', dest='svm_classify',
                        action='store_false',
                        help='Disable SVM classification')
    parser.add_argument('--svm-threshold', type=float, default=0.3,
                        help='SVM probability threshold for LSBG (0-1)')
    parser.add_argument('--svm-checkpoint', default='',
                        help='Path to SVM model checkpoint (.joblib)')

    # Output
    parser.add_argument('--catalog-output',
                        default=os.path.expanduser('~/.ds9/lsbg_catalog.tsv'))

    # Forced photometry
    parser.add_argument('--band-name', default='',
                        help='Band name for forced photometry column suffix')

    # Parallel
    parser.add_argument('--n-workers', type=int, default=0)

    args = parser.parse_args()

    # Ensure ~/.ds9 directory exists
    os.makedirs(os.path.expanduser('~/.ds9'), exist_ok=True)

    if args.mode == 'mask':
        mode_mask(args)
    elif args.mode == 'import-mask':
        mode_import_mask(args)
    elif args.mode == 'clean':
        mode_clean(args)
    elif args.mode == 'detect':
        mode_detect(args)
    elif args.mode == 'photometry':
        mode_photometry(args)
    elif args.mode == 'sersic':
        mode_sersic(args)
    elif args.mode == 'filter':
        mode_filter(args)
    elif args.mode == 'svm-classify':
        mode_svm_classify(args)
    elif args.mode == 'run':
        mode_run(args)
    elif args.mode == 'forced':
        mode_forced(args)


if __name__ == '__main__':
    main()
