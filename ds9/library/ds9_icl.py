#!/usr/bin/env python3
"""
DS9 ICL (Intra-Cluster Light) Detection CLI.

Modes:
    ds9_icl.py FITS --mode mask [--catalog cat.tsv] [--expand-factor 2.5] ...
    ds9_icl.py FITS --mode background --mask mask.fits [--bkg-method polynomial] ...
    ds9_icl.py FITS --mode profile --center 500,400 [--rmax 1000] ...
    ds9_icl.py FITS --mode measure --profile-file profile.tsv [--mu-threshold 26.5] ...
    ds9_icl.py FITS --mode measure-multi --profile-file profile.tsv [--mu-min 25] ...
    ds9_icl.py FITS --mode decompose --profile-file profile.tsv [--pixel-scale 0.06] ...
    ds9_icl.py FITS --mode color --mask mask.fits --center 500,400 --bands F606W:f606w.fits,...
"""

import sys
import os
import argparse
import numpy as np

_script_dir = os.path.dirname(os.path.abspath(__file__))
_project_root = os.path.abspath(os.path.join(_script_dir, '..', '..'))
if _project_root not in sys.path:
    sys.path.insert(0, _project_root)

from icl.config import ICLConfig
from icl.masking import create_source_mask, mask_bright_stars, interpolate_masked
from icl.background import (fit_polynomial_background,
                             fit_chebyshev_background,
                             sep_large_mesh_background)
from icl.profile import find_bcg, measure_sb_profile, measure_sector_profile
from icl.measure import compute_icl_fraction, compute_isophotal_radii, summarize_icl


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


def parse_catalog_tsv(path):
    """Parse TSV catalog from ds9_sextract panel."""
    with open(path, 'r') as f:
        lines = f.readlines()

    if not lines:
        return None

    header = lines[0].strip().split('\t')
    col_idx = {h.strip(): i for i, h in enumerate(header)}

    col_map = {
        'NUMBER': 'number',
        'X_IMAGE': 'x_image',
        'Y_IMAGE': 'y_image',
        'A_IMAGE': 'a_image',
        'B_IMAGE': 'b_image',
        'THETA_IMAGE': 'theta_image',
        'FLUX_AUTO': 'flux_auto',
        'MAG_AUTO': 'mag_auto',
        'KRON_RADIUS': 'kron_radius',
        'FLUX_RADIUS': 'flux_radius',
        'ELLIPTICITY': 'ellipticity',
        'CLASS_STAR': 'class_star',
        'FLAGS': 'flags',
    }

    available = {}
    for tsv_col, key in col_map.items():
        if tsv_col in col_idx:
            available[key] = col_idx[tsv_col]

    if 'number' not in available or 'x_image' not in available:
        print("ERROR: catalog missing NUMBER or X_IMAGE columns", file=sys.stderr)
        return None

    data = {key: [] for key in available}
    for line in lines[1:]:
        line = line.strip()
        if not line:
            continue
        fields = line.split('\t')
        try:
            for key, idx in available.items():
                data[key].append(float(fields[idx]))
        except (IndexError, ValueError):
            continue

    if not data['number']:
        return None

    result = {}
    for key, vals in data.items():
        if key in ('number', 'flags'):
            result[key] = np.array(vals, dtype=int)
        else:
            result[key] = np.array(vals, dtype=float)

    return result


# ---- Mode: mask ----

def mode_mask(args):
    """Create source mask and interpolated image."""
    try:
        import sep
    except ImportError:
        import sep_pjw as sep

    data, header = load_fits_data(args.fits)
    ny, nx = data.shape
    print(f"Image: {nx}x{ny}", file=sys.stderr)

    config = ICLConfig(
        mask_expand_factor=args.expand_factor,
        bright_star_mag_limit=args.bright_star_mag_limit,
        bright_star_radius_scale=args.bright_star_radius_scale,
        interp_method=args.interp_method,
    )

    # Run SEP for segmentation map
    data_c = np.ascontiguousarray(data, dtype=np.float64)
    bkg = sep.Background(data_c)
    data_sub = data_c - bkg

    thresh = args.detect_thresh * bkg.globalrms
    objects, segmap = sep.extract(data_sub, thresh, segmentation_map=True)
    n_sources = len(objects)
    print(f"Detected {n_sources} sources for masking", file=sys.stderr)

    # Create expanded source mask
    mask = create_source_mask(data, segmap, config.mask_expand_factor,
                              max_dilate_radius=args.max_dilate_radius,
                              n_workers=args.n_workers)

    # Add bright star masks if catalog provided
    if args.catalog:
        catalog = parse_catalog_tsv(args.catalog)
        if catalog is not None:
            mask_bright_stars(mask, catalog,
                             config.bright_star_mag_limit,
                             config.bright_star_radius_scale)

    n_masked = int(mask.sum())
    print(f"Masked {n_masked} pixels ({100*n_masked/(ny*nx):.1f}%)", file=sys.stderr)

    # Save mask FITS
    mask_path = args.mask_output
    save_fits(mask.astype(np.int32), header, mask_path)

    # Interpolate masked regions
    print(f"Interpolating masked regions ({config.interp_method})...",
          file=sys.stderr)
    masked_data = interpolate_masked(data, mask, method=config.interp_method)

    masked_path = args.masked_output
    save_fits(masked_data, header, masked_path)

    # Output metadata
    print(f"#ICL_MASK\tN_MASKED={n_masked}\tMASK={mask_path}\tMASKED={masked_path}")


# ---- Mode: background ----

def mode_background(args):
    """Fit large-scale background model."""
    data, header = load_fits_data(args.fits)
    ny, nx = data.shape

    # Load mask
    mask = None
    if args.mask:
        from astropy.io import fits
        with fits.open(args.mask) as hdul:
            mask = hdul[0].data.astype(bool)
        print(f"Loaded mask: {int(mask.sum())} masked pixels", file=sys.stderr)
    else:
        mask = np.zeros((ny, nx), dtype=bool)

    method = args.bkg_method

    if args.iterative:
        from icl.background import iterative_background_icl
        print(f"Iterative background ({method}, {args.bkg_n_iterations} iters)...",
              file=sys.stderr)
        config = ICLConfig(
            bkg_method=method,
            bkg_poly_order=args.bkg_order,
            bkg_sigma_clip=args.bkg_sigma_clip,
            bkg_sep_mesh=args.bkg_sep_mesh,
            bkg_n_iterations=args.bkg_n_iterations,
            bkg_convergence_tol=args.bkg_convergence_tol,
            bkg_refine_thresh=args.bkg_refine_thresh,
            interp_method=args.interp_method,
        )
        background, final_mask, bgsub = iterative_background_icl(data, mask, config)

        # Save updated mask
        mask_out = args.mask_output if args.mask_output else \
            os.path.expanduser('~/.ds9/icl_mask.fits')
        save_fits(final_mask.astype(np.int32), header, mask_out)
    else:
        print(f"Fitting background ({method}, order={args.bkg_order})...",
              file=sys.stderr)

        if method == 'polynomial':
            background = fit_polynomial_background(
                data, mask, order=args.bkg_order, sigma_clip=args.bkg_sigma_clip)
        elif method == 'chebyshev':
            background = fit_chebyshev_background(
                data, mask, order=args.bkg_order, sigma_clip=args.bkg_sigma_clip)
        elif method == 'sep_large':
            background = sep_large_mesh_background(
                data, mask, mesh_size=args.bkg_sep_mesh)
        else:
            print(f"ERROR: Unknown background method: {method}", file=sys.stderr)
            sys.exit(1)
        bgsub = data - background

    # Save background model
    bkg_path = args.bkg_output
    save_fits(background, header, bkg_path)

    # Save background-subtracted image
    sub_path = args.bgsub_output
    save_fits(bgsub, header, sub_path)

    print(f"#ICL_BKG\tMETHOD={method}\tBKG={bkg_path}\tSUB={sub_path}")


# ---- Mode: profile ----

def mode_profile(args):
    """Measure surface brightness profile."""
    data, header = load_fits_data(args.fits)

    config = ICLConfig(
        profile_rmin=args.rmin,
        profile_rmax=args.rmax,
        profile_nsteps=args.nsteps,
        profile_spacing=args.spacing,
        profile_ellipticity=args.ellipticity,
        profile_pa=args.pa,
        profile_sector_width=args.sector_width,
        profile_sector_pa=args.sector_pa,
        mag_zeropoint=args.mag_zeropoint,
        pixel_scale=args.pixel_scale,
    )

    # Determine center
    if args.center:
        parts = args.center.split(',')
        cx = float(parts[0])
        cy = float(parts[1])
    elif args.catalog:
        catalog = parse_catalog_tsv(args.catalog)
        if catalog is None:
            print("ERROR: Failed to parse catalog for BCG finding",
                  file=sys.stderr)
            sys.exit(1)
        cx, cy = find_bcg(data, catalog)
        print(f"Auto-detected BCG at ({cx:.1f}, {cy:.1f})", file=sys.stderr)
    else:
        # Use image center
        ny, nx = data.shape
        cx, cy = nx / 2.0, ny / 2.0

    print(f"Measuring SB profile: center=({cx:.1f},{cy:.1f}), "
          f"rmin={config.profile_rmin}, rmax={config.profile_rmax}",
          file=sys.stderr)

    if args.sector_width < 360.0:
        profile = measure_sector_profile(
            data, cx, cy, args.sector_pa, args.sector_width, config)
    else:
        profile = measure_sb_profile(data, cx, cy, config)

    if not profile:
        print("ERROR: No valid profile points measured", file=sys.stderr)
        sys.exit(1)

    # Save profile TSV
    if args.profile_output:
        with open(args.profile_output, 'w') as f:
            f.write("R_PIX\tR_ARCSEC\tMU\tMU_ERR\tFLUX\tAREA\tNPIX\n")
            for p in profile:
                f.write(f"{p['R_PIX']:.4f}\t{p['R_ARCSEC']:.6f}\t"
                        f"{p['MU']:.4f}\t{p['MU_ERR']:.4f}\t"
                        f"{p['FLUX']:.6e}\t{p['AREA']:.2f}\t{p['NPIX']}\n")
        print(f"Profile saved to {args.profile_output}", file=sys.stderr)

    # Output TSV to stdout
    print("R_PIX\tR_ARCSEC\tMU\tMU_ERR\tFLUX\tAREA\tNPIX")
    for p in profile:
        print(f"{p['R_PIX']:.4f}\t{p['R_ARCSEC']:.6f}\t"
              f"{p['MU']:.4f}\t{p['MU_ERR']:.4f}\t"
              f"{p['FLUX']:.6e}\t{p['AREA']:.2f}\t{p['NPIX']}")


# ---- Mode: measure ----

def mode_measure(args):
    """Compute ICL fraction and isophotal radii from profile."""
    config = ICLConfig(
        icl_mu_threshold=args.mu_threshold,
        icl_mu_levels=args.mu_levels,
        pixel_scale=args.pixel_scale,
    )

    # Load profile from file or stdin
    if args.profile_file:
        with open(args.profile_file, 'r') as f:
            lines = f.readlines()
    else:
        print("ERROR: --profile-file required for measure mode",
              file=sys.stderr)
        sys.exit(1)

    # Parse profile TSV
    header = lines[0].strip().split('\t')
    col_idx = {h: i for i, h in enumerate(header)}

    profile = []
    for line in lines[1:]:
        line = line.strip()
        if not line:
            continue
        fields = line.split('\t')
        try:
            profile.append({
                'R_PIX': float(fields[col_idx['R_PIX']]),
                'R_ARCSEC': float(fields[col_idx['R_ARCSEC']]),
                'MU': float(fields[col_idx['MU']]),
                'MU_ERR': float(fields[col_idx['MU_ERR']]),
                'FLUX': float(fields[col_idx['FLUX']]),
                'AREA': float(fields[col_idx['AREA']]),
                'NPIX': int(float(fields[col_idx['NPIX']])),
            })
        except (IndexError, ValueError, KeyError):
            continue

    if not profile:
        print("ERROR: No valid profile data", file=sys.stderr)
        sys.exit(1)

    summary = summarize_icl(profile, config)

    print(f"ICL Analysis (mu_threshold={config.icl_mu_threshold}):",
          file=sys.stderr)
    print(f"  F_ICL = {summary['F_ICL']:.4f}", file=sys.stderr)
    print(f"  L_ICL = {summary['L_ICL']:.6e}", file=sys.stderr)
    print(f"  L_GAL = {summary['L_GAL']:.6e}", file=sys.stderr)

    # Output TSV
    keys = sorted(summary.keys())
    print("\t".join(keys))
    vals = []
    for k in keys:
        v = summary[k]
        if isinstance(v, float):
            if np.isnan(v):
                vals.append("NaN")
            else:
                vals.append(f"{v:.6e}")
        else:
            vals.append(str(v))
    print("\t".join(vals))


# ---- Mode: measure-multi ----

def _load_profile_tsv(path):
    """Load profile TSV file and return list of dicts."""
    with open(path, 'r') as f:
        lines = f.readlines()
    if not lines:
        return []
    header = lines[0].strip().split('\t')
    col_idx = {h: i for i, h in enumerate(header)}
    profile = []
    for line in lines[1:]:
        line = line.strip()
        if not line:
            continue
        fields = line.split('\t')
        try:
            profile.append({
                'R_PIX': float(fields[col_idx['R_PIX']]),
                'R_ARCSEC': float(fields[col_idx['R_ARCSEC']]),
                'MU': float(fields[col_idx['MU']]),
                'MU_ERR': float(fields[col_idx['MU_ERR']]),
                'FLUX': float(fields[col_idx['FLUX']]),
                'AREA': float(fields[col_idx['AREA']]),
                'NPIX': int(float(fields[col_idx['NPIX']])),
            })
        except (IndexError, ValueError, KeyError):
            continue
    return profile


def mode_measure_multi(args):
    """Compute ICL fraction at multiple SB thresholds with bootstrap errors."""
    from icl.measure import compute_icl_fraction_multi

    if not args.profile_file:
        print("ERROR: --profile-file required for measure-multi mode",
              file=sys.stderr)
        sys.exit(1)

    profile = _load_profile_tsv(args.profile_file)
    if not profile:
        print("ERROR: No valid profile data", file=sys.stderr)
        sys.exit(1)

    print(f"Multi-threshold ICL: mu={args.mu_min}-{args.mu_max}, "
          f"step={args.mu_step}, bootstrap={args.n_bootstrap}",
          file=sys.stderr)

    results = compute_icl_fraction_multi(
        profile,
        mu_min=args.mu_min,
        mu_max=args.mu_max,
        mu_step=args.mu_step,
        n_bootstrap=args.n_bootstrap
    )

    # Output TSV
    print("MU_THRESHOLD\tF_ICL\tF_ICL_ERR\tL_ICL\tL_GAL")
    for r in results:
        print(f"{r['MU_THRESHOLD']:.2f}\t{r['F_ICL']:.6f}\t"
              f"{r['F_ICL_ERR']:.6f}\t{r['L_ICL']:.6e}\t{r['L_GAL']:.6e}")


# ---- Mode: decompose ----

def mode_decompose(args):
    """Fit double Sérsic BCG+ICL decomposition to SB profile."""
    from icl.decompose import fit_double_sersic

    if not args.profile_file:
        print("ERROR: --profile-file required for decompose mode",
              file=sys.stderr)
        sys.exit(1)

    profile = _load_profile_tsv(args.profile_file)
    if not profile:
        print("ERROR: No valid profile data", file=sys.stderr)
        sys.exit(1)

    print(f"Double Sérsic BCG+ICL decomposition...", file=sys.stderr)

    result = fit_double_sersic(
        profile,
        pixel_scale=args.pixel_scale,
        mag_zeropoint=args.mag_zeropoint
    )

    if result is None:
        print("ERROR: Decomposition failed", file=sys.stderr)
        sys.exit(1)

    print(f"  BCG: n={result['n_bcg']:.2f}, re={result['re_bcg']:.1f}px "
          f"({result['re_bcg_arcsec']:.2f}\")", file=sys.stderr)
    print(f"  ICL: n={result['n_icl']:.2f}, re={result['re_icl']:.1f}px "
          f"({result['re_icl_arcsec']:.2f}\")", file=sys.stderr)
    print(f"  R_transition={result['R_transition']:.1f}px, "
          f"f_ICL={result['f_ICL_sersic']:.4f}, chi2={result['chi2']:.4f}",
          file=sys.stderr)

    # Output TSV
    keys = sorted(result.keys())
    print("\t".join(keys))
    vals = []
    for k in keys:
        v = result[k]
        if isinstance(v, float):
            if np.isnan(v):
                vals.append("NaN")
            else:
                vals.append(f"{v:.6e}")
        else:
            vals.append(str(v))
    print("\t".join(vals))


# ---- Mode: color ----

def mode_color(args):
    """Measure multi-band ICL profiles and color indices."""
    from icl.colors import measure_multiband_icl_profiles, compute_color_profiles

    # Parse band spec: "F606W:path1.fits,F814W:path2.fits"
    if not args.bands:
        print("ERROR: --bands required for color mode", file=sys.stderr)
        sys.exit(1)

    band_files = {}
    for spec in args.bands.split(','):
        parts = spec.strip().split(':')
        if len(parts) != 2:
            print(f"ERROR: Invalid band spec '{spec}' (use NAME:path)",
                  file=sys.stderr)
            sys.exit(1)
        band_files[parts[0]] = parts[1]

    if len(band_files) < 2:
        print("ERROR: Need at least 2 bands for color analysis",
              file=sys.stderr)
        sys.exit(1)

    # Load mask
    mask = None
    if args.mask:
        from astropy.io import fits
        with fits.open(args.mask) as hdul:
            mask = hdul[0].data.astype(bool)

    # Center
    if args.center:
        parts = args.center.split(',')
        cx = float(parts[0])
        cy = float(parts[1])
    else:
        print("ERROR: --center required for color mode", file=sys.stderr)
        sys.exit(1)

    config = ICLConfig(
        profile_rmin=args.rmin,
        profile_rmax=args.rmax,
        profile_nsteps=args.nsteps,
        profile_spacing=args.spacing,
        mag_zeropoint=args.mag_zeropoint,
        pixel_scale=args.pixel_scale,
    )

    print(f"Measuring {len(band_files)} bands: {list(band_files.keys())}",
          file=sys.stderr)

    profiles = measure_multiband_icl_profiles(band_files, mask, cx, cy, config)

    # Output per-band profiles + colors
    band_names = sorted(band_files.keys())

    # Header
    cols = ['R_PIX', 'R_ARCSEC']
    for b in band_names:
        cols.extend([f'MU_{b}', f'MU_ERR_{b}'])
    # Add color columns for consecutive pairs
    for i in range(len(band_names) - 1):
        cols.append(f'{band_names[i]}-{band_names[i+1]}')
    print("\t".join(cols))

    # Build radius-indexed data
    # Use first band's radii as reference
    ref_profile = profiles[band_names[0]]
    for j, p in enumerate(ref_profile):
        row = [f"{p['R_PIX']:.4f}", f"{p['R_ARCSEC']:.6f}"]
        mus = {}
        for b in band_names:
            bp = profiles[b]
            if j < len(bp):
                row.extend([f"{bp[j]['MU']:.4f}", f"{bp[j]['MU_ERR']:.4f}"])
                mus[b] = bp[j]['MU']
            else:
                row.extend(["99.0000", "99.0000"])
                mus[b] = 99.0

        # Colors
        for i in range(len(band_names) - 1):
            b1 = band_names[i]
            b2 = band_names[i + 1]
            if mus[b1] < 99.0 and mus[b2] < 99.0:
                row.append(f"{mus[b1] - mus[b2]:.4f}")
            else:
                row.append("NaN")

        print("\t".join(row))


def main():
    parser = argparse.ArgumentParser(
        description='DS9 ICL (Intra-Cluster Light) Detection Pipeline')
    parser.add_argument('fits', help='Input FITS file')
    parser.add_argument('--mode', required=True,
                        choices=['mask', 'background', 'profile',
                                 'measure', 'measure-multi',
                                 'decompose', 'color'],
                        help='Operation mode')

    # Masking args
    parser.add_argument('--catalog', help='TSV catalog file')
    parser.add_argument('--expand-factor', type=float, default=2.5,
                        help='Segmap expansion factor')
    parser.add_argument('--bright-star-mag-limit', type=float, default=18.0)
    parser.add_argument('--bright-star-radius-scale', type=float, default=10.0)
    parser.add_argument('--interp-method', default='linear',
                        choices=['linear', 'cubic', 'nearest'])
    parser.add_argument('--detect-thresh', type=float, default=1.5,
                        help='Detection threshold for segmap (sigma)')
    parser.add_argument('--max-dilate-radius', type=int, default=50,
                        help='Max dilation radius in pixels (prevents over-masking)')
    parser.add_argument('--mask-output',
                        default=os.path.expanduser('~/.ds9/icl_mask.fits'))
    parser.add_argument('--masked-output',
                        default=os.path.expanduser('~/.ds9/icl_masked.fits'))

    # Background args
    parser.add_argument('--mask', help='Mask FITS file for background fitting')
    parser.add_argument('--bkg-method', default='polynomial',
                        choices=['polynomial', 'chebyshev', 'sep_large'])
    parser.add_argument('--bkg-order', type=int, default=3)
    parser.add_argument('--bkg-sigma-clip', type=float, default=3.0)
    parser.add_argument('--bkg-sep-mesh', type=int, default=256)
    parser.add_argument('--bkg-output',
                        default=os.path.expanduser('~/.ds9/icl_background.fits'))
    parser.add_argument('--bgsub-output',
                        default=os.path.expanduser('~/.ds9/icl_bgsub.fits'))
    parser.add_argument('--iterative', action='store_true',
                        help='Use iterative background refinement')
    parser.add_argument('--bkg-n-iterations', type=int, default=3,
                        help='Max iterations for iterative background')
    parser.add_argument('--bkg-convergence-tol', type=float, default=0.01,
                        help='RMS convergence tolerance')
    parser.add_argument('--bkg-refine-thresh', type=float, default=2.0,
                        help='Sigma threshold for residual source detection')

    # Profile args
    parser.add_argument('--center', help='BCG center as x,y (0-indexed)')
    parser.add_argument('--rmin', type=float, default=5.0)
    parser.add_argument('--rmax', type=float, default=1000.0)
    parser.add_argument('--nsteps', type=int, default=80)
    parser.add_argument('--spacing', default='log',
                        choices=['log', 'linear'])
    parser.add_argument('--ellipticity', type=float, default=0.0)
    parser.add_argument('--pa', type=float, default=0.0)
    parser.add_argument('--sector-width', type=float, default=360.0)
    parser.add_argument('--sector-pa', type=float, default=0.0)
    parser.add_argument('--mag-zeropoint', type=float, default=25.0)
    parser.add_argument('--pixel-scale', type=float, default=0.06)
    parser.add_argument('--profile-output',
                        default=os.path.expanduser('~/.ds9/icl_profile.tsv'))

    # Measure args
    parser.add_argument('--profile-file', help='Profile TSV file')
    parser.add_argument('--mu-threshold', type=float, default=26.5)
    parser.add_argument('--mu-levels', default='26.0,27.0,28.0')

    # Multi-threshold args
    parser.add_argument('--mu-min', type=float, default=25.0,
                        help='Min SB threshold for multi-threshold')
    parser.add_argument('--mu-max', type=float, default=28.0,
                        help='Max SB threshold for multi-threshold')
    parser.add_argument('--mu-step', type=float, default=0.5,
                        help='SB threshold step')
    parser.add_argument('--n-bootstrap', type=int, default=200,
                        help='Bootstrap resamples for error estimation')

    # Color args
    parser.add_argument('--bands', help='Band specs: NAME:path,NAME:path,...')

    # Parallel
    parser.add_argument('--n-workers', type=int, default=0)

    args = parser.parse_args()

    if args.mode == 'mask':
        mode_mask(args)
    elif args.mode == 'background':
        mode_background(args)
    elif args.mode == 'profile':
        mode_profile(args)
    elif args.mode == 'measure':
        mode_measure(args)
    elif args.mode == 'measure-multi':
        mode_measure_multi(args)
    elif args.mode == 'decompose':
        mode_decompose(args)
    elif args.mode == 'color':
        mode_color(args)


if __name__ == '__main__':
    main()
