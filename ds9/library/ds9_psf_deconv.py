#!/usr/bin/env python3
"""
DS9 PSF Deconvolution CLI — star finding, PSF building, deconvolution.

Outputs TSV (find_stars) or FITS (build_psf, deconvolve) for Tcl parsing.

Usage:
    ds9_psf_deconv.py FITS --mode find_stars --catalog cat.tsv [--method combined]
    ds9_psf_deconv.py FITS --mode build_psf --star-indices 1,5,12 [--method median]
    ds9_psf_deconv.py FITS --mode deconvolve --psf psf.fits --algorithm rl [--iterations 30]
"""

import sys
import os
import argparse
import numpy as np

_script_dir = os.path.dirname(os.path.abspath(__file__))
_project_root = os.path.abspath(os.path.join(_script_dir, '..', '..'))
if _project_root not in sys.path:
    sys.path.insert(0, _project_root)

from psf_deconv.config import PSFDeconvConfig
from psf_deconv.stars.find_stars import (find_stars_class_star,
                                          find_stars_fwhm,
                                          find_stars_combined)
from psf_deconv.psf.build_psf import (build_psf_median, build_psf_mean,
                                       build_psf_epsf, build_psf_extended)
from psf_deconv.psf.fit_psf import fit_gaussian_psf, fit_moffat_psf
from psf_deconv.deconv.richardson_lucy import (richardson_lucy,
                                                richardson_lucy_accelerated,
                                                richardson_lucy_tv)
from psf_deconv.deconv.wiener import wiener_deconvolve
from psf_deconv.deconv.tikhonov import tikhonov_deconvolve
from psf_deconv.deconv.clean import clean_deconvolve
from psf_deconv.deconv.mem import mem_deconvolve


def parse_catalog_tsv(path):
    """Parse TSV catalog from ds9_sextract panel.

    Returns dict with lowercase keys and numpy arrays.
    Coordinates remain 1-based (for display).
    """
    with open(path, 'r') as f:
        lines = f.readlines()

    if not lines:
        return None

    header = lines[0].strip().split('\t')
    col_idx = {h.strip(): i for i, h in enumerate(header)}

    # Map column names to dict keys
    col_map = {
        'NUMBER': 'number',
        'X_IMAGE': 'x_image',
        'Y_IMAGE': 'y_image',
        'A_IMAGE': 'a_image',
        'B_IMAGE': 'b_image',
        'THETA_IMAGE': 'theta_image',
        'FLUX_AUTO': 'flux_auto',
        'KRON_RADIUS': 'kron_radius',
        'FLUX_RADIUS': 'flux_radius',
        'ELLIPTICITY': 'ellipticity',
        'CLASS_STAR': 'class_star',
        'FLAGS': 'flags',
        'MAG_AUTO': 'mag_auto',
        'FWHM_IMAGE': 'fwhm_image',
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


def load_fits_data(path):
    """Load FITS image data from primary or first image extension."""
    from astropy.io import fits
    with fits.open(path) as hdul:
        # Try primary HDU first
        if hdul[0].data is not None:
            data = hdul[0].data.astype(np.float64)
            header = hdul[0].header
        else:
            # Try first image extension
            for ext in hdul[1:]:
                if ext.data is not None and ext.data.ndim >= 2:
                    data = ext.data.astype(np.float64)
                    header = ext.header
                    break
            else:
                print("ERROR: No image data found in FITS file", file=sys.stderr)
                sys.exit(1)

        # Handle 3D+ data (take first slice)
        while data.ndim > 2:
            data = data[0]

        return data, header


def save_fits(data, header, path):
    """Save data as FITS with WCS header."""
    from astropy.io import fits
    hdu = fits.PrimaryHDU(data.astype(np.float32))
    # Copy WCS keywords
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


def mode_find_stars(args):
    """Find stars in catalog and output TSV."""
    if not args.catalog:
        print("ERROR: --catalog required for find_stars mode", file=sys.stderr)
        sys.exit(1)

    catalog = parse_catalog_tsv(args.catalog)
    if catalog is None:
        print("ERROR: Failed to parse catalog", file=sys.stderr)
        sys.exit(1)

    config = PSFDeconvConfig(
        class_star_thresh=args.class_star_thresh,
        max_ellipticity=args.max_ellipticity,
        fwhm_sigma=args.fwhm_sigma,
        min_flux_snr=args.min_flux_snr,
    )

    method = args.method
    if method == 'combined':
        mask, info = find_stars_combined(catalog, config)
        print(f"Star finding (combined): {info['n_combined']} stars, "
              f"median FWHM={info['median_fwhm']:.2f}", file=sys.stderr)
    elif method == 'class_star':
        mask = find_stars_class_star(catalog, config.class_star_thresh)
        info = {'n_combined': int(mask.sum())}
        print(f"Star finding (CLASS_STAR >= {config.class_star_thresh}): "
              f"{mask.sum()} stars", file=sys.stderr)
    elif method == 'fwhm':
        mask = find_stars_fwhm(catalog, config.fwhm_sigma)
        info = {'n_combined': int(mask.sum())}
        print(f"Star finding (FWHM clustering): {mask.sum()} stars",
              file=sys.stderr)
    else:
        print(f"ERROR: Unknown star finding method: {method}", file=sys.stderr)
        sys.exit(1)

    n_stars = int(mask.sum())
    if n_stars == 0:
        print("ERROR: No stars found", file=sys.stderr)
        sys.exit(1)

    # Output TSV: star indices and positions
    print(f"#PSF_STARS\tN_STARS={n_stars}\tMETHOD={method}")
    print("NUMBER\tX_IMAGE\tY_IMAGE\tA_IMAGE\tB_IMAGE\tCLASS_STAR\tELLIPTICITY")

    indices = np.where(mask)[0]
    for i in indices:
        num = int(catalog['number'][i])
        x = catalog['x_image'][i]
        y = catalog['y_image'][i]
        a = catalog.get('a_image', np.ones(len(mask)))[i]
        b = catalog.get('b_image', np.ones(len(mask)))[i]
        cs = catalog.get('class_star', np.zeros(len(mask)))[i]
        ell = catalog.get('ellipticity', np.zeros(len(mask)))[i]
        print(f"{num}\t{x:.4f}\t{y:.4f}\t{a:.4f}\t{b:.4f}\t{cs:.4f}\t{ell:.4f}")


def mode_build_psf(args):
    """Build PSF from star positions and save to FITS."""
    data, header = load_fits_data(args.fits)

    # Get star positions
    if args.star_indices and args.catalog:
        catalog = parse_catalog_tsv(args.catalog)
        if catalog is None:
            print("ERROR: Failed to parse catalog", file=sys.stderr)
            sys.exit(1)

        indices = [int(x) for x in args.star_indices.split(',')]
        # Match by NUMBER column
        positions = []
        for idx in indices:
            match = np.where(catalog['number'] == idx)[0]
            if len(match) > 0:
                i = match[0]
                # Convert from 1-based to 0-based
                x = catalog['x_image'][i] - 1.0
                y = catalog['y_image'][i] - 1.0
                positions.append((x, y))
            else:
                print(f"WARNING: Star NUMBER {idx} not found in catalog",
                      file=sys.stderr)

        if not positions:
            print("ERROR: No valid star positions found", file=sys.stderr)
            sys.exit(1)
    else:
        print("ERROR: --star-indices and --catalog required", file=sys.stderr)
        sys.exit(1)

    psf_size = args.psf_size
    method = args.psf_method

    print(f"Building PSF ({method}) from {len(positions)} stars, "
          f"size={psf_size}x{psf_size}", file=sys.stderr)

    if method == 'median':
        psf, n_used = build_psf_median(data, positions, psf_size)
        extra_info = {}
    elif method == 'mean':
        psf, n_used = build_psf_mean(data, positions, psf_size)
        extra_info = {}
    elif method == 'epsf':
        psf, n_used = build_psf_epsf(data, positions, psf_size)
        extra_info = {}
    elif method == 'gaussian':
        psf, n_used, extra_info = fit_gaussian_psf(data, positions, psf_size)
    elif method == 'moffat':
        psf, n_used, extra_info = fit_moffat_psf(data, positions, psf_size)
    else:
        print(f"ERROR: Unknown PSF method: {method}", file=sys.stderr)
        sys.exit(1)

    # Save PSF
    output = args.psf_output
    from astropy.io import fits as pf
    hdu = pf.PrimaryHDU(psf.astype(np.float32))
    hdu.header['PSFMETH'] = method
    hdu.header['PSFSIZE'] = psf_size
    hdu.header['NSTARS'] = n_used
    for k, v in extra_info.items():
        hdu.header[k[:8].upper()] = round(v, 6)
    hdu.writeto(output, overwrite=True)

    # Output info TSV
    info_parts = [f"METHOD={method}", f"SIZE={psf_size}", f"N_STARS={n_used}"]
    for k, v in extra_info.items():
        info_parts.append(f"{k.upper()}={v:.4f}")
    print(f"#PSF_BUILT\t" + "\t".join(info_parts))
    print(f"PSF saved to {output}", file=sys.stderr)


def mode_deconvolve(args):
    """Run deconvolution and save result to FITS."""
    data, header = load_fits_data(args.fits)

    # Load PSF
    if not args.psf:
        print("ERROR: --psf required for deconvolve mode", file=sys.stderr)
        sys.exit(1)

    from astropy.io import fits as pf
    with pf.open(args.psf) as hdul:
        psf = hdul[0].data.astype(np.float64)

    # Normalize PSF
    total = psf.sum()
    if total > 0:
        psf /= total

    algorithm = args.algorithm
    print(f"Deconvolving ({algorithm}) image {data.shape}, "
          f"PSF {psf.shape}", file=sys.stderr)

    if algorithm == 'rl':
        result = richardson_lucy(data, psf, iterations=args.iterations)
    elif algorithm == 'rl_accelerated':
        result = richardson_lucy_accelerated(data, psf,
                                             iterations=args.iterations)
    elif algorithm == 'rl_tv':
        result = richardson_lucy_tv(data, psf, iterations=args.iterations,
                                    tv_lambda=args.tv_lambda)
    elif algorithm == 'wiener':
        result = wiener_deconvolve(data, psf, nsr=args.wiener_nsr)
    elif algorithm == 'tikhonov':
        result = tikhonov_deconvolve(data, psf, lam=args.tikhonov_lambda)
    elif algorithm == 'clean':
        result, model, residual = clean_deconvolve(
            data, psf, gain=args.clean_gain, niter=args.clean_niter,
            threshold=args.clean_threshold)
    elif algorithm == 'mem':
        result = mem_deconvolve(data, psf, lam=args.mem_lambda,
                                niter=args.mem_niter)
    else:
        print(f"ERROR: Unknown algorithm: {algorithm}", file=sys.stderr)
        sys.exit(1)

    # Save result
    output = args.output
    save_fits(result, header, output)

    # Output info
    print(f"#DECONV_DONE\tALGORITHM={algorithm}\tOUTPUT={output}")
    print(f"Deconvolution complete: {output}", file=sys.stderr)


def mode_build_psf_extended(args):
    """Build extended PSF from core + wing stars."""
    data, header = load_fits_data(args.fits)

    if not args.catalog:
        print("ERROR: --catalog required for build_psf_extended mode",
              file=sys.stderr)
        sys.exit(1)

    catalog = parse_catalog_tsv(args.catalog)
    if catalog is None:
        print("ERROR: Failed to parse catalog", file=sys.stderr)
        sys.exit(1)

    if 'mag_auto' not in catalog:
        print("ERROR: Catalog must contain MAG_AUTO column", file=sys.stderr)
        sys.exit(1)

    config = PSFDeconvConfig(
        ext_core_mag_min=args.ext_core_mag_min,
        ext_core_mag_max=args.ext_core_mag_max,
        ext_wing_mag_max=args.ext_wing_mag_max,
        ext_core_size=args.ext_core_size,
        ext_wing_size=args.ext_wing_size,
        ext_blend_inner=args.ext_blend_inner,
        ext_blend_outer=args.ext_blend_outer,
        ext_saturation_limit=args.ext_saturation_limit,
    )

    print(f"Building extended PSF: core mag [{args.ext_core_mag_min}, "
          f"{args.ext_core_mag_max}], wing mag < {args.ext_wing_mag_max}",
          file=sys.stderr)

    psf, n_core, n_wing, info = build_psf_extended(data, catalog, config)

    # Save PSF
    output = args.psf_output
    from astropy.io import fits as pf
    hdu = pf.PrimaryHDU(psf.astype(np.float32))
    hdu.header['PSFMETH'] = 'extended'
    hdu.header['PSFSIZE'] = config.ext_wing_size
    hdu.header['NCORE'] = n_core
    hdu.header['NWING'] = n_wing
    hdu.header['CORESIZE'] = config.ext_core_size
    hdu.header['WINGSIZE'] = config.ext_wing_size
    hdu.header['BLENDIN'] = config.ext_blend_inner
    hdu.header['BLENDOUT'] = config.ext_blend_outer
    hdu.writeto(output, overwrite=True)

    print(f"#PSF_EXTENDED\tN_CORE={n_core}\tN_WING={n_wing}\t"
          f"SIZE={config.ext_wing_size}\tOUTPUT={output}")
    print(f"Extended PSF saved to {output}", file=sys.stderr)


def mode_sim_psf(args):
    """Generate simulation PSF (WebbPSF or TinyTim)."""
    from psf_deconv.psf.sim_psf import (
        detect_telescope_instrument, normalize_instrument_name,
        generate_webbpsf, generate_tinytim,
        check_webbpsf_available, check_tinytim_available,
        WEBBPSF_INSTRUMENTS, TINYTIM_INSTRUMENTS)

    telescope = args.sim_telescope.lower()
    instrument = args.sim_instrument
    filter_name = args.sim_filter
    output = args.psf_output

    # Auto-detect from FITS header if needed
    if telescope == 'auto' or instrument == 'auto' or filter_name == 'auto':
        _, header = load_fits_data(args.fits)
        detected = detect_telescope_instrument(header)
        print(f"Auto-detected: {detected}", file=sys.stderr)

        if telescope == 'auto':
            telescope = detected['telescope'].lower()
        if instrument == 'auto':
            instrument = detected['instrument']
        if filter_name == 'auto':
            filter_name = detected['filter']

    if telescope == 'unknown':
        print("ERROR: Cannot determine telescope from FITS header. "
              "Use --sim-telescope to specify.", file=sys.stderr)
        sys.exit(1)

    instrument = normalize_instrument_name(instrument, telescope)

    if telescope == 'jwst':
        if not check_webbpsf_available():
            print("ERROR: webbpsf package not installed. "
                  "Install with: pip install webbpsf", file=sys.stderr)
            sys.exit(1)

        print(f"Generating WebbPSF: {instrument} / {filter_name} "
              f"({args.sim_psf_size}px)", file=sys.stderr)
        psf = generate_webbpsf(
            instrument, filter_name,
            psf_size=args.sim_psf_size,
            oversample=args.sim_oversample,
            jitter_sigma=args.sim_jitter_sigma,
            focus_offset=args.sim_focus_offset,
            output=output)

        print(f"#PSF_SIM\tTELESCOPE=JWST\tINSTRUMENT={instrument}\t"
              f"FILTER={filter_name}\tSIZE={args.sim_psf_size}\tOUTPUT={output}")

    elif telescope == 'hst':
        if not check_tinytim_available():
            print("ERROR: TinyTim executables (tiny1, tiny2, tiny3) not found on PATH.",
                  file=sys.stderr)
            sys.exit(1)

        print(f"Generating TinyTim PSF: {instrument} / {filter_name} "
              f"({args.sim_psf_size}px)", file=sys.stderr)
        psf = generate_tinytim(
            instrument, filter_name,
            psf_size=args.sim_psf_size,
            oversample=args.sim_oversample,
            focus_offset=args.sim_focus_offset,
            output=output)

        print(f"#PSF_SIM\tTELESCOPE=HST\tINSTRUMENT={instrument}\t"
              f"FILTER={filter_name}\tSIZE={args.sim_psf_size}\tOUTPUT={output}")

    else:
        print(f"ERROR: Unsupported telescope: {telescope}. "
              f"Use 'jwst' or 'hst'.", file=sys.stderr)
        sys.exit(1)


def mode_check_sim(args):
    """Check simulation PSF tool availability."""
    from psf_deconv.psf.sim_psf import check_webbpsf_available, check_tinytim_available

    webbpsf_ok = 1 if check_webbpsf_available() else 0
    tinytim_ok = 1 if check_tinytim_available() else 0

    print(f"#SIM_STATUS\tWEBBPSF={webbpsf_ok}\tTINYTIM={tinytim_ok}")

    if webbpsf_ok:
        print("WebbPSF: available", file=sys.stderr)
    else:
        print("WebbPSF: NOT found (pip install webbpsf)", file=sys.stderr)
    if tinytim_ok:
        print("TinyTim: available", file=sys.stderr)
    else:
        print("TinyTim: NOT found (tiny1/tiny2/tiny3 not on PATH)", file=sys.stderr)


def main():
    parser = argparse.ArgumentParser(
        description='DS9 PSF Deconvolution Pipeline')
    parser.add_argument('fits', help='Input FITS file')
    parser.add_argument('--mode', required=True,
                        choices=['find_stars', 'build_psf', 'deconvolve',
                                 'build_psf_extended', 'sim_psf', 'check_sim'],
                        help='Operation mode')

    # Star finding
    parser.add_argument('--catalog', help='TSV catalog file')
    parser.add_argument('--method', default='combined',
                        help='Star finding method: combined, class_star, fwhm')
    parser.add_argument('--class-star-thresh', type=float, default=0.8)
    parser.add_argument('--max-ellipticity', type=float, default=0.2)
    parser.add_argument('--fwhm-sigma', type=float, default=2.0)
    parser.add_argument('--min-flux-snr', type=float, default=10.0)

    # PSF building
    parser.add_argument('--star-indices', help='Comma-separated star NUMBERs')
    parser.add_argument('--psf-method', default='median',
                        help='PSF method: median, mean, epsf, gaussian, moffat')
    parser.add_argument('--psf-size', type=int, default=51)
    parser.add_argument('--psf-output', default=os.path.expanduser(
                        '~/.ds9/psf_current.fits'))

    # Extended PSF
    parser.add_argument('--ext-core-mag-min', type=float, default=18.0)
    parser.add_argument('--ext-core-mag-max', type=float, default=22.0)
    parser.add_argument('--ext-wing-mag-max', type=float, default=16.0)
    parser.add_argument('--ext-core-size', type=int, default=51)
    parser.add_argument('--ext-wing-size', type=int, default=201)
    parser.add_argument('--ext-blend-inner', type=float, default=20.0)
    parser.add_argument('--ext-blend-outer', type=float, default=30.0)
    parser.add_argument('--ext-saturation-limit', type=float, default=60000.0)

    # Simulation PSF
    parser.add_argument('--sim-telescope', default='auto',
                        help='Telescope: auto, jwst, hst')
    parser.add_argument('--sim-instrument', default='auto',
                        help='Instrument: auto, nircam, miri, acs_wfc, ...')
    parser.add_argument('--sim-filter', default='auto',
                        help='Filter: auto, F150W, F814W, ...')
    parser.add_argument('--sim-psf-size', type=int, default=201)
    parser.add_argument('--sim-oversample', type=int, default=1)
    parser.add_argument('--sim-jitter-sigma', type=float, default=0.007)
    parser.add_argument('--sim-focus-offset', type=float, default=0.0)

    # Deconvolution
    parser.add_argument('--psf', help='PSF FITS file')
    parser.add_argument('--algorithm', default='rl',
                        help='Deconv algorithm: rl, rl_accelerated, rl_tv, '
                             'wiener, tikhonov, clean, mem')
    parser.add_argument('--iterations', type=int, default=30)
    parser.add_argument('--wiener-nsr', type=float, default=0.01)
    parser.add_argument('--tikhonov-lambda', type=float, default=0.001)
    parser.add_argument('--tv-lambda', type=float, default=0.001)
    parser.add_argument('--clean-gain', type=float, default=0.1)
    parser.add_argument('--clean-niter', type=int, default=1000)
    parser.add_argument('--clean-threshold', type=float, default=0.0)
    parser.add_argument('--mem-lambda', type=float, default=0.1)
    parser.add_argument('--mem-niter', type=int, default=100)
    parser.add_argument('--output', default=os.path.expanduser(
                        '~/.ds9/deconv_result.fits'))

    args = parser.parse_args()

    if args.mode == 'find_stars':
        mode_find_stars(args)
    elif args.mode == 'build_psf':
        mode_build_psf(args)
    elif args.mode == 'deconvolve':
        mode_deconvolve(args)
    elif args.mode == 'build_psf_extended':
        mode_build_psf_extended(args)
    elif args.mode == 'sim_psf':
        mode_sim_psf(args)
    elif args.mode == 'check_sim':
        mode_check_sim(args)


if __name__ == '__main__':
    main()
