#!/usr/bin/env python3
"""
Simple SEP-based source detection on r+i+z coadd.

Unlike the full LSBG pipeline (which aggressively masks and subtracts background,
removing moderate-brightness extended sources), this script does minimal processing:
- Single-pass SEP background subtraction
- SEP extract + Kron photometry
- Shape parameters and WCS coordinates
- Basic size/ellipticity filtering (no surface brightness cuts - those are for g-band)

The output catalog is compatible with ds9_lsbg.py --mode forced.

Usage: python3 scripts/sep_detect_riz.py <riz_coadd.fits> [options]
"""

import argparse
import math
import sys

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS


def main():
    parser = argparse.ArgumentParser(description='Simple SEP detection on riz coadd')
    parser.add_argument('fits', help='Input FITS (riz coadd)')
    parser.add_argument('--pixel-scale', type=float, default=0.263,
                        help='Pixel scale in arcsec/pixel')
    parser.add_argument('--mag-zeropoint', type=float, default=30.0,
                        help='Magnitude zeropoint')
    parser.add_argument('--detect-thresh', type=float, default=1.5,
                        help='Detection threshold in sigma')
    parser.add_argument('--detect-minarea', type=int, default=50,
                        help='Minimum detection area in pixels')
    parser.add_argument('--deblend-nthresh', type=int, default=32,
                        help='Number of deblending thresholds')
    parser.add_argument('--deblend-mincont', type=float, default=0.005,
                        help='Minimum contrast for deblending')
    parser.add_argument('--filter-kernel', default='gauss',
                        help='Detection filter kernel type')
    parser.add_argument('--filter-fwhm', type=float, default=5.0,
                        help='Filter kernel FWHM in pixels')
    parser.add_argument('--bkg-size', type=int, default=256,
                        help='Background mesh size')
    parser.add_argument('--r-eff-min', type=float, default=2.5,
                        help='Minimum r_eff in arcsec')
    parser.add_argument('--r-eff-max', type=float, default=20.0,
                        help='Maximum r_eff in arcsec')
    parser.add_argument('--ellipticity-max', type=float, default=0.7,
                        help='Maximum ellipticity')
    parser.add_argument('--output', '-o', default=None,
                        help='Output TSV file (default: stdout)')
    args = parser.parse_args()

    try:
        import sep
    except ImportError:
        import sep_pjw as sep

    # Load FITS
    with fits.open(args.fits) as hdul:
        data = hdul[0].data
        header = hdul[0].header
        if data is None and len(hdul) > 1:
            data = hdul[1].data
            header = hdul[1].header

    data = np.ascontiguousarray(data, dtype=np.float64)
    ny, nx = data.shape
    wcs = WCS(header)

    print(f"Image: {nx}x{ny}, pixel_scale={args.pixel_scale}\"/pix, "
          f"ZP={args.mag_zeropoint}", file=sys.stderr)

    # Background subtraction (single pass, large mesh)
    print(f"Background: mesh={args.bkg_size}...", file=sys.stderr)
    bkg = sep.Background(data, bw=args.bkg_size, bh=args.bkg_size)
    data_sub = data - bkg.back()
    rms_global = bkg.globalrms
    print(f"  Global RMS: {rms_global:.4f}", file=sys.stderr)

    # Build convolution filter (Gaussian)
    fwhm = args.filter_fwhm
    sigma = fwhm / 2.3548
    ksize = int(2 * math.ceil(3 * sigma) + 1)
    if ksize < 3:
        ksize = 3
    y_k, x_k = np.mgrid[:ksize, :ksize] - ksize // 2
    kernel = np.exp(-(x_k**2 + y_k**2) / (2 * sigma**2))
    kernel /= kernel.sum()

    # Source extraction
    sep.set_extract_pixstack(10000000)  # riz coadd has high S/N, needs big pixstack
    print(f"Extracting: thresh={args.detect_thresh}sigma, "
          f"minarea={args.detect_minarea}pix...", file=sys.stderr)
    try:
        objects = sep.extract(
            data_sub, args.detect_thresh,
            err=bkg.rms(),
            minarea=args.detect_minarea,
            filter_kernel=kernel,
            deblend_nthresh=args.deblend_nthresh,
            deblend_cont=args.deblend_mincont,
        )
    except Exception as e:
        print(f"ERROR: SEP extract failed: {e}", file=sys.stderr)
        sys.exit(1)

    n_raw = len(objects)
    print(f"  Raw detections: {n_raw}", file=sys.stderr)

    if n_raw == 0:
        print("No sources detected", file=sys.stderr)
        sys.exit(1)

    # Kron photometry
    x = objects['x']
    y = objects['y']
    a = objects['a']
    b = objects['b']
    theta = objects['theta']

    kronrad, krflag = sep.kron_radius(data_sub, x, y, a, b, theta, 6.0)
    # Clamp kronrad
    kronrad = np.where(np.isfinite(kronrad) & (kronrad > 0), kronrad, 2.5)
    kronrad = np.clip(kronrad, 1.0, 10.0)

    flux, fluxerr, flag = sep.sum_ellipse(
        data_sub, x, y, a, b, theta,
        2.5 * kronrad, err=bkg.rms(), subpix=5
    )

    # Magnitudes
    zp = args.mag_zeropoint
    with np.errstate(divide='ignore', invalid='ignore'):
        mag = np.where(flux > 0, -2.5 * np.log10(flux) + zp, 99.0)

    # Two r_eff measures:
    # 1. Second-moment: sqrt(a*b)*ps - used for size filtering (stable, no blending issues)
    # 2. flux_radius: SEP half-light radius - used for mu_eff calculation (more accurate)
    r_eff_moment = np.sqrt(a * b) * args.pixel_scale

    max_radius = np.clip(6.0 * np.sqrt(a * b), 20.0, 200.0)
    r50, r50_flag = sep.flux_radius(
        data_sub, x, y, max_radius, 0.5, subpix=5
    )
    r_eff_flux = r50 * args.pixel_scale
    r_eff_flux = np.where(
        np.isfinite(r_eff_flux) & (r_eff_flux > 0),
        r_eff_flux, r_eff_moment
    )

    # Use second-moment r_eff for filtering (more stable)
    r_eff_arcsec = r_eff_moment

    # Ellipticity
    ellipticity = 1.0 - b / np.where(a > 0, a, 1.0)

    # Surface brightness: mu_eff = mag + 2.5*log10(2*pi*r_eff^2)
    with np.errstate(divide='ignore', invalid='ignore'):
        mu_eff = np.where(
            (r_eff_arcsec > 0) & (mag < 90),
            mag + 2.5 * np.log10(2 * np.pi * r_eff_arcsec**2),
            99.0
        )

    # WCS coordinates
    # X_IMAGE, Y_IMAGE are 1-indexed (DS9/SExtractor convention)
    x_img = x + 1.0
    y_img = y + 1.0
    try:
        coords = wcs.pixel_to_world_values(x, y)
        ra = coords[0]
        dec = coords[1]
    except Exception:
        ra = np.full(n_raw, np.nan)
        dec = np.full(n_raw, np.nan)

    # Filter: size and ellipticity only (NO surface brightness cut on riz)
    # Use moment-based r_eff for filtering (more stable for blended sources)
    ps = args.pixel_scale
    mask = (
        (r_eff_arcsec >= args.r_eff_min) &
        (r_eff_arcsec <= args.r_eff_max) &
        (ellipticity <= args.ellipticity_max) &
        (flux > 0) &
        np.isfinite(mag) &
        (flag < 8)  # allow minor flags, reject only corrupted
    )

    # Cap A/B for output catalog to prevent forced photometry from using
    # absurdly large apertures on blended sources
    MAX_AB_PIX = 30.0  # ~8" at 0.263"/pix
    a_out = np.clip(a, 0, MAX_AB_PIX)
    b_out = np.clip(b, 0, MAX_AB_PIX)

    indices = np.where(mask)[0]
    n_pass = len(indices)
    print(f"  After size/ellipticity filter: {n_pass} sources", file=sys.stderr)
    print(f"    r_eff: {args.r_eff_min}-{args.r_eff_max}\", "
          f"e<{args.ellipticity_max}", file=sys.stderr)

    # Output TSV
    out = sys.stdout
    if args.output:
        os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
        out = open(args.output, 'w')

    columns = [
        'NUMBER', 'X_IMAGE', 'Y_IMAGE', 'RA', 'DEC',
        'A_IMAGE', 'B_IMAGE', 'THETA_IMAGE', 'ELLIPTICITY',
        'FLUX_AUTO', 'FLUXERR_AUTO', 'MAG_AUTO', 'KRON_RADIUS',
        'R_EFF_ARCSEC', 'FLUX_RADIUS_ARCSEC', 'MU_EFF',
    ]
    out.write('\t'.join(columns) + '\n')

    for num, idx in enumerate(indices, 1):
        vals = [
            str(num),
            f"{x_img[idx]:.2f}",
            f"{y_img[idx]:.2f}",
            f"{ra[idx]:.6f}",
            f"{dec[idx]:.6f}",
            f"{a_out[idx]:.3f}",
            f"{b_out[idx]:.3f}",
            f"{theta[idx]:.4f}",
            f"{ellipticity[idx]:.4f}",
            f"{flux[idx]:.6e}",
            f"{fluxerr[idx]:.6e}",
            f"{mag[idx]:.4f}",
            f"{kronrad[idx]:.3f}",
            f"{r_eff_arcsec[idx]:.3f}",
            f"{r_eff_flux[idx]:.3f}",
            f"{mu_eff[idx]:.3f}",
        ]
        out.write('\t'.join(vals) + '\n')

    if args.output:
        out.close()
        print(f"Saved: {args.output} ({n_pass} sources)", file=sys.stderr)
    else:
        print(f"Output: {n_pass} sources to stdout", file=sys.stderr)


if __name__ == '__main__':
    main()
