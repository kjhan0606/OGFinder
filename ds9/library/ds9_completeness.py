#!/usr/bin/env python3
"""Completeness simulation: inject artificial sources and measure recovery rate."""

import sys
import os
import argparse
import numpy as np


def main():
    parser = argparse.ArgumentParser(description='Completeness simulation')
    parser.add_argument('fitsfile', help='Input FITS image')
    parser.add_argument('--n-inject', type=int, default=1000,
                        help='Number of sources to inject per mag bin (default: 1000)')
    parser.add_argument('--mag-min', type=float, default=20.0,
                        help='Minimum magnitude (default: 20)')
    parser.add_argument('--mag-max', type=float, default=28.0,
                        help='Maximum magnitude (default: 28)')
    parser.add_argument('--n-bins', type=int, default=16,
                        help='Number of magnitude bins (default: 16)')
    parser.add_argument('--detect-thresh', type=float, default=1.5)
    parser.add_argument('--detect-minarea', type=int, default=5)
    parser.add_argument('--mag-zeropoint', type=float, default=25.0)
    parser.add_argument('--psf-fwhm', type=float, default=3.0,
                        help='PSF FWHM in pixels for injected sources')
    parser.add_argument('--match-radius', type=float, default=2.0,
                        help='Match radius in pixels (default: 2)')
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

    ny, nx = data.shape
    mask = ~np.isfinite(data) | (data == 0)
    data[~np.isfinite(data)] = 0.0

    # Background
    bkg = sep.Background(data, mask=mask.astype(np.uint8), bw=64, bh=64)
    data_sub = data - bkg
    rms = bkg.rms()

    # Get existing sources (to avoid placing near them)
    objects = sep.extract(data_sub, args.detect_thresh, err=rms,
                           mask=mask.astype(np.uint8), minarea=args.detect_minarea)
    exist_x = objects['x']
    exist_y = objects['y']

    # Create segmentation map to avoid existing sources
    _, segmap = sep.extract(data_sub, args.detect_thresh, err=rms,
                             mask=mask.astype(np.uint8), minarea=args.detect_minarea,
                             segmentation_map=True)

    # Gaussian PSF for injection
    sigma = args.psf_fwhm / 2.3548
    psf_size = int(6 * sigma) | 1  # odd
    psf_half = psf_size // 2
    yy, xx = np.mgrid[-psf_half:psf_half+1, -psf_half:psf_half+1]
    psf = np.exp(-(xx**2 + yy**2) / (2 * sigma**2))
    psf /= np.sum(psf)

    mag_bins = np.linspace(args.mag_min, args.mag_max, args.n_bins + 1)
    mag_centers = 0.5 * (mag_bins[:-1] + mag_bins[1:])

    print("MAG_BIN\tN_INJECTED\tN_RECOVERED\tCOMPLETENESS")
    rng = np.random.RandomState(42)

    for ibin in range(args.n_bins):
        mag_lo = mag_bins[ibin]
        mag_hi = mag_bins[ibin + 1]
        mag_c = mag_centers[ibin]

        # Inject sources
        injected = []
        sim_data = data.copy()
        n_placed = 0

        for _ in range(args.n_inject * 10):  # try more positions
            if n_placed >= args.n_inject:
                break
            ix = rng.randint(psf_half + 10, nx - psf_half - 10)
            iy = rng.randint(psf_half + 10, ny - psf_half - 10)

            # Skip masked/existing
            if mask[iy, ix]:
                continue
            if segmap[iy, ix] > 0:
                continue

            # Random magnitude in bin
            mag = rng.uniform(mag_lo, mag_hi)
            flux = 10.0**((args.mag_zeropoint - mag) / 2.5)

            # Inject
            y0 = iy - psf_half
            x0 = ix - psf_half
            sim_data[y0:y0+psf_size, x0:x0+psf_size] += flux * psf
            injected.append((ix, iy, mag))
            n_placed += 1

        if n_placed == 0:
            print(f"{mag_c:.2f}\t0\t0\t0.000")
            continue

        # Re-extract
        sim_mask = ~np.isfinite(sim_data) | (sim_data == 0)
        sim_data[~np.isfinite(sim_data)] = 0.0
        sim_bkg = sep.Background(sim_data, mask=sim_mask.astype(np.uint8), bw=64, bh=64)
        sim_sub = sim_data - sim_bkg
        sim_rms = sim_bkg.rms()

        rec_objects = sep.extract(sim_sub, args.detect_thresh, err=sim_rms,
                                   mask=sim_mask.astype(np.uint8),
                                   minarea=args.detect_minarea)
        rec_x = rec_objects['x']
        rec_y = rec_objects['y']

        # Match injected to recovered
        n_recovered = 0
        match_r2 = args.match_radius**2
        for inj_x, inj_y, _ in injected:
            dx = rec_x - inj_x
            dy = rec_y - inj_y
            d2 = dx**2 + dy**2
            if len(d2) > 0 and np.min(d2) <= match_r2:
                n_recovered += 1

        completeness = n_recovered / n_placed if n_placed > 0 else 0.0
        print(f"{mag_c:.2f}\t{n_placed}\t{n_recovered}\t{completeness:.4f}")
        print(f"  Mag {mag_c:.1f}: {n_recovered}/{n_placed} = {completeness:.1%}", file=sys.stderr)

    print("Done", file=sys.stderr)


if __name__ == '__main__':
    main()
