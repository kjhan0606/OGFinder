#!/usr/bin/env python3
"""Dual-image mode: detect on one image, measure on another."""

import sys
import os
import argparse
import numpy as np


def main():
    parser = argparse.ArgumentParser(description='Dual-image source extraction')
    parser.add_argument('--detect-image', required=True, help='Detection FITS image')
    parser.add_argument('--measure-image', required=True, help='Measurement FITS image')
    parser.add_argument('--detect-thresh', type=float, default=1.5,
                        help='Detection threshold in sigma')
    parser.add_argument('--detect-minarea', type=int, default=5,
                        help='Minimum detection area')
    parser.add_argument('--deblend-nthresh', type=int, default=32)
    parser.add_argument('--deblend-mincont', type=float, default=0.005)
    parser.add_argument('--phot-aperture', type=float, default=5.0,
                        help='Aperture diameter in pixels')
    parser.add_argument('--mag-zeropoint', type=float, default=25.0)
    parser.add_argument('--catalog', default=None,
                        help='Use existing catalog positions instead of detecting')
    args = parser.parse_args()

    try:
        import sep
    except ImportError:
        print("ERROR: sep is required", file=sys.stderr)
        sys.exit(1)

    try:
        from astropy.io import fits
        from astropy.wcs import WCS
    except ImportError:
        print("ERROR: astropy is required", file=sys.stderr)
        sys.exit(1)

    def load_fits(path):
        with fits.open(path) as hdul:
            for hdu in hdul:
                if hdu.data is not None and hdu.data.ndim >= 2:
                    return hdu.data.astype(np.float64), hdu.header
        return None, None

    # Load images
    det_data, det_header = load_fits(args.detect_image)
    meas_data, meas_header = load_fits(args.measure_image)

    if det_data is None:
        print("ERROR: Cannot read detection image", file=sys.stderr)
        sys.exit(1)
    if meas_data is None:
        print("ERROR: Cannot read measurement image", file=sys.stderr)
        sys.exit(1)

    # WCS from detection image
    try:
        wcs = WCS(det_header)
        has_wcs = wcs.has_celestial
    except Exception:
        has_wcs = False
        wcs = None

    # Background subtraction for both
    def subtract_bg(data):
        mask = ~np.isfinite(data) | (data == 0)
        data[~np.isfinite(data)] = 0.0
        bkg = sep.Background(data, mask=mask.astype(np.uint8), bw=64, bh=64, fw=3, fh=3)
        return data - bkg, bkg, mask

    det_sub, det_bkg, det_mask = subtract_bg(det_data)
    meas_sub, meas_bkg, meas_mask = subtract_bg(meas_data)

    # Detection
    if args.catalog:
        # Use existing catalog
        with open(args.catalog, 'r') as f:
            lines = [l.strip() for l in f if l.strip()]
        headers = lines[0].split('\t')
        col_map = {h.strip(): i for i, h in enumerate(headers)}

        x_arr, y_arr, a_arr, b_arr, theta_arr = [], [], [], [], []
        for line in lines[1:]:
            fields = line.split('\t')
            x_arr.append(float(fields[col_map['X_IMAGE']]) - 1.0)
            y_arr.append(float(fields[col_map['Y_IMAGE']]) - 1.0)
            a_arr.append(float(fields[col_map['A_IMAGE']]))
            b_arr.append(float(fields[col_map['B_IMAGE']]))
            theta_arr.append(float(fields[col_map['THETA_IMAGE']]) * np.pi / 180.0)
        x_arr = np.array(x_arr)
        y_arr = np.array(y_arr)
        a_arr = np.array(a_arr)
        b_arr = np.array(b_arr)
        theta_arr = np.array(theta_arr)
    else:
        # Extract from detection image
        det_rms = det_bkg.rms()
        objects = sep.extract(det_sub, args.detect_thresh, err=det_rms,
                              mask=det_mask.astype(np.uint8),
                              minarea=args.detect_minarea,
                              deblend_nthresh=args.deblend_nthresh,
                              deblend_cont=args.deblend_mincont)
        x_arr = objects['x']
        y_arr = objects['y']
        a_arr = objects['a']
        b_arr = objects['b']
        theta_arr = objects['theta']

    nobj = len(x_arr)
    print(f"Measuring {nobj} sources on measurement image...", file=sys.stderr)

    # Measure on measurement image
    meas_rms = meas_bkg.rms()

    # Kron photometry on measurement image
    kronrad, krflag = sep.kron_radius(meas_sub, x_arr, y_arr,
                                        a_arr, b_arr, theta_arr, 6.0)
    kronrad = np.maximum(kronrad, 3.5)

    flux_auto, fluxerr_auto, flag_auto = sep.sum_ellipse(
        meas_sub, x_arr, y_arr, a_arr, b_arr, theta_arr,
        2.5 * kronrad, subpix=5)

    # Aperture photometry on measurement image
    r_aper = args.phot_aperture / 2.0
    flux_aper, fluxerr_aper, flag_aper = sep.sum_circle(
        meas_sub, x_arr, y_arr, r_aper, subpix=5)

    # Half-light radius
    rmax = np.maximum(6.0 * a_arr, 10.0)
    rmax = np.minimum(rmax, 200.0)
    flux_radius, fr_flag = sep.flux_radius(meas_sub, x_arr, y_arr,
                                            rmax, 0.5, subpix=5)

    # WCS
    if has_wcs:
        world = wcs.pixel_to_world(x_arr, y_arr)
        ra_arr = world.ra.deg
        dec_arr = world.dec.deg
    else:
        ra_arr = np.zeros(nobj)
        dec_arr = np.zeros(nobj)

    # Magnitudes
    zp = args.mag_zeropoint
    mag_auto = np.where(flux_auto > 0, -2.5 * np.log10(flux_auto) + zp, 99.0)
    mag_aper = np.where(flux_aper > 0, -2.5 * np.log10(flux_aper) + zp, 99.0)
    magerr_auto = np.where(flux_auto > 0, 1.0857 * fluxerr_auto / np.maximum(flux_auto, 1e-30), 99.0)

    # Morphology
    ellip = np.where(a_arr > 0, 1.0 - b_arr / a_arr, 0.0)
    theta_deg = theta_arr * 180.0 / np.pi

    # Sort by brightness
    sort_idx = np.argsort(mag_auto)

    # Output header
    print("NUMBER\tX_IMAGE\tY_IMAGE\tALPHA_J2000\tDELTA_J2000\t"
          "MAG_AUTO\tMAG_APER\tFLUX_AUTO\t"
          "FLUXERR_AUTO\tMAGERR_AUTO\t"
          "FLUX_RADIUS\tA_IMAGE\tB_IMAGE\tTHETA_IMAGE\tELLIPTICITY\t"
          "KRON_RADIUS\tFLAGS")

    for num, idx in enumerate(sort_idx, 1):
        print(f"{num}\t{x_arr[idx]+1:.2f}\t{y_arr[idx]+1:.2f}\t"
              f"{ra_arr[idx]:.6f}\t{dec_arr[idx]:.6f}\t"
              f"{mag_auto[idx]:.3f}\t{mag_aper[idx]:.3f}\t{flux_auto[idx]:.2f}\t"
              f"{fluxerr_auto[idx]:.2f}\t{magerr_auto[idx]:.4f}\t"
              f"{flux_radius[idx]:.2f}\t{a_arr[idx]:.2f}\t{b_arr[idx]:.2f}\t"
              f"{theta_deg[idx]:.1f}\t{ellip[idx]:.3f}\t"
              f"{kronrad[idx]:.2f}\t{int(flag_auto[idx])}")

    print(f"Done: {nobj} sources measured", file=sys.stderr)


if __name__ == '__main__':
    main()
