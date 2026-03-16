#!/usr/bin/env python3
"""Multi-band photometry: measure fluxes across multiple band images."""

import sys
import os
import argparse
import numpy as np


def main():
    parser = argparse.ArgumentParser(description='Multi-band photometry')
    parser.add_argument('--detect-image', required=True, help='Detection FITS image')
    parser.add_argument('--bands', required=True,
                        help='Comma-separated band:file pairs (e.g., F435W:f435w.fits,F606W:f606w.fits)')
    parser.add_argument('--catalog', default=None,
                        help='Use existing catalog (TSV) instead of detecting')
    parser.add_argument('--detect-thresh', type=float, default=1.5)
    parser.add_argument('--detect-minarea', type=int, default=5)
    parser.add_argument('--phot-aperture', type=float, default=5.0)
    parser.add_argument('--mag-zeropoint', type=float, default=25.0)
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

    def load_fits(path):
        with fits.open(path) as hdul:
            for hdu in hdul:
                if hdu.data is not None and hdu.data.ndim >= 2:
                    return hdu.data.astype(np.float64)
        return None

    # Parse band list
    bands = []
    for item in args.bands.split(','):
        item = item.strip()
        if ':' in item:
            name, fpath = item.split(':', 1)
        else:
            name = os.path.splitext(os.path.basename(item))[0]
            fpath = item
        bands.append((name.strip(), fpath.strip()))

    # Detection
    det_data = load_fits(args.detect_image)
    if det_data is None:
        print("ERROR: Cannot read detection image", file=sys.stderr)
        sys.exit(1)

    mask = ~np.isfinite(det_data) | (det_data == 0)
    det_data[~np.isfinite(det_data)] = 0.0
    det_bkg = sep.Background(det_data, mask=mask.astype(np.uint8), bw=64, bh=64)
    det_sub = det_data - det_bkg

    if args.catalog:
        # Parse existing catalog
        with open(args.catalog, 'r') as f:
            lines = [l.strip() for l in f if l.strip()]
        headers = [h.strip() for h in lines[0].split('\t')]
        col_map = {h: i for i, h in enumerate(headers)}
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
        objects = sep.extract(det_sub, args.detect_thresh, err=det_bkg.rms(),
                              mask=mask.astype(np.uint8), minarea=args.detect_minarea)
        x_arr = objects['x']
        y_arr = objects['y']
        a_arr = objects['a']
        b_arr = objects['b']
        theta_arr = objects['theta']

    nobj = len(x_arr)
    zp = args.mag_zeropoint
    r_aper = args.phot_aperture / 2.0

    # Kron radius from detection image
    kronrad, _ = sep.kron_radius(det_sub, x_arr, y_arr,
                                   a_arr, b_arr, theta_arr, 6.0)
    kronrad = np.maximum(kronrad, 3.5)

    # Header
    hdr_parts = ["NUMBER", "X_IMAGE", "Y_IMAGE"]
    for bname, _ in bands:
        hdr_parts.extend([f"FLUX_AUTO_{bname}", f"MAG_AUTO_{bname}",
                          f"FLUX_APER_{bname}", f"MAG_APER_{bname}"])
    # Color indices
    for i in range(len(bands) - 1):
        hdr_parts.append(f"COLOR_{bands[i][0]}_{bands[i+1][0]}")

    print('\t'.join(hdr_parts))

    # Measure each band
    band_results = {}
    for bname, bpath in bands:
        bdata = load_fits(bpath)
        if bdata is None:
            print(f"WARNING: Cannot load {bpath}, skipping", file=sys.stderr)
            band_results[bname] = {
                'flux_auto': np.full(nobj, np.nan),
                'mag_auto': np.full(nobj, 99.0),
                'flux_aper': np.full(nobj, np.nan),
                'mag_aper': np.full(nobj, 99.0),
            }
            continue

        bmask = ~np.isfinite(bdata) | (bdata == 0)
        bdata[~np.isfinite(bdata)] = 0.0
        bbkg = sep.Background(bdata, mask=bmask.astype(np.uint8), bw=64, bh=64)
        bsub = bdata - bbkg

        # Kron photometry
        flux_auto, _, _ = sep.sum_ellipse(bsub, x_arr, y_arr,
                                           a_arr, b_arr, theta_arr,
                                           2.5 * kronrad, subpix=5)
        mag_auto = np.where(flux_auto > 0, -2.5 * np.log10(flux_auto) + zp, 99.0)

        # Aperture photometry
        flux_aper, _, _ = sep.sum_circle(bsub, x_arr, y_arr, r_aper, subpix=5)
        mag_aper = np.where(flux_aper > 0, -2.5 * np.log10(flux_aper) + zp, 99.0)

        band_results[bname] = {
            'flux_auto': flux_auto,
            'mag_auto': mag_auto,
            'flux_aper': flux_aper,
            'mag_aper': mag_aper,
        }
        print(f"  Band {bname}: measured {nobj} sources", file=sys.stderr)

    # Sort by detection magnitude (first band)
    first_band = bands[0][0]
    sort_idx = np.argsort(band_results[first_band]['mag_auto'])

    # Output data
    for num, idx in enumerate(sort_idx, 1):
        parts = [str(num), f"{x_arr[idx]+1:.2f}", f"{y_arr[idx]+1:.2f}"]
        for bname, _ in bands:
            br = band_results[bname]
            parts.append(f"{br['flux_auto'][idx]:.2f}")
            parts.append(f"{br['mag_auto'][idx]:.3f}")
            parts.append(f"{br['flux_aper'][idx]:.2f}")
            parts.append(f"{br['mag_aper'][idx]:.3f}")
        # Colors
        for i in range(len(bands) - 1):
            m1 = band_results[bands[i][0]]['mag_auto'][idx]
            m2 = band_results[bands[i+1][0]]['mag_auto'][idx]
            color = m1 - m2 if m1 < 90 and m2 < 90 else 99.0
            parts.append(f"{color:.3f}")
        print('\t'.join(parts))

    print(f"Done: {nobj} sources, {len(bands)} bands", file=sys.stderr)


if __name__ == '__main__':
    main()
