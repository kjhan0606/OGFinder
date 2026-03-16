#!/usr/bin/env python3
"""Sérsic profile fitting for each source."""

import sys
import os
import argparse
import numpy as np

# Add project root to path
_script_dir = os.path.dirname(os.path.abspath(__file__))
_project_root = os.path.abspath(os.path.join(_script_dir, '..', '..'))
if _project_root not in sys.path:
    sys.path.insert(0, _project_root)


def main():
    parser = argparse.ArgumentParser(description='Sérsic profile fitting')
    parser.add_argument('fitsfile', help='Input FITS image')
    parser.add_argument('--catalog', required=True, help='Input TSV catalog')
    parser.add_argument('--max-sources', type=int, default=200)
    parser.add_argument('--mag-zeropoint', type=float, default=25.0)
    args = parser.parse_args()

    try:
        import sep
        from astropy.io import fits
    except ImportError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)

    from sersic_fit.config import SersicConfig
    from sersic_fit.fitter import fit_sersic

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

    # Background subtract
    mask = ~np.isfinite(data) | (data == 0)
    data[~np.isfinite(data)] = 0.0
    bkg = sep.Background(data, mask=mask.astype(np.uint8), bw=64, bh=64)
    data_sub = data - bkg

    ny, nx = data_sub.shape

    # Parse catalog
    with open(args.catalog, 'r') as f:
        lines = [l.strip() for l in f if l.strip()]
    headers = [h.strip() for h in lines[0].split('\t')]
    col_map = {h: i for i, h in enumerate(headers)}

    cfg = SersicConfig(mag_zeropoint=args.mag_zeropoint, max_sources=args.max_sources)

    print("NUMBER\tSERSIC_N\tSERSIC_RE\tSERSIC_IE\tSERSIC_ELLIP\tSERSIC_THETA\tSERSIC_CHI2")

    n_fit = 0
    for line in lines[1:]:
        if n_fit >= cfg.max_sources:
            break
        fields = line.split('\t')
        try:
            num = int(fields[col_map['NUMBER']].strip())
            x = float(fields[col_map['X_IMAGE']].strip()) - 1.0
            y = float(fields[col_map['Y_IMAGE']].strip()) - 1.0
            a = float(fields[col_map['A_IMAGE']].strip())
            b = float(fields[col_map['B_IMAGE']].strip())
            theta = float(fields[col_map['THETA_IMAGE']].strip()) * np.pi / 180.0
            kron_r = float(fields[col_map['KRON_RADIUS']].strip()) if 'KRON_RADIUS' in col_map else 3.5
            flux = float(fields[col_map['FLUX_AUTO']].strip()) if 'FLUX_AUTO' in col_map else 1000.0
        except (ValueError, IndexError, KeyError):
            continue

        # Extract cutout
        halfsize = int(cfg.cutout_scale * kron_r * max(a, 3.0))
        halfsize = max(halfsize, 15)
        halfsize = min(halfsize, 150)

        ix, iy = int(x), int(y)
        x0 = max(0, ix - halfsize)
        y0 = max(0, iy - halfsize)
        x1 = min(nx, ix + halfsize + 1)
        y1 = min(ny, iy + halfsize + 1)

        if x1 - x0 < 10 or y1 - y0 < 10:
            print(f"{num}\tnan\tnan\tnan\tnan\tnan\tnan")
            continue

        cutout = data_sub[y0:y1, x0:x1].copy()
        cx = x - x0
        cy = y - y0

        result = fit_sersic(cutout, cx, cy, a, b, theta, flux, cfg)

        if result is None:
            print(f"{num}\tnan\tnan\tnan\tnan\tnan\tnan")
        else:
            print(f"{num}\t{result['n']:.3f}\t{result['re']:.3f}\t{result['Ie']:.3f}\t"
                  f"{result['ellip']:.3f}\t{result['theta']:.1f}\t{result['chi2']:.4f}")

        n_fit += 1
        if n_fit % 20 == 0:
            print(f"  Fitted {n_fit} sources...", file=sys.stderr)

    print(f"Done: {n_fit} sources fitted", file=sys.stderr)


if __name__ == '__main__':
    main()
