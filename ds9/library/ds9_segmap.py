#!/usr/bin/env python3
"""Generate a segmentation map from a FITS image using SEP."""

import sys
import os
import argparse
import numpy as np

def main():
    parser = argparse.ArgumentParser(description='Generate segmentation map')
    parser.add_argument('fitsfile', help='Input FITS image')
    parser.add_argument('--detect-thresh', type=float, default=1.5,
                        help='Detection threshold in sigma (default: 1.5)')
    parser.add_argument('--detect-minarea', type=int, default=5,
                        help='Minimum detection area in pixels (default: 5)')
    parser.add_argument('--output', default=None,
                        help='Output segmentation map FITS (default: ~/.ds9/segmap.fits)')
    args = parser.parse_args()

    try:
        import sep
    except ImportError:
        try:
            import sep_pjw as sep
        except ImportError:
            print("ERROR: sep is required", file=sys.stderr)
            sys.exit(1)

    try:
        from astropy.io import fits
    except ImportError:
        print("ERROR: astropy is required", file=sys.stderr)
        sys.exit(1)

    if not os.path.exists(args.fitsfile):
        print(f"ERROR: File not found: {args.fitsfile}", file=sys.stderr)
        sys.exit(1)

    # Output path
    if args.output is None:
        outdir = os.path.join(os.path.expanduser('~'), '.ds9')
        os.makedirs(outdir, exist_ok=True)
        args.output = os.path.join(outdir, 'segmap.fits')

    # Load FITS
    with fits.open(args.fitsfile) as hdul:
        data = None
        header = None
        for hdu in hdul:
            if hdu.data is not None and hdu.data.ndim >= 2:
                data = hdu.data.astype(np.float64)
                header = hdu.header
                break
        if data is None:
            print("ERROR: No 2D image data found", file=sys.stderr)
            sys.exit(1)

    # Handle NaN/Inf
    mask = ~np.isfinite(data) | (data == 0)
    data[~np.isfinite(data)] = 0.0

    # Background subtraction
    bkg = sep.Background(data, mask=mask.astype(np.uint8),
                         bw=64, bh=64, fw=3, fh=3)
    data_sub = data - bkg

    # Extract with segmentation map
    objects, segmap = sep.extract(data_sub, args.detect_thresh,
                                  err=bkg.rms(),
                                  mask=mask.astype(np.uint8),
                                  minarea=args.detect_minarea,
                                  segmentation_map=True)

    # Save segmentation map
    hdu_out = fits.PrimaryHDU(segmap.astype(np.int32))
    # Copy WCS if available
    for key in ['CRPIX1', 'CRPIX2', 'CRVAL1', 'CRVAL2',
                'CD1_1', 'CD1_2', 'CD2_1', 'CD2_2',
                'CDELT1', 'CDELT2', 'CROTA2', 'CTYPE1', 'CTYPE2',
                'NAXIS1', 'NAXIS2']:
        if key in header:
            hdu_out.header[key] = header[key]
    hdu_out.header['NSOURCES'] = len(objects)
    hdu_out.writeto(args.output, overwrite=True)

    nobj = len(objects)
    print(f"OK {nobj} {args.output}")

if __name__ == '__main__':
    main()
