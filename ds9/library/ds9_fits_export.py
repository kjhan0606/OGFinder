#!/usr/bin/env python3
"""Export a TSV catalog to FITS binary table format."""

import sys
import os
import argparse

def main():
    parser = argparse.ArgumentParser(description='Export TSV catalog to FITS binary table')
    parser.add_argument('--input', required=True, help='Input TSV catalog file')
    parser.add_argument('--output', required=True, help='Output FITS file')
    args = parser.parse_args()

    try:
        from astropy.io import fits
        from astropy.table import Table
    except ImportError:
        print("ERROR: astropy is required for FITS export", file=sys.stderr)
        sys.exit(1)

    if not os.path.exists(args.input):
        print(f"ERROR: Input file not found: {args.input}", file=sys.stderr)
        sys.exit(1)

    # Read TSV
    with open(args.input, 'r') as f:
        lines = [l.strip() for l in f if l.strip()]

    if len(lines) < 2:
        print("ERROR: No data rows in input", file=sys.stderr)
        sys.exit(1)

    headers = lines[0].split('\t')
    ncols = len(headers)

    # Parse data rows
    columns = {h: [] for h in headers}
    for line in lines[1:]:
        fields = line.split('\t')
        for i, h in enumerate(headers):
            val = fields[i].strip() if i < len(fields) else ''
            columns[h].append(val)

    # Auto-detect column types and build FITS columns
    fits_cols = []
    for h in headers:
        vals = columns[h]
        # Try int first, then float, then string
        try:
            int_vals = [int(v) for v in vals]
            fits_cols.append(fits.Column(name=h, format='J', array=int_vals))
            continue
        except (ValueError, TypeError):
            pass
        try:
            float_vals = [float(v) for v in vals]
            fits_cols.append(fits.Column(name=h, format='D', array=float_vals))
            continue
        except (ValueError, TypeError):
            pass
        # String
        maxlen = max((len(v) for v in vals), default=1)
        if maxlen < 1:
            maxlen = 1
        fits_cols.append(fits.Column(name=h, format=f'{maxlen}A', array=vals))

    # Create FITS binary table
    coldefs = fits.ColDefs(fits_cols)
    hdu = fits.BinTableHDU.from_columns(coldefs)
    hdu.header['EXTNAME'] = 'CATALOG'

    hdulist = fits.HDUList([fits.PrimaryHDU(), hdu])
    hdulist.writeto(args.output, overwrite=True)

    nrows = len(lines) - 1
    print(f"OK {nrows} {ncols}")

if __name__ == '__main__':
    main()
