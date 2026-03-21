#!/usr/bin/env python3
"""
Create inverse-variance weighted r+i+z coadd from DES DR2 tiles.

DES tiles share identical pixel grids (WCS), so no reprojection is needed.
Weight maps (HDU 3, EXTNAME=WGT) are per-pixel inverse-variance.
HDU layout: 0=Primary, 1=SCI, 2=MSK, 3=WGT.

Optimal coadd: coadd = sum(data_i * wgt_i) / sum(wgt_i)

Produces a cutout matching the g-band cutout footprint.

Usage: cd OGFinder && python3 scripts/make_riz_coadd.py
"""

import os
import sys
import numpy as np
from astropy.io import fits

DATADIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', '..', 'data')

BANDS = ['r', 'i', 'z']
TILE_TEMPLATE = os.path.join(DATADIR, 'des_tanoglidis2021_{band}.fits.fz')
G_CUTOUT = os.path.join(DATADIR, 'des_tanoglidis2021_g_cutout.fits')
OUTPUT = os.path.join(DATADIR, 'des_tanoglidis2021_riz_coadd.fits')


def load_tile(band):
    """Load science data (HDU 1/SCI) and weight map (HDU 3/WGT) from DES tile."""
    path = TILE_TEMPLATE.format(band=band)
    if not os.path.exists(path):
        print(f"ERROR: {band}-band tile not found: {path}", file=sys.stderr)
        print("Run: bash scripts/fetch_des_riz.sh", file=sys.stderr)
        sys.exit(1)

    with fits.open(path) as hdul:
        data = hdul[1].data.astype(np.float64)  # SCI
        wgt = hdul[3].data.astype(np.float64)   # WGT (HDU 2 is MSK)
        header = hdul[1].header.copy()
    return data, wgt, header


def get_cutout_slice(full_header, cutout_header):
    """Compute array slices to extract cutout region from full tile.

    Uses CRPIX shift to determine the offset.
    """
    # Full tile CRPIX
    crpix1_full = full_header['CRPIX1']
    crpix2_full = full_header['CRPIX2']
    # Cutout CRPIX
    crpix1_cut = cutout_header['CRPIX1']
    crpix2_cut = cutout_header['CRPIX2']
    # Cutout size
    nx_cut = cutout_header['NAXIS1']
    ny_cut = cutout_header['NAXIS2']

    # Offset: full_pixel = cutout_pixel + offset
    x_off = int(round(crpix1_full - crpix1_cut))
    y_off = int(round(crpix2_full - crpix2_cut))

    return slice(y_off, y_off + ny_cut), slice(x_off, x_off + nx_cut)


def verify_wcs(headers):
    """Verify that all band WCS are identical."""
    ref = headers[0]
    for key in ['CRPIX1', 'CRPIX2', 'CRVAL1', 'CRVAL2', 'CD1_1', 'CD2_2']:
        ref_val = ref[key]
        for i, h in enumerate(headers[1:], 1):
            if abs(h[key] - ref_val) > 1e-10:
                print(f"WARNING: WCS mismatch in {key}: "
                      f"band 0 = {ref_val}, band {i} = {h[key]}",
                      file=sys.stderr)
                return False
    return True


def main():
    if os.path.exists(OUTPUT):
        size_mb = os.path.getsize(OUTPUT) / 1e6
        print(f"Already exists: {OUTPUT} ({size_mb:.1f} MB)")
        print("Delete to regenerate.")
        return

    print("=== r+i+z Inverse-Variance Weighted Coadd ===", file=sys.stderr)

    # Load g-band cutout header for footprint
    cutout_header = None
    if os.path.exists(G_CUTOUT):
        with fits.open(G_CUTOUT) as h:
            cutout_header = h[0].header.copy()
        print(f"Cutout footprint from: {G_CUTOUT}", file=sys.stderr)
        print(f"  Size: {cutout_header['NAXIS1']} x {cutout_header['NAXIS2']}",
              file=sys.stderr)

    # Load bands
    data_list = []
    wgt_list = []
    hdr_list = []
    for band in BANDS:
        print(f"Loading {band}-band...", file=sys.stderr)
        data, wgt, hdr = load_tile(band)
        print(f"  Shape: {data.shape}, "
              f"data range: [{np.nanmin(data):.2f}, {np.nanmax(data):.2f}]",
              file=sys.stderr)
        data_list.append(data)
        wgt_list.append(wgt)
        hdr_list.append(hdr)

    # Verify WCS alignment
    if verify_wcs(hdr_list):
        print("WCS alignment verified: all bands share identical grid",
              file=sys.stderr)

    # Extract cutout region if g-band cutout exists
    if cutout_header is not None:
        sy, sx = get_cutout_slice(hdr_list[0], cutout_header)
        print(f"Extracting cutout: y[{sy.start}:{sy.stop}], "
              f"x[{sx.start}:{sx.stop}]", file=sys.stderr)
        for i in range(len(BANDS)):
            data_list[i] = data_list[i][sy, sx]
            wgt_list[i] = wgt_list[i][sy, sx]

    # Inverse-variance weighted coadd
    # wgt = inverse variance, so coadd = sum(data * wgt) / sum(wgt)
    print("Computing coadd...", file=sys.stderr)
    sum_dw = np.zeros_like(data_list[0])
    sum_w = np.zeros_like(data_list[0])

    for i, band in enumerate(BANDS):
        # Mask bad pixels (wgt <= 0 or non-finite)
        good = np.isfinite(data_list[i]) & np.isfinite(wgt_list[i]) & (wgt_list[i] > 0)
        d = np.where(good, data_list[i], 0.0)
        w = np.where(good, wgt_list[i], 0.0)
        sum_dw += d * w
        sum_w += w
        n_good = np.count_nonzero(good)
        print(f"  {band}: {n_good}/{good.size} good pixels "
              f"({100*n_good/good.size:.1f}%)", file=sys.stderr)

    # Avoid division by zero
    valid = sum_w > 0
    coadd = np.where(valid, sum_dw / sum_w, 0.0)

    print(f"Coadd shape: {coadd.shape}", file=sys.stderr)
    print(f"Coadd range: [{np.nanmin(coadd):.4f}, {np.nanmax(coadd):.4f}]",
          file=sys.stderr)
    print(f"Valid pixels: {np.count_nonzero(valid)}/{valid.size} "
          f"({100*np.count_nonzero(valid)/valid.size:.1f}%)", file=sys.stderr)

    # Write output with r-band WCS (or cutout WCS if available)
    out_header = cutout_header if cutout_header is not None else hdr_list[0]
    # Add coadd metadata
    out_header['COMMENT'] = 'DES DR2 r+i+z inverse-variance weighted coadd'
    out_header['COMMENT'] = 'Tile DES0125-0124, bands: r, i, z'
    out_header['NCOADD'] = (len(BANDS), 'Number of bands coadded')

    hdu = fits.PrimaryHDU(data=coadd.astype(np.float32), header=out_header)
    hdu.writeto(OUTPUT, overwrite=True)

    size_mb = os.path.getsize(OUTPUT) / 1e6
    print(f"\nSaved: {OUTPUT} ({size_mb:.1f} MB)", file=sys.stderr)
    print(f"Next: bash scripts/reproduce_tanoglidis2021_riz.sh", file=sys.stderr)


if __name__ == '__main__':
    main()
