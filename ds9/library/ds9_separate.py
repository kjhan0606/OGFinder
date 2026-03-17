#!/usr/bin/env python3
"""Separate (deblend) a single source into sub-components using SEP.

Given a FITS image and a source's position/size, extract a cutout,
run SEP with aggressive deblending parameters, and output the
sub-sources as TSV to stdout.

Usage:
    python3 ds9_separate.py image.fits \
        --x 150.3 --y 200.1 --a 12.5 --b 8.3 --theta 45.0 \
        --radius-factor 3.0 --deblend-nthresh 64 --deblend-mincont 0.0001 \
        --detect-thresh 0.5 --parent-number 42
"""

import argparse
import sys
import numpy as np

try:
    import sep
except ImportError:
    try:
        import sep_pjw as sep
    except ImportError:
        sep = None


def main():
    parser = argparse.ArgumentParser(description="Separate a source into sub-components")
    parser.add_argument("fits", help="Input FITS file")
    parser.add_argument("--x", type=float, required=True, help="Source X (1-indexed image)")
    parser.add_argument("--y", type=float, required=True, help="Source Y (1-indexed image)")
    parser.add_argument("--a", type=float, default=10.0, help="Semi-major axis (pixels)")
    parser.add_argument("--b", type=float, default=10.0, help="Semi-minor axis (pixels)")
    parser.add_argument("--theta", type=float, default=0.0, help="Position angle (degrees)")
    parser.add_argument("--radius-factor", type=float, default=3.0,
                        help="Cutout radius = factor * max(a,b)")
    parser.add_argument("--deblend-nthresh", type=int, default=64,
                        help="Deblending sub-thresholds")
    parser.add_argument("--deblend-mincont", type=float, default=0.0001,
                        help="Deblending minimum contrast")
    parser.add_argument("--detect-thresh", type=float, default=0.8,
                        help="Detection threshold (sigma)")
    parser.add_argument("--detect-minarea", type=int, default=3,
                        help="Minimum detection area (pixels)")
    parser.add_argument("--parent-number", type=int, default=0,
                        help="Parent source NUMBER for sub-source numbering")
    parser.add_argument("--back-size", type=int, default=32,
                        help="Background mesh size")
    parser.add_argument("--iso-radius", type=float, default=0.0,
                        help="Parent ISO_RADIUS (0 = use A_IMAGE)")
    args = parser.parse_args()

    if sep is None:
        print("ERROR: sep not installed. pip install sep", file=sys.stderr)
        sys.exit(1)

    from astropy.io import fits

    # Load FITS
    with fits.open(args.fits) as hdul:
        for hdu in hdul:
            if hdu.data is not None and hdu.data.ndim >= 2:
                data = hdu.data.astype(np.float64)
                break
        else:
            print("ERROR: no image HDU found", file=sys.stderr)
            sys.exit(1)

    if data.ndim > 2:
        data = data[0] if data.ndim == 3 else data[0, 0]

    ny, nx = data.shape

    # Convert from 1-indexed to 0-indexed
    cx = args.x - 1.0
    cy = args.y - 1.0

    # Compute cutout region
    radius = args.radius_factor * max(args.a, args.b)
    radius = max(radius, 20.0)  # minimum 20 pixels

    x0 = max(0, int(cx - radius))
    y0 = max(0, int(cy - radius))
    x1 = min(nx, int(cx + radius + 1))
    y1 = min(ny, int(cy + radius + 1))

    if x1 - x0 < 5 or y1 - y0 < 5:
        print("ERROR: cutout too small", file=sys.stderr)
        sys.exit(1)

    cutout = data[y0:y1, x0:x1].copy()

    # Ensure C-contiguous
    if not cutout.flags['C_CONTIGUOUS']:
        cutout = np.ascontiguousarray(cutout)

    # Background subtraction on cutout
    bkg_size = min(args.back_size, min(cutout.shape) // 2)
    if bkg_size < 4:
        bkg_size = 4
    try:
        bkg = sep.Background(cutout, bw=bkg_size, bh=bkg_size)
        cutout_sub = cutout - bkg.back()
        globalrms = bkg.globalrms
    except Exception:
        cutout_sub = cutout - np.median(cutout)
        globalrms = np.std(cutout)

    if globalrms <= 0:
        globalrms = 1.0

    # Run SEP extract with aggressive deblending
    try:
        objects = sep.extract(
            cutout_sub,
            thresh=args.detect_thresh,
            err=globalrms,
            minarea=args.detect_minarea,
            deblend_nthresh=args.deblend_nthresh,
            deblend_cont=args.deblend_mincont,
            segmentation_map=False
        )
    except Exception as e:
        print(f"ERROR: sep.extract failed: {e}", file=sys.stderr)
        sys.exit(1)

    if len(objects) < 2:
        print("#SEPARATE\tN_SUB=0\tMESSAGE=No sub-components found")
        sys.exit(0)

    # Convert cutout coordinates back to full-image 1-indexed coordinates
    # Filter: only keep sub-sources that fall within the parent ellipse region
    parent_cx_cut = cx - x0
    parent_cy_cut = cy - y0
    parent_a = args.a
    parent_b = args.b
    theta_rad = np.radians(args.theta)
    cos_t = np.cos(theta_rad)
    sin_t = np.sin(theta_rad)

    sub_sources = []
    for obj in objects:
        # Check if sub-source center is within expanded parent ellipse
        dx = obj['x'] - parent_cx_cut
        dy = obj['y'] - parent_cy_cut
        # Rotate to ellipse frame
        dx_rot = dx * cos_t + dy * sin_t
        dy_rot = -dx * sin_t + dy * cos_t
        # Test against expanded ellipse (2x parent size to be inclusive)
        ell_dist = (dx_rot / (parent_a * 2.0 + 1))**2 + (dy_rot / (parent_b * 2.0 + 1))**2
        if ell_dist > 1.0:
            continue

        sub_sources.append(obj)

    if len(sub_sources) < 2:
        print("#SEPARATE\tN_SUB=0\tMESSAGE=No sub-components found")
        sys.exit(0)

    # Find the brightest sub-source — it keeps the parent's original shape
    fluxes = np.array([obj['flux'] for obj in sub_sources])
    brightest_idx = np.argmax(fluxes)

    # Output header
    n_sub = len(sub_sources)
    print(f"#SEPARATE\tN_SUB={n_sub}\tPARENT={args.parent_number}")

    cols = ["NUMBER", "X_IMAGE", "Y_IMAGE", "A_IMAGE", "B_IMAGE", "THETA_IMAGE",
            "ISO_RADIUS", "FLUX_AUTO", "MAG_AUTO", "FLAGS"]
    print("\t".join(cols))

    mag_zp = 25.0  # default
    for i, obj in enumerate(sub_sources):
        num = f"{args.parent_number}.{i+1}"
        x_img = obj['x'] + x0 + 1.0  # back to 1-indexed
        y_img = obj['y'] + y0 + 1.0

        if i == brightest_idx:
            # Brightest sub-source: keep parent's original ellipse shape
            a_img = args.a
            b_img = args.b
            theta_deg = args.theta
            iso_r = args.iso_radius if args.iso_radius > 0 else max(a_img, 2.0)
        else:
            a_img = obj['a']
            b_img = obj['b']
            theta_deg = np.degrees(obj['theta'])
            iso_r = max(a_img, 2.0)

        flux = obj['flux']
        if flux > 0:
            mag = -2.5 * np.log10(flux) + mag_zp
        else:
            mag = 99.0
        flags = int(obj['flag'])

        vals = [num,
                f"{x_img:.3f}", f"{y_img:.3f}",
                f"{a_img:.3f}", f"{b_img:.3f}",
                f"{theta_deg:.2f}", f"{iso_r:.3f}",
                f"{flux:.4f}", f"{mag:.4f}",
                str(flags)]
        print("\t".join(vals))


if __name__ == "__main__":
    main()
