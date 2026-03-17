#!/usr/bin/env python3
"""Detect a single source at a given image position using SEP.

Given a FITS image and (x, y) in 1-indexed image coordinates,
extract a small cutout, run SEP, and output the nearest detected
source as a single TSV row matching the ds9_sextract catalog format.

Usage:
    python3 ds9_add_source.py image.fits --x 150.3 --y 200.1
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
    parser = argparse.ArgumentParser(description="Detect source at position")
    parser.add_argument("fits", help="Input FITS file")
    parser.add_argument("--x", type=float, required=True, help="X position (1-indexed)")
    parser.add_argument("--y", type=float, required=True, help="Y position (1-indexed)")
    parser.add_argument("--radius", type=float, default=25.0,
                        help="Cutout radius in pixels")
    parser.add_argument("--max-dist", type=float, default=10.0,
                        help="Max distance from click to detected center (pixels)")
    parser.add_argument("--detect-thresh", type=float, default=0.8,
                        help="Detection threshold (sigma)")
    parser.add_argument("--detect-minarea", type=int, default=3,
                        help="Minimum detection area (pixels)")
    parser.add_argument("--deblend-nthresh", type=int, default=64,
                        help="Deblending sub-thresholds")
    parser.add_argument("--deblend-mincont", type=float, default=0.0001,
                        help="Deblending min contrast")
    parser.add_argument("--back-size", type=int, default=32,
                        help="Background mesh size")
    parser.add_argument("--number", type=int, default=9999,
                        help="NUMBER to assign to the detected source")
    parser.add_argument("--mag-zeropoint", type=float, default=25.0,
                        help="Magnitude zeropoint")
    parser.add_argument("--phot-aperture", type=float, default=5.0,
                        help="Aperture photometry radius")
    args = parser.parse_args()

    if sep is None:
        print("ERROR: sep not installed. pip install sep", file=sys.stderr)
        sys.exit(1)

    from astropy.io import fits

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

    # Convert 1-indexed to 0-indexed
    cx = args.x - 1.0
    cy = args.y - 1.0

    # Cutout region
    rad = args.radius
    x0 = max(0, int(cx - rad))
    y0 = max(0, int(cy - rad))
    x1 = min(nx, int(cx + rad + 1))
    y1 = min(ny, int(cy + rad + 1))

    if x1 - x0 < 5 or y1 - y0 < 5:
        print("ERROR: cutout too small", file=sys.stderr)
        sys.exit(1)

    cutout = data[y0:y1, x0:x1].copy()
    if not cutout.flags['C_CONTIGUOUS']:
        cutout = np.ascontiguousarray(cutout)

    # Background
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

    # Extract
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

    if len(objects) == 0:
        print("#ADD_SOURCE\tFOUND=0\tMESSAGE=No source detected at this position")
        sys.exit(0)

    # Find the source closest to the target position (in cutout coords)
    target_cx = cx - x0
    target_cy = cy - y0
    dists = np.sqrt((objects['x'] - target_cx)**2 + (objects['y'] - target_cy)**2)
    best_idx = np.argmin(dists)

    # Reject if too far from click position
    if dists[best_idx] > args.max_dist:
        print("#ADD_SOURCE\tFOUND=0\tMESSAGE=No source within "
              f"{args.max_dist:.0f}px of click position")
        sys.exit(0)

    obj = objects[best_idx]

    # Convert back to full-image 1-indexed
    obj_x = obj['x'] + x0 + 1.0
    obj_y = obj['y'] + y0 + 1.0
    obj_a = obj['a']
    obj_b = obj['b']
    obj_theta = np.degrees(obj['theta'])
    iso_r = max(obj_a, 2.0)
    flux_auto = obj['flux']

    # Aperture photometry
    try:
        flux_aper, fluxerr_aper, _ = sep.sum_circle(
            cutout_sub, [obj['x']], [obj['y']], args.phot_aperture,
            err=globalrms, subpix=5)
        flux_aper = flux_aper[0]
        fluxerr_aper = fluxerr_aper[0]
    except Exception:
        flux_aper = flux_auto
        fluxerr_aper = 0.0

    zp = args.mag_zeropoint
    mag_auto = -2.5 * np.log10(flux_auto) + zp if flux_auto > 0 else 99.0
    mag_aper = -2.5 * np.log10(flux_aper) + zp if flux_aper > 0 else 99.0

    # Kron radius and ellipticity
    try:
        kronrad, krflag = sep.kron_radius(
            cutout_sub, [obj['x']], [obj['y']],
            [obj['a']], [obj['b']], [obj['theta']], 6.0)
        kronrad = kronrad[0]
    except Exception:
        kronrad = 2.5

    ellipticity = 1.0 - obj_b / obj_a if obj_a > 0 else 0.0

    # Peak value
    peak = obj['peak']

    # CLASS_STAR placeholder (can't determine without neural net)
    class_star = 0.5

    # FWHM estimate from a and b
    fwhm = 2.0 * np.sqrt(np.log(2.0)) * (obj_a + obj_b)

    # Output
    print(f"#ADD_SOURCE\tFOUND=1\tNUMBER={args.number}")

    # Match ds9_sextract column format
    cols = ["NUMBER", "X_IMAGE", "Y_IMAGE", "A_IMAGE", "B_IMAGE",
            "THETA_IMAGE", "ELLIPTICITY", "KRON_RADIUS", "ISO_RADIUS",
            "FLUX_AUTO", "FLUX_APER", "MAG_AUTO", "MAG_APER",
            "PEAK", "CLASS_STAR", "FLAGS", "FWHM_IMAGE"]
    print("\t".join(cols))

    vals = [str(args.number),
            f"{obj_x:.3f}", f"{obj_y:.3f}",
            f"{obj_a:.3f}", f"{obj_b:.3f}",
            f"{obj_theta:.2f}", f"{ellipticity:.4f}",
            f"{kronrad:.3f}", f"{iso_r:.3f}",
            f"{flux_auto:.4f}", f"{flux_aper:.4f}",
            f"{mag_auto:.4f}", f"{mag_aper:.4f}",
            f"{peak:.4f}", f"{class_star:.2f}",
            str(int(obj['flag'])), f"{fwhm:.3f}"]
    print("\t".join(vals))


if __name__ == "__main__":
    main()
