#!/usr/bin/env python3
"""Completeness simulation: inject artificial sources and measure recovery rate."""

import sys
import os
import argparse
import numpy as np

# Add project root to path
_script_dir = os.path.dirname(os.path.abspath(__file__))
_project_root = os.path.abspath(os.path.join(_script_dir, '..', '..'))
if _project_root not in sys.path:
    sys.path.insert(0, _project_root)


def _worker_completeness_bin(args):
    """Worker for parallel_map: run completeness for one magnitude bin.

    Each worker makes its own data.copy(), loads shared arrays for
    mask/segmap, and uses a deterministic seed for reproducibility.
    """
    (data_spec, mask_spec, segmap_spec, ibin, mag_lo, mag_hi, mag_c,
     n_inject, detect_thresh, detect_minarea, mag_zeropoint,
     psf, psf_half, psf_size, match_radius, seed) = args

    shm_data = None
    shm_mask = None
    shm_segmap = None
    try:
        import sep

        data, shm_data = data_spec.attach()
        mask, shm_mask = mask_spec.attach()
        segmap, shm_segmap = segmap_spec.attach()

        ny, nx = data.shape
        rng = np.random.RandomState(seed)

        # Inject sources into a copy
        injected = []
        sim_data = data.copy()
        n_placed = 0

        for _ in range(n_inject * 10):
            if n_placed >= n_inject:
                break
            ix = rng.randint(psf_half + 10, nx - psf_half - 10)
            iy = rng.randint(psf_half + 10, ny - psf_half - 10)

            if mask[iy, ix]:
                continue
            if segmap[iy, ix] > 0:
                continue

            mag = rng.uniform(mag_lo, mag_hi)
            flux = 10.0**((mag_zeropoint - mag) / 2.5)

            y0 = iy - psf_half
            x0 = ix - psf_half
            sim_data[y0:y0+psf_size, x0:x0+psf_size] += flux * psf
            injected.append((ix, iy, mag))
            n_placed += 1

        if n_placed == 0:
            return {'mag_c': mag_c, 'n_placed': 0, 'n_recovered': 0,
                    'completeness': 0.0}

        # Re-extract
        sim_mask = ~np.isfinite(sim_data) | (sim_data == 0)
        sim_data[~np.isfinite(sim_data)] = 0.0
        sim_bkg = sep.Background(sim_data, mask=sim_mask.astype(np.uint8),
                                  bw=64, bh=64)
        sim_sub = sim_data - sim_bkg
        sim_rms = sim_bkg.rms()

        rec_objects = sep.extract(sim_sub, detect_thresh, err=sim_rms,
                                   mask=sim_mask.astype(np.uint8),
                                   minarea=detect_minarea)
        rec_x = rec_objects['x']
        rec_y = rec_objects['y']

        # Match
        n_recovered = 0
        match_r2 = match_radius**2
        for inj_x, inj_y, _ in injected:
            dx = rec_x - inj_x
            dy = rec_y - inj_y
            d2 = dx**2 + dy**2
            if len(d2) > 0 and np.min(d2) <= match_r2:
                n_recovered += 1

        completeness = n_recovered / n_placed if n_placed > 0 else 0.0
        return {'mag_c': mag_c, 'n_placed': n_placed,
                'n_recovered': n_recovered, 'completeness': completeness}

    except Exception as e:
        return {'mag_c': mag_c if 'mag_c' in dir() else 0.0,
                'n_placed': 0, 'n_recovered': 0, 'completeness': 0.0}
    finally:
        if shm_data is not None:
            shm_data.close()
        if shm_mask is not None:
            shm_mask.close()
        if shm_segmap is not None:
            shm_segmap.close()


def main():
    parser = argparse.ArgumentParser(description='Completeness simulation')
    parser.add_argument('fitsfile', nargs='?', help='Input FITS image')
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

    from parallel import add_batch_args
    add_batch_args(parser)

    args = parser.parse_args()

    if not args.fitsfile:
        print("ERROR: fitsfile required", file=sys.stderr)
        sys.exit(1)

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

    from parallel import parallel_map, resolve_n_workers, SharedArray

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

    # Segmentation map to avoid existing sources
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

    actual = resolve_n_workers(args.n_workers)

    # Memory guard: each worker copies the full image
    image_bytes = data.nbytes
    try:
        import psutil
        avail_mem = psutil.virtual_memory().available
    except ImportError:
        avail_mem = 64 * 1024**3  # assume 64GB
    max_copies = max(1, int(avail_mem * 0.7 / image_bytes))
    effective_workers = min(actual, args.n_bins, max_copies)

    base_seed = 42

    print("MAG_BIN\tN_INJECTED\tN_RECOVERED\tCOMPLETENESS")

    if effective_workers == 1:
        # Sequential fallback
        rng = np.random.RandomState(base_seed)
        for ibin in range(args.n_bins):
            mag_lo = mag_bins[ibin]
            mag_hi = mag_bins[ibin + 1]
            mag_c = mag_centers[ibin]

            injected = []
            sim_data = data.copy()
            n_placed = 0

            for _ in range(args.n_inject * 10):
                if n_placed >= args.n_inject:
                    break
                ix = rng.randint(psf_half + 10, nx - psf_half - 10)
                iy = rng.randint(psf_half + 10, ny - psf_half - 10)
                if mask[iy, ix]:
                    continue
                if segmap[iy, ix] > 0:
                    continue
                mag = rng.uniform(mag_lo, mag_hi)
                flux = 10.0**((args.mag_zeropoint - mag) / 2.5)
                y0 = iy - psf_half
                x0 = ix - psf_half
                sim_data[y0:y0+psf_size, x0:x0+psf_size] += flux * psf
                injected.append((ix, iy, mag))
                n_placed += 1

            if n_placed == 0:
                print(f"{mag_c:.2f}\t0\t0\t0.000")
                continue

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
            print(f"  Mag {mag_c:.1f}: {n_recovered}/{n_placed} = {completeness:.1%}",
                  file=sys.stderr)
    else:
        # Parallel: share data, mask, segmap
        # Convert mask to uint8 for shared memory
        mask_u8 = mask.astype(np.uint8)
        segmap_i32 = segmap.astype(np.int32)

        with SharedArray(data) as data_spec, \
             SharedArray(mask_u8) as mask_spec, \
             SharedArray(segmap_i32) as segmap_spec:

            tasks = []
            for ibin in range(args.n_bins):
                tasks.append((
                    data_spec, mask_spec, segmap_spec,
                    ibin, mag_bins[ibin], mag_bins[ibin + 1], mag_centers[ibin],
                    args.n_inject, args.detect_thresh, args.detect_minarea,
                    args.mag_zeropoint, psf, psf_half, psf_size,
                    args.match_radius, base_seed + ibin,
                ))

            results = parallel_map(_worker_completeness_bin, tasks,
                                    n_workers=effective_workers,
                                    label="Completeness")

        for r in results:
            print(f"{r['mag_c']:.2f}\t{r['n_placed']}\t{r['n_recovered']}\t"
                  f"{r['completeness']:.4f}")

    print("Done", file=sys.stderr)


if __name__ == '__main__':
    main()
