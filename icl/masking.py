"""Source masking pipeline for ICL detection.

Masks detected sources so that ICL signal is not contaminated
by galaxy/star light during background modeling and SB profile
measurement.
"""

import sys
import numpy as np
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
from scipy.ndimage import binary_dilation, distance_transform_edt


def _dilate_group_worker(args):
    """Worker for parallel dilation of one radius group."""
    segmap_spec, labs, r = args
    from parallel.shared_array import SharedArraySpec
    segmap, shm = segmap_spec.attach()
    try:
        submask = np.isin(segmap, labs)
        if r > 30:
            dist = distance_transform_edt(~submask)
            result = (dist <= r)
        else:
            y, x = np.ogrid[-r:r + 1, -r:r + 1]
            struct = (x * x + y * y) <= r * r
            result = binary_dilation(submask, structure=struct)
        return result
    finally:
        shm.close()


def create_source_mask(data, segmap, expand_factor=2.5, max_dilate_radius=50,
                       n_workers=0):
    """Expand segmentation map segments to create a source mask.

    Segments are grouped by dilation radius and batch-dilated for
    efficiency.  Large radii use distance_transform_edt (O(N) in
    image size regardless of radius) instead of binary_dilation.

    Parameters
    ----------
    data : 2D array
        Image data (used only for shape).
    segmap : 2D int array
        Segmentation map from SEP (0 = background).
    expand_factor : float
        Dilation factor relative to segment equivalent radius.
    max_dilate_radius : int
        Maximum dilation radius in pixels.
    n_workers : int
        0 = auto parallel, 1 = sequential, N = N workers.

    Returns
    -------
    mask : 2D bool array
        True where sources are masked.
    """
    mask = np.zeros(data.shape, dtype=bool)
    labels = np.unique(segmap)
    labels = labels[labels > 0]

    if len(labels) == 0:
        return mask

    # O(N) pixel count per label via bincount
    label_counts = np.bincount(segmap.ravel())

    # Compute dilation radii and group by radius
    n_capped = 0
    radius_groups = defaultdict(list)
    for lab in labels:
        npix = label_counts[lab] if lab < len(label_counts) else 0
        if npix == 0:
            continue
        r_eq = np.sqrt(npix / np.pi)
        r_dilate = max(1, int(r_eq * expand_factor))
        if r_dilate > max_dilate_radius:
            n_capped += 1
            r_dilate = max_dilate_radius
        radius_groups[r_dilate].append(lab)

    if n_capped > 0:
        print(f"  dilation radius capped to {max_dilate_radius}px "
              f"for {n_capped} sources", file=sys.stderr)

    # Merge into at most ~15 groups
    sorted_radii = sorted(radius_groups.keys())
    if len(sorted_radii) > 15:
        all_radii = []
        for r, labs in radius_groups.items():
            all_radii.extend([r] * len(labs))
        all_radii = np.array(all_radii)
        bin_edges = np.percentile(all_radii, np.linspace(0, 100, 16))
        bin_edges = np.unique(np.round(bin_edges).astype(int))
        merged = defaultdict(list)
        merged_r = {}
        for r in sorted_radii:
            b = int(np.searchsorted(bin_edges, r, side='right'))
            merged[b].extend(radius_groups[r])
            merged_r[b] = max(merged_r.get(b, 0), r)
        radius_groups = {merged_r[b]: labs for b, labs in merged.items()}

    # Dilate groups — parallel if multiple groups and workers available
    sorted_groups = sorted(radius_groups.items())
    n_groups = len(sorted_groups)

    from parallel.pool import resolve_n_workers
    actual_workers = resolve_n_workers(n_workers)

    if n_groups >= 2 and actual_workers > 1:
        from parallel.shared_array import SharedArray
        print(f"  dilation: {n_groups} groups, {actual_workers} workers",
              file=sys.stderr)
        with SharedArray(segmap) as seg_spec:
            tasks = [(seg_spec, labs, r) for r, labs in sorted_groups]
            from parallel.pool import parallel_map
            results = parallel_map(
                _dilate_group_worker, tasks,
                n_workers=min(actual_workers, n_groups),
                label="  dilation"
            )
        for submask in results:
            mask |= submask
    else:
        # Sequential fallback
        for gi, (r, labs) in enumerate(sorted_groups):
            submask = np.isin(segmap, labs)
            if r > 30:
                dist = distance_transform_edt(~submask)
                mask |= (dist <= r)
            else:
                y, x = np.ogrid[-r:r + 1, -r:r + 1]
                struct = (x * x + y * y) <= r * r
                mask |= binary_dilation(submask, structure=struct)
            print(f"  dilation group {gi + 1}/{n_groups} "
                  f"(r={r}, {len(labs)} sources)", file=sys.stderr)

    return mask


def mask_bright_stars(mask, catalog, mag_limit=18.0, radius_scale=10.0):
    """Add circular masks for bright stars.

    Bright sources (MAG_AUTO < mag_limit) get circular masks with
    radius proportional to their brightness excess.

    Parameters
    ----------
    mask : 2D bool array
        Existing mask (modified in-place).
    catalog : dict
        Must contain 'x_image', 'y_image', 'mag_auto' arrays.
        Coordinates are 1-based (DS9 convention).
    mag_limit : float
        Magnitude threshold; sources brighter get extra masking.
    radius_scale : float
        Base radius scale factor.
    """
    if catalog is None:
        return

    mag = catalog.get('mag_auto')
    if mag is None:
        return

    x_img = catalog['x_image']
    y_img = catalog['y_image']

    ny, nx = mask.shape

    bright = mag < mag_limit
    if not np.any(bright):
        return

    max_radius = min(ny, nx) // 4
    for i in np.where(bright)[0]:
        # Radius scales with brightness excess
        r = radius_scale * 10 ** (-0.2 * (mag[i] - mag_limit))
        r = max(3, min(int(r + 0.5), max_radius))

        # 0-based center
        cx = int(x_img[i] - 1 + 0.5)
        cy = int(y_img[i] - 1 + 0.5)

        # Bounding box
        y0 = max(0, cy - r)
        y1 = min(ny, cy + r + 1)
        x0 = max(0, cx - r)
        x1 = min(nx, cx + r + 1)

        yy, xx = np.ogrid[y0:y1, x0:x1]
        dist2 = (xx - cx) ** 2 + (yy - cy) ** 2
        mask[y0:y1, x0:x1] |= dist2 <= r * r


def reproject_mask(mask_data, mask_header, target_header):
    """Reproject boolean mask from one WCS frame to another.

    Uses astropy.wcs pixel-to-sky transforms + nearest-neighbor sampling.
    Fast path: same shape + same WCS -> direct copy.

    Parameters
    ----------
    mask_data : 2D array
        Source mask (nonzero = masked).
    mask_header : FITS header
        WCS header of mask_data.
    target_header : FITS header
        WCS header of target image.

    Returns
    -------
    mask_out : 2D bool array
        Reprojected mask matching target image dimensions.
    """
    from astropy.wcs import WCS

    # Target shape from header
    tny = int(target_header.get('NAXIS2', 0))
    tnx = int(target_header.get('NAXIS1', 0))
    if tny == 0 or tnx == 0:
        tny, tnx = mask_data.shape

    mny, mnx = mask_data.shape

    # Check if WCS headers match (fast path)
    same_shape = (mny == tny and mnx == tnx)
    wcs_keys = ['CRPIX1', 'CRPIX2', 'CRVAL1', 'CRVAL2',
                'CD1_1', 'CD1_2', 'CD2_1', 'CD2_2']

    has_wcs_mask = any(k in mask_header for k in ['CRVAL1', 'CD1_1', 'CDELT1'])
    has_wcs_target = any(k in target_header for k in ['CRVAL1', 'CD1_1', 'CDELT1'])

    if has_wcs_mask and has_wcs_target:
        same_wcs = same_shape
        if same_wcs:
            for k in wcs_keys:
                v1 = mask_header.get(k, None)
                v2 = target_header.get(k, None)
                if v1 is None and v2 is None:
                    continue
                if v1 is None or v2 is None:
                    same_wcs = False
                    break
                if abs(float(v1) - float(v2)) > 1e-10:
                    same_wcs = False
                    break

        if same_wcs:
            print("reproject_mask: same WCS, direct copy", file=sys.stderr)
            return mask_data.astype(bool)

        # Case B: WCS coordinate transform
        mask_wcs = WCS(mask_header, naxis=2)
        target_wcs = WCS(target_header, naxis=2)

        # Build target pixel grid
        yy, xx = np.mgrid[0:tny, 0:tnx]
        # Target pixels -> sky -> mask pixels
        sky = target_wcs.pixel_to_world(xx.ravel(), yy.ravel())
        mx, my = mask_wcs.world_to_pixel(sky)
        mx = np.round(mx).astype(int)
        my = np.round(my).astype(int)

        # Bounds check
        valid = (mx >= 0) & (mx < mnx) & (my >= 0) & (my < mny)
        mask_bool = mask_data.astype(bool)

        mask_out = np.zeros(tny * tnx, dtype=bool)
        mask_out[valid] = mask_bool[my[valid], mx[valid]]
        mask_out = mask_out.reshape(tny, tnx)

        n_masked = int(mask_out.sum())
        print(f"reproject_mask: ({mny},{mnx}) -> ({tny},{tnx}), "
              f"{n_masked} pixels masked", file=sys.stderr)
        return mask_out

    else:
        # No WCS available
        if same_shape:
            print("reproject_mask: no WCS, same shape, direct copy",
                  file=sys.stderr)
            return mask_data.astype(bool)
        else:
            raise ValueError(
                f"Cannot reproject mask: no WCS headers and shapes differ "
                f"({mny},{mnx}) vs ({tny},{tnx})")


def interpolate_masked(data, mask, method='linear', block_size=2048,
                       overlap=256):
    """Interpolate masked regions using surrounding unmasked pixels.

    Uses Gaussian-weighted fill on the full image first, then
    block-based refinement for remaining gaps.

    Parameters
    ----------
    data : 2D array
        Input image.
    mask : 2D bool array
        True where pixels should be interpolated.
    method : str
        Interpolation method (kept for API compat, Gaussian fill used).
    block_size : int
        Block size for large-image processing.
    overlap : int
        Overlap between blocks (large enough to avoid seams).

    Returns
    -------
    result : 2D array
        Image with masked regions filled by interpolation.
    """
    ny, nx = data.shape
    result = data.astype(np.float32).copy()

    n_masked = int(mask.sum())
    if n_masked == 0:
        return result

    frac = n_masked / mask.size
    print(f"  interpolation: {n_masked} pixels ({100*frac:.1f}%)",
          file=sys.stderr)

    remaining = mask.copy()
    _gaussian_fill(result, remaining)

    n_left = int(remaining.sum())
    if n_left > 0:
        print(f"  interpolation: {n_left} pixels remaining after global fill",
              file=sys.stderr)
        # Fill any remainder with local median
        if n_left < remaining.size:
            result[remaining] = np.median(result[~remaining])

    return result


def _gaussian_fill(data, mask, max_passes=8):
    """Fill masked pixels using iterative Gaussian-weighted averaging.

    Fills from edges inward: each pass uses a Gaussian filter to
    propagate values from unmasked pixels into masked regions.
    Sigma grows each pass to reach deeper into large gaps.

    The two gaussian_filter calls per pass (data + weight) run in
    parallel threads since scipy releases the GIL during filtering.
    """
    from scipy.ndimage import gaussian_filter

    n_bad = int(mask.sum())
    if n_bad == 0:
        return

    n_good = mask.size - n_bad
    if n_good < 10:
        med = np.median(data[~mask]) if n_good > 0 else 0.0
        data[mask] = med
        mask[:] = False
        return

    weight = (~mask).astype(np.float32)
    tmp = data.copy()
    tmp[mask] = 0.0

    # Start sigma proportional to typical gap size, minimum 20
    sigma = max(20.0, np.sqrt(n_bad / np.pi) * 0.5)

    for i in range(max_passes):
        if not np.any(mask):
            return

        # Run the two gaussian_filter calls in parallel threads
        # (scipy releases GIL during the C computation)
        with ThreadPoolExecutor(max_workers=2) as executor:
            f_data = executor.submit(gaussian_filter, tmp, sigma=sigma)
            f_weight = executor.submit(gaussian_filter, weight, sigma=sigma)
            sd = f_data.result()
            sw = f_weight.result()

        # Fill where enough weight has diffused in
        fillable = mask & (sw > 0.002)
        n_fill = int(fillable.sum())
        if n_fill > 0:
            ratio = sd[fillable] / sw[fillable]
            data[fillable] = ratio
            tmp[fillable] = ratio
            weight[fillable] = 1.0
            mask &= ~fillable

        n_left = int(mask.sum())
        print(f"  interpolation pass {i+1}: sigma={sigma:.0f}, "
              f"filled={n_fill}, remaining={n_left}", file=sys.stderr)

        if n_left == 0:
            return

        sigma *= 2.0
