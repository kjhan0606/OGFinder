"""Source masking pipeline for ICL detection.

Masks detected sources so that ICL signal is not contaminated
by galaxy/star light during background modeling and SB profile
measurement.
"""

import sys
import numpy as np
from collections import defaultdict
from scipy.ndimage import binary_dilation, distance_transform_edt


def create_source_mask(data, segmap, expand_factor=2.5):
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
    radius_groups = defaultdict(list)
    for lab in labels:
        npix = label_counts[lab] if lab < len(label_counts) else 0
        if npix == 0:
            continue
        r_eq = np.sqrt(npix / np.pi)
        r_dilate = max(1, int(r_eq * expand_factor))
        radius_groups[r_dilate].append(lab)

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

    # Dilate each group
    n_groups = len(radius_groups)
    for gi, (r, labs) in enumerate(sorted(radius_groups.items())):
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


def interpolate_masked(data, mask, method='linear', block_size=2048,
                       overlap=64):
    """Interpolate masked regions using surrounding unmasked pixels.

    For large images, uses block-based processing with progress logging.

    Parameters
    ----------
    data : 2D array
        Input image.
    mask : 2D bool array
        True where pixels should be interpolated.
    method : str
        Interpolation method: 'linear', 'cubic', 'nearest'.
    block_size : int
        Block size for large-image processing.
    overlap : int
        Overlap between blocks.

    Returns
    -------
    result : 2D array
        Image with masked regions filled by interpolation.
    """
    ny, nx = data.shape
    result = data.copy()

    n_masked = int(mask.sum())
    if n_masked == 0:
        return result

    frac = n_masked / mask.size
    print(f"  interpolation: {n_masked} pixels ({100*frac:.1f}%)",
          file=sys.stderr)

    # Small image: process directly
    if ny <= block_size and nx <= block_size:
        _interpolate_block(result, mask.copy(), method)
        return result

    # Count total blocks
    step = block_size - overlap
    blocks = []
    for y0 in range(0, ny, step):
        y1 = min(ny, y0 + block_size)
        for x0 in range(0, nx, step):
            x1 = min(nx, x0 + block_size)
            if np.any(mask[y0:y1, x0:x1]):
                blocks.append((y0, y1, x0, x1))

    n_blocks = len(blocks)
    print(f"  interpolation: {n_blocks} blocks to process", file=sys.stderr)

    for bi, (y0, y1, x0, x1) in enumerate(blocks):
        block_mask = mask[y0:y1, x0:x1].copy()
        block_data = result[y0:y1, x0:x1].copy()
        _interpolate_block(block_data, block_mask, method)

        # Write only the non-overlap interior
        iy0 = 0 if y0 == 0 else overlap // 2
        iy1 = block_data.shape[0] if y1 == ny else block_data.shape[0] - overlap // 2
        ix0 = 0 if x0 == 0 else overlap // 2
        ix1 = block_data.shape[1] if x1 == nx else block_data.shape[1] - overlap // 2

        result[y0 + iy0:y0 + iy1, x0 + ix0:x0 + ix1] = \
            block_data[iy0:iy1, ix0:ix1]

        if (bi + 1) % 5 == 0 or bi + 1 == n_blocks:
            print(f"  interpolation: block {bi+1}/{n_blocks}",
                  file=sys.stderr)

    return result


def _interpolate_block(data, mask, method):
    """Interpolate masked pixels in a single block in-place.

    Uses two-pass Gaussian-weighted fill: one pass with sigma
    proportional to gap size fills most pixels, a second pass
    with doubled sigma catches the rest.
    """
    from scipy.ndimage import gaussian_filter

    n_bad = int(mask.sum())
    if n_bad == 0:
        return

    n_good = mask.size - n_bad
    if n_good < 10:
        med = np.median(data[~mask]) if n_good > 0 else 0.0
        data[mask] = med
        return

    # If mask fraction > 90%, just fill with local median
    if n_bad > 0.9 * mask.size:
        data[mask] = np.median(data[~mask])
        return

    weight = (~mask).astype(np.float32)
    tmp = data.astype(np.float32)
    tmp[mask] = 0.0

    # Use large sigma to cover gaps in fewer passes
    sigma = max(30.0, np.sqrt(n_bad / np.pi) * 0.7)

    for _ in range(5):
        if not np.any(mask):
            return

        sd = gaussian_filter(tmp, sigma=sigma)
        sw = gaussian_filter(weight, sigma=sigma)

        fillable = mask & (sw > 0.003)
        n_fill = int(fillable.sum())
        if n_fill > 0:
            ratio = sd[fillable] / sw[fillable]
            data[fillable] = ratio
            tmp[fillable] = ratio
            weight[fillable] = 1.0
            mask &= ~fillable
        sigma *= 2.5

    # Fill remaining with median
    if np.any(mask):
        data[mask] = np.median(data[weight > 0])
