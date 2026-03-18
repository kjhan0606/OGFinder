"""Source masking pipeline for ICL detection.

Masks detected sources so that ICL signal is not contaminated
by galaxy/star light during background modeling and SB profile
measurement.
"""

import numpy as np
from scipy.ndimage import binary_dilation, generate_binary_structure


def create_source_mask(data, segmap, expand_factor=2.5):
    """Expand segmentation map segments to create a source mask.

    Each segment is dilated by expand_factor times its equivalent
    circular radius, using per-segment structuring elements.

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

    for lab in labels:
        seg_pix = segmap == lab
        npix = seg_pix.sum()
        if npix == 0:
            continue

        # Equivalent circular radius
        r_eq = np.sqrt(npix / np.pi)
        r_dilate = max(1, int(r_eq * expand_factor))

        # Create circular structuring element
        y, x = np.ogrid[-r_dilate:r_dilate + 1, -r_dilate:r_dilate + 1]
        struct = (x * x + y * y) <= r_dilate * r_dilate

        dilated = binary_dilation(seg_pix, structure=struct)
        mask |= dilated

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


def interpolate_masked(data, mask, method='linear', block_size=4096,
                       overlap=100):
    """Interpolate masked regions using surrounding unmasked pixels.

    For large images, uses block-based processing to limit memory.

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

    n_masked = mask.sum()
    if n_masked == 0:
        return result

    # Small image: process directly
    if ny <= block_size and nx <= block_size:
        _interpolate_block(result, mask, method)
        return result

    # Large image: block-based processing
    for y0 in range(0, ny, block_size - overlap):
        y1 = min(ny, y0 + block_size)
        for x0 in range(0, nx, block_size - overlap):
            x1 = min(nx, x0 + block_size)

            block_mask = mask[y0:y1, x0:x1]
            if not np.any(block_mask):
                continue

            block_data = result[y0:y1, x0:x1].copy()
            _interpolate_block(block_data, block_mask, method)

            # Write only the non-overlap interior
            iy0 = 0 if y0 == 0 else overlap // 2
            iy1 = block_data.shape[0] if y1 == ny else block_data.shape[0] - overlap // 2
            ix0 = 0 if x0 == 0 else overlap // 2
            ix1 = block_data.shape[1] if x1 == nx else block_data.shape[1] - overlap // 2

            gy0 = y0 + iy0
            gy1 = y0 + iy1
            gx0 = x0 + ix0
            gx1 = x0 + ix1

            result[gy0:gy1, gx0:gx1] = block_data[iy0:iy1, ix0:ix1]

    return result


def _interpolate_block(data, mask, method):
    """Interpolate masked pixels in a single block in-place.

    Uses Gaussian-weighted averaging of unmasked pixels for fast,
    smooth interpolation. Falls back to median for edge cases.
    """
    from scipy.ndimage import gaussian_filter

    n_good = (~mask).sum()
    n_bad = mask.sum()

    if n_bad == 0:
        return
    if n_good < 10:
        med = np.median(data[~mask]) if n_good > 0 else 0.0
        data[mask] = med
        return

    # Gaussian-weighted fill: convolve data*(~mask) and (~mask),
    # then divide to get weighted average at masked locations.
    # Choose sigma proportional to typical masked region size.
    sigma = max(5.0, np.sqrt(n_bad / max(1, n_good)) * 10.0)
    sigma = min(sigma, 50.0)

    weight = (~mask).astype(np.float64)
    filled = data.copy().astype(np.float64)
    filled[mask] = 0.0

    smooth_data = gaussian_filter(filled, sigma=sigma)
    smooth_weight = gaussian_filter(weight, sigma=sigma)
    smooth_weight[smooth_weight < 1e-10] = 1.0

    data[mask] = (smooth_data / smooth_weight)[mask]
