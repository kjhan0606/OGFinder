"""FITS cutout extraction and normalization for galaxy morphology."""

import numpy as np


def extract_cutout(data, x, y, size=64, a_image=5.0):
    """Extract a cutout centered on (x, y) with adaptive scaling.

    Larger galaxies get wider cutouts resized to size x size.
    Edge sources are zero-padded.

    Args:
        data: 2D FITS image array
        x, y: center coordinates (0-based)
        size: output cutout size in pixels
        a_image: semi-major axis (pixels) for adaptive scaling

    Returns:
        cutout: (size, size) float32 array, or None if invalid
    """
    h, w = data.shape
    ix, iy = int(round(x)), int(round(y))

    # Adaptive scaling: larger galaxies need wider field of view
    scale = max(1.0, a_image / 5.0)
    half = int(round(size * scale / 2.0))

    # Source region bounds
    x0 = ix - half
    y0 = iy - half
    x1 = ix + half
    y1 = iy + half

    # Clip to image bounds
    cx0 = max(0, x0)
    cy0 = max(0, y0)
    cx1 = min(w, x1)
    cy1 = min(h, y1)

    if cx1 <= cx0 or cy1 <= cy0:
        return None

    # Extract with zero-padding for edge sources
    raw = np.zeros((y1 - y0, x1 - x0), dtype=np.float32)
    raw[cy0 - y0:cy1 - y0, cx0 - x0:cx1 - x0] = data[cy0:cy1, cx0:cx1]

    # Resize to target size if needed
    if raw.shape[0] != size or raw.shape[1] != size:
        raw = _resize_2d(raw, size, size)

    return raw.astype(np.float32)


def _resize_2d(arr, target_h, target_w):
    """Simple bilinear-ish resize using numpy (no PIL dependency)."""
    h, w = arr.shape
    if h == target_h and w == target_w:
        return arr

    # Use numpy fancy indexing for resize
    row_idx = np.linspace(0, h - 1, target_h)
    col_idx = np.linspace(0, w - 1, target_w)

    # Floor indices and fractions
    r0 = np.floor(row_idx).astype(int)
    c0 = np.floor(col_idx).astype(int)
    r1 = np.minimum(r0 + 1, h - 1)
    c1 = np.minimum(c0 + 1, w - 1)
    rf = row_idx - r0
    cf = col_idx - c0

    # Bilinear interpolation
    rf = rf[:, np.newaxis]
    cf = cf[np.newaxis, :]

    result = ((1 - rf) * (1 - cf) * arr[np.ix_(r0, c0)] +
              (1 - rf) * cf * arr[np.ix_(r0, c1)] +
              rf * (1 - cf) * arr[np.ix_(r1, c0)] +
              rf * cf * arr[np.ix_(r1, c1)])

    return result.astype(np.float32)


def normalize_cutout(cutout, method='asinh', asinh_a=0.1):
    """Normalize cutout using asinh stretch.

    Args:
        cutout: 2D float array
        method: normalization method ('asinh')
        asinh_a: softening parameter for asinh

    Returns:
        normalized: 2D float32 array in [0, 1]
    """
    if cutout is None:
        return None

    arr = cutout.astype(np.float64)

    # Median / MAD-based background subtraction
    median = np.median(arr)
    mad = np.median(np.abs(arr - median))
    sigma = mad * 1.4826 if mad > 0 else 1.0

    # Subtract background and scale
    arr = (arr - median) / sigma

    if method == 'asinh':
        arr = np.arcsinh(arr / asinh_a)

    # Min-max normalize to [0, 1]
    vmin, vmax = arr.min(), arr.max()
    if vmax > vmin:
        arr = (arr - vmin) / (vmax - vmin)
    else:
        arr = np.zeros_like(arr)

    return arr.astype(np.float32)


def batch_extract_cutouts(data, catalog, size=64, asinh_a=0.1):
    """Extract and normalize cutouts for all sources in catalog.

    Args:
        data: 2D FITS image array (background-subtracted)
        catalog: dict with 'x', 'y', 'a_image' arrays (0-based coords)
        size: output cutout size
        asinh_a: asinh softening parameter

    Returns:
        cutouts: (N_valid, 1, size, size) float32 array
        valid_mask: boolean array of length N_sources
    """
    n = len(catalog['x'])
    cutouts = []
    valid_mask = np.zeros(n, dtype=bool)

    for i in range(n):
        x = catalog['x'][i]
        y = catalog['y'][i]
        a = catalog.get('a_image', np.full(n, 5.0))[i]

        cut = extract_cutout(data, x, y, size=size, a_image=a)
        if cut is None:
            continue

        cut = normalize_cutout(cut, asinh_a=asinh_a)
        if cut is None:
            continue

        cutouts.append(cut[np.newaxis, :, :])  # (1, H, W) channel dim
        valid_mask[i] = True

    if len(cutouts) == 0:
        return np.zeros((0, 1, size, size), dtype=np.float32), valid_mask

    return np.stack(cutouts, axis=0).astype(np.float32), valid_mask
