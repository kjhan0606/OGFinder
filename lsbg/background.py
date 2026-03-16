"""Iterative background refinement for LSBG detection.

Standard background estimation subtracts LSBG signal. This module
uses iterative masking and refitting to preserve diffuse emission
while accurately modeling the sky background.

Parallelized via tiled processing using SharedArray + parallel_map.
"""

import sys
import numpy as np
import sep

from icl.background import (fit_polynomial_background,
                             fit_chebyshev_background,
                             sep_large_mesh_background)
from lsbg.masking import iterative_mask_refine


def _cosine_weight_2d(shape, overlap):
    """Create 2D cosine blending weights for tile overlap regions."""
    ny, nx = shape
    wy = np.ones(ny, dtype=np.float64)
    wx = np.ones(nx, dtype=np.float64)

    if overlap > 0:
        ramp = 0.5 * (1.0 - np.cos(np.pi * np.arange(overlap) / overlap))
        # Fade in at start
        if overlap <= ny:
            wy[:overlap] = ramp
            wy[-overlap:] = ramp[::-1]
        if overlap <= nx:
            wx[:overlap] = ramp
            wx[-overlap:] = ramp[::-1]

    return wy[:, np.newaxis] * wx[np.newaxis, :]


def _fit_background_tile(args):
    """Worker function for parallel tile-based background fitting.

    Parameters
    ----------
    args : tuple
        (data_spec, mask_spec, tile_slice, method, order, sigma_clip, mesh_size)
    """
    from parallel import SharedArraySpec

    data_spec, mask_spec, y0, y1, x0, x1, method, order, sigma_clip, mesh_size = args

    data_arr, data_shm = data_spec.attach()
    mask_arr, mask_shm = mask_spec.attach()

    tile_data = data_arr[y0:y1, x0:x1].copy()
    tile_mask = mask_arr[y0:y1, x0:x1].copy()

    try:
        if method == 'polynomial':
            bkg = fit_polynomial_background(tile_data, tile_mask,
                                            order=order,
                                            sigma_clip=sigma_clip)
        elif method == 'chebyshev':
            bkg = fit_chebyshev_background(tile_data, tile_mask,
                                           order=order,
                                           sigma_clip=sigma_clip)
        else:  # sep_large
            bkg = sep_large_mesh_background(tile_data, tile_mask,
                                            mesh_size=mesh_size)
    finally:
        data_shm.close()
        mask_shm.close()

    return (y0, y1, x0, x1, bkg)


def fit_background_parallel(data, mask, method='sep_large', order=3,
                             sigma_clip=3.0, mesh_size=256,
                             tile_size=1024, overlap=128, n_workers=0):
    """Fit background in parallel using tiled processing.

    Splits the image into overlapping tiles, fits background on each
    tile independently, then blends results with cosine weighting.

    Parameters
    ----------
    data : 2D array
        Input image.
    mask : 2D bool array
        Source mask (True = masked).
    method : str
        Background method: 'sep_large', 'polynomial', 'chebyshev'.
    order : int
        Polynomial order (for polynomial/chebyshev).
    sigma_clip : float
        Sigma-clipping threshold.
    mesh_size : int
        SEP mesh size.
    tile_size : int
        Tile size in pixels.
    overlap : int
        Overlap between tiles.
    n_workers : int
        Number of parallel workers (0=auto, 1=sequential).

    Returns
    -------
    background : 2D array
        Fitted background model.
    """
    from parallel import SharedArray, parallel_map, resolve_n_workers

    ny, nx = data.shape

    # For small images, process directly without tiling
    if ny <= tile_size and nx <= tile_size:
        return _fit_single(data, mask, method, order, sigma_clip, mesh_size)

    actual_workers = resolve_n_workers(n_workers)

    # Build tile list
    tiles = []
    step = tile_size - overlap
    for y0 in range(0, ny, step):
        y1 = min(ny, y0 + tile_size)
        for x0 in range(0, nx, step):
            x1 = min(nx, x0 + tile_size)
            tiles.append((y0, y1, x0, x1))
        if y1 == ny:
            break

    # Sequential fallback
    if actual_workers <= 1 or len(tiles) <= 1:
        return _fit_single(data, mask, method, order, sigma_clip, mesh_size)

    print(f"LSBG background: {len(tiles)} tiles, {actual_workers} workers",
          file=sys.stderr)

    data_c = np.ascontiguousarray(data, dtype=np.float64)
    mask_c = mask.astype(np.uint8)

    with SharedArray(data_c) as data_spec, SharedArray(mask_c) as mask_spec:
        tasks = [
            (data_spec, mask_spec, y0, y1, x0, x1,
             method, order, sigma_clip, mesh_size)
            for y0, y1, x0, x1 in tiles
        ]
        results = parallel_map(_fit_background_tile, tasks,
                               n_workers=n_workers,
                               label="Background tiles")

    # Blend tiles with cosine weighting
    background = np.zeros_like(data, dtype=np.float64)
    weights = np.zeros_like(data, dtype=np.float64)

    for y0, y1, x0, x1, tile_bkg in results:
        tile_h = y1 - y0
        tile_w = x1 - x0
        w = _cosine_weight_2d((tile_h, tile_w), overlap)
        background[y0:y1, x0:x1] += tile_bkg * w
        weights[y0:y1, x0:x1] += w

    # Normalize
    good = weights > 0
    background[good] /= weights[good]
    background[~good] = np.median(data[~mask]) if (~mask).any() else 0.0

    return background


def _fit_single(data, mask, method, order, sigma_clip, mesh_size):
    """Fit background on entire image (no tiling)."""
    if method == 'polynomial':
        return fit_polynomial_background(data, mask, order=order,
                                          sigma_clip=sigma_clip)
    elif method == 'chebyshev':
        return fit_chebyshev_background(data, mask, order=order,
                                         sigma_clip=sigma_clip)
    else:
        return sep_large_mesh_background(data, mask, mesh_size=mesh_size)


def iterative_background(data, initial_mask, config):
    """Iterative background refinement.

    Core loop:
    1. Fit background on masked image
    2. Detect residual sources in background-subtracted image
    3. Expand mask to include new detections
    4. Re-interpolate and repeat

    Parameters
    ----------
    data : 2D array
        Original image.
    initial_mask : 2D bool array
        Initial source mask.
    config : LSBGConfig
        Pipeline configuration.

    Returns
    -------
    background : 2D array
        Final refined background model.
    final_mask : 2D bool array
        Final expanded mask.
    cleaned : 2D array
        Background-subtracted image.
    """
    from icl.masking import interpolate_masked

    mask = initial_mask.copy()

    for iteration in range(config.bkg_n_iterations):
        print(f"LSBG background: iteration {iteration + 1}/"
              f"{config.bkg_n_iterations}", file=sys.stderr)

        # Interpolate masked regions before fitting
        filled = interpolate_masked(data, mask, method=config.interp_method)

        # Fit background
        background = fit_background_parallel(
            filled, mask,
            method=config.bkg_method,
            order=config.bkg_poly_order,
            sigma_clip=config.bkg_sigma_clip,
            mesh_size=config.bkg_mesh_size,
            n_workers=config.n_workers
        )

        # Detect residual sources and expand mask
        mask, n_new = iterative_mask_refine(
            data, mask, background,
            sigma_thresh=config.bkg_refine_thresh,
            minarea=config.mask_detect_minarea
        )
        print(f"  -> {n_new} new residual sources detected",
              file=sys.stderr)

        if n_new == 0:
            print(f"  -> converged at iteration {iteration + 1}",
                  file=sys.stderr)
            break

    cleaned = data - background
    return background, mask, cleaned
