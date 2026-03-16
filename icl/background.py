"""Large-scale background modeling for ICL detection.

Standard SExtractor background estimation uses small mesh sizes
(~64 px) that subtract away the ICL signal.  These methods use
much larger scales or global polynomial fits to preserve ICL.
"""

import numpy as np


def fit_polynomial_background(data, mask, order=3, sigma_clip=3.0):
    """Fit a 2D polynomial to unmasked pixels.

    Downsamples to ~100x100 grid for efficiency, fits a 2D polynomial
    of the given order, then evaluates at full resolution.

    Parameters
    ----------
    data : 2D array
        Input image.
    mask : 2D bool array
        True = masked (source) pixels to exclude from fit.
    order : int
        Polynomial order (e.g., 3 = up to x^3*y^0, x^2*y^1, ...).
    sigma_clip : float
        Sigma-clipping iterations on fit residuals.

    Returns
    -------
    background : 2D array
        Fitted background model at full resolution.
    """
    ny, nx = data.shape
    unmask = ~mask

    # Downsample to ~100x100 grid
    step_y = max(1, ny // 100)
    step_x = max(1, nx // 100)
    yy_grid = np.arange(0, ny, step_y)
    xx_grid = np.arange(0, nx, step_x)
    xx_mg, yy_mg = np.meshgrid(xx_grid, yy_grid)

    # Sample unmasked pixels on grid
    sample_mask = unmask[yy_mg, xx_mg]
    sx = xx_mg[sample_mask].astype(float)
    sy = yy_mg[sample_mask].astype(float)
    sz = data[yy_mg[sample_mask], xx_mg[sample_mask]]

    if len(sz) < (order + 1) * (order + 2) // 2:
        # Not enough points: return median
        med = np.median(data[unmask]) if unmask.any() else 0.0
        return np.full_like(data, med)

    # Normalize coordinates to [-1, 1]
    sx_n = 2.0 * sx / max(nx - 1, 1) - 1.0
    sy_n = 2.0 * sy / max(ny - 1, 1) - 1.0

    # Build Vandermonde matrix
    def build_vander(xn, yn, order):
        cols = []
        for i in range(order + 1):
            for j in range(order + 1 - i):
                cols.append(xn ** i * yn ** j)
        return np.column_stack(cols)

    # Sigma-clipping fit
    good = np.ones(len(sz), dtype=bool)
    coeffs = None
    for _ in range(3):
        V = build_vander(sx_n[good], sy_n[good], order)
        coeffs, _, _, _ = np.linalg.lstsq(V, sz[good], rcond=None)
        resid = sz - build_vander(sx_n, sy_n, order) @ coeffs
        sigma = np.std(resid[good])
        if sigma <= 0:
            break
        good = np.abs(resid) < sigma_clip * sigma

    # Evaluate at full resolution
    all_y, all_x = np.mgrid[0:ny, 0:nx]
    all_xn = 2.0 * all_x.astype(float) / max(nx - 1, 1) - 1.0
    all_yn = 2.0 * all_y.astype(float) / max(ny - 1, 1) - 1.0

    # Evaluate row by row to save memory
    background = np.empty_like(data, dtype=float)
    for row in range(ny):
        xn_row = all_xn[row]
        yn_row = all_yn[row]
        V_row = build_vander(xn_row, yn_row, order)
        background[row] = V_row @ coeffs

    return background


def fit_chebyshev_background(data, mask, order=3, sigma_clip=3.0):
    """Fit a 2D Chebyshev polynomial to unmasked pixels.

    Uses numpy.polynomial.chebyshev for numerical stability
    at higher orders.

    Parameters
    ----------
    data : 2D array
        Input image.
    mask : 2D bool array
        True = masked pixels.
    order : int
        Chebyshev polynomial order.
    sigma_clip : float
        Sigma-clipping threshold.

    Returns
    -------
    background : 2D array
        Fitted background model.
    """
    from numpy.polynomial import chebyshev

    ny, nx = data.shape
    unmask = ~mask

    # Downsample
    step_y = max(1, ny // 100)
    step_x = max(1, nx // 100)
    yy_grid = np.arange(0, ny, step_y)
    xx_grid = np.arange(0, nx, step_x)
    xx_mg, yy_mg = np.meshgrid(xx_grid, yy_grid)

    sample_mask = unmask[yy_mg, xx_mg]
    sx = xx_mg[sample_mask].astype(float)
    sy = yy_mg[sample_mask].astype(float)
    sz = data[yy_mg[sample_mask], xx_mg[sample_mask]]

    n_terms = (order + 1) * (order + 2) // 2
    if len(sz) < n_terms:
        med = np.median(data[unmask]) if unmask.any() else 0.0
        return np.full_like(data, med)

    # Normalize to [-1, 1]
    sx_n = 2.0 * sx / max(nx - 1, 1) - 1.0
    sy_n = 2.0 * sy / max(ny - 1, 1) - 1.0

    # Build Chebyshev Vandermonde
    def build_cheb_vander(xn, yn, order):
        cols = []
        for i in range(order + 1):
            Ti = chebyshev.chebval(xn, [0] * i + [1])
            for j in range(order + 1 - i):
                Tj = chebyshev.chebval(yn, [0] * j + [1])
                cols.append(Ti * Tj)
        return np.column_stack(cols)

    # Sigma-clipping fit
    good = np.ones(len(sz), dtype=bool)
    coeffs = None
    for _ in range(3):
        V = build_cheb_vander(sx_n[good], sy_n[good], order)
        coeffs, _, _, _ = np.linalg.lstsq(V, sz[good], rcond=None)
        resid = sz - build_cheb_vander(sx_n, sy_n, order) @ coeffs
        sigma = np.std(resid[good])
        if sigma <= 0:
            break
        good = np.abs(resid) < sigma_clip * sigma

    # Evaluate row by row
    all_y, all_x = np.mgrid[0:ny, 0:nx]
    all_xn = 2.0 * all_x.astype(float) / max(nx - 1, 1) - 1.0
    all_yn = 2.0 * all_y.astype(float) / max(ny - 1, 1) - 1.0

    background = np.empty_like(data, dtype=float)
    for row in range(ny):
        V_row = build_cheb_vander(all_xn[row], all_yn[row], order)
        background[row] = V_row @ coeffs

    return background


def sep_large_mesh_background(data, mask, mesh_size=256):
    """Estimate background using SEP with large mesh size.

    Simplest method: uses sep.Background with a very large mesh
    to avoid subtracting ICL signal.

    Parameters
    ----------
    data : 2D array
        Input image.
    mask : 2D bool array
        True = masked pixels.
    mesh_size : int
        Mesh size in pixels (default 256).

    Returns
    -------
    background : 2D array
        SEP background model.
    """
    import sep

    data_c = np.ascontiguousarray(data, dtype=np.float64)
    mask_byte = mask.astype(np.uint8) if mask is not None else None

    bkg = sep.Background(data_c, mask=mask_byte,
                         bw=mesh_size, bh=mesh_size)
    return bkg.back()
