"""ctypes wrapper for libsersic_fast.so — fast Sérsic model + analytical Jacobian.

Auto-loads the shared library at import time.  Falls back gracefully:
    from sersic_fit.fast import is_available
    if is_available():
        ...  # use fast path
"""

import os
import ctypes
import numpy as np

_lib = None
_FAST_AVAILABLE = False


def _load_library():
    """Try to load libsersic_fast.so from ../lib/ relative to this file."""
    global _lib, _FAST_AVAILABLE

    here = os.path.dirname(os.path.abspath(__file__))
    so_path = os.path.join(here, '..', 'lib', 'libsersic_fast.so')
    so_path = os.path.normpath(so_path)

    if not os.path.isfile(so_path):
        return

    try:
        lib = ctypes.CDLL(so_path)
    except OSError:
        return

    # void sersic_2d_eval(const double *params, const double *x, const double *y,
    #                     double *out, int npix)
    lib.sersic_2d_eval.restype = None
    lib.sersic_2d_eval.argtypes = [
        ctypes.POINTER(ctypes.c_double),  # params
        ctypes.POINTER(ctypes.c_double),  # x
        ctypes.POINTER(ctypes.c_double),  # y
        ctypes.POINTER(ctypes.c_double),  # out
        ctypes.c_int,                     # npix
    ]

    # void sersic_2d_residual(const double *params, const double *data,
    #                         const double *x, const double *y,
    #                         double *residual, int npix)
    lib.sersic_2d_residual.restype = None
    lib.sersic_2d_residual.argtypes = [
        ctypes.POINTER(ctypes.c_double),  # params
        ctypes.POINTER(ctypes.c_double),  # data
        ctypes.POINTER(ctypes.c_double),  # x
        ctypes.POINTER(ctypes.c_double),  # y
        ctypes.POINTER(ctypes.c_double),  # residual
        ctypes.c_int,                     # npix
    ]

    # void sersic_2d_jacobian(const double *params, const double *x, const double *y,
    #                         double *jac, int npix)
    lib.sersic_2d_jacobian.restype = None
    lib.sersic_2d_jacobian.argtypes = [
        ctypes.POINTER(ctypes.c_double),  # params
        ctypes.POINTER(ctypes.c_double),  # x
        ctypes.POINTER(ctypes.c_double),  # y
        ctypes.POINTER(ctypes.c_double),  # jac
        ctypes.c_int,                     # npix
    ]

    # void sersic_2d_cost_and_grad(const double *params, const double *data,
    #     const double *x, const double *y,
    #     double *cost, double *grad, int npix)
    lib.sersic_2d_cost_and_grad.restype = None
    lib.sersic_2d_cost_and_grad.argtypes = [
        ctypes.POINTER(ctypes.c_double),  # params
        ctypes.POINTER(ctypes.c_double),  # data
        ctypes.POINTER(ctypes.c_double),  # x
        ctypes.POINTER(ctypes.c_double),  # y
        ctypes.POINTER(ctypes.c_double),  # cost (scalar output)
        ctypes.POINTER(ctypes.c_double),  # grad (8 output)
        ctypes.c_int,                     # npix
    ]

    # double sersic_2d_fit(const double *data, const double *x, const double *y,
    #     int npix, const double *params_init,
    #     const double *bounds_lo, const double *bounds_hi,
    #     int max_iter, double *params_out)
    lib.sersic_2d_fit.restype = ctypes.c_double
    lib.sersic_2d_fit.argtypes = [
        ctypes.POINTER(ctypes.c_double),  # data
        ctypes.POINTER(ctypes.c_double),  # x
        ctypes.POINTER(ctypes.c_double),  # y
        ctypes.c_int,                     # npix
        ctypes.POINTER(ctypes.c_double),  # params_init
        ctypes.POINTER(ctypes.c_double),  # bounds_lo
        ctypes.POINTER(ctypes.c_double),  # bounds_hi
        ctypes.c_int,                     # max_iter
        ctypes.POINTER(ctypes.c_double),  # params_out
    ]

    _lib = lib
    _FAST_AVAILABLE = True


# Auto-load at import
_load_library()


def is_available():
    """Return True if the C library was loaded successfully."""
    return _FAST_AVAILABLE


def _as_c_double_ptr(arr):
    """Get ctypes double pointer from a contiguous float64 array."""
    return arr.ctypes.data_as(ctypes.POINTER(ctypes.c_double))


def _ensure_f64_contiguous(arr):
    """Ensure array is contiguous float64."""
    return np.ascontiguousarray(arr, dtype=np.float64)


def sersic_2d_fast(params, x_flat, y_flat):
    """Evaluate 2D Sérsic model using C library.

    Parameters
    ----------
    params : array-like, shape (8,)
        [Ie, re, n, xc, yc, ellip, theta, bg]
    x_flat, y_flat : 1D float64 arrays (flattened pixel coords)

    Returns
    -------
    out : 1D float64 array, shape (npix,)
    """
    params = _ensure_f64_contiguous(np.asarray(params, dtype=np.float64))
    x_flat = _ensure_f64_contiguous(x_flat)
    y_flat = _ensure_f64_contiguous(y_flat)
    npix = len(x_flat)

    out = np.empty(npix, dtype=np.float64)

    _lib.sersic_2d_eval(
        _as_c_double_ptr(params),
        _as_c_double_ptr(x_flat),
        _as_c_double_ptr(y_flat),
        _as_c_double_ptr(out),
        ctypes.c_int(npix),
    )
    return out


def sersic_2d_residual_fast(params, data_flat, x_flat, y_flat):
    """Compute residuals (data - model) using C library.

    Parameters
    ----------
    params : array-like, shape (8,)
    data_flat, x_flat, y_flat : 1D float64 arrays

    Returns
    -------
    residual : 1D float64 array, shape (npix,)
    """
    params = _ensure_f64_contiguous(np.asarray(params, dtype=np.float64))
    data_flat = _ensure_f64_contiguous(data_flat)
    x_flat = _ensure_f64_contiguous(x_flat)
    y_flat = _ensure_f64_contiguous(y_flat)
    npix = len(x_flat)

    residual = np.empty(npix, dtype=np.float64)

    _lib.sersic_2d_residual(
        _as_c_double_ptr(params),
        _as_c_double_ptr(data_flat),
        _as_c_double_ptr(x_flat),
        _as_c_double_ptr(y_flat),
        _as_c_double_ptr(residual),
        ctypes.c_int(npix),
    )
    return residual


def sersic_2d_jacobian_fast(params, x_flat, y_flat):
    """Compute analytical Jacobian d(residual)/d(params) using C library.

    Parameters
    ----------
    params : array-like, shape (8,)
    x_flat, y_flat : 1D float64 arrays

    Returns
    -------
    jac : 2D float64 array, shape (npix, 8)
    """
    params = _ensure_f64_contiguous(np.asarray(params, dtype=np.float64))
    x_flat = _ensure_f64_contiguous(x_flat)
    y_flat = _ensure_f64_contiguous(y_flat)
    npix = len(x_flat)

    jac = np.empty((npix, 8), dtype=np.float64)

    _lib.sersic_2d_jacobian(
        _as_c_double_ptr(params),
        _as_c_double_ptr(x_flat),
        _as_c_double_ptr(y_flat),
        _as_c_double_ptr(jac),
        ctypes.c_int(npix),
    )
    return jac


def make_residual_func(data_flat, x_flat, y_flat):
    """Create an optimized residual function with pre-cached pointers.

    Returns a callable residual(params) -> 1D float64 array that avoids
    per-call array conversion overhead for the fixed arrays (data, x, y).
    """
    # Pre-compute everything that doesn't change per-call
    data_flat = _ensure_f64_contiguous(data_flat)
    x_flat = _ensure_f64_contiguous(x_flat)
    y_flat = _ensure_f64_contiguous(y_flat)
    npix = len(x_flat)
    c_npix = ctypes.c_int(npix)

    data_ptr = _as_c_double_ptr(data_flat)
    x_ptr = _as_c_double_ptr(x_flat)
    y_ptr = _as_c_double_ptr(y_flat)

    # Pre-allocate output buffer and params buffer
    residual_buf = np.empty(npix, dtype=np.float64)
    res_ptr = _as_c_double_ptr(residual_buf)
    params_buf = np.empty(8, dtype=np.float64)
    params_ptr = _as_c_double_ptr(params_buf)

    lib_func = _lib.sersic_2d_residual

    def residual(params):
        params_buf[:] = params
        lib_func(params_ptr, data_ptr, x_ptr, y_ptr, res_ptr, c_npix)
        return residual_buf.copy()  # copy needed: scipy stores references

    return residual


def make_cost_grad_func(data_flat, x_flat, y_flat):
    """Create an optimized cost+gradient function for L-BFGS-B.

    Returns a callable cost_and_grad(params) -> (float, ndarray(8,)) that
    computes sum-of-squared-residuals and its gradient in a single C call.
    Only 9 doubles are returned (1 cost + 8 gradient), eliminating per-call
    array allocation overhead.
    """
    data_flat = _ensure_f64_contiguous(data_flat)
    x_flat = _ensure_f64_contiguous(x_flat)
    y_flat = _ensure_f64_contiguous(y_flat)
    npix = len(x_flat)
    c_npix = ctypes.c_int(npix)

    data_ptr = _as_c_double_ptr(data_flat)
    x_ptr = _as_c_double_ptr(x_flat)
    y_ptr = _as_c_double_ptr(y_flat)

    params_buf = np.empty(8, dtype=np.float64)
    params_ptr = _as_c_double_ptr(params_buf)
    cost_buf = np.empty(1, dtype=np.float64)
    cost_ptr = _as_c_double_ptr(cost_buf)
    grad_buf = np.empty(8, dtype=np.float64)
    grad_ptr = _as_c_double_ptr(grad_buf)

    lib_func = _lib.sersic_2d_cost_and_grad

    def cost_and_grad(params):
        params_buf[:] = params
        lib_func(params_ptr, data_ptr, x_ptr, y_ptr, cost_ptr, grad_ptr, c_npix)
        return float(cost_buf[0]), grad_buf.copy()

    return cost_and_grad


def fit_sersic_2d_c(data_flat, x_flat, y_flat, params_init,
                     bounds_lo, bounds_hi, max_iter=500):
    """Run complete Levenberg-Marquardt fitting in C.

    Parameters
    ----------
    data_flat, x_flat, y_flat : 1D float64 arrays
    params_init : array-like, shape (8,)
    bounds_lo, bounds_hi : array-like, shape (8,)
    max_iter : int

    Returns
    -------
    params_out : 1D float64 array, shape (8,)
    chi2 : float (reduced chi2 = cost / (npix - 8))
    """
    data_flat = _ensure_f64_contiguous(data_flat)
    x_flat = _ensure_f64_contiguous(x_flat)
    y_flat = _ensure_f64_contiguous(y_flat)
    params_init = _ensure_f64_contiguous(np.asarray(params_init, dtype=np.float64))
    bounds_lo = np.asarray(bounds_lo, dtype=np.float64).copy()
    bounds_hi = np.asarray(bounds_hi, dtype=np.float64).copy()
    bounds_lo[bounds_lo == -np.inf] = -1e30
    bounds_hi[bounds_hi == np.inf] = 1e30
    bounds_lo = _ensure_f64_contiguous(bounds_lo)
    bounds_hi = _ensure_f64_contiguous(bounds_hi)
    npix = len(x_flat)

    params_out = np.empty(8, dtype=np.float64)

    chi2 = _lib.sersic_2d_fit(
        _as_c_double_ptr(data_flat),
        _as_c_double_ptr(x_flat),
        _as_c_double_ptr(y_flat),
        ctypes.c_int(npix),
        _as_c_double_ptr(params_init),
        _as_c_double_ptr(bounds_lo),
        _as_c_double_ptr(bounds_hi),
        ctypes.c_int(max_iter),
        _as_c_double_ptr(params_out),
    )

    return params_out, float(chi2)
