"""Tests for sersic_fit.fast — C library accuracy and performance.

Usage:
    cd OGFinder && python3 -m sersic_fit.test_fast
"""

import sys
import time
import numpy as np
from scipy.optimize import least_squares

from sersic_fit.utils import sersic_bn, ellipse_radius
from sersic_fit.models import sersic_2d


def _make_grid(size):
    """Create flattened coordinate grids."""
    yy, xx = np.mgrid[0:size, 0:size]
    xx_flat = np.ascontiguousarray(xx.ravel(), dtype=np.float64)
    yy_flat = np.ascontiguousarray(yy.ravel(), dtype=np.float64)
    return xx, yy, xx_flat, yy_flat


def test_eval_accuracy():
    """Compare C sersic_2d_eval vs Python sersic_2d."""
    from sersic_fit.fast import sersic_2d_fast

    print("=== Eval accuracy test ===")
    rng = np.random.RandomState(42)

    test_cases = [
        # Standard case
        [100.0, 10.0, 1.5, 50.0, 50.0, 0.3, 0.5, 1.0],
        # Low n
        [50.0, 5.0, 0.3, 25.0, 25.0, 0.0, 0.0, 0.5],
        # High n, high ellip
        [200.0, 20.0, 6.0, 50.0, 50.0, 0.85, -1.0, 2.0],
        # Small re
        [500.0, 0.5, 1.0, 50.0, 50.0, 0.1, 1.5, 0.1],
        # Large re
        [10.0, 100.0, 4.0, 50.0, 50.0, 0.5, -0.3, 0.0],
    ]

    # Add random cases
    for _ in range(10):
        test_cases.append([
            rng.uniform(1, 500),       # Ie
            rng.uniform(0.5, 50),      # re
            rng.uniform(0.2, 8.0),     # n
            rng.uniform(20, 80),       # xc
            rng.uniform(20, 80),       # yc
            rng.uniform(0, 0.9),       # ellip
            rng.uniform(-np.pi, np.pi),# theta
            rng.uniform(-1, 5),        # bg
        ])

    max_err = 0.0
    for params in test_cases:
        xx, yy, xx_flat, yy_flat = _make_grid(100)
        py_result = sersic_2d(params, xx, yy).ravel()
        c_result = sersic_2d_fast(params, xx_flat, yy_flat)

        diff = np.max(np.abs(py_result - c_result))
        max_err = max(max_err, diff)
        if diff > 1e-8:
            print(f"  WARN: diff={diff:.2e} for params={params}")

    print(f"  Max absolute error: {max_err:.2e}")
    assert max_err < 1e-8, f"Eval accuracy failed: max_err={max_err}"
    print("  PASS")


def test_residual_accuracy():
    """Compare C residual vs Python residual."""
    from sersic_fit.fast import sersic_2d_residual_fast

    print("=== Residual accuracy test ===")
    rng = np.random.RandomState(123)

    params = [100.0, 15.0, 2.0, 50.0, 50.0, 0.3, 0.4, 1.5]
    xx, yy, xx_flat, yy_flat = _make_grid(100)

    # Synthetic data with noise
    data = sersic_2d(params, xx, yy) + rng.normal(0, 0.5, (100, 100))
    data_flat = np.ascontiguousarray(data.ravel(), dtype=np.float64)

    # Test with perturbed params
    params_test = [95.0, 14.0, 2.2, 50.5, 49.5, 0.35, 0.45, 1.3]
    py_resid = (data - sersic_2d(params_test, xx, yy)).ravel()
    c_resid = sersic_2d_residual_fast(params_test, data_flat, xx_flat, yy_flat)

    max_err = np.max(np.abs(py_resid - c_resid))
    print(f"  Max absolute error: {max_err:.2e}")
    assert max_err < 1e-8, f"Residual accuracy failed: max_err={max_err}"
    print("  PASS")


def test_jacobian_accuracy():
    """Compare C analytical Jacobian vs numerical finite differences."""
    from sersic_fit.fast import sersic_2d_fast, sersic_2d_jacobian_fast

    print("=== Jacobian accuracy test ===")
    param_names = ['Ie', 're', 'n', 'xc', 'yc', 'ellip', 'theta', 'bg']

    test_cases = [
        [100.0, 10.0, 1.5, 25.0, 25.0, 0.3, 0.5, 1.0],
        [50.0, 5.0, 0.5, 25.0, 25.0, 0.0, 0.0, 0.5],
        [200.0, 20.0, 4.0, 25.0, 25.0, 0.8, -1.0, 2.0],
        [500.0, 2.0, 1.0, 25.0, 25.0, 0.1, 1.5, 0.1],
    ]

    max_rel_err = 0.0

    for params in test_cases:
        xx, yy, xx_flat, yy_flat = _make_grid(50)
        npix = len(xx_flat)

        # Analytical Jacobian (of residual = -df/d*)
        jac_c = sersic_2d_jacobian_fast(params, xx_flat, yy_flat)

        # Numerical Jacobian via finite differences
        eps = 1e-7
        jac_num = np.zeros((npix, 8), dtype=np.float64)
        p = np.array(params, dtype=np.float64)

        for j in range(8):
            p_plus = p.copy()
            p_minus = p.copy()
            h = max(eps * abs(p[j]), eps)
            p_plus[j] += h
            p_minus[j] -= h

            f_plus = sersic_2d_fast(p_plus, xx_flat, yy_flat)
            f_minus = sersic_2d_fast(p_minus, xx_flat, yy_flat)
            # d(residual)/d(param) = d(data-f)/d(param) = -df/d(param)
            jac_num[:, j] = -(f_plus - f_minus) / (2 * h)

        # Relative error: |a-b| / max(|b|, col_scale * 1e-8)
        # col_scale avoids spurious failures when derivative is ~0
        for j in range(8):
            col_c = jac_c[:, j]
            col_n = jac_num[:, j]
            col_scale = np.max(np.abs(col_n))
            if col_scale < 1e-6:
                # Both analytical and numerical are effectively zero
                abs_err = np.max(np.abs(col_c - col_n))
                if abs_err > 1e-4:
                    print(f"  WARN: param={param_names[j]} abs_err={abs_err:.2e} "
                          f"(near-zero column) for params={params}")
                    max_rel_err = max(max_rel_err, abs_err)
                continue
            scale = np.maximum(np.abs(col_n), col_scale * 1e-8)
            rel_err = np.max(np.abs(col_c - col_n) / scale)
            max_rel_err = max(max_rel_err, rel_err)
            if rel_err > 1e-3:
                print(f"  WARN: param={param_names[j]} rel_err={rel_err:.2e} "
                      f"for params={params}")

    print(f"  Max relative error: {max_rel_err:.2e}")
    assert max_rel_err < 1e-3, f"Jacobian accuracy failed: max_rel_err={max_rel_err}"
    print("  PASS")


def test_edge_cases():
    """Test edge cases: center pixel, extreme parameters."""
    from sersic_fit.fast import sersic_2d_fast, sersic_2d_jacobian_fast

    print("=== Edge case tests ===")

    # Center pixel (r=0)
    params = [100.0, 10.0, 1.5, 5.0, 5.0, 0.3, 0.5, 1.0]
    xx_flat = np.array([5.0], dtype=np.float64)
    yy_flat = np.array([5.0], dtype=np.float64)
    out = sersic_2d_fast(params, xx_flat, yy_flat)
    jac = sersic_2d_jacobian_fast(params, xx_flat, yy_flat)
    assert np.all(np.isfinite(out)), "Non-finite at center pixel"
    assert np.all(np.isfinite(jac)), "Non-finite Jacobian at center pixel"
    print("  Center pixel: OK")

    # Minimum n
    params = [100.0, 10.0, 0.2, 25.0, 25.0, 0.0, 0.0, 0.0]
    xx, yy, xx_flat, yy_flat = _make_grid(50)
    out = sersic_2d_fast(params, xx_flat, yy_flat)
    jac = sersic_2d_jacobian_fast(params, xx_flat, yy_flat)
    assert np.all(np.isfinite(out)), "Non-finite with n=0.2"
    assert np.all(np.isfinite(jac)), "Non-finite Jacobian with n=0.2"
    print("  n=0.2: OK")

    # High ellipticity
    params = [100.0, 10.0, 1.0, 25.0, 25.0, 0.9, 0.0, 0.0]
    out = sersic_2d_fast(params, xx_flat, yy_flat)
    jac = sersic_2d_jacobian_fast(params, xx_flat, yy_flat)
    assert np.all(np.isfinite(out)), "Non-finite with ellip=0.9"
    assert np.all(np.isfinite(jac)), "Non-finite Jacobian with ellip=0.9"
    print("  ellip=0.9: OK")

    # Small re
    params = [500.0, 0.5, 1.0, 25.0, 25.0, 0.1, 0.0, 0.0]
    out = sersic_2d_fast(params, xx_flat, yy_flat)
    jac = sersic_2d_jacobian_fast(params, xx_flat, yy_flat)
    assert np.all(np.isfinite(out)), "Non-finite with re=0.5"
    assert np.all(np.isfinite(jac)), "Non-finite Jacobian with re=0.5"
    print("  re=0.5: OK")

    # Various grid sizes
    for sz in [5, 10, 50, 100, 200]:
        params = [100.0, 10.0, 1.5, sz/2, sz/2, 0.3, 0.5, 1.0]
        xx, yy, xx_flat, yy_flat = _make_grid(sz)
        out = sersic_2d_fast(params, xx_flat, yy_flat)
        jac = sersic_2d_jacobian_fast(params, xx_flat, yy_flat)
        assert np.all(np.isfinite(out)), f"Non-finite at {sz}x{sz}"
        assert np.all(np.isfinite(jac)), f"Non-finite Jacobian at {sz}x{sz}"
    print("  Grid sizes 5-200: OK")

    print("  PASS")


def test_fitting_convergence():
    """Test that C LM fitter converges to same result as Python least_squares."""
    from sersic_fit.fast import fit_sersic_2d_c

    print("=== Fitting convergence test ===")
    rng = np.random.RandomState(99)

    true_params = [150.0, 12.0, 2.5, 50.0, 50.0, 0.3, 0.4, 1.0]
    size = 100
    xx, yy, xx_flat, yy_flat = _make_grid(size)
    data = sersic_2d(true_params, xx, yy) + rng.normal(0, 0.3, (size, size))
    data_flat = np.ascontiguousarray(data.ravel(), dtype=np.float64)

    # Perturbed initial guess
    p0 = [120.0, 10.0, 2.0, 51.0, 49.0, 0.25, 0.3, 1.2]
    bounds_lo = [0, 0.5, 0.2, 40, 40, 0.0, -np.pi, -np.inf]
    bounds_hi = [np.inf, 60, 10, 60, 60, 0.95, np.pi, np.inf]

    # Python least_squares path
    def py_resid(p):
        return (data - sersic_2d(p, xx, yy)).ravel()

    result_py = least_squares(py_resid, p0, bounds=(bounds_lo, bounds_hi),
                              method='trf', max_nfev=500)

    # C LM fitter (entire optimization in C)
    # Replace inf bounds with large finite values for C
    blo = np.array(bounds_lo, dtype=np.float64)
    bhi = np.array(bounds_hi, dtype=np.float64)
    blo[blo == -np.inf] = -1e30
    bhi[bhi == np.inf] = 1e30

    params_c, chi2_c = fit_sersic_2d_c(data_flat, xx_flat, yy_flat, p0,
                                         blo, bhi, max_iter=500)

    # Compare fitted parameters
    param_names = ['Ie', 're', 'n', 'xc', 'yc', 'ellip', 'theta', 'bg']
    print(f"  {'Param':>8s}  {'True':>10s}  {'Python':>10s}  {'C LM':>10s}")
    for i, name in enumerate(param_names):
        print(f"  {name:>8s}  {true_params[i]:10.4f}  "
              f"{result_py.x[i]:10.4f}  {params_c[i]:10.4f}")

    chi2_py = np.sum(result_py.fun**2) / max(1, len(result_py.fun) - 8)
    print(f"  chi2: Python={chi2_py:.6f}  C_LM={chi2_c:.6f}")

    # Both should converge close to true parameters
    for i in range(8):
        if abs(true_params[i]) > 0.01:
            rel = abs(params_c[i] - true_params[i]) / abs(true_params[i])
            assert rel < 0.3, (f"C LM diverged for {param_names[i]}: "
                               f"rel_err={rel:.3f}")

    # chi2 should be reasonable
    assert chi2_c < 1.0, f"C LM chi2 too high: {chi2_c:.4f}"

    print("  PASS")


def benchmark():
    """Benchmark Python least_squares vs C LM fitter at various cutout sizes."""
    from sersic_fit.fast import fit_sersic_2d_c

    print("\n=== Performance benchmark ===")
    rng = np.random.RandomState(42)

    for size in [50, 100, 200, 400]:
        xx, yy, xx_flat, yy_flat = _make_grid(size)
        ctr = size / 2
        params = [100.0, 15.0, 2.0, ctr, ctr, 0.3, 0.5, 1.0]
        data = sersic_2d(params, xx, yy) + rng.normal(0, 0.5, (size, size))
        data_flat = np.ascontiguousarray(data.ravel(), dtype=np.float64)
        p0 = [90.0, 13.0, 1.8, ctr + 0.5, ctr - 0.5, 0.25, 0.45, 1.2]
        bounds_lo = [0, 0.5, 0.2, ctr - 10, ctr - 10, 0.0, -np.pi, -np.inf]
        bounds_hi = [np.inf, 60, 10, ctr + 10, ctr + 10, 0.95, np.pi, np.inf]

        # Replace inf bounds for C
        blo = np.array(bounds_lo, dtype=np.float64)
        bhi = np.array(bounds_hi, dtype=np.float64)
        blo[blo == -np.inf] = -1e30
        bhi[bhi == np.inf] = 1e30

        n_iter = max(5, 50 // (size // 50))

        # Python least_squares (baseline)
        def py_resid(p):
            return (data - sersic_2d(p, xx, yy)).ravel()

        t0 = time.perf_counter()
        for _ in range(n_iter):
            least_squares(py_resid, p0, bounds=(bounds_lo, bounds_hi),
                          method='trf', max_nfev=500)
        t_py = (time.perf_counter() - t0) / n_iter

        # C LM fitter (entire optimization in C)
        t0 = time.perf_counter()
        for _ in range(n_iter):
            fit_sersic_2d_c(data_flat, xx_flat, yy_flat, p0,
                            blo, bhi, max_iter=500)
        t_c = (time.perf_counter() - t0) / n_iter

        print(f"  {size}x{size}: Python={t_py*1000:8.1f}ms  "
              f"C(LM)={t_c*1000:8.1f}ms  speedup={t_py/t_c:.1f}x")


def main():
    from sersic_fit.fast import is_available
    if not is_available():
        print("ERROR: libsersic_fast.so not found. Run build_sersic_fast.sh first.",
              file=sys.stderr)
        sys.exit(1)

    print(f"libsersic_fast.so loaded successfully\n")

    test_eval_accuracy()
    test_residual_accuracy()
    test_jacobian_accuracy()
    test_edge_cases()
    test_fitting_convergence()
    benchmark()

    print("\n=== All tests passed ===")


if __name__ == '__main__':
    main()
