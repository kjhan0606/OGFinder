import numpy as np

from audit_satellite_icl_shared_pixel_covariance import (
    diagonal_covariance,
    lag_weighted_covariance,
    supports_outer_annulus,
)


def test_per_bin_support_ignores_failures_in_other_annuli():
    profile = {"coverage": np.array([0.2, 0.4, 0.6, 0.59, 0.61, 0.7, 0.8]),
               "passes_coverage": False}

    assert supports_outer_annulus(profile)
    assert not supports_outer_annulus({"coverage": np.array([1, 1, 1, 1, 1, 1, 0.59])})
    assert not supports_outer_annulus(None)


def test_diagonal_covariance_includes_signed_target_control_weights():
    variance = np.array([4.0, 9.0, 16.0])
    weights = {0: 0.5, 1: -0.25, 2: 0.25}

    result = diagonal_covariance(weights, np.array([0, 1]), variance, conversion=2.0)

    assert result["n_satellite_weight_pixels"] == 3
    assert result["n_shared_nonzero_weight_pixels"] == 2
    assert result["diagonal_var_icl_target_njy2_arcsec4"] == 13.0
    assert result["diagonal_var_satellite_unweighted_njy2_arcsec4"] == 10.25
    assert result["diagonal_covariance_njy2_arcsec4"] == -0.5
    assert result["diagonal_pixel_noise_correlation"] < 0


def test_diagonal_covariance_is_zero_for_disjoint_support():
    result = diagonal_covariance({2: 1.0}, np.array([0, 1]),
                                 np.array([1.0, 1.0, 3.0]), conversion=1.0)

    assert result["n_shared_nonzero_weight_pixels"] == 0
    assert result["diagonal_covariance_njy2_arcsec4"] == 0.0
    assert result["diagonal_pixel_noise_correlation"] == 0.0


def test_lag_weighted_covariance_matches_direct_signed_pixel_sum():
    variance = np.arange(1, 50, dtype=float).reshape(7, 7)
    wa = np.zeros((7, 7)); wb = np.zeros((7, 7))
    wa[2, 2] = 0.5; wa[3, 4] = 0.5
    wb[2, 3] = 1.0; wb[4, 4] = -0.5
    kernel = np.array([[0.05, 0.1, 0.05],
                       [0.1, 1.0, 0.1],
                       [0.05, 0.1, 0.05]])
    expected = 0.0
    for ay, ax in zip(*np.nonzero(wa)):
        for by, bx in zip(*np.nonzero(wb)):
            dy, dx = by-ay, bx-ax
            if abs(dy) <= 1 and abs(dx) <= 1:
                expected += (wa[ay, ax] * wb[by, bx] * np.sqrt(
                    variance[ay, ax] * variance[by, bx]) * kernel[dy+1, dx+1])

    actual = lag_weighted_covariance(wa, wb, variance, kernel, conversion=2.0)

    assert np.isclose(actual, expected * 4.0, rtol=1e-12, atol=1e-12)


def test_zero_lag_kernel_reproduces_diagonal_variance():
    variance = np.arange(1, 50, dtype=float).reshape(7, 7)
    weights = np.zeros((7, 7)); weights[2, 2] = 0.25; weights[3, 4] = -0.75
    kernel = np.zeros((3, 3)); kernel[1, 1] = 1.0
    actual = lag_weighted_covariance(weights, weights, variance, kernel, conversion=3.0)
    expected = np.sum(weights**2 * variance) * 9.0
    assert np.isclose(actual, expected, rtol=1e-12, atol=1e-12)
