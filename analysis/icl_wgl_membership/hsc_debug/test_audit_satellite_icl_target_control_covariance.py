import numpy as np
import pytest

from audit_satellite_icl_target_control_covariance import (
    annulus_weights, combine_variances, diagonal_weighted_variance,
)
from verify_satellite_icl_target_control_covariance import max_lag_distance_arcsec


def test_annulus_weights_are_uniform_and_match_valid_coverage():
    valid = np.ones((11, 11), dtype=bool)
    valid[5, 6] = False
    weights, coverage, count = annulus_weights(
        valid.shape, (5, 5), 1.0, valid, inner_kpc=1.0, outer_kpc=4.0)
    yy, xx = np.indices(valid.shape)
    radius = np.hypot(xx - 5, yy - 5)
    geometric = (radius >= 1) & (radius < 4)
    assert count == int(geometric.sum()) - 1
    assert np.isclose(weights.sum(), 1.0)
    assert np.allclose(weights[weights > 0], 1.0 / count)
    assert np.isclose(coverage, count / geometric.sum())


def test_target_control_variances_add_and_cross_covariance_is_retained():
    result = combine_variances(2.0, 3.0, 5.0, 7.0,
                               11.0, 13.0, 0.5, 0.75)
    assert result["icl_delta_diagonal_variance"] == 7.0
    assert result["icl_delta_lag_variance"] == 10.0
    assert np.isclose(result["icl_delta_sigma_inflation"], np.sqrt(10 / 7))
    assert result["control_to_target_lag_variance_ratio"] == 7 / 3
    assert np.isclose(result["icl_satellite_lag_noise_correlation"],
                      0.75 / np.sqrt(10 * 13))


def test_combine_variances_rejects_nonpositive_variance():
    with pytest.raises(ValueError):
        combine_variances(2, 3, 5, 0, 11, 13, 0.5, 0.75)


def test_diagonal_variance_ignores_nonfinite_unweighted_pixels():
    weights = np.array([[0.5, 0.5], [0.0, 0.0]])
    variance = np.array([[4.0, 8.0], [np.inf, np.nan]])
    assert np.isclose(diagonal_weighted_variance(weights, variance, 2.0), 12.0)


def test_diagonal_variance_rejects_nonfinite_weighted_pixels():
    weights = np.array([[0.5, 0.5]])
    variance = np.array([[4.0, np.inf]])
    with pytest.raises(ValueError):
        diagonal_weighted_variance(weights, variance, 1.0)


def test_finite_lag_kernel_support_uses_corner_pixel_separation():
    assert max_lag_distance_arcsec((9, 9), 0.168) == pytest.approx(
        np.sqrt(32) * 0.168)
    with pytest.raises(ValueError):
        max_lag_distance_arcsec((8, 9), 0.168)
