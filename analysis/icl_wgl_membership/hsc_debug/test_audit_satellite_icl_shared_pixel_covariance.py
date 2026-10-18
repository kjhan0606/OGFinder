import numpy as np

from audit_satellite_icl_shared_pixel_covariance import diagonal_covariance, supports_outer_annulus


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
