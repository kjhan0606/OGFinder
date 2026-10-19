import numpy as np

from estimate_hsc_coadd_noise_correlation import lag_correlations


def test_lag_correlation_recovers_correlated_noise_and_unit_diagonal():
    rng = np.random.default_rng(20261019)
    white = rng.normal(size=(128, 128))
    image = (4 * white + np.roll(white, 1, axis=1) + np.roll(white, -1, axis=1)) / np.sqrt(18)
    variance = np.ones_like(image)
    valid = np.ones_like(image, dtype=bool)
    result = lag_correlations(image, variance, valid, max_lag=2)
    assert abs(result[(0, 0)][0] - 1.0) < 0.10
    assert result[(0, 1)][0] > 0.15
    assert result[(0, -1)][0] > 0.15


def test_lag_correlation_ignores_masked_extreme_pixels():
    rng = np.random.default_rng(22)
    image = rng.normal(size=(64, 64))
    variance = np.ones_like(image)
    valid = np.ones_like(image, dtype=bool)
    valid[:8] = False
    image[:8] = 1e6
    result = lag_correlations(image, variance, valid, max_lag=1)
    assert abs(result[(0, 0)][0] - 1.0) < 0.15
    assert abs(result[(0, 1)][0]) < 0.15

