import numpy as np
import pytest

from propagate_satellite_icl_measurement_noise import (
    covariance_cholesky,
    partial_pearson,
    partial_pearson_batch,
    simulate_noise_r,
)


def rows(n=12):
    z = np.linspace(0.31, 0.54, n)
    richness = np.linspace(20, 70, n) ** 1.1
    x = np.sin(np.arange(n) * 0.7) + 0.2 * z
    y = 0.55 * x + 0.03 * richness / np.max(richness) + np.cos(np.arange(n) * 0.41)
    association = [{"icl_delta_200_300": str(a), "satellite_delta_50_75": str(b),
                    "redshift": str(c), "richness": str(d)}
                   for a, b, c, d in zip(x, y, z, richness)]
    covariance = [{"icl_delta_diagonal_variance": "0.5",
                   "target_satellite_diagonal_variance_njy2_arcsec4": "0.8",
                   "target_icl_satellite_diagonal_covariance_njy2_arcsec4": "0.05",
                   "icl_delta_lag_variance": "0.7",
                   "target_satellite_lag_variance_njy2_arcsec4": "1.0",
                   "target_icl_satellite_lag_covariance_njy2_arcsec4": "0.08"}
                  for _ in range(n)]
    return association, covariance


def test_batch_partial_pearson_matches_scalar():
    association, _ = rows()
    x = np.array([float(row["icl_delta_200_300"]) for row in association])
    y = np.array([float(row["satellite_delta_50_75"]) for row in association])
    z = np.array([float(row["redshift"]) for row in association])
    rich = np.array([float(row["richness"]) for row in association])
    batch = partial_pearson_batch(np.stack((x, x + 0.1)), np.stack((y, y - 0.1)), z, rich)
    assert batch[0] == pytest.approx(partial_pearson(x, y, z, rich))
    assert batch[1] == pytest.approx(partial_pearson(x + 0.1, y - 0.1, z, rich))


def test_covariance_cholesky_rejects_non_psd_matrix():
    _, covariance = rows(1)
    covariance[0]["target_icl_satellite_lag_covariance_njy2_arcsec4"] = "2"
    with pytest.raises(ValueError, match="positive semidefinite"):
        covariance_cholesky(covariance, "finite_lag")


def test_noise_monte_carlo_is_reproducible_and_uses_both_models():
    association, covariance = rows()
    diagonal_a = simulate_noise_r(association, covariance, "diagonal", 500, 72)
    diagonal_b = simulate_noise_r(association, covariance, "diagonal", 500, 72)
    lag = simulate_noise_r(association, covariance, "finite_lag", 500, 72)
    np.testing.assert_array_equal(diagonal_a, diagonal_b)
    assert np.isfinite(diagonal_a).all()
    assert np.isfinite(lag).all()
    assert not np.array_equal(diagonal_a, lag)
