import numpy as np
import pytest

from extended_psf_injection import satellite_halos_truncated


def test_truncated_satellite_model_preserves_center_and_fifty_kpc_normalization():
    result = satellite_halos_truncated((201, 201), [[100., 100.]], 1.0,
                                       2.5, 15.0)
    assert result[100, 100] == pytest.approx(2.5 * np.exp(50.0 / 15.0), rel=2e-5)
    assert result[100, 150] == pytest.approx(2.5, rel=2e-5)
    assert np.isfinite(result).all()


def test_truncated_satellite_models_add_and_keep_subpixel_centers():
    centers = [[50.0, 50.0], [130.25, 130.0]]
    both = satellite_halos_truncated((181, 181), centers, 1.0, 1.0, 15.0)
    first = satellite_halos_truncated((181, 181), centers[:1], 1.0, 1.0, 15.0)
    second = satellite_halos_truncated((181, 181), centers[1:], 1.0, 1.0, 15.0)
    np.testing.assert_allclose(both, first + second, rtol=1e-6, atol=1e-6)
    assert second[130, 130] != pytest.approx(first[130, 130], rel=1e-2)
