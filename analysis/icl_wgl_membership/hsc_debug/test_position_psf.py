import numpy as np
import pytest

from position_psf import convolve_surface_brightness, satellite_halos_with_psfs


def test_signed_centered_psf_convolution_preserves_flux_and_negative_lobes():
    image = np.zeros((21, 21), dtype=np.float32)
    image[10, 10] = 7.0
    psf = np.array([[0.0, -0.01, 0.0],
                    [-0.01, 1.04, -0.01],
                    [0.0, -0.01, 0.0]], dtype=np.float32)

    result = convolve_surface_brightness(image, psf)

    assert result.sum() == pytest.approx(7.0, rel=1e-6)
    assert result[10, 9] < 0
    assert result[10, 10] == pytest.approx(7.28)


def test_position_specific_satellite_psfs_are_applied_one_per_center():
    psf_a = np.ones((1, 1), dtype=np.float32)
    psf_b = np.array([[0.0, 0.0, 0.0],
                      [0.0, 0.8, 0.2],
                      [0.0, 0.0, 0.0]], dtype=np.float32)
    centers = np.array([[40.0, 40.0], [160.0, 160.0]])
    result = satellite_halos_with_psfs((201, 201), centers, 1.0, 1.0,
                                       15.0, [psf_a, psf_b])

    assert result[40, 40] == pytest.approx(np.exp(50.0 / 15.0), rel=2e-5)
    assert result[160, 161] != pytest.approx(result[160, 159], rel=1e-3)
    assert np.isfinite(result).all()


def test_position_specific_satellite_psf_count_must_match_centers():
    with pytest.raises(ValueError, match="one local PSF"):
        satellite_halos_with_psfs((32, 32), [[16.0, 16.0]], 1.0, 1.0,
                                  15.0, [])
