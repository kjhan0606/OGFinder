import numpy as np
import pytest

from position_psf import (build_core_wing_hybrid_psf,
                          convolve_surface_brightness,
                          satellite_halos_with_psfs)


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


def test_core_wing_hybrid_is_signed_centered_and_unit_normalized():
    core = np.zeros((5, 5), dtype=float)
    core[2, 2] = 1.0
    core[2, 1] = -0.01
    core[2, 2] += 0.01
    wing = np.zeros((21, 21), dtype=float)
    wing[10, 10] = 1.0
    wing[10, 9] = -0.02
    wing[10, 11] = 0.02
    wing[10, 10] -= 0.01
    wing[10, 12] = 0.01

    hybrid, metadata = build_core_wing_hybrid_psf(
        core, wing, splice_radius_arcsec=0.25,
        blend_width_arcsec=0.1, wing_amplitude=0.7,
        pixel_scale_arcsec=0.168)

    assert hybrid.shape == wing.shape
    assert hybrid[10, 10] > 0.99
    assert hybrid[10, 9] < 0
    assert hybrid.sum() == pytest.approx(1.0, abs=2e-7)
    assert hybrid[10, 12] == pytest.approx(
        0.7 * wing[10, 12] / metadata["hybrid_signed_sum_before_normalization"])
    assert metadata["status"] == "exploratory_core_wing_hybrid_sensitivity_only"
    assert metadata["psf_approved"] is False
    assert metadata["science_measurement"] is False


def test_core_wing_hybrid_centers_rectangular_local_stamp():
    core = np.zeros((5, 7), dtype=float)
    core[2, 3] = 1.01
    core[2, 2] = -0.01
    wing = np.zeros((21, 21), dtype=float)
    wing[10, 10] = 1.0

    hybrid, metadata = build_core_wing_hybrid_psf(
        core, wing, splice_radius_arcsec=0.25,
        blend_width_arcsec=0.1, wing_amplitude=0.7,
        pixel_scale_arcsec=0.168)

    assert metadata["local_core_shape"] == [5, 7]
    assert hybrid[10, 10] == pytest.approx(1.01, abs=2e-6)
    assert hybrid[10, 9] == pytest.approx(-0.01, abs=2e-6)


@pytest.mark.parametrize("params", [
    {"splice_radius_arcsec": -1.0, "blend_width_arcsec": 0.5, "wing_amplitude": 1.0},
    {"splice_radius_arcsec": 1.0, "blend_width_arcsec": 0.0, "wing_amplitude": 1.0},
    {"splice_radius_arcsec": 1.0, "blend_width_arcsec": 0.5, "wing_amplitude": 0.0},
])
def test_core_wing_hybrid_rejects_invalid_sensitivity_grid_values(params):
    core = np.zeros((5, 5), dtype=float)
    core[2, 2] = 1.0
    wing = np.zeros((21, 21), dtype=float)
    wing[10, 10] = 1.0
    with pytest.raises(ValueError):
        build_core_wing_hybrid_psf(core, wing, **params, pixel_scale_arcsec=0.168)
