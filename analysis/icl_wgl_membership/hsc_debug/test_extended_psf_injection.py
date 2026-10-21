import numpy as np
import pytest

from extended_psf_injection import (
    calculate_transfer_response, measure_point_source_wing_profile,
    satellite_halos_truncated, select_reference_frames,
    unit_point_source_surface_brightness,
)


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


def test_extended_reference_selection_keeps_i_and_i2_support_separate():
    pairs = ([(f"i-{index}", "i") for index in range(13)]
             + [(f"i2-{index}", "i2") for index in range(44)])
    supported = {(ident, filt) for ident, filt in pairs}

    i_rows = select_reference_frames(pairs, supported, "i")
    i2_rows = select_reference_frames(pairs, supported, "i2")

    assert len(i_rows) == 13 and all(filt == "i" for _, filt in i_rows)
    assert len(i2_rows) == 44 and all(filt == "i2" for _, filt in i2_rows)
    with pytest.raises(ValueError, match="expected 44 i2-band frames"):
        select_reference_frames(pairs[:13], supported, "i2")
    with pytest.raises(ValueError, match="must be i or i2"):
        select_reference_frames(pairs, supported, "z")


def test_transfer_response_calculation_matches_scenario_differences():
    measurements = {
        "baseline": {"icl_raw_200_300_njy_arcsec2": 10.,
                     "satellite_delta_50_75_njy_arcsec2": 20.},
        "icl_only": {"icl_raw_200_300_njy_arcsec2": 13.,
                     "satellite_delta_50_75_njy_arcsec2": 21.},
        "satellite_only": {"icl_raw_200_300_njy_arcsec2": 18.,
                           "satellite_delta_50_75_njy_arcsec2": 25.},
        "both": {"icl_raw_200_300_njy_arcsec2": 22.,
                 "satellite_delta_50_75_njy_arcsec2": 29.},
    }
    result = calculate_transfer_response(measurements, 1., 2.)
    assert result["icl_response_to_diffuse_unit_injection"] == pytest.approx(3.)
    assert result["icl_response_to_all_satellite_envelopes"] == pytest.approx(4.)
    assert result["satellite_response_to_smooth_diffuse"] == pytest.approx(1.)
    assert result["satellite_response_to_satellite_unit_injection"] == pytest.approx(2.5)
    assert result["both_minus_sum_residual_icl_njy_arcsec2"] == pytest.approx(1.)
    assert result["both_minus_sum_residual_satellite_njy_arcsec2"] == pytest.approx(3.)


def test_transfer_response_rejects_nonpositive_amplitudes():
    with pytest.raises(ValueError, match="finite and positive"):
        calculate_transfer_response({}, 0., 1.)


def test_unit_point_source_model_conserves_fractional_pixel_flux():
    pixel_scale = 0.168
    source = unit_point_source_surface_brightness((24, 24), (10.25, 11.5), pixel_scale)
    assert source.sum() * pixel_scale**2 == pytest.approx(1.0, abs=1e-7)
    assert np.count_nonzero(source) == 4


def test_point_source_wing_profile_reports_mask_coverage_and_annular_flux():
    size = 701
    yy, xx = np.indices((size, size))
    radius = np.hypot(xx - 350, yy - 350)
    image = np.full((size, size), 2.0, dtype=float)
    excluded = np.zeros((size, size), dtype=bool)
    excluded[350, 410:510] = True
    variance = np.ones((size, size), dtype=float)
    result = measure_point_source_wing_profile(
        image, excluded, variance, radius, 1.0, 0.168, 1.0)

    inner = result["annuli"]["50_100_kpc"]
    outer = result["annuli"]["200_300_kpc"]
    assert 0 < inner["mask_coverage"] < 1
    assert 0.999 < outer["mask_coverage"] <= 1.0
    assert inner["mean_surface_brightness_njy_arcsec2_per_njy"] == pytest.approx(2.0)
    assert outer["mean_surface_brightness_njy_arcsec2_per_njy"] == pytest.approx(2.0)
