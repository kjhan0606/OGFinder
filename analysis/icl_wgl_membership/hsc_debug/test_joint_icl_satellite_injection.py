import numpy as np
import pytest

from joint_icl_satellite_injection import (
    add_surface_brightness_injection,
    redetect_masks_and_targets,
    satellite_outer_measurement,
    satellite_halo_template,
)


def test_satellite_profile_amplitude_is_normalized_at_50_kpc():
    template = satellite_halo_template((101, 101), np.array([[50.0, 50.0]]),
                                      kpc_per_pixel=1.0, amplitude_at_50=3.0,
                                      scale_kpc=15.0)
    assert template[50, 100] == pytest.approx(3.0)
    assert template[50, 50] == pytest.approx(3.0 * np.exp(50.0 / 15.0))
    assert template[50, 1] == pytest.approx(3.0 * np.exp(1.0 / 15.0))


def test_satellite_profiles_add_linearly_and_validate_inputs():
    centers = np.array([[10., 10.], [20., 20.]])
    total = satellite_halo_template((32, 32), centers, 2.0, 1.0)
    first = satellite_halo_template((32, 32), centers[:1], 2.0, 1.0)
    second = satellite_halo_template((32, 32), centers[1:], 2.0, 1.0)
    np.testing.assert_allclose(total, first + second, rtol=1e-6, atol=1e-6)
    with pytest.raises(ValueError):
        satellite_halo_template((32, 32), centers, 0.0, 1.0)
    with pytest.raises(ValueError):
        satellite_halo_template((32, 32), centers, 1.0, -1.0)


def test_surface_brightness_template_is_converted_to_native_units():
    native = np.zeros((4, 4), dtype=np.float32)
    sb = np.full((4, 4), 10.0, dtype=np.float32)
    converted = add_surface_brightness_injection(native, sb, 2000.0)
    np.testing.assert_allclose(converted * 2000.0, sb, rtol=1e-6)


def test_surface_brightness_injection_rejects_bad_units_or_shape():
    native = np.zeros((2, 2))
    with pytest.raises(ValueError):
        add_surface_brightness_injection(native, native, 0.0)
    with pytest.raises(ValueError):
        add_surface_brightness_injection(native, np.zeros((3, 3)), 1.0)


def test_satellite_measurement_records_missing_support_without_losing_icl_frame(monkeypatch):
    import joint_icl_satellite_injection as joint

    monkeypatch.setattr(joint.satellite, "profile_at", lambda *args, **kwargs: None)
    result = satellite_outer_measurement(
        np.zeros((8, 8)), np.ones((8, 8)), np.zeros((8, 8), dtype=bool),
        np.zeros((8, 8), dtype=np.int32), {}, None, None,
        [({}, {"x": 4., "y": 4., "NUMBER": 1}, 100.)],
        np.array([4., 4.]), 0.168, 1.0, 1.0)

    assert result == {
        "delta_sb_njy_arcsec2": None,
        "n_satellites": 0,
        "n_control_profiles": 0,
        "status": "no_satellite_passed_minimum_support",
    }


def test_redetection_rebuilds_masks_and_recovers_injected_member():
    yy, xx = np.mgrid[:128, :128]
    rng = np.random.default_rng(27)
    variance = np.ones((128, 128), dtype=np.float32)
    bits = np.zeros((128, 128), dtype=np.int32)
    image = rng.normal(0.0, 0.2, size=variance.shape).astype(np.float32)
    image += (500 * np.exp(-((xx - 64) ** 2 + (yy - 64) ** 2) / (2 * 2.0 ** 2))).astype(np.float32)
    image += (30 * np.exp(-((xx - 90) ** 2 + (yy - 64) ** 2) / (2 * 2.0 ** 2))).astype(np.float32)
    members = [{"INDEX": "0", "IN_FRAME": "True", "X_IMAGE": "91", "Y_IMAGE": "65"}]

    result = redetect_masks_and_targets(
        image, variance, bits, members, np.array([64.0, 64.0]), 0.168,
    )

    assert result["diagnostics"]["n_detections"] >= 2
    assert result["diagnostics"]["n_matched_members"] == 1
    assert len(result["targets_pixel_radius"]) == 1
    assert result["icl_mask"].shape == image.shape
    assert result["icl_mask"].any()
    assert result["diagnostics"]["n_bcg_candidates"] == 1
    assert result["diagnostics"]["bcg_mask_policy"] == "preserve_unique_bcg"
    assert result["segmentation"].max() <= result["diagnostics"]["n_detections"]


def test_missing_bcg_can_use_explicit_all_source_mask_policy():
    yy, xx = np.mgrid[:128, :128]
    rng = np.random.default_rng(91)
    variance = np.ones((128, 128), dtype=np.float32)
    bits = np.zeros((128, 128), dtype=np.int32)
    image = rng.normal(0.0, 0.2, size=variance.shape).astype(np.float32)
    image += (80 * np.exp(-((xx - 90) ** 2 + (yy - 64) ** 2) / (2 * 2.0 ** 2))).astype(np.float32)
    members = [{"INDEX": "0", "IN_FRAME": "True", "X_IMAGE": "91", "Y_IMAGE": "65"}]

    with pytest.raises(ValueError, match="found 0 BCG candidates"):
        redetect_masks_and_targets(image, variance, bits, members,
                                   np.array([64.0, 64.0]), 0.168)
    result = redetect_masks_and_targets(image, variance, bits, members,
                                       np.array([64.0, 64.0]), 0.168,
                                       unresolved_bcg_policy="mask_all")

    assert result["diagnostics"]["n_bcg_candidates"] == 0
    assert result["diagnostics"]["bcg_mask_policy"] == "mask_all_sources"
    assert len(result["targets_pixel_radius"]) == 1
