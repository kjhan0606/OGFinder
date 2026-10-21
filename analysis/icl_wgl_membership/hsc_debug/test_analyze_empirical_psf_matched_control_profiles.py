import numpy as np

from analyze_empirical_psf_matched_control_profiles import control_profile, detect_control_masks, frame_summary


def test_constant_control_stamp_has_zero_sky_subtracted_profile():
    image = np.full((401, 401), 7.5)
    variance = np.ones_like(image)
    mask = np.zeros_like(image, dtype=bool)

    rows = control_profile(image, variance, mask, 200.25, 199.75, 100.0)

    assert len(rows) == 22
    assert all(row["status"] == "ok" for row in rows)
    assert all(abs(row["profile_per_target_core_per_pixel"]) < 1e-15 for row in rows)


def test_detected_central_source_blocks_only_inner_profile_bins():
    image = np.zeros((401, 401))
    variance = np.ones_like(image)
    mask = np.zeros_like(image, dtype=bool)
    mask[194:207, 194:207] = True

    rows = control_profile(image, variance, mask, 200, 200, 100.0)

    inner = [row for row in rows if row["r_min_pix"] < 12]
    outer = [row for row in rows if row["r_min_pix"] >= 12]
    assert all(row["status"] == "blocked_control_center_source_or_mask" for row in inner)
    assert all(row["profile_per_target_core_per_pixel"] is None for row in inner)
    assert any(row["status"] == "ok" for row in outer)


def test_frame_summary_averages_stars_before_giving_frames_equal_weight():
    def row(frame, target, control):
        return {
            "pair_status": "complete_profile_pair", "filter": "i2", "kron_factor": 2.5,
            "sky_r_min_pix": 120, "sky_r_max_pix": 150, "r_min_pix": 48.0, "r_max_pix": 64.0,
            "target_cluster_id": frame, "target_profile": target, "control_profile": control,
            "target_minus_control": target - control,
        }

    result = frame_summary([row("A", 1.0, 0.0), row("A", 3.0, 0.0), row("B", 5.0, 0.0)])

    assert len(result) == 1
    assert result[0]["n_frames"] == 2
    assert result[0]["target_frame_equal_mean"] == 3.5
    assert result[0]["paired_difference_frame_equal_mean"] == 3.5


def test_too_few_quality_pixels_is_a_recorded_support_failure():
    image = np.zeros((401, 401))
    variance = np.ones_like(image)
    base_mask = np.ones_like(image, dtype=bool)

    masks, status = detect_control_masks(image, variance, base_mask, 200, 200, 0.168)

    assert status == "blocked_control_too_few_quality_pixels"
    assert len(masks) == 3
    assert all(np.all(mask) for mask in masks)
