import importlib.util
from pathlib import Path

import numpy as np
import pytest


SCRIPT = Path(__file__).with_name("stack_hsc_icl_target_control_2d.py")
SPEC = importlib.util.spec_from_file_location("stack_hsc_icl_target_control_2d", SCRIPT)
stack = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(stack)


def test_pair_common_support_is_strict_and_difference_is_target_minus_control():
    target = np.array([[4.0, 5.0], [np.nan, 8.0]])
    control = np.array([[1.0, 2.0], [3.0, np.nan]])
    target_weight = np.array([[1, 1], [1, 0]])
    control_weight = np.array([[1, 0], [1, 1]])
    t, c, d, valid = stack.pair_common_support(target, control, target_weight, control_weight)
    np.testing.assert_array_equal(valid, [[True, False], [False, False]])
    assert t[0, 0] == 4.0 and c[0, 0] == 1.0 and d[0, 0] == 3.0
    assert np.isnan(d[0, 1])


def test_online_map_uses_each_pair_once_and_reports_between_pair_sem():
    acc = stack.OnlineMap((2, 2))
    acc.add(np.array([[1.0, 2.0], [3.0, 4.0]]), np.array([[1, 1], [0, 1]], bool))
    acc.add(np.array([[3.0, 8.0], [9.0, 10.0]]), np.array([[1, 0], [1, 1]], bool))
    mean, sem = acc.result()
    np.testing.assert_allclose(mean[0, 0], 2.0)
    np.testing.assert_allclose(sem[0, 0], 1.0)
    np.testing.assert_allclose(acc.count, [[2, 1], [1, 2]])
    assert np.isnan(sem[0, 1])


def test_pair_radial_bins_record_common_area_coverage():
    axis = np.array([-1.0, 0.0, 1.0])
    target = np.full((3, 3), 4.0)
    control = np.full((3, 3), 1.0)
    valid = np.ones((3, 3), bool)
    valid[0, 0] = False
    delta = target - control
    rows = stack.pair_radial_profile(target, control, delta, valid, axis,
                                     edges_kpc=np.array([0.0, 1.0, 3.0]))
    assert rows[0]["n_valid_pixels"] == 1
    assert rows[0]["delta_mean_njy_arcsec2"] == 3.0
    assert rows[1]["n_valid_pixels"] == 7
    assert np.isclose(rows[1]["coverage_fraction"], 7 / 8)


def test_radial_stack_weights_pairs_equally_and_records_radial_n():
    rows = [
        {"filter": "i", "profile": [
            {"target_mean_njy_arcsec2": 1.0, "control_mean_njy_arcsec2": 0.0,
             "delta_mean_njy_arcsec2": 1.0, "coverage_fraction": 1.0} for _ in range(20)]},
        {"filter": "i", "profile": [
            {"target_mean_njy_arcsec2": 3.0, "control_mean_njy_arcsec2": 0.0,
             "delta_mean_njy_arcsec2": 3.0, "coverage_fraction": 0.5} for _ in range(20)]},
    ]
    out = stack.summarize_profiles(rows, filters=("i",))
    assert len(out) == 20
    assert out[0]["n_pairs_delta"] == 2
    assert out[0]["mean_delta_njy_arcsec2"] == 2.0
    assert out[0]["pair_sem_delta_njy_arcsec2"] == 1.0
    assert out[0]["median_pair_common_coverage"] == 0.75


def test_native_flux_factor_is_njy_per_pixel():
    info = {"input_unit": "instrumental_flux/pixel", "fluxmag0": 3.631e12,
            "pixel_scale_arcsec": 0.168}
    assert np.isclose(stack.pixel_flux_factor(info), 1.0)


def test_json_safe_preserves_boolean_type():
    assert stack.json_safe(False) is False
    assert stack.json_safe(np.bool_(False)) is False


def test_mask_factor_selects_matching_target_mask_extension():
    assert [stack.mask_hdu_name(factor) for factor in stack.MASK_FACTORS] == [
        "MASK_K2P5", "MASK_K3P5", "MASK_K5P0"]
    with pytest.raises(ValueError):
        stack.mask_hdu_name(4.0)


def test_candidate_control_audit_resolves_only_exact_verified_pairs():
    audit = {
        "schema": "astrafek-s23b-candidate-existing-cutout-support-v1",
        "status": "header_wcs_audit_complete", "science_measurement": False,
        "hsc_query_submitted": False, "data_downloaded": False,
        "image_pixels_read": False, "n_candidates": 1333,
        "candidate_table_sha256": "1b59ca9941849fa7970c43c6e0374db6587716fc4ae5d71b11f84a245ebe6b2c",
        "target": {"n_images": 800, "n_images_with_candidate_centers": 0,
                   "n_unique_candidates_inside": 0, "matches": []},
        "matched_controls": {
            "inputs_sha256": {"manifest": "m" * 64, "download_status": "s" * 64},
            "n_images": 798, "n_images_with_candidate_centers": 1,
            "n_unique_candidates_inside": 1,
            "matches": [{"target_cluster_id": "t1", "frame_id": "c1", "filter": "i",
                         "n_candidate_centers_inside": 1,
                         "candidates": [{"candidate_id": "42", "z": 0.4, "richness": 17.0}]}],
        },
    }
    controls = {"t1": {"id": "c1", "coadd_header_filter": "i"}}
    result = stack.candidate_control_exclusions(audit, controls, "m" * 64, "s" * 64)
    assert result == [{
        "target_id": "t1", "control_id": "c1",
        "reason": "s23b_candidate_center_inside_matched_control",
        "filter": "i", "candidate_ids": ["42"],
        "candidate_redshifts": [0.4], "candidate_richness": [17.0],
    }]

    audit["matched_controls"]["matches"][0]["frame_id"] = "wrong"
    with pytest.raises(ValueError, match="unknown or duplicate"):
        stack.candidate_control_exclusions(audit, controls, "m" * 64, "s" * 64)


def test_candidate_control_pairs_already_rejected_for_filter_mismatch_are_not_double_counted():
    matches = [
        {"target_id": "eligible", "control_id": "c1", "filter": "i",
         "candidate_ids": ["101"]},
        {"target_id": "mismatch", "control_id": "c2", "filter": "i",
         "candidate_ids": ["102"]},
    ]
    newly_excluded, already_excluded = stack.split_candidate_exclusions(
        matches, {"eligible": "i", "mismatch": "i2"})
    assert [row["target_id"] for row in newly_excluded] == ["eligible"]
    assert already_excluded[0]["target_id"] == "mismatch"
    assert already_excluded[0]["reason"] == "already_excluded_target_control_filter_mismatch"
    with pytest.raises(ValueError, match="unknown target"):
        stack.split_candidate_exclusions(matches, {"eligible": "i"})
