import numpy as np
import pytest

from stack_empirical_psf_outer_field_wings_frame_bootstrap import stack_frame_pairs


def _rows():
    rows = []
    for candidate, frame, value in (("a", "f1", 1.0), ("b", "f1", 3.0), ("c", "f2", 5.0)):
        for role, profile in (("target", value + 2.0), ("control", 2.0)):
            rows.append({
                "status": "ok", "profile_per_target_core_per_pixel": profile,
                "cluster_id": frame, "candidate_id": candidate, "role": role,
                "filter": "i2",
                "sky_r_min_pix": 120, "sky_r_max_pix": 150,
                "r_min_kpc": 200.0, "r_max_kpc": 225.0,
            })
    return rows


def test_stars_within_a_frame_are_averaged_before_frame_bootstrap():
    result = stack_frame_pairs(_rows(), n_bootstrap=1000, seed=2)[0]
    assert result["n_complete_star_pairs"] == 3
    assert result["n_independent_target_frames"] == 2
    assert result["n_frames_with_multiple_complete_stars"] == 1
    # f1 contributes mean(1, 3)=2; f2 contributes 5; frames are equally weighted.
    assert result["frame_equal_paired_difference_mean"] == pytest.approx(3.5)
    assert result["bootstrap_unit"].startswith("target/control HSC coadd frame")


def test_blocked_profiles_do_not_enter_the_stack():
    rows = _rows()
    rows[0]["status"] = "blocked_radial_coverage"
    result = stack_frame_pairs(rows, n_bootstrap=1000, seed=2)[0]
    assert result["n_complete_star_pairs"] == 2
    assert result["n_independent_target_frames"] == 2


def test_i_and_i2_profiles_are_kept_in_separate_stack_cells():
    i2 = _rows()
    i = []
    for row in i2:
        copied = dict(row, filter="i")
        if copied["role"] == "target":
            copied["profile_per_target_core_per_pixel"] += 10.0
        i.append(copied)
    result = stack_frame_pairs(i2 + i, n_bootstrap=1000, seed=2)
    by_filter = {row["filter"]: row for row in result}
    assert set(by_filter) == {"i", "i2"}
    assert by_filter["i"]["frame_equal_paired_difference_mean"] == pytest.approx(13.5)
    assert by_filter["i2"]["frame_equal_paired_difference_mean"] == pytest.approx(3.5)
