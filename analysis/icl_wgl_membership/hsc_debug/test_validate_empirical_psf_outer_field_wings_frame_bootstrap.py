import csv

import pytest

from validate_empirical_psf_outer_field_wings_frame_bootstrap import partial_support, recompute_stack


def _profiles():
    rows = []
    for candidate, frame, value, area, sectors in (
            ("a", "f1", 1.0, 0.30, 6), ("b", "f1", 3.0, 0.25, 6), ("c", "f2", 5.0, 0.55, 8)):
        for role, profile in (("target", value + 2), ("control", 2.0)):
            rows.append({"candidate_id": candidate, "cluster_id": frame, "role": role,
                "filter": "i2",
                "status": "ok",
                "sky_r_min_pix": "120", "sky_r_max_pix": "150", "r_min_kpc": "200", "r_max_kpc": "225",
                "profile_per_target_core_per_pixel": str(profile), "area_fraction": str(area),
                "n_covered_sectors": str(sectors), "background": "0.0", "target_core_flux": "10.0",
                "raw_mean": "1.0"})
    return rows


def test_independent_recomputation_averages_stars_before_frame_bootstrap():
    stack = recompute_stack(_profiles(), n_bootstrap=1000, seed=5)
    assert len(stack) == 1
    assert stack[0]["n_complete_star_pairs"] == 3
    assert stack[0]["n_independent_target_frames"] == 2
    assert stack[0]["frame_equal_paired_difference_mean"] == pytest.approx(3.5)


def test_partial_support_gate_counts_frames_not_pixels_as_independent():
    support = partial_support(_profiles(), area_cut=0.2, min_sectors=6)
    assert len(support) == 1
    assert support[0]["n_complete_star_pairs"] == 3
    assert support[0]["n_independent_frames"] == 2
    assert partial_support(_profiles(), area_cut=0.5, min_sectors=8)[0]["n_complete_star_pairs"] == 1


def test_recomputation_keeps_subfilters_separate():
    profiles = _profiles()
    i_profiles = []
    for row in profiles:
        copied = dict(row, filter="i")
        if copied["role"] == "target":
            copied["profile_per_target_core_per_pixel"] = str(float(copied["profile_per_target_core_per_pixel"]) + 10)
        i_profiles.append(copied)
    result = recompute_stack(profiles + i_profiles, n_bootstrap=1000, seed=5)
    assert {row["filter"] for row in result} == {"i", "i2"}
    means = {row["filter"]: row["frame_equal_paired_difference_mean"] for row in result}
    assert means["i"] == pytest.approx(means["i2"] + 10)
