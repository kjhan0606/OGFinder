import pytest

from audit_empirical_psf_outer_star_support import summarize_thresholds


def test_threshold_sweep_counts_stars_and_independent_frames_separately():
    candidates = [
        {"candidate_id": "a", "cluster_id": "f1", "cluster_centric_radius_kpc": 650.0, "eligible": True},
        {"candidate_id": "b", "cluster_id": "f1", "cluster_centric_radius_kpc": 900.0, "eligible": True},
        {"candidate_id": "c", "cluster_id": "f2", "cluster_centric_radius_kpc": 750.0, "eligible": True},
        {"candidate_id": "d", "cluster_id": "f3", "cluster_centric_radius_kpc": 850.0, "eligible": False},
    ]
    rows = summarize_thresholds(candidates, (600, 800))
    assert rows[0]["n_star_control_pairs"] == 3
    assert rows[0]["n_independent_target_frames"] == 2
    assert rows[0]["n_frames_with_multiple_stars"] == 1
    assert rows[1]["n_star_control_pairs"] == 1
    assert rows[1]["n_independent_target_frames"] == 1


def test_threshold_sweep_rejects_unsafe_or_duplicate_thresholds():
    with pytest.raises(ValueError, match="at least 300"):
        summarize_thresholds([], (250,))
    with pytest.raises(ValueError, match="unique"):
        summarize_thresholds([], (600, 600))
