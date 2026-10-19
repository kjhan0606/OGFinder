import pytest

from prepare_wl_random_centers import deterministic_select


def test_selection_is_balanced_bounded_and_reproducible():
    rows = [
        {"lens_id": "L1", "ra_deg": 10.0, "dec_deg": 1.0},
        {"lens_id": "L1", "ra_deg": 11.0, "dec_deg": 1.0},
        {"lens_id": "L1", "ra_deg": 12.0, "dec_deg": 1.0},
        {"lens_id": "L2", "ra_deg": 20.0, "dec_deg": 2.0},
    ]
    redshifts = {"L1": 0.3, "L2": 0.4, "L3": 0.5}

    selected, counts = deterministic_select(rows, redshifts, 2, "seed-A")
    repeat, repeat_counts = deterministic_select(rows, redshifts, 2, "seed-A")

    assert selected == repeat
    assert counts == repeat_counts
    assert len(selected) == 3
    assert counts == {
        "L1": {"available": 3, "selected": 2},
        "L2": {"available": 1, "selected": 1},
        "L3": {"available": 0, "selected": 0},
    }
    assert {row["z_lens"] for row in selected if row["matched_lens_id"] == "L1"} == {0.3}
    assert len({row["random_center_id"] for row in selected}) == len(selected)


def test_selection_rejects_unknown_lenses_invalid_coordinates_and_duplicates():
    redshifts = {"L1": 0.3}
    with pytest.raises(ValueError, match="unknown matched lens"):
        deterministic_select([{"lens_id": "L2", "ra_deg": 10, "dec_deg": 1}], redshifts, 1, "x")
    with pytest.raises(ValueError, match="invalid candidate coordinates"):
        deterministic_select([{"lens_id": "L1", "ra_deg": 360, "dec_deg": 1}], redshifts, 1, "x")
    with pytest.raises(ValueError, match="duplicate coordinates"):
        deterministic_select([
            {"lens_id": "L1", "ra_deg": 10, "dec_deg": 1},
            {"lens_id": "L1", "ra_deg": 10, "dec_deg": 1},
        ], redshifts, 2, "x")


def test_selection_requires_explicit_positive_quota_and_seed():
    with pytest.raises(ValueError, match="positive integer"):
        deterministic_select([], {}, 0, "seed")
    with pytest.raises(ValueError, match="non-empty"):
        deterministic_select([], {}, 1, "")
