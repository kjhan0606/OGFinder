import pytest

from audit_icl_outer_baseline_regions import leave_one_region_out


def test_leave_one_region_out_tracks_full_and_region_means():
    result = leave_one_region_out([1.0, 3.0, 10.0, 14.0], [0, 0, 1, 1])

    assert result["full_mean"] == pytest.approx(7.0)
    assert result["per_region"] == {
        "0": {"n_pairs": 2, "mean": 2.0},
        "1": {"n_pairs": 2, "mean": 12.0},
    }
    assert result["leave_one_region_out"] == {
        "0": {"n_pairs": 2, "mean": 12.0},
        "1": {"n_pairs": 2, "mean": 2.0},
    }


def test_leave_one_region_out_rejects_single_region_and_nonfinite_values():
    with pytest.raises(ValueError, match="at least two occupied regions"):
        leave_one_region_out([1.0, 2.0], [0, 0])
    with pytest.raises(ValueError, match="finite"):
        leave_one_region_out([1.0, float("nan")], [0, 1])
