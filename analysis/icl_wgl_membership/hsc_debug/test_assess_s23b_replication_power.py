import pytest

from assess_s23b_replication_power import evaluate


def ledger():
    return {
        "schema": "astrafek-s23b-spatial-candidate-ledger-v1",
        "status": "offline_coordinate_and_footprint_ledger_complete",
        "science_measurement": False,
        "statistics": {
            "candidate_count": 1333,
            "association_unique_target_cells": {"i2": 33},
            "candidate_cells_new_vs_association": {"i2": 67},
            "candidate_cell_count_is_only_a_target_footprint_upper_bound": True,
        },
    }


def test_candidate_pool_is_below_single_endpoint_power_for_r293():
    result = evaluate(ledger())
    row = next(row for row in result["scenarios"] if row["assumed_partial_r"] == 0.293)
    assert result["s23b_new_i2_target_cell_upper_bound"] == 67
    assert result["combined_target_cell_upper_bound"] == 100
    assert row["required_blocks_80pct_one_frozen_endpoint_alpha_0p05"] == 92
    assert row["additional_blocks_needed_beyond_s23b_cap"] == 25
    assert row["nominal_power_at_s23b_new_target_cell_cap"] == pytest.approx(0.662, abs=0.001)
    assert row["additional_blocks_needed_beyond_current_plus_s23b_target_cell_cap"] == 80
    assert result["science_measurement"] is False


def test_smaller_effect_requires_more_independent_blocks():
    result = evaluate(ledger())
    row = next(row for row in result["scenarios"] if row["assumed_partial_r"] == 0.20)
    assert row["required_blocks_80pct_one_frozen_endpoint_alpha_0p05"] == 196
    assert row["additional_blocks_needed_beyond_s23b_cap"] == 129


@pytest.mark.parametrize("changes", [
    {"status": "bad"},
    {"science_measurement": True},
    {"statistics": {"candidate_count": 1333, "association_unique_target_cells": {"i2": 33},
                    "candidate_cells_new_vs_association": {"i2": 68},
                    "candidate_cell_count_is_only_a_target_footprint_upper_bound": True}},
])
def test_rejects_unreviewed_or_changed_input(changes):
    value = ledger()
    value.update(changes)
    with pytest.raises(ValueError):
        evaluate(value)
