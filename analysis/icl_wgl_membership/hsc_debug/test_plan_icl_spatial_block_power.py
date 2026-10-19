import pytest

from plan_icl_spatial_block_power import plan, required_units, two_sided_power


def test_required_units_is_minimal_for_multiplicity_adjusted_power():
    alpha = 0.05 / 30
    n = required_units(0.293, alpha, 0.8, n_controls=2)
    assert two_sided_power(n, 0.293, alpha, n_controls=2) >= 0.8
    assert two_sided_power(n - 1, 0.293, alpha, n_controls=2) < 0.8


def test_power_increases_with_blocks_and_effect_and_decreases_with_correction():
    base = two_sided_power(33, 0.293, 0.05 / 30, n_controls=2)
    assert two_sided_power(50, 0.293, 0.05 / 30, n_controls=2) > base
    assert two_sided_power(33, 0.4, 0.05 / 30, n_controls=2) > base
    assert two_sided_power(33, 0.293, 0.05, n_controls=2) > base


def test_plan_distinguishes_independent_blocks_from_cluster_rows():
    result = plan(33, 0.293, 30)
    assert result["status"] == "planning_approximation"
    assert result["science_measurement"] is False
    assert result["observed_blocks"] == 33
    observed = next(row for row in result["rows"]
                    if row["effect_r_scenario"] == 0.293)
    assert observed["required_blocks"]["power_0.8"]["adjusted_alpha"] > 33
    assert "not clusters" in result["limitations"][2]


@pytest.mark.parametrize("args", [
    (4, 0.3, 0.05, 2),
    (33, 0.0, 0.05, 2),
    (33, 1.0, 0.05, 2),
    (33, 0.3, 0.0, 2),
])
def test_invalid_inputs_are_rejected(args):
    with pytest.raises(ValueError):
        two_sided_power(*args)
