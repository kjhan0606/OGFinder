import pytest

from plan_icl_sample_size import build_plan, required_independent_rows


def test_required_rows_increase_for_smaller_effect_and_stricter_alpha():
    modest = required_independent_rows(0.20, 0.80, 0.05, covariates=2)
    smaller = required_independent_rows(0.10, 0.80, 0.05, covariates=2)
    corrected = required_independent_rows(0.20, 0.80, 0.05 / 168, covariates=2)
    assert smaller > modest
    assert corrected > modest


def test_plan_translates_supported_rows_to_parent_sample_and_marks_assumptions():
    result = build_plan(parent_lenses=800, usable_i=26, usable_i2=68,
                       effects=(0.20,), powers=(0.80,), design_effects=(1.0, 1.5,))
    assert result["science_measurement"] is False
    assert result["input"]["observed_support_yield"] == {"i": 26 / 800, "i2": 68 / 800}
    rows = result["scenarios"]
    assert rows[1]["inflated_usable_rows"] >= rows[0]["inflated_usable_rows"]
    assert rows[0]["projected_parent_lenses_i"] > rows[0]["projected_parent_lenses_i2"]
    assert result["input"]["effect_sizes_are_assumptions_not_confirmatory_estimates"] is True


def test_rejects_invalid_design_inputs():
    with pytest.raises(ValueError):
        required_independent_rows(1.0, 0.8, 0.05)
    with pytest.raises(ValueError):
        build_plan(parent_lenses=800, usable_i=801, usable_i2=68)
