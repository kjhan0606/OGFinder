import pytest

from assess_icl_power_family_sensitivity import assess, load_inputs


def inputs():
    power = {
        "schema": "astrafek-icl-spatial-block-power-plan-v1",
        "status": "planning_approximation",
        "science_measurement": False,
        "observed_blocks": 33,
        "observed_effect_r_scenario": 0.293,
        "n_controls": 2,
    }
    feasibility = {
        "schema": "astrafek-icl-camira-expansion-feasibility-v1",
        "status": "complete_design_diagnostic",
        "science_measurement": False,
        "cas_query_submitted": False,
        "coordinates_exported": False,
        "images_downloaded": False,
        "observed_independent_blocks": 33,
        "broadest_pool_target_cell_upper_bound_on_connected_blocks": 84,
    }
    return power, feasibility


def test_family_and_effect_sensitivity_changes_power_requirement_without_claiming_detection():
    result = assess(*inputs(), effect_sizes=(0.2, 0.293), family_sizes=(1, 30))
    assert result["science_measurement"] is False
    assert result["observed_connected_blocks"] == 33
    assert result["broadest_existing_pool_connected_block_upper_bound"] == 84
    scenarios = {
        (row["effect_r_scenario"], row["test_family_size_scenario"]): row
        for row in result["scenarios"]
    }
    assert scenarios[(0.293, 1)]["required_independent_blocks"] == 92
    assert scenarios[(0.293, 30)]["required_independent_blocks"] == 180
    assert scenarios[(0.2, 1)]["required_independent_blocks"] == 196
    assert scenarios[(0.293, 1)]["shortfall_even_at_upper_bound"] == 8
    assert scenarios[(0.293, 30)]["shortfall_even_at_upper_bound"] == 96


def test_inputs_must_be_offline_and_match_observed_blocks(tmp_path):
    power, feasibility = inputs()
    feasibility["coordinates_exported"] = True
    with pytest.raises(ValueError, match="offline"):
        assess(power, feasibility)

    power, feasibility = inputs()
    feasibility["observed_independent_blocks"] = 34
    power_path = tmp_path / "power.json"
    feasibility_path = tmp_path / "feasibility.json"
    import json
    power_path.write_text(json.dumps(power))
    feasibility_path.write_text(json.dumps(feasibility))
    with pytest.raises(ValueError, match="disagree"):
        load_inputs(power_path, feasibility_path)
