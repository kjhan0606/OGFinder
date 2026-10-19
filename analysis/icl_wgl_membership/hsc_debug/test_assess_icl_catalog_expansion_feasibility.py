from assess_icl_catalog_expansion_feasibility import assess


def power_summary():
    return {
        "schema": "astrafek-icl-spatial-block-power-plan-v1",
        "status": "planning_approximation",
        "science_measurement": False,
        "observed_blocks": 33,
        "n_comparisons": 30,
        "rows": [{"effect_r_scenario": 0.293,
                  "required_blocks": {"power_0.8": {"adjusted_alpha": 180}}}],
    }


def expansion_summary():
    return {
        "schema": "camira-offline-expansion-pool-audit-v1",
        "status": "complete_catalog_pool_diagnostic",
        "science_measurement": False,
        "cas_query_submitted": False,
        "coordinates_exported": False,
        "images_downloaded": False,
        "current_occupied_5deg_cells": 78,
        "selection_scenarios": {
            "current_selection": {"greedy_additions_after_6arcmin_separation": 57,
                                  "additional_5deg_cells_occupied": 0},
            "wide_z_cluster_redshift": {"greedy_additions_after_6arcmin_separation": 2596,
                                         "additional_5deg_cells_occupied": 6,
                                         "total_5deg_cells_occupied": 84},
        },
    }


def test_target_cells_upper_bound_connected_pair_blocks():
    result = assess(power_summary(), expansion_summary())
    assert result["required_independent_blocks"] == 180
    assert result["observed_independent_blocks"] == 33
    assert result["broadest_pool_new_target_5deg_cells"] == 6
    assert result["broadest_pool_total_5deg_cells"] == 84
    assert result["broadest_pool_target_cell_upper_bound_on_connected_blocks"] == 84
    assert result["minimum_block_gap_even_at_target_cell_upper_bound"] == 96
    assert result["science_measurement"] is False
    assert result["images_downloaded"] is False


def test_rejects_candidate_input_that_already_submitted_a_query():
    expansion = expansion_summary()
    expansion["cas_query_submitted"] = True
    try:
        assess(power_summary(), expansion)
    except ValueError as error:
        assert "offline" in str(error)
    else:
        raise AssertionError("query-submitted candidate pool was incorrectly accepted")
