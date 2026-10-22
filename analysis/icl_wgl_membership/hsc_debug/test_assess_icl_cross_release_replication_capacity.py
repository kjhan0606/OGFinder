from assess_icl_cross_release_replication_capacity import summarize_cells


def test_s20a_pool_only_adds_two_cells_beyond_s23b_and_union_stays_short():
    discovery = {f"d{i}" for i in range(33)}
    s23b = discovery | {f"b{i}" for i in range(67)}
    s20a = discovery | {f"b{i}" for i in range(49)} | {"c0", "c1"}
    result = summarize_cells(discovery, s20a, s23b)
    assert result["s20a_unique_cells_added_beyond_s23b"] == 2
    assert result["s20a_incremental_new_cells_vs_discovery_beyond_s23b"] == 2
    rows = {row["pool"]: row for row in result["pools"]}
    assert rows["S20A_broad"]["n_new_target_cells_vs_i2_discovery"] == 51
    assert rows["S23B_selection_compatible"]["n_new_target_cells_vs_i2_discovery"] == 67
    assert rows["union_catalog_upper_bound"]["n_new_target_cells_vs_i2_discovery"] == 69
    assert rows["union_catalog_upper_bound"]["shortfall_r293_single_endpoint"] == 23
    assert result["union_new_cell_count_is_only_an_upper_bound"] is True


def test_familywise_target_still_exceeds_all_candidate_cells():
    discovery = {f"d{i}" for i in range(33)}
    s23b = discovery | {f"b{i}" for i in range(67)}
    s20a = discovery | {f"b{i}" for i in range(49)} | {"c0", "c1"}
    result = summarize_cells(discovery, s20a, s23b)
    row = result["pools"][2]
    assert row["n_blocks_80pct_r293_bonferroni_30"] == 180
    assert row["n_new_target_cells_vs_i2_discovery"] + len(discovery) == 102
    assert 180 - 102 == 78
