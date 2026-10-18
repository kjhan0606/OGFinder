from analyze_icl_satellite_cut_sensitivity import (
    connected_pair_blocks,
    summarize_association,
)


def test_connected_blocks_link_pairs_sharing_control_sky_cell():
    rows = [
        {"target_ra_deg": 1, "target_dec_deg": 1,
         "control_ra_deg": 100, "control_dec_deg": 1},
        {"target_ra_deg": 11, "target_dec_deg": 1,
         "control_ra_deg": 101, "control_dec_deg": 1},
    ]
    labels, target_cells, control_cells, blocks, largest = connected_pair_blocks(rows)
    assert labels[0] == labels[1]
    assert target_cells == 2
    assert control_cells == 1
    assert blocks == 1
    assert largest == 2


def test_spatial_block_association_is_reproducible():
    rows = []
    for index in range(15):
        rows.append({
            "icl_delta_200_300": float(index % 5 + index * 0.17),
            "satellite_delta_50_75": float((index * 7) % 11 + index * 0.09),
            "z": 0.2 + 0.013 * (index % 7),
            "richness": 10 + (index * 3) % 9,
            "target_ra_deg": 2 + index * 7,
            "target_dec_deg": -30 + index * 3,
            "control_ra_deg": 140 + index * 6,
            "control_dec_deg": 20 - index * 2,
        })
    first = summarize_association(rows, 2.5, 0.4, "i", "all", 200, 42)
    second = summarize_association(rows, 2.5, 0.4, "i", "all", 200, 42)
    assert first == second
    assert first["status"] == "exploratory"
    assert first["n_clusters"] == 15
    assert first["n_joint_blocks"] >= 2
    assert first["valid_bootstrap_replicas"] > 0
