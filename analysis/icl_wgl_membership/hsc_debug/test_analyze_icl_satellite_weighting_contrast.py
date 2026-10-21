import numpy as np
import pytest

from analyze_icl_satellite_weighting_contrast import paired_bootstrap


def _paired_rows(seed=23):
    rng = np.random.default_rng(seed)
    equal, weighted, blocks = [], [], []
    for block in range(8):
        for index in range(4):
            icl = float(rng.normal())
            row = {
                "icl_delta_200_300": icl,
                "satellite_delta_50_75": 0.35 * icl + float(rng.normal()),
                "z": float(rng.normal()),
                "richness": float(rng.normal()),
                "cluster_id": f"{block}-{index}",
                "filter": "i2",
                "icl_radial_bin_kpc": "200-300",
                "satellite_radial_bin_kpc": "50-75",
            }
            equal.append(row)
            weighted.append({
                **row,
                "satellite_delta_50_75": row["satellite_delta_50_75"] + 0.1 * icl,
            })
            blocks.append(str(block))
    return equal, weighted, blocks


def test_paired_block_bootstrap_reproducible_and_reports_delta_interval():
    equal, weighted, blocks = _paired_rows()
    first = paired_bootstrap(equal, weighted, blocks, replicas=1000, seed=91)
    second = paired_bootstrap(equal, weighted, blocks, replicas=1000, seed=91)
    assert first == second
    assert first["n_clusters"] == 32
    assert first["n_spatial_blocks"] == 8
    assert first["valid_replicates"] >= 800
    assert first["paired_delta_ci95_block"][0] <= first["weighted_minus_equal_partial_r"]
    assert first["weighted_minus_equal_partial_r"] <= first["paired_delta_ci95_block"][1]


def test_paired_block_bootstrap_rejects_mismatched_pair_keys():
    equal, weighted, blocks = _paired_rows()
    weighted[0]["cluster_id"] = "wrong-cluster"
    with pytest.raises(ValueError, match="pair key differs"):
        paired_bootstrap(equal, weighted, blocks, replicas=1000, seed=91)


def test_paired_block_bootstrap_requires_enough_replicates():
    equal, weighted, blocks = _paired_rows()
    with pytest.raises(ValueError, match="at least 1000"):
        paired_bootstrap(equal, weighted, blocks, replicas=999, seed=91)
