import numpy as np

from analyze_satellite_common_support import (
    association_result,
    block_bootstrap_mean,
    sky_cells,
)


def test_sky_cells_use_five_degree_ra_dec_grid():
    labels = sky_cells([0.1, 4.9, 5.1], [-0.1, -0.1, 0.1], width=5)
    assert labels.tolist() == ["0_17", "0_17", "1_18"]


def test_block_bootstrap_is_deterministic_and_preserves_bin_covariance():
    values = np.asarray([[1, 3], [2, 5], [4, 9], [8, 17]], dtype=float)
    labels = np.asarray(["a", "a", "b", "c"])
    first = block_bootstrap_mean(values, labels, replicas=500, seed=51)
    second = block_bootstrap_mean(values, labels, replicas=500, seed=51)
    assert all(np.allclose(left, right) if isinstance(left, np.ndarray) else left == right
               for left, right in zip(first, second))
    means, low, high, covariance, n_blocks, largest = first
    assert means.shape == low.shape == high.shape == (2,)
    assert covariance.shape == (2, 2)
    assert covariance[0, 1] > 0
    assert n_blocks == 3
    assert largest == 2


def test_common_support_association_uses_sky_block_bootstrap():
    rows = []
    for index in range(12):
        rows.append({
            "icl_delta_200_300": float(index + (index % 3)),
            "satellite_delta_50_75": float(index * 0.4 + (index % 2)),
            "redshift": 0.2 + 0.03 * index,
            "richness": 15 + index,
            "sky_cell": f"{index % 6}_17",
        })
    first = association_result(rows, "all", "i", replicas=300, seed=91)
    second = association_result(rows, "all", "i", replicas=300, seed=91)
    assert first == second
    assert first["status"] == "exploratory"
    assert first["n_clusters"] == 12
    assert first["n_spatial_blocks"] == 6
    assert first["valid_spearman_replicas"] == 300
    assert len(first["spearman_block_bootstrap_ci95"]) == 2
