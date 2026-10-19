import importlib.util
from pathlib import Path

import numpy as np


SCRIPT = Path(__file__).with_name("bootstrap_icl_2d_spatial_blocks.py")
SPEC = importlib.util.spec_from_file_location("bootstrap_icl_2d_spatial_blocks", SCRIPT)
bootstrap = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(bootstrap)


def test_pairs_sharing_a_sky_cell_across_target_control_roles_are_connected():
    targets = [
        {"ra_deg": 1.0, "dec_deg": 1.0},
        {"ra_deg": 11.0, "dec_deg": 11.0},
        {"ra_deg": 21.0, "dec_deg": 21.0},
    ]
    controls = [
        {"ra_deg": 11.1, "dec_deg": 11.1},
        {"ra_deg": 31.0, "dec_deg": 31.0},
        {"ra_deg": 41.0, "dec_deg": 41.0},
    ]
    labels = bootstrap.connected_pair_blocks(targets, controls, width=5.0)
    assert labels[0] == labels[1]
    assert labels[2] != labels[0]


def test_block_bootstrap_is_reproducible_and_returns_profile_covariance():
    values = np.asarray([[i + radius for radius in range(5)] for i in range(12)], dtype=float)
    labels = np.asarray([f"b{i // 2}" for i in range(12)])
    first = bootstrap.summarize_resamples(values, labels, replicas=500, seed=123)
    second = bootstrap.summarize_resamples(values, labels, replicas=500, seed=123)
    assert first == second
    assert first["n_connected_blocks"] == 6
    assert first["n_pairs_by_annulus"] == [12] * 5
    assert first["point_mean_njy_arcsec2"] == [5.5, 6.5, 7.5, 8.5, 9.5]
    assert len(first["spatial_block_bootstrap"]["covariance_njy2_arcsec4"]) == 5
    assert first["spatial_block_bootstrap"]["complete_profile_replicas"] == 500
