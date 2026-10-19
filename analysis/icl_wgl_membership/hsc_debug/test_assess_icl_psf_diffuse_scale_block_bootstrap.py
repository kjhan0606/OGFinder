import numpy as np
import pytest

from assess_icl_psf_diffuse_scale_block_bootstrap import connected_block_bootstrap


def test_connected_block_bootstrap_keeps_scales_paired_and_is_reproducible():
    labels = np.repeat([f"block-{index}" for index in range(8)], 2)
    base = np.linspace(0.8, 1.0, len(labels))
    values = np.column_stack([base, base - 0.1, base - 0.2, base - 0.3])

    losses = values - 1.0
    first = connected_block_bootstrap(values, losses, labels, replicas=1000, seed=27)
    second = connected_block_bootstrap(values, losses, labels, replicas=1000, seed=27)

    assert first == second
    assert first["n_connected_blocks"] == 8
    assert first["paired_median_retention_change_from_100kpc"]["200"]["point"] == pytest.approx(-0.1)
    assert first["paired_median_retention_change_from_100kpc"]["500"]["bootstrap_95_interval"] == pytest.approx([-0.3, -0.3])


def test_connected_block_bootstrap_requires_connected_spatial_support():
    values = np.ones((10, 4))
    labels = np.repeat(["block-a", "block-b", "block-c", "block-d", "block-e"], 2)

    with pytest.raises(ValueError, match="at least six connected blocks"):
        connected_block_bootstrap(values, values - 1, labels, replicas=1000, seed=1)


def test_connected_block_bootstrap_rejects_nonfinite_or_misaligned_values():
    labels = np.asarray([f"block-{i}" for i in range(8)])
    values = np.ones((8, 4))
    values[0, 2] = np.nan

    with pytest.raises(ValueError, match="values must all be finite"):
        connected_block_bootstrap(values, values - 1, labels, replicas=1000, seed=1)
    with pytest.raises(ValueError, match="one block label per row"):
        connected_block_bootstrap(values, values - 1, labels[:-1], replicas=1000, seed=1)
