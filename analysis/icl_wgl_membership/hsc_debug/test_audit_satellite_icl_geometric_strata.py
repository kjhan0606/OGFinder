import numpy as np
import pytest

from audit_satellite_icl_geometric_strata import (
    block_bootstrap,
    geometric_class,
    paired_difference_bootstrap,
)


def test_geometric_annulus_overlap_boundaries():
    assert geometric_class(124.999) == "no_geometric_overlap"
    assert geometric_class(125.0) == "overlap_possible"
    assert geometric_class(200.0) == "overlap_possible"
    assert geometric_class(375.0) == "overlap_possible"
    assert geometric_class(375.001) == "no_geometric_overlap"
    with pytest.raises(ValueError, match="finite and non-negative"):
        geometric_class(float("nan"))


def synthetic_pairs():
    rng = np.random.default_rng(12)
    rows = []
    for block in range(8):
        for offset in range(4):
            z = 0.30 + 0.015 * block + 0.001 * offset
            richness = 18.0 + 4.5 * block + offset
            icl = rng.normal() + 0.5 * z + 0.01 * richness
            overlap = 0.55 * icl + rng.normal(scale=0.35)
            no_overlap = 0.1 * icl + rng.normal(scale=0.6)
            rows.append({
                "cluster_id": f"{block}-{offset}", "filter": "i2",
                "icl_delta_200_300": str(icl), "z": str(z),
                "richness": str(richness), "satellite_delta": str(overlap),
                "satellite_overlap_possible": str(overlap),
                "satellite_no_overlap": str(no_overlap),
                "target_ra_deg": str(block * 10.0 + 0.01 * offset),
                "target_dec_deg": "0.0",
                "control_ra_deg": str(block * 10.0 + 0.02 * offset),
                "control_dec_deg": "1.0",
            })
    return rows


def test_spatial_block_bootstrap_is_reproducible():
    rows = synthetic_pairs()
    first = block_bootstrap(rows, "satellite_delta", 1000, np.random.default_rng(5))
    second = block_bootstrap(rows, "satellite_delta", 1000, np.random.default_rng(5))
    assert first == second
    assert first["n_blocks"] == 8
    assert first["replicas_valid"] == 1000
    assert first["status"] == "complete_exploratory"


def test_paired_difference_uses_identical_block_resamples():
    rows = synthetic_pairs()
    first = paired_difference_bootstrap(rows, 1000, np.random.default_rng(17))
    second = paired_difference_bootstrap(rows, 1000, np.random.default_rng(17))
    assert first == second
    assert first["n_blocks"] == 8
    assert first["replicas_valid"] == 1000
    assert first["ci95"][0] <= first["point_delta_r"] <= first["ci95"][1]
