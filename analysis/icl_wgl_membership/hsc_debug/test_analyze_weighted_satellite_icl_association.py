import numpy as np
import pytest

from analyze_weighted_satellite_icl_association import (
    block_bootstrap_correlations,
    cluster_satellite_summaries,
    partial_pearson,
)


def member(cluster, delta, weight, band="i2", radial="50-75", index=0):
    return {"cluster_id": cluster, "filter": band, "radial_bin_kpc": radial,
            "member_index": str(index), "delta": str(delta), "camira_weight": str(weight)}


def test_cluster_satellite_summary_uses_member_weight_and_effective_count():
    result = cluster_satellite_summaries([
        member("A", 1.0, 0.25, index=0), member("A", 3.0, 0.75, index=1),
        member("B", 8.0, 1.0), member("A", 99.0, 1.0, radial="30-40", index=2),
    ])
    assert result[("A", "i2")]["satellite_unweighted"] == pytest.approx(2.0)
    assert result[("A", "i2")]["satellite_weighted"] == pytest.approx(2.5)
    assert result[("A", "i2")]["n_satellites"] == 2
    assert result[("A", "i2")]["n_eff_weight"] == pytest.approx(1.6)


def test_cluster_satellite_summary_rejects_duplicate_satellite_profile():
    with pytest.raises(ValueError, match="duplicate satellite"):
        cluster_satellite_summaries([member("A", 1, 0.5), member("A", 2, 0.6)])


def test_partial_pearson_removes_linear_redshift_and_richness_trends():
    z = np.linspace(0.3, 0.55, 20)
    richness = np.linspace(15, 60, 20) ** 1.1
    x = np.sin(np.arange(20) * 0.8) + 0.7 * z + 0.002 * richness
    y = 0.6 * x + 0.4 * np.cos(np.arange(20) * 0.37) - 0.2 * z + 0.001 * richness
    r = partial_pearson(x, y, z, richness)
    assert -1.0 <= r <= 1.0
    with pytest.raises(ValueError, match="rank deficient"):
        partial_pearson(x, y, z, z,)


def test_paired_block_bootstrap_is_reproducible_and_reports_delta_interval():
    rng = np.random.default_rng(31)
    rows = []
    labels = []
    for block in range(8):
        for j in range(4):
            z = 0.30 + 0.02 * block + 0.001 * j
            richness = 15.0 + 5.0 * block + j
            icl = rng.normal() + 0.7 * z + 0.01 * richness
            satellite = 0.4 * icl + rng.normal(scale=0.3) + 0.2 * z - 0.005 * richness
            rows.append({"icl_delta_200_300": icl,
                         "satellite_unweighted": satellite,
                         "satellite_weighted": satellite + rng.normal(scale=0.1),
                         "z": z, "richness": richness})
            labels.append(str(block))
    first = block_bootstrap_correlations(rows, labels, replicas=1000, seed=7)
    second = block_bootstrap_correlations(rows, labels, replicas=1000, seed=7)
    assert first == second
    assert first["n_blocks"] == 8
    assert first["science_measurement"] is False
    assert len(first["bootstrap_ci95_delta_r"]) == 2
