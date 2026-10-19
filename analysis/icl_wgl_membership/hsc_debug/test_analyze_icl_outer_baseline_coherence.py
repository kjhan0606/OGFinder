import numpy as np
import pytest

from analyze_icl_outer_baseline_coherence import bootstrap_coherence, regression_metrics


def test_regression_metrics_recovers_outer_slope_after_covariate_adjustment():
    z = np.linspace(-1.0, 1.0, 40)
    richness = np.sin(np.arange(40) * 0.37)
    outer = 3.0 * z + np.cos(np.arange(40) * 0.21)
    inner = 1.25 * outer + 0.8 * z - 0.4 * richness

    metrics = regression_metrics(inner, outer, np.column_stack((z, richness)))

    assert metrics["partial_slope"] == pytest.approx(1.25, abs=1e-10)
    assert metrics["partial_pearson_r"] == pytest.approx(1.0, abs=1e-10)


def test_regression_metrics_rejects_rank_deficiency_and_nonfinite_values():
    x = np.arange(12.0)
    with pytest.raises(ValueError, match="rank-deficient"):
        regression_metrics(x, x + 1.0, np.column_stack((x, 2 * x)))
    with pytest.raises(ValueError, match="finite"):
        regression_metrics(x, x + 1.0, np.column_stack((x, [0.0] * 11 + [np.nan])))


def test_connected_block_bootstrap_preserves_pair_groups_and_returns_intervals():
    rows = []
    for block in range(15):
        for within in range(4):
            z = 0.1 * within + 0.003 * block**2 + 0.01 * np.sin(block + within)
            richness = 20 + block * 1.3 + within * 0.35 + 0.11 * block * within
            outer = np.sin(block * 0.7 + within * 0.2) + 0.1 * richness
            inner = 1.1 * outer + 0.4 * z + 0.03 * within
            rows.append({
                "ra_deg": block * 10.0 + 0.1 * within,
                "dec_deg": 0.1,
                "control_ra_deg": block * 10.0 + 1.0 + 0.1 * within,
                "control_dec_deg": 0.2,
                "z": z,
                "richness": richness,
                "inner": inner,
                "outer": outer,
                "inner_coverage": 0.3 + 0.01 * within,
                "outer_coverage": 0.2 + 0.02 * block,
            })

    result = bootstrap_coherence(rows, replicas=300, seed=7)

    assert result["n_pairs"] == 60
    assert result["n_connected_blocks"] == 15
    assert result["block_bootstrap_95_intervals"]["partial_slope"]["valid_replicates"] >= 240
    assert result["block_bootstrap_95_intervals"]["partial_pearson_r"]["q025"] <= result["block_bootstrap_95_intervals"]["partial_pearson_r"]["q975"]
