import numpy as np
import pytest

from audit_satellite_icl_coverage_adjustment import (
    block_bootstrap, join_coverage, partial_correlation,
)


def sample_rows(n=12):
    rows = []
    coverage = []
    for i in range(n):
        z = i / n
        richness = float((i * 7 + i * i * 3) % n)
        target = 0.45 + 0.5 * ((i * 5) % n) / n
        control = 0.47 + 0.5 * ((i * 7) % n) / n
        common = min(target, control)
        icl = 2.0 * z + 0.1 * richness + np.sin(i)
        satellite = 1.5 * z - 0.2 * richness + 0.7 * icl + np.cos(i)
        rows.append({"cluster_id": str(i), "filter": "i2",
                     "redshift": str(z), "richness": str(richness),
                     "icl_delta_200_300": str(icl),
                     "satellite_delta_50_75": str(satellite),
                     "target_ra_deg": str(30 + i), "target_dec_deg": "1",
                     "control_ra_deg": str(70 + i), "control_dec_deg": "1"})
        coverage.append({"id": str(i), "filter": "i2", "radial_bin_kpc": "200-300",
                         "target_coverage": str(target), "control_coverage": str(control),
                         "common_coverage": str(common)})
    return rows, coverage


def test_join_coverage_validates_minimum_and_keeps_fixed_rows():
    rows, coverage = sample_rows()
    joined = join_coverage(rows, coverage)
    assert len(joined) == len(rows)
    assert joined[0]["icl_common_coverage"] == float(coverage[0]["common_coverage"])
    bad = [dict(row) for row in coverage]
    bad[0]["common_coverage"] = "0.99"
    with pytest.raises(ValueError, match="inconsistent"):
        join_coverage(rows, bad)


def test_partial_correlation_rejects_collinear_covariates():
    rows, _ = sample_rows()
    with pytest.raises(ValueError, match="rank deficient"):
        partial_correlation(rows, ("redshift", "redshift"))


def test_spatial_block_bootstrap_is_paired_and_reproducible():
    rows, coverage = sample_rows()
    rows = join_coverage(rows, coverage)
    labels = np.repeat(np.arange(6), 2)
    first = block_bootstrap(rows, labels, replicas=100, seed=11)
    second = block_bootstrap(rows, labels, replicas=100, seed=11)
    assert first == second
    assert first["replicas_valid"] + first["replicas_invalid"] == 100
    assert set(first["partial_r_percentiles_2_5_50_97_5"]) == {
        "redshift_richness", "redshift_richness_common_coverage",
        "redshift_richness_target_control_coverage"}
