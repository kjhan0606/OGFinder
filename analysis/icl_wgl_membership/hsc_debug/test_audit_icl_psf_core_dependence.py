import numpy as np
import pytest

from audit_icl_psf_core_dependence import block_bootstrap, partial_correlation


def sample_rows(n=24):
    rows = []
    for i in range(n):
        z = i / n
        richness = float((i * 7 + 3 * i * i) % n)
        coverage = 0.45 + 0.5 * ((i * 5 + i * i) % n) / n
        psf = 1.0 + 0.1 * z + 0.02 * np.sin(i)
        icl = 2.0 * z + 0.2 * richness + 3.0 * psf + np.cos(i)
        rows.append({"icl_delta": float(icl), "psf_fwhm_major_arcsec": float(psf),
                     "redshift": z, "richness": richness,
                     "common_coverage": coverage,
                     "target_coverage": coverage,
                     "control_coverage": min(0.99, coverage + 0.01)})
    return rows


def test_unadjusted_partial_correlation_is_pearson_correlation():
    rows = sample_rows()
    result = partial_correlation(rows, ())
    x = [row["icl_delta"] for row in rows]
    y = [row["psf_fwhm_major_arcsec"] for row in rows]
    assert result["partial_r"] == pytest.approx(np.corrcoef(x, y)[0, 1])


def test_partial_correlation_rejects_rank_deficient_adjustment():
    rows = sample_rows()
    for row in rows:
        row["duplicate_z"] = row["redshift"]
    with pytest.raises(ValueError, match="rank deficient"):
        partial_correlation(rows, ("redshift", "duplicate_z"))


def test_block_bootstrap_is_reproducible_and_reports_invalid_draws():
    rows = sample_rows()
    labels = np.repeat(np.arange(12), 2)
    first = block_bootstrap(rows, labels, ("redshift", "richness", "common_coverage"),
                            replicas=100, seed=17)
    second = block_bootstrap(rows, labels, ("redshift", "richness", "common_coverage"),
                             replicas=100, seed=17)
    assert first == second
    if first["status"] == "complete_descriptive_bootstrap":
        assert first["replicas_valid"] + first["replicas_invalid"] == 100
        assert len(first["percentiles_2_5_50_97_5"]) == 3
