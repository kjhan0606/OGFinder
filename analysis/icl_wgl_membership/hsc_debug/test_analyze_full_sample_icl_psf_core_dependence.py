from __future__ import annotations

import pytest

from analyze_full_sample_icl_psf_core_dependence import join_rows, partial_correlation


def test_join_requires_exact_target_control_psf_pair():
    profiles = [{"id": "T1", "filter": "i", "ra_deg": "10", "dec_deg": "1"}]
    icl = [{"id": "T1", "filter": "i", "control_id": "C1",
            "delta_sb_200_300": "2.0", "coverage_200_300": "0.5",
            "z": "0.3", "richness": "20", "center_masked_k2p5": "False"}]
    psf = [{"target_cluster_id": "T1", "control_id": "C1",
            "target_minus_control_moment_fwhm_major_arcsec": "0.02",
            "target_moment_fwhm_major_arcsec": "1.0",
            "control_moment_fwhm_major_arcsec": "0.98"}]
    targets = {"T1": {"id": "T1", "ra_deg": 10.0, "dec_deg": 1.0}}
    controls = {"C1": {"id": "C1", "target_cluster_id": "T1",
                        "ra_deg": 10.2, "dec_deg": 1.1}}
    row, = join_rows(profiles, icl, psf, targets, controls)
    assert row["psf_core_delta_fwhm_major_arcsec"] == pytest.approx(0.02)
    assert row["icl_delta_njy_arcsec2"] == pytest.approx(2.0)
    assert row["center_masked"] == 0


def test_join_rejects_wrong_matched_control_link():
    profiles = [{"id": "T1", "filter": "i"}]
    icl = [{"id": "T1", "filter": "i", "control_id": "C1",
            "delta_sb_200_300": "2", "coverage_200_300": ".5",
            "z": ".3", "richness": "20", "center_masked_k2p5": "True"}]
    psf = [{"target_cluster_id": "T1", "control_id": "C1",
            "target_minus_control_moment_fwhm_major_arcsec": ".02",
            "target_moment_fwhm_major_arcsec": "1",
            "control_moment_fwhm_major_arcsec": ".98"}]
    with pytest.raises(ValueError, match="control-to-target manifest link"):
        join_rows(profiles, icl, psf, {"T1": {"id": "T1"}},
                  {"C1": {"id": "C1", "target_cluster_id": "T2"}})


def test_partial_correlation_recovers_linear_relation_after_adjustment():
    rows = []
    for i in range(20):
        z = i / 10
        x = i * 0.3 + (i % 3)
        rows.append({"icl_delta_njy_arcsec2": 3 * x + 2 * z,
                     "psf_core_delta_fwhm_major_arcsec": x,
                     "redshift": z})
    result = partial_correlation(rows, ("redshift",))
    assert result["partial_r"] == pytest.approx(1.0)
