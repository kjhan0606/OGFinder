import numpy as np
import pytest

import audit_wl_tomographic_redshift_interaction as audit


def _summary(high_shift=0.0):
    edges = [0.1, 0.2, 0.4]
    profiles = []
    for cohort, baseline in (("lower_z", 100.0), ("upper_z", 60.0)):
        for source_bin in range(1, 5):
            for radial, (r_min, r_max) in enumerate(zip(edges[:-1], edges[1:])):
                offset = (source_bin - 1) * 3.0
                if cohort == "upper_z" and source_bin == 3:
                    offset += high_shift
                profiles.append({
                    "lens_redshift_cohort": cohort,
                    "source_zbin": source_bin,
                    "r_min_mpc": r_min,
                    "r_max_mpc": r_max,
                    "delta_sigma_weak_msun_mpc2": baseline + radial + offset,
                })
    return {
        "schema": "astrafek-hsc-y3-lens-redshift-tomographic-jackknife-v1",
        "status": "exploratory",
        "science_measurement": False,
        "radial_edges_mpc": edges,
        "profile_order": "lens_redshift_cohort lower_z then upper_z; source_zbin 1..4; radial bin; tangential and cross",
        "profile": profiles,
        "covariance": np.eye(32).tolist(),
        "lens_redshift_split": {"n_lenses_total": 20, "n_lenses_lower_z": 10,
                                 "n_lenses_upper_z": 10},
    }


def test_interaction_is_zero_when_source_bin_offsets_match():
    result = audit.analyze_interaction(_summary())
    assert result["global_interaction"]["chi2"] == pytest.approx(0.0)
    assert result["global_interaction"]["degrees_of_freedom"] == 6
    assert result["global_interaction"]["p_value"] == pytest.approx(1.0)
    assert all(row["p_value"] == pytest.approx(1.0)
               for row in result["per_radial_bin_interaction"])
    assert result["science_measurement"] is False


def test_interaction_detects_changed_relative_source_bin_offset():
    result = audit.analyze_interaction(_summary(high_shift=12.0))
    assert result["global_interaction"]["chi2"] > 0
    assert result["global_interaction"]["p_value"] < 1.0
    assert result["per_radial_bin_interaction"][0]["difference_in_differences_msun_mpc2"][1] == pytest.approx(-12.0)


def test_interaction_rejects_wrong_profile_schema():
    summary = _summary()
    summary["science_measurement"] = True
    with pytest.raises(ValueError, match="not a completed"):
        audit.analyze_interaction(summary)
