import numpy as np
import pytest
from scipy.stats import chi2

from audit_wl_cross_shear_null_decomposition import audit, holm


def _summary():
    profiles = []
    # Four source bins, two radii, with profile order: bin, radius, t/cross.
    for source_bin in (1, 2, 3, 4):
        for radial_index, (r_min, r_max) in enumerate(((0.1, 0.2), (0.2, 0.4))):
            profiles.append({
                "source_zbin": source_bin,
                "r_min_mpc": r_min,
                "r_max_mpc": r_max,
                "cross_delta_sigma_weak_msun_mpc2": float(source_bin + radial_index),
                "cross_jackknife_error_msun_mpc2": 1.0,
            })
    covariance = np.eye(16).tolist()
    cross = np.asarray([1.0, 2.0, 2.0, 3.0, 3.0, 4.0, 4.0, 5.0])
    return {
        "schema": "astrafek-hsc-y3-tomographic-consistency-jackknife-v1",
        "science_measurement": False,
        "profile": profiles,
        "covariance": covariance,
        "n_occupied_cells": 5,
        "joint_cross_null": {"chi2": float(cross @ cross), "degrees_of_freedom": 8,
                             "p_value": float(chi2.sf(cross @ cross, 8))},
    }


def test_holm_preserves_order_and_is_monotonic_in_sorted_pvalues():
    adjusted = holm([0.03, 0.01, 0.2])
    assert adjusted == pytest.approx([0.06, 0.03, 0.2])


def test_audit_reproduces_global_null_and_adjusts_radial_family():
    result = audit(_summary(), "a" * 64)
    assert result["science_measurement"] is False
    assert result["global_cross_null"]["chi2"] == pytest.approx(84.0)
    assert len(result["radial_bin_joint_tests"]) == 2
    assert all("p_holm_radial_bins" in row for row in result["radial_bin_joint_tests"])


def test_audit_rejects_inconsistent_covariance_or_science_flag():
    summary = _summary()
    summary["science_measurement"] = True
    with pytest.raises(ValueError, match="non-science"):
        audit(summary, "a" * 64)

    summary = _summary()
    summary["joint_cross_null"]["chi2"] = 99
    with pytest.raises(ValueError, match="disagrees"):
        audit(summary, "a" * 64)
