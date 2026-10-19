import numpy as np
import pytest

from validate_wl_source_reuse_covariance import validate_summary


def _payload():
    edges = [0.1, 0.2, 0.4, 0.7, 1.0, 1.5, 2.0, 3.0]
    dim = 56
    modes = {}
    for mode_name in ("lens_only", "source_only", "joint_endpoints"):
        covariance = np.eye(dim).tolist()
        profile = []
        vector = np.arange(dim, dtype=float)
        for source_bin in range(4):
            for radial_bin in range(7):
                pos = source_bin * 14 + radial_bin * 2
                profile.append({
                    "source_zbin": source_bin + 1,
                    "r_min_mpc": edges[radial_bin],
                    "r_max_mpc": edges[radial_bin + 1],
                    "delta_sigma_weak_msun_mpc2": float(vector[pos]),
                    "cross_delta_sigma_weak_msun_mpc2": float(vector[pos + 1]),
                    "delta_sigma_jackknife_error_msun_mpc2": 1.0,
                    "cross_jackknife_error_msun_mpc2": 1.0,
                })
        consistency = []
        for radial_bin in range(7):
            values = [float(vector[source_bin * 14 + radial_bin * 2])
                      for source_bin in range(4)]
            # Identity covariance gives common mean and chi-square with 3 dof.
            common = float(np.mean(values))
            chi_square = float(sum((value - common) ** 2 for value in values))
            from scipy.stats import chi2
            consistency.append({
                "radial_bin": radial_bin,
                "r_min_mpc": edges[radial_bin],
                "r_max_mpc": edges[radial_bin + 1],
                "source_bin_delta_sigma_msun_mpc2": values,
                "common_value": common,
                "chi2": chi_square,
                "degrees_of_freedom": 3,
                "p_value": float(chi2.sf(chi_square, 3)),
                "covariance_rank": 4,
                "holm_p_value_7_radial_tests": 0.0,
            })
        raw = np.asarray([row["p_value"] for row in consistency])
        order = np.argsort(raw)
        adjusted = np.empty_like(raw)
        running = 0.0
        for rank, index in enumerate(order):
            running = max(running, (len(raw) - rank) * raw[index])
            adjusted[index] = min(1.0, running)
        for row, value in zip(consistency, adjusted):
            row["holm_p_value_7_radial_tests"] = float(value)
        modes[mode_name] = {
            "n_deletion_blocks": 3,
            "covariance_dimension": dim,
            "covariance_rank": dim,
            "relative_covariance_symmetry_error": 0.0,
            "profile": profile,
            "tomographic_consistency": consistency,
            "covariance": covariance,
            "block_deletions": [
                {"block_id": index, "deleted_pair_rows": 1,
                 "used_pairs_remaining_by_source_bin": [10, 10, 10, 10]}
                for index in range(3)
            ],
        }
    return {
        "schema": "astrafek-hsc-y3-source-reuse-covariance-audit-v1",
        "status": "verified_covariance_sensitivity",
        "science_measurement": False,
        "random_lens_subtraction_applied": False,
        "input": {"radial_edges_mpc": edges},
        "source_reuse": {
            "selected_pair_rows": 6,
            "unique_source_objects": 3,
            "source_objects_reused_across_lens_pairs": 2,
            "fraction_unique_sources_reused": 2 / 3,
            "pair_rows_from_reused_sources": 5,
            "fraction_pair_rows_from_reused_sources": 5 / 6,
            "maximum_lens_pair_multiplicity_per_source": 3,
            "reused_sources_spanning_multiple_lens_endpoint_cells": 1,
            "reused_sources_spanning_multiple_source_endpoint_cells": 0,
            "source_endpoint_cell_assignment_is_unique_per_object_id": True,
        },
        "covariance_modes": modes,
    }


def test_validates_covariance_summary_internal_invariants():
    result = validate_summary(_payload())
    assert result["status"] == "verified"
    assert result["science_measurement"] is False
    assert result["covariance_modes"]["joint_endpoints"]["profile_rows_checked"] == 28
    assert result["covariance_modes"]["joint_endpoints"][
        "minimum_covariance_eigenvalue_relative_to_entry_scale"] == 1.0
    assert result["covariance_modes"]["joint_endpoints"][
        "psd_tolerance_relative_to_entry_scale"] == 1e-12


def test_rejects_profile_error_that_disagrees_with_covariance():
    payload = _payload()
    payload["covariance_modes"]["lens_only"]["profile"][0][
        "delta_sigma_jackknife_error_msun_mpc2"] = 2.0
    with pytest.raises(ValueError, match="jackknife_error"):
        validate_summary(payload)


def test_rejects_non_psd_covariance():
    payload = _payload()
    payload["covariance_modes"]["source_only"]["covariance"][0][0] = -1.0
    with pytest.raises(ValueError, match="positive semidefinite"):
        validate_summary(payload)
