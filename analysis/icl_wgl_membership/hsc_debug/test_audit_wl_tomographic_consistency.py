import numpy as np
import pytest

import audit_wl_tomographic_consistency as audit


def test_common_value_gls_handles_consistent_bin_values():
    result = audit.fit_common_value([2.0, 2.0, 2.0, 2.0], np.eye(4))
    assert result["common_value"] == 2.0
    assert result["chi2"] == 0.0
    assert result["degrees_of_freedom"] == 3
    assert result["p_value"] == 1.0


def test_tomographic_jackknife_reuses_common_sky_deletions(monkeypatch):
    pairs = {
        "source_zbin": np.repeat(np.arange(5), 3),
        "lens_ra_deg": np.tile([10.0, 20.0, 30.0], 5),
        "lens_dec_deg": np.zeros(15),
        "source_ra_deg": np.tile([10.01, 20.01, 30.01], 5),
        "source_dec_deg": np.zeros(15),
    }

    def fake_stack(selected, *args, **kwargs):
        n = len(selected["source_zbin"])
        mean_ra = float(np.mean(selected["lens_ra_deg"]))
        return {
            "profile": [
                {"delta_sigma_weak_msun_mpc2": mean_ra, "cross_delta_sigma_weak_msun_mpc2": 0.0},
                {"delta_sigma_weak_msun_mpc2": mean_ra / 2, "cross_delta_sigma_weak_msun_mpc2": 0.0},
            ],
            "counts": {"used_pairs": n},
        }

    monkeypatch.setattr(audit, "stack_s19a_tomographic_pairs", fake_stack)
    result = audit.tomographic_spatial_jackknife(
        pairs, np.array([0.1, 0.2]), np.ones((4, 2)), [0.1, 0.2, 0.4],
        cosmology=object(), cell_width_deg=10,
    )
    assert result["n_occupied_cells"] == 3
    assert len(result["profile"]) == 8
    assert len(result["tomographic_delta_sigma_consistency"]) == 2
    assert len(result["covariance"]) == 16
    assert result["source_bin_counts"][0]["used_pairs"] == 3
    assert result["n_input_unassigned_source_bin_zero_excluded"] == 3
    assert result["n_input_pairs_in_tomographic_bins"] == 12
    assert result["science_measurement"] is False


def test_tomographic_jackknife_rejects_bin_above_four():
    pairs = {
        "source_zbin": np.array([1, 2, 3, 4, 5]),
        "lens_ra_deg": np.array([10.0, 20.0, 30.0, 10.0, 20.0]),
        "lens_dec_deg": np.zeros(5),
        "source_ra_deg": np.array([10.01, 20.01, 30.01, 10.01, 20.01]),
        "source_dec_deg": np.zeros(5),
    }
    with pytest.raises(ValueError, match="integer values 0..4"):
        audit.tomographic_spatial_jackknife(
            pairs, np.array([0.1, 0.2]), np.ones((4, 2)), [0.1, 0.2],
            cosmology=object(), cell_width_deg=10,
        )
