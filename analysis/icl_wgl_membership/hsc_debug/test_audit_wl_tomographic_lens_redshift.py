import numpy as np
import pytest

import audit_wl_tomographic_lens_redshift as audit


def test_equal_count_split_is_deterministic_and_redshift_only():
    ids = np.array(["d", "b", "a", "c", "f", "e"])
    redshifts = np.array([0.5, 0.2, 0.1, 0.3, 0.6, 0.4])
    labels, details = audit._lens_half_labels(ids, redshifts)
    by_id = dict(zip(ids, labels))
    assert [key for key in sorted(ids, key=lambda key: (redshifts[list(ids).index(key)], key))
            if by_id[key] == 0] == ["a", "b", "c"]
    assert details["n_lenses_lower_z"] == details["n_lenses_upper_z"] == 3
    assert details["ties_at_split_boundary"] is False


def test_holm_adjustment_is_monotone_and_retains_input_order():
    adjusted = audit.holm_adjust([0.02, 0.001, 0.04, 0.5])
    assert np.allclose(adjusted, [0.06, 0.004, 0.08, 0.5])


def test_lens_redshift_jackknife_uses_shared_cells_for_both_halves(monkeypatch):
    lens_ids = np.repeat(np.array([f"L{i:02d}" for i in range(16)]), 5)
    lens_z = np.repeat(np.linspace(0.3, 0.55, 16), 5)
    source_bin = np.tile(np.arange(5), 16)
    ra_by_lens = np.repeat(np.array([10.0, 40.0, 80.0, 120.0] * 4), 5)
    pairs = {
        "lens_id": lens_ids,
        "lens_z": lens_z,
        "source_zbin": source_bin,
        "lens_ra_deg": ra_by_lens,
        "lens_dec_deg": np.zeros(80),
        "source_ra_deg": ra_by_lens + 0.01,
        "source_dec_deg": np.zeros(80),
        "object_id": np.arange(80),
    }

    def fake_stack(selected, *args, **kwargs):
        zbin = int(np.asarray(selected["source_zbin"])[0])
        mean_ra = float(np.mean(selected["lens_ra_deg"]))
        mean_z = float(np.mean(selected["lens_z"]))
        return {
            "profile": [
                {"delta_sigma_weak_msun_mpc2": mean_ra + zbin + mean_z,
                 "cross_delta_sigma_weak_msun_mpc2": 0.0},
                {"delta_sigma_weak_msun_mpc2": mean_ra / 2 + zbin + mean_z,
                 "cross_delta_sigma_weak_msun_mpc2": 0.0},
            ],
            "counts": {"used_pairs": len(selected["source_zbin"])},
        }

    monkeypatch.setattr(audit, "stack_s19a_tomographic_pairs", fake_stack)
    result = audit.lens_redshift_tomographic_jackknife(
        pairs, np.array([0.1, 0.2]), np.ones((4, 2)), [0.1, 0.2, 0.4],
        cosmology=object(), cell_width_deg=30,
    )
    assert result["n_occupied_cells"] == 4
    assert result["n_input_unassigned_source_bin_zero_excluded"] == 16
    assert result["n_input_pairs_in_tomographic_bins"] == 64
    assert result["lens_redshift_split"]["n_lenses_lower_z"] == 8
    assert result["lens_redshift_split"]["n_lenses_upper_z"] == 8
    assert len(result["profile"]) == 16
    assert len(result["tomographic_delta_sigma_consistency"]) == 4
    assert result["covariance_shape"] == [32, 32]
    assert result["science_measurement"] is False
    assert result["random_lens_subtraction_applied"] is False
    assert all(row["holm_p_value_14_tests"] is None or row["holm_p_value_14_tests"] >= row["p_value"]
               for row in result["tomographic_delta_sigma_consistency"])


def test_redshift_jackknife_rejects_inconsistent_lens_redshift():
    pairs = {
        "lens_id": np.array(["a", "a", "b", "b", "c", "c"]),
        "lens_z": np.array([0.2, 0.3, 0.4, 0.4, 0.5, 0.5]),
        "source_zbin": np.array([1, 2, 1, 2, 1, 2]),
        "lens_ra_deg": np.array([10, 10, 30, 30, 50, 50], dtype=float),
        "lens_dec_deg": np.zeros(6),
        "source_ra_deg": np.array([10.1, 10.1, 30.1, 30.1, 50.1, 50.1]),
        "source_dec_deg": np.zeros(6),
    }
    with pytest.raises(ValueError, match="inconsistent"):
        audit.lens_redshift_tomographic_jackknife(
            pairs, np.array([0.1, 0.2]), np.ones((4, 2)), [0.1, 0.2],
            cosmology=object(), cell_width_deg=10,
        )
