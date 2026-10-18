import numpy as np

from audit_wl_raw_shape_profiles import radial_bin_index, summarize_profiles


def test_radial_bins_include_last_edge_and_exclude_outside():
    got = radial_bin_index([0.099, 0.1, 0.2, 2.999, 3.0, 3.001])
    assert got.tolist() == [-1, 0, 1, 6, 6, -1]


def test_raw_profile_keeps_polarity_undecoded_and_rotates_components():
    pairs = {
        "lens_id": np.array(["L1"], dtype=object),
        "lens_ra_deg": np.array([150.0]),
        "lens_dec_deg": np.array([0.0]),
        "lens_z": np.array([0.3]),
        "source_ra_deg": np.array([150.01]),
        "source_dec_deg": np.array([0.0]),
        "e1": np.array([-0.1]),
        "e2": np.array([0.0]),
        "sigma_e": np.array([0.2]),
        "weight": np.array([2.0]),
        "mag_a10": np.array([24.0]),
        "resolution": np.array([0.5]),
        "source_zbin": np.array([1]),
        "b_mode_mask_raw": np.array([0]),
    }
    result = summarize_profiles(pairs)
    assert result["science_measurement"] is False
    assert result["b_mode_polarity_decoded"] is False
    assert result["delta_sigma_calculated"] is False

    def profile(mode, group="source_zbin_1"):
        return next(row for row in result["profiles"]
                    if row["b_mode_scenario"] == mode
                    and row["source_bin"] == group
                    and row["radial_bin"] == 0)

    assert profile("no_b_mode_exclusion")["n_lens_source_pairs"] == 1
    assert profile("exclude_raw_zero_hypothesis")["n_lens_source_pairs"] == 0
    kept = profile("exclude_raw_one_hypothesis")
    assert kept["n_lens_source_pairs"] == 1
    assert np.isclose(kept["mean_raw_e_t"], 0.1)
    assert np.isclose(kept["mean_raw_e_x"], 0.0)
