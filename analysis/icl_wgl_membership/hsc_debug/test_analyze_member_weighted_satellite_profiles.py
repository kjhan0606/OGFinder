import hashlib
import json

import pytest

from analyze_member_weighted_satellite_profiles import aggregate, validate_member_metadata


def profile(cluster, member, delta, radial="50-75"):
    return {
        "cluster_id": cluster,
        "filter": "i2",
        "member_index": str(member),
        "radial_bin_kpc": radial,
        "member_detection_separation_arcsec": "0.2",
        "delta_sb_njy_arcsec2": str(delta),
    }


def test_mass_bins_and_membership_weight_are_cluster_equal_weighted():
    rows = [profile("A", 0, 1), profile("A", 1, 3), profile("B", 0, 5)]
    crosswalk = {
        ("A", 0): {"candidate_log10_mstar": 10.7, "camira_weight": 0.25},
        ("A", 1): {"candidate_log10_mstar": 10.8, "camira_weight": 0.75},
        ("B", 0): {"candidate_log10_mstar": 10.9, "camira_weight": 1.0},
    }
    _, result = aggregate(rows, crosswalk)
    all_row = next(r for r in result if r["mass_selection"] == "all")
    mid_row = next(r for r in result if r["mass_selection"] == "10.5-11.0")
    assert all_row["n_clusters"] == 2
    assert all_row["n_satellites"] == 3
    assert float(all_row["equal_cluster_mean_delta_njy_arcsec2"]) == pytest.approx(3.5)
    # Cluster A's weighted mean is 2.5; cluster B's is 5. Equal cluster weighting -> 3.75.
    assert float(all_row["camira_weighted_equal_cluster_mean_delta_njy_arcsec2"]) == pytest.approx(3.75)
    assert mid_row["n_clusters"] == 2
    assert mid_row["n_satellites"] == 3


def test_duplicate_satellite_annulus_row_is_rejected():
    row = profile("A", 0, 1)
    crosswalk = {("A", 0): {"candidate_log10_mstar": 10.7,
                             "camira_weight": 0.5}}
    with pytest.raises(ValueError, match="Duplicate profile row"):
        aggregate([row, row], crosswalk)


def test_profile_without_member_crosswalk_is_rejected():
    with pytest.raises(ValueError, match="No member properties"):
        aggregate([profile("A", 0, 1)], {})


def test_member_metadata_must_match_source_hash_and_column_map(tmp_path):
    member_file = tmp_path / "camira_s20a_wide_sm_member.dat.gz"
    member_file.write_bytes(b"catalog bytes")
    summary_file = tmp_path / "member_summary.json"
    summary = {
        "source_file": member_file.name,
        "sha256": hashlib.sha256(member_file.read_bytes()).hexdigest(),
        "columns": ["ra_cl_deg", "dec_cl_deg", "N_mem", "z_cl",
                    "ra_deg", "dec_deg", "field7", "weight"],
        "field7_note": "undocumented field",
    }
    summary_file.write_text(json.dumps(summary))

    assert validate_member_metadata(summary_file, member_file)["columns"][-1] == "weight"

    summary["sha256"] = "0" * 64
    summary_file.write_text(json.dumps(summary))
    with pytest.raises(ValueError, match="does not match the member-file SHA-256"):
        validate_member_metadata(summary_file, member_file)
