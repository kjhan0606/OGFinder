import pytest

from audit_wl_random_coordinate_reuse import summarize_rows


def test_reuse_is_allowed_between_lenses_but_duplicates_within_lens_are_excluded():
    rows = [
        {"lens_id": "L1", "ra_deg": "10", "dec_deg": "1",
         "camira_3mpc_overlap_count": "0", "geometry_candidate": "False"},
        {"lens_id": "L2", "ra_deg": "10", "dec_deg": "1",
         "camira_3mpc_overlap_count": "0", "geometry_candidate": "False"},
        {"lens_id": "L3", "ra_deg": "10", "dec_deg": "1",
         "camira_3mpc_overlap_count": "1", "geometry_candidate": "False"},
        {"lens_id": "L1", "ra_deg": "11", "dec_deg": "1",
         "camira_3mpc_overlap_count": "0", "geometry_candidate": "True"},
        {"lens_id": "L1", "ra_deg": "11", "dec_deg": "1",
         "camira_3mpc_overlap_count": "0", "geometry_candidate": "False"},
    ]

    result = summarize_rows(rows)

    assert result["n_coordinate_rows"] == 5
    assert result["strict_global_unique"]["n_candidate_rows"] == 1
    assert result["strict_global_unique"]["n_lenses_with_candidates"] == 1
    assert result["lens_scoped_reuse"]["n_candidate_rows"] == 3
    assert result["lens_scoped_reuse"]["n_lenses_with_candidates"] == 2
    assert result["lens_scoped_reuse"]["n_unique_sky_positions"] == 2
    assert result["lens_scoped_reuse"]["n_positions_shared_across_lenses"] == 1
    assert result["lens_scoped_reuse"]["n_duplicate_rows_within_lens_excluded"] == 1


def test_reuse_audit_rejects_bad_overlap_values():
    rows = [{"lens_id": "L1", "ra_deg": "10", "dec_deg": "1",
             "camira_3mpc_overlap_count": "bad", "geometry_candidate": "False"}]
    with pytest.raises(ValueError):
        summarize_rows(rows)
