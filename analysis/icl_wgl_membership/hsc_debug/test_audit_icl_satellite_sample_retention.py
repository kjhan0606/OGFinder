import pytest

from audit_icl_satellite_sample_retention import (
    summarize_candidate_diagnostics, summarize_control_gate_attrition,
    summarize_mask_support, summarize_matched_association,
)


def test_mask_support_tracks_center_and_annulus_gates():
    rows = []
    for i in range(800):
        rows.append({
            "id": str(i), "primary_filter": "i" if i < 309 else "i2",
            "auto_unique_center_candidate": str(i % 2 == 0),
            "center_base_masked": str(i % 4 == 0),
            "center_masked_by_kron": "{'2.5': False, '3.5': True, '5.0': True}",
            "coverage_by_kron": "{'2.5': [0.6,0.6,0.6,0.5,0.6,0.6], '3.5': [0.6,0.6,0.6,0.6,0.6,0.6], '5.0': [0.3,0.3,0.3,0.3,0.3,0.3]}",
        })
    result = summarize_mask_support(rows)
    assert result["i"]["n_frames"] == 309
    assert result["i"]["n_unique_center_candidates"] == 155
    assert result["i"]["n_center_and_unique_unmasked_by_kron"]["2.5"] == 155
    assert result["i"]["coverage_support"]["2.5"]["0.5"]["center_eligible"] == 155
    assert result["i"]["coverage_support"]["3.5"]["0.5"]["all_six_annuli"] == 0


def test_mask_support_rejects_duplicate_or_missing_parent_clusters():
    with pytest.raises(ValueError, match="800 unique"):
        summarize_mask_support([])


def test_candidate_diagnostics_counts_control_gate_attrition():
    rows = [
        {"filter": "i", "cohort": "full_sample", "status": "target_coverage_below_floor"},
        {"filter": "i", "cohort": "full_sample", "status": "insufficient_valid_controls"},
        {"filter": "i", "cohort": "full_sample", "status": "retained"},
    ]
    result = summarize_candidate_diagnostics(rows)
    assert result["i"]["n_target_coverage_passed"] == 2
    assert result["i"]["n_target_coverage_and_control_passed"] == 1
    assert result["i"]["fraction_retained_of_candidates"] == pytest.approx(1 / 3)


def test_control_gate_audit_matches_bin_specific_outer_recovery():
    key = {"cluster_id": "12", "filter": "i2", "member_index": "4", "detection_number": "81"}
    full = {**key, "status": "insufficient_valid_controls", "candidate_control_positions": "8",
            "valid_control_profiles": "1", "valid_controls_by_bin": "4,4,3,2,2,3,3",
            "target_coverage_by_bin": "1,1,1,.9,.8,.7,.6"}
    per_bin = {**key, "status": "retained", "candidate_control_positions": "8",
               "valid_controls_by_bin": "4,4,3,2,2,3,3"}
    profile = {**key, "radial_bin_kpc": "50-75"}

    result = summarize_control_gate_attrition([full], [per_bin], [profile])

    assert result["all_annuli_valid_control_count_distribution"] == {"1": 1}
    assert result["n_with_at_least_minimum_controls_by_annulus"] == [1, 1, 1, 0, 0, 1, 1]
    assert result["outer_50_75_kpc_crosscheck"]["matched_exactly"] is True
    assert result["outer_50_75_kpc_crosscheck"]["n_matching_saved_per_bin_profiles"] == 1


def test_control_gate_audit_rejects_changed_control_geometry():
    key = {"cluster_id": "12", "filter": "i", "member_index": "4", "detection_number": "81"}
    full = {**key, "status": "insufficient_valid_controls", "candidate_control_positions": "8",
            "valid_control_profiles": "1", "valid_controls_by_bin": "3,3,3,3,3,3,2",
            "target_coverage_by_bin": "1,1,1,1,1,1,1"}
    per_bin = {**key, "status": "retained", "candidate_control_positions": "8",
               "valid_controls_by_bin": "3,3,3,3,3,3,3"}
    with pytest.raises(ValueError, match="control geometry changed"):
        summarize_control_gate_attrition([full], [per_bin], [])


def test_matched_association_statistics_are_derived_from_results():
    def row(band, n, blocks, r, spatial_ci, ordinary_ci, rho, rho_ci):
        return {
            "filter": band, "selection": "all_supported_pairs", "status": "complete_exploratory",
            "n_pairs": n, "n_connected_5deg_blocks": blocks,
            "partial_pearson_r_control_redshift_richness": r, "spearman_rho": rho,
            "spearman_spatial_block_ci95": rho_ci,
            "partial_spatial_block_bootstrap": {
                "status": "complete_exploratory", "n_pairs": n, "n_blocks": blocks,
                "partial_r_point": r, "ci95_partial_r": spatial_ci,
            },
            "partial_ordinary_pair_bootstrap": {
                "status": "complete_exploratory", "partial_r_point": r,
                "ci95_partial_r": ordinary_ci,
            },
        }

    result = summarize_matched_association([
        row("i", 26, 13, 0.08, [-0.19, 0.48], [-0.32, 0.56], 0.27, [-0.13, 0.54]),
        row("i2", 68, 33, 0.29, [0.04, 0.50], [-0.01, 0.53], 0.11, [-0.09, 0.30]),
    ])

    assert result["n_cluster_filter_pairs_by_filter"] == {"i": 26, "i2": 68}
    assert result["by_filter"]["i"]["partial_r_control_redshift_richness"] == pytest.approx(0.08)
    assert result["partial_r_i2"] == pytest.approx(0.29)
    assert result["ordinary_pair_ci95_i2"] == pytest.approx([-0.01, 0.53])


def test_matched_association_rejects_incomplete_or_inconsistent_result_rows():
    with pytest.raises(ValueError, match="exactly one"):
        summarize_matched_association([])
    invalid = {
        "filter": "i", "selection": "all_supported_pairs", "status": "complete_exploratory",
        "n_pairs": 26, "n_connected_5deg_blocks": 13,
        "partial_pearson_r_control_redshift_richness": 0.2, "spearman_rho": 0.1,
        "spearman_spatial_block_ci95": [-0.1, 0.2],
        "partial_spatial_block_bootstrap": {
            "status": "complete_exploratory", "n_pairs": 26, "n_blocks": 13,
            "partial_r_point": 0.2, "ci95_partial_r": [-0.1, 0.2],
        },
        "partial_ordinary_pair_bootstrap": {
            "status": "complete_exploratory", "partial_r_point": 0.3,
            "ci95_partial_r": [-0.2, 0.4],
        },
    }
    with pytest.raises(ValueError, match="disagree"):
        summarize_matched_association([invalid, {**invalid, "filter": "i2"}])
