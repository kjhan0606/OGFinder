import pytest

from audit_icl_satellite_sample_retention import summarize_candidate_diagnostics, summarize_mask_support


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
