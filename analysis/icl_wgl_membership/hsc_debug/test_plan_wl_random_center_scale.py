import json

import pytest

from plan_wl_random_center_scale import estimate


def test_estimate_uses_measured_ledger_and_labels_assumptions(tmp_path):
    summary = tmp_path / "row_summary.json"
    status = tmp_path / "row_status.json"
    summary.write_text(json.dumps({
        "schema": "hsc-s19a-shape-row-summary-v1", "status": "complete",
        "science_measurement": False, "n_lenses": 10, "n_pair_rows": 1000,
    }))
    status.write_text(json.dumps({
        "status": "complete", "source_rows_downloaded": True,
        "credentials_saved": False, "total_pair_rows": 1000,
        "jobs": {"a": {"status": "complete", "rows": 400, "bytes": 120000},
                 "b": {"status": "complete", "rows": 600, "bytes": 180000}},
    }))

    result = estimate(summary, status, ratios=[1, 8], count_batch_size=25)

    assert result["status"] == "planning_estimate_only"
    assert result["network_used"] is False
    assert result["lens_export"]["mean_bytes_per_pair"] == 300
    assert result["scenarios"] == [
        {"random_centers_per_lens": 1, "random_centers_total": 10,
         "count_only_batches": 1, "expected_shape_pair_rows": 1000,
         "estimated_csv_bytes": 300000, "estimated_csv_gib": 300000 / 1024**3,
         "random_to_lens_shape_variance_ratio_if_equal_independent_centers": 1.0,
         "random_share_of_lens_plus_random_shape_variance_if_equal_independent_centers": 0.5},
        {"random_centers_per_lens": 8, "random_centers_total": 80,
         "count_only_batches": 4, "expected_shape_pair_rows": 8000,
         "estimated_csv_bytes": 2400000, "estimated_csv_gib": 2400000 / 1024**3,
         "random_to_lens_shape_variance_ratio_if_equal_independent_centers": 0.125,
         "random_share_of_lens_plus_random_shape_variance_if_equal_independent_centers": 1 / 9},
    ]
    assert "assumed equal" in " ".join(result["assumptions_and_limits"])


def test_estimate_rejects_incomplete_or_credentialled_export(tmp_path):
    summary = tmp_path / "row_summary.json"
    status = tmp_path / "row_status.json"
    summary.write_text(json.dumps({
        "schema": "hsc-s19a-shape-row-summary-v1", "status": "complete",
        "science_measurement": False, "n_lenses": 1, "n_pair_rows": 1,
    }))
    status.write_text(json.dumps({
        "status": "complete", "source_rows_downloaded": True,
        "credentials_saved": True, "total_pair_rows": 1,
        "jobs": {"a": {"status": "complete", "rows": 1, "bytes": 10}},
    }))
    with pytest.raises(ValueError, match="unexpected provenance"):
        estimate(summary, status)
