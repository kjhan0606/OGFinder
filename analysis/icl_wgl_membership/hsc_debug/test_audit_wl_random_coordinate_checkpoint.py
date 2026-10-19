import hashlib
import json
from pathlib import Path

import pytest

import audit_wl_random_coordinate_checkpoint as audit_module


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _fixture(tmp_path, monkeypatch, *, status="running", wrong_rows=False):
    plan_dir = tmp_path / "plan"
    results = plan_dir / "cas_coordinates"
    results.mkdir(parents=True)
    plan_path = plan_dir / "random_control_coordinate_plan.json"
    plan_path.write_text("{}\n", encoding="utf-8")
    first_sql = plan_dir / "random_coordinates_000.sql"
    second_sql = plan_dir / "random_coordinates_001.sql"
    first_sql.write_text("SELECT 1\n", encoding="utf-8")
    second_sql.write_text("SELECT 2\n", encoding="utf-8")
    batches = [
        {"path": str(first_sql), "expected_rows": 2},
        {"path": str(second_sql), "expected_rows": 1},
    ]
    manifest = {
        "batches": batches,
        "batch_size": 1,
        "candidate_pair_rows_expected": 3,
    }
    monkeypatch.setattr(
        audit_module,
        "validate_coordinate_plan",
        lambda _plan: (manifest, batches, [{"lens_id": "a"}, {"lens_id": "b"}], ["a", "b"],
                       {"a": 2, "b": 1}),
    )
    output = results / "random_coordinates_000.csv"
    csv_rows = "lens_id,ra_deg,dec_deg\na,10,1\na,11,2\n"
    if wrong_rows:
        csv_rows = "lens_id,ra_deg,dec_deg\na,10,1\n"
    output.write_text(csv_rows, encoding="utf-8")
    checkpoint = {
        "schema": "hsc-pdr3-random-control-coordinate-run-v1",
        "status": status,
        "plan_sha256": _sha(plan_path),
        "jobs": {
            "random_coordinates_000": {
                "status": "complete", "job_id": "100", "rows": 2,
                "sha256": _sha(output),
            },
            "random_coordinates_001": {"status": "submitted", "job_id": "101"},
        },
        "credentials_saved": False,
        "source_rows_returned": False,
        "shape_rows_returned": False,
        "photoz_rows_returned": False,
        "images_returned": False,
        "science_measurement": False,
    }
    (results / "coordinate_status.json").write_text(json.dumps(checkpoint), encoding="utf-8")
    return plan_dir


def test_audit_reports_verified_partial_checkpoint_without_claiming_liveness(tmp_path, monkeypatch):
    plan_dir = _fixture(tmp_path, monkeypatch)

    result = audit_module.audit(plan_dir)

    assert result["status"] == "verified_partial"
    assert result["checkpoint_status"] == "running"
    assert result["n_batches_verified"] == 1
    assert result["n_batches_remaining"] == 1
    assert result["verified_coordinate_rows"] == 2
    assert result["process_liveness_verified"] is False
    assert result["science_measurement"] is False


def test_audit_rejects_batch_rows_that_disagree_with_approved_counts(tmp_path, monkeypatch):
    plan_dir = _fixture(tmp_path, monkeypatch, wrong_rows=True)

    with pytest.raises(ValueError, match="per-lens candidate counts"):
        audit_module.audit(plan_dir)


def test_audit_rejects_complete_checkpoint_with_incomplete_batch(tmp_path, monkeypatch):
    plan_dir = _fixture(tmp_path, monkeypatch, status="complete")

    with pytest.raises(ValueError, match="claims completion"):
        audit_module.audit(plan_dir)
