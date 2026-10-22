import hashlib
from pathlib import Path

import pytest

from audit_wl_random_coordinate_partial_export import audit_checkpoint


def fixture(tmp_path):
    results = tmp_path / "cas_coordinates"
    results.mkdir(parents=True)
    csv_path = results / "random_coordinates_000.csv"
    csv_path.write_text(
        "lens_id,ra_deg,dec_deg\nL1,10.0,-1.0\nL1,10.1,-1.1\n",
        encoding="utf-8",
    )
    digest = hashlib.sha256(csv_path.read_bytes()).hexdigest()
    manifest = {"batch_size": 1, "candidate_pair_rows_expected": 3}
    batches = [
        {"path": "/plan/random_coordinates_000.sql", "expected_rows": 2},
        {"path": "/plan/random_coordinates_001.sql", "expected_rows": 1},
    ]
    state = {
        "schema": "hsc-pdr3-random-control-coordinate-run-v1",
        "status": "running", "plan_sha256": "plan-hash",
        "credentials_saved": False, "source_rows_returned": False,
        "shape_rows_returned": False, "photoz_rows_returned": False,
        "images_returned": False, "science_measurement": False,
        "jobs": {
            "random_coordinates_000": {
                "status": "complete", "job_id": "100", "rows": 2,
                "sha256": digest,
            },
            "random_coordinates_001": {"status": "submitted", "job_id": "101"},
        },
    }
    counts = {"L1": 2, "L2": 1}
    return manifest, batches, ["L1", "L2"], counts, state, results


def test_audit_verifies_saved_batches_and_leaves_submitted_job_unresolved(tmp_path):
    manifest, batches, lens_ids, counts, state, results = fixture(tmp_path)

    result = audit_checkpoint(manifest, batches, lens_ids, counts, state,
                              results, "plan-hash")

    assert result["status"] == "partial_verified"
    assert result["verified_coordinate_rows"] == 2
    assert result["unresolved_batches"] == [{
        "batch": "random_coordinates_001", "status": "submitted",
        "job_id": "101", "saved_csv_present": False,
    }]
    assert result["network_used"] is False
    assert result["rows_downloaded"] is False


def test_audit_rejects_corrupt_csv_or_changed_checkpoint_scope(tmp_path):
    manifest, batches, lens_ids, counts, state, results = fixture(tmp_path)
    (results / "random_coordinates_000.csv").write_text("corrupt\n", encoding="utf-8")
    with pytest.raises(ValueError, match="checksum mismatch"):
        audit_checkpoint(manifest, batches, lens_ids, counts, state,
                         results, "plan-hash")

    manifest, batches, lens_ids, counts, state, results = fixture(tmp_path / "other")
    state["source_rows_returned"] = True
    with pytest.raises(ValueError, match="outside the coordinate-only scope"):
        audit_checkpoint(manifest, batches, lens_ids, counts, state,
                         results, "plan-hash")


def test_audit_rejects_checkpoint_bound_to_another_plan(tmp_path):
    manifest, batches, lens_ids, counts, state, results = fixture(tmp_path)

    with pytest.raises(ValueError, match="different coordinate plan"):
        audit_checkpoint(manifest, batches, lens_ids, counts, state,
                         results, "new-plan-hash")
