import json
from pathlib import Path

import pytest

import prepare_hsc_psf_picker_target_control_plan as planner


def _write_json(path: Path, value):
    path.write_text(json.dumps(value) + "\n", encoding="utf-8")
    return path


def test_plan_deduplicates_centers_and_excludes_previously_approved_positions(tmp_path):
    target_manifest = _write_json(tmp_path / "targets.json", {"clusters": [
        {"id": "T1", "ra_deg": 10.0, "dec_deg": 1.0},
        {"id": "T2", "ra_deg": 20.0, "dec_deg": 2.0},
    ]})
    control_manifest = _write_json(tmp_path / "controls.json", {"clusters": [
        {"id": "C1", "target_cluster_id": "T1", "ra_deg": 30.0, "dec_deg": 3.0},
        {"id": "C2", "target_cluster_id": "T2", "ra_deg": 20.0, "dec_deg": 2.0},
    ]})
    approved_mapping = _write_json(tmp_path / "approved.json", [
        {"ra_deg": 10.0, "dec_deg": 1.0},
    ])
    prior_audit = _write_json(tmp_path / "audit.json", {"n_returned": 2, "total_bytes": 30})
    prior_psf_dir = tmp_path / "psfs"
    prior_psf_dir.mkdir()
    (prior_psf_dir / "a.fits").write_bytes(b"a" * 10)
    (prior_psf_dir / "b.fits").write_bytes(b"b" * 20)
    output = tmp_path / "plan"

    plan = planner.prepare(
        target_manifest, control_manifest, approved_mapping, prior_audit,
        prior_psf_dir, output, expected_targets=2, expected_controls=2,
    )

    assert plan["status"] == "prepared_not_submitted"
    assert plan["n_target_centers"] == 2
    assert plan["n_control_centers"] == 2
    assert plan["n_unique_combined_centers"] == 3
    assert plan["n_already_covered_by_approved_psfs"] == 1
    assert plan["n_requests_remaining"] == 2
    assert plan["estimated_total_bytes_median"] == 30
    assert plan["estimated_total_bytes_mean"] == 30
    assert plan["network_used"] is False
    assert plan["request_submitted"] is False
    assert plan["images_requested"] is False
    assert plan["shape_rows_requested"] is False
    assert plan["photoz_requested"] is False
    rows = json.loads((output / "position_mapping.json").read_text())
    overlap = next(row for row in rows if row["ra_deg"] == 20.0)
    assert {role["role"] for role in overlap["roles"]} == {
        "target_cluster_center", "matched_control_center",
    }
    coordinates = (output / "psf_picker_coordinates.txt").read_text().splitlines()
    assert len(coordinates) == 3  # header plus two uncovered positions
    with pytest.raises(FileExistsError, match="refusing to overwrite"):
        planner.prepare(
            target_manifest, control_manifest, approved_mapping, prior_audit,
            prior_psf_dir, output, expected_targets=2, expected_controls=2,
        )


def test_plan_rejects_audit_mismatch_and_bad_coordinates(tmp_path):
    target_manifest = _write_json(tmp_path / "targets.json", {"clusters": [
        {"id": "T1", "ra_deg": 400.0, "dec_deg": 0.0},
    ]})
    control_manifest = _write_json(tmp_path / "controls.json", {"clusters": []})
    approved_mapping = _write_json(tmp_path / "approved.json", [])
    prior_audit = _write_json(tmp_path / "audit.json", {"n_returned": 1, "total_bytes": 1})
    prior_psf_dir = tmp_path / "psfs"
    prior_psf_dir.mkdir()
    (prior_psf_dir / "a.fits").write_bytes(b"a")

    with pytest.raises(ValueError, match="unexpected center counts"):
        planner.prepare(
            target_manifest, control_manifest, approved_mapping, prior_audit,
            prior_psf_dir, tmp_path / "wrong-counts", expected_targets=2, expected_controls=0,
        )
    with pytest.raises(ValueError, match="celestial coordinates"):
        planner.prepare(
            target_manifest, control_manifest, approved_mapping, prior_audit,
            prior_psf_dir, tmp_path / "bad-coordinate", expected_targets=1, expected_controls=0,
        )
    valid_targets = _write_json(tmp_path / "valid_targets.json", {"clusters": [
        {"id": "T1", "ra_deg": 10.0, "dec_deg": 0.0},
    ]})
    mismatched_audit = _write_json(tmp_path / "mismatched_audit.json", {
        "n_returned": 1, "total_bytes": 2,
    })
    with pytest.raises(ValueError, match="disagree with the saved audit"):
        planner.prepare(
            valid_targets, control_manifest, approved_mapping, mismatched_audit,
            prior_psf_dir, tmp_path / "bad-audit", expected_targets=1, expected_controls=0,
        )
