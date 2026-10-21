import json
from pathlib import Path

import pytest

import audit_wl_random_coordinate_candidate_counts as audit


def test_annulus_area_uses_arcmin_squared():
    area = audit.annulus_area_arcmin2(0.0, 60.0)
    assert area == pytest.approx(3.141592653589793)
    with pytest.raises(ValueError, match="radii"):
        audit.annulus_area_arcmin2(60.0, 60.0)


def test_run_summarizes_all_800_aggregate_counts_without_cas(tmp_path, monkeypatch):
    coordinate_plan = tmp_path / "coordinate_plan"
    coordinate_plan.mkdir()
    coordinate_manifest_path = coordinate_plan / "random_control_coordinate_plan.json"
    coordinate_manifest_path.write_text("{}\n")
    count_plan = tmp_path / "count_plan"
    (count_plan / "cas_counts").mkdir(parents=True)
    (count_plan / "random_counts.csv").write_text("approved aggregate counts\n")
    (count_plan / "random_control_count_plan.json").write_text(json.dumps({
        "candidate_query": {"candidate_target_per_lens": 200},
        "batches": [{"path": str(count_plan / "random_counts.csv")}],
    }))
    (count_plan / "cas_counts/count_status.json").write_text("{}\n")
    (count_plan / "cas_counts/count_summary.json").write_text(json.dumps({
        "status": "complete", "n_lenses": 800, "sum_candidate_counts": sum(range(200, 1000)),
    }))
    coordinate_sql = coordinate_plan / "random_coordinates_000.sql"
    coordinate_sql.write_text("coordinate-only SQL plan\n")
    lenses = [{
        "lens_id": str(index), "ra_deg": 10.0 + index / 1000,
        "dec_deg": -1.0, "z_lens": 0.25 + index * 0.0004,
        "radius_arcsec": 500.0,
    } for index in range(800)]
    ids = [row["lens_id"] for row in lenses]
    counts = dict(zip(ids, range(200, 1000)))
    manifest = {
        "count_plan": str(count_plan), "count_plan_sha256": "count-plan-sha",
        "candidate_pair_rows_expected": sum(counts.values()),
        "candidate_query": {"outer_radius_arcsec": 1800.0},
    }
    monkeypatch.setattr(audit, "validate_coordinate_plan", lambda _path: (
        manifest, [{"path": str(coordinate_sql)}], lenses, ids, counts,
    ))

    result = audit.run(coordinate_plan, tmp_path / "audit")

    assert result["n_lenses"] == 800
    assert result["sum_candidate_counts"] == sum(range(200, 1000))
    assert result["lenses_meeting_planned_candidate_count"] == 800
    assert result["zero_candidate_lenses"] == 0
    assert result["random_coordinates_read"] is False
    assert result["cas_requests_made"] is False
