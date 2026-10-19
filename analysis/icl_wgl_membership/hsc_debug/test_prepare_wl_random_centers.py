import csv
import json

import pytest

from fetch_ancillary_catalogues import sha
from prepare_wl_random_centers import deterministic_select, prepare


def test_selection_is_balanced_bounded_and_reproducible():
    rows = [
        {"lens_id": "L1", "ra_deg": 10.0, "dec_deg": 1.0},
        {"lens_id": "L1", "ra_deg": 11.0, "dec_deg": 1.0},
        {"lens_id": "L1", "ra_deg": 12.0, "dec_deg": 1.0},
        {"lens_id": "L2", "ra_deg": 20.0, "dec_deg": 2.0},
    ]
    redshifts = {"L1": 0.3, "L2": 0.4, "L3": 0.5}

    selected, counts = deterministic_select(rows, redshifts, 2, "seed-A")
    repeat, repeat_counts = deterministic_select(rows, redshifts, 2, "seed-A")

    assert selected == repeat
    assert counts == repeat_counts
    assert len(selected) == 3
    assert counts == {
        "L1": {"available": 3, "selected": 2},
        "L2": {"available": 1, "selected": 1},
        "L3": {"available": 0, "selected": 0},
    }
    assert {row["z_lens"] for row in selected if row["matched_lens_id"] == "L1"} == {0.3}
    assert len({row["random_center_id"] for row in selected}) == len(selected)


def test_selection_rejects_unknown_lenses_invalid_coordinates_and_duplicates():
    redshifts = {"L1": 0.3}
    with pytest.raises(ValueError, match="unknown matched lens"):
        deterministic_select([{"lens_id": "L2", "ra_deg": 10, "dec_deg": 1}], redshifts, 1, "x")
    with pytest.raises(ValueError, match="invalid candidate coordinates"):
        deterministic_select([{"lens_id": "L1", "ra_deg": 360, "dec_deg": 1}], redshifts, 1, "x")
    with pytest.raises(ValueError, match="duplicate coordinates"):
        deterministic_select([
            {"lens_id": "L1", "ra_deg": 10, "dec_deg": 1},
            {"lens_id": "L1", "ra_deg": 10, "dec_deg": 1},
        ], redshifts, 2, "x")


def test_selection_requires_explicit_positive_quota_and_seed():
    with pytest.raises(ValueError, match="positive integer"):
        deterministic_select([], {}, 0, "seed")
    with pytest.raises(ValueError, match="non-empty"):
        deterministic_select([], {}, 1, "")


def test_prepare_writes_downstream_manifest_only_from_hashed_geometry(monkeypatch, tmp_path):
    import prepare_wl_random_centers as planner

    lens_ids = [f"L{i:03d}" for i in range(800)]
    redshifts = {lens_id: 0.3 + i / 10000 for i, lens_id in enumerate(lens_ids)}
    lenses = [{"lens_id": lens_id, "z_lens": redshifts[lens_id]} for lens_id in lens_ids]
    expected_counts = {lens_id: (2 if lens_id == "L000" else 0) for lens_id in lens_ids}

    coordinate_plan = tmp_path / "coordinates"
    cas_dir = coordinate_plan / "cas_coordinates"
    cas_dir.mkdir(parents=True)
    count_plan = tmp_path / "counts"
    count_plan.mkdir()
    preflight = tmp_path / "preflight"
    preflight.mkdir()
    source_manifest = tmp_path / "sample_manifest.json"
    source_manifest.write_text(json.dumps({"h0": 70.0, "omega_m": 0.3}))
    (count_plan / "random_control_count_plan.json").write_text(json.dumps({
        "source_preflight": str(preflight),
    }))
    (preflight / "preflight_manifest.json").write_text(json.dumps({
        "source_manifest": str(source_manifest),
        "source_sha256": sha(source_manifest),
    }))
    plan_doc = {"count_plan": str(count_plan)}
    plan_path = coordinate_plan / "random_control_coordinate_plan.json"
    plan_path.write_text(json.dumps(plan_doc))
    coordinate_summary = {
        "status": "complete", "coordinate_plan_sha256": sha(plan_path),
        "n_coordinate_rows": 2, "output_columns": ["lens_id", "ra_deg", "dec_deg"],
    }
    coordinate_summary_path = cas_dir / "coordinate_summary.json"
    coordinate_summary_path.write_text(json.dumps(coordinate_summary))
    coordinate_status = {
        "status": "complete", "credentials_saved": False,
        "source_rows_returned": False, "shape_rows_returned": False,
        "photoz_rows_returned": False, "images_returned": False,
        "science_measurement": False, "total_coordinate_rows": 2,
    }
    (cas_dir / "coordinate_status.json").write_text(json.dumps(coordinate_status))

    candidate_path = tmp_path / "candidate_geometry_audit.tsv"
    with candidate_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, ["lens_id", "ra_deg", "dec_deg",
            "geometry_candidate", "science_ready"], delimiter="\t")
        writer.writeheader()
        writer.writerow({"lens_id": "L000", "ra_deg": "10", "dec_deg": "1",
                         "geometry_candidate": True, "science_ready": False})
        writer.writerow({"lens_id": "L000", "ra_deg": "11", "dec_deg": "1",
                         "geometry_candidate": False, "science_ready": False})
    per_lens_path = tmp_path / "per_lens_candidate_counts.tsv"
    with per_lens_path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, ["lens_id", "n_candidates", "n_geometry_candidates"], delimiter="\t")
        writer.writeheader()
        writer.writerow({"lens_id": "L000", "n_candidates": 2, "n_geometry_candidates": 1})
    geometry_summary = {
        "schema": "astrafek-wl-random-coordinate-geometry-audit-v1",
        "status": "complete_geometry_screen_only", "science_measurement": False,
        "random_controls_selected": False, "coordinate_plan": str(coordinate_plan),
        "coordinate_plan_sha256": sha(plan_path),
        "coordinate_summary_sha256": sha(coordinate_summary_path),
        "candidate_geometry_audit_sha256": sha(candidate_path),
        "per_lens_candidate_counts_sha256": sha(per_lens_path),
        "n_lenses": 800, "n_coordinates": 2, "n_geometry_candidates": 1,
    }
    geometry_summary_path = tmp_path / "summary.json"
    geometry_summary_path.write_text(json.dumps(geometry_summary))
    geometry_dir = tmp_path / "geometry"
    geometry_dir.mkdir()
    for source in (candidate_path, per_lens_path, geometry_summary_path):
        (geometry_dir / source.name).write_bytes(source.read_bytes())
    geometry_summary_path = geometry_dir / "summary.json"
    geometry_summary["candidate_geometry_audit_sha256"] = sha(geometry_dir / candidate_path.name)
    geometry_summary["per_lens_candidate_counts_sha256"] = sha(geometry_dir / per_lens_path.name)
    geometry_summary_path.write_text(json.dumps(geometry_summary))

    monkeypatch.setattr(planner, "validate_coordinate_plan", lambda _: (
        plan_doc, [], lenses, lens_ids, expected_counts,
    ))
    output = tmp_path / "selected"
    result = prepare(geometry_dir, output, per_lens=1, seed="test-seed")

    manifest = json.loads((output / "random_centers.json").read_text())
    assert result["hsc_query_submitted"] is False
    assert manifest["status"] == "selected_for_shape_count_preflight"
    assert manifest["science_measurement"] is False
    assert manifest["h0_km_s_mpc"] == 70.0
    assert manifest["omega_m"] == 0.3
    assert manifest["selection"]["n_random_centers"] == 1
    assert manifest["centers"][0]["matched_lens_id"] == "L000"
