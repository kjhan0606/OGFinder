import csv
import json

import pytest
import numpy as np

from prepare_wl_random_shape_count_plan import prepare as prepare_count_plan
from prepare_wl_random_shape_rows import (
    OUTPUT_COLUMNS, prepare, shape_rows_query,
)
from run_wl_random_shape_stack import verify_and_load_pairs
from submit_wl_random_shape_counts import sha
from submit_wl_random_shape_rows import read_shape_rows, run, validate_plan


def _centers():
    return {
        "schema": "astrafek-wl-random-centers-v1",
        "status": "selected_for_shape_count_preflight",
        "science_measurement": False,
        "h0_km_s_mpc": 70.0,
        "omega_m": 0.3,
        "centers": [
            {"random_center_id": "R-001", "matched_lens_id": "L-1",
             "ra_deg": 12.0, "dec_deg": -2.0, "z_lens": 0.31},
            {"random_center_id": "R-002", "matched_lens_id": "L-1",
             "ra_deg": 13.0, "dec_deg": -1.0, "z_lens": 0.31},
            {"random_center_id": "R-003", "matched_lens_id": "L-2",
             "ra_deg": 30.0, "dec_deg": 1.0, "z_lens": 0.42},
        ],
    }


def _make_count_inputs(tmp_path):
    centers_path = tmp_path / "centers.json"
    centers_path.write_text(json.dumps(_centers()))
    count_plan_dir = tmp_path / "count-plan"
    count_manifest = prepare_count_plan(centers_path, count_plan_dir, batch_size=2)
    plan_path = count_plan_dir / "random_shape_count_plan.json"
    plan_sha = sha(plan_path)
    results = tmp_path / "count-results"
    results.mkdir()
    center_rows = _centers()["centers"]
    values = {"R-001": 2, "R-002": 1, "R-003": 1}
    jobs = {}
    for index, start in enumerate(range(0, len(center_rows), 2)):
        name = f"random_shape_counts_{index:03d}"
        path = results / f"{name}.csv"
        chunk = center_rows[start:start + 2]
        with path.open("w", encoding="utf-8", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=[
                "random_center_id", "matched_lens_id", "n_shape_sources"])
            writer.writeheader()
            for center in chunk:
                writer.writerow({"random_center_id": center["random_center_id"],
                                 "matched_lens_id": center["matched_lens_id"],
                                 "n_shape_sources": values[center["random_center_id"]]})
        jobs[name] = {"status": "complete", "rows": len(chunk), "sha256": sha(path)}
    summary = {
        "status": "complete", "plan_sha256": plan_sha,
        "n_random_centers": 3, "n_matched_lenses": 2,
        "sum_shape_source_counts": 4,
        "source_shape_rows_returned": False, "photoz_rows_returned": False,
        "images_returned": False, "science_measurement": False,
    }
    (results / "count_summary.json").write_text(json.dumps(summary))
    state = {
        "status": "complete", "plan_sha256": plan_sha, "jobs": jobs,
        "credentials_saved": False, "source_shape_rows_returned": False,
        "photoz_rows_returned": False, "images_returned": False,
        "science_measurement": False,
    }
    (results / "count_status.json").write_text(json.dumps(state))
    return count_plan_dir, results, count_manifest


def _shape_row(center, object_id="O-1"):
    row = {name: "0" for name in OUTPUT_COLUMNS}
    row.update({
        "random_center_id": center["random_center_id"],
        "matched_lens_id": center["matched_lens_id"],
        "random_ra_deg": str(center["ra_deg"]),
        "random_dec_deg": str(center["dec_deg"]),
        "lens_z": str(center["z_lens"]),
        "object_id": object_id,
        "source_ra_deg": "45.0", "source_dec_deg": "10.0",
        "weight": "1.0", "rms_e": "0.2", "b_mode_mask_raw": "1",
        "source_zbin": "2",
    })
    return row


def test_random_shape_row_sql_matches_count_query_aperture_and_explicit_columns():
    center = {"random_center_id": "R'1", "matched_lens_id": "L'2",
              "ra_deg": 12.125, "dec_deg": -2.5, "z_lens": 0.31,
              "radius_arcsec": 900.0}
    sql = shape_rows_query(center)
    assert "SELECT DISTINCT" in sql
    assert "'R''1' AS random_center_id" in sql
    assert "'L''2' AS matched_lens_id" in sql
    assert "coneSearch(c.i_coord, 12.125, -2.5, 900)" in sql
    assert "b.b_mode_mask AS b_mode_mask_raw" in sql
    assert "b.hsc_y3_zbin AS source_zbin" in sql
    assert "photoz" not in sql.lower()
    assert "SELECT *" not in sql


def test_prepare_and_validate_plan_are_bound_to_verified_count_rows(tmp_path):
    count_plan, count_results, _ = _make_count_inputs(tmp_path)
    row_plan_dir = tmp_path / "row-plan"
    manifest = prepare(count_plan, count_results, row_plan_dir)
    validated, centers, batches, counts, count_plan_sha, row_plan_sha = validate_plan(row_plan_dir)

    assert manifest["status"] == "prepared_not_submitted"
    assert manifest["n_random_centers"] == 3
    assert manifest["n_query_batches"] == 2
    assert manifest["expected_random_source_pair_rows"] == 4
    assert manifest["output_columns"] == OUTPUT_COLUMNS
    assert manifest["source_shape_rows_returned"] is False
    assert manifest["photoz_pdf_rows_included"] is False
    assert manifest["images_returned"] is False
    assert manifest["execution_requires_explicit_user_approval"] is True
    assert len(centers) == 3 and len(batches) == 2
    assert counts["R-001"]["n_shape_sources"] == 2
    assert count_plan_sha == sha(count_plan / "random_shape_count_plan.json")
    assert len(row_plan_sha) == 64
    assert validated["expected_random_source_pair_rows"] == 4


def test_row_submitter_requires_separate_explicit_approval(tmp_path):
    count_plan, count_results, _ = _make_count_inputs(tmp_path)
    row_plan_dir = tmp_path / "row-plan"
    prepare(count_plan, count_results, row_plan_dir)
    results_dir = tmp_path / "row-results"

    with pytest.raises(ValueError, match="explicit user approval"):
        run(row_plan_dir, results_dir, "test-user", submit=True, approved=False)

    assert not results_dir.exists()


def test_shape_row_response_checks_center_mapping_and_pair_uniqueness(tmp_path):
    center = {"random_center_id": "R-001", "matched_lens_id": "L-1",
              "ra_deg": 12.0, "dec_deg": -2.0, "z_lens": 0.31}
    path = tmp_path / "rows.csv"
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=OUTPUT_COLUMNS)
        writer.writeheader()
        writer.writerow(_shape_row(center, "O-1"))
        writer.writerow(_shape_row(center, "O-2"))
    assert read_shape_rows(path, [center], {"R-001": {"n_shape_sources": 2}}) == {"R-001": 2}

    with path.open("a", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=OUTPUT_COLUMNS)
        writer.writerow(_shape_row(center, "O-2"))
    with pytest.raises(ValueError, match="duplicated"):
        read_shape_rows(path, [center], {"R-001": {"n_shape_sources": 3}})


def test_row_response_rejects_changed_lens_mapping(tmp_path):
    center = {"random_center_id": "R-001", "matched_lens_id": "L-1",
              "ra_deg": 12.0, "dec_deg": -2.0, "z_lens": 0.31}
    row = _shape_row(center)
    row["matched_lens_id"] = "L-2"
    path = tmp_path / "rows.csv"
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=OUTPUT_COLUMNS)
        writer.writeheader()
        writer.writerow(row)
    with pytest.raises(ValueError, match="association"):
        read_shape_rows(path, [center], {"R-001": {"n_shape_sources": 1}})


def test_stack_loader_requires_complete_hashed_rows_and_preserves_center_mapping(tmp_path):
    count_plan, count_results, _ = _make_count_inputs(tmp_path)
    row_plan_dir = tmp_path / "row-plan"
    manifest = prepare(count_plan, count_results, row_plan_dir)
    row_manifest_path = row_plan_dir / "random_shape_row_plan.json"
    rows_dir = tmp_path / "row-results"
    rows_dir.mkdir()
    centers = _centers()["centers"]
    expected = {"R-001": 2, "R-002": 1, "R-003": 1}
    jobs = {}
    for index, start in enumerate(range(0, len(centers), manifest["batch_size"])):
        name = f"random_shape_rows_{index:03d}"
        path = rows_dir / f"{name}.csv"
        chunk = centers[start:start + manifest["batch_size"]]
        with path.open("w", encoding="utf-8", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=OUTPUT_COLUMNS)
            writer.writeheader()
            for center in chunk:
                for source_index in range(expected[center["random_center_id"]]):
                    writer.writerow(_shape_row(center, f"{center['random_center_id']}-O{source_index}"))
        jobs[name] = {
            "status": "complete", "rows": sum(expected[c["random_center_id"]] for c in chunk),
            "sha256": sha(path),
        }
    plan_sha = sha(row_manifest_path)
    count_plan_sha = manifest["source_count_plan_sha256"]
    summary = {
        "status": "complete", "plan_sha256": plan_sha,
        "n_source_pairs": sum(expected.values()), "source_shape_rows_downloaded": True,
        "photoz_pdf_rows_included": False, "images_returned": False,
        "science_measurement": False,
    }
    state = {
        "schema": "hsc-s19a-random-shape-row-run-v1", "status": "complete",
        "plan_sha256": plan_sha, "count_plan_sha256": count_plan_sha,
        "source_shape_rows_downloaded": True, "credentials_saved": False,
        "photoz_pdf_rows_included": False, "images_returned": False,
        "science_measurement": False, "jobs": jobs,
    }
    (rows_dir / "row_status.json").write_text(json.dumps(state))
    (rows_dir / "row_summary.json").write_text(json.dumps(summary))

    pairs, provenance = verify_and_load_pairs(rows_dir, row_plan_dir, masked_value=0)

    assert len(pairs["object_id"]) == 4
    assert set(pairs["lens_id"]) == {"R-001", "R-002", "R-003"}
    assert np.all(pairs["b_mode_pass"])
    assert provenance["n_source_pairs"] == 4
    assert provenance["n_random_centers"] == 3
