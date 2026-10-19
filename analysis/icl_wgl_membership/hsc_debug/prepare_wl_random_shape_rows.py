#!/usr/bin/env python3
"""Prepare (without submission) minimal HSC shape-row SQL for matched random centers."""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

from prepare_wl_random_shape_count_plan import sha
from submit_wl_random_shape_counts import read_count_rows, validate_plan as validate_count_plan


OUTPUT_COLUMNS = [
    "random_center_id", "matched_lens_id", "random_ra_deg", "random_dec_deg", "lens_z",
    "object_id", "source_ra_deg", "source_dec_deg", "e1", "e2", "sigma_e", "rms_e",
    "weight", "m", "c1", "c2", "mag_a10", "resolution", "psf_shape11",
    "psf_shape22", "psf_shape12", "b_mode_mask_raw", "source_zbin",
]
DEFAULT_ROW_BYTES = {"low": 256, "high": 480}


def sql_literal(value: str) -> str:
    value = str(value)
    if not value or any(ord(char) < 32 for char in value):
        raise ValueError("IDs must be non-empty printable strings")
    return "'" + value.replace("'", "''") + "'"


def sql_number(value: float) -> str:
    number = float(value)
    if not math.isfinite(number):
        raise ValueError("random center coordinates and redshifts must be finite")
    return format(number, ".12g")


def shape_rows_query(center: dict[str, str | float]) -> str:
    """Return the same minimal 23-column shape schema used by the lens stack."""
    center_id = sql_literal(center["random_center_id"])
    matched_id = sql_literal(center["matched_lens_id"])
    ra, dec = sql_number(center["ra_deg"]), sql_number(center["dec_deg"])
    redshift, radius = sql_number(center["z_lens"]), sql_number(center["radius_arcsec"])
    projection = [
        f"{center_id} AS random_center_id", f"{matched_id} AS matched_lens_id",
        f"{ra} AS random_ra_deg", f"{dec} AS random_dec_deg", f"{redshift} AS lens_z",
        "b.object_id AS object_id", "c.i_ra AS source_ra_deg", "c.i_dec AS source_dec_deg",
        "a.i_hsmshaperegauss_e1 AS e1", "a.i_hsmshaperegauss_e2 AS e2",
        "b.i_hsmshaperegauss_derived_sigma_e AS sigma_e",
        "b.i_hsmshaperegauss_derived_rms_e AS rms_e",
        "b.i_hsmshaperegauss_derived_weight AS weight",
        "b.i_hsmshaperegauss_derived_shear_bias_m AS m",
        "b.i_hsmshaperegauss_derived_shear_bias_c1 AS c1",
        "b.i_hsmshaperegauss_derived_shear_bias_c2 AS c2",
        "b.i_apertureflux_10_mag AS mag_a10",
        "b.i_hsmshaperegauss_resolution AS resolution",
        "a.i_sdssshape_psf_shape11 AS psf_shape11",
        "a.i_sdssshape_psf_shape22 AS psf_shape22",
        "a.i_sdssshape_psf_shape12 AS psf_shape12",
        "b.b_mode_mask AS b_mode_mask_raw", "b.hsc_y3_zbin AS source_zbin",
    ]
    return (
        "SELECT DISTINCT\n  " + ",\n  ".join(projection) +
        "\nFROM s19a_wide.meas2 a\n"
        "INNER JOIN s19a_wide.weaklensing_hsm_regauss b USING (object_id)\n"
        "INNER JOIN s19a_wide.meas c USING (object_id)\n"
        f"WHERE coneSearch(c.i_coord, {ra}, {dec}, {radius})"
    )


def load_verified_counts(plan_dir: Path, results_dir: Path):
    count_manifest, centers, batches, plan_sha = validate_count_plan(plan_dir)
    plan_path = Path(plan_dir).resolve() / "random_shape_count_plan.json"
    results_dir = Path(results_dir).resolve()
    checkpoint_path = results_dir / "count_status.json"
    summary_path = results_dir / "count_summary.json"
    if not checkpoint_path.is_file() or not summary_path.is_file():
        raise ValueError("complete approved random shape-count results are required")
    checkpoint, summary = (json.loads(checkpoint_path.read_text(encoding="utf-8")),
                           json.loads(summary_path.read_text(encoding="utf-8")))
    if (checkpoint.get("status") != "complete" or checkpoint.get("plan_sha256") != plan_sha
            or summary.get("status") != "complete" or summary.get("plan_sha256") != plan_sha):
        raise ValueError("random shape-count results do not match the validated plan")
    if any(checkpoint.get(key) is not False for key in (
            "credentials_saved", "source_shape_rows_returned", "photoz_rows_returned",
            "images_returned", "science_measurement")):
        raise ValueError("count-only checkpoint contains data outside its declared scope")
    if any(summary.get(key) is not False for key in (
            "source_shape_rows_returned", "photoz_rows_returned", "images_returned",
            "science_measurement")):
        raise ValueError("count-only summary contains data outside its declared scope")
    if (summary.get("n_random_centers") != len(centers)
            or summary.get("n_matched_lenses") != count_manifest["n_matched_lenses"]):
        raise ValueError("count-only summary does not cover the prepared center set")
    jobs = checkpoint.get("jobs", {})
    size = int(count_manifest["batch_size"])
    all_counts = {}
    batch_records = []
    for index, batch in enumerate(batches):
        name = f"random_shape_counts_{index:03d}"
        item = jobs.get(name, {})
        path = results_dir / f"{name}.csv"
        if item.get("status") != "complete" or not path.is_file() or sha(path) != item.get("sha256"):
            raise ValueError(f"missing, incomplete, or changed count batch: {name}")
        chunk = centers[index * size:(index + 1) * size]
        rows = read_count_rows(path, chunk)
        if item.get("rows") != len(rows):
            raise ValueError(f"saved count-row total changed for {name}")
        all_counts.update(rows)
        batch_records.append({
            "name": name, "sha256": item["sha256"],
            "n_random_centers": len(chunk),
            "n_shape_sources": sum(row["n_shape_sources"] for row in rows.values()),
        })
    expected_total = int(summary["sum_shape_source_counts"])
    observed_total = sum(row["n_shape_sources"] for row in all_counts.values())
    if len(all_counts) != len(centers) or observed_total != expected_total:
        raise ValueError("count-only rows disagree with the completed summary")
    return count_manifest, centers, batches, all_counts, batch_records, plan_sha, summary


def prepare(count_plan_dir: Path, count_results_dir: Path, output: Path) -> dict:
    output = Path(output).resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite row plan: {output}")
    (count_manifest, centers, count_batches, counts, count_records,
     count_plan_sha, count_summary) = load_verified_counts(count_plan_dir, count_results_dir)
    batch_size = int(count_manifest["batch_size"])
    output.mkdir(parents=True, exist_ok=False)
    row_batches = []
    for index, start in enumerate(range(0, len(centers), batch_size)):
        chunk = centers[start:start + batch_size]
        sql = "\nUNION ALL\n".join(shape_rows_query(row) for row in chunk) + "\n"
        path = output / f"random_shape_rows_{index:03d}.sql"
        path.write_text(sql, encoding="utf-8")
        expected_rows = sum(counts[row["random_center_id"]]["n_shape_sources"] for row in chunk)
        row_batches.append({
            "path": str(path), "n_random_centers": len(chunk),
            "first_random_center_id": chunk[0]["random_center_id"],
            "last_random_center_id": chunk[-1]["random_center_id"],
            "expected_pair_rows_from_count_only": expected_rows,
            "sha256": sha(path),
        })
    total_rows = sum(row["n_shape_sources"] for row in counts.values())
    manifest = {
        "schema": "hsc-s19a-random-shape-row-export-plan-v1",
        "status": "prepared_not_submitted",
        "source_count_plan": str((Path(count_plan_dir).resolve() / "random_shape_count_plan.json")),
        "source_count_plan_sha256": count_plan_sha,
        "source_count_summary": str((Path(count_results_dir).resolve() / "count_summary.json")),
        "source_count_summary_sha256": sha(Path(count_results_dir).resolve() / "count_summary.json"),
        "source_count_results_dir": str(Path(count_results_dir).resolve()),
        "source_centers": count_manifest["source_centers"],
        "source_centers_sha256": count_manifest["source_centers_sha256"],
        "release": "pdr3",
        "query_tables": ["s19a_wide.meas2", "s19a_wide.weaklensing_hsm_regauss", "s19a_wide.meas"],
        "n_random_centers": len(centers),
        "n_matched_lenses": count_manifest["n_matched_lenses"],
        "n_query_batches": len(row_batches), "batch_size": batch_size,
        "expected_random_source_pair_rows": total_rows,
        "estimated_uncompressed_csv_bytes": {
            "low_assumption_bytes_per_row": DEFAULT_ROW_BYTES["low"],
            "high_assumption_bytes_per_row": DEFAULT_ROW_BYTES["high"],
            "low": total_rows * DEFAULT_ROW_BYTES["low"],
            "high": total_rows * DEFAULT_ROW_BYTES["high"],
            "low_gib": total_rows * DEFAULT_ROW_BYTES["low"] / (1024 ** 3),
            "high_gib": total_rows * DEFAULT_ROW_BYTES["high"] / (1024 ** 3),
            "basis": "planning range only; measured size must be computed from returned rows",
        },
        "output_columns": OUTPUT_COLUMNS,
        "row_semantics": "one row per (random_center_id, object_id); source IDs may repeat across centers and matched lenses",
        "photoz_pdf_rows_included": False,
        "source_shape_rows_returned": False,
        "images_returned": False,
        "network_used": False,
        "query_submitted": False,
        "science_measurement": False,
        "credentials_saved": False,
        "execution_requires_explicit_user_approval": True,
        "count_batches": count_records,
        "batches": row_batches,
        "notes": [
            "The source query uses the same joins and coneSearch as the count-only preflight.",
            "Each returned batch must match exact per-center distinct-object counts before it is accepted.",
            "The exported rows contain shape-calibration and HSC tomographic-bin fields, not per-object photo-z PDFs.",
            "This plan does not approve random-center validity or a science measurement.",
        ],
    }
    plan_path = output / "random_shape_row_plan.json"
    plan_path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2, allow_nan=False) + "\n",
                         encoding="utf-8")
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--count-plan-dir", type=Path, required=True)
    parser.add_argument("--count-results-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = prepare(args.count_plan_dir, args.count_results_dir, args.output)
    print(json.dumps({key: result[key] for key in (
        "status", "n_random_centers", "n_matched_lenses", "n_query_batches",
        "expected_random_source_pair_rows", "estimated_uncompressed_csv_bytes",
        "output_columns", "query_submitted", "science_measurement")},
        ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
