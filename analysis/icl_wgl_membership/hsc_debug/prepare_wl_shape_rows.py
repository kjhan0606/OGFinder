#!/usr/bin/env python3
"""Prepare (but never submit) minimal HSC S19A lens-source row queries.

The plan is tied to the verified count-only result. It requests only the
explicit columns needed for calibrated shape stacking and retains one row per
(lens_id, object_id) pair. No network calls or source-row retrieval occur here.
"""
import argparse
import json
import math
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from submit_wl_shape_counts import load_plan, read_count_csv

DEFAULT_PREFLIGHT = Path(
    "/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/"
    "wl_shape_preflight_2026-10-18_01"
)
DEFAULT_OUTPUT = Path(
    "/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/"
    "wl_shape_rows_plan_2026-10-18_04"
)

OUTPUT_COLUMNS = [
    "lens_id", "lens_ra_deg", "lens_dec_deg", "lens_z", "object_id",
    "source_ra_deg", "source_dec_deg", "e1", "e2", "sigma_e", "rms_e",
    "weight", "m", "c1", "c2", "mag_a10", "resolution",
    "psf_shape11", "psf_shape22", "psf_shape12", "b_mode_mask_raw",
    "source_zbin",
]


def sql_literal(value):
    return "'" + str(value).replace("'", "''") + "'"


def number(value):
    value = float(value)
    if not math.isfinite(value):
        raise ValueError("lens coordinates and redshifts must be finite")
    return format(value, ".12g")


def lens_query(lens):
    lens_id = sql_literal(lens["lens_id"])
    ra = number(lens["ra_deg"])
    dec = number(lens["dec_deg"])
    redshift = number(lens["z_lens"])
    radius = number(lens["radius_arcsec"])
    projection = [
        f"{lens_id} AS lens_id", f"{ra} AS lens_ra_deg",
        f"{dec} AS lens_dec_deg", f"{redshift} AS lens_z",
        "b.object_id AS object_id", "c.i_ra AS source_ra_deg",
        "c.i_dec AS source_dec_deg",
        "a.i_hsmshaperegauss_e1 AS e1",
        "a.i_hsmshaperegauss_e2 AS e2",
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


def load_counts(preflight, query_plan):
    result_dir = preflight / "cas_counts"
    status_path = result_dir / "count_status.json"
    summary_path = result_dir / "count_summary.json"
    if not status_path.is_file() or not summary_path.is_file():
        raise ValueError("verified count-only results are required before planning row export")
    state = json.loads(status_path.read_text(encoding="utf-8"))
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    if state.get("status") != "complete" or summary.get("status") != "complete":
        raise ValueError("count-only result is not complete")
    if state.get("source_rows_downloaded") is not False or summary.get("source_rows_downloaded") is not False:
        raise ValueError("count-only checkpoint does not certify that source rows were absent")
    counts = {}
    for index in range(len(query_plan["sql_files"])):
        path = result_dir / f"shape_counts_{index:03d}.csv"
        start = index * int(query_plan["batch_size"])
        lens_ids = [str(row["lens_id"]) for row in query_plan["lenses"][start:start + int(query_plan["batch_size"])] ]
        rows = read_count_csv(path, lens_ids)
        for row in rows:
            counts[row["lens_id"]] = row["n_shape_sources"]
    if len(counts) != 800 or len(set(counts)) != 800:
        raise ValueError("count table must cover all 800 unique lenses")
    if sum(counts.values()) != int(summary["sum_lens_cone_shape_counts"]):
        raise ValueError("count row sum disagrees with the verified summary")
    return counts, state, summary


def prepare(preflight, output):
    preflight = preflight.resolve()
    output = output.resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite existing plan directory: {output}")
    query_plan = load_plan(preflight)
    counts, state, summary = load_counts(preflight, query_plan)
    batch_size = int(query_plan["batch_size"])
    output.mkdir(parents=True, exist_ok=False)
    batches = []
    lenses = query_plan["lenses"]
    for index, start in enumerate(range(0, len(lenses), batch_size)):
        batch = lenses[start:start + batch_size]
        sql = "\nUNION ALL\n".join(lens_query(lens) for lens in batch) + "\n"
        path = output / f"shape_rows_{index:03d}.sql"
        path.write_text(sql, encoding="utf-8")
        expected_rows = sum(counts[str(lens["lens_id"])] for lens in batch)
        batches.append({
            "path": str(path), "n_lenses": len(batch),
            "first_lens_id": str(batch[0]["lens_id"]),
            "last_lens_id": str(batch[-1]["lens_id"]),
            "expected_pair_rows_from_count_only": expected_rows,
            "sha256": __import__("hashlib").sha256(path.read_bytes()).hexdigest(),
        })
    total = int(summary["sum_lens_cone_shape_counts"])
    manifest = {
        "schema": "hsc-s19a-shape-row-export-plan-v1",
        "status": "prepared_not_submitted",
        "source_count_summary": str((preflight / "cas_counts" / "count_summary.json")),
        "source_count_plan_sha256": state["plan_sha256"],
        "release": "pdr3", "table_join": [
            "s19a_wide.meas2", "s19a_wide.weaklensing_hsm_regauss", "s19a_wide.meas",
        ],
        "n_lenses": 800, "n_query_batches": len(batches),
        "lens_source_pair_rows_expected": total,
        "unique_source_rows_expected": None,
        "output_columns": OUTPUT_COLUMNS,
        "row_semantics": "one row per lens-source pair; the same object_id may occur for multiple lenses",
        "redshift_input": "HSC hsc_y3_zbin (0 means no tomographic bin); calibrated n(z) must be supplied separately",
        "photoz_pdf_rows_included": False,
        "estimated_uncompressed_csv_bytes": {
            "low_assumption_bytes_per_row": 256,
            "high_assumption_bytes_per_row": 480,
            "low": total * 256, "high": total * 480,
            "low_gib": total * 256 / (1024**3),
            "high_gib": total * 480 / (1024**3),
            "basis": "planning range only; 22 explicit CSV columns, not measured from returned rows",
        },
        "batches": batches,
        "source_rows_downloaded": False,
        "science_measurement": False,
        "credentials_saved": False,
        "execution_requires_explicit_user_approval": True,
        "notes": [
            "SELECT DISTINCT and the same joins/coneSearch as the count query should return exactly the prior per-lens counts.",
            "Each returned batch must be checked against exact per-lens counts and unique (lens_id, object_id) pairs before it is accepted.",
            "The total of 3,069,544 is a lens-source pair upper bound for unique source transfer; overlapping source IDs are intentionally retained as lens pairs here.",
            "Full DNNz P(z) archives are not part of this plan; HSC Y3 source n(z) calibration is a separate dependency and download approval.",
        ],
    }
    path = output / "row_export_plan.json"
    path.write_text(json.dumps(manifest, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--preflight", type=Path, default=DEFAULT_PREFLIGHT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    manifest = prepare(args.preflight, args.output)
    estimate = manifest["estimated_uncompressed_csv_bytes"]
    print(json.dumps({
        "status": manifest["status"], "network_used": False,
        "lenses": manifest["n_lenses"], "batches": manifest["n_query_batches"],
        "lens_source_pair_rows_expected": manifest["lens_source_pair_rows_expected"],
        "estimated_uncompressed_gib": [round(estimate["low_gib"], 2), round(estimate["high_gib"], 2)],
        "output_columns": len(manifest["output_columns"]),
        "shape_rows_retrieved": False, "output": str(args.output.resolve()),
    }, ensure_ascii=False))


if __name__ == "__main__":
    main()
