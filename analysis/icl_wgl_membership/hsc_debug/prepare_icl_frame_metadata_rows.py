#!/usr/bin/env python3
"""Prepare an offline HSC PDR3 frame-metadata row export for ICL noise work.

The plan is derived from the completed and audited count-only query. It emits
SQL only; submission requires a separate, exact user approval.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

from submit_icl_frame_metadata_counts import BATCH_SIZE, validate_plan


DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
COUNT_PLAN = DATA / "icl_frame_metadata_count_plan_2026-11-05_03"
COUNT_RESULTS = DATA / "icl_frame_metadata_count_results_2026-11-05_03"
COUNT_AUDIT = DATA / "icl_frame_metadata_count_audit_2026-11-05_02"
DEFAULT_OUTPUT = DATA / "icl_frame_metadata_row_plan_2026-11-05_01"
EXPECTED_COUNT_PLAN_SHA256 = "3dac4db631d0dcbee5c5f4509e639f2bd90ca2af99fa7cf5f66c2a9f60bd9014"
EXPECTED_COUNT_AUDIT_SHA256 = "63e7424cb6abe3f89b7138bb99ea1072fbb367ab1e18b014974981334522810b"

OUTPUT_COLUMNS = [
    "center_id", "sample_role", "ra_deg", "dec_deg", "coadd_band",
    "skymap_id", "tract", "patch", "patch_s", "frame_num", "visit", "ccd",
    "filter", "exptime", "gain1", "gain2", "gain3", "gain4", "fluxmag0",
    "sigma_sky", "mjd",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def sql_quote(value: str) -> str:
    return "'" + value.replace("'", "''") + "'"


def sql_for_batch(centers: list[dict]) -> str:
    values = ",\n        ".join(
        f"({sql_quote(row['center_id'])}, {sql_quote(row['role'])}, "
        f"{float(row['ra_deg']):.10f}, {float(row['dec_deg']):.10f}, "
        f"{sql_quote(row['band'])})"
        for row in centers
    )
    columns = ",\n    ".join(f"f.{name}" for name in (
        "filter", "exptime", "gain1", "gain2", "gain3", "gain4",
        "fluxmag0", "sigma_sky", "mjd",
    ))
    return f"""-- Approved-scope candidate: per-center contributing frame metadata only.
-- Do not submit without explicit approval of this exact plan SHA-256.
WITH requested(center_id, sample_role, ra_deg, dec_deg, band) AS (
    VALUES
        {values}
), candidates AS (
    SELECT r.center_id, r.sample_role, r.ra_deg, r.dec_deg, r.band,
           m.skymap_id, m.tract, m.patch, m.patch_s, m.ra2000, m.dec2000
    FROM requested AS r
    JOIN pdr3_wide.mosaic AS m ON m.band = r.band
    WHERE m.areacube && radec_to_coord(r.ra_deg, r.dec_deg)
), chosen_patch AS (
    SELECT DISTINCT ON (center_id)
           center_id, ra_deg, dec_deg, band, skymap_id, tract, patch,
           patch_s, ra2000, dec2000
    FROM candidates
    ORDER BY center_id,
        2.0 * acos(LEAST(1.0, GREATEST(-1.0,
            sin(radians(dec_deg)) * sin(radians(dec2000))
            + cos(radians(dec_deg)) * cos(radians(dec2000))
              * cos(radians(ra2000 - ra_deg))
        ))), skymap_id
)
SELECT r.center_id, r.sample_role, r.ra_deg, r.dec_deg,
       r.band AS coadd_band, p.skymap_id, p.tract, p.patch, p.patch_s,
       mf.frame_num, mf.visit, mf.ccd,
       {columns}
FROM requested AS r
JOIN chosen_patch AS p USING (center_id)
JOIN pdr3_wide.mosaicframe AS mf
  ON mf.skymap_id = p.skymap_id AND mf.band = p.band
JOIN pdr3_wide.frame AS f ON f.frame_num = mf.frame_num
ORDER BY r.center_id, mf.frame_num;
"""


def load_audited_centers() -> tuple[list[dict], dict, dict]:
    plan, centers, batches, plan_hash = validate_plan(COUNT_PLAN)
    if plan_hash != EXPECTED_COUNT_PLAN_SHA256:
        raise ValueError("count-only parent plan hash changed")
    status = json.loads((COUNT_RESULTS / "status.json").read_text(encoding="utf-8"))
    cas_summary = json.loads((COUNT_RESULTS / "summary.json").read_text(encoding="utf-8"))
    if (status.get("status") != "complete" or status.get("query_submitted") is not True
            or status.get("plan_sha256") != plan_hash or len(status.get("jobs", {})) != 16):
        raise ValueError("approved count-only parent query is not complete")
    for batch in batches:
        job = status["jobs"].get(batch["name"], {})
        path = COUNT_RESULTS / f"{batch['name']}.csv"
        if (job.get("status") != "complete" or not path.is_file()
                or sha256(path) != job.get("sha256")):
            raise ValueError(f"count-only parent batch is missing or changed: {batch['name']}")

    audit_path = COUNT_AUDIT / "summary.json"
    audit = json.loads(audit_path.read_text(encoding="utf-8"))
    if (sha256(audit_path) != EXPECTED_COUNT_AUDIT_SHA256
            or audit.get("status") != "complete_exploratory"
            or audit.get("query_plan_sha256") != plan_hash
            or audit.get("n_centers") != 1598
            or audit.get("n_target_centers") != 800
            or audit.get("n_control_centers") != 798
            or audit.get("centers_without_patch") != 0
            or audit.get("science_measurement") is not False):
        raise ValueError("count-only audit is incomplete or differs from the approved parent")
    rows_path = Path(audit["verified_rows_tsv"])
    if sha256(rows_path) != audit.get("verified_rows_tsv_sha256"):
        raise ValueError("audited center table hash does not match its summary")
    counts: dict[str, dict] = {}
    with rows_path.open(encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        expected_header = {
            "center_id", "sample_role", "ra_deg", "dec_deg", "band", "skymap_id",
            "tract", "patch", "patch_s", "n_patch_frame_memberships",
            "n_unique_frames", "n_unique_visits", "n_unique_ccds",
        }
        if set(reader.fieldnames or []) != expected_header:
            raise ValueError("audited frame-count table columns changed")
        for row in reader:
            if row["center_id"] in counts or not row["skymap_id"]:
                raise ValueError("duplicate center or missing mapped patch in audited rows")
            counts[row["center_id"]] = row
    if len(counts) != 1598 or set(counts) != {row["center_id"] for row in centers}:
        raise ValueError("audited center table does not match all requested count-plan centers")

    ordered = []
    for center in centers:
        row = counts[center["center_id"]]
        if (row["sample_role"] != center["role"] or row["band"] != center["band"]
                or abs(float(row["ra_deg"]) - center["ra_deg"]) > 1e-8
                or abs(float(row["dec_deg"]) - center["dec_deg"]) > 1e-8):
            raise ValueError(f"audited center identity changed: {center['center_id']}")
        for field in ("n_patch_frame_memberships", "n_unique_frames", "n_unique_visits", "n_unique_ccds"):
            row[field] = int(row[field])
        if row["n_patch_frame_memberships"] != row["n_unique_frames"]:
            raise ValueError("parent counts show duplicated frame membership for a center")
        ordered.append({**center, **row})
    expected_total = sum(row["n_patch_frame_memberships"] for row in ordered)
    if expected_total != 67071 or expected_total != audit["total_counts"]["n_patch_frame_memberships"]:
        raise ValueError("per-center counts do not reproduce the 67,071-row metadata estimate")
    return ordered, audit, cas_summary


def build(output_dir: Path) -> dict:
    output_dir = Path(output_dir).resolve()
    if output_dir.exists():
        raise FileExistsError(f"refusing to overwrite metadata row plan: {output_dir}")
    centers, audit, cas_summary = load_audited_centers()
    output_dir.mkdir(parents=True)
    batches = []
    expected_rows_total = 0
    for start in range(0, len(centers), BATCH_SIZE):
        index = start // BATCH_SIZE
        batch = centers[start:start + BATCH_SIZE]
        sql = sql_for_batch(batch)
        path = output_dir / f"frame_metadata_rows_{index:03d}.sql"
        path.write_text(sql, encoding="utf-8")
        expected_rows = sum(row["n_patch_frame_memberships"] for row in batch)
        expected_rows_total += expected_rows
        batches.append({
            "name": path.stem,
            "path": str(path),
            "sha256": sha256(path),
            "n_centers": len(batch),
            "expected_rows_from_approved_count_query": expected_rows,
            "first_center_id": batch[0]["center_id"],
            "last_center_id": batch[-1]["center_id"],
            "sql_bytes": path.stat().st_size,
        })
    if len(batches) != 16 or expected_rows_total != 67071:
        raise ValueError("unexpected metadata batch or row estimate")
    sql_aggregate = hashlib.sha256("\n".join(item["sha256"] for item in batches).encode()).hexdigest()
    input_paths = [
        COUNT_PLAN / "plan.json", COUNT_RESULTS / "status.json",
        COUNT_RESULTS / "summary.json", COUNT_AUDIT / "summary.json",
        Path(audit["verified_rows_tsv"]),
    ]
    plan = {
        "schema": "astrafek-icl-frame-metadata-row-plan-v1",
        "status": "prepared_not_submitted",
        "parent_count_plan_sha256": EXPECTED_COUNT_PLAN_SHA256,
        "parent_count_audit_sha256": EXPECTED_COUNT_AUDIT_SHA256,
        "release": "pdr3",
        "tables": ["pdr3_wide.mosaic", "pdr3_wide.mosaicframe", "pdr3_wide.frame"],
        "purpose": "retrieve contributing CCD-frame calibration metadata for existing HSC-I target/control coadd centers, to assess whether synthetic-source Poisson variance can be propagated",
        "n_target_centers": 800,
        "n_control_centers": 798,
        "n_centers": len(centers),
        "n_query_batches": len(batches),
        "batch_size_centers": BATCH_SIZE,
        "estimated_rows": expected_rows_total,
        "row_estimate_source": "exact sum of per-center mosaicframe memberships from the approved and audited count-only result",
        "output_columns": OUTPUT_COLUMNS,
        "columns_not_returned": ["image pixels", "shape measurements", "photo-z", "object catalogs", "per-pixel coadd weights", "credentials"],
        "images_returned": False,
        "shape_rows_returned": False,
        "photoz_rows_returned": False,
        "network_used": False,
        "query_submitted": False,
        "credentials_saved": False,
        "science_measurement": False,
        "execution_requires_explicit_user_approval": True,
        "count_parent_cas_summary_sha256": sha256(COUNT_RESULTS / "summary.json"),
        "input_sha256": {str(path): sha256(path) for path in input_paths},
        "batches": batches,
        "sql_aggregate_sha256": sql_aggregate,
        "estimated_csv_size": "approximately 13–27 MB, using a conservative 200–400 bytes per 21-column row; actual size depends on returned numeric formatting",
        "scientific_limit": "Frame-level exposure, amplifier gains, zeropoint, sky RMS, and MJD do not provide pixel-specific coadd weights or map each pixel to an amplifier. Do not use this export alone to set shot-noise variance or promote science_measurement.",
        "notes": [
            "HSC PDR3 mosaicframe explicitly maps (skymap_id, band) to frame_num, visit, and CCD; frame_num is the primary key in pdr3_wide.frame.",
            "Original frame.filter is returned because HSC-I2 is merged into HSC-I in the coadd schema.",
            "The output is center-frame membership metadata, not object-level catalog data.",
            "A separate approval is required before any CAS submission; the previous count-only approval does not authorize these rows.",
        ],
    }
    plan_path = output_dir / "plan.json"
    plan_path.write_text(json.dumps(plan, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    result = {
        "status": plan["status"], "plan": str(plan_path),
        "plan_sha256": sha256(plan_path), "sql_aggregate_sha256": sql_aggregate,
        "centers": len(centers), "batches": len(batches),
        "estimated_rows": expected_rows_total,
        "columns": len(OUTPUT_COLUMNS), "sql_bytes": sum(x["sql_bytes"] for x in batches),
        "network_used": False, "query_submitted": False,
        "images": False, "shape_rows": False, "photoz": False,
        "science_measurement": False,
    }
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return plan


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    build(args.output_dir)


if __name__ == "__main__":
    main()
