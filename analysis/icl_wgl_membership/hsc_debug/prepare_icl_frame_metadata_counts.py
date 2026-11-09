#!/usr/bin/env python3
"""Prepare an offline, count-only HSC PDR3 frame-metadata plan for ICL coadds.

Reads existing target/control download manifests and emits SQL only. It does
not access HSC/CAS, request image products, or save credentials. Execution
requires a separate explicit user approval.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
TARGET_MANIFEST = DATA / "download_01/download_manifest.json"
TARGET_STATUS = DATA / "download_01/download_status.json"
CONTROL_MANIFEST = DATA / "control_download_01/download_manifest.json"
CONTROL_STATUS = DATA / "control_download_01/download_status.json"
DEFAULT_OUTPUT = DATA / "icl_frame_metadata_count_plan_2026-11-04_01"
BATCH_SIZE = 100


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def quote(value: str) -> str:
    return "'" + value.replace("'", "''") + "'"


def load_group(manifest_path: Path, status_path: Path, role: str) -> list[dict]:
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    status = json.loads(status_path.read_text(encoding="utf-8"))
    if status.get("status") != "complete" or len(status.get("completed", [])) != status.get("n_expected"):
        raise ValueError(f"{role} image download is not a complete, auditable sample")
    clusters = manifest.get("clusters")
    if not isinstance(clusters, list) or len(clusters) != status["n_expected"]:
        raise ValueError(f"{role} manifest rows do not match the completed download")
    completed_ids = {str(item.get("id")) for item in status["completed"]}
    centers = []
    for row in clusters:
        center_id = str(row.get("id", ""))
        ra, dec = row.get("ra_deg"), row.get("dec_deg")
        band = str(row.get("band", "HSC-I"))
        if (not center_id or center_id not in completed_ids
                or not isinstance(ra, (int, float)) or not isinstance(dec, (int, float))
                or not 0.0 <= float(ra) < 360.0 or not -90.0 <= float(dec) <= 90.0):
            raise ValueError(f"invalid or un-downloaded {role} center row: {center_id!r}")
        if band != "HSC-I":
            raise ValueError(f"unexpected HSC coadd band for {role} row {center_id}: {band}")
        centers.append({
            "center_id": f"{role}:{center_id}",
            "role": role,
            "ra_deg": float(ra),
            "dec_deg": float(dec),
            "band": band,
        })
    if len(completed_ids) != len(centers):
        raise ValueError(f"{role} completed image IDs are not a one-to-one manifest match")
    return centers


def sql_for_batch(centers: list[dict]) -> str:
    values = ",\n        ".join(
        f"({quote(x['center_id'])}, {quote(x['role'])}, {x['ra_deg']:.10f}, "
        f"{x['dec_deg']:.10f}, {quote(x['band'])})"
        for x in centers
    )
    return f"""-- Offline-prepared count-only query; submit only after user approval.
-- One returned row per image center, with counts of candidate coadd-frame membership.
WITH requested(center_id, sample_role, ra_deg, dec_deg, band) AS (
    VALUES
        {values}
), candidates AS (
    SELECT
        r.center_id, r.sample_role, r.ra_deg, r.dec_deg, r.band,
        m.skymap_id, m.tract, m.patch, m.patch_s, m.ra2000, m.dec2000
    FROM requested AS r
    JOIN pdr3_wide.mosaic AS m ON m.band = r.band
    -- PDR3 mosaic exposes areacube directly. The schema browser documents
    -- radec_to_coord() and the areacube overlap operator for spatial searches.
    WHERE m.areacube && radec_to_coord(r.ra_deg, r.dec_deg)
), chosen_patch AS (
    SELECT DISTINCT ON (center_id)
        center_id, ra_deg, dec_deg, band, skymap_id, tract, patch, patch_s, ra2000, dec2000
    FROM candidates
    ORDER BY center_id,
        2.0 * acos(LEAST(1.0, GREATEST(-1.0,
            sin(radians(dec_deg))
            * sin(radians(dec2000))
            + cos(radians(dec_deg))
            * cos(radians(dec2000))
            * cos(radians(ra2000 - ra_deg))
        ))),
        skymap_id
), counts AS (
    SELECT
        p.center_id,
        COUNT(mf.frame_num) AS n_patch_frame_memberships,
        COUNT(DISTINCT mf.frame_num) AS n_unique_frames,
        COUNT(DISTINCT mf.visit) AS n_unique_visits,
        COUNT(DISTINCT mf.ccd) AS n_unique_ccds
    FROM chosen_patch AS p
    LEFT JOIN pdr3_wide.mosaicframe AS mf
      ON mf.skymap_id = p.skymap_id AND mf.band = p.band
    GROUP BY p.center_id
)
SELECT
    r.center_id, r.sample_role, r.ra_deg, r.dec_deg, r.band,
    p.skymap_id, p.tract, p.patch, p.patch_s,
    COALESCE(c.n_patch_frame_memberships, 0) AS n_patch_frame_memberships,
    COALESCE(c.n_unique_frames, 0) AS n_unique_frames,
    COALESCE(c.n_unique_visits, 0) AS n_unique_visits,
    COALESCE(c.n_unique_ccds, 0) AS n_unique_ccds
FROM requested AS r
LEFT JOIN chosen_patch AS p USING (center_id)
LEFT JOIN counts AS c USING (center_id)
ORDER BY r.center_id;
"""


def build(output_dir: Path) -> dict:
    if output_dir.exists():
        raise FileExistsError(f"refusing to overwrite metadata count plan: {output_dir}")
    centers = load_group(TARGET_MANIFEST, TARGET_STATUS, "target")
    centers += load_group(CONTROL_MANIFEST, CONTROL_STATUS, "control")
    centers.sort(key=lambda row: row["center_id"])
    if len(centers) != 1598 or len({row["center_id"] for row in centers}) != len(centers):
        raise ValueError("expected 800 target plus 798 unique matched-control centers")

    output_dir.mkdir(parents=True)
    batches = []
    for start in range(0, len(centers), BATCH_SIZE):
        index = start // BATCH_SIZE
        batch = centers[start:start + BATCH_SIZE]
        sql = sql_for_batch(batch)
        path = output_dir / f"frame_metadata_counts_{index:03d}.sql"
        path.write_text(sql, encoding="utf-8")
        batches.append({
            "name": path.stem,
            "path": str(path),
            "sha256": sha256(path),
            "n_centers": len(batch),
            "first_center_id": batch[0]["center_id"],
            "last_center_id": batch[-1]["center_id"],
            "sql_bytes": path.stat().st_size,
        })

    input_paths = (TARGET_MANIFEST, TARGET_STATUS, CONTROL_MANIFEST, CONTROL_STATUS)
    plan = {
        "schema": "astrafek-icl-frame-metadata-count-plan-v1",
        "status": "prepared_not_submitted",
        "release": "pdr3",
        "tables": ["pdr3_wide.mosaic", "pdr3_wide.mosaicframe"],
        "purpose": "map existing target/control coadd cutout centers to the nearest HSC-I mosaic areacube and count its contributing CCD frames",
        "n_target_images": 800,
        "n_control_images": 798,
        "n_centers": len(centers),
        "n_query_batches": len(batches),
        "batch_size_centers": BATCH_SIZE,
        "rows_returned_per_complete_query": "one count-only row per input image center (1,598 expected)",
        "output_columns": ["center_id", "sample_role", "ra_deg", "dec_deg", "band", "skymap_id", "tract", "patch", "patch_s", "n_patch_frame_memberships", "n_unique_frames", "n_unique_visits", "n_unique_ccds"],
        "source_metadata_columns_not_returned_in_this_count_stage": ["frame.exptime", "frame.gain1", "frame.gain2", "frame.gain3", "frame.gain4", "frame.fluxmag0", "frame.sigma_sky"],
        "images_returned": False,
        "shape_rows_returned": False,
        "photoz_rows_returned": False,
        "network_used": False,
        "query_submitted": False,
        "credentials_saved": False,
        "science_measurement": False,
        "execution_requires_explicit_user_approval": True,
        "input_sha256": {str(path): sha256(path) for path in input_paths},
        "batches": batches,
        "notes": [
            "The official HSC PDR3 mosaicframe table maps each (skymap_id, band) mosaic to contributing frame_num, visit, and CCD rows.",
            "The official pdr3_wide.frame schema exposes exptime and gain1..gain4; this first plan retrieves counts only. A later row-export plan will be sized from observed counts and require separate approval.",
            "HSC-I coadd band normalization is used; original HSC-I versus HSC-I2 per-frame filter will be retained in any later metadata row export.",
            "The nearest-center tie-break among mosaic areacube matches is deterministic; returned patch IDs must be reconciled with the downloaded cutout WCS before using frame gains scientifically.",
        ],
    }
    path = output_dir / "plan.json"
    path.write_text(json.dumps(plan, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": plan["status"],
        "science_measurement": False,
        "centers": len(centers),
        "batches": len(batches),
        "output": str(path.resolve()),
        "sql_sha256": hashlib.sha256("\n".join(b["sha256"] for b in batches).encode()).hexdigest(),
    }, ensure_ascii=False))
    return plan


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    build(args.output_dir.resolve())


if __name__ == "__main__":
    main()
