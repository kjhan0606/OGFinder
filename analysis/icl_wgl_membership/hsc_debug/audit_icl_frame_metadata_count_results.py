#!/usr/bin/env python3
"""Verify the approved PDR3 count-only frame metadata query and summarize depth."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from submit_icl_frame_metadata_counts import (
    load_expected_centers,
    read_result,
    sha256,
    validate_plan,
)


DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
DEFAULT_PLAN = DATA / "icl_frame_metadata_count_plan_2026-11-05_03"
DEFAULT_RESULTS = DATA / "icl_frame_metadata_count_results_2026-11-05_03"
DEFAULT_OUTPUT = DATA / "icl_frame_metadata_count_audit_2026-11-05_01"
CONTROL_MANIFEST = DATA / "control_download_01/download_manifest.json"
COUNT_FIELDS = (
    "n_patch_frame_memberships", "n_unique_frames", "n_unique_visits", "n_unique_ccds",
)


def distribution(values: list[int]) -> dict:
    x = np.asarray(values, dtype=float)
    return {
        "n": int(x.size),
        "min": int(x.min()),
        "p10": float(np.quantile(x, 0.10)),
        "median": float(np.median(x)),
        "mean": float(x.mean()),
        "p90": float(np.quantile(x, 0.90)),
        "max": int(x.max()),
        "sum": int(x.sum()),
    }


def run(plan_dir: Path, results_dir: Path, output: Path) -> dict:
    if output.exists():
        raise FileExistsError(f"refusing to overwrite frame-count audit: {output}")
    plan, centers, batches, plan_hash = validate_plan(plan_dir)
    status_path = results_dir / "status.json"
    summary_path = results_dir / "summary.json"
    status = json.loads(status_path.read_text(encoding="utf-8"))
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    if (plan_hash != "3dac4db631d0dcbee5c5f4509e639f2bd90ca2af99fa7cf5f66c2a9f60bd9014"
            or status.get("plan_sha256") != plan_hash
            or status.get("status") != "complete"
            or status.get("query_submitted") is not True
            or status.get("credentials_saved") is not False
            or status.get("images_returned") is not False
            or status.get("shape_rows_returned") is not False
            or status.get("photoz_rows_returned") is not False
            or status.get("science_measurement") is not False):
        raise ValueError("CAS checkpoint is incomplete, not the approved v3, or outside its scope")
    if (summary.get("status") != "complete" or summary.get("plan_sha256") != plan_hash
            or summary.get("n_centers") != 1598
            or summary.get("n_target_centers") != 800
            or summary.get("n_control_centers") != 798
            or summary.get("images_returned") is not False
            or summary.get("shape_rows_returned") is not False
            or summary.get("photoz_rows_returned") is not False
            or summary.get("credentials_saved") is not False
            or summary.get("science_measurement") is not False):
        raise ValueError("CAS result summary differs from the approved count-only contract")
    expected_centers = {row["center_id"]: row for row in centers}
    rows: dict[str, dict] = {}
    batch_hashes = {}
    jobs = status.get("jobs", {})
    if len(jobs) != 16:
        raise ValueError("expected exactly 16 completed CAS batches")
    for index, batch in enumerate(batches):
        name = batch["name"]
        job = jobs.get(name, {})
        path = results_dir / f"{name}.csv"
        if (job.get("status") != "complete" or not path.is_file()
                or sha256(path) != job.get("sha256")):
            raise ValueError(f"missing or changed count batch: {name}")
        start = index * 100
        expected = centers[start:start + batch["n_centers"]]
        observed = read_result(path, expected)
        if len(observed) != batch["n_centers"] or job.get("rows") != len(observed):
            raise ValueError(f"count-row total mismatch: {name}")
        if set(rows) & set(observed):
            raise ValueError(f"duplicate center across batches: {name}")
        rows.update(observed)
        batch_hashes[name] = job["sha256"]
    if set(rows) != set(expected_centers):
        raise ValueError("verified CAS rows do not cover all requested centers exactly once")
    if summary.get("centers_without_patch") != sum(not row["skymap_id"] for row in rows.values()):
        raise ValueError("summary patch-coverage count does not reproduce")

    by_role = {}
    for role in ("target", "control"):
        selected = [row for row in rows.values() if row["sample_role"] == role]
        by_role[role] = {
            "n_centers": len(selected),
            "centers_without_patch": sum(not row["skymap_id"] for row in selected),
            "count_distributions": {
                field: distribution([int(row[field]) for row in selected])
                for field in COUNT_FIELDS
            },
        }
    if by_role["target"]["n_centers"] != 800 or by_role["control"]["n_centers"] != 798:
        raise ValueError("target/control role counts disagree with the request")
    control_doc = json.loads(CONTROL_MANIFEST.read_text(encoding="utf-8"))
    paired_rows = []
    for control in control_doc["clusters"]:
        target_id = str(control["target_cluster_id"])
        target_center_id = f"target:{target_id}"
        control_center_id = f"control:{control['id']}"
        if target_center_id not in rows or control_center_id not in rows:
            raise ValueError(f"missing target/control count for matched pair {target_id}")
        target_row, control_row = rows[target_center_id], rows[control_center_id]
        if target_row["sample_role"] != "target" or control_row["sample_role"] != "control":
            raise ValueError(f"role mismatch in matched pair {target_id}")
        paired = {"target_cluster_id": target_id, "target_center_id": target_center_id,
                  "control_center_id": control_center_id}
        for field in COUNT_FIELDS:
            paired[f"target_{field}"] = int(target_row[field])
            paired[f"control_{field}"] = int(control_row[field])
            paired[f"delta_{field}"] = int(target_row[field]) - int(control_row[field])
        paired_rows.append(paired)
    if len(paired_rows) != 798 or len({row["target_cluster_id"] for row in paired_rows}) != 798:
        raise ValueError("expected 798 unique target-control pairs")
    paired_differences = {}
    for field in COUNT_FIELDS:
        delta = [row[f"delta_{field}"] for row in paired_rows]
        paired_differences[field] = {
            **distribution(delta),
            "fraction_equal": float(np.mean(np.asarray(delta) == 0)),
            "fraction_target_higher": float(np.mean(np.asarray(delta) > 0)),
            "fraction_control_higher": float(np.mean(np.asarray(delta) < 0)),
        }
    totals = {field: sum(int(row[field]) for row in rows.values()) for field in COUNT_FIELDS}
    summary_fields = {
        "n_patch_frame_memberships": "sum_patch_frame_memberships",
        "n_unique_frames": "sum_unique_frames_per_center",
        "n_unique_visits": "sum_unique_visits_per_center",
    }
    for field, summary_field in summary_fields.items():
        if int(summary[summary_field]) != totals[field]:
            raise ValueError(f"total {field} does not match CAS summary")
    output.mkdir(parents=True, exist_ok=False)
    tsv = output / "frame_counts_verified.tsv"
    fields = ["center_id", "sample_role", "ra_deg", "dec_deg", "band", "skymap_id",
              "tract", "patch", "patch_s", *COUNT_FIELDS]
    with tsv.open("w", encoding="utf-8", newline="") as stream:
        import csv
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(sorted(rows.values(), key=lambda row: row["center_id"]))
    paired_tsv = output / "target_control_frame_count_differences.tsv"
    paired_fields = list(paired_rows[0])
    with paired_tsv.open("w", encoding="utf-8", newline="") as stream:
        import csv
        writer = csv.DictWriter(stream, fieldnames=paired_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(sorted(paired_rows, key=lambda row: int(row["target_cluster_id"])))
    result = {
        "schema": "astrafek-icl-frame-metadata-count-audit-v1",
        "status": "complete_exploratory",
        "science_measurement": False,
        "network_used_by_auditor": False,
        "query_was_approved_and_submitted": True,
        "query_plan_sha256": plan_hash,
        "query_submitted": True,
        "query_result_scope": "1,598 target/control centers; HSC-I patch and integer frame/visit/CCD counts only",
        "n_centers": len(rows),
        "n_target_centers": by_role["target"]["n_centers"],
        "n_control_centers": by_role["control"]["n_centers"],
        "centers_without_patch": sum(not row["skymap_id"] for row in rows.values()),
        "by_role": by_role,
        "n_paired_target_control_centers": len(paired_rows),
        "paired_target_minus_control_differences": paired_differences,
        "total_counts": totals,
        "not_returned": {
            "images": True, "shape_rows": True, "photoz_rows": True,
            "gain_exposure_sky_values": True, "credentials": True,
        },
        "input_hashes": {
            "plan_json_sha256": plan_hash,
            "cas_status_sha256": sha256(status_path),
            "cas_summary_sha256": sha256(summary_path),
            "control_manifest_sha256": sha256(CONTROL_MANIFEST),
            "batch_hashes_sha256": hashlib.sha256("\n".join(
                f"{name}:{digest}" for name, digest in sorted(batch_hashes.items())
            ).encode()).hexdigest(),
        },
        "verified_rows_tsv": str(tsv),
        "verified_rows_tsv_sha256": sha256(tsv),
        "paired_differences_tsv": str(paired_tsv),
        "paired_differences_tsv_sha256": sha256(paired_tsv),
        "scientific_use_limit": "Frame-count differences can diagnose heterogeneous coadd coverage, but counts alone do not provide pixel-level weights or Poisson variance.",
    }
    report = output / "summary.json"
    report.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan-dir", type=Path, default=DEFAULT_PLAN)
    parser.add_argument("--results", type=Path, default=DEFAULT_RESULTS)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    run(args.plan_dir, args.results, args.output)
