#!/usr/bin/env python3
"""Verify saved random-control coordinate batches without contacting HSC CAS."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from submit_wl_random_control_coordinates import (
    DEFAULT_PLAN,
    read_coordinates,
    validate_coordinate_plan,
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def audit(plan_dir: Path) -> dict:
    """Check checkpoint and completed files against the approved SQL plan."""
    manifest, batches, _lenses, lens_ids, expected_counts = validate_coordinate_plan(plan_dir)
    results_dir = plan_dir / "cas_coordinates"
    checkpoint_path = results_dir / "coordinate_status.json"
    if not checkpoint_path.is_file():
        raise FileNotFoundError(f"coordinate checkpoint is missing: {checkpoint_path}")
    checkpoint = json.loads(checkpoint_path.read_text(encoding="utf-8"))
    if checkpoint.get("schema") != "hsc-pdr3-random-control-coordinate-run-v1":
        raise ValueError("unexpected coordinate checkpoint schema")
    plan_path = plan_dir / "random_control_coordinate_plan.json"
    plan_hash = sha256(plan_path)
    if checkpoint.get("plan_sha256") != plan_hash:
        raise ValueError("coordinate checkpoint belongs to a different plan")
    if checkpoint.get("status") not in {
        "pending", "running", "complete", "interrupted", "stopped_on_error"
    }:
        raise ValueError("unknown coordinate checkpoint status")
    if any(checkpoint.get(key) is not False for key in (
        "credentials_saved", "source_rows_returned", "shape_rows_returned",
        "photoz_rows_returned", "images_returned", "science_measurement"
    )):
        raise ValueError("checkpoint reports data outside the approved coordinate-only scope")

    expected_names = {Path(batch["path"]).stem for batch in batches}
    jobs = checkpoint.get("jobs", {})
    if not isinstance(jobs, dict) or set(jobs) - expected_names:
        raise ValueError("checkpoint contains an unplanned batch")
    records = []
    total_rows = 0
    for index, batch in enumerate(batches):
        name = Path(batch["path"]).stem
        item = jobs.get(name, {"status": "pending"})
        status = item.get("status")
        if status not in {"pending", "submitting", "submitted", "complete", "error"}:
            raise ValueError(f"unknown batch status for {name}: {status}")
        csv_path = results_dir / f"{name}.csv"
        part_path = results_dir / f"{name}.csv.part"
        if status == "complete":
            if not csv_path.is_file():
                raise FileNotFoundError(f"completed coordinate batch is missing: {csv_path}")
            digest = sha256(csv_path)
            if digest != item.get("sha256"):
                raise ValueError(f"coordinate batch checksum mismatch: {name}")
            start = index * int(manifest["batch_size"])
            ids = lens_ids[start:start + int(manifest["batch_size"])]
            counts = {lens_id: expected_counts[lens_id] for lens_id in ids}
            rows = read_coordinates(csv_path, counts)
            if rows != item.get("rows") or rows != batch.get("expected_rows"):
                raise ValueError(f"coordinate batch row count differs from approved plan: {name}")
            total_rows += rows
        else:
            digest = None
            rows = None
            if csv_path.exists():
                raise ValueError(f"uncheckpointed coordinate file must not be counted: {csv_path}")
        records.append({
            "batch": name,
            "status": status,
            "job_id": str(item["job_id"]) if item.get("job_id") is not None else None,
            "rows_verified": rows,
            "bytes": csv_path.stat().st_size if status == "complete" else None,
            "sha256": digest,
            "partial_download_present": part_path.exists(),
        })

    if checkpoint["status"] == "complete" and any(row["status"] != "complete" for row in records):
        raise ValueError("checkpoint claims completion while planned batches remain incomplete")
    completed = sum(row["status"] == "complete" for row in records)
    return {
        "schema": "astrafek-wl-random-coordinate-checkpoint-audit-v1",
        "status": "verified_complete" if completed == len(records) else "verified_partial",
        "checkpoint_status": checkpoint["status"],
        "checkpoint_error_type": checkpoint.get("error_type"),
        "n_lenses": len(lens_ids),
        "n_batches": len(records),
        "n_batches_verified": completed,
        "n_batches_remaining": len(records) - completed,
        "expected_coordinate_rows": int(manifest["candidate_pair_rows_expected"]),
        "verified_coordinate_rows": total_rows,
        "results_dir": str(results_dir.resolve()),
        "plan_sha256": plan_hash,
        "checkpoint_sha256": sha256(checkpoint_path),
        "credentials_saved": False,
        "source_rows_returned": False,
        "shape_rows_returned": False,
        "photoz_rows_returned": False,
        "images_returned": False,
        "science_measurement": False,
        "process_liveness_verified": False,
        "batches": records,
    }


def write_report(output_dir: Path, result: dict, script_path: Path) -> None:
    if output_dir.exists():
        raise FileExistsError(f"refusing to overwrite audit output: {output_dir}")
    result["audit_script_sha256"] = sha256(script_path)
    output_dir.mkdir(parents=True, exist_ok=False)
    json_path = output_dir / "checkpoint_audit.json"
    json_path.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    complete_rows = [row for row in result["batches"] if row["status"] == "complete"]
    lines = [
        "# HSC random-control coordinate checkpoint audit",
        "",
        f"- Result: `{result['status']}`; checkpoint: `{result['checkpoint_status']}`",
        f"- Verified batches: {result['n_batches_verified']}/{result['n_batches']}",
        f"- Verified coordinate rows: {result['verified_coordinate_rows']}/"
        f"{result['expected_coordinate_rows']}",
        f"- Saved complete-batch files checked: {len(complete_rows)}",
        f"- Checkpoint SHA-256: `{result['checkpoint_sha256']}`",
        f"- Plan SHA-256: `{result['plan_sha256']}`",
        "- This audit makes no HSC request and verifies no live process state.",
        "- No source rows, shapes, photo-z values, or images were returned; "
        "`science_measurement=false`.",
        "",
        "## Batch state",
        "",
        "| Batch | State | Job | Verified rows | SHA-256 |",
        "|---|---|---:|---:|---|",
    ]
    for row in result["batches"]:
        lines.append(
            f"| {row['batch']} | {row['status']} | {row['job_id'] or '—'} | "
            f"{row['rows_verified'] if row['rows_verified'] is not None else '—'} | "
            f"{row['sha256'] or '—'} |"
        )
    (output_dir / "assessment.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan-dir", type=Path, default=DEFAULT_PLAN)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    result = audit(args.plan_dir.resolve())
    write_report(args.output_dir.resolve(), result, Path(__file__).resolve())
    print(json.dumps({key: result[key] for key in (
        "status", "checkpoint_status", "n_batches_verified", "n_batches_remaining",
        "expected_coordinate_rows", "verified_coordinate_rows", "science_measurement",
    )}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
