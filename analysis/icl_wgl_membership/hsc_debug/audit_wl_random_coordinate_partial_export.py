#!/usr/bin/env python3
"""Audit saved batches from the approved random-coordinate export offline.

This tool never contacts HSC CAS and never submits a query or downloads data.
It verifies only batches already marked complete in the local checkpoint.
Submitted jobs are reported as unresolved and are not retried or fetched.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from submit_wl_random_control_coordinates import (
    DEFAULT_PLAN,
    read_coordinates,
    sha,
    validate_coordinate_plan,
)


def audit_checkpoint(plan_manifest: dict, batches: list[dict], lens_ids: list[str],
                     counts: dict[str, int], state: dict, results_dir: Path,
                     plan_sha256: str) -> dict:
    if state.get("schema") != "hsc-pdr3-random-control-coordinate-run-v1":
        raise ValueError("unexpected coordinate-export checkpoint schema")
    if state.get("plan_sha256") != plan_sha256:
        raise ValueError("checkpoint belongs to a different coordinate plan")
    restricted_false = (
        "credentials_saved", "source_rows_returned", "shape_rows_returned",
        "photoz_rows_returned", "images_returned", "science_measurement",
    )
    if any(state.get(key) is not False for key in restricted_false):
        raise ValueError("checkpoint contains state outside the coordinate-only scope")

    jobs = state.get("jobs")
    if not isinstance(jobs, dict):
        raise ValueError("checkpoint jobs must be an object")
    expected_names = {Path(batch["path"]).stem for batch in batches}
    if not set(jobs).issubset(expected_names):
        raise ValueError("checkpoint contains a job outside the approved plan")

    verified, unresolved = [], []
    for index, batch in enumerate(batches):
        name = Path(batch["path"]).stem
        item = jobs.get(name, {"status": "pending"})
        status = item.get("status")
        result_path = results_dir / f"{name}.csv"
        if status == "complete":
            if not result_path.is_file():
                raise ValueError(f"completed coordinate batch is missing: {name}")
            actual_sha = sha(result_path)
            if actual_sha != item.get("sha256"):
                raise ValueError(f"completed coordinate batch checksum mismatch: {name}")
            start = index * int(plan_manifest["batch_size"])
            batch_lens = lens_ids[start:start + int(plan_manifest["batch_size"])]
            expected_counts = {lens_id: counts[lens_id] for lens_id in batch_lens}
            rows = read_coordinates(result_path, expected_counts)
            if rows != int(batch["expected_rows"]) or rows != item.get("rows"):
                raise ValueError(f"completed coordinate batch row count mismatch: {name}")
            verified.append({
                "batch": name, "job_id": str(item.get("job_id", "")),
                "rows": rows, "bytes": result_path.stat().st_size,
                "sha256": actual_sha,
            })
        else:
            if status not in (None, "pending", "submitting", "submitted", "error"):
                raise ValueError(f"unrecognized checkpoint state for {name}: {status}")
            unresolved.append({
                "batch": name, "status": status or "pending",
                "job_id": str(item["job_id"]) if item.get("job_id") is not None else None,
                "saved_csv_present": result_path.exists(),
            })

    if state.get("status") == "complete" and unresolved:
        raise ValueError("checkpoint claims completion but has unresolved batches")
    status = "complete_verified" if not unresolved else "partial_verified"
    return {
        "schema": "astrafek-wl-random-coordinate-partial-export-audit-v1",
        "status": status,
        "plan_sha256": plan_sha256,
        "checkpoint_status": state.get("status", "unknown"),
        "n_plan_lenses": len(lens_ids),
        "n_plan_batches": len(batches),
        "expected_coordinate_rows": int(plan_manifest["candidate_pair_rows_expected"]),
        "n_verified_batches": len(verified),
        "verified_coordinate_rows": sum(row["rows"] for row in verified),
        "verified_batches": verified,
        "unresolved_batches": unresolved,
        "network_used": False,
        "query_submitted": False,
        "rows_downloaded": False,
        "credentials_saved": False,
        "source_rows_returned": False,
        "shape_rows_returned": False,
        "photoz_rows_returned": False,
        "images_returned": False,
        "science_measurement": False,
    }


def audit(plan_dir: Path, results_dir: Path) -> dict:
    plan_dir, results_dir = Path(plan_dir).resolve(), Path(results_dir).resolve()
    plan_manifest, batches, _, lens_ids, counts = validate_coordinate_plan(plan_dir)
    checkpoint_path = results_dir / "coordinate_status.json"
    if not checkpoint_path.is_file():
        raise FileNotFoundError(f"coordinate checkpoint not found: {checkpoint_path}")
    state = json.loads(checkpoint_path.read_text(encoding="utf-8"))
    plan_path = plan_dir / "random_control_coordinate_plan.json"
    return audit_checkpoint(plan_manifest, batches, lens_ids, counts, state,
                            results_dir, sha(plan_path))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", type=Path, default=DEFAULT_PLAN)
    parser.add_argument("--results", type=Path,
                        help="defaults to <plan>/cas_coordinates")
    parser.add_argument("--output", type=Path,
                        help="optional new JSON file; refuses to overwrite")
    args = parser.parse_args()
    results = args.results or (args.plan / "cas_coordinates")
    result = audit(args.plan, results)
    rendered = json.dumps(result, ensure_ascii=False, indent=2, sort_keys=True) + "\n"
    if args.output:
        if args.output.exists():
            raise FileExistsError(f"refusing to overwrite {args.output}")
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    print(rendered, end="")


if __name__ == "__main__":
    main()
