#!/usr/bin/env python3
"""Reconcile completed HSC random-shape exports with saved stack provenance.

This local integrity audit checks immutable plan hashes, all 124 CSV hashes and
byte counts, completed job/row totals, and the random-row provenance recorded by
the two existing stack runs. It makes no HSC/CAS requests and reads no
credentials. It does not inspect or recompute the scientific estimator.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter
from pathlib import Path


DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
PLAN_DIR = DATA / "wl_random_shape_rows_lens_reuse_2026-11-03_01"
RESULTS = PLAN_DIR / "results"
EXPECTED_PLAN_SHA256 = "fc883444d2b23663d5f6034d773b76db948add29f1cda0c23384d17cfed2996a"
STACKS = (
    DATA / "wl_random_shape_stack_lens_reuse_2026-11-04_01/summary.json",
    DATA / "wl_random_shape_stack_lens_reuse_2026-11-04_02p5/summary.json",
)
DEFAULT_OUTPUT = DATA / "wl_random_shape_export_integrity_2026-11-05_01/summary.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def audit(output: Path) -> dict:
    if output.exists():
        raise FileExistsError(f"refusing to overwrite {output}")
    plan_path = PLAN_DIR / "random_shape_row_plan.json"
    state_path = RESULTS / "row_status.json"
    summary_path = RESULTS / "row_summary.json"
    plan = read_json(plan_path)
    state = read_json(state_path)
    summary = read_json(summary_path)
    plan_sha = sha256(plan_path)

    problems: list[str] = []
    if plan_sha != EXPECTED_PLAN_SHA256:
        problems.append("immutable row-plan hash changed")
    if state.get("status") != "complete" or state.get("plan_sha256") != plan_sha:
        problems.append("row checkpoint is incomplete or belongs to another plan")
    if summary.get("status") != "complete" or summary.get("plan_sha256") != plan_sha:
        problems.append("row summary is incomplete or belongs to another plan")
    for key in ("credentials_saved", "images_returned", "photoz_pdf_rows_included", "science_measurement"):
        if state.get(key) is not False:
            problems.append(f"row checkpoint scope flag {key} is not false")
    if state.get("source_shape_rows_downloaded") is not True:
        problems.append("row checkpoint does not record the approved shape-row export")

    jobs = state.get("jobs", {})
    expected_batches = int(plan.get("n_query_batches", -1))
    batch_records = []
    states = Counter(job.get("status") for job in jobs.values())
    total_rows = 0
    total_bytes = 0
    if len(jobs) != expected_batches or expected_batches != 124:
        problems.append(f"expected 124 jobs; found {len(jobs)}")
    for index in range(expected_batches):
        name = f"random_shape_rows_{index:03d}"
        job = jobs.get(name)
        path = RESULTS / f"{name}.csv"
        if not isinstance(job, dict) or job.get("status") != "complete":
            problems.append(f"{name} is not marked complete")
            continue
        if not path.is_file():
            problems.append(f"{name}.csv is missing")
            continue
        actual_bytes = path.stat().st_size
        actual_sha = sha256(path)
        expected_bytes = int(job.get("bytes", -1))
        expected_sha = str(job.get("sha256", ""))
        if actual_bytes != expected_bytes:
            problems.append(f"{name}.csv byte count differs from checkpoint")
        if actual_sha != expected_sha:
            problems.append(f"{name}.csv hash differs from checkpoint")
        total_rows += int(job.get("rows", -1))
        total_bytes += actual_bytes
        batch_records.append({
            "name": name,
            "job_id": str(job.get("job_id", "")),
            "status": job.get("status"),
            "rows": int(job.get("rows", -1)),
            "bytes": actual_bytes,
            "sha256": actual_sha,
        })

    expected_rows = int(summary.get("n_source_pairs", -1))
    expected_bytes = int(summary.get("n_csv_bytes", -1))
    if total_rows != expected_rows or expected_rows != int(plan.get("expected_random_source_pair_rows", -2)):
        problems.append("batch row total does not match the completed summary and immutable plan")
    if total_bytes != expected_bytes:
        problems.append("batch byte total does not match the completed summary")
    if summary.get("count_plan_sha256") != state.get("count_plan_sha256"):
        problems.append("count-plan hash differs between checkpoint and summary")

    stack_records = []
    batch_hashes = {row["name"]: row["sha256"] for row in batch_records}
    for path in STACKS:
        if not path.is_file():
            problems.append(f"saved stack summary is missing: {path}")
            continue
        stack = read_json(path)
        recorded_hashes = stack.get("input", {}).get("random_row_batch_sha256", {})
        agrees = recorded_hashes == batch_hashes
        if stack.get("input", {}).get("random_row_plan_sha256") != plan_sha:
            problems.append(f"saved stack references a different random-row plan: {path}")
        if not agrees:
            problems.append(f"saved stack batch hashes differ from current CSVs: {path}")
        stack_records.append({
            "path": str(path),
            "sha256": sha256(path),
            "status": stack.get("status"),
            "science_measurement": stack.get("science_measurement"),
            "random_row_plan_sha256": stack.get("input", {}).get("random_row_plan_sha256"),
            "n_random_centers": stack.get("input", {}).get("n_random_centers"),
            "n_source_pairs": stack.get("input", {}).get("n_source_pairs"),
            "all_batch_hashes_match": agrees,
            "block_width_deg": stack.get("input", {}).get("block_width_deg"),
        })
        if stack.get("science_measurement") is not False:
            problems.append(f"saved stack changed its science gate: {path}")

    result = {
        "schema": "astrafek-wl-random-shape-export-integrity-v1",
        "status": "complete" if not problems else "audit_incomplete",
        "science_measurement": False,
        "network_used": False,
        "credentials_read_or_saved": False,
        "images_returned": False,
        "photoz_pdf_rows_included": False,
        "plan": {"path": str(plan_path), "sha256": plan_sha},
        "checkpoint": {
            "path": str(state_path),
            "status": state.get("status"),
            "n_jobs": len(jobs),
            "job_states": dict(states),
            "stale_error_type_field": state.get("error_type"),
            "scope_flags_consistent": not any("scope flag" in p for p in problems),
        },
        "summary": {
            "path": str(summary_path),
            "status": summary.get("status"),
            "n_random_centers": summary.get("n_random_centers"),
            "n_matched_lenses": summary.get("n_matched_lenses"),
            "n_source_pairs": total_rows,
            "n_csv_bytes": total_bytes,
            "count_plan_sha256": summary.get("count_plan_sha256"),
        },
        "batches": batch_records,
        "stacks": stack_records,
        "problems": problems,
        "interpretation": (
            "All exported source-row files are hash-bound to both stack runs. "
            "The audit verifies file/provenance integrity only; it does not "
            "validate shape calibration, random-aperture matching, or the science claim."
        ),
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": result["status"],
        "science_measurement": False,
        "n_jobs": len(jobs),
        "job_states": dict(states),
        "n_source_pairs": total_rows,
        "n_csv_bytes": total_bytes,
        "stale_error_type_field": state.get("error_type"),
        "problems": problems,
        "output": str(output),
    }, ensure_ascii=False))
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    audit(args.output.resolve())


if __name__ == "__main__":
    main()
