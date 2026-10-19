#!/usr/bin/env python3
"""Estimate random-lens count and row-export scale from the approved lens export.

The calculation assumes a random-center shape-row rate and bytes-per-row equal
to the measured 800-lens export. It is a planning approximation only; count-only
results at selected random centers must replace it before any row export.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path


def estimate(row_summary_path: Path, row_status_path: Path,
             ratios=(1, 2, 4, 8, 10), count_batch_size=25):
    summary = json.loads(Path(row_summary_path).read_text(encoding="utf-8"))
    status = json.loads(Path(row_status_path).read_text(encoding="utf-8"))
    if (summary.get("schema") != "hsc-s19a-shape-row-summary-v1"
            or summary.get("status") != "complete"
            or summary.get("science_measurement") is not False
            or status.get("status") != "complete"
            or status.get("source_rows_downloaded") is not True
            or status.get("credentials_saved") is not False):
        raise ValueError("approved lens row export is incomplete or has unexpected provenance")
    n_lenses = int(summary["n_lenses"])
    n_pairs = int(summary["n_pair_rows"])
    if n_lenses < 1 or n_pairs < 1 or status.get("total_pair_rows") != n_pairs:
        raise ValueError("row summary and export checkpoint totals disagree")
    jobs = status.get("jobs", {})
    completed = [job for job in jobs.values() if job.get("status") == "complete"]
    if (not jobs or len(completed) != len(jobs)
            or sum(int(job.get("rows", -1)) for job in completed) != n_pairs):
        raise ValueError("row export batch ledger is incomplete or inconsistent")
    byte_total = sum(int(job.get("bytes", -1)) for job in completed)
    if byte_total <= 0 or any(int(job.get("bytes", 0)) <= 0 for job in completed):
        raise ValueError("row export byte ledger is invalid")
    if count_batch_size < 1:
        raise ValueError("count batch size must be positive")

    unique_ratios = sorted(set(int(ratio) for ratio in ratios))
    if not unique_ratios or any(ratio < 1 for ratio in unique_ratios):
        raise ValueError("ratios must be positive integers")
    rows = []
    for ratio in unique_ratios:
        n_centers = n_lenses * ratio
        expected_rows = n_pairs * ratio
        expected_bytes = byte_total * ratio
        rows.append({
            "random_centers_per_lens": ratio,
            "random_centers_total": n_centers,
            "count_only_batches": math.ceil(n_centers / count_batch_size),
            "expected_shape_pair_rows": expected_rows,
            "estimated_csv_bytes": expected_bytes,
            "estimated_csv_gib": expected_bytes / 1024**3,
            "random_to_lens_shape_variance_ratio_if_equal_independent_centers": 1.0 / ratio,
            "random_share_of_lens_plus_random_shape_variance_if_equal_independent_centers": 1.0 / (ratio + 1.0),
        })
    return {
        "schema": "astrafek-wl-random-center-scale-plan-v1",
        "status": "planning_estimate_only",
        "science_measurement": False,
        "network_used": False,
        "query_submitted": False,
        "lens_export": {
            "row_summary": str(Path(row_summary_path).resolve()),
            "row_status": str(Path(row_status_path).resolve()),
            "n_lenses": n_lenses, "n_pair_rows": n_pairs,
            "csv_bytes": byte_total,
            "mean_bytes_per_pair": byte_total / n_pairs,
            "mean_pairs_per_lens": n_pairs / n_lenses,
        },
        "count_batch_size": count_batch_size,
        "scenarios": rows,
        "assumptions_and_limits": [
            "Random-center source counts and row bytes are assumed equal to the measured cluster-lens export; verify with count-only queries before row export.",
            "Variance ratios assume equal independent per-center shape variance; repeated sources, overlapping apertures, depth, and spatial covariance can change them.",
            "A ratio is a planning scenario, not a selected or science-approved random catalog size.",
        ],
        "next_gate": "geometry-screen the approved random coordinates, choose a provisional ratio, prepare exact count-only SQL, and obtain approval before submitting it",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--row-summary", type=Path, required=True)
    parser.add_argument("--row-status", type=Path, required=True)
    parser.add_argument("--ratios", type=int, nargs="+", default=[1, 2, 4, 8, 10])
    parser.add_argument("--count-batch-size", type=int, default=25)
    parser.add_argument("--output", type=Path, help="optional new JSON output; refuses overwrite")
    args = parser.parse_args()
    result = estimate(args.row_summary, args.row_status, args.ratios, args.count_batch_size)
    rendered = json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False) + "\n"
    if args.output:
        output = args.output.resolve()
        if output.exists():
            parser.error(f"refusing to overwrite {output}")
        output.parent.mkdir(parents=True, exist_ok=True)
        output.write_text(rendered, encoding="utf-8")
    print(rendered, end="")


if __name__ == "__main__":
    main()
