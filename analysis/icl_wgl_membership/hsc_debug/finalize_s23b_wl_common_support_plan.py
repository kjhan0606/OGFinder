#!/usr/bin/env python3
"""Finalize the local S23B WL source-row plan from approved count-only results.

This script validates the completed 400-center count-only result, updates a
new post-count row-plan manifest with observed control counts, and writes an
approval summary. It never contacts HSC CAS or retrieves shape rows.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import statistics
import sys

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from prepare_s23b_wl_common_support_plan import canonical_sha, sha, write_json
from prepare_wl_random_shape_rows import load_verified_counts

DEFAULT_PLAN_DIR = Path(
    "/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/"
    "wl_s23b_wl_common_support_plan_2026-11-10_01"
)


def quantile(values: list[float], probability: float) -> float:
    ordered = sorted(float(value) for value in values)
    if not ordered:
        raise ValueError("cannot calculate a quantile of an empty list")
    position = (len(ordered) - 1) * probability
    lower = int(position)
    upper = min(lower + 1, len(ordered) - 1)
    return ordered[lower] + (position - lower) * (ordered[upper] - ordered[lower])


def fmt_gib(value: int) -> str:
    return f"{int(value) / 1024**3:.2f} GiB"


def load_lens_counts(plan_dir: Path) -> dict[str, int]:
    with (plan_dir / "selected_lenses.tsv").open(encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    counts = {str(row["lens_id"]): int(row["n_raw_shape_sources"]) for row in rows}
    if len(counts) != 80 or any(value < 1000 for value in counts.values()):
        raise ValueError("selected lens table no longer matches the 80-lens raw-count screen")
    return counts


def finalize(plan_dir: Path) -> dict:
    plan_dir = Path(plan_dir).resolve()
    outputs = [
        "row_export_plan_post_count.json", "ROW_EXPORT_PLAN_POST_COUNT_SHA256",
        "post_count_assessment.json", "post_count_assessment.md",
        "post_count_approval_summary.json", "random_count_execution_summary.json",
    ]
    existing = [name for name in outputs if (plan_dir / name).exists()]
    if existing:
        raise FileExistsError(f"refusing to overwrite post-count outputs: {existing}")

    pre_row_path = plan_dir / "row_export_plan.json"
    pre_row = json.loads(pre_row_path.read_text(encoding="utf-8"))
    if sha(pre_row_path) != (plan_dir / "ROW_EXPORT_PLAN_SHA256").read_text().strip():
        raise ValueError("pre-count row-plan file changed after the previous approval stage")
    if pre_row.get("query_plan_sha256") != canonical_sha(pre_row):
        raise ValueError("pre-count row-plan canonical hash is invalid")

    count_manifest, centers, count_batches, counts, count_records, count_sha, count_summary = (
        load_verified_counts(plan_dir, plan_dir / "results")
    )
    if count_sha != pre_row.get("count_only_preflight_plan_sha256"):
        raise ValueError("completed counts do not match the plan used to build the row SQL")
    if (len(centers) != 400 or count_summary.get("n_matched_lenses") != 80
            or count_summary.get("status") != "complete"):
        raise ValueError("the approved count result does not cover the expected 400 controls")
    lens_counts = load_lens_counts(plan_dir)
    by_lens: dict[str, list[int]] = {lens_id: [] for lens_id in lens_counts}
    for center_id, row in counts.items():
        lens_id = row["matched_lens_id"]
        if lens_id not in by_lens:
            raise ValueError(f"count result references an unselected lens: {lens_id}")
        by_lens[lens_id].append(int(row["n_shape_sources"]))
    if any(len(values) != 5 for values in by_lens.values()):
        raise ValueError("the five-control assignment is incomplete for one or more lenses")

    lens_total = sum(lens_counts.values())
    control_total = sum(row["n_shape_sources"] for row in counts.values())
    if control_total != int(count_summary["sum_shape_source_counts"]):
        raise ValueError("verified control rows do not sum to the saved count summary")
    ratios = [(sum(values) / 5) / lens_counts[lens_id]
              for lens_id, values in by_lens.items()]
    control_sums_per_lens = [sum(values) for values in by_lens.values()]
    per_center_counts = [int(row["n_shape_sources"]) for row in counts.values()]
    n_zero = sum(value == 0 for value in per_center_counts)
    n_lenses_any_zero = sum(any(value == 0 for value in values) for values in by_lens.values())
    n_lenses_all_zero = sum(all(value == 0 for value in values) for values in by_lens.values())
    support = {
        "n_random_centers": len(centers), "n_matched_lenses": len(lens_counts),
        "controls_per_lens": 5,
        "random_source_pair_rows_from_count_only": control_total,
        "lens_source_pair_rows_from_count_only": lens_total,
        "combined_source_pair_rows_expected": lens_total + control_total,
        "zero_count_controls": n_zero,
        "positive_count_controls": len(per_center_counts) - n_zero,
        "control_count_min_median_max": [
            min(per_center_counts), statistics.median(per_center_counts), max(per_center_counts),
        ],
        "lens_sum_of_five_control_counts_min_median_max": [
            min(control_sums_per_lens), statistics.median(control_sums_per_lens),
            max(control_sums_per_lens),
        ],
        "lenses_with_at_least_one_zero_control": n_lenses_any_zero,
        "lenses_with_all_five_zero_controls": n_lenses_all_zero,
        "lenses_with_all_five_positive_controls": len(lens_counts) - n_lenses_any_zero,
        "per_lens_mean_control_to_lens_raw_count_ratio_q16_q50_q84": {
            "q16": quantile(ratios, 0.16), "q50": quantile(ratios, 0.50),
            "q84": quantile(ratios, 0.84),
        },
        "aggregate_mean_control_to_lens_raw_count_ratio": control_total / (5 * lens_total),
        "interpretation": "raw released-catalog availability only; project shape/redshift/B-mode/radial cuts have not been applied",
    }

    random_row_batches = [dict(batch) for batch in pre_row["random_batches"]]
    batch_size = int(count_manifest["batch_size"])
    for index, batch in enumerate(random_row_batches):
        chunk = centers[index * batch_size:(index + 1) * batch_size]
        expected = sum(counts[row["random_center_id"]]["n_shape_sources"] for row in chunk)
        if batch.get("n_random_centers") != len(chunk):
            raise ValueError(f"source-row SQL batch has a changed center count: {batch.get('name')}")
        batch["expected_rows_from_count_preflight"] = expected
        batch["count_preflight_batch_name"] = f"random_shape_counts_{index:03d}"
        batch["count_preflight_batch_sha256"] = count_records[index]["sha256"]
    if len(random_row_batches) != len(count_batches) or len(random_row_batches) != 16:
        raise ValueError("random row batches do not align with the completed count batches")

    old_estimates = pre_row["estimated_uncompressed_csv_bytes"]
    lens_bpr = float(old_estimates["point_basis_bytes_per_row"]["lens"])
    control_bpr = float(old_estimates["point_basis_bytes_per_row"]["random_control"])
    combined_total = lens_total + control_total
    estimated_bytes = {
        "low": combined_total * 256,
        "point": int(round(lens_total * lens_bpr + control_total * control_bpr)),
        "high": combined_total * 480,
        "low_gib": combined_total * 256 / 1024**3,
        "point_gib": (lens_total * lens_bpr + control_total * control_bpr) / 1024**3,
        "high_gib": combined_total * 480 / 1024**3,
        "planning_bytes_per_row_range": [256, 480],
        "point_basis_bytes_per_row": {"lens": lens_bpr, "random_control": control_bpr},
        "unit": "uncompressed CSV bytes",
        "basis": "observed count-only pair totals multiplied by historical measured bytes per pair; exact output bytes depend on returned field formatting",
    }

    row_plan = dict(pre_row)
    row_plan.pop("query_plan_sha256", None)
    row_plan["schema"] = "astrafek-s23b-wl-common-support-shape-row-plan-v2"
    row_plan["count_only_preflight_status"] = "complete"
    row_plan["count_only_preflight_plan_sha256"] = count_sha
    row_plan["count_only_preflight_summary_sha256"] = sha(plan_dir / "results" / "count_summary.json")
    row_plan["count_only_preflight_results_sha256"] = {
        record["name"]: record["sha256"] for record in count_records
    }
    row_plan["expected_random_source_pair_rows"] = control_total
    row_plan["expected_combined_source_pair_rows_from_count_only"] = combined_total
    row_plan["pre_count_random_source_pair_row_proxies_from_other_sample"] = (
        row_plan.pop("random_source_pair_rows_planning_proxies")
    )
    row_plan["pre_count_combined_pair_proxies_from_other_sample"] = (
        row_plan.pop("combined_source_pair_rows_planning_proxies")
    )
    row_plan["control_count_support_summary"] = support
    row_plan["estimated_uncompressed_csv_bytes"] = estimated_bytes
    row_plan["random_batches"] = random_row_batches
    row_plan["count_preflight_required_before_exact_control_row_total"] = False
    row_plan["count_preflight_completed"] = True
    row_plan["query_submitted"] = False
    row_plan["source_rows_downloaded"] = False
    row_plan["science_measurement"] = False
    historical_notes = [
        note for note in row_plan.get("notes", [])
        if not note.startswith("No source-row or count-only query is submitted")
        and not note.startswith("The 400-control count-only plan is a separate data request")
    ]
    row_plan["notes"] = historical_notes + [
        "The separately approved 400-center count-only query is complete and hash-pinned in this manifest.",
        "Expected shape-row pair counts are the verified count-only values; each returned source-row batch must still match every per-center count and unique (center_id, object_id) pair.",
        "This manifest does not authorize or submit the source-row query or transfer.",
    ]
    row_plan["query_plan_sha256"] = canonical_sha(row_plan)
    post_row_path = plan_dir / "row_export_plan_post_count.json"
    write_json(post_row_path, row_plan)
    post_row_sha = sha(post_row_path)
    (plan_dir / "ROW_EXPORT_PLAN_POST_COUNT_SHA256").write_text(post_row_sha + "\n", encoding="ascii")

    results_summary_path = plan_dir / "results" / "count_summary.json"
    assessment = {
        "schema": "astrafek-s23b-wl-common-support-post-count-assessment-v1",
        "status": "complete_count_only_audit_source_row_plan_prepared",
        "count_plan_sha256": count_sha,
        "count_summary_sha256": sha(results_summary_path),
        "count_plan_status": "complete",
        "count_batches_verified": len(count_records),
        "count_centers_verified": len(counts),
        "count_scope": count_manifest["output_columns"],
        "control_count_support": support,
        "expected_lens_source_pair_rows": lens_total,
        "expected_random_source_pair_rows": control_total,
        "expected_combined_source_pair_rows": combined_total,
        "estimated_uncompressed_csv_bytes": estimated_bytes,
        "lens_output_columns": row_plan["lens_output_columns"],
        "random_output_columns": row_plan["random_output_columns"],
        "row_export_plan_file_sha256": post_row_sha,
        "count_only_query_submitted": True,
        "source_row_query_submitted": False,
        "source_rows_downloaded": False,
        "shape_values_returned": False,
        "coordinates_returned": False,
        "photoz_returned": False,
        "images_returned": False,
        "network_used_for_count_only_query": True,
        "credentials_saved": False,
        "science_measurement": False,
        "execution_requires_explicit_user_approval": True,
        "limitations": [
            "Counts are raw released-catalog availability and precede project shape, source-redshift, B-mode, and radial cuts.",
            "The source-row query is not submitted. Returned row counts must still match per-center count-only values.",
            "The byte estimate uses historical bytes per pair; actual bytes are not known until a separately approved transfer.",
            "B-mode polarity, calibrated source n(z), lensing efficiency, final common support, and spatial covariance remain unresolved.",
        ],
    }
    write_json(plan_dir / "post_count_assessment.json", assessment)
    write_json(plan_dir / "random_count_execution_summary.json", {
        "schema": "astrafek-s23b-random-count-execution-summary-v1",
        "status": "complete", "plan_sha256": count_sha,
        "n_random_centers": len(counts), "n_matched_lenses": len(lens_counts),
        "n_verified_batches": len(count_records),
        "sum_shape_source_counts": control_total,
        "zero_source_centers": n_zero,
        "source_shape_rows_returned": False, "photoz_rows_returned": False,
        "images_returned": False, "science_measurement": False,
    })
    approval = {
        "schema": "astrafek-s23b-wl-source-row-approval-summary-v1",
        "status": "prepared_not_submitted",
        "action": "HSC PDR3 S19A weak-lensing source-row export for selected S23B lenses and matched random controls",
        "row_export_plan_file_sha256": post_row_sha,
        "count_only_plan_sha256": count_sha,
        "n_lenses": len(lens_counts), "n_random_centers": len(counts),
        "lens_source_pair_rows_expected": lens_total,
        "random_source_pair_rows_expected": control_total,
        "combined_source_pair_rows_expected": combined_total,
        "lens_row_sql_batches": len(row_plan["lens_batches"]),
        "random_row_sql_batches": len(random_row_batches),
        "lens_output_columns": row_plan["lens_output_columns"],
        "random_output_columns": row_plan["random_output_columns"],
        "estimated_uncompressed_csv_bytes": estimated_bytes,
        "photoz_pdf_rows_included": False, "images_included": False,
        "query_submitted": False, "source_rows_downloaded": False,
        "credentials_saved": False, "science_measurement": False,
        "requires_explicit_user_approval": True,
    }
    write_json(plan_dir / "post_count_approval_summary.json", approval)

    ratio = support["per_lens_mean_control_to_lens_raw_count_ratio_q16_q50_q84"]
    report = [
        "# S23B weak-lensing control count audit and source-row plan",
        "",
        "The separately approved count-only query completed 16/16 batches and all 400 center counts passed ID, integer, row-count, and SHA-256 checks. Its output contained only the random-center ID, paired lens ID, and integer raw source count. No shape values, source coordinates, photo-z, or images were returned. Credentials were not saved.",
        "",
        "## Observed raw support",
        "",
        f"The 400 controls contain {control_total:,} raw count pairs; the 80 selected lenses contain {lens_total:,}. The combined expected source-row count is {combined_total:,} lens/control pairs. Eighteen control centers have zero objects; {len(per_center_counts) - n_zero} are positive. No lens has all five controls at zero, while {n_lenses_any_zero} lenses have at least one zero-count control.",
        f"Per-center control counts have min/median/max {min(per_center_counts):,}/{statistics.median(per_center_counts):,.1f}/{max(per_center_counts):,}. Across the 80 lenses, the mean of each lens's five-control counts divided by its lens count has q16/median/q84 {ratio['q16']:.3f}/{ratio['q50']:.3f}/{ratio['q84']:.3f}. Aggregate random-per-center to lens count ratio is {support['aggregate_mean_control_to_lens_raw_count_ratio']:.3f}. These are pre-cut catalog availability diagnostics, not final matched source support.",
        "",
        "## Prepared source-row export (not submitted)",
        "",
        f"The post-count plan expects {lens_total:,} lens rows and {control_total:,} control rows, totaling {combined_total:,} lens/control pairs. Lens SQL is split into {len(row_plan['lens_batches'])} batches; control SQL into {len(random_row_batches)} batches. On return, every center's row total and each unique (center ID, object ID) pair must be checked against the count preflight.",
        f"The estimated uncompressed CSV volume is {fmt_gib(estimated_bytes['low'])} low, {fmt_gib(estimated_bytes['point'])} point, and {fmt_gib(estimated_bytes['high'])} high. Low/high assume 256/480 bytes per row; the point estimate applies historical measured 296.4-byte lens and 328.8-byte control widths to the observed counts. Actual file size is not yet known.",
        "",
        "Lens rows have 22 columns and control rows 23. Both preserve the shape components, measurement weights and shear-bias terms, source magnitudes/resolution, PSF moments, raw B-mode mask, and HSC Y3 source-redshift bin. Full photo-z PDFs and images are excluded.",
        "",
        "The row plan remains `prepared_not_submitted`; separate explicit approval tied to its exact file SHA-256 is required before any row query or transfer. Apply the same shape, source-redshift, B-mode, and radial cuts to lenses and controls. Verify B-mode polarity and calibrated n(z)/lensing efficiency and audit footprint, PSF, and covariance before calling this a measurement. Keep `science_measurement=false`.",
        "",
        f"Count-only plan SHA-256: `{count_sha}`  ",
        f"Source-row plan file SHA-256: `{post_row_sha}`  ",
        f"Plan directory: `{plan_dir}`  ",
        "No source-row query was submitted and no source-shape values were downloaded.",
        "",
    ]
    (plan_dir / "post_count_assessment.md").write_text("\n".join(report), encoding="utf-8")
    return {
        "status": assessment["status"], "count_plan_sha256": count_sha,
        "row_export_plan_file_sha256": post_row_sha,
        "n_lenses": len(lens_counts), "n_random_centers": len(counts),
        "expected_combined_source_pair_rows": combined_total,
        "estimated_uncompressed_csv_bytes": estimated_bytes,
        "count_only_query_submitted": True, "source_row_query_submitted": False,
        "source_rows_downloaded": False, "science_measurement": False,
        "output": str(plan_dir),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan-dir", type=Path, default=DEFAULT_PLAN_DIR)
    args = parser.parse_args()
    print(json.dumps(finalize(args.plan_dir), ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
