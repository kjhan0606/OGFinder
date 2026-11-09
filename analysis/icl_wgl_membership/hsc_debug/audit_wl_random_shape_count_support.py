#!/usr/bin/env python3
"""Audit completed count-only random-center batches without querying HSC."""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import json
from pathlib import Path

from prepare_wl_random_shape_count_plan import sha
from submit_wl_random_shape_counts import read_count_rows, validate_plan
from submit_wl_shape_counts import load_plan as load_lens_count_plan, read_count_csv


def quantile(values: list[float], probability: float) -> float | None:
    if not values:
        return None
    ordered = sorted(values)
    position = (len(ordered) - 1) * probability
    lower = int(position)
    upper = min(lower + 1, len(ordered) - 1)
    fraction = position - lower
    return ordered[lower] * (1 - fraction) + ordered[upper] * fraction


def audit(plan_dir: Path, results_dir: Path, lens_preflight_dir: Path | None = None) -> dict:
    plan_dir, results_dir = Path(plan_dir).resolve(), Path(results_dir).resolve()
    manifest, centers, batches, plan_sha = validate_plan(plan_dir)
    checkpoint_path = results_dir / "count_status.json"
    if not checkpoint_path.is_file():
        raise FileNotFoundError(f"missing count checkpoint: {checkpoint_path}")
    checkpoint = json.loads(checkpoint_path.read_text(encoding="utf-8"))
    if checkpoint.get("plan_sha256") != plan_sha:
        raise ValueError("count checkpoint does not match the validated plan")
    for key in ("credentials_saved", "source_shape_rows_returned", "photoz_rows_returned",
                "images_returned", "science_measurement"):
        if checkpoint.get(key) is not False:
            raise ValueError(f"checkpoint violates count-only scope: {key}")

    batch_size = int(manifest["batch_size"])
    jobs = checkpoint.get("jobs", {})
    counts: dict[str, dict] = {}
    verified_batches = []
    unverified_batches = []
    for index, batch in enumerate(batches):
        name = f"random_shape_counts_{index:03d}"
        item = jobs.get(name, {})
        path = results_dir / f"{name}.csv"
        start = index * batch_size
        chunk = centers[start:start + batch_size]
        if item.get("status") != "complete" or not path.is_file():
            unverified_batches.append({"name": name, "status": item.get("status", "missing"),
                                      "job_id": item.get("job_id")})
            continue
        digest = sha(path)
        if digest != item.get("sha256"):
            raise ValueError(f"completed CSV checksum mismatch: {name}")
        rows = read_count_rows(path, chunk)
        if item.get("rows") != len(rows) or len(rows) != len(chunk):
            raise ValueError(f"completed CSV row count mismatch: {name}")
        counts.update(rows)
        verified_batches.append({"name": name, "job_id": item["job_id"],
                                 "n_centers": len(rows), "sha256": digest,
                                 "sum_shape_source_counts": sum(
                                     row["n_shape_sources"] for row in rows.values())})

    expected_by_lens: dict[str, list[str]] = defaultdict(list)
    for center in centers:
        expected_by_lens[center["matched_lens_id"]].append(center["random_center_id"])
    observed_by_lens: dict[str, list[int]] = defaultdict(list)
    zero_centers = 0
    for center_id, row in counts.items():
        value = row["n_shape_sources"]
        observed_by_lens[row["matched_lens_id"]].append(value)
        zero_centers += value == 0
    complete_lens_counts = {
        lens_id: values for lens_id, values in observed_by_lens.items()
        if len(values) == len(expected_by_lens[lens_id])
    }
    positive_centers_per_complete_lens = Counter(
        sum(value > 0 for value in values) for values in complete_lens_counts.values())
    verified_count_total = sum(row["n_shape_sources"] for row in counts.values())
    run_complete = (checkpoint.get("status") == "complete"
                    and len(verified_batches) == len(batches)
                    and len(counts) == len(centers))
    if checkpoint.get("status") == "complete" and not run_complete:
        raise ValueError("checkpoint says complete but not all planned count rows verified")
    result = {
        "schema": "astrafek-wl-random-shape-count-support-audit-v1",
        "status": "complete" if run_complete else "partial",
        "science_measurement": False,
        "count_only": True,
        "plan_sha256": plan_sha,
        "checkpoint_path": str(checkpoint_path),
        "n_expected_random_centers": len(centers),
        "n_expected_matched_lenses": len(expected_by_lens),
        "n_expected_batches": len(batches),
        "n_verified_batches": len(verified_batches),
        "n_unverified_batches": len(unverified_batches),
        "n_verified_centers": len(counts),
        "n_complete_lenses": len(complete_lens_counts),
        "n_verified_count_rows": len(counts),
        "sum_shape_source_counts_verified": verified_count_total,
        "n_zero_count_centers_verified": zero_centers,
        "positive_centers_per_complete_lens": {
            str(key): positive_centers_per_complete_lens[key]
            for key in sorted(positive_centers_per_complete_lens)},
        "verified_batches": verified_batches,
        "unverified_batches": unverified_batches,
        "source_shape_rows_returned": False,
        "photoz_rows_returned": False,
        "images_returned": False,
    }
    if lens_preflight_dir is not None:
        lens_preflight_dir = Path(lens_preflight_dir).resolve()
        lens_manifest = load_lens_count_plan(lens_preflight_dir)
        lens_results = lens_preflight_dir / "cas_counts"
        lens_checkpoint_path = lens_results / "count_status.json"
        lens_summary_path = lens_results / "count_summary.json"
        if not lens_checkpoint_path.is_file() or not lens_summary_path.is_file():
            raise FileNotFoundError("lens count preflight checkpoint/summary is missing")
        lens_checkpoint = json.loads(lens_checkpoint_path.read_text(encoding="utf-8"))
        lens_summary = json.loads(lens_summary_path.read_text(encoding="utf-8"))
        if lens_checkpoint.get("status") != "complete" or lens_summary.get("status") != "complete":
            raise ValueError("lens count preflight is not complete")
        if (lens_checkpoint.get("source_rows_downloaded") is not False
                or lens_checkpoint.get("science_measurement") is not False
                or lens_summary.get("source_rows_downloaded") is not False
                or lens_summary.get("science_measurement") is not False):
            raise ValueError("lens preflight declares data outside its count-only scope")
        lens_ids = [str(row["lens_id"]) for row in lens_manifest["lenses"]]
        lens_batch_size = int(lens_manifest["batch_size"])
        lens_counts = {}
        lens_jobs = lens_checkpoint.get("jobs", {})
        for index, _ in enumerate(lens_manifest["sql_files"]):
            name = f"shape_counts_{index:03d}"
            item = lens_jobs.get(name, {})
            path = lens_results / f"{name}.csv"
            chunk = lens_ids[index * lens_batch_size:(index + 1) * lens_batch_size]
            if item.get("status") != "complete" or not path.is_file() or sha(path) != item.get("sha256"):
                raise ValueError(f"lens count batch missing or checksum mismatch: {name}")
            lens_counts.update({row["lens_id"]: row["n_shape_sources"]
                                for row in read_count_csv(path, chunk)})
        if (len(lens_counts) != len(lens_ids)
                or len(lens_counts) != int(lens_summary["n_count_rows"])
                or sum(lens_counts.values()) != int(lens_summary["sum_lens_cone_shape_counts"])
                or lens_checkpoint.get("plan_sha256") != lens_summary.get("query_plan_sha256")):
            raise ValueError("verified lens counts do not match the complete summary")

        lens_support: dict[str, list[int]] = defaultdict(list)
        positive_lens_random_to_lens_ratios = []
        for lens_id, values in complete_lens_counts.items():
            if lens_id in lens_counts:
                lens_support["lens_zero" if lens_counts[lens_id] == 0 else "lens_positive"].append(
                    sum(value > 0 for value in values))
                if lens_counts[lens_id] > 0:
                    positive_lens_random_to_lens_ratios.append(
                        (sum(values) / len(values)) / lens_counts[lens_id])
        result["lens_support_comparison"] = {
            "preflight_dir": str(lens_preflight_dir),
            "preflight_plan_sha256": lens_checkpoint["plan_sha256"],
            "n_verified_lens_counts": len(lens_counts),
            "n_complete_random_matched_lenses_with_lens_count": sum(map(len, lens_support.values())),
            "lens_preflight_zero_counts": sum(value == 0 for value in lens_counts.values()),
            "complete_random_lenses_by_lens_preflight_and_positive_random_centers": {
                key: {str(k): v for k, v in sorted(Counter(values).items())}
                for key, values in sorted(lens_support.items())},
            "positive_lens_random_mean_to_lens_count_ratio_q16_q50_q84": {
                key: quantile(positive_lens_random_to_lens_ratios, probability)
                for key, probability in (("q16", 0.16), ("q50", 0.50), ("q84", 0.84))},
        }
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan-dir", required=True, type=Path)
    parser.add_argument("--results-dir", required=True, type=Path)
    parser.add_argument("--lens-preflight-dir", type=Path,
                        help="optional verified 800-lens count-preflight directory")
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    output = args.output.resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite audit output: {output}")
    result = audit(args.plan_dir, args.results_dir, args.lens_preflight_dir)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({key: result[key] for key in (
        "status", "n_verified_batches", "n_verified_centers", "n_complete_lenses",
        "sum_shape_source_counts_verified", "n_zero_count_centers_verified",
        "positive_centers_per_complete_lens", "science_measurement")}, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
