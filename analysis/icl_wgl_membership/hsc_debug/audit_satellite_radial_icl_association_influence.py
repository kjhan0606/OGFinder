#!/usr/bin/env python3
"""Delete-one-cluster and delete-one-connected-block sensitivities for radial ICL links."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from analyze_satellite_outer_icl_baseline_association import metrics, read_tsv, sha256
from evaluate_prestack_pair_support import pair_block_labels


DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
DEFAULT_INPUT = DATA / "satellite_outer_icl_baseline_assoc_2026-10-20_07"
RADIAL_BINS = ("5-10", "10-15", "15-20", "20-30", "30-40", "40-50", "50-75")
METRICS = (
    "spearman_rho",
    "partial_pearson_r_control_z_richness",
    "partial_pearson_r_control_z_richness_satellite_count_area_support",
)


def finite(value):
    return value is not None and math.isfinite(float(value))


def summarize_deletions(point: dict, deletion_records: list[dict], scheme: str) -> dict:
    result = {}
    for label in ("icl_raw", "icl_outer_subtracted"):
        result[label] = {}
        for metric in METRICS:
            center = point[label].get(metric)
            values = np.asarray([r[label].get(metric) for r in deletion_records
                                 if finite(r[label].get(metric))], dtype=float)
            if not len(values) or not finite(center):
                result[label][metric] = {"n_valid": int(len(values)), "range": None,
                                         "n_sign_reversals": None}
                continue
            result[label][metric] = {
                "n_valid": int(len(values)),
                "range": [float(np.min(values)), float(np.max(values))],
                "n_sign_reversals": int(np.sum(np.sign(values) != np.sign(float(center)))),
                "all_same_sign_as_full_sample": bool(np.all(np.sign(values) == np.sign(float(center)))),
            }
    return {"scheme": scheme, "metrics": result}


def analyze_group(rows: list[dict], band: str, radial_bin: str) -> tuple[dict, list[dict], list[dict]]:
    if len(rows) < 12:
        raise ValueError(f"insufficient common-support rows for {band}/{radial_bin}: {len(rows)}")
    if len({r["cluster_id"] for r in rows}) != len(rows):
        raise ValueError(f"duplicate cluster in {band}/{radial_bin}")
    point = metrics(rows)
    block_targets = [{"ra_deg": r["ra_deg"], "dec_deg": r["dec_deg"]} for r in rows]
    block_controls = [{"ra_deg": r["control_ra_deg"], "dec_deg": r["control_dec_deg"]} for r in rows]
    _, labels, n_target, n_control, n_blocks, largest = pair_block_labels(
        block_targets, block_controls, np.ones(len(rows), dtype=bool))
    if n_blocks < 3:
        raise ValueError(f"at least three connected blocks required for {band}/{radial_bin}")
    block_deletions = []
    for label in sorted(set(labels)):
        deleted = [row for row, row_label in zip(rows, labels) if row_label != label]
        block_deletions.append({"deleted_block": label, "n_rows_remaining": len(deleted),
                                **metrics(deleted)})
    cluster_deletions = []
    for row in rows:
        deleted = [item for item in rows if item["cluster_id"] != row["cluster_id"]]
        cluster_deletions.append({"deleted_cluster_id": row["cluster_id"],
                                  "n_rows_remaining": len(deleted), **metrics(deleted)})
    result = {
        "filter": band, "radial_bin_kpc": radial_bin, "n_clusters": len(rows),
        "n_target_cells": int(n_target), "n_control_cells": int(n_control),
        "n_connected_blocks": int(n_blocks), "largest_block_clusters": int(largest),
        "point_estimates": point,
        "delete_one_block": summarize_deletions(point, block_deletions, "connected_5_degree_target_control_block"),
        "delete_one_cluster": summarize_deletions(point, cluster_deletions, "cluster"),
    }
    for deletion in block_deletions:
        deletion["filter"], deletion["radial_bin_kpc"] = band, radial_bin
    for deletion in cluster_deletions:
        deletion["filter"], deletion["radial_bin_kpc"] = band, radial_bin
    return result, block_deletions, cluster_deletions


def write_tsv(path: Path, rows: list[dict]) -> None:
    fields = ["filter", "radial_bin_kpc", "deleted_unit", "n_rows_remaining",
              "icl_raw_spearman_rho", "icl_outer_subtracted_spearman_rho",
              "icl_raw_partial_z_richness", "icl_outer_subtracted_partial_z_richness",
              "icl_raw_partial_extended", "icl_outer_subtracted_partial_extended"]
    with path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", extrasaction="raise")
        writer.writeheader()
        for row in rows:
            writer.writerow({
                "filter": row["filter"], "radial_bin_kpc": row["radial_bin_kpc"],
                "deleted_unit": row.get("deleted_block", row.get("deleted_cluster_id")),
                "n_rows_remaining": row["n_rows_remaining"],
                "icl_raw_spearman_rho": row["icl_raw"].get("spearman_rho"),
                "icl_outer_subtracted_spearman_rho": row["icl_outer_subtracted"].get("spearman_rho"),
                "icl_raw_partial_z_richness": row["icl_raw"].get("partial_pearson_r_control_z_richness"),
                "icl_outer_subtracted_partial_z_richness": row["icl_outer_subtracted"].get("partial_pearson_r_control_z_richness"),
                "icl_raw_partial_extended": row["icl_raw"].get("partial_pearson_r_control_z_richness_satellite_count_area_support"),
                "icl_outer_subtracted_partial_extended": row["icl_outer_subtracted"].get("partial_pearson_r_control_z_richness_satellite_count_area_support"),
            })


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    if args.output_dir.exists():
        parser.error(f"refusing to overwrite existing output: {args.output_dir}")
    summary_path = args.input_dir / "summary.json"
    detail_path = args.input_dir / "matched_radial_cluster_values.tsv"
    if not summary_path.is_file() or not detail_path.is_file():
        parser.error("radial association summary/detail input is missing")
    source = json.loads(summary_path.read_text(encoding="utf-8"))
    if (source.get("schema") != "astrafek-satellite-outer-icl-baseline-association-v1"
            or source.get("status") != "complete_exploratory"
            or source.get("science_measurement") is not False
            or len(source.get("radial_common_support_results", [])) != 14):
        parser.error("input is not the expected complete exploratory common-support analysis")
    expected_hash = source.get("outputs_sha256", {}).get("matched_radial_cluster_values.tsv")
    if not expected_hash or sha256(detail_path) != expected_hash:
        parser.error("radial detail checksum differs from its source summary")
    rows = read_tsv(detail_path)
    float_fields = ("icl_raw", "outer_baseline", "icl_outer_subtracted", "satellite",
                    "redshift", "richness", "ra_deg", "dec_deg", "control_ra_deg",
                    "control_dec_deg", "inner_coverage", "outer_coverage")
    for row in rows:
        for field in float_fields:
            row[field] = float(row[field])
        row["n_satellites"] = int(row["n_satellites"])
    if len(rows) != sum(item["n_clusters"] for item in source["results"]) * len(RADIAL_BINS):
        parser.error("radial detail row count is inconsistent with fixed samples")
    groups = {}
    for row in rows:
        groups.setdefault((row["filter"], row["radial_bin_kpc"]), []).append(row)
    results, block_rows, cluster_rows = [], [], []
    for band in ("i", "i2"):
        sample_by_radius = []
        expected_n = 77 if band == "i" else 183
        for radial_bin in RADIAL_BINS:
            group = groups.get((band, radial_bin), [])
            if len(group) != expected_n:
                parser.error(f"unexpected common sample size for {band}/{radial_bin}: {len(group)}")
            ids = {r["cluster_id"] for r in group}
            sample_by_radius.append(ids)
            result, blocks, clusters = analyze_group(group, band, radial_bin)
            results.append(result)
            block_rows.extend(blocks)
            cluster_rows.extend(clusters)
        if not all(ids == sample_by_radius[0] for ids in sample_by_radius[1:]):
            parser.error(f"cluster membership changes across radial bins for {band}")

    args.output_dir.mkdir(parents=True, exist_ok=False)
    block_path = args.output_dir / "leave_one_block_out.tsv"
    cluster_path = args.output_dir / "leave_one_cluster_out.tsv"
    write_tsv(block_path, block_rows)
    write_tsv(cluster_path, cluster_rows)
    result = {
        "schema": "astrafek-satellite-radial-icl-influence-audit-v1",
        "status": "complete_exploratory_sensitivity", "science_measurement": False,
        "method": {
            "sample": "same fixed cluster/filter IDs across all seven satellite radial bins",
            "deletions": "one cluster or one connected target/control 5-degree block removed at a time",
            "reported": "ranges and sign reversals of raw, outer-subtracted, and covariate-adjusted correlations; no p-values",
        },
        "limitations": [
            "leave-one-out ranges are influence diagnostics, not confidence intervals or hypothesis tests",
            "the 400-500 kpc annulus is not verified blank sky",
            "radial satellite-member lists can differ even though cluster IDs are fixed",
            "PSF wings, sky transfer, source completeness, extinction, and full coadd covariance remain unresolved",
            "the 14 filter-radius comparisons are exploratory and have no multiplicity correction",
        ],
        "source_summary": {"path": str(summary_path.resolve()), "sha256": sha256(summary_path)},
        "source_radial_table": {"path": str(detail_path.resolve()), "sha256": sha256(detail_path)},
        "analysis_script_sha256": sha256(Path(__file__).resolve()),
        "results": results,
        "outputs_sha256": {block_path.name: sha256(block_path), cluster_path.name: sha256(cluster_path)},
    }
    result_path = args.output_dir / "summary.json"
    with result_path.open("x", encoding="utf-8") as stream:
        json.dump(result, stream, ensure_ascii=False, indent=2, allow_nan=False)
        stream.write("\n")
    print(json.dumps({
        "status": result["status"], "science_measurement": False,
        "n_filter_radius_groups": len(results),
        "block_deletion_rows": len(block_rows), "cluster_deletion_rows": len(cluster_rows),
        "output": str(args.output_dir),
    }, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
