#!/usr/bin/env python3
"""Audit connected HSC tract partition sensitivity for the paired WL stack.

This small offline audit checks whether plausible tract-center adjacency
thresholds change the connected-component membership used by the saved
random-subtracted weak-lensing influence analysis. It reads no shape rows and
submits no HSC query.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from fetch_ancillary_catalogues import sha
from run_wl_hsc_tract_region_influence import (
    nearest_region_labels,
    tract_centers_from_patch_qa,
    tract_region_labels,
)


ROOT = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
DEFAULT_PATCH_QA = ROOT / "catalogues_2026-10-13_01/patch_qa.csv"
DEFAULT_REGION_SUMMARY = ROOT / "wl_random_tract_region_influence_2026-11-06_01/summary.json"
DEFAULT_RANDOM_CENTERS = ROOT / "wl_random_centers_lens_reuse_2026-11-03_01/random_centers.json"
DEFAULT_LENS_MANIFEST = ROOT / "wl_shape_preflight_2026-10-18_01/preflight_manifest.json"
DEFAULT_OUTPUT = ROOT / "wl_random_tract_partition_sensitivity_2026-11-06_01"
THRESHOLDS_DEG = (1.6, 2.0, 2.4, 3.0)


def partition_signature(keys, labels):
    groups = {}
    for key, label in zip(keys, labels):
        groups.setdefault(int(label), []).append(str(key))
    return tuple(sorted(tuple(sorted(group)) for group in groups.values()))


def run(args):
    output = args.output.resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite output directory: {output}")
    region_summary = json.loads(args.region_summary.read_text(encoding="utf-8"))
    random_doc = json.loads(args.random_centers.read_text(encoding="utf-8"))
    lens_manifest = json.loads(args.lens_manifest.read_text(encoding="utf-8"))
    common_lens_ids = {
        str(value) for deletion in region_summary["region_deletions"]
        for value in deletion["deleted_lens_ids"]
    }
    lens_rows = [row for row in lens_manifest["lenses"]
                 if str(row["lens_id"]) in common_lens_ids]
    centers = [row for row in random_doc["centers"]
               if str(row["matched_lens_id"]) in common_lens_ids]
    if (len(common_lens_ids) != 435 or len(lens_rows) != 435 or len(centers) != 2175
            or len({row["random_center_id"] for row in centers}) != 2175):
        raise ValueError("saved common-support lens/random-center manifests are incomplete")

    keys, xyz = tract_centers_from_patch_qa(args.patch_qa)
    rows = []
    reference_signature = None
    for threshold in THRESHOLDS_DEG:
        labels = tract_region_labels(xyz, threshold)
        signature = partition_signature(keys, labels)
        if reference_signature is None:
            reference_signature = signature
        same_partition = signature == reference_signature
        lens_regions = nearest_region_labels(
            [row["ra_deg"] for row in lens_rows],
            [row["dec_deg"] for row in lens_rows], xyz, labels,
        )
        center_regions = nearest_region_labels(
            [row["ra_deg"] for row in centers],
            [row["dec_deg"] for row in centers], xyz, labels,
        )
        lens_counts = {str(int(region)): int((lens_regions == region).sum())
                       for region in sorted(set(lens_regions.tolist()))}
        random_counts = {str(int(region)): int((center_regions == region).sum())
                         for region in sorted(set(center_regions.tolist()))}
        rows.append({
            "adjacency_deg": threshold,
            "n_tract_components": len(signature),
            "same_tract_partition_as_lowest_threshold": same_partition,
            "common_lens_counts_by_component": lens_counts,
            "random_center_counts_by_component": random_counts,
        })
    if not all(row["same_tract_partition_as_lowest_threshold"] for row in rows):
        raise ValueError("tract topology changes over the reviewed adjacency thresholds")
    if len({json.dumps(row["common_lens_counts_by_component"], sort_keys=True)
            for row in rows}) != 1:
        raise ValueError("lens component assignment changes despite an invariant partition")
    if len({json.dumps(row["random_center_counts_by_component"], sort_keys=True)
            for row in rows}) != 1:
        raise ValueError("random-center component assignment changes despite an invariant partition")

    script_hash = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    result = {
        "schema": "astrafek-wl-random-tract-partition-sensitivity-v1",
        "status": "complete_exploratory",
        "science_measurement": False,
        "network_used": False,
        "shape_rows_read": False,
        "query_submitted": False,
        "method": "compare exact tract-ID connected-component partitions and nearest-tract lens/random-center labels over the specified center-adjacency thresholds",
        "thresholds_deg": list(THRESHOLDS_DEG),
        "partition_invariant": True,
        "n_tracts": len(keys),
        "n_tract_components": len(reference_signature),
        "common_support_lenses": len(common_lens_ids),
        "common_support_random_centers": len(centers),
        "threshold_results": rows,
        "interpretation_limit": "topology stability does not validate tract polygons, nearest-center boundaries, official field labels, or a field covariance",
        "inputs": {
            "region_influence_summary_sha256": sha(args.region_summary),
            "patch_qa_sha256": sha(args.patch_qa),
            "random_center_manifest_sha256": sha(args.random_centers),
            "lens_preflight_manifest_sha256": sha(args.lens_manifest),
            "runner_sha256": script_hash,
        },
    }
    output.mkdir(parents=True, exist_ok=False)
    with (output / "summary.json").open("x", encoding="utf-8") as stream:
        json.dump(result, stream, ensure_ascii=False, indent=2, allow_nan=False)
        stream.write("\n")
    print(json.dumps({
        "status": result["status"], "partition_invariant": True,
        "thresholds_deg": list(THRESHOLDS_DEG), "tract_components": len(reference_signature),
        "common_lenses": len(common_lens_ids), "random_centers": len(centers),
        "output": str(output / "summary.json"),
    }, ensure_ascii=False))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--patch-qa", type=Path, default=DEFAULT_PATCH_QA)
    parser.add_argument("--region-summary", type=Path, default=DEFAULT_REGION_SUMMARY)
    parser.add_argument("--random-centers", type=Path, default=DEFAULT_RANDOM_CENTERS)
    parser.add_argument("--lens-manifest", type=Path, default=DEFAULT_LENS_MANIFEST)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    run(parser.parse_args())


if __name__ == "__main__":
    main()
