#!/usr/bin/env python3
"""Audit HSC tract-region influence on matched target-control ICL profiles."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from run_wl_hsc_tract_region_influence import (
    nearest_region_labels,
    tract_centers_from_patch_qa,
    tract_region_labels,
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path):
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def _key(row):
    return (str(row["target_id"]), str(row["filter"]),
            int(float(row["r_min_kpc"])), int(float(row["r_max_kpc"])))


def leave_one_region_out_means(values, regions):
    """Return the full equal-row mean and per-region deletion means."""
    values = np.asarray(values, dtype=float)
    regions = np.asarray(regions)
    if (values.ndim != 1 or regions.shape != values.shape or len(values) < 2
            or not np.isfinite(values).all()):
        raise ValueError("values and region labels must be matching finite vectors")
    unique = np.unique(regions)
    if len(unique) < 2:
        raise ValueError("leave-one-region-out requires at least two regions")
    full_mean = float(values.mean())
    deletions = []
    for region in unique:
        keep = regions != region
        if not np.any(keep):
            raise ValueError("region deletion leaves no rows")
        deletions.append((region.item() if hasattr(region, "item") else region,
                          int(np.count_nonzero(keep)), float(values[keep].mean())))
    return full_mean, deletions


def audit_influence(stack_dir: Path, lens_manifest_path: Path, patch_qa_path: Path,
                    *, adjacency_deg=2.0):
    stack_summary_path = stack_dir / "summary.json"
    evaluation_path = stack_dir / "evaluation.json"
    stack_summary = json.loads(stack_summary_path.read_text(encoding="utf-8"))
    evaluation = json.loads(evaluation_path.read_text(encoding="utf-8"))
    if stack_summary.get("status") != "complete_exploratory" or stack_summary.get("science_measurement") is not False:
        raise ValueError("source ICL stack must be a completed exploratory product")
    paired_path = Path(stack_summary["outputs"]["paired_profiles"])
    expected_paired_sha = stack_summary["outputs"]["paired_profiles_sha256"]
    if paired_path.resolve().parent != stack_dir.resolve() or sha256(paired_path) != expected_paired_sha:
        raise ValueError("paired profile table is outside the stack or its checksum changed")
    broad_pair_path = Path(evaluation["broad_radial_pair_profiles"])
    broad_profile_path = Path(evaluation["broad_radial_profiles"])
    if any(path.resolve().parent != stack_dir.resolve() for path in (broad_pair_path, broad_profile_path)):
        raise ValueError("broad profile inputs must be inside the source stack directory")

    manifest = json.loads(lens_manifest_path.read_text(encoding="utf-8"))
    lenses = manifest.get("clusters", [])
    if len(lenses) != 800:
        raise ValueError(f"expected 800 target clusters, found {len(lenses)}")
    lens_by_id = {str(row["id"]): row for row in lenses}
    if len(lens_by_id) != 800:
        raise ValueError("target manifest has duplicate IDs")

    pair_rows = read_tsv(paired_path)
    if len(pair_rows) != 783:
        raise ValueError(f"expected 783 matched target-control pairs, found {len(pair_rows)}")
    pair_by_target_filter = {}
    target_coords = {}
    for row in pair_rows:
        target_id, filter_name = str(row["target_id"]), str(row["filter"])
        if target_id not in lens_by_id:
            raise ValueError(f"unknown target ID in paired profiles: {target_id}")
        key = (target_id, filter_name)
        if key in pair_by_target_filter:
            raise ValueError(f"duplicate target/filter pair: {key}")
        ra, dec = float(row["ra_deg"]), float(row["dec_deg"])
        expected = lens_by_id[target_id]
        ra_delta = (ra - float(expected["ra_deg"]) + 180.0) % 360.0 - 180.0
        dec_delta = dec - float(expected["dec_deg"])
        if max(abs(ra_delta), abs(dec_delta)) > 1e-7:
            raise ValueError(f"target coordinates disagree with manifest: {target_id}")
        pair_by_target_filter[key] = row
        target_coords[target_id] = (ra, dec)

    keys, tract_centers = tract_centers_from_patch_qa(patch_qa_path)
    base_regions = tract_region_labels(tract_centers, adjacency_deg)
    manifest_region = nearest_region_labels(
        [row["ra_deg"] for row in lenses], [row["dec_deg"] for row in lenses],
        tract_centers, base_regions,
    )
    manifest_region_by_id = {str(row["id"]): int(region)
                             for row, region in zip(lenses, manifest_region)}
    pair_region = {target_id: manifest_region_by_id[target_id]
                   for target_id in target_coords}

    stability = []
    for threshold in (1.6, 1.8, 2.0, 2.2, 2.5, 3.0):
        regions = tract_region_labels(tract_centers, threshold)
        assigned = nearest_region_labels(
            [row["ra_deg"] for row in lenses], [row["dec_deg"] for row in lenses],
            tract_centers, regions,
        )
        _, counts = np.unique(assigned, return_counts=True)
        stability.append({
            "adjacency_deg": threshold,
            "n_tract_components": int(len(np.unique(regions))),
            "occupied_lens_component_counts": sorted(map(int, counts), reverse=True),
        })
    expected_counts = stability[0]["occupied_lens_component_counts"]
    if any(item["n_tract_components"] != stability[0]["n_tract_components"]
           or item["occupied_lens_component_counts"] != expected_counts for item in stability):
        raise ValueError("tract-region partition is not stable across the configured threshold grid")

    broad_rows = read_tsv(broad_pair_path)
    reference_rows = read_tsv(broad_profile_path)
    if len(broad_rows) != 4698:
        raise ValueError(f"expected 4698 pair-by-annulus rows, found {len(broad_rows)}")
    broad_keys = [_key(row) for row in broad_rows]
    if len(set(broad_keys)) != len(broad_keys):
        raise ValueError("duplicate target/filter/annulus row in broad paired profiles")
    expected_pairs = set(pair_by_target_filter)
    if {key[:2] for key in broad_keys} != expected_pairs:
        raise ValueError("broad profile rows do not match the 783 paired-profile entries")

    metrics = []
    for ref in reference_rows:
        filter_name = str(ref["filter"])
        r_min = int(float(ref["r_min_kpc"]))
        r_max = int(float(ref["r_max_kpc"]))
        subset = [row for row in broad_rows
                  if row["filter"] == filter_name
                  and int(float(row["r_min_kpc"])) == r_min
                  and int(float(row["r_max_kpc"])) == r_max
                  and row["delta_njy_arcsec2"] not in ("", "nan", "NaN", "None")]
        values = np.asarray([float(row["delta_njy_arcsec2"]) for row in subset], dtype=float)
        coverage = np.asarray([float(row["common_coverage_fraction"]) for row in subset], dtype=float)
        if not len(values) or not np.isfinite(values).all() or not np.isfinite(coverage).all():
            raise ValueError(f"invalid or empty pair profile in {filter_name} {r_min}-{r_max} kpc")
        row_regions = np.asarray([pair_region[row["target_id"]] for row in subset])
        full_mean, deletion_means = leave_one_region_out_means(values, row_regions)
        if len(subset) != int(ref["n_pairs"]) or not math.isclose(
                full_mean, float(ref["mean_delta_njy_arcsec2"]), rel_tol=1e-11, abs_tol=1e-10):
            raise ValueError(f"pairwise reconstruction disagrees with stacked profile: {filter_name} {r_min}-{r_max}")
        loo = []
        for region, n_retained, loo_mean in deletion_means:
            shift = loo_mean - full_mean
            loo.append({
                "region_id": int(region),
                "n_pairs_retained": n_retained,
                "mean_delta_njy_arcsec2": loo_mean,
                "shift_njy_arcsec2": shift,
                "fractional_shift": shift / full_mean if full_mean else None,
                "sign_flip": bool(full_mean != 0 and loo_mean != 0 and np.sign(full_mean) != np.sign(loo_mean)),
            })
        metrics.append({
            "filter": filter_name,
            "r_min_kpc": r_min,
            "r_max_kpc": r_max,
            "n_pairs": len(values),
            "mean_delta_njy_arcsec2": full_mean,
            "pair_sd_njy_arcsec2": float(values.std(ddof=1)) if len(values) > 1 else None,
            "pair_sem_njy_arcsec2": float(values.std(ddof=1) / np.sqrt(len(values))) if len(values) > 1 else None,
            "median_common_coverage": float(np.median(coverage)),
            "leave_one_region_out": loo,
        })

    pair_region_counts = {}
    for target_id, region in pair_region.items():
        pair_region_counts[str(region)] = pair_region_counts.get(str(region), 0) + 1
    return {
        "schema": "astrafek-icl-tract-region-influence-v1",
        "status": "exploratory_influence",
        "science_measurement": False,
        "covariance_estimated": False,
        "method": "equal-pair-weight means reconstructed from saved broad target-control annular profiles; delete all pairs in one approximate connected HSC tract region",
        "adjacency_deg": float(adjacency_deg),
        "n_tracts": len(keys),
        "n_tract_components": int(len(np.unique(base_regions))),
        "n_occupied_pair_regions": int(len(set(pair_region.values()))),
        "manifest_lens_counts_by_region": {
            str(region): int(count) for region, count in zip(*np.unique(manifest_region, return_counts=True))
        },
        "matched_pair_target_counts_by_region": pair_region_counts,
        "partition_threshold_sensitivity": stability,
        "n_matched_pairs": len(pair_by_target_filter),
        "profile_metrics": metrics,
        "limitations": [
            "tract-center graph is an approximate footprint grouping, not a tract-polygon overlap or official HSC field map",
            "only three regions are occupied, so leave-one-region shifts are not a stable covariance, uncertainty interval, or significance test",
            "equal-pair-weight annular profiles inherit the original low common-pixel support and mask/sky/PSF limitations",
            "same-filter i and i2 results are kept separate; their response functions are not pooled",
            "no new HSC query or download was performed",
        ],
        "input": {
            "stack_dir": str(stack_dir.resolve()),
            "stack_summary_sha256": sha256(stack_summary_path),
            "evaluation_sha256": sha256(evaluation_path),
            "paired_profiles_sha256": sha256(paired_path),
            "broad_pair_profiles_sha256": sha256(broad_pair_path),
            "broad_profiles_sha256": sha256(broad_profile_path),
            "lens_manifest_path": str(lens_manifest_path.resolve()),
            "lens_manifest_sha256": sha256(lens_manifest_path),
            "patch_qa_path": str(patch_qa_path.resolve()),
            "patch_qa_sha256": sha256(patch_qa_path),
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stack-dir", type=Path, required=True)
    parser.add_argument("--lens-manifest", type=Path, required=True)
    parser.add_argument("--patch-qa", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--adjacency-deg", type=float, default=2.0)
    args = parser.parse_args()
    if args.output.exists():
        parser.error(f"refusing to overwrite existing output: {args.output}")
    result = audit_influence(
        args.stack_dir.resolve(), args.lens_manifest.resolve(), args.patch_qa.resolve(),
        adjacency_deg=args.adjacency_deg,
    )
    output = args.output.resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("x", encoding="utf-8") as stream:
        json.dump(result, stream, indent=2, ensure_ascii=False, allow_nan=False)
        stream.write("\n")
    print(json.dumps({
        "status": result["status"], "science_measurement": False,
        "n_matched_pairs": result["n_matched_pairs"],
        "manifest_lens_counts_by_region": result["manifest_lens_counts_by_region"],
        "matched_pair_target_counts_by_region": result["matched_pair_target_counts_by_region"],
        "profile_metrics": len(result["profile_metrics"]), "output": str(output),
    }, ensure_ascii=False))


if __name__ == "__main__":
    main()
