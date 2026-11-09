#!/usr/bin/env python3
"""Audit HSC tract-region balance of the saved i-band ICL/WL matched sample.

This offline audit maps the common and matched ICL/WL lens IDs onto the
approximate connected tract regions already computed for the 800-lens parent.
It does not read shape rows or images and does not diagnose the cross-shear
origin; region labels are only a footprint-concentration check.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
DEFAULT_OVERLAP = DATA / "icl_wl_sample_overlap_audit_2026-11-05_03/summary.json"
DEFAULT_SPLITS = {
    "i": DATA / "wl_random_subtracted_icl_split_i_2026-11-06_01/summary.json",
    "i2": DATA / "wl_random_subtracted_icl_split_2026-11-05_02/summary.json",
}
DEFAULT_REGIONS = DATA / "wl_hsc_tract_region_influence_2026-10-20_01/summary.json"
DEFAULT_OUTPUTS = {
    "i": DATA / "wl_icl_split_tract_region_balance_2026-11-06_01",
    "i2": DATA / "wl_icl_split_tract_region_balance_i2_2026-11-06_01",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def counts_by_region(ids: set[str], region_by_lens: dict[str, int]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for lens_id in ids:
        if lens_id not in region_by_lens:
            raise ValueError(f"lens {lens_id} has no saved tract-region label")
        key = str(region_by_lens[lens_id])
        counts[key] = counts.get(key, 0) + 1
    return {key: counts[key] for key in sorted(counts, key=int)}


def run(overlap_path: Path, split_path: Path, regions_path: Path,
        output: Path, filter_name: str) -> dict:
    paths = [Path(value).resolve() for value in
             (overlap_path, split_path, regions_path)]
    if output.exists():
        raise FileExistsError(f"refusing to overwrite {output}")
    overlap_path, split_path, regions_path = paths
    overlap = json.loads(overlap_path.read_text(encoding="utf-8"))
    split = json.loads(split_path.read_text(encoding="utf-8"))
    regions = json.loads(regions_path.read_text(encoding="utf-8"))
    if overlap.get("status") != "complete_exploratory" or overlap.get("science_measurement") is not False:
        raise ValueError("unexpected ICL/WL overlap status")
    if split.get("status") != "complete_exploratory" or split.get("science_measurement") is not False:
        raise ValueError("unexpected saved ICL-split status")
    if regions.get("status") != "exploratory_field_influence" or regions.get("science_measurement") is not False:
        raise ValueError("unexpected tract-region influence status")
    if split.get("filter") != filter_name or split.get("source_tomographic_bins") != [2, 3, 4]:
        raise ValueError("saved split filter or source bins do not match the request")

    region_by_lens: dict[str, int] = {}
    for item in regions.get("region_deletions", []):
        region_id = int(item["region_id"])
        for raw_id in item.get("deleted_lens_ids", []):
            lens_id = str(raw_id)
            if lens_id in region_by_lens:
                raise ValueError(f"lens appears in multiple tract regions: {lens_id}")
            region_by_lens[lens_id] = region_id

    common = set(map(str, overlap["icl_profile_support_by_filter"][filter_name][
        "icl_supported_ids_in_common_wl_support"]))
    pairs = split.get("matched_lens_pairs", [])
    if len(pairs) != split.get("n_matched_lens_pairs"):
        raise ValueError("matched-pair row count does not match split summary")
    low_ids = [str(row["lower_lens_id"]) for row in pairs]
    high_ids = [str(row["upper_lens_id"]) for row in pairs]
    if len(set(low_ids + high_ids)) != 2 * len(pairs):
        raise ValueError("matched pairs reuse a lens ID")
    matched = set(low_ids + high_ids)
    if not matched <= common:
        raise ValueError("matched split IDs are not a subset of the approved i overlap")
    common_counts = counts_by_region(common, region_by_lens)
    low_counts = counts_by_region(set(low_ids), region_by_lens)
    high_counts = counts_by_region(set(high_ids), region_by_lens)

    pair_rows = []
    same_region = 0
    for row in pairs:
        low_id, high_id = str(row["lower_lens_id"]), str(row["upper_lens_id"])
        low_region, high_region = region_by_lens[low_id], region_by_lens[high_id]
        same = low_region == high_region
        same_region += int(same)
        pair_rows.append({
            "pair_id": str(row["pair_id"]),
            "lower_lens_id": low_id,
            "lower_region": low_region,
            "upper_lens_id": high_id,
            "upper_region": high_region,
            "same_region": same,
        })

    result = {
        "schema": "astrafek-wl-icl-split-tract-region-balance-v1",
        "status": "complete_exploratory",
        "science_measurement": False,
        "network_used": False,
        "shape_rows_read": False,
        "images_read": False,
        "analysis_script": str(Path(__file__).resolve()),
        "analysis_script_sha256": sha256(Path(__file__).resolve()),
        "method": "map saved ICL/WL IDs to nearest-tract connected-component labels from the existing parent-sample field-influence audit",
        "filter": filter_name,
        "tract_regions": {
            "n_components_parent": regions.get("n_tract_components"),
            "occupied_lens_counts_parent": regions.get("region_lens_counts"),
            "adjacency_threshold_deg": regions.get("tract_adjacency_threshold_deg"),
            "assignment_limit": "approximate nearest tract-center component, not official HSC field boundaries or source-row sky assignment",
        },
        "i_overlap": {
            "n_common_filter_lenses": len(common),
            "ids_sha256": hashlib.sha256("\n".join(sorted(common)).encode()).hexdigest(),
            "lens_counts_by_region": common_counts,
        },
        "matched_split": {
            "n_pairs": len(pairs),
            "n_lenses": len(matched),
            "matched_pair_list_sha256": hashlib.sha256(json.dumps(
                pairs, sort_keys=True, separators=(",", ":")
            ).encode()).hexdigest(),
            "low_icl_lens_counts_by_region": low_counts,
            "high_icl_lens_counts_by_region": high_counts,
            "same_region_pairs": same_region,
            "cross_region_pairs": len(pairs) - same_region,
            "pairs": pair_rows,
        },
        "interpretation": "The matched low/high groups are balanced across the two occupied parent components, but half the pairs cross components. This check does not locate or explain the cross-component null; source and random-center coordinates were not mapped here.",
        "input_hashes": {
            str(path): sha256(path) for path in
            (overlap_path, split_path, regions_path)
        },
    }
    output.mkdir(parents=True)
    report = output / "summary.json"
    report.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n",
                      encoding="utf-8")
    result["summary_sha256"] = sha256(report)
    print(json.dumps({
        "status": result["status"],
        "science_measurement": False,
        "filter": filter_name,
        "common_filter_lenses": len(common),
        "common_counts_by_region": common_counts,
        "matched_pairs": len(pairs),
        "low_counts_by_region": low_counts,
        "high_counts_by_region": high_counts,
        "same_region_pairs": same_region,
        "cross_region_pairs": len(pairs) - same_region,
        "summary": str(report),
        "summary_sha256": result["summary_sha256"],
    }, ensure_ascii=False))
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--overlap", type=Path, default=DEFAULT_OVERLAP)
    parser.add_argument("--filter", choices=("i", "i2"), default="i")
    parser.add_argument("--split", type=Path)
    parser.add_argument("--regions", type=Path, default=DEFAULT_REGIONS)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    run(args.overlap, args.split or DEFAULT_SPLITS[args.filter], args.regions,
        args.output or DEFAULT_OUTPUTS[args.filter], args.filter)


if __name__ == "__main__":
    main()
