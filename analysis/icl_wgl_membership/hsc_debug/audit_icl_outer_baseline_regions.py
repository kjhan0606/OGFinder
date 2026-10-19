#!/usr/bin/env python3
"""Test whether the ICL outer-annulus baseline is concentrated in one HSC region."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

from run_wl_hsc_tract_region_influence import (
    nearest_region_labels,
    tract_centers_from_patch_qa,
    tract_region_labels,
)


ROOT = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
DEFAULT_OUTER = ROOT / "icl_outer_annulus_baseline_2026-10-20_02"
DEFAULT_STACK = ROOT / "icl_2d_matched_stack_2026-10-19_02"
DEFAULT_PATCH_QA = ROOT / "catalogues_2026-10-13_01/patch_qa.csv"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def leave_one_region_out(values: list[float], regions: list[int]) -> dict:
    if len(values) != len(regions) or len(values) < 2:
        raise ValueError("matching values and at least two rows are required")
    if not all(math.isfinite(float(value)) for value in values):
        raise ValueError("regional values must be finite")
    unique = sorted(set(regions))
    if len(unique) < 2:
        raise ValueError("at least two occupied regions are required")
    full = sum(values) / len(values)
    per_region = {}
    deletions = {}
    for region in unique:
        subset = [value for value, label in zip(values, regions) if label == region]
        keep = [value for value, label in zip(values, regions) if label != region]
        per_region[str(region)] = {"n_pairs": len(subset), "mean": sum(subset) / len(subset)}
        deletions[str(region)] = {"n_pairs": len(keep), "mean": sum(keep) / len(keep)}
    return {"full_mean": full, "per_region": per_region, "leave_one_region_out": deletions}


def analyze(outer_dir: Path, stack_dir: Path, patch_qa: Path, output: Path,
            adjacency_deg: float = 2.0) -> dict:
    if output.exists():
        raise FileExistsError(f"refusing to overwrite {output}")
    summary_path = outer_dir / "summary.json"
    outer_table = outer_dir / "paired_outer_subtracted_profiles.tsv"
    stack_summary_path = stack_dir / "summary.json"
    stack_table = stack_dir / "paired_radial_profiles.tsv"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    stack_summary = json.loads(stack_summary_path.read_text(encoding="utf-8"))
    if summary.get("status") != "complete_exploratory_sensitivity" or summary.get("science_measurement") is not False:
        raise ValueError("outer-annulus input must be complete and exploratory")
    if (summary.get("output_sha256") != sha256(outer_table)
            or summary.get("inputs", {}).get("paired_profiles_sha256") != sha256(stack_table)):
        raise ValueError("outer or source-profile checksum mismatch")
    if stack_summary.get("science_measurement") is not False:
        raise ValueError("source stack must remain exploratory")

    source_rows = read_tsv(stack_table)
    source_by_key = {}
    for row in source_rows:
        key = (str(row["target_id"]), str(row["filter"]))
        if key in source_by_key:
            raise ValueError(f"duplicate source pair {key}")
        source_by_key[key] = row
    keys, centers = tract_centers_from_patch_qa(patch_qa)
    components = tract_region_labels(centers, adjacency_deg)
    if len(set(map(int, components))) < 2:
        raise ValueError("patch-QA tract graph has fewer than two connected regions")

    output_rows = read_tsv(outer_table)
    grouped: dict[tuple[str, str], list[dict]] = {}
    for row in output_rows:
        key = (str(row["target_id"]), str(row["filter"]))
        source = source_by_key.get(key)
        if source is None:
            raise ValueError(f"outer-baseline row missing from source stack: {key}")
        if str(source["control_id"]) != str(row["control_id"]):
            raise ValueError(f"matched control changed for {key}")
        grouped.setdefault((key[1], f"{float(row['r_min_kpc']):g}-{float(row['r_max_kpc']):g}"), []).append({
            "target_id": key[0],
            "ra_deg": float(source["ra_deg"]),
            "dec_deg": float(source["dec_deg"]),
            "raw": float(row["inner_delta_njy_arcsec2"]),
            "outer": float(row["outer_400_500_delta_njy_arcsec2"]),
            "corrected": float(row["outer_subtracted_delta_njy_arcsec2"]),
        })

    metrics = {}
    for (band, annulus), rows in sorted(grouped.items()):
        labels = nearest_region_labels([row["ra_deg"] for row in rows],
                                       [row["dec_deg"] for row in rows], centers, components)
        labels = [int(label) for label in labels]
        if len(set(labels)) < 2:
            raise ValueError(f"{band} {annulus} kpc sample occupies fewer than two regions")
        metrics.setdefault(band, {})[annulus] = {
            "n_pairs": len(rows),
            "n_occupied_regions": len(set(labels)),
            "raw_profile": leave_one_region_out([r["raw"] for r in rows], labels),
            "outer_baseline": leave_one_region_out([r["outer"] for r in rows], labels),
            "outer_subtracted_profile": leave_one_region_out([r["corrected"] for r in rows], labels),
        }

    result = {
        "schema": "astrafek-icl-outer-baseline-region-influence-v1",
        "status": "complete_exploratory_influence",
        "science_measurement": False,
        "covariance_estimated": False,
        "method": "equal-pair means and leave-one-approximate-tract-region-out influence for the same 5%-support-selected rows",
        "adjacency_deg": adjacency_deg,
        "n_patch_qa_tracts": len(keys),
        "n_tract_components": len(set(map(int, components))),
        "by_filter": metrics,
        "inputs": {
            "outer_summary": str(summary_path.resolve()),
            "outer_summary_sha256": sha256(summary_path),
            "outer_rows": str(outer_table.resolve()),
            "outer_rows_sha256": sha256(outer_table),
            "stack_summary": str(stack_summary_path.resolve()),
            "stack_summary_sha256": sha256(stack_summary_path),
            "stack_profiles": str(stack_table.resolve()),
            "stack_profiles_sha256": sha256(stack_table),
            "patch_qa": str(patch_qa.resolve()),
            "patch_qa_sha256": sha256(patch_qa),
        },
        "limitations": [
            "tract-center regions are approximate groups, not official HSC footprint polygons",
            "three occupied regions permit influence checks only, not a stable regional covariance or confidence interval",
            "the per-annulus 5% support threshold yields different selected rows across annuli",
            "outer-annulus subtraction can remove real extended ICL and this analysis cannot identify the physical origin of an offset",
            "no HSC query or download was made; science_measurement=false",
        ],
    }
    result["script_sha256"] = sha256(Path(__file__).resolve())
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outer-dir", type=Path, default=DEFAULT_OUTER)
    parser.add_argument("--stack-dir", type=Path, default=DEFAULT_STACK)
    parser.add_argument("--patch-qa", type=Path, default=DEFAULT_PATCH_QA)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--adjacency-deg", type=float, default=2.0)
    args = parser.parse_args()
    result = analyze(args.outer_dir.resolve(), args.stack_dir.resolve(),
                     args.patch_qa.resolve(), args.output.resolve(), args.adjacency_deg)
    print(json.dumps({"status": result["status"], "science_measurement": result["science_measurement"],
                      "n_tract_components": result["n_tract_components"],
                      "by_filter": result["by_filter"], "output": str(args.output.resolve())},
                     ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
