#!/usr/bin/env python3
"""Leave-one-cluster influence audit for the selected ICL--satellite relation."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

from audit_satellite_icl_region_fixed_effects import (
    DEFAULT_ASSOC, DEFAULT_BASELINE, DEFAULT_EXTRACTION, DEFAULT_ICL,
    partial_correlation,
)
from audit_satellite_icl_spatial_blocks import load_verified_rows


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def leave_one_cluster_out(rows):
    rows = list(rows)
    ids = [str(row["cluster_id"]) for row in rows]
    if len(rows) < 6 or len(ids) != len(set(ids)):
        raise ValueError("need at least six rows with unique cluster IDs")
    base = partial_correlation(rows)["partial_r"]
    records = []
    for removed_id in ids:
        retained = [row for row in rows if str(row["cluster_id"]) != removed_id]
        result = partial_correlation(retained)
        records.append({
            "removed_cluster_id": removed_id,
            "partial_r": result["partial_r"],
            "delta_r_from_full": result["partial_r"] - base,
            "n_remaining": result["n_rows"],
            "residual_degrees_of_freedom": result["residual_degrees_of_freedom"],
        })
    values = [item["partial_r"] for item in records]
    return {
        "n_clusters": len(rows),
        "full_sample_partial_r": base,
        "leave_one_out_min_r": min(values),
        "leave_one_out_max_r": max(values),
        "leave_one_out_median_r": sorted(values)[len(values) // 2],
        "n_sign_changes": sum(value * base < 0 for value in values),
        "largest_absolute_change": max(records, key=lambda item: abs(item["delta_r_from_full"])),
        "records": records,
    }


def analyze(association: Path, icl: Path, extraction: Path,
            baseline_path: Path, output_path: Path) -> dict:
    output_path = Path(output_path).resolve()
    if output_path.exists():
        raise FileExistsError(f"refusing to overwrite existing output: {output_path}")
    baseline_path = Path(baseline_path).resolve()
    baseline = json.loads(baseline_path.read_text(encoding="utf-8"))
    if baseline.get("schema") != "hsc-satellite-icl-spatial-block-bootstrap-v1" or baseline.get("science_measurement") is not False:
        raise ValueError("unexpected exploratory baseline")
    rows, provenance = load_verified_rows(association, icl, extraction)
    baseline_r = {
        (item["selection"], item["filter"]): float(
            item["partial_spatial_block_bootstrap"]["partial_r_point"])
        for item in baseline.get("results", [])
    }
    for key, value in provenance.items():
        if baseline.get("provenance", {}).get(key) != value:
            raise ValueError(f"baseline provenance mismatch for {key}")
    results = []
    sample_definitions = (
        ("all_supported_pairs", 68, 26,
         lambda row: True),
        ("center_eligible", 15, 8,
         lambda row: str(row["center_eligible"]).lower() == "true"),
    )
    for selection, expected_i2, expected_i, keep in sample_definitions:
        for band, expected_count in (("i2", expected_i2), ("i", expected_i)):
            selected = [row for row in rows
                        if row["filter"] == band and keep(row)]
            if len(selected) != expected_count:
                raise ValueError(f"fixed {selection}/{band} sample count changed")
            result = leave_one_cluster_out(selected)
            if not math.isclose(result["full_sample_partial_r"],
                                baseline_r[(selection, band)],
                                rel_tol=0, abs_tol=1e-12):
                raise ValueError(f"baseline partial correlation mismatch for {selection}/{band}")
            result["filter"] = band
            result["selection"] = selection
            result["baseline_reproduction"] = "exact_match_within_1e-12"
            results.append(result)
    output = {
        "schema": "astrafek-satellite-icl-cluster-influence-v1",
        "status": "complete_exploratory_influence_audit",
        "science_measurement": False,
        "covariance_estimated": False,
        "method": {
            "sample": "fixed validated all-supported and stricter center-eligible subsets; one cluster row removed at a time within each filter",
            "statistic": "partial Pearson correlation of 200-300 kpc ICL and 50-75 kpc satellite light, controlling redshift and richness",
            "interpretation": "single-cluster influence diagnostic, not a confidence interval or detection test",
        },
        "results": results,
        "provenance": provenance,
        "inputs": {
            "baseline": str(baseline_path),
            "baseline_sha256": sha256(baseline_path),
        },
        "limitations": [
            "the association sample is selected by prior support and extraction gates",
            "leave-one-cluster ranges do not quantify sampling uncertainty",
            "PSF wings, sky transfer, mask leakage, and measurement covariance remain unresolved",
            "no new HSC query or data transfer was performed",
        ],
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("x", encoding="utf-8") as stream:
        json.dump(output, stream, ensure_ascii=False, indent=2, allow_nan=False)
        stream.write("\n")
    return output


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--association", type=Path, default=DEFAULT_ASSOC)
    parser.add_argument("--icl", type=Path, default=DEFAULT_ICL)
    parser.add_argument("--extraction-summary", type=Path, default=DEFAULT_EXTRACTION)
    parser.add_argument("--baseline", type=Path, default=DEFAULT_BASELINE)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = analyze(args.association, args.icl, args.extraction_summary,
                     args.baseline, args.output)
    print(json.dumps({key: [{k: item[k] for k in (
        "selection", "filter", "n_clusters", "full_sample_partial_r", "leave_one_out_min_r",
        "leave_one_out_max_r", "n_sign_changes", "largest_absolute_change")}
        for item in result["results"]] if key == "results" else result[key]
        for key in ("status", "science_measurement", "results")},
        ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
