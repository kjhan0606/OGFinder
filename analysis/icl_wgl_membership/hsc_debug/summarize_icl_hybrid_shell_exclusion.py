#!/usr/bin/env python3
"""Independently audit and summarize hybrid-PSF shell-exclusion batches."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np


DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
ASSOCIATION_SHA256 = "483f81540543fc46ec1cf4333bbd3e238799b493e48fe03235da3038f5714917"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def quantiles(values: list[float]) -> dict:
    array = np.asarray(values, dtype=float)
    if not len(array) or not np.isfinite(array).all():
        raise ValueError("quantile inputs must be nonempty and finite")
    q16, median, q84 = np.quantile(array, [0.16, 0.5, 0.84])
    return {"q16": float(q16), "median": float(median), "q84": float(q84),
            "min": float(array.min()), "max": float(array.max())}


def audit_filter(summary_path: Path, expected_filter: str,
                 association_keys: set[tuple[str, str]]) -> dict:
    summary_path = summary_path.resolve()
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    if (summary.get("status") != "complete" or summary.get("n_errors") != 0
            or summary.get("science_measurement") is not False
            or summary.get("observational_validation_passed") is not False
            or summary.get("flux_covariance_available") is not False
            or summary.get("exclude_supported_satellite_shells_from_icl") is not True):
        raise ValueError(f"{summary_path}: batch status or guard flags are unexpected")
    records = summary.get("records", [])
    if (summary.get("frame_filter") != expected_filter
            or len(records) != {"i": 12, "i2": 36}[expected_filter]
            or len(records) != summary.get("n_selected")
            or len(records) != summary.get("n_completed")):
        raise ValueError(f"{summary_path}: filter or record counts are inconsistent")
    pair_keys = [(str(row["cluster_id"]), str(row["filter"]))
                 for row in summary.get("selected_pairs", [])]
    record_keys = [(str(row["cluster_id"]), str(row["filter"])) for row in records]
    if (pair_keys != record_keys or len(set(record_keys)) != len(record_keys)
            or any(key not in association_keys for key in record_keys)):
        raise ValueError(f"{summary_path}: selected frame keys do not match the fixed association table")

    response_names = ("icl_response_to_diffuse_unit_injection",
                      "icl_response_to_all_satellite_envelopes")
    raw, excluded, paired = {key: [] for key in response_names}, {key: [] for key in response_names}, {key: [] for key in response_names}
    point_raw, point_excluded, coverage_raw, coverage_excluded = [], [], [], []
    shell_count = overlap_pixels = 0
    frame_hashes = []
    for record in records:
        frame_path = Path(record["path"]).resolve()
        actual_hash = sha256(frame_path)
        if actual_hash != record.get("sha256"):
            raise ValueError(f"{frame_path}: frame JSON checksum mismatch")
        frame = json.loads(frame_path.read_text(encoding="utf-8"))
        if (frame.get("science_measurement") is not False
                or frame.get("observational_validation_passed") is not False
                or frame.get("flux_covariance_available") is not False
                or frame.get("hybrid_psf", {}).get("psf_approved") is not False
                or frame.get("icl_satellite_shell_exclusion", {}).get("applied") is not True):
            raise ValueError(f"{frame_path}: frame guard flags are unexpected")
        for key in response_names:
            a = float(frame["fixed_mask_transfer_response"][key])
            b = float(frame["satellite_shell_excluded_transfer_response"][key])
            raw[key].append(a)
            excluded[key].append(b)
            paired[key].append(b - a)
        annulus = "200_300_kpc"
        p0 = frame["point_source_wing_profile"]["annuli"][annulus]
        p1 = frame["point_source_wing_profile_satellite_shells_excluded"]["annuli"][annulus]
        point_raw.append(float(p0["masked_annulus_flux_fraction"]))
        point_excluded.append(float(p1["masked_annulus_flux_fraction"]))
        coverage_raw.append(float(p0["mask_coverage"]))
        coverage_excluded.append(float(p1["mask_coverage"]))
        shell_count += len(frame["icl_satellite_shell_exclusion"]["supported_profiles"])
        overlap_pixels += int(frame["icl_satellite_shell_exclusion"]["overlap_pixels_removed"])
        frame_hashes.append({"cluster_id": frame["cluster_id"],
                             "filter": frame["filter"], "sha256": actual_hash})

    return {
        "filter": expected_filter,
        "summary_path": str(summary_path),
        "summary_sha256": sha256(summary_path),
        "status": summary["status"],
        "n_selected": len(records),
        "n_completed": summary["n_completed"],
        "n_errors": summary["n_errors"],
        "hybrid_parameters": summary["hybrid_parameters"],
        "psf_model_mode": summary["psf_model_mode"],
        "public_reference_transfer_mode": summary["public_reference_transfer_mode"],
        "science_measurement": summary["science_measurement"],
        "psf_approved": summary["reference"]["psf_approved"],
        "responses": {
            key: {"fixed_mask": quantiles(raw[key]),
                  "shell_excluded": quantiles(excluded[key]),
                  "paired_shell_excluded_minus_fixed": quantiles(paired[key])}
            for key in response_names
        },
        "point_source_unmasked_annulus_flux_fraction_200_300_kpc": {
            "fixed_mask": quantiles(point_raw),
            "shell_excluded": quantiles(point_excluded),
            "paired_shell_excluded_minus_fixed": quantiles(
                [b - a for a, b in zip(point_raw, point_excluded)])
        },
        "mask_coverage_200_300_kpc": {
            "fixed_mask": quantiles(coverage_raw),
            "shell_excluded": quantiles(coverage_excluded),
            "paired_shell_excluded_minus_fixed": quantiles(
                [b - a for a, b in zip(coverage_raw, coverage_excluded)])
        },
        "supported_shell_count": shell_count,
        "overlap_pixels_removed": overlap_pixels,
        "frame_hashes": frame_hashes,
        "limitations": summary["limitations"],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--association-table", type=Path, required=True)
    parser.add_argument("--i-summary", type=Path, required=True)
    parser.add_argument("--i2-summary", type=Path, required=True)
    args = parser.parse_args()
    if sha256(args.association_table) != ASSOCIATION_SHA256:
        raise ValueError("association table checksum differs from the approved fixed sample")
    with args.association_table.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    if len(rows) != 74:
        raise ValueError(f"expected the audited 74-row association table, got {len(rows)}")
    association_keys = {(str(row["target_id"]), str(row["filter"])) for row in rows}
    results = [audit_filter(args.i_summary, "i", association_keys),
               audit_filter(args.i2_summary, "i2", association_keys)]
    combined = {(row["cluster_id"], row["filter"])
                for result in results for row in result["frame_hashes"]}
    if (len(combined) != 48
            or results[0]["n_completed"] != 12
            or results[1]["n_completed"] != 36):
        raise ValueError(f"expected 48 unique association frames across both filters, got {len(combined)}")
    print(json.dumps({
        "status": "complete",
        "association_table": str(args.association_table.resolve()),
        "association_sha256": ASSOCIATION_SHA256,
        "association_rows": len(rows),
        "n_association_psf_supported_frames": len(combined),
        "filters": results,
        "science_measurement": False,
        "observational_validation_passed": False,
        "flux_covariance_available": False,
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
