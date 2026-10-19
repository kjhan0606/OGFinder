#!/usr/bin/env python3
"""Audit where the 800-cluster ICL/satellite analysis loses sample support.

This is a provenance-checked descriptive audit of existing summaries and TSVs.
It submits no queries, reads no image pixels, and makes no science measurement.
"""
from __future__ import annotations

import argparse
import ast
import csv
import hashlib
import json
import math
from collections import Counter
from pathlib import Path


ROOT = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
DEFAULT_MASKS = ROOT / "prestack_cluster_masks_2026-10-16_01"
DEFAULT_SATELLITES = ROOT / "satellite_outer_profile_full_2026-10-18_03"
DEFAULT_PER_BIN = ROOT / "satellite_outer_profile_perbin_full_2026-10-18_01"
DEFAULT_ASSOCIATION = ROOT / "satellite_icl_spatial_blocks_2026-10-20_02"
DEFAULT_OUTPUT = ROOT / "icl_satellite_sample_retention_2026-10-19_01/summary.json"
RING_INDEX = 3  # 200--300 kpc in the saved six-annulus coverage vector.
FACTORS = ("2.5", "3.5", "5.0")
THRESHOLDS = (0.4, 0.5)
SATELLITE_ANNULI_KPC = ("5-10", "10-15", "15-20", "20-30", "30-40", "40-50", "50-75")
MIN_CONTROLS = 3


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_literal(value, field):
    try:
        return ast.literal_eval(value)
    except (ValueError, SyntaxError) as error:
        raise ValueError(f"invalid saved Python literal in {field}") from error


def parse_bool(value, field):
    value = str(value).strip().lower()
    if value in ("true", "1"):
        return True
    if value in ("false", "0"):
        return False
    raise ValueError(f"invalid boolean in {field}: {value!r}")


def summarize_mask_support(rows):
    if len(rows) != 800 or len({str(row["id"]) for row in rows}) != 800:
        raise ValueError("expected 800 unique prepared cluster/filter records")
    result = {}
    for band in ("i", "i2"):
        selected = [row for row in rows if row["primary_filter"] == band]
        metrics = {
            "n_frames": len(selected),
            "n_unique_center_candidates": 0,
            "n_center_unmasked_base": 0,
            "n_center_and_unique_unmasked_by_kron": {},
            "coverage_support": {},
        }
        for factor in FACTORS:
            metrics["n_center_and_unique_unmasked_by_kron"][factor] = 0
            metrics["coverage_support"][factor] = {
                str(t): {"center_eligible": 0, "all_six_annuli": 0} for t in THRESHOLDS
            }
        for row in selected:
            unique = parse_bool(row["auto_unique_center_candidate"], "auto_unique_center_candidate")
            base_unmasked = not parse_bool(row["center_base_masked"], "center_base_masked")
            metrics["n_unique_center_candidates"] += int(unique)
            metrics["n_center_unmasked_base"] += int(base_unmasked)
            centers = parse_literal(row["center_masked_by_kron"], "center_masked_by_kron")
            coverages = parse_literal(row["coverage_by_kron"], "coverage_by_kron")
            for factor in FACTORS:
                if factor not in centers or factor not in coverages or len(coverages[factor]) != 6:
                    raise ValueError(f"missing expected coverage vector for Kron {factor}")
                center_ok = unique and not bool(centers[factor])
                metrics["n_center_and_unique_unmasked_by_kron"][factor] += int(center_ok)
                values = [float(v) for v in coverages[factor]]
                if not all(math.isfinite(v) and 0 <= v <= 1 for v in values):
                    raise ValueError("coverage fractions must be finite values in [0, 1]")
                for threshold in THRESHOLDS:
                    entry = metrics["coverage_support"][factor][str(threshold)]
                    if center_ok and values[RING_INDEX] >= threshold:
                        entry["center_eligible"] += 1
                    if center_ok and all(value >= threshold for value in values):
                        entry["all_six_annuli"] += 1
        if not selected:
            raise ValueError(f"no records for filter {band}")
        result[band] = metrics
    if sum(v["n_frames"] for v in result.values()) != 800:
        raise ValueError("filter counts do not sum to the parent sample")
    return result


def summarize_candidate_diagnostics(rows):
    by_filter = {}
    for band in ("i", "i2"):
        selected = [row for row in rows if row["filter"] == band and row["cohort"] == "full_sample"]
        statuses = Counter(row["status"] for row in selected)
        eligible_target_coverage = statuses.get("insufficient_valid_controls", 0) + statuses.get("retained", 0)
        if sum(statuses.values()) != len(selected):
            raise ValueError("candidate diagnostic status tally mismatch")
        by_filter[band] = {
            "n_candidate_satellite_positions": len(selected),
            "status_counts": dict(sorted(statuses.items())),
            "n_target_coverage_passed": eligible_target_coverage,
            "n_target_coverage_and_control_passed": statuses.get("retained", 0),
            "fraction_retained_of_candidates": statuses.get("retained", 0) / len(selected) if selected else None,
        }
    return by_filter


def _satellite_key(row):
    return (str(row["cluster_id"]), str(row["filter"]), str(row["member_index"]),
            str(row["detection_number"]))


def summarize_control_gate_attrition(full_rows, per_bin_rows, per_bin_profiles):
    """Explain full-profile control losses and cross-check independent-bin recovery."""
    full_bad = {_satellite_key(row): row for row in full_rows
                if row["status"] == "insufficient_valid_controls"}
    if len(full_bad) != sum(row["status"] == "insufficient_valid_controls" for row in full_rows):
        raise ValueError("duplicate full-profile insufficient-control candidate")
    per_bin_by_key = {_satellite_key(row): row for row in per_bin_rows}
    if len(per_bin_by_key) != len(per_bin_rows) or not set(full_bad).issubset(per_bin_by_key):
        raise ValueError("full and per-annulus candidate diagnostics do not join one-to-one")

    outer_profile_keys = {_satellite_key(row) for row in per_bin_profiles
                          if row["radial_bin_kpc"] == SATELLITE_ANNULI_KPC[-1]}
    if len(outer_profile_keys) != sum(row["radial_bin_kpc"] == SATELLITE_ANNULI_KPC[-1]
                                      for row in per_bin_profiles):
        raise ValueError("duplicate outer-annulus satellite profile")

    result = {
        "minimum_valid_controls": MIN_CONTROLS,
        "annuli_kpc": list(SATELLITE_ANNULI_KPC),
        "n_full_profile_candidates_with_insufficient_controls": len(full_bad),
        "by_filter": {},
        "all_annuli_valid_control_count_distribution": {},
        "n_with_at_least_minimum_controls_by_annulus": [0] * len(SATELLITE_ANNULI_KPC),
        "outer_50_75_kpc_crosscheck": {},
        "interpretation": [
            "control positions are generated before per-annulus coverage cuts; most candidates existed but failed the complete-profile gate",
            "independent-annulus extraction can retain a satellite in a supported bin even when the complete seven-bin profile fails",
        ],
    }
    control_distribution = Counter()
    outer_eligible = set()
    for key, row in full_bad.items():
        per_bin = per_bin_by_key[key]
        if (row["valid_controls_by_bin"] != per_bin["valid_controls_by_bin"]
                or row["candidate_control_positions"] != per_bin["candidate_control_positions"]):
            raise ValueError(f"control geometry changed between extraction runs: {key}")
        counts = [int(value) for value in row["valid_controls_by_bin"].split(",")]
        target_coverage = [float(value) for value in row["target_coverage_by_bin"].split(",")]
        if len(counts) != len(SATELLITE_ANNULI_KPC) or len(target_coverage) != len(counts):
            raise ValueError("unexpected satellite radial-bin diagnostic width")
        valid_full = int(row["valid_control_profiles"])
        if valid_full not in range(MIN_CONTROLS):
            raise ValueError("candidate is not below the minimum complete-profile control count")
        candidate_n = int(row["candidate_control_positions"])
        if not 0 <= valid_full <= candidate_n <= 8:
            raise ValueError("invalid candidate/valid-control counts")
        control_distribution[str(valid_full)] += 1
        eligible_bins = [n >= MIN_CONTROLS for n in counts]
        for index, eligible in enumerate(eligible_bins):
            result["n_with_at_least_minimum_controls_by_annulus"][index] += int(eligible)
        if target_coverage[-1] >= 0.60 and eligible_bins[-1]:
            outer_eligible.add(key)

    for band in ("i", "i2"):
        selected = [(key, row) for key, row in full_bad.items() if row["filter"] == band]
        result["by_filter"][band] = {
            "n_candidates": len(selected),
            "candidate_control_position_distribution": dict(sorted(Counter(
                int(row["candidate_control_positions"]) for _, row in selected
            ).items())),
            "complete_profile_valid_control_distribution": dict(sorted(Counter(
                int(row["valid_control_profiles"]) for _, row in selected
            ).items())),
            "n_with_at_least_minimum_controls_by_annulus": [
                sum(int(row["valid_controls_by_bin"].split(",")[index]) >= MIN_CONTROLS
                    for _, row in selected)
                for index in range(len(SATELLITE_ANNULI_KPC))
            ],
        }
    result["all_annuli_valid_control_count_distribution"] = dict(sorted(control_distribution.items()))
    if outer_eligible != (outer_profile_keys & set(full_bad)):
        raise ValueError("independent 50-75 kpc extraction does not match eligible full-profile failures")
    result["outer_50_75_kpc_crosscheck"] = {
        "n_full_profile_failures_with_target_and_control_support": len(outer_eligible),
        "n_matching_saved_per_bin_profiles": len(outer_profile_keys & set(full_bad)),
        "matched_exactly": True,
    }
    return result


def audit(masks_dir, satellites_dir, per_bin_dir, association_dir):
    masks_dir, satellites_dir, per_bin_dir, association_dir = map(
        Path, (masks_dir, satellites_dir, per_bin_dir, association_dir)
    )
    mask_summary_path = masks_dir / "summary.json"
    records_path = masks_dir / "records.tsv"
    mask_summary = json.loads(mask_summary_path.read_text())
    if (mask_summary.get("status") != "complete" or mask_summary.get("n_expected") != 800
            or mask_summary.get("n_success") != 800 or mask_summary.get("n_error") != 0
            or mask_summary.get("science_measurement") is not False):
        raise ValueError("the 800-frame mask-preparation summary is not complete")
    expected_hash = mask_summary["output_sha256"][str(records_path)]
    if sha256(records_path) != expected_hash:
        raise ValueError("prepared mask records checksum mismatch")
    with records_path.open(encoding="utf-8", newline="") as stream:
        mask_rows = list(csv.DictReader(stream, delimiter="\t"))

    satellite_summary_path = satellites_dir / "summary.json"
    candidate_path = satellites_dir / "candidate_diagnostics.tsv"
    satellite_summary = json.loads(satellite_summary_path.read_text())
    if (satellite_summary.get("status") != "complete"
            or satellite_summary.get("sample_scope") != "full_manifest"
            or satellite_summary.get("n_frames_requested") != 800
            or satellite_summary.get("n_frames_completed") != 800
            or satellite_summary.get("n_frame_errors") != 0
            or satellite_summary.get("science_measurement") is not False):
        raise ValueError("authoritative full-sample satellite profile run is incomplete")
    candidate_hash = satellite_summary["output_sha256"][str(candidate_path)]
    if sha256(candidate_path) != candidate_hash:
        raise ValueError("candidate-diagnostics table checksum mismatch")
    with candidate_path.open(encoding="utf-8", newline="") as stream:
        candidate_rows = list(csv.DictReader(stream, delimiter="\t"))

    per_bin_summary_path = per_bin_dir / "summary.json"
    per_bin_summary = json.loads(per_bin_summary_path.read_text())
    if (per_bin_summary.get("status") != "complete"
            or per_bin_summary.get("n_frames_requested") != 800
            or per_bin_summary.get("n_frames_completed") != 800
            or per_bin_summary.get("n_frame_errors") != 0
            or per_bin_summary.get("science_measurement") is not False):
        raise ValueError("per-annulus satellite profile run is incomplete")
    per_bin_candidates_path = per_bin_dir / "candidate_diagnostics.tsv"
    per_bin_profiles_path = per_bin_dir / "satellite_profiles.tsv"
    for path in (per_bin_candidates_path, per_bin_profiles_path):
        expected = per_bin_summary.get("output_sha256", {}).get(str(path))
        if not expected or sha256(path) != expected:
            raise ValueError(f"per-annulus output checksum mismatch: {path.name}")
    with per_bin_candidates_path.open(encoding="utf-8", newline="") as stream:
        per_bin_candidate_rows = list(csv.DictReader(stream, delimiter="\t"))
    with per_bin_profiles_path.open(encoding="utf-8", newline="") as stream:
        per_bin_profile_rows = list(csv.DictReader(stream, delimiter="\t"))
    assoc_summary_path = association_dir / "summary.json"
    assoc = json.loads(assoc_summary_path.read_text())
    if assoc.get("status") != "complete_exploratory" or assoc.get("science_measurement") is not False:
        raise ValueError("spatial-block association summary is incomplete or has unexpected flags")

    per_bin_outer = {}
    for band in ("i", "i2"):
        values = per_bin_summary["filters"]["full_sample"][band]["bins"]["50-75"]
        per_bin_outer[band] = {key: values[key] for key in ("n_clusters", "n_satellites", "mean_delta", "bootstrap_low", "bootstrap_high")}
    matched = {row["filter"]: row["n_pairs"] for row in assoc["results"] if row["selection"] == "all_supported_pairs"}
    spatial = {row["filter"]: row["n_connected_5deg_blocks"] for row in assoc["results"] if row["selection"] == "all_supported_pairs"}
    if matched != {"i": 26, "i2": 68}:
        raise ValueError(f"unexpected matched ICL-satellite support: {matched}")
    if spatial != {"i": 13, "i2": 33}:
        raise ValueError(f"unexpected spatial block counts: {spatial}")
    mask_metrics = summarize_mask_support(mask_rows)
    candidate_metrics = summarize_candidate_diagnostics(candidate_rows)
    return {
        "schema": "hsc-icl-satellite-sample-retention-audit-v1",
        "status": "complete_descriptive",
        "science_measurement": False,
        "parent_sample_n_frames": 800,
        "prepared_mask_support_by_filter": mask_metrics,
        "strict_all_radial_bins_satellite_candidate_retention": candidate_metrics,
        "complete_profile_control_gate_attrition": summarize_control_gate_attrition(
            candidate_rows, per_bin_candidate_rows, per_bin_profile_rows),
        "per_annulus_satellite_50_75_kpc_stack": per_bin_outer,
        "matched_icl_satellite_association": {
            "n_cluster_filter_pairs_by_filter": matched,
            "n_connected_5deg_blocks_by_filter": spatial,
            "partial_r_i2": 0.2926385303703714,
            "spatial_block_ci95_i2": [0.04234471282373385, 0.4972908249681033],
            "ordinary_pair_ci95_i2": [-0.00709738518280587, 0.5260703386161849],
            "spearman_spatial_block_ci95_i2": [-0.08829931741011786, 0.30059075612694064],
        },
        "sample_size_interpretation": [
            "800 input frames are not 800 usable ICL-satellite pairs; the audited join retains 26 i and 68 i2 pairs.",
            "Per-annulus stacks recover more support than requiring all six satellite radial bins simultaneously.",
            "The i2 partial-Pearson interval is positive under the spatial-block bootstrap but includes zero under ordinary-pair and Spearman diagnostics; this is exploratory, not a confirmed association.",
            "The prepared support and candidate-control gates are diagnostic; no threshold was relaxed or sample reselected here.",
        ],
        "inputs": {
            "mask_summary_sha256": sha256(mask_summary_path),
            "mask_records_sha256": sha256(records_path),
            "satellite_summary_sha256": sha256(satellite_summary_path),
            "candidate_diagnostics_sha256": candidate_hash,
            "per_bin_summary_sha256": sha256(per_bin_summary_path),
            "association_summary_sha256": sha256(assoc_summary_path),
        },
        "limitations": [
            "the input 800 cluster frames are the fixed selected parent, not a complete HSC CAMIRA population",
            "cluster/filter rows and 5-degree blocks are not interchangeable independent samples",
            "PSF wings, sky transfer, extinction, correlated coadd noise, and member completeness remain uncalibrated",
            "weak-lensing random-lens subtraction and total spatial covariance remain pending",
        ],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--masks", type=Path, default=DEFAULT_MASKS)
    parser.add_argument("--satellites", type=Path, default=DEFAULT_SATELLITES)
    parser.add_argument("--per-bin", type=Path, default=DEFAULT_PER_BIN)
    parser.add_argument("--association", type=Path, default=DEFAULT_ASSOCIATION)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    output = args.output.resolve()
    if output.exists():
        parser.error(f"refusing to overwrite {output}")
    result = audit(args.masks, args.satellites, args.per_bin, args.association)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False) + "\n")
    print(json.dumps({
        "status": result["status"],
        "science_measurement": False,
        "parent_frames": result["parent_sample_n_frames"],
        "matched_pairs": result["matched_icl_satellite_association"]["n_cluster_filter_pairs_by_filter"],
        "output": str(output),
        "output_sha256": sha256(output),
    }, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
