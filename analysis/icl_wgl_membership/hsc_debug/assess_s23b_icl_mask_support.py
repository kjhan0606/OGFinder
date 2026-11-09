#!/usr/bin/env python3
"""Summarize the engineering mask and radial-area support for S23B targets."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np


DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
TARGET_PLAN = DATA / "s23b_icl_target_image_plan_2026-11-07_01"
TARGET_WORK = TARGET_PLAN / "target_download_01"
TARGET_AUDIT = DATA / "s23b_icl_target_image_audit_2026-11-07_01"
MASK_WORK = DATA / "s23b_icl_mask_support_2026-11-07_01/target"
ZERO_SUPPORT = DATA / "s23b_icl_zero_support_audit_2026-11-07_02/summary.json"
DEFAULT_OUTPUT = DATA / "s23b_icl_mask_support_assessment_2026-11-07_01"
EXPECTED_TARGET_MANIFEST_SHA256 = "e487582ba40cd3e58da67704242db4d7219e381fb08f4b33fed47ce4ff268b44"
EXPECTED_TARGET_AUDIT_SHA256 = "8580d64ba76e1c1f64331f7ce29b8d6f9dde4b56914b699e9e09e57f339e4b25"
EXPECTED_ZERO_IDS = {"S23B_10647", "S23B_11786", "S23B_432", "S23B_8226"}
ANNULI = ((50.0, 100.0), (100.0, 150.0), (150.0, 200.0),
          (200.0, 300.0), (300.0, 400.0), (400.0, 500.0))


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def summarize(values: list[float]) -> dict:
    x = np.asarray(values, dtype=float)
    if not len(x) or not np.isfinite(x).all() or np.any((x < 0) | (x > 1)):
        raise ValueError("coverage values must be finite fractions in [0, 1]")
    return {
        "n_frames": int(len(x)),
        "median": float(np.median(x)),
        "p10": float(np.percentile(x, 10)),
        "p25": float(np.percentile(x, 25)),
        "p75": float(np.percentile(x, 75)),
        "p90": float(np.percentile(x, 90)),
        "n_zero": int(np.count_nonzero(x == 0)),
        "n_ge_0p5": int(np.count_nonzero(x >= 0.5)),
        "n_ge_0p8": int(np.count_nonzero(x >= 0.8)),
    }


def assess(output_dir: Path) -> dict:
    output_dir = output_dir.resolve()
    if output_dir.exists():
        raise FileExistsError(f"refusing to overwrite {output_dir}")
    manifest_path = TARGET_PLAN / "target_manifest.json"
    audit_path = TARGET_AUDIT / "summary.json"
    support_summary_path = MASK_WORK / "summary.json"
    records_path = MASK_WORK / "records.tsv"
    coverage_path = MASK_WORK / "coverage.tsv"
    for path in (manifest_path, audit_path, support_summary_path, records_path,
                 coverage_path, ZERO_SUPPORT):
        if not path.is_file():
            raise FileNotFoundError(path)

    manifest_sha = sha256(manifest_path)
    audit_sha = sha256(audit_path)
    if manifest_sha != EXPECTED_TARGET_MANIFEST_SHA256 or audit_sha != EXPECTED_TARGET_AUDIT_SHA256:
        raise ValueError("target manifest or independent image audit differs from the pinned input")
    manifest, image_audit = read_json(manifest_path), read_json(audit_path)
    support, zero = read_json(support_summary_path), read_json(ZERO_SUPPORT)
    if (image_audit.get("status") != "complete_header_wcs_file_integrity_audit"
            or image_audit.get("inputs_sha256", {}).get(str(manifest_path)) != manifest_sha
            or image_audit.get("n_expected_frames") != 1043
            or image_audit.get("n_verified_frames") != 1043
            or image_audit.get("science_measurement") is not False
            or support.get("n_expected") != 1043
            or support.get("n_success") != 1039
            or support.get("n_error") != 4
            or support.get("science_measurement") is not False
            or support.get("mask_ready_for_science") is not False
            or zero.get("status") != "complete_zero_pixel_support_audit"):
        raise ValueError("image audit, mask-support preparation, or zero-support audit is incomplete")
    if len({str(x["id"]) for x in manifest.get("clusters", [])}) != 1043:
        raise ValueError("target manifest does not contain 1,043 unique IDs")
    zero_ids = {str(x["id"]) for x in zero.get("frames", [])}
    error_ids = {str(x["id"]) for x in support.get("errors", [])}
    if zero_ids != EXPECTED_ZERO_IDS or error_ids != EXPECTED_ZERO_IDS:
        raise ValueError("mask preparation errors do not match the audited zero-support frames")
    if any(int(x.get("n_unflagged_finite_positive_variance", -1)) != 0 for x in zero["frames"]):
        raise ValueError("one or more rejected frames has unflagged pixel support")

    records = read_tsv(records_path)
    coverage_rows = read_tsv(coverage_path)
    record_by_id = {row["id"]: row for row in records}
    expected_ids = {str(row["id"]) for row in manifest["clusters"]} - EXPECTED_ZERO_IDS
    if len(records) != 1039 or set(record_by_id) != expected_ids:
        raise ValueError("mask-support record identities do not match the 1,039 usable-frame subset")
    if any(row.get("science_measurement") != "False"
           or row.get("mask_ready_for_science") != "False" for row in records):
        raise ValueError("mask-support records contain an invalid science-readiness flag")
    if any(int(row.get("n_members", -1)) != 0 for row in records):
        raise ValueError("target manifest unexpectedly contains member positions")
    if len(coverage_rows) != len(records) * 3 * len(ANNULI):
        raise ValueError("radial-coverage rows are incomplete")
    if {row["id"] for row in coverage_rows} != expected_ids:
        raise ValueError("radial-coverage IDs differ from mask records")

    group_specs = {"all": expected_ids, "i": set(), "i2": set()}
    for ident, row in record_by_id.items():
        band = row["primary_filter"].lower()
        if band not in {"i", "i2"}:
            raise ValueError(f"unexpected HSC filter family: {band}")
        group_specs[band].add(ident)
    coverage = {}
    for factor in ("2.5", "3.5", "5.0"):
        coverage[factor] = {}
        for group, selected_ids in group_specs.items():
            coverage[factor][group] = {}
            for lo, hi in ANNULI:
                values = [float(row["coverage"]) for row in coverage_rows
                          if row["id"] in selected_ids and row["kron_factor"] == factor
                          and float(row["r_min_kpc"]) == lo and float(row["r_max_kpc"]) == hi]
                coverage[factor][group][f"{int(lo)}-{int(hi)}"] = summarize(values)

    filter_counts = {band: len(ids) for band, ids in group_specs.items() if band != "all"}
    if filter_counts != {"i": 251, "i2": 788}:
        raise ValueError(f"unexpected post-exclusion filter counts: {filter_counts}")
    result = {
        "schema": "astrafek-s23b-icl-mask-support-assessment-v1",
        "status": "support_assessment_complete_science_gate_closed",
        "n_expected_targets": 1043, "n_mask_support_success": 1039,
        "n_zero_pixel_exclusions": 4, "zero_pixel_exclusion_ids": sorted(EXPECTED_ZERO_IDS),
        "post_exclusion_filter_counts": filter_counts,
        "auto_unique_center_candidates": int(support["n_auto_unique_center_candidates"]),
        "centers_with_base_quality_mask": int(support["n_center_base_masked"]),
        "centers_masked_by_kron": support["n_center_masked_by_kron"],
        "target_member_positions_in_manifest": 0,
        "radial_unmasked_area_fraction": coverage,
        "science_measurement": False, "images_ready_for_stacking": False,
        "mask_ready_for_science": False,
        "interpretation_limits": [
            "Coverage is an engineering diagnostic from automated source masks; no mask threshold or pair-selection cut is approved here.",
            "The centroid source was not manually confirmed as the BCG; unique central detections are candidates only.",
            "No CAMIRA member positions were present, so these masks do not include a member-catalog completeness pass.",
            "Target-control common support, sky residuals, extinction, PSF wings, coadd covariance, and image-derived measurement uncertainties remain unevaluated.",
            "Do not pool the i and i2 filters; keep filter-specific stacks and calibrations separate.",
        ],
        "inputs_sha256": {
            str(manifest_path): manifest_sha, str(audit_path): audit_sha,
            str(support_summary_path): sha256(support_summary_path),
            str(records_path): sha256(records_path), str(coverage_path): sha256(coverage_path),
            str(ZERO_SUPPORT): sha256(ZERO_SUPPORT),
        },
    }
    output_dir.mkdir(parents=True, exist_ok=False)
    summary_path = output_dir / "summary.json"
    summary_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": result["status"], "n_mask_support_success": result["n_mask_support_success"],
        "n_zero_pixel_exclusions": result["n_zero_pixel_exclusions"],
        "post_exclusion_filter_counts": filter_counts,
        "k2p5_all": coverage["2.5"]["all"],
        "science_measurement": False, "images_ready_for_stacking": False,
    }, indent=2))
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    assess(args.output_dir)


if __name__ == "__main__":
    main()
