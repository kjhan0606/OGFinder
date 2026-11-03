#!/usr/bin/env python3
"""Independently audit ICL/satellite covariance pair selection and scope.

Recomputes the annulus-support selection from the saved ICL support and
satellite-profile tables, then compares exact target/filter keys with both
pixel-covariance products and the fixed ICL--satellite association sample.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def run(args):
    icl_rows = read_tsv(args.icl_support)
    satellite_rows = read_tsv(args.satellite_profiles)
    association_rows = read_tsv(args.association_rows)
    if not icl_rows or not satellite_rows or not association_rows:
        raise ValueError("support inputs must be non-empty")

    icl_keys = {
        (str(row["id"]), str(row["filter"]))
        for row in icl_rows if float(row["coverage_200_300"]) >= 0.40
    }
    eligible_profiles = [
        row for row in satellite_rows
        if row["cohort"] == "full_sample"
        and row["radial_bin_kpc"] == "50-75"
        and float(row["target_coverage"]) >= 0.60
        and float(row["control_coverage"]) >= 0.60
        and int(row["n_controls"]) >= 3
        and (str(row["cluster_id"]), str(row["filter"])) in icl_keys
    ]
    profile_keys = [
        (str(row["cluster_id"]), str(row["filter"]), str(row["member_index"]),
         str(row["detection_number"]))
        for row in eligible_profiles
    ]
    if len(set(profile_keys)) != len(profile_keys):
        raise ValueError("eligible satellite profile identities are duplicated")
    expected_pairs = {(row[0], row[1]) for row in profile_keys}
    association_keys = {(str(row["target_id"]), str(row["filter"])) for row in association_rows}
    if len(association_keys) != len(association_rows):
        raise ValueError("association target/filter keys are duplicated")

    target_summary_path = args.target_covariance / "summary.json"
    target_table_path = args.target_covariance / "pair_pixel_covariance.tsv"
    combined_summary_path = args.target_control_covariance / "summary.json"
    combined_table_path = args.target_control_covariance / "target_control_pixel_covariance.tsv"
    target_summary = json.loads(target_summary_path.read_text(encoding="utf-8"))
    combined_summary = json.loads(combined_summary_path.read_text(encoding="utf-8"))
    for label, summary, table, table_name in (
        ("target", target_summary, target_table_path, "pair_pixel_covariance.tsv"),
        ("target_control", combined_summary, combined_table_path, "target_control_pixel_covariance.tsv"),
    ):
        if summary.get("status") not in {"complete_exploratory", "complete_diagnostic"}:
            raise ValueError(f"{label} covariance product is not complete")
        if summary.get("science_measurement") is not False:
            raise ValueError(f"{label} covariance product has an unexpected science gate")
        if sha256(table) != summary.get("output_sha256", {}).get(str(table)):
            raise ValueError(f"{label} covariance table hash differs from its summary")

    target_inputs = target_summary.get("input_sha256", {})
    for source in (args.icl_support, args.satellite_profiles):
        expected_hash = target_inputs.get(str(source.resolve()))
        if expected_hash is None or expected_hash != sha256(source):
            raise ValueError(f"target covariance summary does not verify support input {source}")
    if target_summary.get("selection", {}).get("n_cluster_filter_pairs") != len(expected_pairs):
        raise ValueError("target covariance summary pair count differs from recomputed support")
    if target_summary.get("selection", {}).get("n_satellite_profiles") != len(eligible_profiles):
        raise ValueError("target covariance summary profile count differs from recomputed support")

    target_rows = read_tsv(target_table_path)
    target_output_keys = {(str(row["cluster_id"]), str(row["filter"])) for row in target_rows}
    combined_rows = read_tsv(combined_table_path)
    combined_output_keys = {(str(row["cluster_id"]), str(row["filter"])) for row in combined_rows}
    if len(target_output_keys) != len(target_rows) or len(combined_output_keys) != len(combined_rows):
        raise ValueError("covariance output pair keys are duplicated")
    if target_output_keys != expected_pairs or combined_output_keys != expected_pairs:
        raise ValueError("covariance outputs do not exactly match independently selected support")
    missing_association = sorted(association_keys - combined_output_keys)
    if missing_association:
        raise ValueError(f"association pairs absent from covariance support: {missing_association}")
    if combined_summary.get("selection", {}).get("n_completed_pairs") != len(expected_pairs):
        raise ValueError("target/control summary completed-pair count differs from support")
    if combined_summary.get("selection", {}).get("n_errors") != 0:
        raise ValueError("target/control covariance product contains pair errors")

    return {
        "schema": "astrafek-icl-satellite-covariance-support-scope-verification-v1",
        "status": "verified",
        "science_measurement": False,
        "icl_support_sha256": sha256(args.icl_support),
        "satellite_profiles_sha256": sha256(args.satellite_profiles),
        "association_rows_sha256": sha256(args.association_rows),
        "target_covariance_summary_sha256": sha256(target_summary_path),
        "target_control_covariance_summary_sha256": sha256(combined_summary_path),
        "n_icl_support_keys_at_coverage_floor": len(icl_keys),
        "n_satellite_profiles_selected": len(eligible_profiles),
        "n_covariance_pairs_recomputed": len(expected_pairs),
        "n_covariance_output_pairs_verified": len(combined_output_keys),
        "n_association_pairs": len(association_keys),
        "n_association_pairs_covered": len(association_keys & combined_output_keys),
        "missing_association_pairs": missing_association,
        "n_supported_pairs_outside_association_sample": len(combined_output_keys - association_keys),
        "pairs_by_filter": {
            band: sum(key[1] == band for key in expected_pairs) for band in ("i", "i2")
        },
        "association_pairs_by_filter": {
            band: sum(key[1] == band for key in association_keys) for band in ("i", "i2")
        },
        "selection": {
            "icl_annulus_kpc": [200, 300],
            "icl_coverage_floor": 0.40,
            "satellite_annulus_kpc": [50, 75],
            "target_and_control_satellite_coverage_floor": 0.60,
            "minimum_azimuthal_controls": 3,
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--icl-support", type=Path, required=True)
    parser.add_argument("--satellite-profiles", type=Path, required=True)
    parser.add_argument("--association-rows", type=Path, required=True)
    parser.add_argument("--target-covariance", type=Path, required=True)
    parser.add_argument("--target-control-covariance", type=Path, required=True)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = run(args)
    if args.output:
        if args.output.exists():
            raise FileExistsError(f"refusing to overwrite {args.output}")
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
