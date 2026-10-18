#!/usr/bin/env python3
"""Validate per-annulus ICL support outputs and optional strict-sample parity."""
import argparse
import csv
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path):
    with Path(path).open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def require(condition, message):
    if not condition:
        raise ValueError(message)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--strict-pair-table", type=Path,
                        help="Optional existing 50-300 kpc all-annuli pair table")
    parser.add_argument("--program", type=Path,
                        help="Optional local copy of the measurement script to hash-check")
    args = parser.parse_args()
    output = args.output.resolve()
    summary_path = output / "summary.json"
    require(summary_path.is_file(), f"Missing summary: {summary_path}")
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    require(summary.get("schema") == "hsc-icl-per-annulus-support-v1",
            "Unexpected per-annulus ICL schema")
    require(summary.get("status") == "complete"
            and summary.get("science_gates", {}).get("science_measurement") is False,
            "Run is incomplete or science_measurement is not false")
    require(summary.get("n_target_profiles") == 800
            and summary.get("n_control_profiles") == 798,
            "Expected the existing 800-target/798-control sample")
    require(summary.get("coverage_floor") == 0.4 and summary.get("kron_factor") == 2.5,
            "Unexpected support threshold or mask")

    checked_hashes = {}
    for name, expected in summary.get("output_sha256", {}).items():
        path = Path(name)
        require(path.is_file() and sha256(path) == expected,
                f"Output checksum mismatch: {path}")
        checked_hashes[path.name] = expected
    required_outputs = {"annulus_pair_profiles.tsv", "annulus_stack_summary.tsv",
                        "pair_metrics_k2p5_c40_r50_300.tsv"}
    require(required_outputs.issubset(checked_hashes),
            "Output checksum manifest is incomplete")
    if args.program:
        require(sha256(args.program) == summary.get("program_sha256"),
                "Measurement program checksum mismatch")

    profiles = read_tsv(output / "annulus_pair_profiles.tsv")
    stack_rows = read_tsv(output / "annulus_stack_summary.tsv")
    metrics = read_tsv(output / "pair_metrics_k2p5_c40_r50_300.tsv")
    require(len(metrics) == summary.get("n_200_300_pairs"),
            "200-300 kpc pair count differs from summary")
    require(len({(row["id"], row["filter"], row["radial_bin_kpc"])
                 for row in profiles}) == len(profiles),
            "Duplicate cluster/filter/annulus rows")
    require(len({(row["id"], row["filter"]) for row in metrics}) == len(metrics),
            "Duplicate 200-300 kpc cluster/filter rows")

    grouped = defaultdict(list)
    floor = float(summary["coverage_floor"])
    for row in profiles:
        target = float(row["target_sb_njy_arcsec2"])
        control = float(row["control_sb_njy_arcsec2"])
        delta = float(row["delta_sb_njy_arcsec2"])
        target_coverage = float(row["target_coverage"])
        control_coverage = float(row["control_coverage"])
        common_coverage = float(row["common_coverage"])
        require(all(math.isfinite(value) for value in
                    (target, control, delta, target_coverage, control_coverage,
                     common_coverage)), "Non-finite retained pair value")
        require(target_coverage >= floor and control_coverage >= floor
                and math.isclose(common_coverage, min(target_coverage, control_coverage),
                                 rel_tol=0, abs_tol=1e-12),
                "Retained pair violates annulus coverage rule")
        require(math.isclose(target - control, delta, rel_tol=1e-10, abs_tol=1e-10),
                "Target-control difference is inconsistent")
        require(row["filter"] in ("i", "i2"), "Unexpected filter")
        grouped[(row["filter"], row["radial_bin_kpc"])].append(row)

    summary_by_key = {(row["filter"], row["radial_bin_kpc"]): row for row in stack_rows}
    require(len(summary_by_key) == 12, "Expected two filters across six radial annuli")
    for key, rows in grouped.items():
        require(key in summary_by_key, f"Missing support summary for {key}")
        record = summary_by_key[key]
        require(int(record["n_pairs"]) == len(rows), f"Pair count mismatch for {key}")
        mean = sum(float(row["delta_sb_njy_arcsec2"]) for row in rows) / len(rows)
        require(math.isclose(float(record["mean_pair_delta_njy_arcsec2"]), mean,
                             rel_tol=1e-10, abs_tol=1e-10), f"Mean mismatch for {key}")

    metric_by_key = {(row["id"], row["filter"]): row for row in metrics}
    profile_outer = {(row["id"], row["filter"]): row for row in profiles
                     if row["radial_bin_kpc"] == "200-300"}
    require(set(metric_by_key) == set(profile_outer),
            "Pair metrics table does not match 200-300 kpc profile support")
    for key, row in metric_by_key.items():
        source = profile_outer[key]
        require(math.isclose(float(row["delta_sb_200_300"]),
                             float(source["delta_sb_njy_arcsec2"]),
                             rel_tol=1e-10, abs_tol=1e-10),
                f"Pair metric delta mismatch for {key}")
        require(math.isclose(float(row["coverage_200_300"]),
                             float(source["common_coverage"]),
                             rel_tol=1e-10, abs_tol=1e-10),
                f"Pair metric coverage mismatch for {key}")

    strict_matches = 0
    if args.strict_pair_table:
        strict_rows = read_tsv(args.strict_pair_table)
        for strict in strict_rows:
            key = (strict["id"], strict["filter"])
            require(key in metric_by_key,
                    f"Strict all-annuli pair missing from annulus support: {key}")
            current = metric_by_key[key]
            require(math.isclose(float(strict["delta_sb_200_300"]),
                                 float(current["delta_sb_200_300"]),
                                 rel_tol=1e-10, abs_tol=1e-10),
                    f"Strict pair delta changed for {key}")
            strict_matches += 1

    result = {
        "status": "verified", "science_measurement": False,
        "target_profiles": summary["n_target_profiles"],
        "control_profiles": summary["n_control_profiles"],
        "annulus_pair_rows": len(profiles), "pair_metrics_200_300": len(metrics),
        "strict_pair_rows_verified": strict_matches,
        "verified_output_hashes": len(checked_hashes),
        "per_filter_bin_support": summary["per_filter_bin_support"],
    }
    print(json.dumps(result, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except ValueError as error:
        raise SystemExit(f"Validation failed: {error}")
