#!/usr/bin/env python3
"""Validate ICL annulus mask/coverage sensitivity and reference-cut parity."""
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
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def read_tsv(path):
    with Path(path).open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def require(condition, message):
    if not condition:
        raise ValueError(message)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--sample", type=Path, required=True)
    parser.add_argument("--profiles", type=Path, required=True)
    parser.add_argument("--controls", type=Path, required=True)
    parser.add_argument("--per-annulus-reference", type=Path, required=True,
                        help="Verified K2.5/40% per-annulus ICL output")
    parser.add_argument("--program", type=Path)
    args = parser.parse_args()
    output, sample, profiles, controls, reference = (
        args.output.resolve(), args.sample.resolve(), args.profiles.resolve(),
        args.controls.resolve(), args.per_annulus_reference.resolve())
    summary_path = output / "summary.json"
    require(summary_path.is_file(), f"Missing summary: {summary_path}")
    summary = read_json(summary_path)
    require(summary.get("schema") == "hsc-icl-annulus-selection-sensitivity-v1"
            and summary.get("status") == "complete"
            and summary.get("science_measurement") is False,
            "Unexpected schema, incomplete run, or science_measurement is not false")
    require(summary.get("n_targets") == 800 and summary.get("n_matched_controls") == 798
            and summary.get("n_configurations") == 24,
            "Expected full 800-target/798-control sensitivity grid")

    target_table = profiles / "cluster_profiles.tsv"
    target_archive = profiles / "cluster_profiles.npz"
    profiles_summary_path = profiles / "summary.json"
    extraction_summary_path = controls / "extraction_summary.json"
    manifest_path = sample / "sample_manifest.json"
    for path in (manifest_path, target_table, target_archive, profiles_summary_path,
                 extraction_summary_path):
        require(summary.get("input_sha256", {}).get(str(path)) == sha256(path),
                f"Input checksum mismatch: {path}")

    profile_summary = read_json(profiles_summary_path)
    extraction_summary = read_json(extraction_summary_path)
    require(profile_summary.get("status") == "complete"
            and profile_summary.get("science_gates", {}).get("science_measurement") is False
            and extraction_summary.get("status") == "complete"
            and extraction_summary.get("science_measurement") is False,
            "Source profile/control products are incomplete or non-exploratory")

    control_hashes = []
    for record in extraction_summary.get("records", []):
        control_id = str(record["control_id"])
        baseline = controls / control_id / "baseline.npz"
        digest = sha256(baseline)
        require(record.get("baseline_sha256") == digest,
                f"Control baseline hash mismatch: {control_id}")
        require(profile_summary.get("source_sha256", {}).get(str(baseline)) == digest,
                f"Control baseline not linked to target stack: {control_id}")
        control_hashes.append((str(baseline), digest))
    aggregate = hashlib.sha256(json.dumps(
        sorted(control_hashes), separators=(",", ":")).encode()).hexdigest()
    require(summary.get("input_sha256", {}).get("control_baseline_hashes_sha256") == aggregate,
            "Aggregated control-baseline checksum mismatch")

    checked = {}
    for name, expected in summary.get("output_sha256", {}).items():
        path = Path(name)
        require(path.is_file() and sha256(path) == expected,
                f"Output checksum mismatch: {path}")
        checked[path.name] = expected
    require({"icl_annulus_sensitivity.tsv", "selected_pair_profiles.tsv"}.issubset(checked),
            "Output hash manifest is incomplete")
    if args.program:
        require(summary.get("program_sha256") == sha256(args.program),
                "Analysis program checksum mismatch")

    rows = read_tsv(output / "icl_annulus_sensitivity.tsv")
    details = read_tsv(output / "selected_pair_profiles.tsv")
    require(len(rows) == 24, "Expected 24 filter/factor/coverage configurations")
    keys = [(float(row["kron_factor"]), float(row["coverage_floor"]), row["filter"])
            for row in rows]
    require(len(set(keys)) == len(rows), "Duplicate sensitivity-grid summary row")
    grouped = defaultdict(list)
    for row in details:
        key = (float(row["kron_factor"]), float(row["coverage_floor"]), row["filter"])
        grouped[key].append(row)
        target, control = float(row["target_flux"]), float(row["control_flux"])
        delta = float(row["delta_flux"])
        target_cov = float(row["target_coverage"])
        control_cov = float(row["control_coverage"])
        require(all(math.isfinite(value) for value in
                    (target, control, delta, target_cov, control_cov)),
                "Non-finite selected pair")
        require(target_cov >= key[1] and control_cov >= key[1],
                f"Selected pair violates coverage threshold: {key}/{row['id']}")
        require(math.isclose(delta, target - control, rel_tol=1e-10, abs_tol=1e-10),
                f"Target-control flux mismatch: {key}/{row['id']}")
    require(len({(float(row["kron_factor"]), float(row["coverage_floor"]),
                  row["filter"], row["id"]) for row in details}) == len(details),
            "Duplicate selected pair profile")

    for row in rows:
        key = (float(row["kron_factor"]), float(row["coverage_floor"]), row["filter"])
        selected = grouped[key]
        require(int(row["n_pairs"]) == len(selected), f"Pair count mismatch: {key}")
        expected_mean = (sum(float(item["delta_flux"]) for item in selected) / len(selected)
                         if selected else None)
        if expected_mean is not None:
            require(math.isclose(float(row["mean_pair_delta_njy_arcsec2"]), expected_mean,
                                 rel_tol=1e-10, abs_tol=1e-10),
                    f"Mean mismatch: {key}")
            n_blocks = int(row["n_joint_spatial_blocks"])
            require(n_blocks >= 1, f"Selected pairs lack spatial-block labels: {key}")
            if n_blocks >= 2:
                require(math.isfinite(float(row["spatial_block_bootstrap_95_low"]))
                        and math.isfinite(float(row["spatial_block_bootstrap_95_high"])),
                        f"Invalid bootstrap interval: {key}")
            else:
                require(not row["spatial_block_bootstrap_95_low"]
                        and not row["spatial_block_bootstrap_95_high"],
                        f"One-block configuration must report an unavailable interval: {key}")

    reference_summary = read_json(reference / "summary.json")
    reference_table = reference / "annulus_pair_profiles.tsv"
    require(reference_summary.get("status") == "complete"
            and reference_summary.get("science_gates", {}).get("science_measurement") is False
            and reference_summary.get("output_sha256", {}).get(str(reference_table)) == sha256(reference_table),
            "Per-annulus reference product is incomplete or changed")
    ref_rows = [row for row in read_tsv(reference_table)
                if row["radial_bin_kpc"] == "200-300"]
    ref_by_key = {(row["id"], row["filter"]): row for row in ref_rows}
    baseline_rows = [row for row in details
                     if float(row["kron_factor"]) == 2.5
                     and float(row["coverage_floor"]) == 0.4]
    baseline_by_key = {(row["id"], row["filter"]): row for row in baseline_rows}
    require(set(ref_by_key) == set(baseline_by_key),
            "K2.5/40% sensitivity selection differs from per-annulus reference")
    for key, reference_row in ref_by_key.items():
        current = baseline_by_key[key]
        for field, ref_field in (("target_flux", "target_sb_njy_arcsec2"),
                                 ("control_flux", "control_sb_njy_arcsec2"),
                                 ("delta_flux", "delta_sb_njy_arcsec2"),
                                 ("target_coverage", "target_coverage"),
                                 ("control_coverage", "control_coverage")):
            require(math.isclose(float(current[field]), float(reference_row[ref_field]),
                                 rel_tol=1e-10, abs_tol=1e-10),
                    f"Reference parity failed for {key}/{field}")

    result = {
        "status": "verified", "science_measurement": False,
        "sensitivity_configurations": len(rows), "selected_pair_rows": len(details),
        "verified_output_hashes": len(checked),
        "reference_pairs_verified": len(ref_by_key),
        "reference_filter_counts": {
            band: sum(row["filter"] == band for row in baseline_rows) for band in ("i", "i2")
        },
    }
    print(json.dumps(result, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except ValueError as error:
        raise SystemExit(f"Validation failed: {error}")
