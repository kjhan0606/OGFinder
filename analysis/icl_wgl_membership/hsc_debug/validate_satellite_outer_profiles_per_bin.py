#!/usr/bin/env python3
"""Verify per-annulus satellite-profile provenance, support and table invariants."""
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
    parser.add_argument("--expected-frames", type=int)
    args = parser.parse_args()
    output = args.output.resolve()
    summary_path = output / "summary.json"
    require(summary_path.is_file(), f"Missing summary: {summary_path}")
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    require(summary.get("schema") == "hsc-satellite-outer-profile-per-annulus-v1",
            "Unexpected per-annulus profile schema")
    require(summary.get("status") == "complete"
            and summary.get("science_measurement") is False,
            "Run is incomplete or science_measurement is not false")

    n_frames = int(summary.get("n_frames_requested", -1))
    frames = summary.get("frame_records", [])
    require(n_frames > 0 and summary.get("n_frames_completed") == n_frames
            and summary.get("n_frame_errors") == 0 and not summary.get("errors")
            and len(frames) == n_frames,
            "Frame provenance is incomplete or inconsistent")
    if args.expected_frames is not None:
        require(n_frames == args.expected_frames,
                f"Expected {args.expected_frames} frames, found {n_frames}")
    frame_ids = [str(row.get("cluster_id")) for row in frames]
    require(len(set(frame_ids)) == n_frames, "Frame provenance has duplicate cluster IDs")
    require(all(row.get("image_sha256") and row.get("mask_sha256") for row in frames),
            "At least one frame lacks image/mask checksums")

    audit_path = summary.get("artifact_audit")
    audit_hash = summary.get("artifact_audit_sha256")
    if audit_hash:
        audit_path = Path(audit_path)
        require(audit_path.is_file() and sha256(audit_path) == audit_hash,
                "Preparation artifact-audit checksum mismatch")
        audit = json.loads(audit_path.read_text(encoding="utf-8"))
        require(audit.get("passed") is True and audit.get("science_measurement") is False,
                "Preparation audit is not a passing exploratory audit")

    checked_hashes = {}
    for name, expected in summary.get("output_sha256", {}).items():
        path = Path(name)
        require(path.is_file() and sha256(path) == expected,
                f"Output checksum mismatch: {path}")
        checked_hashes[path.name] = expected
    required = {"satellite_profiles.tsv", "cluster_profiles.tsv", "filter_profiles.tsv",
                "candidate_diagnostics.tsv", "satellite_profile_review.png"}
    require(required.issubset(checked_hashes), "Output checksum manifest is incomplete")

    profiles = read_tsv(output / "satellite_profiles.tsv")
    cluster_rows = read_tsv(output / "cluster_profiles.tsv")
    filter_rows = read_tsv(output / "filter_profiles.tsv")
    diagnostics = read_tsv(output / "candidate_diagnostics.tsv")
    require(len(profiles) == int(summary.get("n_satellite_annulus_profiles", -1)),
            "Satellite annulus row count differs from summary")
    require(len(diagnostics) == int(summary.get("n_candidate_diagnostics", -1)),
            "Candidate diagnostics count differs from summary")
    require(len({(row["cluster_id"], row["filter"], row["member_index"],
                  row["radial_bin_kpc"]) for row in profiles}) == len(profiles),
            "Satellite profile table has duplicate cluster/member/bin rows")
    require(len({(row["cluster_id"], row["filter"], row["member_index"])
                 for row in diagnostics}) == len(diagnostics),
            "Candidate diagnostics has duplicate cluster/member rows")

    floor = float(summary["method"]["coverage_floor"])
    min_controls = int(summary["method"]["minimum_controls_per_annulus"])
    profile_groups = defaultdict(list)
    for row in profiles:
        numeric = {key: float(row[key]) for key in (
            "target_sb_njy_arcsec2", "control_sb_njy_arcsec2", "delta_sb_njy_arcsec2",
            "target_err_nominal_njy_arcsec2", "target_coverage", "control_coverage")}
        require(all(math.isfinite(value) for value in numeric.values()),
                "Non-finite value in retained annulus")
        require(numeric["target_coverage"] >= floor
                and numeric["control_coverage"] >= floor,
                "Retained annulus is below its coverage floor")
        require(int(row["n_controls"]) >= min_controls,
                "Retained annulus has too few controls")
        require(math.isclose(numeric["target_sb_njy_arcsec2"]
                             - numeric["control_sb_njy_arcsec2"],
                             numeric["delta_sb_njy_arcsec2"],
                             rel_tol=1e-10, abs_tol=1e-10),
                "Target-control difference is inconsistent")
        key = (row["cluster_id"], row["filter"], row["cohort"], row["radial_bin_kpc"])
        profile_groups[key].append(row)

    expected_cluster_keys = set()
    for key, rows in profile_groups.items():
        cluster_id, band, cohort, radial_bin = key
        expected_cluster_keys.add(key)
        matched = [row for row in cluster_rows
                   if (row["cluster_id"], row["filter"], row["cohort"],
                       row["radial_bin_kpc"]) == key]
        require(len(matched) == 1, f"Missing/duplicate cluster summary for {key}")
        record = matched[0]
        require(int(record["n_satellites"]) == len(rows),
                f"Satellite count mismatch for {key}")
        for field in ("target_sb_njy_arcsec2", "control_sb_njy_arcsec2", "delta_sb_njy_arcsec2"):
            output_field = "mean_" + field.replace("_njy_arcsec2", "_njy_arcsec2")
            expected_mean = sum(float(row[field]) for row in rows) / len(rows)
            require(math.isclose(float(record[output_field]), expected_mean,
                                 rel_tol=1e-10, abs_tol=1e-10),
                    f"Cluster mean mismatch for {key}/{field}")
    require(len(cluster_rows) == len(expected_cluster_keys),
            "Cluster summary contains unsupported or missing annuli")

    filter_lookup = {(row["cohort"], row["filter"], row["radial_bin_kpc"]): row
                     for row in filter_rows}
    require(len(filter_lookup) == len(filter_rows), "Duplicate filter/bin summary rows")
    expected_filter_keys = set()
    for cohort, band, radial_bin in {(key[2], key[1], key[3]) for key in profile_groups}:
        selected = [row for row in cluster_rows
                    if row["cohort"] == cohort and row["filter"] == band
                    and row["radial_bin_kpc"] == radial_bin]
        key = (cohort, band, radial_bin)
        expected_filter_keys.add(key)
        require(key in filter_lookup, f"Missing filter summary for {key}")
        record = filter_lookup[key]
        require(int(record["n_clusters"]) == len(selected),
                f"Cluster count mismatch for {key}")
        expected_satellites = sum(int(row["n_satellites"]) for row in selected)
        require(int(record["n_satellites"]) == expected_satellites,
                f"Satellite count mismatch for {key}")
        mean = sum(float(row["mean_delta_sb_njy_arcsec2"]) for row in selected) / len(selected)
        require(math.isclose(float(record["equal_cluster_mean_delta_sb_njy_arcsec2"]),
                             mean, rel_tol=1e-10, abs_tol=1e-10),
                f"Equal-cluster mean mismatch for {key}")
        low, high = float(record["cluster_bootstrap_95_low"]), float(record["cluster_bootstrap_95_high"])
        require(math.isfinite(low) and math.isfinite(high) and low <= high,
                f"Invalid bootstrap interval for {key}")
    require(set(filter_lookup) == expected_filter_keys,
            "Filter table contains unsupported or missing annuli")

    groups = {}
    for key in sorted(expected_filter_keys):
        cohort, band, radial_bin = key
        record = filter_lookup[key]
        groups[f"{cohort}/{band}/{radial_bin}"] = {
            "clusters": int(record["n_clusters"]),
            "satellites": int(record["n_satellites"]),
        }
    result = {
        "status": "verified", "science_measurement": False,
        "frames": n_frames, "satellites_with_any_bin": len({
            (row["cluster_id"], row["filter"], row["member_index"]) for row in profiles}),
        "satellite_annulus_profiles": len(profiles), "diagnostic_rows": len(diagnostics),
        "verified_output_hashes": len(checked_hashes),
        "artifact_audit_verified": bool(audit_hash), "groups": groups,
    }
    print(json.dumps(result, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except ValueError as error:
        raise SystemExit(f"Validation failed: {error}")
