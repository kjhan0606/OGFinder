#!/usr/bin/env python3
"""Verify a saved satellite-centered pilot's provenance and table invariants."""
import argparse
import csv
import hashlib
import json
import math
from collections import Counter
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


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    output = args.output.resolve()
    summary_path = output / "summary.json"
    if not summary_path.is_file():
        parser.error(f"Missing summary: {summary_path}")
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    if summary.get("schema") != "hsc-satellite-outer-profile-pilot-v1":
        parser.error("Unexpected pilot schema")
    if summary.get("status") != "complete" or summary.get("science_measurement") is not False:
        parser.error("Pilot is incomplete or its science gate is not false")
    if summary.get("n_frames_requested") != 50 or summary.get("n_frames_completed") != 50:
        parser.error("Expected all 50 validation frames")
    if summary.get("n_frame_errors") != 0 or summary.get("errors"):
        parser.error("Pilot reports frame errors")
    if len(summary.get("frame_records", [])) != 50:
        parser.error("Frame provenance count differs from 50")
    if any(not record.get("image_sha256") or not record.get("mask_sha256")
           for record in summary["frame_records"]):
        parser.error("At least one frame lacks image/mask checksums")

    checked_hashes = {}
    for name, expected in summary.get("output_sha256", {}).items():
        path = Path(name)
        if not path.is_file() or sha256(path) != expected:
            parser.error(f"Output checksum mismatch: {path}")
        checked_hashes[path.name] = expected
    required_outputs = {"satellite_profiles.tsv", "cluster_profiles.tsv",
                        "filter_profiles.tsv", "candidate_diagnostics.tsv",
                        "satellite_profile_review.png"}
    if not required_outputs.issubset(checked_hashes):
        parser.error("Output checksum manifest omits a required table")

    profiles = read_tsv(output / "satellite_profiles.tsv")
    cluster_profiles = read_tsv(output / "cluster_profiles.tsv")
    filter_profiles = read_tsv(output / "filter_profiles.tsv")
    diagnostics = read_tsv(output / "candidate_diagnostics.tsv")
    edges = summary["method"]["satellite_radius_bins_kpc"]
    n_bins = len(edges) - 1
    n_satellites = int(summary["n_satellite_profiles"])
    if len(profiles) != n_satellites * n_bins:
        parser.error("Satellite profile row count does not equal satellites × radial bins")
    if Counter(row["status"] for row in diagnostics).get("retained", 0) != n_satellites:
        parser.error("Candidate diagnostics retained count differs from summary")
    floor = float(summary["method"]["coverage_floor"])
    min_controls = int(summary["method"]["minimum_controls"])
    fields = ("target_sb_njy_arcsec2", "control_sb_njy_arcsec2",
              "delta_sb_njy_arcsec2", "target_coverage", "control_coverage")
    for row in profiles:
        values = [float(row[key]) for key in fields]
        if not all(math.isfinite(value) for value in values):
            parser.error("Non-finite value in satellite profile table")
        if float(row["target_coverage"]) < floor or float(row["control_coverage"]) < floor:
            parser.error("A retained radial bin falls below its coverage floor")
        if int(row["n_controls"]) < min_controls:
            parser.error("A retained satellite has too few valid controls")
        if not math.isclose(values[0] - values[1], values[2], rel_tol=1e-10, abs_tol=1e-10):
            parser.error("Saved target-control profile difference is inconsistent")

    actual_groups = {}
    cluster_members = {}
    for cohort, band in sorted({(row["cohort"], row["filter"]) for row in profiles}):
        selected = [row for row in profiles if row["cohort"] == cohort and row["filter"] == band]
        galaxies = {(row["cluster_id"], row["member_index"]) for row in selected}
        clusters = {row["cluster_id"] for row in selected}
        actual_groups[(cohort, band)] = len(galaxies)
        cluster_members[(cohort, band)] = len(clusters)
    for cohort, bands in summary["filters"].items():
        for band, record in bands.items():
            key = (cohort, band)
            if record["n_satellites"] != actual_groups.get(key, 0):
                parser.error(f"Satellite count mismatch for {cohort}/{band}")
            if record["n_clusters"] != cluster_members.get(key, 0):
                parser.error(f"Cluster count mismatch for {cohort}/{band}")
    if len(filter_profiles) != sum(len(bands) for bands in summary["filters"].values()) * n_bins:
        parser.error("Equal-cluster filter table row count differs from summary groups")
    expected_cluster_rows = sum(record["n_clusters"] for bands in summary["filters"].values()
                                for record in bands.values()) * n_bins
    if len(cluster_profiles) != expected_cluster_rows:
        parser.error("Per-cluster profile table row count differs from summary groups")

    result = {
        "status": "verified", "science_measurement": False,
        "frames": len(summary["frame_records"]), "satellites": n_satellites,
        "radial_bins_per_satellite": n_bins, "profile_rows": len(profiles),
        "diagnostic_rows": len(diagnostics), "verified_output_hashes": len(checked_hashes),
        "qualitative_plot": "satellite_profile_review.png",
        "groups": {f"{cohort}/{band}": {"satellites": n, "clusters": cluster_members[(cohort, band)]}
                   for (cohort, band), n in actual_groups.items()},
    }
    print(json.dumps(result, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
