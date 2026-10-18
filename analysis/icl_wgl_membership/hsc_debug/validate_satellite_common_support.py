#!/usr/bin/env python3
"""Validate common-cluster satellite profile and ICL association outputs."""
import argparse
import csv
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path

import numpy as np


BINS = ("5-10", "10-15", "15-20", "20-30", "30-40", "40-50", "50-75")


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
    parser.add_argument("--satellite", type=Path, required=True)
    parser.add_argument("--icl", type=Path, required=True)
    parser.add_argument("--program", type=Path)
    args = parser.parse_args()
    output, satellite, icl = (args.output.resolve(), args.satellite.resolve(), args.icl.resolve())
    summary_path = output / "summary.json"
    require(summary_path.is_file(), f"Missing summary: {summary_path}")
    summary = read_json(summary_path)
    require(summary.get("schema") == "hsc-satellite-common-cluster-support-v1"
            and summary.get("status") == "complete"
            and summary.get("science_measurement") is False,
            "Unexpected schema, incomplete run, or science_measurement is not false")
    require(summary.get("block_width_deg") == 5.0 and summary.get("replicas", 0) >= 100,
            "Unexpected spatial bootstrap configuration")

    sat_summary_path = satellite / "summary.json"
    sat_profiles_path = satellite / "cluster_profiles.tsv"
    icl_summary_path = icl / "summary.json"
    icl_pairs_path = icl / "pair_metrics_k2p5_c40_r50_300.tsv"
    manifest_path = Path(read_json(sat_summary_path)["input_manifest"])
    input_paths = (sat_summary_path, sat_profiles_path, manifest_path,
                   icl_summary_path, icl_pairs_path)
    for path in input_paths:
        require(summary.get("input_sha256", {}).get(str(path)) == sha256(path),
                f"Input checksum mismatch: {path}")

    output_hashes = summary.get("output_sha256", {})
    checked = {}
    for name, expected in output_hashes.items():
        path = Path(name)
        require(path.is_file() and sha256(path) == expected,
                f"Output checksum mismatch: {path}")
        checked[path.name] = expected
    required_outputs = {"common_support_profile.tsv", "profile_covariance.tsv",
                        "common_cluster_ids.tsv", "icl_satellite_common_support.tsv",
                        "association_results.tsv", "profile_bootstrap_covariance.npz"}
    require(required_outputs.issubset(checked), "Output checksum manifest is incomplete")
    if args.program:
        require(summary.get("program_sha256") == sha256(args.program),
                "Analysis program checksum mismatch")

    sat_summary = read_json(sat_summary_path)
    icl_summary = read_json(icl_summary_path)
    require(sat_summary.get("status") == "complete"
            and sat_summary.get("science_measurement") is False
            and icl_summary.get("status") == "complete"
            and icl_summary.get("science_gates", {}).get("science_measurement") is False,
            "An input is incomplete or not exploratory")
    source = read_tsv(sat_profiles_path)
    source_by_band_bin = defaultdict(dict)
    for row in source:
        key = (row["filter"], row["radial_bin_kpc"])
        ident = str(row["cluster_id"])
        require(ident not in source_by_band_bin[key], f"Duplicate source profile: {key}/{ident}")
        source_by_band_bin[key][ident] = row

    members = read_tsv(output / "common_cluster_ids.tsv")
    member_by_band = defaultdict(dict)
    for row in members:
        key = (row["filter"], row["cluster_id"])
        require(row["cluster_id"] not in member_by_band[row["filter"]],
                f"Duplicate common cluster ID: {key}")
        member_by_band[row["filter"]][row["cluster_id"]] = row

    profile_rows = read_tsv(output / "common_support_profile.tsv")
    require(len(profile_rows) == 14, "Expected two filters across seven radial bins")
    profile_by_key = {(row["filter"], row["radial_bin_kpc"]): row for row in profile_rows}
    require(len(profile_by_key) == len(profile_rows), "Duplicate common profile summary row")
    covariance_rows = read_tsv(output / "profile_covariance.tsv")
    require(len(covariance_rows) == 98, "Expected two 7x7 profile covariance matrices")
    covariance = defaultdict(dict)
    for row in covariance_rows:
        value = float(row["bootstrap_covariance_njy2_arcsec4"])
        require(math.isfinite(value), "Non-finite profile covariance")
        covariance[row["filter"]][(row["radial_bin_kpc"], row["other_radial_bin_kpc"])] = value
    covariance_quality = {}
    for band in ("i", "i2"):
        expected_common = set.intersection(*(set(source_by_band_bin[(band, radial_bin)])
                                             for radial_bin in BINS))
        actual_common = set(member_by_band[band])
        require(actual_common == expected_common,
                f"Common cluster intersection mismatch for {band}")
        require(len(actual_common) == summary["n_common_clusters"][band],
                f"Common cluster count mismatch for {band}")
        for radial_bin in BINS:
            record = profile_by_key.get((band, radial_bin))
            require(record is not None, f"Missing profile summary: {band}/{radial_bin}")
            values = [float(source_by_band_bin[(band, radial_bin)][ident][
                "mean_delta_sb_njy_arcsec2"]) for ident in sorted(expected_common)]
            expected_mean = float(np.mean(values))
            require(int(record["n_common_clusters"]) == len(expected_common),
                    f"Per-bin common count mismatch: {band}/{radial_bin}")
            require(math.isclose(float(record["mean_equal_cluster_delta_sb_njy_arcsec2"]),
                                 expected_mean, rel_tol=1e-10, abs_tol=1e-10),
                    f"Equal-cluster mean mismatch: {band}/{radial_bin}")
            require(math.isfinite(float(record["spatial_block_bootstrap_95_low"]))
                    and math.isfinite(float(record["spatial_block_bootstrap_95_high"]))
                    and float(record["spatial_block_bootstrap_95_low"])
                    <= float(record["spatial_block_bootstrap_95_high"]),
                    f"Invalid block-bootstrap interval: {band}/{radial_bin}")
        for first in BINS:
            for second in BINS:
                a = covariance[band].get((first, second))
                b = covariance[band].get((second, first))
                require(a is not None and b is not None
                        and math.isclose(a, b, rel_tol=1e-8, abs_tol=1e-10),
                        f"Covariance is incomplete or asymmetric: {band}/{first}/{second}")
        matrix = np.asarray([[covariance[band][(first, second)] for second in BINS]
                             for first in BINS], dtype=float)
        eigenvalues = np.linalg.eigvalsh(matrix)
        rank = int(np.linalg.matrix_rank(matrix))
        require(rank == len(BINS) and float(eigenvalues.min()) > 0,
                f"Profile covariance is not positive definite/full rank for {band}")
        covariance_quality[band] = {
            "rank": rank, "minimum_eigenvalue": float(eigenvalues.min()),
            "condition_number": float(np.linalg.cond(matrix)),
        }

    association = read_tsv(output / "icl_satellite_common_support.tsv")
    require(len({(row["cluster_id"], row["filter"]) for row in association}) == len(association),
            "Duplicate ICL/satellite association row")
    pair_rows = read_tsv(icl_pairs_path)
    pair_keys = {(str(row["id"]), row["filter"]) for row in pair_rows}
    association_keys = {(row["cluster_id"], row["filter"]) for row in association}
    expected_keys = {(ident, band) for band in ("i", "i2")
                     for ident in member_by_band[band] if (ident, band) in pair_keys}
    require(association_keys == expected_keys, "ICL join does not match common cluster support")
    require(summary["n_matched_icl_clusters"] == {
        band: sum(row["filter"] == band for row in association) for band in ("i", "i2")
    }, "ICL matched counts differ from summary")

    results = read_tsv(output / "association_results.tsv")
    require(len(results) == 4, "Expected all/center-eligible rows for both filters")
    result_by_key = {(row["filter"], row["selection"]): row for row in results}
    require(len(result_by_key) == 4, "Duplicate association result groups")
    for band in ("i", "i2"):
        rows = [row for row in association if row["filter"] == band]
        eligible = [row for row in rows if row["center_eligible"] == "True"]
        require(int(result_by_key[(band, "all_common_annulus_clusters")]["n_clusters"])
                == len(rows), f"Association count mismatch for {band}")
        require(int(result_by_key[(band, "center_eligible")]["n_clusters"])
                == len(eligible), f"Center-eligible count mismatch for {band}")

    result = {
        "status": "verified", "science_measurement": False,
        "common_clusters": {band: len(member_by_band[band]) for band in ("i", "i2")},
        "matched_icl_clusters": summary["n_matched_icl_clusters"],
        "profile_bins_per_filter": 7, "verified_covariance_elements": len(covariance_rows),
        "covariance_quality": covariance_quality,
        "verified_output_hashes": len(checked),
    }
    print(json.dumps(result, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except ValueError as error:
        raise SystemExit(f"Validation failed: {error}")
