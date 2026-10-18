#!/usr/bin/env python3
"""Validate ICL/satellite association sensitivity across the full cut grid."""
import argparse
import csv
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path


FACTORS = (2.5, 3.5, 5.0)
THRESHOLDS = (0.25, 0.35, 0.40, 0.50)


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
    parser.add_argument("--icl-sensitivity", type=Path, required=True)
    parser.add_argument("--satellite", type=Path, required=True)
    parser.add_argument("--common-support", type=Path, required=True)
    parser.add_argument("--program", type=Path)
    args = parser.parse_args()
    output, icl, satellite, common = (
        args.output.resolve(), args.icl_sensitivity.resolve(),
        args.satellite.resolve(), args.common_support.resolve())
    summary_path = output / "summary.json"
    require(summary_path.is_file(), f"Missing summary: {summary_path}")
    summary = read_json(summary_path)
    require(summary.get("schema") == "hsc-icl-satellite-cut-sensitivity-v1"
            and summary.get("status") == "complete"
            and summary.get("science_measurement") is False,
            "Unexpected schema, incomplete run, or science_measurement is not false")
    require(summary.get("n_configurations") == 24 and summary.get("n_result_groups") == 48,
            "Expected 24 configurations and 48 all/center-eligible results")

    icl_summary_path = icl / "summary.json"
    icl_detail_path = icl / "selected_pair_profiles.tsv"
    sat_summary_path = satellite / "summary.json"
    sat_profiles_path = satellite / "cluster_profiles.tsv"
    common_summary_path = common / "summary.json"
    common_ids_path = common / "common_cluster_ids.tsv"
    for path in (icl_summary_path, icl_detail_path, sat_summary_path, sat_profiles_path,
                 common_summary_path, common_ids_path):
        require(summary.get("input_sha256", {}).get(str(path)) == sha256(path),
                f"Input checksum mismatch: {path}")

    checked = {}
    for name, expected in summary.get("output_sha256", {}).items():
        path = Path(name)
        require(path.is_file() and sha256(path) == expected,
                f"Output checksum mismatch: {path}")
        checked[path.name] = expected
    require({"association_sensitivity.tsv", "matched_cluster_values.tsv"}.issubset(checked),
            "Output checksum manifest is incomplete")
    if args.program:
        require(summary.get("program_sha256") == sha256(args.program),
                "Analysis program checksum mismatch")

    icl_summary = read_json(icl_summary_path)
    sat_summary = read_json(sat_summary_path)
    common_summary = read_json(common_summary_path)
    require(icl_summary.get("science_measurement") is False
            and sat_summary.get("science_measurement") is False
            and common_summary.get("science_measurement") is False,
            "An input unexpectedly claims a science measurement")
    common_ids = defaultdict(set)
    common_rows = read_tsv(common_ids_path)
    for row in common_rows:
        key = (row["cluster_id"], row["filter"])
        require(row["filter"] in ("i", "i2") and row["cluster_id"] not in common_ids[row["filter"]],
                f"Duplicate or unexpected common cluster: {key}")
        common_ids[row["filter"]].add(row["cluster_id"])
    require({band: len(common_ids[band]) for band in ("i", "i2")}
            == summary.get("n_common_satellite_clusters"),
            "Common satellite sample counts differ from summary")

    satellite_outer = {}
    for row in read_tsv(sat_profiles_path):
        if row["radial_bin_kpc"] == "50-75":
            key = (row["cluster_id"], row["filter"])
            require(key not in satellite_outer, f"Duplicate satellite outer row: {key}")
            satellite_outer[key] = row

    icl_details = read_tsv(icl_detail_path)
    expected_by_config = defaultdict(set)
    for row in icl_details:
        factor, threshold, band = float(row["kron_factor"]), float(row["coverage_floor"]), row["filter"]
        ident = row["id"]
        if ident in common_ids[band]:
            expected_by_config[(factor, threshold, band)].add(ident)
    result_rows = read_tsv(output / "association_sensitivity.tsv")
    detail_rows = read_tsv(output / "matched_cluster_values.tsv")
    require(len(result_rows) == 48, "Association result table must contain 48 groups")
    grouped = defaultdict(list)
    for row in detail_rows:
        key = (float(row["kron_factor"]), float(row["coverage_floor"]), row["filter"])
        ident = row["cluster_id"]
        require(ident in common_ids[row["filter"]], f"Non-common cluster in association: {key}/{ident}")
        require((ident, row["filter"]) in satellite_outer,
                f"Missing satellite profile for association: {key}/{ident}")
        expected_sat = float(satellite_outer[(ident, row["filter"])]["mean_delta_sb_njy_arcsec2"])
        require(math.isclose(float(row["satellite_delta_50_75"]), expected_sat,
                             rel_tol=1e-10, abs_tol=1e-10),
                f"Satellite outer value mismatch: {key}/{ident}")
        grouped[key].append(row)
    require(len({(float(row["kron_factor"]), float(row["coverage_floor"]),
                  row["filter"], row["cluster_id"]) for row in detail_rows}) == len(detail_rows),
            "Duplicate matched association cluster/configuration rows")

    results_by_key = {(float(row["kron_factor"]), float(row["coverage_floor"]),
                       row["filter"], row["selection"]): row for row in result_rows}
    require(len(results_by_key) == 48, "Duplicate result group")
    for factor in FACTORS:
        for threshold in THRESHOLDS:
            for band in ("i", "i2"):
                key = (factor, threshold, band)
                rows = grouped[key]
                require({row["cluster_id"] for row in rows} == expected_by_config[key],
                        f"Association join sample mismatch: {key}")
                all_result = results_by_key[(factor, threshold, band, "all_common_clusters")]
                center_result = results_by_key[(factor, threshold, band, "center_eligible")]
                require(int(all_result["n_clusters"]) == len(rows),
                        f"All-sample association count mismatch: {key}")
                eligible = [row for row in rows if row["center_eligible"] == "True"]
                require(int(center_result["n_clusters"]) == len(eligible),
                        f"Center-eligible association count mismatch: {key}")
                if len(rows) >= 5 and all_result["status"] == "exploratory":
                    rho = all_result["spearman_rho"]
                    require(rho != "" and math.isfinite(float(rho)),
                            f"Missing Spearman estimate: {key}")
                    require(int(all_result["valid_bootstrap_replicas"]) > 0,
                            f"No valid block bootstrap draws: {key}")

    result = {
        "status": "verified", "science_measurement": False,
        "configurations": len({(float(row["kron_factor"]), float(row["coverage_floor"]))
                                for row in result_rows}),
        "result_groups": len(result_rows), "matched_cluster_config_rows": len(detail_rows),
        "common_clusters": {band: len(common_ids[band]) for band in ("i", "i2")},
        "verified_output_hashes": len(checked),
    }
    print(json.dumps(result, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except ValueError as error:
        raise SystemExit(f"Validation failed: {error}")
