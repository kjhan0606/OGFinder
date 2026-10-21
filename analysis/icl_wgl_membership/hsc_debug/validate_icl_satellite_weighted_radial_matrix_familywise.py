#!/usr/bin/env python3
"""Independently validate a saved 168-cell CAMIRA-weight sensitivity run."""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

from audit_icl_satellite_radial_support_matrix import ICL_BINS, SATELLITE_BINS, sha256


ESTIMATORS = ("equal_satellite", "camira_weighted")
FILTERS = ("i", "i2")


def read_tsv(path):
    with Path(path).open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def holm_adjust(values):
    values = [float(x) for x in values]
    order = sorted(range(len(values)), key=values.__getitem__)
    adjusted = [0.0] * len(values)
    running = 0.0
    for rank, index in enumerate(order):
        running = max(running, (len(values) - rank) * values[index])
        adjusted[index] = min(1.0, running)
    return adjusted


def validate(run_dir: Path):
    run_dir = Path(run_dir).resolve()
    summary_path = run_dir / "summary.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    if (summary.get("schema") != "astrafek-icl-satellite-weighted-radial-matrix-familywise-v1"
            or summary.get("status") != "complete_exploratory_familywise_sensitivity"
            or summary.get("science_measurement") is not False
            or summary.get("n_expected_cells") != 168
            or summary.get("n_eligible_cells") != 168):
        raise ValueError("unexpected weighted radial-matrix result schema/status")
    associations = run_dir / "weighted_radial_matrix_association.tsv"
    support = run_dir / "weighted_radial_matrix_support.tsv"
    null_path = run_dir / "joint_maxT_null.tsv"
    for path in (associations, support, null_path):
        expected = summary.get("output_sha256", {}).get(str(path))
        if expected is None or sha256(path) != expected:
            raise ValueError(f"result hash mismatch: {path}")
    provenance = summary.get("provenance", {})
    for section in provenance.values():
        for item in section.values():
            if isinstance(item, dict) and "path" in item:
                path = Path(item["path"])
                if not path.is_file() or sha256(path) != item["sha256"]:
                    raise ValueError(f"source provenance hash mismatch: {path}")
    result_rows = read_tsv(associations)
    support_rows = read_tsv(support)
    null_rows = read_tsv(null_path)
    if len(result_rows) != 168 or len(support_rows) != 168 or len(null_rows) != 5000:
        raise ValueError("result/support/null table dimensions differ from the run contract")
    keys = [(r["satellite_estimator"], r["filter"], r["icl_radial_bin_kpc"],
             r["satellite_radial_bin_kpc"]) for r in result_rows]
    expected = {(estimator, band, icl, sat) for estimator in ESTIMATORS for band in FILTERS
                for icl in ICL_BINS for sat in SATELLITE_BINS}
    if len(keys) != len(set(keys)) or set(keys) != expected:
        raise ValueError("weighted association table does not contain the exact 168-cell family")
    support_keys = [(r["estimator"], r["filter"], r["icl_radial_bin_kpc"],
                     r["satellite_radial_bin_kpc"]) for r in support_rows]
    if len(support_keys) != len(set(support_keys)) or set(support_keys) != expected:
        raise ValueError("support table does not match the 168-cell family")
    support_map = {key: row for key, row in zip(support_keys, support_rows)}
    for estimator, band, icl, sat in expected:
        row = support_map[(estimator, band, icl, sat)]
        if int(row["n_common"]) < 10 or int(row["n_spatial_blocks"]) < 8:
            raise ValueError("an output cell violates the predeclared support gate")
    for band in FILTERS:
        for icl in ICL_BINS:
            for sat in SATELLITE_BINS:
                equal = support_map[("equal_satellite", band, icl, sat)]
                weighted = support_map[("camira_weighted", band, icl, sat)]
                if (equal["n_common"] != weighted["n_common"]
                        or equal["n_spatial_blocks"] != weighted["n_spatial_blocks"]):
                    raise ValueError("equal/weighted estimator support unexpectedly differs")
    max_null = [float(r["max_abs_t"]) for r in null_rows]
    if not all(math.isfinite(x) and x >= 0 for x in max_null):
        raise ValueError("max-T null contains invalid values")
    raw_p = [float(r["p_wild_cluster_raw"]) for r in result_rows]
    holm = holm_adjust(raw_p)
    for row, adjusted in zip(result_rows, holm):
        observed = abs(float(row["t_cluster_cr1"]))
        expected_p = (1 + sum(value >= observed for value in max_null)) / (len(max_null) + 1)
        if not math.isclose(float(row["p_holm_all_eligible_cells"]), adjusted,
                            rel_tol=0, abs_tol=1e-15):
            raise ValueError("saved Holm p-value does not reproduce")
        if not math.isclose(float(row["p_maxT_single_step_complete_null"]), expected_p,
                            rel_tol=0, abs_tol=1e-15):
            raise ValueError("saved max-T p-value does not reproduce")
        if float(row["p_holm_all_eligible_cells"]) < float(row["p_wild_cluster_raw"]):
            raise ValueError("Holm adjustment is smaller than the raw p-value")
    return {
        "status": "pass", "science_measurement": False,
        "n_cells": len(result_rows), "n_maxT_replicates": len(max_null),
        "n_raw_p_lt_0p05": sum(p < 0.05 for p in raw_p),
        "n_holm_p_lt_0p05": sum(float(r["p_holm_all_eligible_cells"]) < 0.05 for r in result_rows),
        "minimum_holm_p": min(float(r["p_holm_all_eligible_cells"]) for r in result_rows),
        "minimum_maxT_p": min(float(r["p_maxT_single_step_complete_null"]) for r in result_rows),
        "summary_sha256": sha256(summary_path),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("run_dir", type=Path)
    args = parser.parse_args()
    print(json.dumps(validate(args.run_dir), ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
