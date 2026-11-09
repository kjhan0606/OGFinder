#!/usr/bin/env python3
"""Familywise block-bootstrap intervals for full-sample radial variance ratios.

This uses only the pinned output table from the completed 783-pair radial
variance analysis. It jointly resamples spatial blocks across filters and
annuli so the maximum-t interval accounts for the 20 correlated comparisons.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
from evaluate_prestack_pair_support import pair_block_labels


DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
PARENT = DATA / "icl_variance_radial_balance_full783_2026-11-06_01"
PARENT_SUMMARY = PARENT / "summary.json"
PARENT_TABLE = PARENT / "target_control_radial_comparison.tsv"
OUTPUT = DATA / "icl_variance_radial_balance_familywise_2026-11-06_02"
EXPECTED_SUMMARY_SHA256 = "f155075dfc5bfdc978fdcff9d3252d369ee81d5ef1af23ae0371aef1d4d04632"
EXPECTED_TABLE_SHA256 = "db77c9a32797aafa79efaac7f7e18d7e14a8f79c96b18e8fe296aceb6f35830d"
N_BOOT = 10000
SEED = 20261107
RADII = ((0, 50), (50, 75), (75, 100), (100, 150), (150, 200),
         (200, 300), (300, 400), (400, 500), (500, 600), (600, 900))


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def run(output: Path, replicas: int = N_BOOT) -> dict:
    output = output.resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite {output}")
    if replicas < 5000:
        raise ValueError("use at least 5,000 block bootstrap replicates")
    if sha256(PARENT_SUMMARY) != EXPECTED_SUMMARY_SHA256:
        raise ValueError("pinned full-sample variance summary changed")
    if sha256(PARENT_TABLE) != EXPECTED_TABLE_SHA256:
        raise ValueError("pinned full-sample pair table changed")
    parent = json.loads(PARENT_SUMMARY.read_text(encoding="utf-8"))
    if (parent.get("status") != "complete_diagnostic"
            or parent.get("sample_scope") != "all783"
            or parent.get("n_pairs") != 783
            or parent.get("science_measurement") is not False
            or parent.get("network_used") is not False):
        raise ValueError("parent result does not match expected offline exploratory scope")
    expected_table_hash = parent.get("output_sha256", {}).get("target_control_radial_comparison.tsv")
    if sha256(PARENT_TABLE) != expected_table_hash:
        raise ValueError("parent table hash does not match its summary")

    rows = read_rows(PARENT_TABLE)
    pairs = {}
    pair_coordinates = {}
    for row in rows:
        key = (str(row["cluster_id"]), row["filter"])
        pairs.setdefault(key, {})[(int(row["r_min_kpc"]), int(row["r_max_kpc"]))] = row
        pair_coordinates.setdefault(key, {
            "target": {"ra_deg": float(row["target_ra_deg"]), "dec_deg": float(row["target_dec_deg"])},
            "control": {"ra_deg": float(row["control_ra_deg"]), "dec_deg": float(row["control_dec_deg"])},
        })
    keys = sorted(pairs, key=lambda k: (k[1], int(k[0])))
    if len(keys) != 783:
        raise ValueError(f"expected 783 target/filter pairs, got {len(keys)}")
    target_coords = [pair_coordinates[key]["target"] for key in keys]
    control_coords = [pair_coordinates[key]["control"] for key in keys]
    _, labels, n_target_cells, n_control_cells, n_blocks, largest = pair_block_labels(
        target_coords, control_coords, np.ones(len(keys), dtype=bool)
    )
    if n_blocks < 6:
        raise ValueError(f"only {n_blocks} connected blocks across both filters")
    block_values = np.unique(labels)
    members = {label: np.flatnonzero(labels == label) for label in block_values}
    columns = []
    metrics = []
    for band in ("i", "i2"):
        for lo, hi in RADII:
            vector = np.full(len(keys), np.nan)
            for index, key in enumerate(keys):
                if key[1] != band:
                    continue
                row = pairs[key].get((lo, hi))
                if row is not None:
                    value = float(row["log_target_control_variance_ratio"])
                    vector[index] = value
            if np.isfinite(vector).sum() < 20:
                raise ValueError(f"insufficient support for {band} {lo}-{hi} kpc")
            columns.append(vector)
            metrics.append({"filter": band, "r_min_kpc": lo, "r_max_kpc": hi,
                            "n_pairs": int(np.isfinite(vector).sum()),
                            "point_log_ratio": float(np.nanmedian(vector))})
    matrix = np.column_stack(columns)
    point = np.asarray([metric["point_log_ratio"] for metric in metrics])
    rng = np.random.default_rng(SEED)
    boot = np.empty((replicas, len(metrics)), dtype=np.float64)
    for replicate in range(replicas):
        draw = rng.choice(block_values, size=len(block_values), replace=True)
        indices = np.concatenate([members[label] for label in draw])
        boot[replicate] = np.nanmedian(matrix[indices], axis=0)
    se = np.std(boot, axis=0, ddof=1)
    if np.any(~np.isfinite(se)) or np.any(se <= 0):
        raise ValueError("invalid bootstrap standard error for one or more metrics")
    max_t = np.nanmax(np.abs((boot - point) / se), axis=1)
    critical = float(np.quantile(max_t, 0.95))
    for i, metric in enumerate(metrics):
        pointwise = np.quantile(boot[:, i], [0.025, 0.975])
        low, high = point[i] - critical * se[i], point[i] + critical * se[i]
        metric.update({
            "bootstrap_se_log_ratio": float(se[i]),
            "pointwise_95_ratio_ci": [float(np.exp(pointwise[0])), float(np.exp(pointwise[1]))],
            "familywise_95_log_ratio_ci": [float(low), float(high)],
            "familywise_95_ratio_ci": [float(np.exp(low)), float(np.exp(high))],
            "familywise_interval_includes_equal_variance": bool(low <= 0 <= high),
        })
    output.mkdir(parents=True, exist_ok=False)
    metrics_path = output / "familywise_intervals.tsv"
    with metrics_path.open("x", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(metrics[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(metrics)
    summary = {
        "schema": "astrafek-icl-radial-variance-familywise-v1",
        "status": "complete_diagnostic",
        "science_measurement": False,
        "network_used": False,
        "sample_scope": "all 783 same-filter target/control pairs (302 i, 481 i2)",
        "n_spatial_blocks_joint_across_filters": n_blocks,
        "largest_block_pairs": largest,
        "n_target_sky_cells": n_target_cells,
        "n_control_sky_cells": n_control_cells,
        "n_metrics": len(metrics),
        "bootstrap": {"replicates": replicas, "seed": SEED,
                      "method": "joint connected-5-degree-block resampling across both filters and 20 correlated radial metrics; centered max-absolute-t simultaneous intervals"},
        "max_t_critical_95": critical,
        "n_familywise_intervals_excluding_equal_variance": sum(
            not metric["familywise_interval_includes_equal_variance"] for metric in metrics),
        "inputs_sha256": {str(PARENT_SUMMARY): sha256(PARENT_SUMMARY),
                          str(PARENT_TABLE): sha256(PARENT_TABLE)},
        "output_sha256": {metrics_path.name: sha256(metrics_path)},
        "limitations": [
            "connected 5-degree blocks are approximate equal-area cells, not validated HSC survey regions",
            "the max-t bands adjust the 20 radial variance contrasts, not every downstream scientific analysis",
            "the delivered variance plane does not supply effective gain for synthetic-source Poisson noise",
            "these diagnostic intervals do not define the covariance of the masked radial ICL estimator",
            "science_measurement remains false",
        ],
    }
    summary_path = output / "summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True, allow_nan=False) + "\n",
                            encoding="utf-8")
    print(json.dumps({"status": summary["status"], "n_metrics": len(metrics),
                      "familywise_exclusions": summary["n_familywise_intervals_excluding_equal_variance"],
                      "max_t_critical_95": critical, "science_measurement": False}, indent=2))
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--bootstrap", type=int, default=N_BOOT)
    args = parser.parse_args()
    run(args.output, args.bootstrap)


if __name__ == "__main__":
    main()
