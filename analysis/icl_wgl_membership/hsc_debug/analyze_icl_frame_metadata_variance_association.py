#!/usr/bin/env python3
"""Relate approved HSC frame metadata to delivered coadd variance balance.

This is an offline, exploratory matched-pair analysis. It does not convert
frame metadata into pixel weights, infer source-Poisson noise, or change the
science-measurement gate.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
from evaluate_prestack_pair_support import pair_block_labels

METADATA_AUDIT = DATA / "icl_frame_metadata_row_audit_2026-11-06_02"
VARIANCE_DIR = DATA / "icl_variance_radial_balance_full783_2026-11-06_01"
DEFAULT_OUTPUT = DATA / "icl_frame_metadata_variance_association_2026-11-06_02"
EXPECTED_PLAN_SHA256 = "49aa4b7601dd42243209d7b3b394a4dd3b4cbf00d87d28cbb9467553c29f4ce3"
EXPECTED_METADATA_SHA256 = "2560e8ece327bdd5ae02108f5187b40efa6bf3fa8b602b929287da9b64dffaf5"
N_BOOT = 5000
SEED = 20261106

PREDICTORS = (
    "delta_log_frame_count",
    "delta_log_visit_count",
    "delta_log_sky_rms_flux_proxy",
    "delta_log_fluxmag0",
    "delta_log_gain_geomean",
)


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty table: {path}")
    with path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t",
                                lineterminator="\n", extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def summarize_center(row: dict[str, str]) -> dict:
    numeric = ("n_unique_frames", "n_unique_visits", "n_unique_ccds",
               "median_sigma_sky", "median_fluxmag0", "median_gain1",
               "median_gain2", "median_gain3", "median_gain4",
               "n_filter_i", "n_filter_i2")
    result = {key: float(row[key]) for key in numeric}
    result["sample_role"] = row["sample_role"]
    if any(not math.isfinite(result[key]) or result[key] <= 0
           for key in ("n_unique_frames", "n_unique_visits", "median_sigma_sky",
                       "median_fluxmag0", "median_gain1", "median_gain2",
                       "median_gain3", "median_gain4")):
        raise ValueError(f"invalid positive metadata value for {row['center_id']}")
    # A coarse per-frame sky RMS proxy in calibrated flux units. It is not a
    # coadd variance because frame weights, pixel/CCD mapping, and correlations
    # are not available in this export.
    result["sky_rms_flux_proxy"] = result["median_sigma_sky"] * 3.631e12 / result["median_fluxmag0"]
    result["gain_geomean"] = math.exp(float(np.mean(np.log([
        result[f"median_gain{i}"] for i in range(1, 5)]))))
    return result


def partial_corr(x: np.ndarray, y: np.ndarray, covariates: np.ndarray) -> float:
    design = np.column_stack([np.ones(len(x)), covariates])
    rx = x - design @ np.linalg.lstsq(design, x, rcond=None)[0]
    ry = y - design @ np.linalg.lstsq(design, y, rcond=None)[0]
    if np.std(rx) == 0 or np.std(ry) == 0:
        return float("nan")
    return float(np.corrcoef(rx, ry)[0, 1])


def block_interval(rows: list[dict], predictor: str, replicas: int, seed: int) -> dict:
    x = np.asarray([float(row[predictor]) for row in rows])
    y = np.asarray([float(row["log_target_control_variance_ratio"]) for row in rows])
    cov = np.column_stack([
        np.asarray([float(row["redshift"]) for row in rows]),
        np.asarray([float(row["log_richness"]) for row in rows]),
    ])
    labels = np.asarray([row["pair_block_5deg"] for row in rows], dtype=str)
    groups = np.unique(labels)
    if len(groups) < 6:
        raise ValueError(f"at least six connected spatial blocks required; got {len(groups)}")
    by_group = {group: np.flatnonzero(labels == group) for group in groups}
    rng = np.random.default_rng(seed)
    draws = []
    for _ in range(replicas):
        selected = rng.choice(groups, size=len(groups), replace=True)
        indices = np.concatenate([by_group[group] for group in selected])
        value = partial_corr(x[indices], y[indices], cov[indices])
        if np.isfinite(value):
            draws.append(value)
    if len(draws) < replicas * 0.9:
        raise ValueError(f"too few finite bootstrap draws: {len(draws)}/{replicas}")
    return {
        "n_pairs": len(rows),
        "n_connected_5deg_blocks": int(len(groups)),
        "largest_block_pairs": int(max(map(len, by_group.values()))),
        "partial_r_z_logrichness": partial_corr(x, y, cov),
        "block_bootstrap_ci95": np.percentile(draws, [2.5, 97.5]).tolist(),
        "bootstrap_fraction_positive": float(np.mean(np.asarray(draws) > 0)),
        "replicas_requested": replicas,
        "replicas_valid": len(draws),
        "seed": seed,
    }


def joint_familywise(rows: list[dict], replicas: int, seed: int) -> dict:
    """Single-step max-T intervals for five predictors in both filter strata."""
    target_coords = [{"ra_deg": float(row["target_ra_deg"]),
                      "dec_deg": float(row["target_dec_deg"])} for row in rows]
    control_coords = [{"ra_deg": float(row["control_ra_deg"]),
                       "dec_deg": float(row["control_dec_deg"])} for row in rows]
    _, labels, n_target_cells, n_control_cells, n_blocks, largest = pair_block_labels(
        target_coords, control_coords, np.ones(len(rows), dtype=bool))
    groups = np.unique(labels)
    by_group = {group: np.flatnonzero(labels == group) for group in groups}
    metrics = [(band, predictor) for band in ("i", "i2") for predictor in PREDICTORS]
    observed = []
    for band, predictor in metrics:
        subset = [row for row in rows if row["filter"] == band]
        x = np.asarray([float(row[predictor]) for row in subset])
        y = np.asarray([float(row["log_target_control_variance_ratio"]) for row in subset])
        cov = np.column_stack([
            np.asarray([float(row["redshift"]) for row in subset]),
            np.asarray([float(row["log_richness"]) for row in subset]),
        ])
        observed.append(partial_corr(x, y, cov))
    observed = np.asarray(observed)

    rng = np.random.default_rng(seed)
    draws = np.full((replicas, len(metrics)), np.nan, dtype=float)
    for replica in range(replicas):
        chosen = rng.choice(groups, size=len(groups), replace=True)
        indices = np.concatenate([by_group[group] for group in chosen])
        for metric_index, (band, predictor) in enumerate(metrics):
            sample_rows = [rows[int(index)] for index in indices
                           if rows[int(index)]["filter"] == band]
            if len(sample_rows) < 12:
                continue
            x = np.asarray([float(row[predictor]) for row in sample_rows])
            y = np.asarray([float(row["log_target_control_variance_ratio"])
                            for row in sample_rows])
            cov = np.column_stack([
                np.asarray([float(row["redshift"]) for row in sample_rows]),
                np.asarray([float(row["log_richness"]) for row in sample_rows]),
            ])
            draws[replica, metric_index] = partial_corr(x, y, cov)
    valid = np.all(np.isfinite(draws), axis=1)
    draws = draws[valid]
    if len(draws) < replicas * 0.9:
        raise ValueError(f"too few complete joint bootstrap draws: {len(draws)}/{replicas}")
    se = np.std(draws, axis=0, ddof=1)
    if np.any(se <= 0):
        raise ValueError("degenerate max-T bootstrap standard error")
    centered_t = np.abs((draws - observed) / se)
    max_t = np.max(centered_t, axis=1)
    critical = float(np.quantile(max_t, 0.95))
    inference = {}
    for index, (band, predictor) in enumerate(metrics):
        adjusted_p = (1 + int(np.count_nonzero(max_t >= abs(observed[index] / se[index])))) / (len(max_t) + 1)
        inference.setdefault(band, {})[predictor] = {
            "partial_r_z_logrichness": float(observed[index]),
            "maxT_simultaneous_ci95": [
                float(observed[index] - critical * se[index]),
                float(observed[index] + critical * se[index]),
            ],
            "maxT_adjusted_p": float(adjusted_p),
            "bootstrap_se": float(se[index]),
        }
    return {
        "status": "exploratory_single_step_maxT",
        "family": "all 10 filter-by-metadata associations jointly",
        "n_joint_pairs": len(rows),
        "n_target_cells": n_target_cells,
        "n_control_cells": n_control_cells,
        "n_connected_5deg_blocks": n_blocks,
        "largest_block_pairs": largest,
        "critical_maxT_95": critical,
        "replicas_requested": replicas,
        "replicas_valid": len(draws),
        "seed": seed,
        "intervals_by_filter": inference,
    }


def run(output: Path, replicas: int, seed: int) -> dict:
    if output.exists():
        raise FileExistsError(f"refusing to overwrite {output}")
    audit_path = METADATA_AUDIT / "summary.json"
    metadata_path = Path(json.loads(audit_path.read_text())["aggregate_tables"]["center_aggregates"])
    variance_path = VARIANCE_DIR / "target_control_radial_comparison.tsv"
    variance_summary_path = VARIANCE_DIR / "summary.json"
    audit = json.loads(audit_path.read_text(encoding="utf-8"))
    variance_summary = json.loads(variance_summary_path.read_text(encoding="utf-8"))
    if audit.get("status") != "complete_exploratory" or audit.get("plan_sha256") != EXPECTED_PLAN_SHA256:
        raise ValueError("frame metadata audit is missing, incomplete, or from an unexpected plan")
    if audit.get("result_files_sha256", {}).get("center_aggregates.tsv") != EXPECTED_METADATA_SHA256:
        raise ValueError("approved frame metadata aggregates differ from the audited file hash")
    if audit.get("n_frame_rows") != 67071 or audit.get("n_centers") != 1598:
        raise ValueError("frame metadata scope differs from the approved 67,071-row plan")
    if variance_summary.get("status") != "complete_diagnostic" or variance_summary.get("n_pairs") != 783:
        raise ValueError("full-sample coadd variance product is missing or has unexpected support")

    metadata_rows = read_tsv(metadata_path)
    metadata_by_id = {row["center_id"]: summarize_center(row) for row in metadata_rows}
    if len(metadata_by_id) != 1598:
        raise ValueError(f"expected 1,598 unique frame metadata centers; got {len(metadata_by_id)}")
    variance_rows = [row for row in read_tsv(variance_path)
                     if row["r_min_kpc"] == "200" and row["r_max_kpc"] == "300"]
    pair_keys = set()
    joined = []
    for row in variance_rows:
        key = (row["cluster_id"], row["filter"])
        if key in pair_keys:
            raise ValueError(f"duplicate target/filter pair in variance table: {key}")
        pair_keys.add(key)
        band = row["filter"]
        target_key = f"target:{row['cluster_id']}"
        control_key = f"control:{row['control_id']}"
        if target_key not in metadata_by_id or control_key not in metadata_by_id:
            raise ValueError(f"missing frame metadata for paired endpoints: {key}")
        target, control = metadata_by_id[target_key], metadata_by_id[control_key]
        record = {field: row[field] for field in (
            "cluster_id", "control_id", "filter", "redshift", "log_richness",
            "target_variance", "control_variance", "log_target_control_variance_ratio",
            "pair_block_5deg", "target_ra_deg", "target_dec_deg",
            "control_ra_deg", "control_dec_deg")}
        record["target_frame_count"] = target["n_unique_frames"]
        record["control_frame_count"] = control["n_unique_frames"]
        record["target_visit_count"] = target["n_unique_visits"]
        record["control_visit_count"] = control["n_unique_visits"]
        record["target_ccd_count"] = target["n_unique_ccds"]
        record["control_ccd_count"] = control["n_unique_ccds"]
        record["target_median_sigma_sky"] = target["median_sigma_sky"]
        record["control_median_sigma_sky"] = control["median_sigma_sky"]
        record["target_sky_rms_flux_proxy"] = target["sky_rms_flux_proxy"]
        record["control_sky_rms_flux_proxy"] = control["sky_rms_flux_proxy"]
        record["target_median_fluxmag0"] = target["median_fluxmag0"]
        record["control_median_fluxmag0"] = control["median_fluxmag0"]
        record["target_gain_geomean"] = target["gain_geomean"]
        record["control_gain_geomean"] = control["gain_geomean"]
        record["delta_log_frame_count"] = math.log(target["n_unique_frames"] / control["n_unique_frames"])
        record["delta_log_visit_count"] = math.log(target["n_unique_visits"] / control["n_unique_visits"])
        record["delta_log_sky_rms_flux_proxy"] = math.log(target["sky_rms_flux_proxy"] / control["sky_rms_flux_proxy"])
        record["delta_log_fluxmag0"] = math.log(target["median_fluxmag0"] / control["median_fluxmag0"])
        record["delta_log_gain_geomean"] = math.log(target["gain_geomean"] / control["gain_geomean"])
        joined.append(record)
    if len(joined) != 774:
        raise ValueError(f"expected 774 supported paired 200-300 kpc rows; got {len(joined)}")

    results = {"filters": {}, "predictors": list(PREDICTORS)}
    for band_index, band in enumerate(("i", "i2")):
        band_rows = [row for row in joined if row["filter"] == band]
        if len(band_rows) != (299 if band == "i" else 475):
            raise ValueError(f"unexpected {band}-band paired support: {len(band_rows)}")
        results["filters"][band] = {}
        for predictor_index, predictor in enumerate(PREDICTORS):
            result = block_interval(band_rows, predictor, replicas,
                                    seed + 100 * band_index + predictor_index)
            results["filters"][band][predictor] = result
    results["familywise_maxT"] = joint_familywise(joined, replicas, seed + 1000)
    results.update({
        "schema": "icl-frame-metadata-vs-coadd-variance-v1",
        "status": "complete_exploratory",
        "science_measurement": False,
        "n_pairs_by_filter": {band: sum(row["filter"] == band for row in joined)
                              for band in ("i", "i2")},
        "n_paired_rows": len(joined),
        "radius_kpc": [200, 300],
        "outcome": "log(target/control median delivered coadd VARIANCE)",
        "adjustment": ["redshift", "log richness"],
        "resampling": "5-degree connected target/control spatial blocks within filter",
        "limitations": [
            "Frame membership does not identify pixel-to-frame weights or which amplifier covers each pixel.",
            "The normalized frame sky RMS is only a coarse proxy, not a coadd noise model.",
            "Five metadata predictors are exploratory; intervals are pointwise and not familywise-adjusted.",
            "Delivered coadd variance omits source-Poisson and correlated-noise validation.",
            "The analysis cannot update ICL flux errors or the science-measurement gate.",
        ],
        "input_hashes": {
            "metadata_audit": sha256(audit_path),
            "metadata_aggregates": sha256(metadata_path),
            "variance_summary": sha256(variance_summary_path),
            "variance_comparison": sha256(variance_path),
        },
        "script_sha256": sha256(Path(__file__).resolve()),
    })
    output.mkdir(parents=True, exist_ok=False)
    write_tsv(output / "joined_pairs.tsv", joined)
    summary_path = output / "summary.json"
    summary_path.write_text(json.dumps(results, indent=2, ensure_ascii=False) + "\n",
                            encoding="utf-8")
    results["output_sha256"] = {"joined_pairs.tsv": sha256(output / "joined_pairs.tsv")}
    summary_path.write_text(json.dumps(results, indent=2, ensure_ascii=False) + "\n",
                            encoding="utf-8")
    return results


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--replicas", type=int, default=N_BOOT)
    parser.add_argument("--seed", type=int, default=SEED)
    args = parser.parse_args()
    result = run(args.output, args.replicas, args.seed)
    print(json.dumps({"status": result["status"], "output": str(args.output),
                      "n_paired_rows": result["n_paired_rows"],
                      "n_pairs_by_filter": result["n_pairs_by_filter"],
                      "science_measurement": result["science_measurement"]}, indent=2))


if __name__ == "__main__":
    main()
