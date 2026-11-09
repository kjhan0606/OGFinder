#!/usr/bin/env python3
"""Audit how truncating saved coadd residual-lag kernels changes noise factors.

Uses the independently verified 9x9 lag-covariance measurements for 783
existing same-filter target/control pairs. For each square cutoff radius 0..4
pixels it sums the symmetric variance-normalized kernel, compares target and
control values, and estimates paired connected-5-degree-block intervals.
This is a finite-kernel diagnostic, not a full aperture covariance correction.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
import sys
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from analyze_icl_coadd_covariance_paired_blocks import (
    DEFAULT_AUDIT, DEFAULT_SAMPLE, audit_inputs, read_tsv, sha256,
)
from evaluate_prestack_pair_support import pair_block_labels

DEFAULT_OUTPUT = DEFAULT_SAMPLE / "icl_lag_kernel_truncation_2026-11-06_01"
LAG_RADIUS = 4
EXPECTED_PAIRS = {"i": 302, "i2": 481}


def write_tsv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty table: {path}")
    with path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t",
                                lineterminator="\n", extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def load_kernels(audit_dir: Path) -> tuple[dict, Path]:
    lag_path = audit_dir / "lag_correlation.tsv"
    rows = read_tsv(lag_path)
    kernels = {}
    for row in rows:
        key = (str(row["cluster_id"]), str(row["role"]), str(row["filter"]))
        lag = (int(row["dy"]), int(row["dx"]))
        value = float(row["median_variance_normalized_covariance"])
        if not np.isfinite(value) or lag in kernels.setdefault(key, {}):
            raise ValueError(f"invalid or duplicate lag value for {key}/{lag}")
        kernels[key][lag] = value
    expected_lags = {(dy, dx) for dy in range(-LAG_RADIUS, LAG_RADIUS + 1)
                     for dx in range(-LAG_RADIUS, LAG_RADIUS + 1)}
    if len(kernels) != 2 * sum(EXPECTED_PAIRS.values()):
        raise ValueError(f"unexpected kernel count: {len(kernels)}")
    symmetric = {}
    for key, values in kernels.items():
        if set(values) != expected_lags:
            raise ValueError(f"incomplete 9x9 kernel for {key}")
        matrix = np.empty((9, 9), dtype=float)
        for dy, dx in expected_lags:
            matrix[dy + 4, dx + 4] = 0.5 * (values[(dy, dx)] + values[(-dy, -dx)])
        symmetric[key] = matrix
    return symmetric, lag_path


def connected_block_bootstrap(target: np.ndarray, control: np.ndarray,
                              labels: np.ndarray, *, replicas: int,
                              seed: int) -> dict:
    groups = np.unique(labels)
    members = [np.flatnonzero(labels == label) for label in groups]
    if len(members) < 3 or sum(map(len, members)) != len(target):
        raise ValueError("invalid connected-block partition")
    delta = target - control
    rng = np.random.default_rng(seed)
    draws_target = np.empty(replicas)
    draws_control = np.empty(replicas)
    draws_paired_delta = np.empty(replicas)
    for index in range(replicas):
        draw = rng.integers(0, len(members), size=len(members))
        selected = np.concatenate([members[int(i)] for i in draw])
        draws_target[index] = np.median(target[selected])
        draws_control[index] = np.median(control[selected])
        draws_paired_delta[index] = np.median(delta[selected])
    return {
        "n_blocks": len(members),
        "largest_block_pairs": max(map(len, members)),
        "target_median": float(np.median(target)),
        "control_median": float(np.median(control)),
        "paired_target_minus_control_median": float(np.median(delta)),
        "paired_delta_block_bootstrap_ci95": [
            float(np.quantile(draws_paired_delta, .025)),
            float(np.quantile(draws_paired_delta, .975)),
        ],
        "difference_of_medians_block_bootstrap_ci95": [
            float(np.quantile(draws_target - draws_control, .025)),
            float(np.quantile(draws_target - draws_control, .975)),
        ],
        "replicas": replicas,
        "seed": seed,
    }


def analyze(sample: Path, audit_dir: Path, output: Path,
            replicas: int, seed: int) -> dict:
    sample, audit_dir, output = map(Path, (sample, audit_dir, output))
    output = output.resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite output: {output}")
    upstream, frame_rows, manifest, control_manifest = audit_inputs(sample, audit_dir)
    kernels, lag_path = load_kernels(audit_dir)
    target_by_id = {str(row["id"]): row for row in manifest["clusters"]}
    control_by_target = {
        str(row["target_cluster_id"]): row for row in control_manifest["clusters"]
    }
    frame_lookup = {
        (str(row["cluster_id"]), str(row["role"]), str(row["filter"])): row
        for row in frame_rows
    }
    by_filter, pair_rows, lag_rows = {}, [], []
    for filter_index, band in enumerate(("i", "i2")):
        ids = sorted(ident for ident, role, filt in frame_lookup
                     if role == "target" and filt == band)
        if len(ids) != EXPECTED_PAIRS[band]:
            raise ValueError(f"unexpected {band}-band pair support: {len(ids)}")
        target_positions, control_positions = [], []
        per_pair = []
        for ident in ids:
            target = target_by_id.get(ident)
            control = control_by_target.get(ident)
            tframe = frame_lookup.get((ident, "target", band))
            cframe = frame_lookup.get((ident, "matched_control", band))
            tkernel = kernels.get((ident, "target", band))
            ckernel = kernels.get((ident, "matched_control", band))
            if any(value is None for value in (target, control, tframe, cframe, tkernel, ckernel)):
                raise ValueError(f"incomplete matched kernel pair: {ident}/{band}")
            target_positions.append({"ra_deg": float(target["ra_deg"]),
                                     "dec_deg": float(target["dec_deg"])})
            control_positions.append({"ra_deg": float(control["ra_deg"]),
                                      "dec_deg": float(control["dec_deg"])})
            row = {"cluster_id": ident, "filter": band}
            for radius in range(LAG_RADIUS + 1):
                region = np.zeros((9, 9), dtype=bool)
                for iy, dy in enumerate(range(-4, 5)):
                    for ix, dx in enumerate(range(-4, 5)):
                        region[iy, ix] = max(abs(dy), abs(dx)) <= radius
                target_factor = float(tkernel[region].sum())
                control_factor = float(ckernel[region].sum())
                row[f"target_factor_r{radius}"] = target_factor
                row[f"control_factor_r{radius}"] = control_factor
                row[f"paired_delta_r{radius}"] = target_factor - control_factor
                if radius > 0:
                    row[f"target_shell_r{radius}"] = target_factor - float(tkernel[
                        np.fromfunction(lambda y, x: np.maximum(np.abs(y - 4), np.abs(x - 4)) <= radius - 1,
                                        (9, 9), dtype=int).astype(bool)].sum())
                    row[f"control_shell_r{radius}"] = control_factor - float(ckernel[
                        np.fromfunction(lambda y, x: np.maximum(np.abs(y - 4), np.abs(x - 4)) <= radius - 1,
                                        (9, 9), dtype=int).astype(bool)].sum())
            per_pair.append(row)

        _, labels, n_target_cells, n_control_cells, n_blocks, largest = pair_block_labels(
            target_positions, control_positions, np.ones(len(ids), dtype=bool),
        )
        labels = np.asarray(labels, dtype=str)
        metrics = {}
        for radius in range(LAG_RADIUS + 1):
            target = np.asarray([row[f"target_factor_r{radius}"] for row in per_pair])
            control = np.asarray([row[f"control_factor_r{radius}"] for row in per_pair])
            metrics[f"square_lag_radius_{radius}_pixels"] = connected_block_bootstrap(
                target, control, labels, replicas=replicas,
                seed=seed + 10000 * filter_index + 100 * radius,
            )
        shells = {}
        for radius in range(1, LAG_RADIUS + 1):
            target = np.asarray([row[f"target_shell_r{radius}"] for row in per_pair])
            control = np.asarray([row[f"control_shell_r{radius}"] for row in per_pair])
            shells[f"shell_{radius}_pixels"] = connected_block_bootstrap(
                target, control, labels, replicas=replicas,
                seed=seed + 50000 + 10000 * filter_index + 100 * radius,
            )
        by_filter[band] = {
            "n_pairs": len(ids), "n_target_sky_cells": n_target_cells,
            "n_control_sky_cells": n_control_cells, "n_connected_5deg_blocks": n_blocks,
            "largest_block_pairs": largest, "cumulative_square_kernel_factors": metrics,
            "incremental_shell_factors": shells,
        }
        pair_rows.extend(per_pair)
        for radius in range(LAG_RADIUS + 1):
            lag_rows.append({
                "filter": band, "max_abs_lag_pixels": radius,
                **metrics[f"square_lag_radius_{radius}_pixels"],
            })

    input_paths = [audit_dir / "summary.json", audit_dir / "independent_audit.json",
                   lag_path, sample / "download_01/download_manifest.json",
                   sample / "control_download_01/download_manifest.json",
                   sample / "download_01/download_status.json",
                   sample / "control_download_01/download_status.json"]
    result = {
        "schema": "astrafek-icl-coadd-lag-kernel-truncation-v1",
        "status": "complete_exploratory_diagnostic",
        "science_measurement": False,
        "hsc_network_used": False,
        "images_read": False,
        "data_requested": False,
        "method": "sum symmetric per-frame median normalized residual-covariance kernels inside nested square lag cutoffs 0..4 pixels; compare matched target/control frames and bootstrap connected 5-degree sky blocks within filter",
        "n_pairs": sum(EXPECTED_PAIRS.values()),
        "filters": by_filter,
        "limitations": [
            "The source kernels were measured from a finite set of non-overlapping 64-pixel patches and each lag value is a per-frame median, not an aperture-specific covariance field.",
            "The existing input ends at 4 pixels; this audit tests truncation within the available 9x9 kernel and cannot bound longer-range correlation.",
            "Summed square-kernel factors approximate stationary top-hat noise and are not exact covariance of the masked annular ICL estimator.",
            "Source Poisson noise, sky subtraction transfer, PSF wings and full field covariance are not included; keep science_measurement=false.",
        ],
        "input_sha256": {str(path): sha256(path) for path in input_paths},
        "upstream_covariance_summary_sha256": sha256(audit_dir / "summary.json"),
        "upstream_independent_audit_sha256": sha256(audit_dir / "independent_audit.json"),
        "analysis_script_sha256": sha256(Path(__file__).resolve()),
    }
    output.mkdir(parents=True, exist_ok=False)
    write_tsv(output / "per_pair_kernel_factors.tsv", pair_rows)
    write_tsv(output / "lag_cutoff_summary.tsv", lag_rows)
    summary_path = output / "summary.json"
    summary_path.write_text(json.dumps(result, ensure_ascii=False, indent=2,
                                      allow_nan=False, sort_keys=True) + "\n",
                            encoding="utf-8")
    result["output_sha256"] = {
        path.name: sha256(path) for path in
        (summary_path, output / "per_pair_kernel_factors.tsv", output / "lag_cutoff_summary.tsv")
    }
    print(json.dumps({
        "status": result["status"], "science_measurement": False,
        "filters": {band: {
            "pairs": section["n_pairs"],
            "blocks": section["n_connected_5deg_blocks"],
            "cumulative_factors": {
                radius: values["target_median"] for radius, values in
                ((str(i), section["cumulative_square_kernel_factors"][f"square_lag_radius_{i}_pixels"])
                 for i in range(LAG_RADIUS + 1))
            },
            "paired_differences": {
                radius: values["paired_target_minus_control_median"] for radius, values in
                ((str(i), section["cumulative_square_kernel_factors"][f"square_lag_radius_{i}_pixels"])
                 for i in range(LAG_RADIUS + 1))
            },
        } for band, section in by_filter.items()},
        "output": str(output), "output_sha256": result["output_sha256"],
    }, indent=2, ensure_ascii=False))
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", type=Path, default=DEFAULT_SAMPLE)
    parser.add_argument("--audit", type=Path, default=DEFAULT_AUDIT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--replicas", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=20261106)
    args = parser.parse_args()
    if args.replicas < 1000:
        parser.error("use at least 1,000 spatial-block bootstrap replicas")
    analyze(args.sample, args.audit, args.output, args.replicas, args.seed)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
