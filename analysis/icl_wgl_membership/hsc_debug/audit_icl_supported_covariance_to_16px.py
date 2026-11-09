#!/usr/bin/env python3
"""Measure longer residual-covariance lags for the supported ICL pairs.

This focuses the existing HSC coadd noise-patch estimator on the 224
target/control pairs that contribute to the 200--300 proper-kpc ICL contrast.
It extends the local normalized lag kernel to +/-16 pixels using only
already acquired images, archived target masks, and existing control masks.
Per-frame checkpoints allow a Lageunha CPU run to resume without overwriting
or losing completed work. Results remain diagnostic and science-gated false.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import time

import numpy as np

HERE = Path(__file__).resolve().parent
import sys
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from audit_icl_full_sample_coadd_covariance import (
    DEFAULT_SAMPLE, build_records, frame_correlations, read_json, sha256, write_json,
)
from analyze_icl_annulus_aperture_covariance import load_profiles
from evaluate_prestack_pair_support import pair_block_labels

DEFAULT_OUTPUT = DEFAULT_SAMPLE / "icl_supported_covariance_to_16px_2026-11-06_01"
AUDIT_DIR = DEFAULT_SAMPLE / "icl_coadd_noise_correlation_full_2026-11-03_01"
PATCH_SIZE = 64
MAX_PATCHES = 32
MAX_LAG = 16
MIN_VALID_FRACTION = 0.80
BOOTSTRAP_REPLICAS = 5000
BOOTSTRAP_SEED = 20261106


def write_tsv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty table: {path}")
    with path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t",
                                lineterminator="\n", extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def block_summary(target: np.ndarray, control: np.ndarray,
                  labels: np.ndarray, replicas: int, seed: int) -> dict:
    groups = np.unique(labels)
    members = [np.flatnonzero(labels == group) for group in groups]
    if len(members) < 3 or sum(map(len, members)) != len(target):
        raise ValueError("invalid filter-specific connected-block assignment")
    delta = target - control
    rng = np.random.default_rng(seed)
    boot_target = np.empty(replicas)
    boot_control = np.empty(replicas)
    boot_delta = np.empty(replicas)
    for index in range(replicas):
        draw = rng.integers(0, len(members), size=len(members))
        selected = np.concatenate([members[int(j)] for j in draw])
        boot_target[index] = np.median(target[selected])
        boot_control[index] = np.median(control[selected])
        boot_delta[index] = np.median(delta[selected])
    return {
        "n_connected_blocks": len(members),
        "largest_block_pairs": max(map(len, members)),
        "target_median_factor": float(np.median(target)),
        "control_median_factor": float(np.median(control)),
        "paired_median_target_minus_control": float(np.median(delta)),
        "paired_delta_ci95": [float(np.quantile(boot_delta, .025)),
                               float(np.quantile(boot_delta, .975))],
        "difference_of_medians_ci95": [float(np.quantile(boot_target-boot_control, .025)),
                                       float(np.quantile(boot_target-boot_control, .975))],
        "replicas": replicas, "seed": seed,
    }


def symmetric_kernel(rows: list[dict]) -> np.ndarray:
    raw = {(int(row["dy"]), int(row["dx"])):
           float(row["median_variance_normalized_covariance"]) for row in rows}
    expected = {(dy, dx) for dy in range(-MAX_LAG, MAX_LAG+1)
                for dx in range(-MAX_LAG, MAX_LAG+1)}
    if set(raw) != expected:
        raise ValueError("frame checkpoint does not contain the complete 33x33 lag kernel")
    kernel = np.empty((2*MAX_LAG+1, 2*MAX_LAG+1), dtype=float)
    for dy, dx in expected:
        kernel[dy+MAX_LAG, dx+MAX_LAG] = 0.5 * (raw[(dy, dx)] + raw[(-dy, -dx)])
    if not np.all(np.isfinite(kernel)):
        raise ValueError("non-finite kernel values")
    return kernel


def run(output: Path, *, resume: bool = False, replicas: int = BOOTSTRAP_REPLICAS) -> dict:
    output = Path(output).resolve()
    profiles, profile_info = load_profiles()
    records, setup = build_records(DEFAULT_SAMPLE.resolve())
    pair_by_id = {str(row["cluster_id"]): row for row in setup["pair_rows"]}
    supported = set(profiles)
    selected = [row for row in records
                if (str(row["cluster_id"]), str(row["filter"])) in supported]
    expected_keys = {(ident, band, role) for ident, band in supported
                     for role in ("target", "matched_control")}
    actual_keys = {(str(r["cluster_id"]), str(r["filter"]), str(r["role"])) for r in selected}
    if len(supported) != 224 or actual_keys != expected_keys or len(selected) != 448:
        raise ValueError("supported ICL image list is not exactly 224 paired target/control frames")
    selected.sort(key=lambda row: (row["filter"], str(row["cluster_id"]), row["role"]))
    input_hashes = {
        str(path): sha256(path) for path in (
            AUDIT_DIR / "summary.json", AUDIT_DIR / "independent_audit.json",
            AUDIT_DIR / "lag_correlation.tsv", *profile_info["paths"],
        )
    }
    params = {
        "schema": "astrafek-icl-long-lag-run-v1", "supported_pairs": 224,
        "patch_size": PATCH_SIZE, "max_patches": MAX_PATCHES, "max_lag": MAX_LAG,
        "min_valid_fraction": MIN_VALID_FRACTION,
        "input_sha256": input_hashes,
    }
    fingerprint = hashlib.sha256(json.dumps(params, sort_keys=True,
                                             separators=(",", ":")).encode()).hexdigest()
    if output.exists():
        if not resume:
            raise FileExistsError(f"refusing existing output without --resume: {output}")
        state_path = output / "progress.json"
        if not state_path.is_file():
            raise ValueError("resume directory has no progress checkpoint")
        state = read_json(state_path)
        if state.get("run_fingerprint") != fingerprint:
            raise ValueError("resume checkpoint input hashes or parameters differ")
    else:
        if resume:
            raise FileNotFoundError("--resume requested but output directory does not exist")
        output.mkdir(parents=True, exist_ok=False)
        state = {"schema": "astrafek-icl-long-lag-checkpoint-v1", "status": "running",
                 "run_fingerprint": fingerprint, "params": params,
                 "completed_frames": {}, "science_measurement": False}
        write_json(output / "progress.json", state)

    for index, record in enumerate(selected, 1):
        frame_key = f"{record['filter']}:{record['cluster_id']}:{record['role']}"
        checkpoint_path = output / "frames" / f"{index:04d}.json"
        old = state["completed_frames"].get(frame_key)
        if old:
            if not checkpoint_path.is_file() or sha256(checkpoint_path) != old.get("sha256"):
                raise ValueError(f"completed frame checkpoint is missing or changed: {frame_key}")
            continue
        checkpoint_path.parent.mkdir(parents=True, exist_ok=True)
        print(f"[{index}/{len(selected)}] long-lag covariance {frame_key}", flush=True)
        detail, lag_rows = frame_correlations(
            record, patch_size=PATCH_SIZE, max_patches=MAX_PATCHES,
            max_lag=MAX_LAG, min_valid_fraction=MIN_VALID_FRACTION,
        )
        detail["connected_5deg_block"] = record["connected_5deg_block"]
        payload = {"detail": detail, "lag_rows": lag_rows}
        temporary = checkpoint_path.with_suffix(".json.tmp")
        temporary.write_text(json.dumps(payload, ensure_ascii=False, allow_nan=False) + "\n",
                             encoding="utf-8")
        temporary.replace(checkpoint_path)
        state["completed_frames"][frame_key] = {
            "checkpoint": str(checkpoint_path), "sha256": sha256(checkpoint_path),
        }
        state["n_frames_complete"] = len(state["completed_frames"])
        state["last_frame"] = frame_key
        write_json(output / "progress.json", state)

    frame_details, all_lags = [], []
    frame_by_key = {}
    for index, record in enumerate(selected, 1):
        frame_key = f"{record['filter']}:{record['cluster_id']}:{record['role']}"
        checkpoint = output / "frames" / f"{index:04d}.json"
        payload = read_json(checkpoint)
        detail, lag_rows = payload["detail"], payload["lag_rows"]
        if ((detail.get("cluster_id"), detail.get("filter"), detail.get("role"))
                != (str(record["cluster_id"]), str(record["filter"]), str(record["role"]))):
            raise ValueError(f"frame checkpoint identity mismatch: {frame_key}")
        frame_details.append(detail)
        all_lags.extend(lag_rows)
        frame_by_key[(str(detail["cluster_id"]), str(detail["filter"]), str(detail["role"]))] = detail

    by_filter = {}
    for filter_index, band in enumerate(("i", "i2")):
        ids = sorted(ident for ident, filt in supported if filt == band)
        pair_target, pair_control = [], []
        for ident in ids:
            target = frame_by_key[(ident, band, "target")]
            control = frame_by_key[(ident, band, "matched_control")]
            pair_target.append(target); pair_control.append(control)
        target_positions = [{"ra_deg": float(pair_by_id[ident]["target_ra_deg"]),
                             "dec_deg": float(pair_by_id[ident]["target_dec_deg"])}
                            for ident in ids]
        control_positions = [{"ra_deg": float(pair_by_id[ident]["control_ra_deg"]),
                              "dec_deg": float(pair_by_id[ident]["control_dec_deg"])}
                             for ident in ids]
        _, labels, n_target_cells, n_control_cells, n_blocks, largest = pair_block_labels(
            target_positions, control_positions, np.ones(len(ids), dtype=bool),
        )
        labels = np.asarray(labels, dtype=str)
        target_kernels = []
        control_kernels = []
        for target, control in zip(pair_target, pair_control):
            trows = [row for row in all_lags if row["cluster_id"] == target["cluster_id"]
                     and row["filter"] == band and row.get("role") == "target"]
            crows = [row for row in all_lags if row["cluster_id"] == control["cluster_id"]
                     and row["filter"] == band and row.get("role") == "matched_control"]
            target_kernels.append(symmetric_kernel(trows))
            control_kernels.append(symmetric_kernel(crows))
        target_kernels = np.stack(target_kernels)
        control_kernels = np.stack(control_kernels)
        cumulative, incremental = {}, {}
        previous_target = np.zeros(len(ids)); previous_control = np.zeros(len(ids))
        for radius in range(MAX_LAG + 1):
            region = np.zeros((2*MAX_LAG+1, 2*MAX_LAG+1), dtype=bool)
            for iy, dy in enumerate(range(-MAX_LAG, MAX_LAG+1)):
                for ix, dx in enumerate(range(-MAX_LAG, MAX_LAG+1)):
                    region[iy, ix] = max(abs(dy), abs(dx)) <= radius
            target = target_kernels[:, region].sum(axis=1)
            control = control_kernels[:, region].sum(axis=1)
            cumulative[f"square_lag_radius_{radius}_pixels"] = block_summary(
                target, control, labels, replicas,
                BOOTSTRAP_SEED + 10000*filter_index + 100*radius,
            )
            if radius:
                incremental[f"shell_{radius}_pixels"] = block_summary(
                    target-previous_target, control-previous_control, labels, replicas,
                    BOOTSTRAP_SEED + 50000 + 10000*filter_index + 100*radius,
                )
            previous_target, previous_control = target, control
        by_filter[band] = {
            "n_pairs": len(ids), "n_target_sky_cells": n_target_cells,
            "n_control_sky_cells": n_control_cells, "n_connected_5deg_blocks": n_blocks,
            "largest_block_pairs": largest, "cumulative_square_kernel_factors": cumulative,
            "incremental_shell_factors": incremental,
        }

    frame_path, lag_path = output / "frame_noise_correlation.tsv", output / "lag_correlation.tsv"
    write_tsv(frame_path, frame_details)
    write_tsv(lag_path, all_lags)
    summary = {
        "schema": "astrafek-icl-supported-long-lag-noise-v1",
        "status": "complete_exploratory_diagnostic", "science_measurement": False,
        "hsc_network_used": False, "images_read": True, "data_requested": False,
        "sample": {"supported_200_300_kpc_pairs": len(supported), "image_frames": len(selected),
                   "filter_pairs": {band: sum(ident_band == band for _, ident_band in supported)
                                    for band in ("i", "i2")}},
        "method": {
            "patch_size": PATCH_SIZE, "max_patches_per_frame": MAX_PATCHES,
            "max_lag_pixels": MAX_LAG, "min_valid_fraction": MIN_VALID_FRACTION,
            "lag_kernel": "median variance-normalized residual covariance from 32 deterministic non-overlapping 64x64 patches per image; affine plane removed, 4-sigma clipping",
            "aperture_relevance": "only 224 exact pairs that support the saved 200-300 proper-kpc ICL contrast",
        },
        "filters": by_filter,
        "input_sha256": input_hashes,
        "run_fingerprint": fingerprint,
        "analysis_script_sha256": sha256(Path(__file__).resolve()),
        "limitations": [
            "The 33x33 kernel is still truncated at 16 pixels and its lag values are patch medians rather than a spatially varying covariance map.",
            "The 32-patch per-image estimate is less precise than the earlier 96-patch 9x9 measurement.",
            "Square stationary kernel sums approximate aperture noise; exact radial mask covariance, source Poisson, sky transfer and PSF wings remain unresolved.",
            "5-degree blocks are an approximate footprint partition. Keep science_measurement=false.",
        ],
    }
    summary_path = output / "summary.json"
    summary_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2,
                                      allow_nan=False, sort_keys=True) + "\n",
                            encoding="utf-8")
    hashes = {path.name: sha256(path) for path in (summary_path, frame_path, lag_path)}
    state.update(status="complete_exploratory_diagnostic", n_frames_complete=len(selected),
                 output_sha256=hashes, summary_sha256=hashes["summary.json"])
    write_json(output / "progress.json", state)
    print(json.dumps({"status": summary["status"], "science_measurement": False,
                      "filters": by_filter, "output": str(output),
                      "output_sha256": hashes}, ensure_ascii=False, indent=2), flush=True)
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--replicas", type=int, default=BOOTSTRAP_REPLICAS)
    args = parser.parse_args()
    if args.replicas < 1000:
        parser.error("at least 1000 spatial-block bootstrap replicas are required")
    run(args.output, resume=args.resume, replicas=args.replicas)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
