#!/usr/bin/env python3
"""Estimate connected-sky-block uncertainty in paired ICL transfer vs scale."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

from audit_icl_psf_core_dependence import connected_block_labels, load_verified_inputs
from summarize_icl_psf_diffuse_scale_grid import load_run, sha256


DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
DEFAULT_GRID = DATA / "icl_psf_diffuse_scale_grid_2026-10-20_01/summary.json"
DEFAULT_ICL = DATA / "icl_annulus_support_2026-10-18_01"
DEFAULT_EXTRACTION = DATA / "prestack_extraction_2026-10-15_01/extraction_summary.json"
DEFAULT_PSFS = DATA / "psf_core_variation_2026-10-20_01/summary.json"
SCALES = (100.0, 200.0, 300.0, 500.0)


def summary(values: np.ndarray) -> dict:
    return {"mean": float(np.mean(values)), "median": float(np.median(values)),
            "q16": float(np.quantile(values, 0.16)),
            "q84": float(np.quantile(values, 0.84))}


def connected_block_bootstrap(retention: np.ndarray, loss: np.ndarray,
                              labels: np.ndarray, replicas: int, seed: int) -> dict:
    retention, loss = np.asarray(retention, dtype=float), np.asarray(loss, dtype=float)
    labels = np.asarray(labels, dtype=str)
    if (retention.ndim != 2 or loss.shape != retention.shape
            or len(retention) != len(labels) or retention.shape[1] < 2):
        raise ValueError("retention/loss must be aligned rows-by-scales with one block label per row")
    if not np.isfinite(retention).all() or not np.isfinite(loss).all():
        raise ValueError("bootstrap values must all be finite")
    blocks = np.unique(labels)
    if len(blocks) < 6:
        raise ValueError(f"at least six connected blocks are required, found {len(blocks)}")
    indices_by_block = [np.flatnonzero(labels == block) for block in blocks]
    rng = np.random.default_rng(seed)
    n_scales = retention.shape[1]
    median_retention = np.empty((replicas, n_scales), dtype=float)
    mean_retention = np.empty((replicas, n_scales), dtype=float)
    median_loss = np.empty((replicas, n_scales), dtype=float)
    median_paired_delta = np.empty((replicas, n_scales - 1), dtype=float)
    mean_paired_delta = np.empty((replicas, n_scales - 1), dtype=float)
    median_paired_loss_delta = np.empty((replicas, n_scales - 1), dtype=float)
    for replica in range(replicas):
        chosen = rng.integers(0, len(blocks), size=len(blocks))
        indices = np.concatenate([indices_by_block[int(index)] for index in chosen])
        sample_retention = retention[indices]
        sample_loss = loss[indices]
        median_retention[replica] = np.median(sample_retention, axis=0)
        mean_retention[replica] = sample_retention.mean(axis=0)
        median_loss[replica] = np.median(sample_loss, axis=0)
        median_paired_delta[replica] = np.median(
            sample_retention[:, 1:] - sample_retention[:, [0]], axis=0)
        mean_paired_delta[replica] = np.mean(
            sample_retention[:, 1:] - sample_retention[:, [0]], axis=0)
        median_paired_loss_delta[replica] = np.median(
            sample_loss[:, 1:] - sample_loss[:, [0]], axis=0)
    point_median = np.median(retention, axis=0)
    point_mean = np.mean(retention, axis=0)
    point_loss_median = np.median(loss, axis=0)
    point_median_delta = np.median(retention[:, 1:] - retention[:, [0]], axis=0)
    point_mean_delta = np.mean(retention[:, 1:] - retention[:, [0]], axis=0)
    point_loss_delta = np.median(loss[:, 1:] - loss[:, [0]], axis=0)
    return {
        "n_connected_blocks": int(len(blocks)),
        "largest_block_rows": int(max(map(len, indices_by_block))),
        "replicas_requested": int(replicas),
        "replicas_valid": int(len(median_retention)),
        "seed": int(seed),
        "median_retention_by_scale": {
            str(int(scale)): {"point": float(point_median[index]),
                              "bootstrap_95_interval": np.quantile(median_retention[:, index], [0.025, 0.975]).tolist()}
            for index, scale in enumerate(SCALES)
        },
        "mean_retention_sensitivity_by_scale": {
            str(int(scale)): {"point": float(point_mean[index]),
                              "bootstrap_95_interval": np.quantile(mean_retention[:, index], [0.025, 0.975]).tolist()}
            for index, scale in enumerate(SCALES)
        },
        "paired_median_retention_change_from_100kpc": {
            str(int(scale)): {"point": float(point_median_delta[index - 1]),
                              "bootstrap_95_interval": np.quantile(median_paired_delta[:, index - 1], [0.025, 0.975]).tolist()}
            for index, scale in enumerate(SCALES[1:], start=1)
        },
        "paired_mean_retention_change_sensitivity_from_100kpc": {
            str(int(scale)): {"point": float(point_mean_delta[index - 1]),
                              "bootstrap_95_interval": np.quantile(mean_paired_delta[:, index - 1], [0.025, 0.975]).tolist()}
            for index, scale in enumerate(SCALES[1:], start=1)
        },
        "median_plane_transfer_loss_by_scale": {
            str(int(scale)): {"point": float(point_loss_median[index]),
                              "bootstrap_95_interval": np.quantile(median_loss[:, index], [0.025, 0.975]).tolist()}
            for index, scale in enumerate(SCALES)
        },
        "paired_median_plane_loss_change_from_100kpc": {
            str(int(scale)): {"point": float(point_loss_delta[index - 1]),
                              "bootstrap_95_interval": np.quantile(median_paired_loss_delta[:, index - 1], [0.025, 0.975]).tolist()}
            for index, scale in enumerate(SCALES[1:], start=1)
        },
    }


def analyze(grid_path: Path, icl_dir: Path, extraction_path: Path,
            psf_summary_path: Path, replicas: int = 5000,
            seed: int = 20261020) -> dict:
    grid_path = Path(grid_path).resolve()
    grid = json.loads(grid_path.read_text(encoding="utf-8"))
    if (grid.get("status") != "verified_exploratory_scale_grid"
            or grid.get("science_measurement") is not False
            or grid.get("scale_lengths_kpc") != list(SCALES)):
        raise ValueError("scale-grid summary is incomplete or has unexpected scales")

    runs = {}
    for scale in SCALES:
        record = grid["runs"][str(int(scale))]
        directory = Path(record["directory"])
        run = load_run(scale, directory)
        if (run["summary_sha256"] != record["summary_sha256"]
                or run["audit_sha256"] != record["audit_sha256"]):
            raise ValueError(f"scale-grid provenance mismatch at {scale:g} kpc")
        runs[scale] = run

    # This join independently verifies the approved 57-frame PSF sample and
    # provides the matched target/control sky coordinates for block labels.
    spatial_rows, spatial_provenance = load_verified_inputs(
        icl_dir, extraction_path, Path(grid["runs"]["200"]["directory"]), psf_summary_path)
    by_filter = {}
    for filter_index, band in enumerate(("i", "i2")):
        selected = [row for row in spatial_rows if row["filter"] == band]
        expected = 13 if band == "i" else 44
        if len(selected) != expected:
            raise ValueError(f"fixed {band} PSF support changed: {len(selected)} != {expected}")
        keys = [(str(row["cluster_id"]), band) for row in selected]
        if len(set(keys)) != len(keys):
            raise ValueError(f"duplicate cluster/filter row in spatial support: {band}")
        scale_rows = {scale: runs[scale]["records"] for scale in SCALES}
        if any(any(key not in scale_rows[scale] for key in keys) for scale in SCALES):
            raise ValueError(f"spatial support does not match scale runs for {band}")

        raw = np.asarray([
            [float(scale_rows[scale][key]["transfer_response"][
                "icl_response_to_diffuse_unit_injection"]) for scale in SCALES]
            for key in keys
        ], dtype=float)
        corrected = np.asarray([
            [float(scale_rows[scale][key]["transfer_response"][
                "icl_plane_corrected_response_to_diffuse_unit_injection"]) for scale in SCALES]
            for key in keys
        ], dtype=float)
        positive = np.isfinite(raw).all(axis=1) & np.isfinite(corrected).all(axis=1) & (raw > 0).all(axis=1)
        if not positive.any():
            raise ValueError(f"no frames have positive raw response across every scale for {band}")
        retention = corrected[positive] / raw[positive]
        transfer_loss = corrected[positive] - raw[positive]
        retained_rows = [row for row, keep in zip(selected, positive) if keep]
        labels = connected_block_labels(retained_rows)
        block_ids, block_counts = np.unique(labels, return_counts=True)
        if len(block_ids) < 6:
            raise ValueError(f"positive-raw support has too few connected blocks for {band}: {len(block_ids)}")
        result = connected_block_bootstrap(retention, transfer_loss, labels,
                                           replicas, seed + 10000 * filter_index)
        raw_minima = {}
        for index, scale in enumerate(SCALES):
            selected_raw = raw[positive, index]
            row_index = int(np.argmin(selected_raw))
            raw_minima[str(int(scale))] = {
                "minimum_positive_raw_response": float(selected_raw[row_index]),
                "cluster_id": retained_rows[row_index]["cluster_id"],
                "retention_ratio_at_minimum": float(retention[row_index, index]),
            }
        result.update({
            "n_frames_before_positive_raw_selection": int(len(selected)),
            "n_frames_used": int(int(positive.sum())),
            "n_frames_excluded_for_nonpositive_or_invalid_raw_response": int((~positive).sum()),
            "descriptive_frame_scatter_by_scale": {
                str(int(scale)): summary(retention[:, index])
                for index, scale in enumerate(SCALES)
            },
            "paired_change_frame_scatter_from_100kpc": {
                str(int(scale)): summary(retention[:, index] - retention[:, 0])
                for index, scale in enumerate(SCALES[1:], start=1)
            },
            "minimum_positive_raw_response_diagnostics": raw_minima,
            "median_absolute_plane_transfer_loss_by_scale": {
                str(int(scale)): summary(transfer_loss[:, index])
                for index, scale in enumerate(SCALES)
            },
        })
        by_filter[band] = result

    return {
        "schema": "astrafek-icl-psf-diffuse-scale-connected-block-bootstrap-v1",
        "status": "complete_exploratory_spatial_block_bootstrap",
        "science_measurement": False,
        "estimand": "median per-frame plane-corrected/raw diffuse ICL response ratio and median absolute corrected-minus-raw response among frames with positive raw response at every scale",
        "scales_kpc": list(SCALES),
        "method": "paired same-frame connected target/control 5-degree sky-block bootstrap; the same selected blocks resampled jointly across all four profile scales",
        "interval_note": "percentile 95% spatial-block bootstrap intervals for paired medians; means are secondary sensitivity estimates and no p-values are computed",
        "filters": by_filter,
        "inputs": {
            "scale_grid_summary": {"path": str(grid_path), "sha256": sha256(grid_path)},
            "icl_support_summary_sha256": spatial_provenance["icl_summary_sha256"],
            "icl_annulus_table_sha256": spatial_provenance["icl_annulus_table_sha256"],
            "matched_control_extraction_sha256": spatial_provenance["extraction_summary_sha256"],
            "position_psf_batch_summary_sha256": spatial_provenance["joint_summary_sha256"],
            "position_psf_audit_sha256": spatial_provenance["joint_audit_sha256"],
            "psf_core_summary_sha256": spatial_provenance["psf_core_summary_sha256"],
        },
        "analysis_script_sha256": sha256(Path(__file__).resolve()),
        "limitations": [
            "The inference is conditional on the selected 57-frame support and positive raw responses across all scales.",
            "A very small positive raw response can make the mean retention ratio unstable; minimum denominators are recorded and the paired median is primary.",
            "The 13-frame i-band interval is based on 10 blocks before response filtering and may remain imprecise.",
            "Connected 5-degree blocks are a spatial dependence approximation; pixel/coadd covariance is not included here.",
            "This is an injection-transfer sensitivity, not an observed ICL correction or detection.",
            "PSF wings, HSC upstream sky processing, photon noise, extinction and passband uncertainty remain unresolved.",
        ],
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--grid-summary", type=Path, default=DEFAULT_GRID)
    parser.add_argument("--icl-dir", type=Path, default=DEFAULT_ICL)
    parser.add_argument("--extraction-summary", type=Path, default=DEFAULT_EXTRACTION)
    parser.add_argument("--psf-core-summary", type=Path, default=DEFAULT_PSFS)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--replicas", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=20261020)
    args = parser.parse_args()
    if args.replicas < 1000:
        parser.error("at least 1000 spatial-block bootstrap replicates are required")
    if args.output.exists():
        parser.error(f"refusing to overwrite existing output: {args.output}")
    result = analyze(args.grid_summary, args.icl_dir, args.extraction_summary,
                     args.psf_core_summary, args.replicas, args.seed)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n",
                           encoding="utf-8")
    print(json.dumps({"status": result["status"], "science_measurement": False,
                      "filters": result["filters"], "output": str(args.output)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
