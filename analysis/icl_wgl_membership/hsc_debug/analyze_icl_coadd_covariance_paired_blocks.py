#!/usr/bin/env python3
"""Block-bootstrap target/control differences in measured coadd noise covariance.

This uses only the completed full-sample local-lag audit and the existing
target/control acquisition manifests. It is a diagnostic comparison, not an
ICL pixel-estimator covariance correction or a science measurement.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from evaluate_prestack_pair_support import pair_block_labels


DEFAULT_SAMPLE = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
DEFAULT_AUDIT = DEFAULT_SAMPLE / "icl_coadd_noise_correlation_full_2026-11-03_01"
DEFAULT_OUTPUT = DEFAULT_SAMPLE / "icl_coadd_covariance_paired_blocks_2026-11-03_01"
METRICS = {
    "variance_factor": "top_hat_variance_inflation_approx_square_lag_window",
    "sigma_factor": "top_hat_sigma_inflation_approx_square_lag_window",
    "zero_lag_ratio": "variance_plane_ratio_at_zero_lag",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def finite(row: dict[str, str], key: str) -> float:
    value = float(row[key])
    if not math.isfinite(value):
        raise ValueError(f"non-finite {key}")
    return value


def audit_inputs(sample: Path, audit_dir: Path) -> tuple[dict, list[dict[str, str]], dict, dict]:
    summary_path = audit_dir / "summary.json"
    verification_path = audit_dir / "independent_audit.json"
    frame_path = audit_dir / "frame_noise_correlation.tsv"
    required = (summary_path, verification_path, frame_path,
                sample / "download_01/download_manifest.json",
                sample / "control_download_01/download_manifest.json",
                sample / "download_01/download_status.json",
                sample / "control_download_01/download_status.json")
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise FileNotFoundError(", ".join(missing))
    summary, verification = read_json(summary_path), read_json(verification_path)
    if (summary.get("status") != "complete_exploratory"
            or summary.get("science_measurement") is not False
            or summary.get("sample", {}).get("same_filter_pairs") != 783
            or verification.get("status") != "verified"
            or verification.get("science_measurement") is not False
            or verification.get("input_acquisition_hashes_matched") is not True
            or verification.get("frame_and_lag_arithmetic_reproduced") is not True
            or verification.get("n_frame_rows") != 1566):
        raise ValueError("full-sample covariance audit is incomplete or unexpected")
    if sha256(summary_path) != verification.get("summary_sha256"):
        raise ValueError("independently verified summary hash mismatch")
    manifest = read_json(required[3])
    control_manifest = read_json(required[4])
    target_status, control_status = read_json(required[5]), read_json(required[6])
    if (target_status.get("status") != "complete" or target_status.get("downloaded") != 800
            or control_status.get("status") != "complete" or control_status.get("downloaded") != 798):
        raise ValueError("target/control acquisition records are incomplete")
    for path in required[3:]:
        expected = summary.get("inputs", {}).get(str(path))
        if expected is not None and sha256(path) != expected:
            raise ValueError(f"acquisition record hash changed: {path}")
    return summary, read_tsv(frame_path), manifest, control_manifest


def analyze(sample: Path, audit_dir: Path, output: Path, replicas: int, seed: int) -> dict:
    sample, audit_dir, output = sample.resolve(), audit_dir.resolve(), output.resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite {output}")
    if replicas < 1000:
        raise ValueError("replicas must be at least 1000")
    summary, frames, manifest, control_manifest = audit_inputs(sample, audit_dir)
    target_by_id = {str(row["id"]): row for row in manifest["clusters"]}
    control_by_target = {str(row["target_cluster_id"]): row for row in control_manifest["clusters"]}
    frame_by_key = {}
    for row in frames:
        key = (str(row["cluster_id"]), str(row["role"]), str(row["filter"]))
        if key in frame_by_key:
            raise ValueError(f"duplicate frame record: {key}")
        frame_by_key[key] = row

    filters = {}
    for band, expected_n in (("i", 302), ("i2", 481)):
        target_ids = sorted(ident for ident, role, filt in frame_by_key
                            if role == "target" and filt == band)
        if len(target_ids) != expected_n:
            raise ValueError(f"unexpected {band}-band target support: {len(target_ids)}")
        target_rows, control_rows = [], []
        paired_frames = []
        for ident in target_ids:
            target = target_by_id.get(ident)
            control = control_by_target.get(ident)
            trow = frame_by_key.get((ident, "target", band))
            crow = frame_by_key.get((ident, "matched_control", band))
            if target is None or control is None or trow is None or crow is None:
                raise ValueError(f"incomplete pair for {ident}/{band}")
            target_band = (target.get("coadd_header_filter") or target.get("filter")
                           or target.get("band"))
            if target_band and target_band != band:
                raise ValueError(f"target manifest/filter mismatch for {ident}")
            if control.get("coadd_header_filter") != band:
                raise ValueError(f"control manifest/filter mismatch for {ident}")
            target_rows.append({"ra_deg": float(target["ra_deg"]),
                                "dec_deg": float(target["dec_deg"])})
            control_rows.append({"ra_deg": float(control["ra_deg"]),
                                 "dec_deg": float(control["dec_deg"])})
            paired_frames.append((trow, crow))

        _, labels, n_target_cells, n_control_cells, n_blocks, largest = pair_block_labels(
            target_rows, control_rows, np.ones(expected_n, dtype=bool))
        label_values = np.asarray(labels, dtype=str)
        groups = np.unique(label_values)
        members = [np.flatnonzero(label_values == group) for group in groups]
        if len(members) != n_blocks or sum(map(len, members)) != expected_n:
            raise ValueError(f"filter-specific block partition is incomplete for {band}")
        rng = np.random.default_rng(seed + (0 if band == "i" else 10000))
        result_by_metric = {}
        for metric, field in METRICS.items():
            target = np.asarray([finite(pair[0], field) for pair in paired_frames])
            control = np.asarray([finite(pair[1], field) for pair in paired_frames])
            delta = target - control
            boot_target = np.empty(replicas, dtype=float)
            boot_control = np.empty(replicas, dtype=float)
            boot_paired_delta = np.empty(replicas, dtype=float)
            for replica in range(replicas):
                draw = rng.integers(0, len(members), size=len(members))
                indices = np.concatenate([members[int(index)] for index in draw])
                boot_target[replica] = np.median(target[indices])
                boot_control[replica] = np.median(control[indices])
                boot_paired_delta[replica] = np.median(delta[indices])
            result_by_metric[metric] = {
                "target_median": float(np.median(target)),
                "control_median": float(np.median(control)),
                "paired_target_minus_control_median": float(np.median(delta)),
                "paired_delta_block_bootstrap_ci95": [
                    float(np.quantile(boot_paired_delta, .025)),
                    float(np.quantile(boot_paired_delta, .975)),
                ],
                "difference_of_medians_block_bootstrap_ci95": [
                    float(np.quantile(boot_target - boot_control, .025)),
                    float(np.quantile(boot_target - boot_control, .975)),
                ],
                "paired_delta_q16_q84_across_pairs": [
                    float(np.quantile(delta, .16)), float(np.quantile(delta, .84))
                ],
            }
        filters[band] = {
            "n_same_filter_pairs": expected_n,
            "filter_specific_connected_5deg_blocks": n_blocks,
            "largest_block_pairs": largest,
            "n_target_sky_cells": n_target_cells,
            "n_control_sky_cells": n_control_cells,
            "bootstrap_replicas": replicas,
            "seed": seed + (0 if band == "i" else 10000),
            "metrics": result_by_metric,
        }

    input_paths = [audit_dir / "summary.json", audit_dir / "independent_audit.json",
                   audit_dir / "frame_noise_correlation.tsv",
                   sample / "download_01/download_manifest.json",
                   sample / "control_download_01/download_manifest.json",
                   sample / "download_01/download_status.json",
                   sample / "control_download_01/download_status.json"]
    result = {
        "schema": "astrafek-icl-coadd-covariance-paired-block-bootstrap-v1",
        "status": "complete_diagnostic",
        "science_measurement": False,
        "method": "paired target-minus-control per-frame local-lag metric; resample connected 5-degree sky blocks independently within each filter",
        "n_total_pairs": 783,
        "filters": filters,
        "interpretation": "A small paired difference is evidence only that the measured local residual lag statistic is similar in targets and controls; it does not establish equal radial ICL-estimator covariance.",
        "limitations": [
            "The 9x9 finite-lag sum is an approximate square top-hat variance factor, not the exact radial ICL profile covariance.",
            "Per-frame residual kernels are medians over sampled patches and are not spatially resolved across each aperture.",
            "This comparison does not propagate source Poisson noise, upstream sky transfer, PSF wings, or mask-growth response.",
            "The bootstrap intervals describe spatial-block stability for the selected existing sample and do not promote the ICL result to a science measurement.",
        ],
        "inputs": {str(path): sha256(path) for path in input_paths},
        "parent_covariance_summary_sha256": sha256(audit_dir / "summary.json"),
        "parent_independent_audit_sha256": sha256(audit_dir / "independent_audit.json"),
        "analysis_script_sha256": sha256(Path(__file__).resolve()),
    }
    output.mkdir(parents=True, exist_ok=False)
    out_path = output / "summary.json"
    out_path.write_text(json.dumps(result, indent=2, ensure_ascii=False,
                                   allow_nan=False, sort_keys=True) + "\n",
                        encoding="utf-8")
    result["output_sha256"] = sha256(out_path)
    print(json.dumps({"status": result["status"], "science_measurement": False,
                      "filters": {band: {"pairs": values["n_same_filter_pairs"],
                                         "blocks": values["filter_specific_connected_5deg_blocks"],
                                         "metrics": values["metrics"]}
                                  for band, values in filters.items()},
                      "summary": str(out_path), "summary_sha256": result["output_sha256"]},
                     indent=2, ensure_ascii=False))
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", type=Path, default=DEFAULT_SAMPLE)
    parser.add_argument("--audit", type=Path, default=DEFAULT_AUDIT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--replicas", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=20261103)
    args = parser.parse_args()
    analyze(args.sample, args.audit, args.output, args.replicas, args.seed)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
