#!/usr/bin/env python3
"""Audit and spatially summarize fixed-mask joint-injection batch results."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from run_joint_icl_satellite_batch import (
    DEFAULT_DATA,
    DEFAULT_ICL_PAIRS,
    DEFAULT_SATELLITE_TABLE,
    select_jointly_supported_pairs,
)


DEFAULT_MANIFEST = DEFAULT_DATA / "sample_manifest.json"
METRICS = (
    "icl_response_to_diffuse_unit_injection",
    "icl_response_to_all_satellite_envelopes",
    "satellite_response_to_smooth_diffuse",
    "satellite_response_to_satellite_unit_injection",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def spatial_bootstrap(values: np.ndarray, labels: np.ndarray,
                      replicas: int, seed: int) -> dict:
    values = np.asarray(values, dtype=float)
    labels = np.asarray(labels)
    groups = np.unique(labels)
    if len(values) != len(labels) or not len(values):
        raise ValueError("spatial bootstrap requires paired values and labels")
    group_indices = [np.flatnonzero(labels == group) for group in groups]
    rng = np.random.default_rng(seed)
    means = np.empty(replicas, dtype=float)
    for i in range(replicas):
        chosen = rng.integers(0, len(groups), size=len(groups))
        indices = np.concatenate([group_indices[j] for j in chosen])
        means[i] = values[indices].mean()
    return {
        "n_clusters": int(len(values)),
        "n_spatial_blocks_5deg": int(len(groups)),
        "largest_block_clusters": int(max(map(len, group_indices))),
        "mean": float(values.mean()),
        "median": float(np.median(values)),
        "cluster_q16_q84": np.quantile(values, [0.16, 0.84]).tolist(),
        "spatial_block_bootstrap_95": np.quantile(means, [0.025, 0.975]).tolist(),
        "bootstrap_replicas": int(replicas),
    }


def audit(batch_summary_path: Path, manifest_path: Path,
          icl_pairs_path: Path, satellite_profiles_path: Path,
          output_path: Path, replicas: int = 5000,
          seed: int = 20261018) -> dict:
    if output_path.exists():
        raise FileExistsError(f"refusing to overwrite {output_path}")
    summary = json.loads(batch_summary_path.read_text(encoding="utf-8"))
    expected = select_jointly_supported_pairs(icl_pairs_path, satellite_profiles_path)
    observed = [(str(r["cluster_id"]), str(r["filter"])) for r in summary.get("records", [])]
    if summary.get("status") != "complete" or summary.get("n_errors") != 0:
        raise ValueError("batch is incomplete; no inferential summary will be produced")
    if observed != expected:
        raise ValueError("batch output IDs/order do not match the frozen joint-support selection")
    if summary.get("n_completed") != len(expected):
        raise ValueError("batch record count differs from the selected support sample")
    if summary.get("science_measurement") is not False:
        raise ValueError("unexpected science-measurement gate")
    if summary.get("inputs", {}).get(str(icl_pairs_path)) != sha256(icl_pairs_path):
        raise ValueError("ICL support input hash mismatch")
    if summary.get("inputs", {}).get(str(satellite_profiles_path)) != sha256(satellite_profiles_path):
        raise ValueError("satellite support input hash mismatch")

    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    clusters = {str(c["id"]): c for c in manifest["clusters"]}
    rows = []
    max_additivity_residual = 0.0
    for record in summary["records"]:
        path = Path(record["path"])
        if not path.is_file() or sha256(path) != record["sha256"]:
            raise ValueError(f"frame output missing or hash mismatch: {path}")
        frame = json.loads(path.read_text(encoding="utf-8"))
        if (frame.get("cluster_id") != str(record["cluster_id"])
                or frame.get("filter") != record["filter"]
                or frame.get("science_measurement") is not False
                or frame.get("observational_validation_passed") is not False
                or frame.get("flux_covariance_available") is not False):
            raise ValueError(f"frame identity/gate mismatch: {path}")
        if frame.get("script_sha256") != summary.get("single_frame_script_sha256"):
            raise ValueError(f"single-frame script hash mismatch: {path}")
        if str(record["cluster_id"]) not in clusters:
            raise ValueError(f"cluster ID absent from sample manifest: {record['cluster_id']}")
        response = frame["transfer_response"]
        if not all(math.isfinite(float(response[key])) for key in METRICS):
            raise ValueError(f"non-finite transfer response: {path}")
        residual = max(abs(float(response["both_minus_sum_residual_icl_njy_arcsec2"])),
                       abs(float(response["both_minus_sum_residual_satellite_njy_arcsec2"])))
        max_additivity_residual = max(max_additivity_residual, residual)
        c = clusters[str(record["cluster_id"])]
        ra, dec = float(c["ra_deg"]), float(c["dec_deg"])
        rows.append({"cluster_id": str(record["cluster_id"]), "filter": record["filter"],
                     "ra_deg": ra, "dec_deg": dec,
                     "spatial_block_5deg": f"{int(ra//5)}_{int((dec+90)//5)}",
                     **{key: float(response[key]) for key in METRICS}})
    if max_additivity_residual > 1e-4:
        raise ValueError(f"fixed-mask response is not additive: max residual={max_additivity_residual}")

    by_filter = {}
    for band in sorted({r["filter"] for r in rows}):
        subset = [r for r in rows if r["filter"] == band]
        labels = np.asarray([r["spatial_block_5deg"] for r in subset])
        metrics = {}
        for metric_index, metric in enumerate(METRICS):
            vals = np.asarray([r[metric] for r in subset], dtype=float)
            metrics[metric] = spatial_bootstrap(vals, labels, replicas,
                                                seed + 100 * metric_index + (0 if band == "i" else 1))
        by_filter[band] = {"n_clusters": len(subset),
                           "n_spatial_blocks_5deg": len(set(labels)),
                           "metrics": metrics}

    result = {
        "schema": "hsc-joint-icl-satellite-batch-audit-v1",
        "status": "verified_exploratory",
        "science_measurement": False,
        "observational_validation_passed": False,
        "flux_covariance_available": False,
        "batch_summary": str(batch_summary_path),
        "batch_summary_sha256": sha256(batch_summary_path),
        "manifest_sha256": sha256(manifest_path),
        "icl_pair_input_sha256": sha256(icl_pairs_path),
        "satellite_profile_input_sha256": sha256(satellite_profiles_path),
        "n_clusters": len(rows),
        "n_spatial_blocks_5deg_total": len({r["spatial_block_5deg"] for r in rows}),
        "max_additivity_residual_njy_arcsec2": max_additivity_residual,
        "by_filter": by_filter,
        "records": rows,
        "bootstrap": {"unit": "5-degree sky cells, with all clusters in each draw retained together",
                      "replicas": replicas, "seed": seed},
        "limitations": [
            "The input sample is selected by existing K2.5/40% ICL and 50--75 kpc satellite support.",
            "Spatial intervals quantify variation in fixed-mask injection response only; they are not confidence intervals for an astrophysical signal.",
            "Masks/detections are fixed, and the injection does not include PSF convolution, sky-processing transfer or source Poisson noise.",
            "The run does not estimate full image-level measurement covariance or the ICL--satellite population relation.",
        ],
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n",
                           encoding="utf-8")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--batch-summary", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--icl-pairs", type=Path, default=DEFAULT_ICL_PAIRS)
    parser.add_argument("--satellite-profiles", type=Path, default=DEFAULT_SATELLITE_TABLE)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--replicas", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=20261018)
    args = parser.parse_args()
    if args.replicas < 100:
        parser.error("at least 100 bootstrap replicas are required")
    result = audit(args.batch_summary.resolve(), args.manifest.resolve(),
                   args.icl_pairs.resolve(), args.satellite_profiles.resolve(),
                   args.output.resolve(), args.replicas, args.seed)
    print(json.dumps({k: result[k] for k in
                      ("status", "n_clusters", "n_spatial_blocks_5deg_total",
                       "max_additivity_residual_njy_arcsec2", "by_filter")}, indent=2))


if __name__ == "__main__":
    main()
