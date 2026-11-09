#!/usr/bin/env python3
"""Scale the approved HSC visit-split covariance diagnostic to complete pairs.

This is an offline exploratory noise audit over already downloaded HSC-I
visit warps. It does not correct the ICL profile or promote science status.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np

import analyze_icl_visit_split_covariance_pilot as pilot

ROOT = pilot.ROOT
PLAN_DIR = pilot.PLAN_DIR
DOWNLOAD_DIR = pilot.DOWNLOAD_DIR
TARGET_MANIFEST = pilot.TARGET_MANIFEST
CONTROL_MANIFEST = pilot.CONTROL_MANIFEST
PAIR_TABLE = pilot.PAIR_TABLE
DEFAULT_OUTPUT = ROOT / "icl_visit_split_covariance_sample_2026-11-06_02"
N_BOOT = 5000
SEED = 20261106


def read_tsv(path: Path) -> list[dict]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def paired_block_bootstrap(paired: list[dict], metric: str) -> dict:
    result = {}
    for band in ("i", "i2"):
        rows = [r for r in paired if r["filter"] == band
                and r.get(f"{metric}_target") is not None
                and r.get(f"{metric}_control") is not None]
        deltas = np.asarray([r[f"{metric}_target"] - r[f"{metric}_control"] for r in rows], dtype=float)
        if not rows:
            result[band] = {"n_pairs": 0, "n_blocks": 0, "median_target_minus_control": None,
                            "block_bootstrap_ci95": None}
            continue
        blocks = sorted({r["block"] for r in rows})
        index = {block: np.asarray([i for i, row in enumerate(rows) if row["block"] == block], dtype=int)
                 for block in blocks}
        rng = np.random.default_rng(SEED + (0 if band == "i" else 1) + sum(ord(c) for c in metric))
        estimates = np.empty(N_BOOT, dtype=float)
        for j in range(N_BOOT):
            chosen = rng.choice(blocks, size=len(blocks), replace=True)
            selected = np.concatenate([index[b] for b in chosen])
            estimates[j] = float(np.median(deltas[selected]))
        result[band] = {
            "n_pairs": len(rows), "n_blocks": len(blocks),
            "median_target_minus_control": float(np.median(deltas)),
            "block_bootstrap_ci95": [float(x) for x in np.quantile(estimates, [0.025, 0.975])],
            "n_replicates": N_BOOT,
        }
    return result


def run(output_dir: Path) -> dict:
    output_dir = Path(output_dir).resolve()
    if output_dir.exists():
        raise FileExistsError(f"refusing to overwrite {output_dir}")
    input_paths = [PLAN_DIR / "plan.json", DOWNLOAD_DIR / "summary.json",
                   DOWNLOAD_DIR / "checkpoint.json", TARGET_MANIFEST,
                   CONTROL_MANIFEST, PAIR_TABLE]
    for path in input_paths:
        if not path.is_file():
            raise FileNotFoundError(path)
    input_hashes = {str(path): pilot.sha256(path) for path in input_paths}
    if (input_hashes[str(PLAN_DIR / "plan.json")] != pilot.PLAN_SHA256
            or input_hashes[str(DOWNLOAD_DIR / "summary.json")] != pilot.DOWNLOAD_SHA256):
        raise ValueError("approved plan/download summary hash mismatch")
    plan = json.loads((PLAN_DIR / "plan.json").read_text(encoding="utf-8"))
    download = json.loads((DOWNLOAD_DIR / "summary.json").read_text(encoding="utf-8"))
    checkpoint = json.loads((DOWNLOAD_DIR / "checkpoint.json").read_text(encoding="utf-8"))
    if (download.get("completed_archive_count") != 184
            or download.get("http_404_request_count") != 4
            or download.get("science_measurement") is not False):
        raise ValueError("approved acquisition checkpoint changed")
    pairs = read_tsv(PAIR_TABLE)
    target_manifest = json.loads(TARGET_MANIFEST.read_text(encoding="utf-8"))
    control_manifest = json.loads(CONTROL_MANIFEST.read_text(encoding="utf-8"))
    coadd_rows = {
        "target": {str(r["id"]): r for r in target_manifest["clusters"]},
        "control": {str(r["id"]): r for r in control_manifest["clusters"]},
    }
    plan_index = {r["request_id"]: i for i, r in enumerate(plan["requests"])}
    record_map = checkpoint["requests"]
    center_rows = []
    kernel_rows = []
    paired_rows = []
    pair_errors = []
    expected_pair_total = len(pairs)
    for pair_no, pair in enumerate(pairs, 1):
        pair_key = (pair["filter"], str(pair["cluster_id"]))
        ids = {"target": str(pair["cluster_id"]), "control": str(pair["control_id"])}
        request_ids = {role: f"{pair['filter']}:{role}:{ids[role]}" for role in ids}
        entries = {}
        missing = []
        for role, request_id in request_ids.items():
            if request_id not in plan_index:
                missing.append({"role": role, "request_id": request_id, "status": "not_in_approved_plan"})
                continue
            idx = plan_index[request_id]
            record = record_map.get(f"request_{idx:03d}", {})
            if record.get("status") == "http_404_permanent" or not record.get("archive_sha256"):
                missing.append({"role": role, "request_id": request_id,
                                "status": record.get("status", "missing_checkpoint")})
                continue
            entries[role] = (idx, plan["requests"][idx], record)
        if missing:
            pair_errors.append({"filter": pair["filter"], "cluster_id": ids["target"],
                                "status": "approved_archive_gap", "details": missing})
            continue
        current = {}
        try:
            for role in ("target", "control"):
                idx, request, record = entries[role]
                center = pilot.center_covariance(
                    request, idx, pair, role, coadd_rows[role][ids[role]], record, output_dir)
                center.update({"filter": pair["filter"], "cluster_id": ids["target"],
                               "control_id": ids["control"], "block": pair["connected_5deg_block"]})
                lag8_counts = np.asarray(center["kernel_pair_counts"], dtype=np.int64)
                lag8_values = center["kernel"]
                center["lag_sum_8_full_support"] = bool(
                    np.all(lag8_counts > 0)
                    and all(value is not None for row in lag8_values for value in row))
                current[role] = center
        except Exception as exc:
            pair_errors.append({"filter": pair["filter"], "cluster_id": ids["target"],
                                "status": "analysis_error", "error": f"{type(exc).__name__}: {exc}"})
            continue
        center_rows.extend(current.values())
        for role, center in current.items():
            slim = {k: v for k, v in center.items() if k not in ("kernel", "kernel_pair_counts")}
            slim["split_a_visits"] = ",".join(map(str, center["split_a_visits"]))
            slim["split_b_visits"] = ",".join(map(str, center["split_b_visits"]))
            center_rows[-2 if role == "target" else -1] = center
            # A separate table row keeps structured visit lists compact and legible.
            center_rows[-2 if role == "target" else -1]["_slim"] = slim
            for dy in range(-8, 9):
                for dx in range(-8, 9):
                    cov = center["kernel"][dy + 8][dx + 8]
                    count = center["kernel_pair_counts"][dy + 8][dx + 8]
                    kernel_rows.append({"filter": pair["filter"], "cluster_id": ids["target"],
                                        "role": role, "block": pair["connected_5deg_block"],
                                        "dy": dy, "dx": dx, "covariance": cov,
                                        "pair_count": count})
        row = {"filter": pair["filter"], "cluster_id": ids["target"],
               "control_id": ids["control"], "block": pair["connected_5deg_block"]}
        for metric, key in (("rms", "split_difference_rms_normalized_before_plane"),
                            ("zero_lag", "normalized_covariance_zero_lag"),
                            ("coverage", "coverage_in_geometric_annulus")):
            row[f"{metric}_target"] = current["target"].get(key)
            row[f"{metric}_control"] = current["control"].get(key)
        for role in ("target", "control"):
            center = current[role]
            row[f"lag_sum_8_{role}"] = (center["symmetric_cumulative_lag_sum"]["8"]
                                         if center["lag_sum_8_full_support"] else None)
            row[f"lag_sum_8_full_support_{role}"] = center["lag_sum_8_full_support"]
        paired_rows.append(row)
        print(f"processed {pair_no}/{expected_pair_total}: {pair['filter']} {ids['target']} "
              f"({len(paired_rows)} complete pairs)", flush=True)

    metrics = {name: paired_block_bootstrap(paired_rows, name)
               for name in ("rms", "zero_lag", "lag_sum_8", "coverage")}
    output_dir.mkdir(parents=True, exist_ok=False)
    center_table = output_dir / "per_center.tsv"
    center_fields = sorted({k for row in center_rows for k in row.get("_slim", row).keys()})
    with center_table.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=center_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in center_rows:
            writer.writerow(row["_slim"])
    kernel_table = output_dir / "covariance_kernels.tsv"
    kernel_fields = ["filter", "cluster_id", "role", "block", "dy", "dx", "covariance", "pair_count"]
    with kernel_table.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=kernel_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(kernel_rows)
    paired_table = output_dir / "paired_comparisons.tsv"
    paired_fields = list(paired_rows[0]) if paired_rows else ["filter", "cluster_id", "control_id", "block"]
    with paired_table.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=paired_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(paired_rows)
    summary = {
        "schema": "astrafek-icl-visit-split-covariance-sample-v1",
        "status": "complete_exploratory" if not pair_errors else "complete_with_data_gaps",
        "science_measurement": False, "hsc_network_used": False, "gpu_used": False,
        "input_pair_count": expected_pair_total,
        "complete_pair_count": len(paired_rows),
        "complete_pairs_by_filter": {band: sum(r["filter"] == band for r in paired_rows) for band in ("i", "i2")},
        "analyzed_center_count": len(center_rows),
        "failed_or_missing_pairs": pair_errors,
        "paired_block_bootstrap": metrics,
        "method": "alternating sorted visit-ID halves, inverse-variance stacks, nearest-neighbor affine WCS mapping, normalized A-B difference, local affine-plane/4-sigma clipped covariance in the 200-300 proper-kpc annulus; paired 5-degree block bootstrap",
        "limitations": [
            "Split visit sets are not identical between each target and its control; the pair table records both exact partitions.",
            "Nearest-neighbor mapping does not reproduce the HSC coadd resampling kernel.",
            "The split difference retains residual sky, source variability/structure, visit-dependent PSF and mask effects, and variance-plane calibration errors.",
            "The local plane fit, sigma clipping, overlapping patch windows, and selected annulus mask define an estimator-conditional covariance diagnostic.",
            "No correction is applied to ICL or weak-lensing science profiles; science_measurement remains false.",
        ],
        "input_sha256": input_hashes,
        "per_center_sha256": pilot.sha256(center_table),
        "kernel_table_sha256": pilot.sha256(kernel_table),
        "paired_table_sha256": pilot.sha256(paired_table),
    }
    summary_path = output_dir / "summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True, allow_nan=False) + "\n",
                            encoding="utf-8")
    print(json.dumps({"status": summary["status"], "science_measurement": False,
                      "complete_pairs": summary["complete_pair_count"],
                      "complete_pairs_by_filter": summary["complete_pairs_by_filter"],
                      "output": str(summary_path), "paired_block_bootstrap": metrics}, indent=2), flush=True)
    return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    run(parser.parse_args().output_dir)
