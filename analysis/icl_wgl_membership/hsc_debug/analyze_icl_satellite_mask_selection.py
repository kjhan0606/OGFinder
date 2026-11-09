#!/usr/bin/env python3
"""Measure how annulus-mask coverage selects ICL-supported satellite candidates.

Uses saved mask-component audit rows, the matching SEP detection catalogues,
and cluster metadata. No image pixels or HSC services are accessed. Reports
descriptive quartile contrasts with cluster-block bootstrap intervals.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np


DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
AUDIT = DATA / "icl_satellite_outer_coverage_loss_2026-11-05_03"
MANIFEST = DATA / "download_01/download_manifest.json"
PREP = DATA / "satellite_detection_psf_2026-10-18_01"
DEFAULT_OUTPUT = DATA / "icl_satellite_mask_selection_2026-11-05_01"
PREDICTORS = ("cluster_radius_kpc", "detection_snr", "redshift", "richness")
FILTERS = ("i", "i2")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with Path(path).open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def numeric(row: dict, key: str) -> float:
    value = float(row[key])
    if not np.isfinite(value):
        raise ValueError(f"non-finite {key}")
    return value


def block_contrast(rows: list[dict], predictor: str, outcome: str,
                   replicas: int, seed: int) -> dict:
    ordered = sorted(rows, key=lambda row: (row[predictor], row["cluster_id"]))
    n = len(ordered)
    labels = np.empty(n, dtype=np.int8)
    for index, row in enumerate(ordered):
        labels[index] = min(3, index * 4 // n)
    values = np.asarray([int(row[outcome]) for row in ordered], dtype=float)
    groups: dict[str, list[int]] = {}
    for index, row in enumerate(ordered):
        groups.setdefault(row["cluster_id"], []).append(index)
    cluster_ids = sorted(groups)
    if len(cluster_ids) < 8:
        raise ValueError(f"too few independent cluster blocks: {len(cluster_ids)}")

    def rates(indices: np.ndarray) -> tuple[float, float]:
        selected_labels = labels[indices]
        selected_values = values[indices]
        q1, q4 = selected_labels == 0, selected_labels == 3
        return float(selected_values[q1].mean()), float(selected_values[q4].mean())

    all_indices = np.arange(n)
    low, high = rates(all_indices)
    rng = np.random.default_rng(seed)
    draws = np.empty(replicas, dtype=float)
    for rep in range(replicas):
        sampled = rng.choice(cluster_ids, size=len(cluster_ids), replace=True)
        indices = np.concatenate([np.asarray(groups[cluster], dtype=int) for cluster in sampled])
        low_rate, high_rate = rates(indices)
        draws[rep] = high_rate - low_rate
    return {
        "n_candidates": n, "n_cluster_blocks": len(cluster_ids),
        "q1_predictor_range": [ordered[0][predictor], ordered[max(0, n // 4 - 1)][predictor]],
        "q4_predictor_range": [ordered[3 * n // 4][predictor], ordered[-1][predictor]],
        "q1_outcome_fraction": low, "q4_outcome_fraction": high,
        "q4_minus_q1": high - low,
        "cluster_block_bootstrap_ci95": np.quantile(draws, [0.025, 0.975]).tolist(),
        "bootstrap_replicas": replicas,
    }


def run(output: Path, replicas: int, seed: int) -> dict:
    if output.exists():
        raise FileExistsError(f"refusing to overwrite {output}")
    summary_path = AUDIT / "summary.json"
    candidates_path = AUDIT / "candidate_coverage_components.tsv"
    frame_path = AUDIT / "frame_support_components.tsv"
    for path in (summary_path, candidates_path, frame_path, MANIFEST):
        if not path.is_file():
            raise FileNotFoundError(path)
    audit_summary = json.loads(summary_path.read_text(encoding="utf-8"))
    if (audit_summary.get("status") != "complete_offline_diagnostic"
            or audit_summary.get("science_measurement") is not False
            or audit_summary.get("schema") != "astrafek-icl-satellite-outer-coverage-loss-v3"):
        raise ValueError("coverage-loss input is not the verified exact-pixel audit")
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    cluster_info = {}
    for row in manifest["clusters"]:
        cluster_info[str(row["id"])] = {
            "redshift": float(row["z"]), "richness": float(row["richness"]),
        }

    by_cluster_detection = {}
    detection_hashes = {}
    candidate_rows = read_tsv(candidates_path)
    for ident in sorted({row["cluster_id"] for row in candidate_rows}):
        path = PREP / f"cluster_{ident}" / "detections.tsv"
        detection_hashes[str(path)] = sha256(path)
        detections = read_tsv(path)
        by_cluster_detection[ident] = {str(row["NUMBER"]): row for row in detections}

    enriched = []
    for row in candidate_rows:
        ident = str(row["cluster_id"])
        detection = by_cluster_detection[ident].get(str(row["detection_number"]))
        if detection is None:
            raise ValueError(f"candidate detection missing for {ident}/{row['detection_number']}")
        meta = cluster_info[ident]
        snr = numeric(detection, "SNR_NOMINAL")
        item = {
            "cluster_id": ident, "filter": row["filter"],
            "cluster_radius_kpc": numeric(row, "cluster_radius_kpc"),
            "detection_snr": snr,
            "redshift": meta["redshift"], "richness": meta["richness"],
            "coverage_pass": int(numeric(row, "after_source_neighbor_masks_coverage_50_75") >= 0.40),
            "profile_retained": int(row["outer_annulus_profile_retained"] == "True"),
        }
        enriched.append(item)

    results = {"coverage_pass": {}, "profile_retained": {}}
    for band_index, band in enumerate(FILTERS):
        rows = [row for row in enriched if row["filter"] == band]
        if len(rows) < 100:
            raise ValueError(f"insufficient {band} candidate support: {len(rows)}")
        results["coverage_pass"][band] = {
            predictor: block_contrast(rows, predictor, "coverage_pass", replicas,
                                      seed + band_index * 100 + predictor_index)
            for predictor_index, predictor in enumerate(PREDICTORS)
        }
        results["profile_retained"][band] = {
            predictor: block_contrast(rows, predictor, "profile_retained", replicas,
                                      seed + 1000 + band_index * 100 + predictor_index)
            for predictor_index, predictor in enumerate(PREDICTORS)
        }

    output.mkdir(parents=True, exist_ok=False)
    table_path = output / "candidate_selection_inputs.tsv"
    fields = ["cluster_id", "filter", *PREDICTORS, "coverage_pass", "profile_retained"]
    with table_path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(enriched)
    summary = {
        "schema": "astrafek-icl-satellite-mask-selection-v1",
        "status": "complete_offline_diagnostic", "science_measurement": False,
        "hsc_network_used": False, "science_image_pixels_read": False,
        "n_candidates": len(enriched), "n_cluster_blocks": len({x["cluster_id"] for x in enriched}),
        "coverage_threshold": 0.40, "bootstrap_unit": "cluster_id",
        "bootstrap_replicas": replicas, "seed": seed,
        "results": results,
        "limitations": [
            "Quartile contrasts are descriptive, with no joint adjustment for the correlated predictors.",
            "Cluster-block bootstrap accounts for candidate clustering but not calibration, PSF, sky, or mask-model uncertainty.",
            "Detection SNR is taken from the saved SEP catalogue and is not a completeness-calibrated flux estimate.",
            "This estimates mask-selection patterns only; it is not a flux correction or ICL–satellite measurement.",
        ],
        "input_sha256": {str(path): sha256(path) for path in
                         (summary_path, candidates_path, frame_path, MANIFEST)},
        "detection_table_sha256": detection_hashes,
        "output_sha256": {str(table_path): sha256(table_path)},
    }
    (output / "summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8")
    print(json.dumps({"status": summary["status"], "output": str(output),
                      "n_candidates": len(enriched), "n_clusters": summary["n_cluster_blocks"],
                      "coverage_pass_q4_minus_q1": {
                          band: {key: value["q4_minus_q1"] for key, value in predictors.items()}
                          for band, predictors in results["coverage_pass"].items()}},
                     ensure_ascii=False, indent=2))
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--replicas", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=20261105)
    args = parser.parse_args()
    if args.replicas < 1000:
        parser.error("at least 1,000 cluster-bootstrap replicas are required")
    run(args.output.resolve(), args.replicas, args.seed)


if __name__ == "__main__":
    main()
