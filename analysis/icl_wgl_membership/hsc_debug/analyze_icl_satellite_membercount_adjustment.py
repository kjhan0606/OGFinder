#!/usr/bin/env python3
"""Assess whether the number of detected members drives the 94-pair ICL relation.

Uses existing target-control profiles only. The member count is verified by
recounting the contributing 50--75 kpc profile rows. Inference resamples shared
5-degree target/control sky blocks and keeps the i/i2 measurements together.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from evaluate_prestack_pair_support import pair_block_labels

ROOT = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
PAIR_TABLE = ROOT / "satellite_outer_icl_association_perbin_2026-10-18_03/matched_cluster_values.tsv"
MEMBER_TABLE = ROOT / "satellite_outer_profile_perbin_full_2026-10-18_01/satellite_profiles.tsv"
TARGET_MANIFEST = ROOT / "sample_manifest.json"
CONTROL_MANIFEST = ROOT / "control_download_01/download_manifest.json"
DEFAULT_OUTPUT = ROOT / "icl_satellite_membercount_adjustment_2026-11-06_02"
BANDS = ("i", "i2")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with Path(path).open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def partial_r(rows: list[dict], x_name: str, y_name: str,
              controls: tuple[str, ...]) -> float:
    x = np.asarray([float(row[x_name]) for row in rows], dtype=float)
    y = np.asarray([float(row[y_name]) for row in rows], dtype=float)
    z = np.column_stack([np.ones(len(rows))] + [
        np.asarray([float(row[name]) for row in rows], dtype=float)
        for name in controls
    ])
    if len(rows) <= z.shape[1] + 2 or np.linalg.matrix_rank(z) != z.shape[1]:
        return float("nan")
    rx = x - z @ np.linalg.lstsq(z, x, rcond=None)[0]
    ry = y - z @ np.linalg.lstsq(z, y, rcond=None)[0]
    denom = float(np.linalg.norm(rx) * np.linalg.norm(ry))
    if not np.isfinite(denom) or denom <= 0:
        return float("nan")
    return float(rx @ ry / denom)


METRICS = {
    "icl_vs_satellite_baseline": ("icl_delta_200_300", "satellite_delta_50_75", ("redshift", "richness")),
    "icl_vs_satellite_adjusted_for_member_count": ("icl_delta_200_300", "satellite_delta_50_75", ("redshift", "richness", "log_n_satellites")),
    "member_count_vs_icl": ("log_n_satellites", "icl_delta_200_300", ("redshift", "richness")),
    "member_count_vs_satellite_brightness": ("log_n_satellites", "satellite_delta_50_75", ("redshift", "richness")),
}


def load_inputs() -> tuple[list[dict], dict]:
    for path in (PAIR_TABLE, MEMBER_TABLE, TARGET_MANIFEST, CONTROL_MANIFEST):
        if not path.is_file():
            raise FileNotFoundError(path)
    pair_rows = read_tsv(PAIR_TABLE)
    if len(pair_rows) != 94:
        raise ValueError(f"expected 94 supported cluster/filter rows, got {len(pair_rows)}")
    selected = {(r["cluster_id"], r["filter"]): r for r in pair_rows
                if r["radial_bin_kpc"] == "50-75"} if "radial_bin_kpc" in pair_rows[0] else {
                    (r["cluster_id"], r["filter"]): r for r in pair_rows
                }
    if len(selected) != 94:
        raise ValueError(f"expected 94 ICL/satellite support rows, got {len(selected)}")

    member_rows = [r for r in read_tsv(MEMBER_TABLE) if r["radial_bin_kpc"] == "50-75"]
    member_counts: dict[tuple[str, str], int] = {}
    member_ids: dict[tuple[str, str], set[str]] = {}
    for row in member_rows:
        key = row["cluster_id"], row["filter"]
        member_counts[key] = member_counts.get(key, 0) + 1
        member_ids.setdefault(key, set()).add(row["member_index"])
    if any(len(member_ids[key]) != count for key, count in member_counts.items()):
        raise ValueError("duplicate member index in a 50--75 kpc profile")

    targets_doc = json.loads(TARGET_MANIFEST.read_text(encoding="utf-8"))
    controls_doc = json.loads(CONTROL_MANIFEST.read_text(encoding="utf-8"))
    targets = {str(r["id"]): r for r in targets_doc["clusters"]}
    controls = {str(r["target_cluster_id"]): r for r in controls_doc["clusters"]}
    if len(targets) != 800 or len(controls) != 798:
        raise ValueError("unexpected parent target/control manifest sizes")

    rows = []
    for key, source in selected.items():
        cluster_id, band = key
        if key not in member_counts or int(source["n_satellites"]) != member_counts[key]:
            raise ValueError(f"saved satellite count does not match the contributing member rows: {key}")
        if int(source["n_satellites"]) < 1:
            raise ValueError(f"invalid zero member count in supported row: {key}")
        target, control = targets[cluster_id], controls.get(cluster_id)
        if control is None or str(control["id"]) != str(source.get("control_id", control["id"])):
            raise ValueError(f"target/control identity mismatch for {key}")
        count = int(source["n_satellites"])
        rows.append({
            "cluster_id": cluster_id, "filter": band,
            "icl_delta_200_300": float(source["icl_delta_200_300"]),
            "satellite_delta_50_75": float(source["satellite_delta_50_75"]),
            "n_satellites": count, "log_n_satellites": float(np.log(count)),
            "redshift": float(source["redshift"]),
            "richness": float(source["richness"]),
            "log_richness": float(np.log(float(source["richness"]))),
            "target_ra_deg": float(target["ra_deg"]),
            "target_dec_deg": float(target["dec_deg"]),
            "control_ra_deg": float(control["ra_deg"]),
            "control_dec_deg": float(control["dec_deg"]),
        })
    rows.sort(key=lambda row: (row["filter"], row["cluster_id"]))
    if [sum(r["filter"] == band for r in rows) for band in BANDS] != [26, 68]:
        raise ValueError("supported sample filter counts changed")
    return rows, {str(path): sha256(path) for path in
                  (PAIR_TABLE, MEMBER_TABLE, TARGET_MANIFEST, CONTROL_MANIFEST)}


def summarize(values: list[float], replicas: int) -> dict:
    array = np.asarray(values, dtype=float)
    array = array[np.isfinite(array)]
    if len(array) < replicas * 0.8:
        raise ValueError(f"insufficient valid bootstrap draws: {len(array)}/{replicas}")
    return {"ci95": np.quantile(array, [0.025, 0.975]).tolist(),
            "bootstrap_fraction_positive": float(np.mean(array > 0)),
            "replicas_valid": int(len(array))}


def run(output: Path, replicas: int, seed: int) -> dict:
    output = Path(output).resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite {output}")
    if replicas < 2000:
        raise ValueError("at least 2,000 spatial-block bootstrap replicates are required")
    rows, input_hashes = load_inputs()
    target_coords = [{"ra_deg": r["target_ra_deg"], "dec_deg": r["target_dec_deg"]} for r in rows]
    control_coords = [{"ra_deg": r["control_ra_deg"], "dec_deg": r["control_dec_deg"]} for r in rows]
    _, block_labels, n_target_cells, n_control_cells, n_blocks, largest = pair_block_labels(
        target_coords, control_coords, np.ones(len(rows), dtype=bool))
    block_labels = np.asarray(block_labels, dtype=str)
    groups = np.unique(block_labels)
    index_by_group = {group: np.flatnonzero(block_labels == group) for group in groups}

    point, draws = {}, {band: {**{metric: [] for metric in METRICS},
                               "member_count_adjustment_delta": []} for band in BANDS}
    for band in BANDS:
        subset = [r for r in rows if r["filter"] == band]
        point[band] = {name: partial_r(subset, *spec) for name, spec in METRICS.items()}
    rng = np.random.default_rng(seed)
    valid = 0
    for _ in range(replicas):
        chosen = rng.choice(groups, size=len(groups), replace=True)
        indices = np.concatenate([index_by_group[group] for group in chosen])
        replicate_rows = [rows[int(i)] for i in indices]
        replica_metrics = {}
        for band in BANDS:
            subset = [r for r in replicate_rows if r["filter"] == band]
            replica_metrics[band] = {name: partial_r(subset, *spec) for name, spec in METRICS.items()}
        if any(not np.isfinite(value) for band in BANDS for value in replica_metrics[band].values()):
            continue
        valid += 1
        for band in BANDS:
            for name, value in replica_metrics[band].items():
                draws[band][name].append(value)
            draws[band]["member_count_adjustment_delta"].append(
                replica_metrics[band]["icl_vs_satellite_adjusted_for_member_count"]
                - replica_metrics[band]["icl_vs_satellite_baseline"])
    if valid < replicas * 0.8:
        raise ValueError(f"too few complete joint filter draws: {valid}/{replicas}")

    result_rows = []
    for band in BANDS:
        subset = [r for r in rows if r["filter"] == band]
        counts = np.asarray([r["n_satellites"] for r in subset], dtype=int)
        result = {"filter": band, "n_pairs": len(subset), "n_connected_5deg_blocks": len(set(block_labels[[i for i, r in enumerate(rows) if r['filter'] == band]])),
                  "member_count_min": int(counts.min()), "member_count_median": float(np.median(counts)),
                  "member_count_max": int(counts.max())}
        for name in METRICS:
            result[name + "_partial_r"] = point[band][name]
            result[name + "_bootstrap_ci95"] = json.dumps(summarize(draws[band][name], replicas)["ci95"])
            result[name + "_bootstrap_fraction_positive"] = summarize(draws[band][name], replicas)["bootstrap_fraction_positive"]
        result["member_count_adjustment_delta_partial_r"] = point[band]["icl_vs_satellite_adjusted_for_member_count"] - point[band]["icl_vs_satellite_baseline"]
        result["member_count_adjustment_delta_ci95"] = json.dumps(summarize(draws[band]["member_count_adjustment_delta"], replicas)["ci95"])
        result_rows.append(result)

    output.mkdir(parents=True, exist_ok=False)
    joined_path = output / "membercount_adjusted_pairs.tsv"
    with joined_path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    result_path = output / "membercount_adjustment_by_filter.tsv"
    with result_path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(result_rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(result_rows)
    summary = {
        "schema": "astrafek-icl-satellite-membercount-adjustment-v1",
        "status": "complete_exploratory_sensitivity",
        "science_measurement": False,
        "question": "Does the number of contributing detected satellites account for the supported ICL--satellite association?",
        "method": {
            "sample": "fixed 94-pair common-support sample; i=26, i2=68",
            "outcome": "satellite equal-member mean excess at 50--75 kpc versus ICL target-control excess at 200--300 kpc",
            "baseline_controls": ["redshift", "linear richness"],
            "adjusted_controls": ["redshift", "linear richness", "log contributing satellite count"],
            "uncertainty": "joint resampling of 5-degree connected target/control sky blocks across filters",
            "replicas_requested": replicas, "replicas_valid": valid, "seed": seed,
            "n_connected_blocks_joint": int(n_blocks), "largest_block_pairs_joint": int(largest),
            "n_target_cells": int(n_target_cells), "n_control_cells": int(n_control_cells),
            "confidence_intervals": "pointwise central percentile 95%; exploratory and not multiplicity-adjusted",
        },
        "input_sha256": input_hashes,
        "analysis_script_sha256": sha256(Path(__file__).resolve()),
        "n_pairs": len(rows), "results": result_rows,
        "limitations": [
            "member count is a precision/support proxy and may also be a mediator or collider; adjustment is not causal",
            "the selected 94-pair sample and 5-degree approximate sky blocks remain",
            "member detection completeness, source-Poisson noise, PSF wings, sky-transfer covariance, and pixel coadd weights remain unresolved",
            "this analysis does not change science_measurement=false",
        ],
        "output_sha256": {str(path): sha256(path) for path in (joined_path, result_path)},
    }
    summary_path = output / "summary.json"
    summary_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--replicas", type=int, default=10000)
    parser.add_argument("--seed", type=int, default=20261106)
    args = parser.parse_args()
    summary = run(args.output, args.replicas, args.seed)
    print(json.dumps({"status": summary["status"], "science_measurement": False,
                      "output": str(args.output.resolve()), "results": summary["results"]},
                     ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
