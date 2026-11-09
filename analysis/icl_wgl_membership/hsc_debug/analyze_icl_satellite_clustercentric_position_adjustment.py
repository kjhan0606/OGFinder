#!/usr/bin/env python3
"""Test whether host-centric satellite-position mix explains the ICL relation.

The runner uses the fixed 94-pair common-support table and independently
reconstructs each 50--75 kpc stack's contributing-member count and median
clustercentric radius from the per-member profiles. Inference uses joint
5-degree target/control sky-block resampling across filters.
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
DEFAULT_OUTPUT = ROOT / "icl_satellite_clustercentric_position_adjustment_2026-11-06_01"
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
    design = np.column_stack([np.ones(len(rows))] + [
        np.asarray([float(row[name]) for row in rows], dtype=float)
        for name in controls
    ])
    if len(rows) <= design.shape[1] + 2 or np.linalg.matrix_rank(design) != design.shape[1]:
        return float("nan")
    rx = x - design @ np.linalg.lstsq(design, x, rcond=None)[0]
    ry = y - design @ np.linalg.lstsq(design, y, rcond=None)[0]
    denominator = float(np.linalg.norm(rx) * np.linalg.norm(ry))
    if not np.isfinite(denominator) or denominator <= 0:
        return float("nan")
    return float(rx @ ry / denominator)


METRICS = {
    "icl_vs_satellite_baseline": (
        "icl_delta_200_300", "satellite_delta_50_75", ("redshift", "richness")),
    "icl_vs_satellite_adjusted_for_position_mix": (
        "icl_delta_200_300", "satellite_delta_50_75",
        ("redshift", "richness", "log_median_cluster_radius_kpc")),
    "position_mix_vs_icl": (
        "log_median_cluster_radius_kpc", "icl_delta_200_300", ("redshift", "richness")),
    "position_mix_vs_satellite_brightness": (
        "log_median_cluster_radius_kpc", "satellite_delta_50_75", ("redshift", "richness")),
}


def load_inputs() -> tuple[list[dict], dict]:
    for path in (PAIR_TABLE, MEMBER_TABLE, TARGET_MANIFEST, CONTROL_MANIFEST):
        if not path.is_file():
            raise FileNotFoundError(path)
    pair_rows = read_tsv(PAIR_TABLE)
    if len(pair_rows) != 94:
        raise ValueError(f"expected 94 supported cluster/filter pairs, got {len(pair_rows)}")
    member_rows = [row for row in read_tsv(MEMBER_TABLE)
                   if row["radial_bin_kpc"] == "50-75"]
    members_by_key: dict[tuple[str, str], list[dict[str, str]]] = {}
    for row in member_rows:
        key = str(row["cluster_id"]), row["filter"]
        members_by_key.setdefault(key, []).append(row)

    targets_doc = json.loads(TARGET_MANIFEST.read_text(encoding="utf-8"))
    controls_doc = json.loads(CONTROL_MANIFEST.read_text(encoding="utf-8"))
    targets = {str(row["id"]): row for row in targets_doc["clusters"]}
    controls = {str(row["target_cluster_id"]): row for row in controls_doc["clusters"]}
    if len(targets) != 800 or len(controls) != 798:
        raise ValueError("unexpected target/control manifest sizes")

    result = []
    for source in pair_rows:
        cluster_id, band = str(source["cluster_id"]), source["filter"]
        key = cluster_id, band
        members = members_by_key.get(key, [])
        if not members or len(members) != int(source["n_satellites"]):
            raise ValueError(f"contributing-member count differs from saved 50--75 kpc profile: {key}")
        radii = np.asarray([float(row["cluster_radius_kpc"]) for row in members], dtype=float)
        if not np.all(np.isfinite(radii)) or np.any(radii < 0):
            raise ValueError(f"invalid clustercentric radius in contributing members: {key}")
        if cluster_id not in targets or cluster_id not in controls:
            raise ValueError(f"missing target/control coordinates: {cluster_id}")
        target, control = targets[cluster_id], controls[cluster_id]
        median_radius = float(np.median(radii))
        if median_radius <= 0:
            raise ValueError(f"non-positive median clustercentric radius: {key}")
        result.append({
            "cluster_id": cluster_id, "filter": band,
            "icl_delta_200_300": float(source["icl_delta_200_300"]),
            "satellite_delta_50_75": float(source["satellite_delta_50_75"]),
            "n_satellites": len(members),
            "median_cluster_radius_kpc": median_radius,
            "log_median_cluster_radius_kpc": float(np.log(median_radius)),
            "q16_cluster_radius_kpc": float(np.quantile(radii, 0.16)),
            "q84_cluster_radius_kpc": float(np.quantile(radii, 0.84)),
            "redshift": float(source["redshift"]),
            "richness": float(source["richness"]),
            "target_ra_deg": float(target["ra_deg"]),
            "target_dec_deg": float(target["dec_deg"]),
            "control_ra_deg": float(control["ra_deg"]),
            "control_dec_deg": float(control["dec_deg"]),
        })
    result.sort(key=lambda row: (row["filter"], row["cluster_id"]))
    counts = [sum(row["filter"] == band for row in result) for band in BANDS]
    if counts != [26, 68]:
        raise ValueError(f"supported-sample band counts changed: {counts}")
    return result, {str(path): sha256(path) for path in
                    (PAIR_TABLE, MEMBER_TABLE, TARGET_MANIFEST, CONTROL_MANIFEST)}


def interval(values: list[float], requested: int) -> dict:
    array = np.asarray(values, dtype=float)
    array = array[np.isfinite(array)]
    if array.size < requested * 0.8:
        raise ValueError(f"insufficient complete joint bootstrap draws: {array.size}/{requested}")
    return {"ci95": np.quantile(array, [0.025, 0.975]).tolist(),
            "fraction_positive": float(np.mean(array > 0)),
            "valid_draws": int(array.size)}


def run(output: Path, replicas: int, seed: int) -> dict:
    output = Path(output).resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite {output}")
    if replicas < 2000:
        raise ValueError("at least 2,000 joint spatial-block resamples are required")
    rows, input_hashes = load_inputs()
    targets = [{"ra_deg": row["target_ra_deg"], "dec_deg": row["target_dec_deg"]} for row in rows]
    controls = [{"ra_deg": row["control_ra_deg"], "dec_deg": row["control_dec_deg"]} for row in rows]
    _, labels, n_target_cells, n_control_cells, n_blocks, largest = pair_block_labels(
        targets, controls, np.ones(len(rows), dtype=bool))
    labels = np.asarray(labels, dtype=str)
    groups = np.unique(labels)
    block_rows = {group: np.flatnonzero(labels == group) for group in groups}

    point = {}
    draws = {band: {**{metric: [] for metric in METRICS},
                    "position_adjustment_delta": []} for band in BANDS}
    for band in BANDS:
        subset = [row for row in rows if row["filter"] == band]
        point[band] = {name: partial_r(subset, *spec) for name, spec in METRICS.items()}
    rng = np.random.default_rng(seed)
    valid_draws = 0
    for _ in range(replicas):
        selected_groups = rng.choice(groups, size=len(groups), replace=True)
        indices = np.concatenate([block_rows[group] for group in selected_groups])
        bootstrap_rows = [rows[int(index)] for index in indices]
        draw_by_band = {}
        for band in BANDS:
            subset = [row for row in bootstrap_rows if row["filter"] == band]
            draw_by_band[band] = {name: partial_r(subset, *spec)
                                  for name, spec in METRICS.items()}
        if any(not np.isfinite(value) for band in BANDS
               for value in draw_by_band[band].values()):
            continue
        valid_draws += 1
        for band in BANDS:
            draws[band]["position_adjustment_delta"].append(
                draw_by_band[band]["icl_vs_satellite_adjusted_for_position_mix"]
                - draw_by_band[band]["icl_vs_satellite_baseline"])
            for metric, value in draw_by_band[band].items():
                draws[band][metric].append(value)
    if valid_draws < replicas * 0.8:
        raise ValueError(f"too few complete joint filter draws: {valid_draws}/{replicas}")

    result_rows = []
    for band in BANDS:
        subset = [row for row in rows if row["filter"] == band]
        medians = np.asarray([row["median_cluster_radius_kpc"] for row in subset])
        result = {
            "filter": band, "n_pairs": len(subset),
            "n_connected_5deg_blocks": len(set(labels[[i for i, row in enumerate(rows)
                                                         if row["filter"] == band]])),
            "median_hostcentric_radius_kpc": float(np.median(medians)),
            "q16_hostcentric_radius_kpc": float(np.quantile(medians, 0.16)),
            "q84_hostcentric_radius_kpc": float(np.quantile(medians, 0.84)),
        }
        for metric in METRICS:
            stats = interval(draws[band][metric], replicas)
            result[f"{metric}_partial_r"] = point[band][metric]
            result[f"{metric}_ci95"] = json.dumps(stats["ci95"])
        delta = interval(draws[band]["position_adjustment_delta"], replicas)
        result["position_adjustment_delta_partial_r"] = (
            point[band]["icl_vs_satellite_adjusted_for_position_mix"]
            - point[band]["icl_vs_satellite_baseline"])
        result["position_adjustment_delta_ci95"] = json.dumps(delta["ci95"])
        result_rows.append(result)

    output.mkdir(parents=True, exist_ok=False)
    joined_path, result_path = output / "position_adjusted_pairs.tsv", output / "position_adjustment_by_filter.tsv"
    for path, write_rows in ((joined_path, rows), (result_path, result_rows)):
        with path.open("x", encoding="utf-8", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(write_rows[0]), delimiter="\t",
                                    lineterminator="\n", extrasaction="raise")
            writer.writeheader()
            writer.writerows(write_rows)
    summary = {
        "schema": "astrafek-icl-satellite-clustercentric-position-adjustment-v1",
        "status": "complete_exploratory_sensitivity", "science_measurement": False,
        "method": {
            "sample": "fixed 94-pair common support; i=26, i2=68",
            "member_position_statistic": "median host-centric proper radius of contributing detected members in each 50--75 kpc satellite profile",
            "baseline_controls": ["redshift", "linear CAMIRA richness"],
            "adjusted_controls": ["redshift", "linear CAMIRA richness", "log median host-centric member radius"],
            "resampling": "joint 5-degree connected target/control blocks across filters",
            "replicas_requested": replicas, "replicas_valid": valid_draws, "seed": seed,
            "n_connected_blocks_joint": int(n_blocks), "largest_block_pairs_joint": int(largest),
            "n_target_cells": int(n_target_cells), "n_control_cells": int(n_control_cells),
            "confidence_intervals": "pointwise central percentile 95%; exploratory, not multiplicity-adjusted",
        },
        "input_sha256": input_hashes, "n_pairs": len(rows), "results": result_rows,
        "limitations": [
            "Median host-centric radius summarizes only the detected members entering the selected profile and does not correct membership incompleteness.",
            "The 94-pair common-support sample is selected; position adjustment is not a causal or population correction.",
            "Position summaries are noisy for stacks with few contributing members.",
            "Source-Poisson variance, PSF wings, sky transfer, and total image covariance remain unresolved.",
            "No ICL/satellite profile is corrected and science_measurement remains false.",
        ],
        "output_sha256": {str(path): sha256(path) for path in (joined_path, result_path)},
        "analysis_script_sha256": sha256(Path(__file__).resolve()),
    }
    (output / "summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2, allow_nan=False) + "\n",
        encoding="utf-8")
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--replicas", type=int, default=10000)
    parser.add_argument("--seed", type=int, default=20261106)
    args = parser.parse_args()
    summary = run(args.output, args.replicas, args.seed)
    print(json.dumps({"status": summary["status"], "science_measurement": False,
                      "output": str(args.output.resolve()), "results": summary["results"],
                      "replicas_valid": summary["method"]["replicas_valid"]},
                     ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
