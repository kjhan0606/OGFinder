#!/usr/bin/env python3
"""Describe selection induced by 20% annular coverage in forced satellite stacks.

Uses existing candidate diagnostics and measured profile outputs only. It
reports position- and cluster-level support without treating selected
unmatched CAMIRA positions as confirmed satellites.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import re

import numpy as np


DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
RUNS = {
    "coverage_60": DATA / "satellite_forced_outer_profile_unmatched_camira_unmasked_2026-11-06_01",
    "coverage_20": DATA / "satellite_forced_outer_profile_unmatched_camira_unmasked_cov20_2026-11-06_01",
}
TARGETS = DATA / "download_01/download_manifest.json"
DEFAULT_OUTPUT = DATA / "icl_satellite_forced_position_coverage_selection_2026-11-06_01"
PROFILE_BINS = ["5-10", "10-15", "15-20", "20-30", "30-40", "40-50", "50-75"]
CLUSTER_RADIUS_BINS = [(50, 100), (100, 200), (200, 300), (300, 500)]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def read_tsv(path: Path) -> list[dict[str, str]]:
    with Path(path).open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty table: {path}")
    with path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def summarize(values: list[float]) -> dict:
    arr = np.asarray(values, dtype=float)
    return {"n": int(len(arr)), "median": float(np.median(arr)) if len(arr) else None,
            "q16": float(np.quantile(arr, .16)) if len(arr) else None,
            "q84": float(np.quantile(arr, .84)) if len(arr) else None,
            "p95": float(np.quantile(arr, .95)) if len(arr) else None}


def read_run(run_path: Path, expected_floor: float) -> tuple[dict, list[dict], set[tuple[str, str, str, str]]]:
    summary = read_json(run_path / "summary.json")
    recorded_floor = summary.get("coverage_floor")
    if recorded_floor is None:
        match = re.search(r"target coverage >= ([0-9.]+)", summary.get("method", {}).get("annulus_support", ""))
        recorded_floor = float(match.group(1)) if match else None
    if (summary.get("status") != "complete_exploratory"
            or summary.get("science_measurement") is not False
            or recorded_floor != expected_floor
            or summary.get("method", {}).get("mask_other_unmatched_camira") is not False
            or summary.get("n_candidate_positions") != 3197
            or summary.get("n_target_frames_processed") != 721):
        raise ValueError(f"unexpected or incomplete run: {run_path}")
    diagnostics = read_tsv(run_path / "candidate_diagnostics.tsv")
    profiles = read_tsv(run_path / "forced_satellite_profiles.tsv")
    if len(diagnostics) != 3197:
        raise ValueError(f"incomplete candidate diagnostics: {run_path}")
    keys = {(row["cluster_id"], row["filter"], row["member_index"], row["radial_bin_kpc"])
            for row in profiles}
    return summary, diagnostics, keys


def run(output: Path) -> dict:
    if output.exists():
        raise FileExistsError(f"refusing existing output directory: {output}")
    target_manifest = read_json(TARGETS)
    frame_by_id = {str(row["id"]): row for row in target_manifest["clusters"]}
    summaries, diagnostics_by_floor, profile_keys = {}, {}, {}
    for name, floor in (("coverage_60", 0.6), ("coverage_20", 0.2)):
        summaries[name], diagnostics_by_floor[name], profile_keys[name] = read_run(RUNS[name], floor)

    records = []
    for row20, row60 in zip(diagnostics_by_floor["coverage_20"], diagnostics_by_floor["coverage_60"]):
        key20 = (row20["cluster_id"], row20["filter"], row20["member_index"])
        key60 = (row60["cluster_id"], row60["filter"], row60["member_index"])
        if key20 != key60:
            raise ValueError("coverage-run candidate keys are not aligned")
        ident, band, member = key20
        frame = frame_by_id.get(ident)
        if frame is None:
            raise ValueError(f"missing target metadata for cluster {ident}")
        cov20 = [float(value) for value in row20["target_coverage_by_bin"].split(",")]
        cov60 = [float(value) for value in row60["target_coverage_by_bin"].split(",")]
        if len(cov20) != len(PROFILE_BINS) or len(cov60) != len(PROFILE_BINS):
            raise ValueError("unexpected radial coverage vector length")
        for idx, radial_bin in enumerate(PROFILE_BINS):
            key20_bin = (ident, band, member, radial_bin)
            key60_bin = (ident, band, member, radial_bin)
            records.append({
                "cluster_id": ident, "filter": band, "member_index": member,
                "cluster_radius_kpc": float(row20["cluster_radius_kpc"]),
                "redshift": float(frame["z"]), "richness": float(frame["richness"]),
                "satellite_bin_kpc": radial_bin,
                "target_coverage_20": cov20[idx], "target_coverage_60": cov60[idx],
                "target_pass_20": cov20[idx] >= .2,
                "target_pass_60": cov60[idx] >= .6,
                "retained_20": key20_bin in profile_keys["coverage_20"],
                "retained_60": key60_bin in profile_keys["coverage_60"],
                "n_valid_control_profiles": int(row20["valid_control_profiles"]),
            })

    output_rows = []
    for band in ("i", "i2"):
        band_rows = [row for row in records if row["filter"] == band]
        for satellite_bin in PROFILE_BINS:
            for rlo, rhi in CLUSTER_RADIUS_BINS:
                subset = [row for row in band_rows if row["satellite_bin_kpc"] == satellite_bin
                          and rlo <= row["cluster_radius_kpc"] < rhi]
                if not subset:
                    continue
                output_rows.append({
                    "filter": band, "satellite_bin_kpc": satellite_bin,
                    "clustercentric_radius_kpc": f"{rlo}-{rhi}",
                    "n_candidate_positions": len(subset),
                    "n_target_coverage_ge_20pct": sum(row["target_pass_20"] for row in subset),
                    "n_retained_20pct_with_controls": sum(row["retained_20"] for row in subset),
                    "n_target_coverage_ge_60pct": sum(row["target_pass_60"] for row in subset),
                    "n_retained_60pct_with_controls": sum(row["retained_60"] for row in subset),
                    "fraction_retained_20pct": float(np.mean([row["retained_20"] for row in subset])),
                    "median_target_coverage": float(np.median([row["target_coverage_20"] for row in subset])),
                    "p95_target_coverage": float(np.quantile([row["target_coverage_20"] for row in subset], .95)),
                })

    overall = {}
    for band in ("i", "i2"):
        selected = [row for row in records if row["filter"] == band and row["satellite_bin_kpc"] == "50-75"]
        cluster_ids = sorted({row["cluster_id"] for row in selected}, key=int)
        by_cluster = []
        for ident in cluster_ids:
            cluster = [row for row in selected if row["cluster_id"] == ident]
            frame = frame_by_id[ident]
            by_cluster.append({
                "cluster_id": ident, "filter": band, "redshift": float(frame["z"]),
                "richness": float(frame["richness"]), "n_forced_candidates": len(cluster),
                "n_retained_20pct": sum(row["retained_20"] for row in cluster),
                "any_retained_20pct": any(row["retained_20"] for row in cluster),
                "max_target_coverage_50_75": max(row["target_coverage_20"] for row in cluster),
            })
        passing = [row for row in by_cluster if row["any_retained_20pct"]]
        overall[band] = {
            "n_candidate_positions_50_75": len(selected),
            "candidate_target_coverage_20pct": summarize([row["target_coverage_20"] for row in selected]),
            "n_clusters_with_candidate": len(by_cluster),
            "n_clusters_with_retained_profile": len(passing),
            "cluster_retention_fraction": len(passing) / len(by_cluster) if by_cluster else None,
            "candidate_cluster_radius_distribution_kpc": summarize([row["cluster_radius_kpc"] for row in selected]),
            "retained_cluster_radius_distribution_kpc": summarize([row["cluster_radius_kpc"] for row in selected if row["retained_20"]]),
            "candidate_redshift_distribution": summarize([row["redshift"] for row in by_cluster]),
            "retained_redshift_distribution": summarize([row["redshift"] for row in passing]),
            "candidate_richness_distribution": summarize([row["richness"] for row in by_cluster]),
            "retained_richness_distribution": summarize([row["richness"] for row in passing]),
            "cluster_rows": by_cluster,
        }

    output.mkdir(parents=True, exist_ok=False)
    write_tsv(output / "coverage_selection_by_radius.tsv", output_rows)
    cluster_rows = [row for band in ("i", "i2") for row in overall[band]["cluster_rows"]]
    write_tsv(output / "cluster_retention_50_75.tsv", cluster_rows)
    result = {
        "schema": "astrafek-icl-forced-position-coverage-selection-v1",
        "status": "complete_exploratory", "science_measurement": False,
        "n_candidate_positions": 3197, "filters": overall,
        "mask_and_threshold_interpretation": "descriptive support audit; no causal or missing-at-random assumption is tested",
        "limitations": [
            "Selection is measured among CAMIRA positions without nearby SEP detections, not among confirmed satellite galaxies.",
            "The coverage floor is a post-hoc analysis setting; multiple threshold choices are not familywise adjusted.",
            "Only the target coverage vector is available for every candidate; control-coverage acceptance is observed only through retained profile rows.",
            "Clustercentric and survey-region correlations are not modeled in these descriptive rates.",
            "No flux, completeness, PSF, sky-transfer, source-Poisson, or total covariance correction is made.",
        ],
        "input_sha256": {str(path): sha256(path) for run_path in RUNS.values()
                          for path in (run_path / "summary.json", run_path / "candidate_diagnostics.tsv",
                                       run_path / "forced_satellite_profiles.tsv")}
                        | {str(TARGETS): sha256(TARGETS)},
        "output_sha256": {name: sha256(output / name) for name in
                          ("coverage_selection_by_radius.tsv", "cluster_retention_50_75.tsv")},
    }
    (output / "summary.json").write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    result = run(args.output)
    print(json.dumps({"status": result["status"], "filters": {
        band: {key: value for key, value in metrics.items() if key != "cluster_rows"}
        for band, metrics in result["filters"].items()},
        "science_measurement": result["science_measurement"]}, indent=2))


if __name__ == "__main__":
    main()
