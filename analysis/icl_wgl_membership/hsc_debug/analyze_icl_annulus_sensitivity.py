#!/usr/bin/env python3
"""Measure full-sample ICL annulus sensitivity to mask and area cuts.

Uses existing target and matched-control profiles only. It varies Kron mask
factor and common target/control coverage floor independently for the
200--300 kpc annulus, keeping filters separate and recording every selected
pair. This is a sensitivity analysis, not a science-approved measurement.
"""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import platform
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from evaluate_prestack_pair_support import pair_block_labels
from measure_icl_annulus_support import block_bootstrap_mean


FACTORS = (2.5, 3.5, 5.0)
THRESHOLDS = (0.25, 0.35, 0.40, 0.50)
ANNULUS_INDEX = 3
ANNULUS_NAME = "200-300"
PROFILE_BIN_COUNT = 6


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def read_tsv(path):
    with Path(path).open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def write_tsv(path, rows, fields):
    with Path(path).open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t",
                                extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def pair_is_eligible(target_row, target_coverage, target_raw, control_record,
                     control_coverage, control_raw, factor_index, threshold):
    if target_row["filter"] != str(control_record["header_filter"]):
        return False
    target_cov = float(target_coverage[factor_index, ANNULUS_INDEX])
    control_cov = float(control_coverage[factor_index, ANNULUS_INDEX])
    target_flux = float(target_raw[factor_index, ANNULUS_INDEX])
    control_flux = float(control_raw[factor_index, ANNULUS_INDEX])
    return (target_cov >= threshold and control_cov >= threshold
            and math.isfinite(target_flux) and math.isfinite(control_flux))


def evaluate_config(target_rows, target_raw, target_coverage, control_by_target,
                    control_profiles, factor_index, factor, threshold, band,
                    replicas, seed):
    selected = []
    for index, target in enumerate(target_rows):
        ident = str(target["id"])
        control_record = control_by_target.get(ident)
        if control_record is None or target["filter"] != band:
            continue
        control_id = str(control_record["control_id"])
        profile = control_profiles.get(control_id)
        if profile is None:
            raise ValueError(f"Missing verified control profile: {control_id}")
        control_raw, control_coverage = profile
        if not pair_is_eligible(
                target, target_coverage[index], target_raw[index], control_record,
                control_coverage, control_raw, factor_index, threshold):
            continue
        target_flux = float(target_raw[index, factor_index, ANNULUS_INDEX])
        control_flux = float(control_raw[factor_index, ANNULUS_INDEX])
        common_coverage = min(
            float(target_coverage[index, factor_index, ANNULUS_INDEX]),
            float(control_coverage[factor_index, ANNULUS_INDEX]))
        selected.append({
            "id": ident, "filter": band, "control_id": control_id,
            "target_flux": target_flux, "control_flux": control_flux,
            "delta_flux": target_flux - control_flux,
            "target_coverage": float(target_coverage[index, factor_index, ANNULUS_INDEX]),
            "control_coverage": float(control_coverage[factor_index, ANNULUS_INDEX]),
            "common_coverage": common_coverage,
            "target_ra_deg": float(target["ra_deg"]),
            "target_dec_deg": float(target["dec_deg"]),
            "control_ra_deg": float(control_record["ra_deg"]),
            "control_dec_deg": float(control_record["dec_deg"]),
            "z": float(target["z"]), "richness": float(target["richness"]),
            "auto_unique_center_candidate": str(target["auto_unique_center_candidate"]),
            "center_masked_by_kron": str(bool(json.loads(target["center_masked_by_kron"]).get(
                str(factor), False))),
        })

    if selected:
        target_block_rows = [{"ra_deg": row["target_ra_deg"], "dec_deg": row["target_dec_deg"]}
                             for row in selected]
        control_block_rows = [{"ra_deg": row["control_ra_deg"], "dec_deg": row["control_dec_deg"]}
                              for row in selected]
        indices, labels, target_cells, control_cells, n_blocks, largest = pair_block_labels(
            target_block_rows, control_block_rows, np.ones(len(selected), dtype=bool))
        if len(indices) != len(selected):
            raise ValueError("Spatial block labels changed pair order or count")
        values = np.asarray([row["delta_flux"] for row in selected], dtype=float)
        mean_delta = float(values.mean())
        mean_target = float(np.mean([row["target_flux"] for row in selected]))
        mean_control = float(np.mean([row["control_flux"] for row in selected]))
        low, high = block_bootstrap_mean(values, labels, replicas, seed)
    else:
        target_cells = control_cells = n_blocks = largest = 0
        mean_delta = mean_target = mean_control = low = high = None

    summary = {
        "kron_factor": factor, "coverage_floor": threshold, "radial_bin_kpc": ANNULUS_NAME,
        "filter": band, "n_pairs": len(selected),
        "n_target_sky_cells": target_cells, "n_control_sky_cells": control_cells,
        "n_joint_spatial_blocks": n_blocks, "largest_joint_block": largest,
        "mean_target_sb_njy_arcsec2": mean_target,
        "mean_control_sb_njy_arcsec2": mean_control,
        "mean_pair_delta_njy_arcsec2": mean_delta,
        "spatial_block_bootstrap_95_low": low,
        "spatial_block_bootstrap_95_high": high,
    }
    detail = [{
        **row, "kron_factor": factor, "coverage_floor": threshold,
        "radial_bin_kpc": ANNULUS_NAME,
    } for row in selected]
    return summary, detail


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", type=Path,
                        default=Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13"))
    parser.add_argument("--profiles", type=Path)
    parser.add_argument("--controls", type=Path)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--bootstrap-replicas", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=20261018)
    args = parser.parse_args()
    if not platform.node().lower().startswith("lageunha"):
        parser.error("Run data-profile analysis on lageunha")
    if args.bootstrap_replicas < 100:
        parser.error("bootstrap replicas must be at least 100")

    sample = args.sample.resolve()
    profiles_dir = (args.profiles or sample / "prestack_flux_profiles_2026-10-18_02").resolve()
    controls_dir = (args.controls or sample / "prestack_extraction_2026-10-15_01").resolve()
    work = args.work.resolve()
    if work.exists():
        parser.error(f"Refusing to overwrite existing output: {work}")
    target_table = profiles_dir / "cluster_profiles.tsv"
    target_archive = profiles_dir / "cluster_profiles.npz"
    profiles_summary_path = profiles_dir / "summary.json"
    extraction_summary_path = controls_dir / "extraction_summary.json"
    manifest_path = sample / "sample_manifest.json"
    required = (target_table, target_archive, profiles_summary_path,
                extraction_summary_path, manifest_path)
    if any(not path.is_file() for path in required):
        parser.error("A required existing profile, control, or manifest input is missing")

    profiles_summary = read_json(profiles_summary_path)
    extraction_summary = read_json(extraction_summary_path)
    manifest = read_json(manifest_path)
    target_rows = read_tsv(target_table)
    control_records = extraction_summary.get("records", [])
    if (len(manifest.get("clusters", [])) != 800 or len(target_rows) != 800
            or profiles_summary.get("status") != "complete"
            or profiles_summary.get("n_targets_measured") != 800
            or profiles_summary.get("science_gates", {}).get("science_measurement") is not False
            or extraction_summary.get("status") != "complete"
            or extraction_summary.get("science_measurement") is not False
            or len(control_records) != 798):
        parser.error("Expected complete exploratory 800-target/798-control inputs")
    for path in (target_table, target_archive):
        if profiles_summary.get("output_sha256", {}).get(str(path)) != sha256(path):
            parser.error(f"Target profile checksum mismatch: {path}")
    if [str(row["id"]) for row in target_rows] != [str(row["id"])
                                                   for row in manifest["clusters"]]:
        parser.error("Target profile row order differs from sample manifest")

    control_by_target = {}
    for record in control_records:
        target_id = str(record["target_cluster_id"])
        if target_id in control_by_target:
            parser.error(f"Duplicate control assignment for {target_id}")
        control_by_target[target_id] = record
    if len(control_by_target) != 798:
        parser.error("Expected 798 unique matched target IDs")

    target_raw_archive = np.load(target_archive)
    target_raw = target_raw_archive["raw"].copy()
    target_coverage = target_raw_archive["coverage"].copy()
    target_raw_archive.close()
    if target_raw.shape != (800, 3, PROFILE_BIN_COUNT) or target_coverage.shape != target_raw.shape:
        parser.error("Unexpected target profile cube shape")

    extraction_by_id = {str(row["control_id"]): row
                        for row in extraction_summary["records"]}
    control_profiles = {}
    control_hashes = {}
    provenance = profiles_summary.get("source_sha256", {})
    for target_id, record in control_by_target.items():
        if str(record.get("target_cluster_id")) != target_id:
            raise ValueError(f"Control target mismatch: {record.get('control_id')}")
        control_id = str(record["control_id"])
        baseline_path = controls_dir / control_id / "baseline.npz"
        expected_hash = record.get("baseline_sha256")
        if (not baseline_path.is_file() or not expected_hash
                or sha256(baseline_path) != expected_hash
                or provenance.get(str(baseline_path)) != expected_hash):
            raise ValueError(f"Control baseline checksum/provenance mismatch: {control_id}")
        if control_id not in extraction_by_id:
            raise ValueError(f"Control missing from extraction summary: {control_id}")
        with np.load(baseline_path) as baseline:
            raw = np.asarray(baseline["raw"], dtype=float)
            coverage = np.asarray(baseline["coverage"], dtype=float)
        if raw.shape != (3, PROFILE_BIN_COUNT) or coverage.shape != raw.shape:
            raise ValueError(f"Unexpected control profile shape: {control_id}")
        control_profiles[control_id] = (raw, coverage)
        control_hashes[str(baseline_path)] = expected_hash

    summary_rows, detail_rows = [], []
    for factor_index, factor in enumerate(FACTORS):
        for threshold_index, threshold in enumerate(THRESHOLDS):
            for band_index, band in enumerate(("i", "i2")):
                result, detail = evaluate_config(
                    target_rows, target_raw, target_coverage, control_by_target,
                    control_profiles, factor_index, factor, threshold, band,
                    args.bootstrap_replicas,
                    args.seed + factor_index * 10000 + threshold_index * 1000 + band_index)
                summary_rows.append(result)
                detail_rows.extend(detail)

    work.mkdir(parents=True, exist_ok=False)
    summary_out = work / "icl_annulus_sensitivity.tsv"
    details_out = work / "selected_pair_profiles.tsv"
    write_tsv(summary_out, summary_rows, list(summary_rows[0]))
    detail_fields = ["kron_factor", "coverage_floor", "radial_bin_kpc", "id", "filter",
                     "control_id", "z", "richness", "target_flux", "control_flux",
                     "delta_flux", "target_coverage", "control_coverage", "common_coverage",
                     "target_ra_deg", "target_dec_deg", "control_ra_deg", "control_dec_deg",
                     "auto_unique_center_candidate", "center_masked_by_kron"]
    detail_rows = [{key: row[key] for key in detail_fields} for row in detail_rows]
    write_tsv(details_out, detail_rows, detail_fields)
    outputs = (summary_out, details_out)
    summary = {
        "schema": "hsc-icl-annulus-selection-sensitivity-v1", "status": "complete",
        "science_measurement": False,
        "method": {
            "radial_bin_kpc": ANNULUS_NAME,
            "selection": "independent target and matched-control coverage thresholds; exact header filter match; finite raw profiles",
            "mask_factors": list(FACTORS), "coverage_floors": list(THRESHOLDS),
            "weighting": "equal target-control pairs; filters kept separate",
            "uncertainty": "5-degree connected target/control sky-block bootstrap; exploratory",
        },
        "n_targets": 800, "n_matched_controls": 798,
        "n_configurations": len(summary_rows), "bootstrap_replicas": args.bootstrap_replicas,
        "seed": args.seed,
        "input_sha256": {
            str(manifest_path): sha256(manifest_path), str(target_table): sha256(target_table),
            str(target_archive): sha256(target_archive),
            str(profiles_summary_path): sha256(profiles_summary_path),
            str(extraction_summary_path): sha256(extraction_summary_path),
            "control_baseline_hashes_sha256": hashlib.sha256(
                json.dumps(sorted(control_hashes.items()), separators=(",", ":")).encode()
            ).hexdigest(),
        },
        "program_sha256": sha256(Path(__file__).resolve()),
        "output_sha256": {str(path): sha256(path) for path in outputs},
        "limitations": [
            "Coverage and Kron-factor cuts are sensitivity cases; this analysis does not choose a science cut.",
            "The 200-300 kpc annulus is measured independently; other radial bins are not required.",
            "The bootstrap does not include pixel-level coadd covariance, PSF, sky-transfer, extinction, or source-mask systematics.",
            "All outputs remain exploratory and do not change science_measurement=false.",
        ],
    }
    (work / "summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False,
                                                   allow_nan=False) + "\n", encoding="utf-8")
    primary = [row for row in summary_rows if row["kron_factor"] == 2.5
               and row["coverage_floor"] == 0.40]
    print(json.dumps({"status": summary["status"], "science_measurement": False,
                      "n_configurations": len(summary_rows), "k2p5_c40": primary,
                      "output": str(work)}, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
