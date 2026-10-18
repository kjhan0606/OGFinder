#!/usr/bin/env python3
"""Re-stack existing ICL target/control profiles with annulus-specific support.

This reads previously measured 800-cluster profiles and their 798 matched
control baselines. It performs no image I/O, archive query, or download. The
same exploratory K2.5 mask and 40% target/control coverage floor are retained;
each radial annulus is selected independently so that a low-coverage neighbor
annulus does not discard a usable 200--300 kpc pair.
"""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import platform
import sys
import time

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
from evaluate_prestack_pair_support import pair_block_labels


BINS = ((50, 100), (100, 150), (150, 200), (200, 300), (300, 400), (400, 500))
BIN_NAMES = [f"{lo}-{hi}" for lo, hi in BINS]
K2P5_INDEX = 0


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
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


def block_bootstrap_mean(values, labels, replicas, seed):
    values = np.asarray(values, dtype=float)
    labels = np.asarray(labels, dtype=object)
    unique = np.unique(labels)
    if len(unique) < 2:
        return [None, None]
    members = {label: values[labels == label] for label in unique}
    rng = np.random.default_rng(seed)
    draws = np.empty(replicas, dtype=float)
    for index in range(replicas):
        selected_blocks = rng.choice(unique, size=len(unique), replace=True)
        sample = np.concatenate([members[label] for label in selected_blocks])
        draws[index] = np.mean(sample)
    return [float(value) for value in np.percentile(draws, (2.5, 97.5))]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", type=Path,
                        default=Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13"))
    parser.add_argument("--profiles", type=Path,
                        help="Existing full-sample ICL profile output directory")
    parser.add_argument("--controls", type=Path,
                        help="Existing matched-control extraction directory")
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--coverage-floor", type=float, default=.40)
    parser.add_argument("--bootstrap-replicas", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=20261018)
    args = parser.parse_args()
    if not platform.node().lower().startswith("lageunha"):
        parser.error("Run data-profile analysis on lageunha")
    if not 0 < args.coverage_floor <= 1 or args.bootstrap_replicas < 100:
        parser.error("Invalid coverage/bootstrap configuration")

    sample = args.sample.resolve()
    profiles_dir = (args.profiles or sample / "prestack_flux_profiles_2026-10-18_02").resolve()
    controls_dir = (args.controls or sample / "prestack_extraction_2026-10-15_01").resolve()
    work = args.work.resolve()
    if work.exists():
        parser.error(f"Refusing to overwrite existing output: {work}")
    sample_manifest = sample / "sample_manifest.json"
    target_table = profiles_dir / "cluster_profiles.tsv"
    target_archive = profiles_dir / "cluster_profiles.npz"
    profile_summary_path = profiles_dir / "summary.json"
    extraction_summary_path = controls_dir / "extraction_summary.json"
    required = (sample_manifest, target_table, target_archive, profile_summary_path,
                extraction_summary_path)
    if any(not path.is_file() for path in required):
        parser.error("A required manifest, target profile, or control summary is missing")
    manifest = read_json(sample_manifest)
    profile_summary = read_json(profile_summary_path)
    extraction_summary = read_json(extraction_summary_path)
    target_rows = read_tsv(target_table)
    controls = extraction_summary.get("records", [])
    if (len(manifest.get("clusters", [])) != 800 or len(target_rows) != 800
            or profile_summary.get("status") != "complete"
            or profile_summary.get("n_targets_measured") != 800
            or profile_summary.get("science_gates", {}).get("science_measurement") is not False
            or extraction_summary.get("status") != "complete"
            or len(controls) != 798
            or extraction_summary.get("science_measurement") is not False):
        parser.error("Expected complete exploratory 800-target/798-control inputs")
    expected_target_hash = profile_summary.get("output_sha256", {}).get(str(target_table))
    expected_archive_hash = profile_summary.get("output_sha256", {}).get(str(target_archive))
    if (not expected_target_hash or sha256(target_table) != expected_target_hash
            or not expected_archive_hash or sha256(target_archive) != expected_archive_hash):
        parser.error("Target profile table/archive checksum mismatch")
    if [str(row["id"]) for row in target_rows] != [str(row["id"])
                                                      for row in manifest["clusters"]]:
        parser.error("Target profile row order differs from the sample manifest")

    control_by_target = {}
    for record in controls:
        target_id = str(record["target_cluster_id"])
        if target_id in control_by_target:
            parser.error(f"Duplicate control assignment for {target_id}")
        control_by_target[target_id] = record
    if len(control_by_target) != 798:
        parser.error("Control assignments do not have 798 unique target IDs")

    with np.load(target_archive) as archive:
        target_raw = archive["raw"].copy()
        target_coverage = archive["coverage"].copy()
    if target_raw.shape != (800, 3, 6) or target_coverage.shape != target_raw.shape:
        parser.error("Unexpected target profile cube shape")

    work.mkdir(parents=True, exist_ok=False)
    started = time.time()
    pair_rows_by_bin = {name: [] for name in BIN_NAMES}
    control_rows_for_blocks = []
    target_rows_for_blocks = []
    accepted_rows = []
    invalid_reasons = {"missing_control": 0, "filter_mismatch": 0}
    input_baseline_hashes = {}
    for index, target in enumerate(target_rows):
        ident = str(target["id"])
        control_record = control_by_target.get(ident)
        if control_record is None:
            invalid_reasons["missing_control"] += 1
            control_rows_for_blocks.append({"filter": "", "ra_deg": float(target["ra_deg"]),
                                            "dec_deg": float(target["dec_deg"])})
            target_rows_for_blocks.append({**target, "ra_deg": float(target["ra_deg"]),
                                           "dec_deg": float(target["dec_deg"])})
            continue
        control_id = str(control_record["control_id"])
        if str(control_record.get("target_cluster_id")) != ident:
            raise ValueError(f"Control target mismatch for {control_id}: expected {ident}")
        baseline_path = controls_dir / control_id / "baseline.npz"
        expected = control_record.get("baseline_sha256")
        if not baseline_path.is_file() or not expected or sha256(baseline_path) != expected:
            raise ValueError(f"Control baseline checksum mismatch: {control_id}")
        if profile_summary.get("source_sha256", {}).get(str(baseline_path)) != expected:
            raise ValueError(f"Stack provenance does not record this control: {control_id}")
        input_baseline_hashes[str(baseline_path)] = expected
        with np.load(baseline_path) as baseline:
            raw = np.asarray(baseline["raw"], dtype=float)
            coverage = np.asarray(baseline["coverage"], dtype=float)
        if raw.shape != (3, 6) or coverage.shape != (3, 6):
            raise ValueError(f"Unexpected control profile shape: {control_id}")
        control_filter = str(control_record["header_filter"])
        control_row = {"filter": control_filter, "ra_deg": float(control_record["ra_deg"]),
                       "dec_deg": float(control_record["dec_deg"]), "control_id": control_id}
        control_rows_for_blocks.append(control_row)
        target_rows_for_blocks.append({**target, "ra_deg": float(target["ra_deg"]),
                                       "dec_deg": float(target["dec_deg"])})
        if target["filter"] != control_filter:
            invalid_reasons["filter_mismatch"] += 1
            continue
        for bin_index, (lo, hi) in enumerate(BINS):
            common = min(float(target_coverage[index, K2P5_INDEX, bin_index]),
                         float(coverage[K2P5_INDEX, bin_index]))
            target_value = float(target_raw[index, K2P5_INDEX, bin_index])
            control_value = float(raw[K2P5_INDEX, bin_index])
            if common < args.coverage_floor or not (math.isfinite(target_value)
                                                     and math.isfinite(control_value)):
                continue
            delta = target_value - control_value
            frame = manifest["clusters"][index]
            members_key = ("members_%d_%d" % (lo, hi))
            n_members = int(float(target[members_key]))
            area_mpc2 = math.pi * (hi**2 - lo**2) / 1.0e6
            center_mask = json.loads(target["center_masked_by_kron"])
            row = {
                "id": ident, "filter": target["filter"], "control_id": control_id,
                "z": float(target["z"]), "richness": float(target["richness"]),
                "ra_deg": float(target["ra_deg"]), "dec_deg": float(target["dec_deg"]),
                "auto_unique_center_candidate": str(target["auto_unique_center_candidate"]),
                "center_masked_k2p5": str(bool(center_mask.get("2.5", False))),
                "radial_bin_kpc": f"{lo}-{hi}", "target_sb_njy_arcsec2": target_value,
                "control_sb_njy_arcsec2": control_value,
                "delta_sb_njy_arcsec2": delta,
                "target_coverage": float(target_coverage[index, K2P5_INDEX, bin_index]),
                "control_coverage": float(coverage[K2P5_INDEX, bin_index]),
                "common_coverage": common, "n_member_positions": n_members,
                "member_density_mpc2": n_members / area_mpc2,
            }
            pair_rows_by_bin[row["radial_bin_kpc"]].append(row)
            accepted_rows.append(row)

    summary_rows = []
    metrics_200_300 = []
    for bin_index, (lo, hi) in enumerate(BINS):
        name = f"{lo}-{hi}"
        selected_rows = pair_rows_by_bin[name]
        for band in ("i", "i2"):
            valid_ids = {row["id"] for row in selected_rows if row["filter"] == band}
            selected = np.asarray([str(row["id"]) in valid_ids and row["filter"] == band
                                   for row in target_rows], dtype=bool)
            # The pair rows already encode coverage, filter agreement and finite profiles.
            selected_count = int(np.count_nonzero(selected))
            if selected_count:
                selected_indices, labels, n_target_cells, n_control_cells, n_blocks, largest = pair_block_labels(
                    target_rows_for_blocks, control_rows_for_blocks, selected)
                selected_values = np.asarray([
                    row["delta_sb_njy_arcsec2"] for row in selected_rows
                    if row["filter"] == band], dtype=float)
                if len(selected_values) != len(selected_indices):
                    raise ValueError(f"Pair-block row count mismatch for {band}/{name}")
                low, high = block_bootstrap_mean(selected_values, labels,
                                                 args.bootstrap_replicas,
                                                 args.seed + bin_index * 100 + (0 if band == "i" else 1))
                mean = float(np.mean(selected_values))
            else:
                n_target_cells = n_control_cells = n_blocks = largest = 0
                low = high = mean = None
            summary_rows.append({
                "filter": band, "radial_bin_kpc": name, "n_pairs": selected_count,
                "n_target_sky_cells": n_target_cells, "n_control_sky_cells": n_control_cells,
                "n_joint_spatial_blocks": n_blocks, "largest_joint_block": largest,
                "mean_pair_delta_njy_arcsec2": mean,
                "spatial_block_bootstrap_95_low": low,
                "spatial_block_bootstrap_95_high": high,
            })

    for row in pair_rows_by_bin["200-300"]:
        metrics_200_300.append({
            "id": row["id"], "filter": row["filter"], "control_id": row["control_id"],
            "richness": row["richness"], "z": row["z"],
            "auto_unique_center_candidate": row["auto_unique_center_candidate"],
            "center_masked_k2p5": row["center_masked_k2p5"],
            "members_50_100": "", "members_100_150": "", "members_150_200": "",
            "members_200_300": row["n_member_positions"],
            "density_mpc2_50_100": "", "density_mpc2_100_150": "",
            "density_mpc2_150_200": "", "density_mpc2_200_300": row["member_density_mpc2"],
            "delta_sb_50_100": "", "delta_sb_100_150": "", "delta_sb_150_200": "",
            "delta_sb_200_300": row["delta_sb_njy_arcsec2"],
            "coverage_50_100": "", "coverage_100_150": "", "coverage_150_200": "",
            "coverage_200_300": row["common_coverage"],
        })

    profile_fields = ["id", "filter", "control_id", "z", "richness", "ra_deg", "dec_deg",
                      "auto_unique_center_candidate", "center_masked_k2p5", "radial_bin_kpc",
                      "target_sb_njy_arcsec2", "control_sb_njy_arcsec2", "delta_sb_njy_arcsec2",
                      "target_coverage", "control_coverage", "common_coverage",
                      "n_member_positions", "member_density_mpc2"]
    summary_fields = ["filter", "radial_bin_kpc", "n_pairs", "n_target_sky_cells",
                      "n_control_sky_cells", "n_joint_spatial_blocks", "largest_joint_block",
                      "mean_pair_delta_njy_arcsec2", "spatial_block_bootstrap_95_low",
                      "spatial_block_bootstrap_95_high"]
    pair_fields = ["id", "filter", "control_id", "richness", "z",
                   "auto_unique_center_candidate", "center_masked_k2p5",
                   "members_50_100", "members_100_150", "members_150_200", "members_200_300",
                   "density_mpc2_50_100", "density_mpc2_100_150", "density_mpc2_150_200",
                   "density_mpc2_200_300", "delta_sb_50_100", "delta_sb_100_150",
                   "delta_sb_150_200", "delta_sb_200_300", "coverage_50_100",
                   "coverage_100_150", "coverage_150_200", "coverage_200_300"]
    write_tsv(work / "annulus_pair_profiles.tsv", accepted_rows, profile_fields)
    write_tsv(work / "annulus_stack_summary.tsv", summary_rows, summary_fields)
    write_tsv(work / "pair_metrics_k2p5_c40_r50_300.tsv", metrics_200_300, pair_fields)
    outputs = [work / "annulus_pair_profiles.tsv", work / "annulus_stack_summary.tsv",
               work / "pair_metrics_k2p5_c40_r50_300.tsv"]
    summary = {
        "schema": "hsc-icl-per-annulus-support-v1", "status": "complete",
        "science_gates": {"science_measurement": False, "sky_correction_approved": False,
                          "psf_approved": False, "mask_ready_for_science": False},
        "n_target_profiles": len(target_rows), "n_control_profiles": len(controls),
        "n_pairs_by_annulus": {
            name: sum(item["radial_bin_kpc"] == name for item in accepted_rows)
            for name in BIN_NAMES
        },
        "per_filter_bin_support": summary_rows,
        "n_200_300_pairs": len(metrics_200_300),
        "invalid_reasons": invalid_reasons,
        "coverage_floor": args.coverage_floor, "kron_factor": 2.5,
        "extraction": "raw target-control profile; FLUXMAG0-converted nJy/arcsec^2",
        "support_rule": "target and matched control pass coverage independently in the selected radial annulus; filters must match exactly",
        "uncertainty": "5-degree connected target/control sky-block bootstrap for each filter and annulus; exploratory",
        "limitations": [
            "This is an annulus-specific support sensitivity analysis, not a science-approved ICL measurement.",
            "Different annuli use different target/control pairs; the resulting radial means are not a common-sample profile.",
            "The K2.5 mask and 40% coverage floor are exploratory choices and the automatic masks/BCG centers are not science-certified.",
            "PSF wings, upstream sky transfer, Galactic extinction, bandpass response, residual source light and cluster selection remain unvalidated.",
        ],
        "input_sha256": {
            str(sample_manifest): sha256(sample_manifest),
            str(target_table): sha256(target_table), str(target_archive): sha256(target_archive),
            str(profile_summary_path): sha256(profile_summary_path),
            str(extraction_summary_path): sha256(extraction_summary_path),
            "control_baseline_hashes_sha256": hashlib.sha256(
                json.dumps(sorted(input_baseline_hashes.items()), separators=(",", ":")).encode()).hexdigest(),
        },
        "program_sha256": sha256(Path(__file__).resolve()),
        "output_sha256": {str(path): sha256(path) for path in outputs},
        "elapsed_seconds": time.time() - started,
        "bootstrap": {"replicas": args.bootstrap_replicas, "seed": args.seed},
    }
    (work / "summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False,
                                                   allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({"status": summary["status"], "n_200_300_pairs": len(metrics_200_300),
                      "n_pairs_by_annulus": summary["n_pairs_by_annulus"],
                      "filters": {row["filter"] + "/" + row["radial_bin_kpc"]: row["n_pairs"]
                                  for row in summary_rows},
                      "output": str(work)}, ensure_ascii=False), flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
