#!/usr/bin/env python3
"""Block-bootstrap the exploratory CAMIRA-density / diffuse-light association.

This uses the existing K2.5, 40%-coverage, 50--300 kpc full-sample support
selection. It is a spatial robustness check, not a satellite outer-light
measurement or a science-approved ICL result.
"""
import argparse
import csv
import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.stats import pearsonr, spearmanr

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from evaluate_prestack_pair_support import pair_block_labels

BIN_NAMES = ("50-100", "100-150", "150-200", "200-300")
FILTERS = ("i", "i2")


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
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def partial_pearson(x, y, richness, redshift):
    covariates = np.column_stack((np.asarray(richness, float), np.asarray(redshift, float)))
    scale = np.std(covariates, axis=0)
    if np.any(~np.isfinite(covariates)) or np.any(scale <= 0):
        return np.nan
    design = np.column_stack((np.ones(len(x)), (covariates - np.mean(covariates, axis=0)) / scale))
    xr = np.asarray(x, float) - design @ np.linalg.lstsq(design, x, rcond=None)[0]
    yr = np.asarray(y, float) - design @ np.linalg.lstsq(design, y, rcond=None)[0]
    if np.std(xr) <= 0 or np.std(yr) <= 0:
        return np.nan
    return float(pearsonr(xr, yr).statistic)


def percentile(values, q):
    finite = np.asarray(values, float)
    finite = finite[np.isfinite(finite)]
    return float(np.percentile(finite, q)) if len(finite) else None


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13"))
    parser.add_argument("--stack", type=Path,
                        default=Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/prestack_flux_profiles_2026-10-18_02"))
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--replicas", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=20261019)
    args = parser.parse_args()
    data, stack, work = args.data.resolve(), args.stack.resolve(), args.work.resolve()
    if work.exists():
        parser.error(f"Refusing to overwrite existing output: {work}")
    if args.replicas < 100:
        parser.error("replicas must be at least 100")

    sample_path = data / "sample_manifest.json"
    stack_summary_path = stack / "summary.json"
    profiles_path = stack / "cluster_profiles.tsv"
    pairs_path = stack / "pair_metrics_k2p5_c40_r50_300.tsv"
    control_summary_path = data / "prestack_extraction_2026-10-15_01/extraction_summary.json"
    required = (sample_path, stack_summary_path, profiles_path, pairs_path, control_summary_path)
    absent = [str(path) for path in required if not path.is_file()]
    if absent:
        parser.error("Missing input(s): " + ", ".join(absent))

    sample = read_json(sample_path)
    stack_summary = read_json(stack_summary_path)
    control_summary = read_json(control_summary_path)
    if stack_summary.get("status") != "complete" or stack_summary.get("science_gates", {}).get("science_measurement") is not False:
        parser.error("Expected a complete exploratory stack with science_measurement=false")
    output_hashes = stack_summary.get("output_sha256", {})
    for path in (profiles_path, pairs_path):
        if output_hashes.get(str(path)) != sha256(path):
            parser.error(f"Stack output checksum mismatch: {path}")
    source_hashes = stack_summary.get("source_sha256", {})
    if source_hashes.get(str(sample_path)) != sha256(sample_path):
        parser.error("Sample manifest checksum mismatch")
    if stack_summary.get("n_targets_measured") != 800 or stack_summary.get("n_matched_controls") != 798:
        parser.error("Unexpected target/control counts")

    profiles = {str(row["id"]): row for row in read_tsv(profiles_path)}
    pairs = read_tsv(pairs_path)
    controls = {str(row["target_cluster_id"]): row for row in control_summary.get("records", [])}
    if len(profiles) != 800 or len(pairs) != 137:
        parser.error("Unexpected full-profile or selected-pair table size")
    if {band: sum(row["filter"] == band for row in pairs) for band in FILTERS} != {"i": 44, "i2": 93}:
        parser.error("Unexpected exact-filter support in the saved selection")

    result_rows = []
    filter_diagnostics = {}
    for filter_index, band in enumerate(FILTERS):
        selected = [row for row in pairs if row["filter"] == band]
        n = len(selected)
        if n < 8:
            filter_diagnostics[band] = {"n_pairs": n, "status": "insufficient_pairs"}
            continue
        target_rows, control_rows = [], []
        for row in selected:
            target = profiles.get(str(row["id"]))
            control = controls.get(str(row["id"]))
            if target is None or control is None or str(control.get("control_id")) != str(row["control_id"]):
                parser.error(f"Target/control provenance mismatch for {row['id']}")
            target_rows.append({"ra_deg": float(target["ra_deg"]), "dec_deg": float(target["dec_deg"])})
            control_rows.append({"ra_deg": float(control["ra_deg"]), "dec_deg": float(control["dec_deg"])})
        _, labels, n_target_cells, n_control_cells, n_blocks, largest = pair_block_labels(
            target_rows, control_rows, np.ones(n, dtype=bool))
        if band == "i" and n_blocks != 15:
            parser.error(f"Expected 15 joint 5-degree blocks for i, found {n_blocks}")
        if band == "i2" and n_blocks != 37:
            parser.error(f"Expected 37 joint 5-degree blocks for i2, found {n_blocks}")
        group_ids = np.unique(labels)
        members = [np.flatnonzero(labels == group) for group in group_ids]
        rng = np.random.default_rng(args.seed + filter_index * 1000)
        by_bin = []
        for bin_index, bin_name in enumerate(BIN_NAMES):
            # The table stores each annular count in a separate field.
            density_fields = ("density_mpc2_50_100", "density_mpc2_100_150",
                              "density_mpc2_150_200", "density_mpc2_200_300")
            x = np.asarray([float(row[density_fields[bin_index]]) for row in selected], float)
            y = np.asarray([float(row[f"delta_sb_{bin_name.replace('-', '_')}"]) for row in selected], float)
            richness = np.asarray([float(row["richness"]) for row in selected], float)
            redshift = np.asarray([float(row["z"]) for row in selected], float)
            rho = float(spearmanr(x, y).statistic)
            p_nominal = float(spearmanr(x, y).pvalue)
            partial = partial_pearson(x, y, richness, redshift)
            rho_boot = np.full(args.replicas, np.nan)
            partial_boot = np.full(args.replicas, np.nan)
            for replica in range(args.replicas):
                draw = rng.integers(0, len(members), size=len(members))
                indices = np.concatenate([members[int(i)] for i in draw])
                bx, by = x[indices], y[indices]
                if len(indices) >= 5 and np.std(bx) > 0 and np.std(by) > 0:
                    rho_boot[replica] = spearmanr(bx, by).statistic
                    partial_boot[replica] = partial_pearson(
                        bx, by, richness[indices], redshift[indices])
            row_result = {
                "filter": band, "radial_bin_kpc": bin_name, "n_pairs": n,
                "n_joint_5deg_blocks": n_blocks, "largest_joint_block": largest,
                "n_target_cells": n_target_cells, "n_control_cells": n_control_cells,
                "mean_satellite_density_mpc2": float(np.mean(x)),
                "mean_delta_sb_njy_arcsec2": float(np.mean(y)),
                "spearman_rho": rho, "spearman_p_nominal": p_nominal,
                "spearman_block_bootstrap_ci95_low": percentile(rho_boot, 2.5),
                "spearman_block_bootstrap_ci95_high": percentile(rho_boot, 97.5),
                "spearman_block_bootstrap_fraction_positive": float(np.mean(rho_boot[np.isfinite(rho_boot)] > 0)),
                "partial_pearson_r_control_richness_z": partial if math.isfinite(partial) else None,
                "partial_block_bootstrap_ci95_low": percentile(partial_boot, 2.5),
                "partial_block_bootstrap_ci95_high": percentile(partial_boot, 97.5),
                "partial_block_bootstrap_fraction_positive": float(np.mean(partial_boot[np.isfinite(partial_boot)] > 0)),
                "valid_spearman_replicas": int(np.isfinite(rho_boot).sum()),
                "valid_partial_replicas": int(np.isfinite(partial_boot).sum()),
            }
            by_bin.append(row_result)
            result_rows.append(row_result)
        filter_diagnostics[band] = {
            "n_pairs": n, "n_joint_5deg_blocks": n_blocks,
            "largest_joint_block": largest, "n_target_cells": n_target_cells,
            "n_control_cells": n_control_cells,
            "center_masked_k2p5": sum(row["center_masked_k2p5"] == "True" for row in selected),
            "automatic_unique_center_candidates": sum(row["auto_unique_center_candidate"] == "True" for row in selected),
            "status": "exploratory_block_bootstrap",
        }

    work.mkdir(parents=True, exist_ok=False)
    fields = list(result_rows[0].keys()) if result_rows else ["filter", "radial_bin_kpc", "n_pairs"]
    write_tsv(work / "block_bootstrap_metrics.tsv", result_rows, fields)
    summary = {
        "schema": "hsc-full-sample-satellite-association-block-bootstrap-v1",
        "status": "complete", "science_measurement": False,
        "selection": "Kron 2.5; common target/control coverage >=40%; 50-300 kpc; exact filter matches",
        "n_pairs": len(pairs), "filters": filter_diagnostics,
        "bootstrap": {"method": "resample joint connected 5-degree target/control sky blocks with replacement",
                      "replicas": args.replicas, "seed": args.seed,
                      "interval": "central percentile 95%; fraction_positive is stability, not a null-test p-value"},
        "limitations": [
            "CAMIRA member positions are counts only; no member photometry, stellar mass, membership probability or satellite-centered outer-light profile is used.",
            "The target-minus-control diffuse profile is exploratory and does not pass the sky, PSF, extinction or mask science gates.",
            "Block bootstrap addresses spatial grouping at the chosen 5-degree scale but does not propagate measurement/injection uncertainty or correct the four-bin multiple-testing problem.",
            "The data-driven pair selection and low usable sample size can bias the association; results must not be reported as a physical detection.",
        ],
        "source_sha256": {str(path): sha256(path) for path in required},
        "input_stack_summary": str(stack_summary_path),
        "output_sha256": {str(work / "block_bootstrap_metrics.tsv"): sha256(work / "block_bootstrap_metrics.tsv")},
        "script_sha256": {str(Path(__file__).resolve()): sha256(__file__)},
    }
    (work / "summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
                                        encoding="utf-8")
    print(json.dumps({"status": "complete", "pairs": len(pairs), "filters": filter_diagnostics,
                      "output": str(work)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
