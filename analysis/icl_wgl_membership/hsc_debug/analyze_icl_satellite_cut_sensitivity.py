#!/usr/bin/env python3
"""Test ICL--satellite association stability across mask/area choices.

Joins each existing 200--300 kpc ICL sensitivity selection to the same-filter
satellite profile on the seven-bin common-cluster sample. Uncertainty resamples
connected five-degree blocks joining both target and control skies. No image or
catalogue rows are fetched; all outputs remain exploratory.
"""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import sys

import numpy as np
from scipy.stats import pearsonr, spearmanr

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from analyze_satellite_common_support import partial_pearson
from evaluate_prestack_pair_support import pair_block_labels


FACTORS = (2.5, 3.5, 5.0)
THRESHOLDS = (0.25, 0.35, 0.40, 0.50)


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


def connected_pair_blocks(rows):
    target = [{"ra_deg": float(row["target_ra_deg"]),
               "dec_deg": float(row["target_dec_deg"])} for row in rows]
    control = [{"ra_deg": float(row["control_ra_deg"]),
                "dec_deg": float(row["control_dec_deg"])} for row in rows]
    selected = np.ones(len(rows), dtype=bool)
    indices, labels, target_cells, control_cells, n_blocks, largest = pair_block_labels(
        target, control, selected)
    if len(indices) != len(rows):
        raise ValueError("Pair-block helper changed the selected row order")
    return labels, target_cells, control_cells, n_blocks, largest


def summarize_association(rows, factor, threshold, band, selection, replicas, seed):
    n = len(rows)
    base = {"kron_factor": factor, "coverage_floor": threshold, "filter": band,
            "selection": selection, "n_clusters": n}
    if n < 5:
        return {**base, "status": "insufficient_support", "n_joint_blocks": 0,
                "spearman_rho": None, "spearman_p_nominal": None,
                "spearman_block_bootstrap_95_low": None,
                "spearman_block_bootstrap_95_high": None,
                "partial_pearson_r_z_richness": None,
                "partial_block_bootstrap_95_low": None,
                "partial_block_bootstrap_95_high": None,
                "valid_bootstrap_replicas": 0}
    labels, target_cells, control_cells, n_blocks, largest = connected_pair_blocks(rows)
    if n_blocks < 2:
        return {**base, "status": "insufficient_spatial_blocks",
                "n_target_cells": target_cells, "n_control_cells": control_cells,
                "n_joint_blocks": n_blocks, "largest_joint_block": largest,
                "spearman_rho": None, "spearman_p_nominal": None,
                "spearman_block_bootstrap_95_low": None,
                "spearman_block_bootstrap_95_high": None,
                "partial_pearson_r_z_richness": None,
                "partial_block_bootstrap_95_low": None,
                "partial_block_bootstrap_95_high": None,
                "valid_bootstrap_replicas": 0}

    x = np.asarray([float(row["icl_delta_200_300"]) for row in rows])
    y = np.asarray([float(row["satellite_delta_50_75"]) for row in rows])
    z = np.asarray([float(row["z"]) for row in rows])
    richness = np.asarray([float(row["richness"]) for row in rows])
    rho, p_value = spearmanr(x, y)
    partial = partial_pearson(x, y, z, richness)
    unique = np.unique(labels)
    members = {label: np.flatnonzero(labels == label) for label in unique}
    rng = np.random.default_rng(seed)
    rho_draws, partial_draws = [], []
    for _ in range(replicas):
        selected_blocks = rng.choice(unique, size=len(unique), replace=True)
        indices = np.concatenate([members[label] for label in selected_blocks])
        bx, by = x[indices], y[indices]
        if np.std(bx) <= 0 or np.std(by) <= 0:
            continue
        rho_draws.append(float(spearmanr(bx, by).statistic))
        value = partial_pearson(bx, by, z[indices], richness[indices])
        if math.isfinite(value):
            partial_draws.append(value)

    def interval(values):
        return ([float(v) for v in np.percentile(values, (2.5, 97.5))]
                if values else [None, None])

    return {
        **base, "status": "exploratory", "n_target_cells": target_cells,
        "n_control_cells": control_cells, "n_joint_blocks": n_blocks,
        "largest_joint_block": largest, "mean_icl_delta_200_300_njy_arcsec2": float(x.mean()),
        "mean_satellite_delta_50_75_njy_arcsec2": float(y.mean()),
        "spearman_rho": float(rho), "spearman_p_nominal": float(p_value),
        "spearman_block_bootstrap_95_low": interval(rho_draws)[0],
        "spearman_block_bootstrap_95_high": interval(rho_draws)[1],
        "partial_pearson_r_z_richness": partial if math.isfinite(partial) else None,
        "partial_block_bootstrap_95_low": interval(partial_draws)[0],
        "partial_block_bootstrap_95_high": interval(partial_draws)[1],
        "valid_bootstrap_replicas": len(rho_draws),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--icl-sensitivity", type=Path, required=True)
    parser.add_argument("--satellite", type=Path, required=True,
                        help="Existing per-annulus satellite output")
    parser.add_argument("--common-support", type=Path, required=True,
                        help="Verified seven-bin common-cluster output")
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--replicas", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=20261019)
    args = parser.parse_args()
    if args.replicas < 100:
        parser.error("replicas must be at least 100")
    icl_dir, satellite, common, work = (
        args.icl_sensitivity.resolve(), args.satellite.resolve(),
        args.common_support.resolve(), args.work.resolve())
    if work.exists():
        parser.error(f"Refusing to overwrite existing output: {work}")

    icl_summary_path = icl_dir / "summary.json"
    icl_detail_path = icl_dir / "selected_pair_profiles.tsv"
    sat_summary_path = satellite / "summary.json"
    sat_profiles_path = satellite / "cluster_profiles.tsv"
    common_summary_path = common / "summary.json"
    common_ids_path = common / "common_cluster_ids.tsv"
    if any(not path.is_file() for path in (icl_summary_path, icl_detail_path,
                                           sat_summary_path, sat_profiles_path,
                                           common_summary_path, common_ids_path)):
        parser.error("A required sensitivity, satellite, or common-sample input is missing")
    icl_summary, sat_summary, common_summary = (
        read_json(icl_summary_path), read_json(sat_summary_path), read_json(common_summary_path))
    if (icl_summary.get("status") != "complete"
            or icl_summary.get("science_measurement") is not False
            or icl_summary.get("schema") != "hsc-icl-annulus-selection-sensitivity-v1"
            or sat_summary.get("status") != "complete"
            or sat_summary.get("science_measurement") is not False
            or sat_summary.get("n_frames_completed") != 800
            or common_summary.get("status") != "complete"
            or common_summary.get("science_measurement") is not False):
        parser.error("Expected complete exploratory 800-frame inputs")
    for path, summary in ((icl_detail_path, icl_summary),
                          (sat_profiles_path, sat_summary),
                          (common_ids_path, common_summary)):
        if summary.get("output_sha256", {}).get(str(path)) != sha256(path):
            parser.error(f"Input checksum mismatch: {path}")

    common_rows = read_tsv(common_ids_path)
    common_ids = {band: set() for band in ("i", "i2")}
    common_info = {}
    for row in common_rows:
        band, ident = row["filter"], str(row["cluster_id"])
        if band not in common_ids or ident in common_ids[band]:
            parser.error(f"Unexpected or duplicate common-cluster record: {band}/{ident}")
        common_ids[band].add(ident)
        common_info[(ident, band)] = row

    satellite_rows = read_tsv(sat_profiles_path)
    satellite_outer = {}
    for row in satellite_rows:
        if row["radial_bin_kpc"] != "50-75":
            continue
        key = (str(row["cluster_id"]), row["filter"])
        if key in satellite_outer:
            parser.error(f"Duplicate satellite outer profile: {key}")
        satellite_outer[key] = row

    icl_details = read_tsv(icl_detail_path)
    by_config = {}
    for factor in FACTORS:
        for threshold in THRESHOLDS:
            for band in ("i", "i2"):
                by_config[(factor, threshold, band)] = []
    seen = set()
    for row in icl_details:
        factor, threshold, band = (float(row["kron_factor"]),
                                  float(row["coverage_floor"]), row["filter"])
        key = (factor, threshold, band, str(row["id"]))
        if key in seen:
            parser.error(f"Duplicate ICL sensitivity pair: {key}")
        seen.add(key)
        if (factor, threshold, band) not in by_config:
            parser.error(f"Unexpected ICL sensitivity configuration: {key[:3]}")
        ident = str(row["id"])
        if ident not in common_ids[band]:
            continue
        satellite_profile = satellite_outer.get((ident, band))
        if satellite_profile is None:
            parser.error(f"Missing common-sample satellite outer profile: {ident}/{band}")
        common_record = common_info[(ident, band)]
        for field, common_field in (("z", "redshift"), ("richness", "richness")):
            if not math.isclose(float(row[field]), float(common_record[common_field]),
                                rel_tol=0, abs_tol=1e-8):
                parser.error(f"ICL/satellite metadata mismatch: {ident}/{band}/{field}")
        joined = {
            "cluster_id": ident, "filter": band, "kron_factor": factor,
            "coverage_floor": threshold,
            "icl_delta_200_300": float(row["delta_flux"]),
            "satellite_delta_50_75": float(satellite_profile["mean_delta_sb_njy_arcsec2"]),
            "z": float(row["z"]), "richness": float(row["richness"]),
            "target_ra_deg": float(row["target_ra_deg"]),
            "target_dec_deg": float(row["target_dec_deg"]),
            "control_ra_deg": float(row["control_ra_deg"]),
            "control_dec_deg": float(row["control_dec_deg"]),
            "center_eligible": (row["auto_unique_center_candidate"] == "True"
                                and row["center_masked_by_kron"] == "False"),
        }
        by_config[(factor, threshold, band)].append(joined)

    result_rows = []
    matched_rows = []
    for factor_index, factor in enumerate(FACTORS):
        for threshold_index, threshold in enumerate(THRESHOLDS):
            for band_index, band in enumerate(("i", "i2")):
                all_rows = by_config[(factor, threshold, band)]
                matched_rows.extend(all_rows)
                result_rows.append(summarize_association(
                    all_rows, factor, threshold, band, "all_common_clusters", args.replicas,
                    args.seed + factor_index * 100000 + threshold_index * 1000 + band_index * 10))
                eligible = [row for row in all_rows if row["center_eligible"]]
                result_rows.append(summarize_association(
                    eligible, factor, threshold, band, "center_eligible", args.replicas,
                    args.seed + factor_index * 100000 + threshold_index * 1000 + band_index * 10 + 1))

    work.mkdir(parents=True, exist_ok=False)
    result_path = work / "association_sensitivity.tsv"
    detail_path = work / "matched_cluster_values.tsv"
    result_fields = sorted({key for row in result_rows for key in row})
    write_tsv(result_path, result_rows, result_fields)
    detail_fields = ["cluster_id", "filter", "kron_factor", "coverage_floor",
                     "icl_delta_200_300", "satellite_delta_50_75", "z", "richness",
                     "target_ra_deg", "target_dec_deg", "control_ra_deg", "control_dec_deg",
                     "center_eligible"]
    write_tsv(detail_path, matched_rows, detail_fields)
    summary = {
        "schema": "hsc-icl-satellite-cut-sensitivity-v1", "status": "complete",
        "science_measurement": False,
        "method": {
            "satellite_sample": "same-filter cluster intersection with at least one retained profile in every one of seven radial bins",
            "icl_variable": "annulus-specific 200-300 kpc target-control excess",
            "uncertainty": "connected 5-degree blocks link target and control sky cells; block bootstrap; exploratory",
            "covariates": "partial Pearson residualizes redshift and richness linearly",
            "comparisons": "all 24 Kron-factor/coverage configurations reported without selecting a preferred cut",
        },
        "n_common_satellite_clusters": {band: len(common_ids[band]) for band in ("i", "i2")},
        "n_configurations": 24, "n_result_groups": len(result_rows),
        "replicas": args.replicas, "seed": args.seed,
        "input_sha256": {
            str(icl_summary_path): sha256(icl_summary_path),
            str(icl_detail_path): sha256(icl_detail_path),
            str(sat_summary_path): sha256(sat_summary_path),
            str(sat_profiles_path): sha256(sat_profiles_path),
            str(common_summary_path): sha256(common_summary_path),
            str(common_ids_path): sha256(common_ids_path),
        },
        "program_sha256": sha256(Path(__file__).resolve()),
        "output_sha256": {str(path): sha256(path) for path in (result_path, detail_path)},
        "limitations": [
            "The common satellite cluster set is selected and differs by filter.",
            "All configurations share the same images and masks; this is sensitivity analysis, not an independent replication.",
            "Partial-correlation intervals and nominal p-values do not correct for testing multiple configurations.",
            "No calibrated PSF wings, sky-transfer correction, pixel covariance, source-mass weighting, or extinction correction is applied.",
            "science_measurement remains false.",
        ],
    }
    (work / "summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False,
                                                   allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({"status": summary["status"], "science_measurement": False,
                      "n_common_satellite_clusters": summary["n_common_satellite_clusters"],
                      "n_matched_rows": len(matched_rows),
                      "results": result_rows,
                      "output": str(work)}, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
