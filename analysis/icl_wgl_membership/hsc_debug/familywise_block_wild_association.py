#!/usr/bin/env python3
"""Exploratory familywise inference for the ICL--satellite cut sensitivity.

Fits cluster-level satellite excess on ICL excess, redshift, and richness. A
restricted wild-cluster bootstrap imposes the no-association null. The same
Webb six-point weight is used for a spatial block across every cut/filter test,
then the maximum absolute studentized statistic supplies single-step max-T
adjusted p-values. The analysis is diagnostic only; it does not set a science
cut or turn on the science-measurement flag.
"""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from evaluate_prestack_pair_support import pair_block_labels


INPUT = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_cut_sensitivity_2026-10-18_01")
WEBB_WEIGHTS = np.asarray([
    -math.sqrt(1.5), -1.0, -math.sqrt(0.5),
    math.sqrt(0.5), 1.0, math.sqrt(1.5),
])


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_tsv(path):
    with Path(path).open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def write_tsv(path, rows, fields):
    with Path(path).open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t",
                                extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def cluster_t(beta, residual, design, inverse, labels):
    """CR1 cluster-robust t statistic for the ICL coefficient."""
    n, k = design.shape
    groups = np.unique(labels)
    g = len(groups)
    if n <= k or g < 2:
        return math.nan
    scores = np.vstack([
        design[labels == group].T @ residual[labels == group]
        for group in groups
    ])
    meat = scores.T @ scores
    correction = (g / (g - 1)) * ((n - 1) / (n - k))
    variance = correction * (inverse @ meat @ inverse)[-1, -1]
    if not math.isfinite(variance) or variance <= 0:
        return math.nan
    return float(beta[-1] / math.sqrt(variance))


def prepare_group(rows, block_by_row, min_blocks):
    n = len(rows)
    blocks = np.asarray([block_by_row[id(row)] for row in rows], dtype=str)
    unique_blocks = np.unique(blocks)
    output = {"n_clusters": n, "n_spatial_blocks": len(unique_blocks),
              "largest_spatial_block": max(np.bincount(
                  np.searchsorted(unique_blocks, blocks)), default=0)}
    if n < 5:
        return None, {**output, "status": "insufficient_clusters"}
    if len(unique_blocks) < min_blocks:
        return None, {**output, "status": "insufficient_spatial_blocks"}

    x = np.asarray([float(row["icl_delta_200_300"]) for row in rows])
    y = np.asarray([float(row["satellite_delta_50_75"]) for row in rows])
    z = np.asarray([float(row["z"]) for row in rows])
    richness = np.asarray([float(row["richness"]) for row in rows])
    nuisance = np.column_stack((np.ones(n), z, richness))
    design = np.column_stack((nuisance, x))
    if np.linalg.matrix_rank(design) != design.shape[1]:
        return None, {**output, "status": "rank_deficient"}

    bread = np.linalg.inv(design.T @ design)
    null_bread = np.linalg.inv(nuisance.T @ nuisance)
    null_beta = null_bread @ nuisance.T @ y
    null_fit = nuisance @ null_beta
    null_resid = y - null_fit
    # HC2 adjustment within the restricted model before applying block weights.
    leverage = np.sum((nuisance @ null_bread) * nuisance, axis=1)
    adjusted_resid = null_resid / np.sqrt(np.maximum(1.0 - leverage, 0.05))
    beta = bread @ design.T @ y
    residual = y - design @ beta
    t_observed = cluster_t(beta, residual, design, bread, blocks)
    if not math.isfinite(t_observed):
        return None, {**output, "status": "invalid_cluster_variance"}

    # Partial correlation is reported descriptively alongside the regression test.
    x_resid = x - nuisance @ (null_bread @ nuisance.T @ x)
    y_resid = null_resid
    denom = float(np.linalg.norm(x_resid) * np.linalg.norm(y_resid))
    partial_r = float(x_resid @ y_resid / denom) if denom > 0 else math.nan
    return {
        "rows": rows, "x": x, "y": y, "nuisance": nuisance,
        "design": design, "bread": bread, "blocks": blocks,
        "unique_blocks": unique_blocks, "null_fit": null_fit,
        "adjusted_resid": adjusted_resid, "t_observed": t_observed,
        "beta_icl": float(beta[-1]), "partial_r": partial_r,
        "n_clusters": n, "n_spatial_blocks": len(unique_blocks),
        "largest_spatial_block": output["largest_spatial_block"],
    }, {**output, "status": "eligible"}


def bootstrap_t(group, block_weights):
    blocks = group["blocks"]
    # All tests share the same global block weights, preserving cross-cut
    # dependence in the max statistic.
    weight = np.asarray([block_weights[b] for b in blocks])
    y_star = group["null_fit"] + group["adjusted_resid"] * weight
    design, bread = group["design"], group["bread"]
    beta = bread @ design.T @ y_star
    residual = y_star - design @ beta
    return cluster_t(beta, residual, design, bread, blocks)


def holm_adjust(p_values):
    order = np.argsort(p_values)
    adjusted = np.empty(len(p_values), dtype=float)
    running = 0.0
    m = len(p_values)
    for rank, index in enumerate(order):
        running = max(running, (m - rank) * p_values[index])
        adjusted[index] = min(1.0, running)
    return adjusted


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=INPUT)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--replicas", type=int, default=10000)
    parser.add_argument("--seed", type=int, default=20261018)
    parser.add_argument("--min-blocks", type=int, default=6)
    args = parser.parse_args()
    source = args.input.resolve()
    work = args.work.resolve()
    if args.replicas < 1000:
        parser.error("replicas must be at least 1000")
    if args.min_blocks < 2:
        parser.error("min-blocks must be at least 2")
    if work.exists():
        parser.error(f"Refusing to overwrite existing output: {work}")

    summary_path = source / "summary.json"
    detail_path = source / "matched_cluster_values.tsv"
    if not summary_path.is_file() or not detail_path.is_file():
        parser.error("Sensitivity summary or matched cluster data is missing")
    source_summary = json.loads(summary_path.read_text(encoding="utf-8"))
    if source_summary.get("schema") != "hsc-icl-satellite-cut-sensitivity-v1" \
            or source_summary.get("status") != "complete" \
            or source_summary.get("science_measurement") is not False:
        parser.error("Expected the verified exploratory cut-sensitivity output")
    if source_summary.get("output_sha256", {}).get(str(detail_path)) != sha256(detail_path):
        parser.error("Matched-cluster input checksum differs from its summary")

    rows = read_tsv(detail_path)
    if len(rows) != len({(row["cluster_id"], row["filter"], row["kron_factor"],
                          row["coverage_floor"]) for row in rows}):
        parser.error("Duplicate cluster/filter/configuration row")
    target = [{"ra_deg": float(row["target_ra_deg"]),
               "dec_deg": float(row["target_dec_deg"])} for row in rows]
    control = [{"ra_deg": float(row["control_ra_deg"]),
                "dec_deg": float(row["control_dec_deg"])} for row in rows]
    indices, labels, n_target, n_control, n_global_blocks, largest_global = pair_block_labels(
        target, control, np.ones(len(rows), dtype=bool))
    if len(indices) != len(rows):
        raise ValueError("Spatial block construction changed row order")
    block_by_row = {id(row): str(label) for row, label in zip(rows, labels)}

    grouped = {}
    filters = ("i", "i2")
    factors = (2.5, 3.5, 5.0)
    thresholds = (0.25, 0.35, 0.40, 0.50)
    selections = ("all_common_clusters", "center_eligible")
    for band in filters:
        for factor in factors:
            for threshold in thresholds:
                for selection in selections:
                    grouped[(band, factor, threshold, selection)] = []
    for row in rows:
        keybase = (row["filter"], float(row["kron_factor"]),
                   float(row["coverage_floor"]))
        grouped.setdefault((*keybase, "all_common_clusters"), []).append(row)
        if row["center_eligible"] == "True":
            grouped.setdefault((*keybase, "center_eligible"), []).append(row)

    prepared, records = {}, []
    for key in sorted(grouped):
        group, record = prepare_group(grouped[key], block_by_row, args.min_blocks)
        base = {"filter": key[0], "kron_factor": key[1],
                "coverage_floor": key[2], "selection": key[3]}
        if group is None:
            records.append({**base, **record})
        else:
            prepared[key] = group
    if len(prepared) < 2:
        parser.error("Fewer than two groups meet the prespecified support rule")

    global_blocks = sorted(set(str(label) for label in labels))
    rng = np.random.default_rng(args.seed)
    max_t_draws = np.empty(args.replicas, dtype=float)
    group_t_draws = {key: np.empty(args.replicas, dtype=float) for key in prepared}
    invalid_draws = {key: 0 for key in prepared}
    for replica in range(args.replicas):
        draw = rng.choice(WEBB_WEIGHTS, size=len(global_blocks), replace=True)
        weights = dict(zip(global_blocks, draw))
        maximum = 0.0
        for key, group in prepared.items():
            value = bootstrap_t(group, weights)
            group_t_draws[key][replica] = value
            if math.isfinite(value):
                maximum = max(maximum, abs(value))
            else:
                invalid_draws[key] += 1
        max_t_draws[replica] = maximum
    invalid = {key: value for key, value in invalid_draws.items() if value}
    if invalid:
        raise ValueError(f"Non-finite wild-bootstrap t draws: {invalid}")

    raw_p = np.asarray([
        (1 + np.count_nonzero(np.abs(group_t_draws[key]) >= abs(group["t_observed"])))
        / (args.replicas + 1)
        for key, group in prepared.items()
    ])
    max_t_p = np.asarray([
        (1 + np.count_nonzero(max_t_draws >= abs(group["t_observed"])))
        / (args.replicas + 1)
        for group in prepared.values()
    ])
    holm_p = holm_adjust(raw_p)
    for index, (key, group) in enumerate(prepared.items()):
        records.append({
            "filter": key[0], "kron_factor": key[1], "coverage_floor": key[2],
            "selection": key[3], "status": "exploratory_inference",
            "n_clusters": group["n_clusters"],
            "n_spatial_blocks": group["n_spatial_blocks"],
            "largest_spatial_block": group["largest_spatial_block"],
            "partial_pearson_r_z_richness": group["partial_r"],
            "beta_icl": group["beta_icl"], "t_cluster_cr1": group["t_observed"],
            "valid_bootstrap_replicas": args.replicas,
            "p_wild_cluster_raw": raw_p[index],
            "p_holm_across_eligible_groups": holm_p[index],
            "p_maxT_single_step_complete_null": max_t_p[index],
        })

    work.mkdir(parents=True, exist_ok=False)
    results_path = work / "familywise_association.tsv"
    fields = ["filter", "kron_factor", "coverage_floor", "selection", "status",
              "n_clusters", "n_spatial_blocks", "largest_spatial_block",
              "partial_pearson_r_z_richness", "beta_icl", "t_cluster_cr1",
              "valid_bootstrap_replicas",
              "p_wild_cluster_raw", "p_holm_across_eligible_groups",
              "p_maxT_single_step_complete_null"]
    write_tsv(results_path, records, fields)
    null_path = work / "maxT_null.tsv"
    write_tsv(null_path, [{"replica": i + 1, "max_abs_t": float(value)}
                          for i, value in enumerate(max_t_draws)],
              ["replica", "max_abs_t"])
    summary = {
        "schema": "hsc-icl-satellite-familywise-block-wild-v1",
        "status": "complete_exploratory", "science_measurement": False,
        "method": {
            "model": "satellite 50-75 kpc excess ~ ICL 200-300 kpc excess + redshift + richness",
            "null": "restricted model without ICL term; HC2-adjusted residuals",
            "bootstrap": "Webb six-point wild-cluster bootstrap-t; common weights by connected target/control 5-degree sky block",
            "multiplicity": "single-step max absolute t under complete null plus Holm adjusted raw bootstrap p-values",
            "eligible_rule": f"at least 5 clusters and {args.min_blocks} connected spatial blocks; exclusions retained in results",
            "family": "48 expected filter x Kron x coverage x selection groups; empty and support-limited groups retained",
            "replicas": args.replicas, "seed": args.seed,
        },
        "input": str(detail_path), "input_sha256": sha256(detail_path),
        "global_spatial_blocks": n_global_blocks,
        "global_target_cells": n_target, "global_control_cells": n_control,
        "largest_global_block_rows": largest_global,
        "n_expected_groups": len(filters) * len(factors) * len(thresholds) * len(selections),
        "n_eligible_groups": len(prepared), "n_reported_groups": len(records),
        "program_sha256": sha256(Path(__file__).resolve()),
        "output_sha256": {str(path): sha256(path) for path in (results_path, null_path)},
        "limitations": [
            "Max-T calibration is a single-step complete-null diagnostic; strong familywise control is not claimed because subset pivotality is unverified.",
            "The 5-degree connected blocks approximate spatial dependence and do not include pixel/coadd covariance or long-range sky systematics.",
            "The 24 mask/coverage configurations reuse the same images and are not independent replication.",
            "The satellite and ICL measurements remain selected, exploratory, and uncalibrated for PSF wings, sky transfer, extinction, and membership completeness.",
            "No effect is interpreted causally; science_measurement remains false.",
        ],
    }
    summary_path_out = work / "summary.json"
    summary_path_out.write_text(json.dumps(summary, indent=2, ensure_ascii=False,
                                           allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({"status": summary["status"], "science_measurement": False,
                      "global_spatial_blocks": n_global_blocks,
                      "eligible_groups": len(prepared), "reported_groups": len(records),
                      "results": str(results_path), "output": str(work)},
                     indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
