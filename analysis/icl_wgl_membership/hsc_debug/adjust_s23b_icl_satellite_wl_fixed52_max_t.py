#!/usr/bin/env python3
"""Familywise wild-cluster max-T audit for the fixed-52 joint ICL/WL family.

Uses only hash-pinned, already-derived Astrafek tables. No HSC/CAS request,
image access, raw shape-row CSV read, or data download occurs. All inference
remains exploratory and science_measurement stays false.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from evaluate_prestack_pair_support import pair_block_labels


DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
RELATION_DIR = DATA / "s23b_icl_satellite_wl_fixed52_relation_2026-11-11_03"
ICL_PROFILE_PATH = DATA / "s23b_icl_satellite_relation_2026-11-07_03/paired_icl_satellite_profiles.tsv"
DEFAULT_OUTPUT = DATA / "s23b_icl_satellite_wl_fixed52_maxt_2026-11-11_01"

RELATION_SUMMARY = RELATION_DIR / "summary.json"
METRICS_PATH = RELATION_DIR / "metrics.tsv"
JOINED_PATH = RELATION_DIR / "joined_lens_profiles.tsv"

EXPECTED_SHA256 = {
    "relation_summary.json": "a7a9f30107bb800986d5d1aa4bf4e2012497092170e3bc27bf246587f7840e2f",
    "metrics.tsv": "128bd628b15446d18edc23a932c055575264cdf237a747b5b7d57799ee6d0f69",
    "joined_lens_profiles.tsv": "2a9b5b22cba53dd98c2b81ac7266a3b18c0df3059c407031f54b0a639ba277e9",
    "icl_profiles.tsv": "98ed9e0022bfa68d5921ead62dce55be7df1cfa169e6a777d66d1e9760f43994",
}
EXPECTED_RELATION_RUNNER_SHA256 = "27839cf7bfa1683d94e6bf3dd30e2a7eac9a0f71d56c4bedeaa574beb3a41c8a"
EXPECTED_FIXED_IDS_SHA256 = "18312f27139b3750abdbf71543f22d9b7b18b2d3b24cbd29ed686a377524bb7d"
REPLICATES = 10_000
SEED = 20261111
MIN_LENSES = 12
MIN_BLOCKS = 6
WEBB = np.array([
    -math.sqrt(1.5), -1.0, -math.sqrt(0.5),
    math.sqrt(0.5), 1.0, math.sqrt(1.5),
])
SUPPORT_FLOORS = {"all": 0.0, "support_ge_0p10": 0.10}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_json(path: Path) -> dict:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def read_tsv(path: Path) -> list[dict[str, str]]:
    with Path(path).open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f"Refusing to write an empty table: {path}")
    with Path(path).open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(
            stream, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def id_sha256(values: set[str]) -> str:
    return hashlib.sha256("\n".join(sorted(values)).encode("utf-8")).hexdigest()


def partial_r(x: np.ndarray, y: np.ndarray, reduced: np.ndarray) -> float:
    residual_x = x - reduced @ np.linalg.lstsq(reduced, x, rcond=None)[0]
    residual_y = y - reduced @ np.linalg.lstsq(reduced, y, rcond=None)[0]
    denominator = float(np.linalg.norm(residual_x) * np.linalg.norm(residual_y))
    if denominator <= 0 or not math.isfinite(denominator):
        return float("nan")
    return float(residual_x @ residual_y / denominator)


def cr1_t(y: np.ndarray, model: dict) -> tuple[float, float]:
    design = model["design"]
    beta = model["xtx_inverse"] @ design.T @ y
    residual = y - design @ beta
    meat = np.zeros((design.shape[1], design.shape[1]), dtype=float)
    for indices in model["group_indices"]:
        score = design[indices].T @ residual[indices]
        meat += np.outer(score, score)
    correction = model["cr1_correction"]
    covariance = correction * model["xtx_inverse"] @ meat @ model["xtx_inverse"]
    variance = float(covariance[-1, -1])
    if not math.isfinite(variance) or variance <= 0:
        return float(beta[-1]), float("nan")
    return float(beta[-1]), float(beta[-1] / math.sqrt(variance))


def load_and_check() -> tuple[dict, list[dict], list[dict], dict[str, dict]]:
    sources = {
        "relation_summary.json": RELATION_SUMMARY,
        "metrics.tsv": METRICS_PATH,
        "joined_lens_profiles.tsv": JOINED_PATH,
        "icl_profiles.tsv": ICL_PROFILE_PATH,
    }
    observed = {name: sha256(path) for name, path in sources.items()}
    for name, actual in observed.items():
        if actual != EXPECTED_SHA256[name]:
            raise ValueError(f"Pinned input SHA-256 mismatch for {name}: {actual}")

    summary = read_json(RELATION_SUMMARY)
    if (summary.get("status") != "complete_exploratory_fixed52_joint_relation"
            or summary.get("science_measurement") is not False
            or summary.get("network_used") is not False
            or summary.get("hsc_query_submitted") is not False
            or summary.get("new_data_downloaded") is not False
            or summary.get("raw_source_shape_csvs_read") is not False
            or summary.get("images_read") is not False
            or summary.get("runner_sha256") != EXPECTED_RELATION_RUNNER_SHA256
            or summary.get("n_fixed_lenses") != 52
            or summary.get("i_and_i2_lens_sets_disjoint") is not True):
        raise ValueError("Relation summary has unexpected status, scope, runner, or sample")
    if summary.get("input_sha256", {}).get("icl_profiles.tsv") != EXPECTED_SHA256["icl_profiles.tsv"]:
        raise ValueError("Relation run does not reference the pinned ICL profile table")

    metrics = read_tsv(METRICS_PATH)
    joined = read_tsv(JOINED_PATH)
    icl_rows = read_tsv(ICL_PROFILE_PATH)
    if len(metrics) != 116 or len(joined) != 52 * 7:
        raise ValueError("Expected 116 association metrics and 52 by 7 joined lens profiles")
    fixed_ids = {row["cluster_id"] for row in joined}
    if len(fixed_ids) != 52 or id_sha256(fixed_ids) != EXPECTED_FIXED_IDS_SHA256:
        raise ValueError("Joined relation table differs from the fixed 52-lens set")

    by_cluster: dict[str, list[dict[str, str]]] = {}
    for row in joined:
        by_cluster.setdefault(row["cluster_id"], []).append(row)
    for cluster_id, rows in by_cluster.items():
        if len(rows) != 7 or {int(row["radial_bin"]) for row in rows} != set(range(7)):
            raise ValueError(f"Incomplete radial rows for {cluster_id}")
        for key in ("filter", "redshift", "richness", "icl_contrast_njy_arcsec2",
                    "satellite_log10_density_ratio", "inner_common_support_fraction",
                    "outer_common_support_fraction", "control_id"):
            if len({row[key] for row in rows}) != 1:
                raise ValueError(f"Inconsistent joined value {key} for {cluster_id}")

    icl_by_cluster = {}
    for row in icl_rows:
        cluster_id = str(row["target_id"])
        if cluster_id not in fixed_ids or row["filter"] not in {"i", "i2"}:
            continue
        if cluster_id in icl_by_cluster:
            raise ValueError(f"Duplicate ICL/filter row for {cluster_id}")
        icl_by_cluster[cluster_id] = row
    if set(icl_by_cluster) != fixed_ids:
        raise ValueError("ICL profile table does not provide exactly one row per fixed lens")

    for cluster_id, row in icl_by_cluster.items():
        joined_row = by_cluster[cluster_id][0]
        checks = (
            (row["filter"], joined_row["filter"]),
            (float(row["z_cl"]), float(joined_row["redshift"])),
            (float(row["richness"]), float(joined_row["richness"])),
            (float(row["icl_relative_contrast_njy_arcsec2"]),
             float(joined_row["icl_contrast_njy_arcsec2"])),
            (float(row["log10_inner_outer_candidate_density_ratio"]),
             float(joined_row["satellite_log10_density_ratio"])),
        )
        for left, right in checks:
            if isinstance(left, str):
                if left != right:
                    raise ValueError(f"Filter mismatch in joined ICL table for {cluster_id}")
            elif not math.isclose(left, right, rel_tol=1e-10, abs_tol=1e-10):
                raise ValueError(f"Profile mismatch in joined ICL table for {cluster_id}")

    return summary, metrics, joined, icl_by_cluster


def make_models(metrics: list[dict[str, str]], joined: list[dict[str, str]],
                icl_by_cluster: dict[str, dict], block_by_cluster: dict[str, str]) -> list[dict]:
    by_cluster: dict[str, list[dict[str, str]]] = {}
    for row in joined:
        by_cluster.setdefault(row["cluster_id"], []).append(row)
    models = []
    for metric in metrics:
        band = metric["filter"]
        support_group = metric["support_group"]
        floor = SUPPORT_FLOORS[support_group]
        selected_ids = [
            cluster_id for cluster_id in sorted(by_cluster)
            if icl_by_cluster[cluster_id]["filter"] == band
            and min(float(icl_by_cluster[cluster_id]["inner_common_support_fraction"]),
                    float(icl_by_cluster[cluster_id]["outer_common_support_fraction"])) >= floor
        ]

        relationship = metric["relationship"]
        rows_for_model = []
        if relationship == "icl_vs_satellite":
            for cluster_id in selected_ids:
                row = icl_by_cluster[cluster_id]
                rows_for_model.append((
                    cluster_id,
                    float(row["icl_relative_contrast_njy_arcsec2"]),
                    float(row["log10_inner_outer_candidate_density_ratio"]),
                    [float(row["z_cl"]), math.log(float(row["richness"]))],
                ))
        else:
            radial_bin = int(metric["radial_bin"])
            component = metric["component"]
            for cluster_id in selected_ids:
                row = next(r for r in by_cluster[cluster_id]
                           if int(r["radial_bin"]) == radial_bin)
                if int(row["n_lens_pairs"]) <= 0 or int(row["n_random_pairs"]) <= 0:
                    continue
                icl_value = float(row["icl_contrast_njy_arcsec2"])
                satellite_value = float(row["satellite_log10_density_ratio"])
                redshift = float(row["redshift"])
                log_richness = math.log(float(row["richness"]))
                if metric["predictor"] == "icl_contrast_njy_arcsec2":
                    x_value = icl_value
                    covariates = [satellite_value, redshift, log_richness]
                elif metric["predictor"] == "satellite_log10_density_ratio":
                    x_value = satellite_value
                    covariates = [icl_value, redshift, log_richness]
                else:
                    raise ValueError(f"Unknown WL predictor: {metric['predictor']}")
                outcome_key = ("delta_sigma_random_subtracted_msun_mpc2"
                               if component == "tangential"
                               else "cross_random_subtracted_msun_mpc2")
                rows_for_model.append((cluster_id, x_value, float(row[outcome_key]), covariates))

        n = len(rows_for_model)
        if n < MIN_LENSES:
            raise ValueError(f"Metric {metric} has fewer than {MIN_LENSES} lenses")
        x = np.asarray([row[1] for row in rows_for_model], dtype=float)
        y = np.asarray([row[2] for row in rows_for_model], dtype=float)
        covariates = np.asarray([row[3] for row in rows_for_model], dtype=float)
        if not (np.isfinite(x).all() and np.isfinite(y).all()
                and np.isfinite(covariates).all()):
            raise ValueError("Non-finite value in a fixed-sample relation model")
        reduced = np.column_stack((np.ones(n), covariates))
        design = np.column_stack((reduced, x))
        if np.linalg.matrix_rank(design) != design.shape[1]:
            raise ValueError(f"Rank-deficient relation design for {metric}")
        labels = np.asarray([block_by_cluster[row[0]] for row in rows_for_model])
        blocks = np.unique(labels)
        if len(blocks) < MIN_BLOCKS:
            raise ValueError(f"Metric has fewer than {MIN_BLOCKS} spatial blocks: {metric}")
        observed_r = partial_r(x, y, reduced)
        xtx_inverse = np.linalg.inv(design.T @ design)
        model = {
            "filter": band,
            "support_group": support_group,
            "relationship": relationship,
            "component": metric["component"],
            "radial_bin": metric["radial_bin"],
            "predictor": metric["predictor"],
            "outcome": metric["outcome"],
            "covariates": metric["covariates"],
            "cluster_ids": [row[0] for row in rows_for_model],
            "labels": labels,
            "design": design,
            "xtx_inverse": xtx_inverse,
            "group_indices": [np.flatnonzero(labels == block) for block in blocks],
            "cr1_correction": (len(blocks) / (len(blocks) - 1.0))
                              * ((n - 1.0) / (n - design.shape[1])),
            "reduced_fit": reduced @ np.linalg.lstsq(reduced, y, rcond=None)[0],
            "restricted_resid": None,
            "n_lenses": n,
            "n_spatial_blocks": len(blocks),
            "relation_reported_n_spatial_blocks": int(metric["n_spatial_blocks"]),
            "largest_block_n": max(int(np.count_nonzero(labels == block)) for block in blocks),
            "partial_r": observed_r,
        }
        model["restricted_resid"] = y - model["reduced_fit"]
        model["observed_slope"], model["observed_t"] = cr1_t(y, model)
        if not math.isfinite(model["partial_r"]) or not math.isfinite(model["observed_t"]):
            raise ValueError(f"Non-finite observed test statistic for {metric}")
        if (model["n_lenses"] != int(metric["n_lenses"])
                or not math.isclose(model["partial_r"], float(metric["partial_r"]),
                                    rel_tol=1e-9, abs_tol=1e-10)):
            raise ValueError(f"Recreated point estimate differs from pinned metric: {metric}")
        models.append(model)
    return models


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--replicates", type=int, default=REPLICATES)
    args = parser.parse_args()
    if args.output.exists():
        parser.error(f"Refusing to overwrite an existing result: {args.output}")
    if args.replicates < 1000:
        parser.error("At least 1,000 Webb wild-cluster replicates are required")

    summary, metrics, joined, icl_by_cluster = load_and_check()
    ids = sorted(icl_by_cluster)
    target_rows = [{
        "ra_deg": float(icl_by_cluster[cluster_id]["ra_deg"]),
        "dec_deg": float(icl_by_cluster[cluster_id]["dec_deg"]),
    } for cluster_id in ids]
    control_rows = [{
        "ra_deg": float(icl_by_cluster[cluster_id]["control_ra_deg"]),
        "dec_deg": float(icl_by_cluster[cluster_id]["control_dec_deg"]),
    } for cluster_id in ids]
    _, full_labels, n_target_cells, n_control_cells, n_blocks, largest_block = pair_block_labels(
        target_rows, control_rows, np.ones(len(ids), dtype=bool)
    )
    block_by_cluster = dict(zip(ids, map(str, full_labels)))
    if len(block_by_cluster) != 52 or n_blocks < MIN_BLOCKS:
        raise ValueError("The common fixed-sample spatial block partition is invalid")

    models = make_models(metrics, joined, icl_by_cluster, block_by_cluster)
    family_size = len(models)
    if family_size != 116:
        raise ValueError(f"Expected the complete 116-test relation family, found {family_size}")
    direct_n = sum(model["relationship"] == "icl_vs_satellite" for model in models)
    wl_n = family_size - direct_n
    if direct_n != 4 or wl_n != 112:
        raise ValueError("The direct ICL/satellite and WL family counts changed")

    block_ids = sorted(set(block_by_cluster.values()))
    rng = np.random.default_rng(SEED)
    t_draws = np.full((args.replicates, family_size), np.nan, dtype=float)
    for draw in range(args.replicates):
        multipliers = {
            block: float(WEBB[rng.integers(0, len(WEBB))]) for block in block_ids
        }
        for index, model in enumerate(models):
            wild = np.asarray([multipliers[label] for label in model["labels"]])
            y_star = model["reduced_fit"] + model["restricted_resid"] * wild
            _, statistic = cr1_t(y_star, model)
            t_draws[draw, index] = statistic
        if (draw + 1) % 1000 == 0:
            print(f"bootstrap {draw + 1}/{args.replicates}", flush=True)

    valid = np.isfinite(t_draws).all(axis=1)
    valid_draws = t_draws[valid]
    if len(valid_draws) < int(0.99 * args.replicates):
        raise RuntimeError(f"Only {len(valid_draws)}/{args.replicates} complete max-T draws")
    max_abs = np.max(np.abs(valid_draws), axis=1)

    output_rows = []
    for index, model in enumerate(models):
        abs_draw = np.abs(valid_draws[:, index])
        observed = abs(model["observed_t"])
        p_unadjusted = (1 + int(np.count_nonzero(abs_draw >= observed))) / (1 + len(abs_draw))
        p_adjusted = (1 + int(np.count_nonzero(max_abs >= observed))) / (1 + len(max_abs))
        output_rows.append({
            "filter": model["filter"],
            "support_group": model["support_group"],
            "relationship": model["relationship"],
            "predictor": model["predictor"],
            "outcome": model["outcome"],
            "component": model["component"],
            "radial_bin": model["radial_bin"],
            "n_lenses": model["n_lenses"],
            "n_spatial_blocks": model["n_spatial_blocks"],
            "relation_reported_n_spatial_blocks": model["relation_reported_n_spatial_blocks"],
            "largest_block_n": model["largest_block_n"],
            "partial_r": model["partial_r"],
            "ols_slope": model["observed_slope"],
            "cr1_t": model["observed_t"],
            "wild_cluster_p_unadjusted": p_unadjusted,
            "wild_cluster_p_max_t_116": p_adjusted,
            "bootstrap_replicates": len(valid_draws),
            "inferential_claim": False,
        })

    raw_rejections = sum(row["wild_cluster_p_unadjusted"] < 0.05 for row in output_rows)
    adjusted_rejections = sum(row["wild_cluster_p_max_t_116"] < 0.05 for row in output_rows)
    smallest_adjusted = min(output_rows, key=lambda row: row["wild_cluster_p_max_t_116"])
    summary_path = args.output.resolve()
    runner_path = Path(__file__).resolve()
    input_hashes = {
        "relation_summary.json": sha256(RELATION_SUMMARY),
        "metrics.tsv": sha256(METRICS_PATH),
        "joined_lens_profiles.tsv": sha256(JOINED_PATH),
        "icl_profiles.tsv": sha256(ICL_PROFILE_PATH),
    }

    args.output.mkdir(parents=True, exist_ok=False)
    table_path = args.output / "max_t_metrics.tsv"
    write_tsv(table_path, output_rows)
    pointwise_excludes = summary["key_results"]["wl_pointwise_intervals_excluding_zero"]
    pointwise_total = summary["key_results"]["wl_pointwise_intervals_total"]
    report_lines = [
        "# S23B fixed-52 ICL–satellite–weak-lensing familywise audit",
        "",
        "This analysis applies a restricted Webb six-point wild-cluster bootstrap with a shared 5-degree target/control block multiplier to the 116 partial-correlation tests in the saved fixed-52 relation. Four tests are direct ICL–satellite relations; 112 test the WL profile against ICL or satellite structure. The analysis reads only saved tables and makes no HSC/CAS request, image read, raw source-row CSV read, or download. `science_measurement=false`.",
        "",
        f"The fixed 52-lens sample forms {n_blocks} connected 5-degree blocks (largest block: {largest_block} lenses). The full family uses {args.replicates:,} shared-multiplier replicates; {len(valid_draws):,} were valid. CR1 cluster-studentized OLS slopes are tested under the corresponding restricted zero-slope model. This is a complementary null-based familywise procedure; it does not adjust the percentile intervals in the relation report.",
        "",
        f"The earlier percentile-bootstrap report has {pointwise_excludes} of {pointwise_total} WL partial-correlation intervals that exclude zero pointwise. The wild-cluster max-T family has {raw_rejections} unadjusted p-values below 0.05 and {adjusted_rejections} max-T adjusted p-values below 0.05. Its smallest adjusted p-value is {smallest_adjusted['wild_cluster_p_max_t_116']:.4f} for {smallest_adjusted['filter']} / {smallest_adjusted['support_group']} / {smallest_adjusted['relationship']} / {smallest_adjusted['component']} / radial bin {smallest_adjusted['radial_bin']} (partial r = {smallest_adjusted['partial_r']:.3f}).",
        "",
        "Interpretation remains limited to this selected 52-lens sample. The 5-degree connected-cell partition is an approximate spatial dependence model, not a validated HSC survey covariance. The analysis does not resolve source photo-z, shear response, additive/multiplicative residuals, PSF, magnification, member dilution, masks, or source reuse. It does not establish a physical ICL–satellite–halo relation or a detection.",
        "",
        "## Provenance",
        "",
        f"- Runner SHA-256: `{sha256(runner_path)}`",
        f"- Relation input runner SHA-256: `{summary['runner_sha256']}`",
        f"- Fixed lens ID SHA-256: `{EXPECTED_FIXED_IDS_SHA256}`",
        f"- Family size: {family_size} (4 direct ICL–satellite, 112 WL association tests)",
        f"- Result directory: `{summary_path}`",
        "- Full per-test estimates and p-values: `max_t_metrics.tsv`",
    ]
    assessment_path = args.output / "assessment.md"
    assessment_path.write_text("\n".join(report_lines) + "\n", encoding="utf-8")
    result = {
        "schema": "astrafek-s23b-icl-satellite-wl-fixed52-max-t-v1",
        "status": "complete_exploratory_familywise_max_t",
        "science_measurement": False,
        "network_used": False,
        "hsc_query_submitted": False,
        "new_data_downloaded": False,
        "images_read": False,
        "raw_source_shape_csvs_read": False,
        "source_tables_read": [str(RELATION_SUMMARY), str(METRICS_PATH),
                                str(JOINED_PATH), str(ICL_PROFILE_PATH)],
        "family_size": family_size,
        "n_direct_icl_satellite_tests": direct_n,
        "n_wl_association_tests": wl_n,
        "n_fixed_lenses": len(ids),
        "fixed_lens_ids_sha256": EXPECTED_FIXED_IDS_SHA256,
        "spatial_partition": {
            "definition": "connected components of target/control 5-degree spatial cells, built once on all 52 paired centers and shared across the family",
            "n_blocks": n_blocks,
            "n_target_cells": n_target_cells,
            "n_control_cells": n_control_cells,
            "largest_block_n": largest_block,
        },
        "method": "restricted Webb six-point wild-cluster bootstrap under each zero-predictor-slope null; CR1 cluster-studentized OLS; one shared multiplier per 5-degree block per replicate; max absolute t across the full 116-test family",
        "replicates_requested": args.replicates,
        "replicates_valid": int(len(valid_draws)),
        "seed": SEED,
        "key_results": {
            "wild_cluster_p_unadjusted_below_0p05": raw_rejections,
            "wild_cluster_p_max_t_116_below_0p05": adjusted_rejections,
            "minimum_adjusted_p": smallest_adjusted["wild_cluster_p_max_t_116"],
            "minimum_adjusted_metric": {
                key: smallest_adjusted[key] for key in (
                    "filter", "support_group", "relationship", "component",
                    "radial_bin", "partial_r", "wild_cluster_p_unadjusted",
                    "wild_cluster_p_max_t_116",
                )
            },
            "prior_wl_percentile_intervals_excluding_zero": pointwise_excludes,
            "prior_wl_percentile_intervals_total": pointwise_total,
        },
        "inputs": {
            "sha256": input_hashes,
            "relation_runner_sha256": summary["runner_sha256"],
        },
        "runner_path": str(runner_path),
        "runner_sha256": sha256(runner_path),
        "output_directory": str(summary_path),
        "output_sha256": {
            "max_t_metrics.tsv": sha256(table_path),
            "assessment.md": sha256(assessment_path),
        },
        "limitations": [
            "selected fixed 52-lens support sample is not representative of all S23B clusters",
            "connected 5-degree cells are approximate spatial blocks, not validated HSC survey regions",
            "small and uneven block counts limit wild-cluster inference",
            "max-T p-values complement rather than adjust the saved percentile bootstrap intervals",
            "inherits unresolved photo-z, shear calibration, PSF, magnification, mask, member-dilution, source-reuse, and image-support limitations",
        ],
    }
    summary_path_file = args.output / "summary.json"
    summary_path_file.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n",
                                 encoding="utf-8")
    print(json.dumps({
        "status": result["status"],
        "science_measurement": False,
        "family_size": family_size,
        "blocks": n_blocks,
        "valid_replicates": int(len(valid_draws)),
        "raw_rejections": raw_rejections,
        "max_t_rejections": adjusted_rejections,
        "minimum_adjusted_metric": result["key_results"]["minimum_adjusted_metric"],
        "output_directory": str(summary_path),
    }, ensure_ascii=False, indent=2), flush=True)


if __name__ == "__main__":
    main()
