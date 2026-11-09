#!/usr/bin/env python3
"""Assess covariate-based selection sensitivity in the 224-to-94 ICL sample.

This offline diagnostic fits a pooled logistic support model using only
pre-outcome covariates measured for all 224 ICL-supported cluster/filter pairs.
It compares raw and inverse-probability weighted partial correlations for the
94 pairs with a saved 50--75 kpc satellite profile. Weights are exploratory;
the missing-at-random assumption is not established and no science gate is
changed.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
from scipy.special import expit

HERE = Path(__file__).resolve().parent
import sys
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from evaluate_prestack_pair_support import pair_block_labels
from audit_icl_target_control_outer_sky import partial_pearson

DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
ICL_METRICS = DATA / "icl_annulus_support_2026-10-18_01/pair_metrics_k2p5_c40_r50_300.tsv"
SUPPORTED = DATA / "icl_satellite_supported_94_association_2026-11-05_01/supported_pairs.tsv"
TARGET_MANIFEST = DATA / "sample_manifest.json"
CONTROL_MANIFEST = DATA / "control_download_01/download_manifest.json"
DEFAULT_OUTPUT = DATA / "icl_satellite_joint_support_ipw_2026-11-05_06"
BOOTSTRAPS = 2000
SEED = 20261105


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def fit_logit(x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray, dict]:
    if x.ndim != 2 or y.shape != (x.shape[0],) or len(np.unique(y)) != 2:
        raise ValueError("invalid or non-varying logistic selection input")

    def objective(beta: np.ndarray) -> tuple[float, np.ndarray]:
        eta = x @ beta
        probability = expit(eta)
        return (float(np.logaddexp(0.0, eta).sum() - np.dot(y, eta)),
                x.T @ (probability - y))

    result = minimize(lambda beta: objective(beta), np.zeros(x.shape[1]), jac=True,
                      method="BFGS", options={"gtol": 1e-7, "maxiter": 2000})
    gradient_inf_norm = float(np.linalg.norm(result.jac, ord=np.inf))
    if not np.isfinite(gradient_inf_norm) or gradient_inf_norm > 2e-5:
        raise RuntimeError(
            f"selection logistic fit failed stationarity check: {result.message}; "
            f"gradient_inf_norm={gradient_inf_norm:.3g}"
        )
    p = np.clip(expit(x @ result.x), 1e-6, 1.0 - 1e-6)
    return result.x, p, {
        "converged_by_stationarity": True,
        "optimizer_reported_success": bool(result.success),
        "optimizer_message": str(result.message),
        "gradient_inf_norm": gradient_inf_norm,
        "iterations": int(result.nit),
    }


def weighted_partial_r(x: np.ndarray, y: np.ndarray, z: np.ndarray,
                       richness: np.ndarray, weights: np.ndarray) -> float:
    if len(x) < 6 or not (len(x) == len(y) == len(z) == len(richness) == len(weights)):
        return float("nan")
    if not all(np.isfinite(v).all() for v in (x, y, z, richness, weights)):
        return float("nan")
    if np.any(weights <= 0) or np.linalg.matrix_rank(
            np.column_stack((np.ones(len(x)), z, richness))) < 3:
        return float("nan")
    design = np.column_stack((np.ones(len(x)), z, richness))
    root_w = np.sqrt(weights)
    beta_x = np.linalg.lstsq(design * root_w[:, None], x * root_w, rcond=None)[0]
    beta_y = np.linalg.lstsq(design * root_w[:, None], y * root_w, rcond=None)[0]
    rx, ry = x - design @ beta_x, y - design @ beta_y
    mx = np.average(rx, weights=weights)
    my = np.average(ry, weights=weights)
    rx, ry = rx - mx, ry - my
    den = np.sqrt(np.dot(weights, rx * rx) * np.dot(weights, ry * ry))
    return float(np.dot(weights, rx * ry) / den) if den > 0 else float("nan")


def association_under_propensity(rows: list[dict], selected: np.ndarray,
                                 probability: np.ndarray) -> dict:
    result = {}
    for band in ("i", "i2"):
        idx_all = np.asarray([i for i, row in enumerate(rows) if row["filter"] == band])
        idx_sel = idx_all[selected[idx_all]]
        prevalence = float(selected[idx_all].mean())
        weights = prevalence / probability[idx_sel]
        weights /= weights.mean()
        data = [rows[i] for i in idx_sel]
        eff_n = float(weights.sum() ** 2 / np.dot(weights, weights))
        result[band] = {
            "n_supported": len(idx_sel),
            "ipw_partial_r_z_richness": weighted_partial_r(
                np.asarray([r["icl_delta_200_300"] for r in data]),
                np.asarray([r["satellite_delta_50_75"] for r in data]),
                np.asarray([r["redshift"] for r in data]),
                np.asarray([r["richness"] for r in data]), weights),
            "effective_n": eff_n,
            "max_stabilized_weight": float(weights.max()),
            "selected_propensity_min": float(probability[idx_sel].min()),
            "selected_propensity_max": float(probability[idx_sel].max()),
        }
    return result


def run(output: Path, replicas: int, seed: int) -> dict:
    inputs = (ICL_METRICS, SUPPORTED, TARGET_MANIFEST, CONTROL_MANIFEST)
    for path in inputs:
        if not path.is_file():
            raise FileNotFoundError(path)
    if output.exists():
        raise FileExistsError(f"refusing to overwrite {output}")
    metrics = read_tsv(ICL_METRICS)
    selected_rows = read_tsv(SUPPORTED)
    target_doc = json.loads(TARGET_MANIFEST.read_text(encoding="utf-8"))
    control_doc = json.loads(CONTROL_MANIFEST.read_text(encoding="utf-8"))
    targets = {str(row["id"]): row for row in target_doc["clusters"]}
    controls = {str(row["target_cluster_id"]): row for row in control_doc["clusters"]}
    selected_by_key = {(row["cluster_id"], row["filter"]): row for row in selected_rows}
    if len(metrics) != 224 or len({(r["id"], r["filter"]) for r in metrics}) != 224:
        raise ValueError("expected 224 unique ICL-supported cluster/filter pairs")
    if len(selected_rows) != 94 or len(selected_by_key) != 94:
        raise ValueError("expected 94 unique supported 50--75 kpc satellite profiles")
    keys = {(row["id"], row["filter"]) for row in metrics}
    if not set(selected_by_key).issubset(keys):
        raise ValueError("satellite support sample is not nested in the ICL-supported set")

    rows = []
    for source in sorted(metrics, key=lambda r: (int(r["id"]), r["filter"])):
        ident, band = source["id"], source["filter"]
        if ident not in targets or ident not in controls:
            raise ValueError(f"missing target/control coordinates for {ident}")
        sat = selected_by_key.get((ident, band))
        z = float(source["z"])
        richness = float(source["richness"])
        coverage = float(source["coverage_200_300"])
        masked = source["center_masked_k2p5"] == "True"
        row = {
            "cluster_id": ident, "filter": band, "selected_50_75": sat is not None,
            "icl_delta_200_300": float(sat["icl_delta_200_300"]) if sat else float("nan"),
            "satellite_delta_50_75": float(sat["satellite_delta_50_75"]) if sat else float("nan"),
            "redshift": z, "richness": richness, "icl_coverage_200_300": coverage,
            "center_masked_k2p5": masked,
            "target_ra_deg": float(targets[ident]["ra_deg"]),
            "target_dec_deg": float(targets[ident]["dec_deg"]),
            "control_ra_deg": float(controls[ident]["ra_deg"]),
            "control_dec_deg": float(controls[ident]["dec_deg"]),
        }
        if sat and (abs(row["icl_delta_200_300"] - float(source["delta_sb_200_300"])) > 5e-7
                    or abs(z - float(sat["redshift"])) > 1e-8
                    or abs(richness - float(sat["richness"])) > 1e-6):
            raise ValueError(f"outcome/covariate mismatch in saved pair {ident} {band}")
        rows.append(row)

    y = np.asarray([row["selected_50_75"] for row in rows], dtype=float)
    z = np.asarray([row["redshift"] for row in rows])
    richness = np.asarray([row["richness"] for row in rows])
    log_richness = np.log(richness)
    coverage = np.asarray([row["icl_coverage_200_300"] for row in rows])
    masked = np.asarray([row["center_masked_k2p5"] for row in rows], dtype=float)
    is_i2 = np.asarray([row["filter"] == "i2" for row in rows], dtype=float)
    standardize = lambda a: (a - a.mean()) / a.std(ddof=0)
    z_s, rich_s, cov_s = standardize(z), standardize(log_richness), standardize(coverage)
    # Prespecified descriptive model: filter, redshift, richness, ICL radial
    # coverage and center-mask state, all measured for selected and unselected pairs.
    design = np.column_stack((np.ones(len(rows)), is_i2, z_s, rich_s, cov_s, masked))
    names = ["intercept", "filter_i2", "redshift_per_SD", "log_richness_per_SD",
             "icl_coverage_per_SD", "center_masked"]
    beta, propensity, fit_meta = fit_logit(design, y)
    for row, p in zip(rows, propensity):
        row["estimated_support_probability"] = float(p)
        row["unstabilized_ipw_if_selected"] = float(1.0 / p)

    # Model specification sensitivity: remove center-mask state, allow
    # nonlinear coverage, or fit filter-specific slopes. These are diagnostic
    # alternatives, not model-selection steps based on the association result.
    no_mask_design = design[:, :5]
    _, p_no_mask, meta_no_mask = fit_logit(no_mask_design, y)
    quadratic_design = np.column_stack((design, cov_s ** 2))
    _, p_quadratic, meta_quadratic = fit_logit(quadratic_design, y)
    p_filter_specific = np.full(len(y), np.nan)
    filter_specific_meta = {}
    for band in ("i", "i2"):
        take = np.asarray([row["filter"] == band for row in rows])
        band_design = np.column_stack((np.ones(take.sum()), z_s[take], rich_s[take],
                                       cov_s[take], masked[take]))
        _, band_probability, band_meta = fit_logit(band_design, y[take])
        p_filter_specific[take] = band_probability
        filter_specific_meta[band] = band_meta
    model_specs = {
        "primary_pooled_additive": {"probability": propensity, "fit": fit_meta},
        "pooled_without_center_mask": {"probability": p_no_mask, "fit": meta_no_mask},
        "pooled_quadratic_coverage": {"probability": p_quadratic, "fit": meta_quadratic},
        "filter_specific_additive": {"probability": p_filter_specific,
                                      "fit": filter_specific_meta},
    }
    model_sensitivity = {
        name: {"fit": item["fit"], "association_by_filter": association_under_propensity(
            rows, y.astype(bool), item["probability"])}
        for name, item in model_specs.items()
    }

    by_filter = {}
    selected_indices = np.flatnonzero(y == 1)
    for band in ("i", "i2"):
        idx_all = np.asarray([i for i, row in enumerate(rows) if row["filter"] == band])
        idx_sel = idx_all[y[idx_all] == 1]
        idx_not_sel = idx_all[y[idx_all] == 0]
        balance = {}
        for name, values in (("redshift", z), ("richness", richness),
                             ("icl_coverage_200_300", coverage),
                             ("center_masked_k2p5", masked)):
            selected_values, unselected_values = values[idx_sel], values[idx_not_sel]
            pooled_sd = np.sqrt((selected_values.var(ddof=1) + unselected_values.var(ddof=1)) / 2.0)
            balance[name] = {
                "selected_mean": float(selected_values.mean()),
                "unselected_mean": float(unselected_values.mean()),
                "standardized_mean_difference": float(
                    (selected_values.mean() - unselected_values.mean()) / pooled_sd
                    if pooled_sd > 0 else 0.0),
            }
        raw_rows = [selected_by_key[(rows[i]["cluster_id"], band)] for i in idx_sel]
        raw_r = partial_pearson(
            [float(r["icl_delta_200_300"]) for r in raw_rows],
            [float(r["satellite_delta_50_75"]) for r in raw_rows],
            [[float(r["redshift"]) for r in raw_rows],
             [float(r["richness"]) for r in raw_rows]],
        )
        prevalence = float(y[idx_all].mean())
        w = prevalence / propensity[idx_sel]
        w = w / w.mean()
        selected_data = [rows[i] for i in idx_sel]
        weighted_r = weighted_partial_r(
            np.asarray([r["icl_delta_200_300"] for r in selected_data]),
            np.asarray([r["satellite_delta_50_75"] for r in selected_data]),
            np.asarray([r["redshift"] for r in selected_data]),
            np.asarray([r["richness"] for r in selected_data]), w,
        )
        raw_r_logrichness = partial_pearson(
            [float(r["icl_delta_200_300"]) for r in raw_rows],
            [float(r["satellite_delta_50_75"]) for r in raw_rows],
            [[float(r["redshift"]) for r in raw_rows],
             [float(np.log(float(r["richness"]))) for r in raw_rows]],
        )
        weighted_r_logrichness = weighted_partial_r(
            np.asarray([r["icl_delta_200_300"] for r in selected_data]),
            np.asarray([r["satellite_delta_50_75"] for r in selected_data]),
            np.asarray([r["redshift"] for r in selected_data]),
            np.log(np.asarray([r["richness"] for r in selected_data])), w,
        )
        effective_n = float(w.sum() ** 2 / np.dot(w, w))
        by_filter[band] = {
            "n_icl_supported": int(len(idx_all)), "n_satellite_supported": int(len(idx_sel)),
            "support_fraction": prevalence, "selection_balance": balance,
            "raw_partial_r_z_richness": float(raw_r),
            "ipw_partial_r_z_richness": weighted_r,
            "ipw_minus_raw": weighted_r - float(raw_r),
            "raw_partial_r_z_logrichness": float(raw_r_logrichness),
            "ipw_partial_r_z_logrichness": weighted_r_logrichness,
            "ipw_minus_raw_z_logrichness": weighted_r_logrichness - float(raw_r_logrichness),
            "stabilized_weight_min": float(w.min()), "stabilized_weight_median": float(np.median(w)),
            "stabilized_weight_p95": float(np.percentile(w, 95)),
            "stabilized_weight_max": float(w.max()), "effective_n": effective_n,
            "effective_n_fraction_of_selected": effective_n / len(idx_sel),
        }

    target_coords = [{"ra_deg": r["target_ra_deg"], "dec_deg": r["target_dec_deg"]} for r in rows]
    control_coords = [{"ra_deg": r["control_ra_deg"], "dec_deg": r["control_dec_deg"]} for r in rows]
    _, block_labels, n_target_cells, n_control_cells, n_blocks, largest = pair_block_labels(
        target_coords, control_coords, np.ones(len(rows), dtype=bool))
    groups = np.unique(block_labels)
    by_block = {group: np.flatnonzero(block_labels == group) for group in groups}
    rng = np.random.default_rng(seed)
    draws = {band: [] for band in ("i", "i2")}
    fit_failures = 0
    precision_warnings = 0
    max_warning_gradient = 0.0
    max_warning_coefficient = 0.0
    too_few = 0
    for _ in range(replicas):
        sampled_groups = rng.choice(groups, size=len(groups), replace=True)
        idx = np.concatenate([by_block[g] for g in sampled_groups])
        try:
            b_draw, p_draw, draw_fit = fit_logit(design[idx], y[idx])
        except (ValueError, RuntimeError, np.linalg.LinAlgError):
            fit_failures += 1
            continue
        if not draw_fit["optimizer_reported_success"]:
            precision_warnings += 1
            max_warning_gradient = max(max_warning_gradient, draw_fit["gradient_inf_norm"])
            max_warning_coefficient = max(max_warning_coefficient,
                                          float(np.max(np.abs(b_draw))))
        for band in ("i", "i2"):
            use = np.asarray([rows[i]["filter"] == band and y[i] == 1 for i in idx])
            selected_i = idx[use]
            if len(selected_i) < 6:
                too_few += 1
                continue
            d = [rows[i] for i in selected_i]
            probs = p_draw[use]
            weights = np.asarray([y[idx][use].mean()]) / probs
            weights /= weights.mean()
            corr = weighted_partial_r(
                np.asarray([r["icl_delta_200_300"] for r in d]),
                np.asarray([r["satellite_delta_50_75"] for r in d]),
                np.asarray([r["redshift"] for r in d]),
                np.asarray([r["richness"] for r in d]), weights,
            )
            if np.isfinite(corr):
                draws[band].append(corr)
    bootstrap = {}
    for band, values in draws.items():
        if len(values) < replicas * 0.8:
            raise RuntimeError(f"too few valid {band} block-bootstrap estimates: {len(values)}/{replicas}")
        bootstrap[band] = {
            "ipw_partial_r_ci95": np.percentile(values, [2.5, 97.5]).tolist(),
            "fraction_positive": float(np.mean(np.asarray(values) > 0)),
            "valid_replicates": len(values),
        }

    output.mkdir(parents=True, exist_ok=False)
    pair_path = output / "pair_selection_propensities.tsv"
    with pair_path.open("x", encoding="utf-8", newline="") as stream:
        columns = [*rows[0].keys(), "connected_5deg_block"]
        writer = csv.DictWriter(stream, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row, label in zip(rows, block_labels):
            writer.writerow({**row, "connected_5deg_block": label})
    summary = {
        "schema": "astrafek-icl-satellite-joint-support-ipw-sensitivity-v1",
        "status": "complete_offline_sensitivity",
        "science_measurement": False, "hsc_network_used": False,
        "parent_icl_supported_pairs": len(rows), "selected_satellite_profile_pairs": int(y.sum()),
        "selection_model": {
            "model": "pooled unpenalized logistic regression",
            "predictors": names, "coefficients": {n: float(v) for n, v in zip(names, beta)},
            "fit": fit_meta, "propensity_min": float(propensity.min()),
            "propensity_p01": float(np.percentile(propensity, 1)),
            "propensity_median": float(np.median(propensity)),
            "propensity_p99": float(np.percentile(propensity, 99)),
            "propensity_max": float(propensity.max()),
            "selected_propensity_min": float(propensity[selected_indices].min()),
            "selected_propensity_max": float(propensity[selected_indices].max()),
            "missing_not_at_random_assumption_verified": False,
        },
        "association_by_filter": by_filter,
        "model_specification_sensitivity": model_sensitivity,
        "block_bootstrap": {
            "method": "resample connected 5-degree target/control sky blocks; refit pooled selection model in each draw",
            "n_connected_blocks": int(n_blocks), "largest_block_pairs": int(largest),
            "target_cells": int(n_target_cells), "control_cells": int(n_control_cells),
            "replicates_requested": replicas,
            "selection_model_nonconvergent_draws": fit_failures,
            "optimizer_precision_warnings_accepted_by_gradient": precision_warnings,
            "max_accepted_warning_gradient_inf_norm": max_warning_gradient,
            "max_abs_coefficient_in_precision_warning_draws": max_warning_coefficient,
            "too_few_selected_rows": too_few, "by_filter": bootstrap,
            "seed": seed,
        },
        "interpretation": "IPW is a model-dependent sensitivity only. It does not correct selection unless profile availability is conditionally independent of unobserved satellite signal given the chosen covariates; that assumption is not established. No detection or completeness claim is made.",
        "input_sha256": {str(path): sha256(path) for path in inputs},
        "output_sha256": {str(pair_path): sha256(pair_path)},
    }
    summary_path = output / "summary.json"
    summary_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
                            encoding="utf-8")
    summary["summary_sha256"] = sha256(summary_path)
    print(json.dumps({
        "status": summary["status"], "science_measurement": False,
        "selected": summary["selected_satellite_profile_pairs"],
        "blocks": n_blocks, "support_model": summary["selection_model"],
        "association_by_filter": by_filter, "bootstrap": bootstrap,
        "output": str(output.resolve()),
    }, ensure_ascii=False, indent=2))
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--replicas", type=int, default=BOOTSTRAPS)
    parser.add_argument("--seed", type=int, default=SEED)
    args = parser.parse_args()
    if args.replicas < 500:
        parser.error("at least 500 connected-block bootstrap replicates are required")
    run(args.output.resolve(), args.replicas, args.seed)


if __name__ == "__main__":
    main()
