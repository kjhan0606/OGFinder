#!/usr/bin/env python3
"""Test whether delivered variance metrics alter satellite-support IPW results."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
from scipy.optimize import minimize
from scipy.special import expit

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from evaluate_prestack_pair_support import pair_block_labels
from audit_icl_target_control_outer_sky import partial_pearson

DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
ICL_METRICS = DATA / "icl_annulus_support_2026-10-18_01/pair_metrics_k2p5_c40_r50_300.tsv"
SUPPORTED = DATA / "icl_satellite_supported_94_association_2026-11-05_01/supported_pairs.tsv"
VARIANCE = DATA / "icl_variance_radial_balance_full783_2026-11-06_01/target_control_radial_comparison.tsv"
VARIANCE_SUMMARY = DATA / "icl_variance_radial_balance_full783_2026-11-06_01/summary.json"
TARGET_MANIFEST = DATA / "sample_manifest.json"
CONTROL_MANIFEST = DATA / "control_download_01/download_manifest.json"
DEFAULT_OUTPUT = DATA / "icl_satellite_ipw_variance_sensitivity_2026-11-06_01"
EXPECTED_SHA = {
    "support": "79f7c4c6b18b29ad837e91da4a4359433b09a564cad62efaf60b4cce655685da",
    "variance": "db77c9a32797aafa79efaac7f7e18d7e14a8f79c96b18e8fe296aceb6f35830d",
}
MODELS = (
    "primary_support_covariates",
    "plus_log_variance_ratio",
    "plus_variance_ratio_and_mean_depth",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def fit_logit(design: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray, dict]:
    if design.ndim != 2 or y.shape != (design.shape[0],) or len(np.unique(y)) != 2:
        raise ValueError("selection resample has invalid design or only one support class")

    def objective(beta: np.ndarray) -> tuple[float, np.ndarray]:
        eta = design @ beta
        probability = expit(eta)
        return (float(np.logaddexp(0.0, eta).sum() - np.dot(y, eta)),
                design.T @ (probability - y))

    opt = minimize(lambda beta: objective(beta), np.zeros(design.shape[1]), jac=True,
                   method="BFGS", options={"gtol": 1e-7, "maxiter": 2000})
    grad = float(np.linalg.norm(opt.jac, ord=np.inf))
    if not np.isfinite(grad) or grad > 2e-5:
        raise RuntimeError(f"logistic fit failed stationarity: {opt.message}; gradient={grad:.3g}")
    probability = np.clip(expit(design @ opt.x), 1e-6, 1 - 1e-6)
    return opt.x, probability, {
        "converged_by_stationarity": True,
        "optimizer_reported_success": bool(opt.success),
        "optimizer_message": str(opt.message),
        "gradient_inf_norm": grad,
        "iterations": int(opt.nit),
    }


def weighted_partial_r(rows: list[dict], weights: np.ndarray) -> float:
    x = np.asarray([float(row["icl_delta_200_300"]) for row in rows])
    y = np.asarray([float(row["satellite_delta_50_75"]) for row in rows])
    z = np.asarray([float(row["redshift"]) for row in rows])
    richness = np.asarray([float(row["richness"]) for row in rows])
    design = np.column_stack((np.ones(len(rows)), z, richness))
    if (len(rows) < 6 or not np.isfinite(np.column_stack((x, y, z, richness, weights))).all()
            or np.any(weights <= 0) or np.linalg.matrix_rank(design) < design.shape[1]):
        return float("nan")
    root_w = np.sqrt(weights)
    bx = np.linalg.lstsq(design * root_w[:, None], x * root_w, rcond=None)[0]
    by = np.linalg.lstsq(design * root_w[:, None], y * root_w, rcond=None)[0]
    rx, ry = x - design @ bx, y - design @ by
    rx -= np.average(rx, weights=weights)
    ry -= np.average(ry, weights=weights)
    den = np.sqrt(np.dot(weights, rx * rx) * np.dot(weights, ry * ry))
    return float(np.dot(weights, rx * ry) / den) if den > 0 else float("nan")


def make_design(rows: list[dict], model: str, indices: np.ndarray | None = None) -> np.ndarray:
    source = rows if indices is None else [rows[int(i)] for i in indices]
    columns = [np.ones(len(source)),
               np.asarray([row["filter"] == "i2" for row in source], dtype=float),
               np.asarray([float(row["z_s"]) for row in source]),
               np.asarray([float(row["log_richness_s"]) for row in source]),
               np.asarray([float(row["coverage_s"]) for row in source]),
               np.asarray([float(row["center_masked_k2p5"]) for row in source])]
    if model in MODELS[1:]:
        columns.append(np.asarray([float(row["log_variance_ratio_s"]) for row in source]))
    if model == MODELS[2]:
        columns.append(np.asarray([float(row["mean_log_variance_s"]) for row in source]))
    return np.column_stack(columns)


def association(rows: list[dict], selected: np.ndarray, probabilities: np.ndarray) -> dict:
    by_filter = {}
    for band in ("i", "i2"):
        all_idx = np.asarray([i for i, row in enumerate(rows) if row["filter"] == band])
        sel_idx = all_idx[selected[all_idx]]
        prevalence = float(selected[all_idx].mean())
        weights = prevalence / probabilities[sel_idx]
        weights /= weights.mean()
        selected_rows = [rows[i] for i in sel_idx]
        by_filter[band] = {
            "n_icl_supported": int(len(all_idx)),
            "n_satellite_supported": int(len(sel_idx)),
            "support_fraction": prevalence,
            "ipw_partial_r_z_richness": weighted_partial_r(selected_rows, weights),
            "effective_n": float(weights.sum() ** 2 / np.dot(weights, weights)),
            "max_stabilized_weight": float(weights.max()),
            "selected_propensity_min": float(probabilities[sel_idx].min()),
            "selected_propensity_max": float(probabilities[sel_idx].max()),
        }
    return by_filter


def run(output: Path, replicas: int, seed: int) -> dict:
    output = output.resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite {output}")
    inputs = (ICL_METRICS, SUPPORTED, VARIANCE, VARIANCE_SUMMARY, TARGET_MANIFEST, CONTROL_MANIFEST)
    for path in inputs:
        if not path.is_file():
            raise FileNotFoundError(path)
    if sha256(SUPPORTED) != EXPECTED_SHA["support"] or sha256(VARIANCE) != EXPECTED_SHA["variance"]:
        raise ValueError("pinned support or full variance input changed")
    variance_summary = json.loads(VARIANCE_SUMMARY.read_text(encoding="utf-8"))
    if variance_summary.get("science_measurement") is not False:
        raise ValueError("variance input must remain diagnostic-only")

    metrics = read_tsv(ICL_METRICS)
    support_rows = read_tsv(SUPPORTED)
    variance_rows = [row for row in read_tsv(VARIANCE)
                     if float(row["r_min_kpc"]) == 200.0 and float(row["r_max_kpc"]) == 300.0]
    targets_doc = json.loads(TARGET_MANIFEST.read_text(encoding="utf-8"))
    controls_doc = json.loads(CONTROL_MANIFEST.read_text(encoding="utf-8"))
    targets = {str(row["id"]): row for row in targets_doc["clusters"]}
    controls = {str(row["target_cluster_id"]): row for row in controls_doc["clusters"]}
    support_by_key = {(row["cluster_id"], row["filter"]): row for row in support_rows}
    variance_by_key = {(row["cluster_id"], row["filter"]): row for row in variance_rows}
    if (len(metrics) != 224 or len({(row["id"], row["filter"]) for row in metrics}) != 224
            or len(support_by_key) != 94 or len(variance_by_key) != len(variance_rows)):
        raise ValueError("unexpected or duplicate support/variance parent rows")

    rows = []
    for source in sorted(metrics, key=lambda row: (int(row["id"]), row["filter"])):
        ident, band = source["id"], source["filter"]
        if ident not in targets or ident not in controls:
            raise ValueError(f"missing matched target/control manifest entry for {ident}")
        var = variance_by_key.get((ident, band))
        if var is None or var["control_id"] != str(controls[ident]["id"]):
            raise ValueError(f"missing or mismatched variance pair for {ident} {band}")
        ratio = float(var["log_target_control_variance_ratio"])
        tv, cv = float(var["target_variance"]), float(var["control_variance"])
        if not np.isfinite([ratio, tv, cv]).all() or min(tv, cv) <= 0:
            raise ValueError(f"invalid variance metrics for {ident} {band}")
        sat = support_by_key.get((ident, band))
        z, richness = float(source["z"]), float(source["richness"])
        row = {
            "cluster_id": ident, "filter": band,
            "control_id": str(controls[ident]["id"]),
            "selected_50_75": sat is not None,
            "icl_delta_200_300": float(sat["icl_delta_200_300"]) if sat else float("nan"),
            "satellite_delta_50_75": float(sat["satellite_delta_50_75"]) if sat else float("nan"),
            "redshift": z, "richness": richness,
            "icl_coverage_200_300": float(source["coverage_200_300"]),
            "center_masked_k2p5": source["center_masked_k2p5"] == "True",
            "log_variance_ratio": ratio,
            "target_variance_200_300": tv,
            "control_variance_200_300": cv,
            "mean_log_variance": 0.5 * (np.log(tv) + np.log(cv)),
            "target_ra_deg": float(targets[ident]["ra_deg"]),
            "target_dec_deg": float(targets[ident]["dec_deg"]),
            "control_ra_deg": float(controls[ident]["ra_deg"]),
            "control_dec_deg": float(controls[ident]["dec_deg"]),
        }
        if sat and (abs(row["icl_delta_200_300"] - float(source["delta_sb_200_300"])) > 5e-7
                    or abs(z - float(sat["redshift"])) > 1e-8
                    or abs(richness - float(sat["richness"])) > 1e-6):
            raise ValueError(f"selected outcome mismatch for {ident} {band}")
        rows.append(row)

    if sum(bool(row["selected_50_75"]) for row in rows) != 94:
        raise ValueError("expected exact 94 selected satellite profiles")
    for band, expected in (("i", 71), ("i2", 153)):
        if sum(row["filter"] == band for row in rows) != expected:
            raise ValueError(f"unexpected parent count in {band}")

    def std(values: np.ndarray) -> np.ndarray:
        scale = float(values.std(ddof=0))
        if not np.isfinite(scale) or scale <= 0:
            raise ValueError("non-varying support-model covariate")
        return (values - values.mean()) / scale

    z_s = std(np.asarray([row["redshift"] for row in rows]))
    lr_s = std(np.log([row["richness"] for row in rows]))
    coverage_s = std(np.asarray([row["icl_coverage_200_300"] for row in rows]))
    ratio_s = std(np.asarray([row["log_variance_ratio"] for row in rows]))
    depth_s = std(np.asarray([row["mean_log_variance"] for row in rows]))
    for i, row in enumerate(rows):
        row.update({"z_s": z_s[i], "log_richness_s": lr_s[i], "coverage_s": coverage_s[i],
                    "log_variance_ratio_s": ratio_s[i], "mean_log_variance_s": depth_s[i]})
    y = np.asarray([row["selected_50_75"] for row in rows], dtype=float)
    labels_target = [{"ra_deg": row["target_ra_deg"], "dec_deg": row["target_dec_deg"]} for row in rows]
    labels_control = [{"ra_deg": row["control_ra_deg"], "dec_deg": row["control_dec_deg"]} for row in rows]
    _, block_labels, target_cells, control_cells, n_blocks, largest = pair_block_labels(
        labels_target, labels_control, np.ones(len(rows), dtype=bool))

    fits = {}
    probabilities = {}
    design_names = {
        MODELS[0]: ["intercept", "filter_i2", "redshift_per_SD", "log_richness_per_SD",
                    "icl_coverage_per_SD", "center_masked"],
        MODELS[1]: ["intercept", "filter_i2", "redshift_per_SD", "log_richness_per_SD",
                    "icl_coverage_per_SD", "center_masked", "log_variance_ratio_per_SD"],
        MODELS[2]: ["intercept", "filter_i2", "redshift_per_SD", "log_richness_per_SD",
                    "icl_coverage_per_SD", "center_masked", "log_variance_ratio_per_SD",
                    "mean_log_variance_per_SD"],
    }
    for model in MODELS:
        beta, p, meta = fit_logit(make_design(rows, model), y)
        fits[model] = {"coefficients": dict(zip(design_names[model], map(float, beta))), "fit": meta}
        probabilities[model] = p
        for row, value in zip(rows, p):
            row[f"propensity_{model}"] = float(value)

    point_estimates = {model: association(rows, y.astype(bool), probabilities[model]) for model in MODELS}
    # Paired resampling keeps the baseline and variance-augmented IPW estimates
    # on the exact same connected target/control blocks.
    groups = np.unique(block_labels)
    by_block = {group: np.flatnonzero(block_labels == group) for group in groups}
    rng = np.random.default_rng(seed)
    draws = {model: {band: [] for band in ("i", "i2")} for model in MODELS}
    delta_draws = {model: {band: [] for band in ("i", "i2")} for model in MODELS[1:]}
    fit_failures = {model: 0 for model in MODELS}
    for _ in range(replicas):
        indices = np.concatenate([by_block[group] for group in rng.choice(groups, size=len(groups), replace=True)])
        local_y = y[indices]
        if len(np.unique(local_y)) < 2:
            for model in MODELS:
                fit_failures[model] += 1
            continue
        estimates = {model: {} for model in MODELS}
        for model in MODELS:
            try:
                _, p_draw, _ = fit_logit(make_design(rows, model, indices), local_y)
                draw_rows = [rows[int(i)] for i in indices]
                local_selected = local_y.astype(bool)
                estimated = association(draw_rows, local_selected, p_draw)
            except (ValueError, RuntimeError, np.linalg.LinAlgError):
                fit_failures[model] += 1
                continue
            estimates[model] = estimated
            for band in ("i", "i2"):
                value = estimated[band]["ipw_partial_r_z_richness"]
                if np.isfinite(value):
                    draws[model][band].append(float(value))
        for model in MODELS[1:]:
            for band in ("i", "i2"):
                base = estimates[MODELS[0]].get(band, {}).get("ipw_partial_r_z_richness", np.nan)
                current = estimates[model].get(band, {}).get("ipw_partial_r_z_richness", np.nan)
                if np.isfinite(base) and np.isfinite(current):
                    delta_draws[model][band].append(float(current - base))

    results = {}
    for model in MODELS:
        results[model] = {"support_model": fits[model], "association_by_filter": {}}
        for band in ("i", "i2"):
            values = np.asarray(draws[model][band])
            if len(values) < replicas * 0.8:
                raise RuntimeError(f"too few {model}/{band} valid paired draws: {len(values)}/{replicas}")
            record = dict(point_estimates[model][band])
            record.update({"paired_spatial_block_bootstrap_ci95": np.percentile(values, [2.5, 97.5]).tolist(),
                           "bootstrap_valid": len(values)})
            if model != MODELS[0]:
                deltas = np.asarray(delta_draws[model][band])
                if len(deltas) < replicas * 0.8:
                    raise RuntimeError(f"too few paired change draws for {model}/{band}")
                record["change_from_primary"] = (
                    point_estimates[model][band]["ipw_partial_r_z_richness"]
                    - point_estimates[MODELS[0]][band]["ipw_partial_r_z_richness"])
                record["paired_change_ci95"] = np.percentile(deltas, [2.5, 97.5]).tolist()
                record["change_bootstrap_valid"] = len(deltas)
            results[model]["association_by_filter"][band] = record

    output.mkdir(parents=True, exist_ok=False)
    joined_path = output / "support_variance_join.tsv"
    fields = list(rows[0]) + ["connected_5deg_block"]
    with joined_path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows({**row, "connected_5deg_block": str(label)} for row, label in zip(rows, block_labels))
    summary = {
        "schema": "astrafek-icl-satellite-ipw-variance-sensitivity-v1",
        "status": "complete_offline_sensitivity",
        "science_measurement": False,
        "hsc_network_used": False,
        "parent_icl_supported_pairs": len(rows),
        "selected_satellite_profile_pairs": int(y.sum()),
        "filter_counts": {band: sum(row["filter"] == band for row in rows) for band in ("i", "i2")},
        "models": results,
        "block_bootstrap": {
            "method": "same connected 5-degree target/control components resampled for all propensity models; each selection model refit per replicate",
            "n_connected_blocks": int(n_blocks), "largest_block_pairs": int(largest),
            "target_cells": int(target_cells), "control_cells": int(control_cells),
            "replicas_requested": replicas, "seed": seed,
            "selection_model_nonconvergent_draws": fit_failures,
        },
        "selection_model_limits": [
            "IPW is a sensitivity analysis; missing-at-random conditional on measured covariates is not established",
            "the delivered variance plane is not total noise and cannot represent source-Poisson/effective-gain terms",
            "the model does not estimate image-level ICL estimator covariance or control all PSF, sky, mask and completeness effects",
        ],
        "input_sha256": {str(path): sha256(path) for path in inputs},
        "output_sha256": {str(joined_path): sha256(joined_path)},
    }
    summary_path = output / "summary.json"
    summary_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True, allow_nan=False) + "\n",
                            encoding="utf-8")
    summary["summary_sha256"] = sha256(summary_path)
    print(json.dumps({
        "status": summary["status"], "science_measurement": False,
        "n_parent": len(rows), "n_selected": int(y.sum()), "blocks": n_blocks,
        "results": {model: {band: data["ipw_partial_r_z_richness"]
                            for band, data in item["association_by_filter"].items()}
                    for model, item in results.items()},
        "output": str(output), "summary_sha256": summary["summary_sha256"],
    }, indent=2))
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--replicas", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=20261106)
    args = parser.parse_args()
    if args.replicas < 1000:
        parser.error("at least 1,000 spatial-block bootstrap replicates are required")
    run(args.output, args.replicas, args.seed)


if __name__ == "__main__":
    main()
