#!/usr/bin/env python3
"""Jointly propagate satellite-support IPW and finite-lag pixel-noise correction."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from analyze_icl_satellite_ipw_variance_sensitivity import (
    DEFAULT_OUTPUT as IPW_VARIANCE_OUTPUT,
    MODELS,
    fit_logit,
    make_design,
    read_tsv,
    sha256,
)
from assess_icl_satellite_pixel_noise_bias import corrected_partial_correlation
from verify_satellite_icl_covariance_support_scope import run as verify_covariance_scope

DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
DEFAULT_COVARIANCE = DATA / "satellite_icl_target_control_covariance_2026-11-03_02"
DEFAULT_OUTPUT = DATA / "icl_satellite_ipw_pixelnoise_joint_2026-11-06_01"
EXPECTED_COVARIANCE_SHA = "cc85ffd59fc94154a6597d3feca0e4d811682b9fcda7b79169e27ab778765167"


def estimate(rows: list[dict[str, str]], indices: np.ndarray, model: str,
             selected: np.ndarray, probabilities: np.ndarray,
             covariance: dict[tuple[str, str], dict[str, str]],
             band: str) -> dict | None:
    chosen = indices[(selected[indices]) & np.asarray([rows[i]["filter"] == band for i in indices])]
    if len(chosen) < 6:
        return None
    selected_rows = [rows[int(i)] for i in chosen]
    p = probabilities[chosen]
    # Stabilize within filter and resampled parent sample, then normalize to
    # unit mean so the estimator is unchanged by a common scale factor.
    prevalence = float(np.mean([selected[int(i)] for i in indices
                                if rows[int(i)]["filter"] == band]))
    weights = prevalence / p
    weights /= weights.mean()
    values = []
    for row in selected_rows:
        cov = covariance[(row["cluster_id"], row["filter"])]
        values.append((float(cov["icl_delta_lag_variance"]),
                       float(cov["target_satellite_lag_variance_njy2_arcsec4"]),
                       float(cov["target_icl_satellite_lag_covariance_njy2_arcsec4"])))
    values = np.asarray(values, dtype=float)
    corrected = corrected_partial_correlation(
        [float(row["icl_delta_200_300"]) for row in selected_rows],
        [float(row["satellite_delta_50_75"]) for row in selected_rows],
        [float(row["redshift"]) for row in selected_rows],
        [float(row["richness"]) for row in selected_rows],
        values[:, 0], values[:, 1], values[:, 2], weights=weights,
    )
    if corrected is None:
        return None
    return {
        **corrected,
        "n_selected_rows": len(chosen),
        "effective_n": float(weights.sum() ** 2 / np.dot(weights, weights)),
        "max_stabilized_weight": float(weights.max()),
    }


def run(output: Path, ipw_variance_dir: Path, covariance_dir: Path,
        replicas: int, seed: int) -> dict:
    output = output.resolve()
    ipw_variance_dir = ipw_variance_dir.resolve()
    covariance_dir = covariance_dir.resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite {output}")
    propensity_path = ipw_variance_dir / "support_variance_join.tsv"
    ipw_summary_path = ipw_variance_dir / "summary.json"
    covariance_path = covariance_dir / "target_control_pixel_covariance.tsv"
    covariance_summary_path = covariance_dir / "summary.json"
    for path in (propensity_path, ipw_summary_path, covariance_path, covariance_summary_path):
        if not path.is_file():
            raise FileNotFoundError(path)
    ipw_summary = json.loads(ipw_summary_path.read_text(encoding="utf-8"))
    covariance_summary = json.loads(covariance_summary_path.read_text(encoding="utf-8"))
    if (ipw_summary.get("status") != "complete_offline_sensitivity"
            or ipw_summary.get("science_measurement") is not False
            or ipw_summary.get("hsc_network_used") is not False):
        raise ValueError("parent propensity analysis is not an offline diagnostic")
    if (covariance_summary.get("status") != "complete_exploratory"
            or covariance_summary.get("science_measurement") is not False
            or sha256(covariance_path) != EXPECTED_COVARIANCE_SHA
            or covariance_summary.get("output_sha256", {}).get(str(covariance_path)) != EXPECTED_COVARIANCE_SHA):
        raise ValueError("pinned finite-lag covariance input is invalid or changed")

    rows = read_tsv(propensity_path)
    if len(rows) != 224 or sum(row["selected_50_75"] == "True" for row in rows) != 94:
        raise ValueError("expected 224 ICL-supported rows and exactly 94 satellite-supported outcomes")
    for row in rows:
        mask_value = row["center_masked_k2p5"].strip().lower()
        if mask_value not in {"true", "false", "1", "0"}:
            raise ValueError(f"unexpected center-mask state: {mask_value}")
        row["center_masked_k2p5"] = "1" if mask_value in {"true", "1"} else "0"
    selected = np.asarray([row["selected_50_75"] == "True" for row in rows], dtype=bool)
    labels = np.asarray([row["connected_5deg_block"] for row in rows], dtype=str)
    if len(np.unique(labels)) < 6:
        raise ValueError("too few connected spatial blocks")

    covariance_rows = read_tsv(covariance_path)
    covariance = {(row["cluster_id"], row["filter"]): row for row in covariance_rows}
    if len(covariance) != len(covariance_rows):
        raise ValueError("duplicate pair keys in pixel covariance table")
    selected_keys = {(row["cluster_id"], row["filter"]) for row, use in zip(rows, selected) if use}
    if selected_keys != set(covariance):
        raise ValueError("finite-lag covariance table does not exactly cover the 94 selected pairs")
    for key in selected_keys:
        ipw_row = next(row for row in rows if (row["cluster_id"], row["filter"]) == key)
        cov_row = covariance[key]
        if ipw_row["control_id"] != cov_row["control_id"]:
            raise ValueError(f"control identity mismatch in covariance join: {key}")
    scope = verify_covariance_scope(
        type("Args", (), {
            "icl_support": DATA / "icl_annulus_support_2026-10-18_01/pair_metrics_k2p5_c40_r50_300.tsv",
            "satellite_profiles": DATA / "satellite_outer_profile_perbin_full_2026-10-18_01/satellite_profiles.tsv",
            "association_rows": DATA / "icl_satellite_pixel_noise_bias_supported_94_2026-11-06_03/association_rows.tsv",
            "target_covariance": DATA / "satellite_icl_lag_covariance_2026-11-03_02",
            "target_control_covariance": covariance_dir,
        })()
    )
    if (scope.get("status") != "verified" or scope.get("n_covariance_output_pairs_verified") != 94
            or scope.get("n_association_pairs_covered") != 94 or scope.get("missing_association_pairs")):
        raise ValueError("independent pixel-covariance support-scope audit failed")

    all_indices = np.arange(len(rows))
    point = {}
    fit_point = {}
    probabilities_by_model = {}
    for model in MODELS:
        design = make_design(rows, model)
        beta, probability, fit_meta = fit_logit(design, selected.astype(float))
        probabilities_by_model[model] = probability
        fit_point[model] = {"coefficients": [float(value) for value in beta], "fit": fit_meta}
        point[model] = {}
        for band in ("i", "i2"):
            result = estimate(rows, all_indices, model, selected, probability, covariance, band)
            if result is None or result["pixel_noise_corrected_partial_r"] is None:
                raise ValueError(f"invalid joint corrected estimate for {model}/{band}")
            point[model][band] = result

    groups = np.unique(labels)
    by_block = {group: np.flatnonzero(labels == group) for group in groups}
    rng = np.random.default_rng(seed)
    draws = {model: {band: {"raw": [], "corrected": []} for band in ("i", "i2")} for model in MODELS}
    deltas = {model: {band: [] for band in ("i", "i2")} for model in MODELS[1:]}
    invalid_fit_draws = {model: 0 for model in MODELS}
    invalid_estimate_draws = {model: 0 for model in MODELS}
    for _ in range(replicas):
        indices = np.concatenate([by_block[group] for group in rng.choice(groups, size=len(groups), replace=True)])
        local_selected = selected[indices]
        if len(np.unique(local_selected)) < 2:
            for model in MODELS:
                invalid_fit_draws[model] += 1
            continue
        draw_rows = [rows[int(index)] for index in indices]
        local_indices = np.arange(len(draw_rows))
        model_results = {model: {} for model in MODELS}
        model_failed = set()
        for model in MODELS:
            try:
                _, p_draw, _ = fit_logit(make_design(draw_rows, model), local_selected.astype(float))
            except (ValueError, RuntimeError, np.linalg.LinAlgError):
                invalid_fit_draws[model] += 1
                model_failed.add(model)
                continue
            for band in ("i", "i2"):
                result = estimate(draw_rows, local_indices, model, local_selected, p_draw, covariance, band)
                if result is None or result["pixel_noise_corrected_partial_r"] is None:
                    invalid_estimate_draws[model] += 1
                    continue
                model_results[model][band] = result
                for key, name in (("raw_partial_r", "raw"),
                                  ("pixel_noise_corrected_partial_r", "corrected")):
                    value = result.get(key)
                    if value is not None and np.isfinite(value):
                        draws[model][band][name].append(float(value))
        for model in MODELS[1:]:
            for band in ("i", "i2"):
                if (model in model_failed or MODELS[0] in model_failed
                        or band not in model_results.get(model, {})
                        or band not in model_results.get(MODELS[0], {})):
                    continue
                deltas[model][band].append(
                    model_results[model][band]["pixel_noise_corrected_partial_r"]
                    - model_results[MODELS[0]][band]["pixel_noise_corrected_partial_r"])

    results = {}
    for model in MODELS:
        results[model] = {"support_model": fit_point[model], "association_by_filter": {}}
        for band in ("i", "i2"):
            raw = np.asarray(draws[model][band]["raw"])
            corrected = np.asarray(draws[model][band]["corrected"])
            if min(len(raw), len(corrected)) < replicas * 0.8:
                raise RuntimeError(f"too few joint-correction bootstrap draws for {model}/{band}")
            record = dict(point[model][band])
            record["raw_ipw_partial_r_ci95"] = np.percentile(raw, [2.5, 97.5]).tolist()
            record["pixel_noise_corrected_ipw_partial_r_ci95"] = np.percentile(
                corrected, [2.5, 97.5]).tolist()
            record["bootstrap_valid"] = len(corrected)
            if model != MODELS[0]:
                delta = np.asarray(deltas[model][band])
                if len(delta) < replicas * 0.8:
                    raise RuntimeError(f"too few paired model-change draws for {model}/{band}")
                record["change_from_primary_corrected"] = (
                    point[model][band]["pixel_noise_corrected_partial_r"]
                    - point[MODELS[0]][band]["pixel_noise_corrected_partial_r"])
                record["paired_change_ci95"] = np.percentile(delta, [2.5, 97.5]).tolist()
                record["change_bootstrap_valid"] = len(delta)
            results[model]["association_by_filter"][band] = record

    output.mkdir(parents=True, exist_ok=False)
    scope_path = output / "covariance_scope_verification.json"
    scope_path.write_text(json.dumps(scope, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    summary = {
        "schema": "astrafek-icl-satellite-ipw-pixelnoise-joint-v1",
        "status": "complete_offline_sensitivity",
        "science_measurement": False,
        "hsc_network_used": False,
        "parent_pairs": 224,
        "selected_pairs_with_pixel_covariance": 94,
        "filter_counts": {band: {
            "parent": sum(row["filter"] == band for row in rows),
            "selected": sum(bool(use) and row["filter"] == band for row, use in zip(rows, selected)),
        } for band in ("i", "i2")},
        "support_models": MODELS,
        "results_by_model": results,
        "block_bootstrap": {
            "method": "resample connected 5-degree target/control blocks among 224 ICL-supported pairs; refit support model per draw; calculate IPW and finite-lag-noise-corrected partial r on the selected rows",
            "n_connected_blocks": int(len(groups)),
            "largest_block_pairs": int(max(map(len, by_block.values()))),
            "replicas_requested": replicas,
            "seed": seed,
            "invalid_support_model_draws": invalid_fit_draws,
            "invalid_corrected_estimate_draws": invalid_estimate_draws,
        },
        "limitations": [
            "the finite-lag pixel-noise correction conditions on saved stationary covariance estimates and does not propagate their estimation uncertainty",
            "IPW missing-at-random conditional on measured covariates is not established; weights do not prove or guarantee a completeness correction",
            "the noise correction does not include source Poisson/effective gain, sky transfer, PSF wings, calibration, or astrophysical covariance",
            "the 5-degree connected components are approximate survey blocks, not verified contiguous HSC fields",
            "all outcomes remain selected-sample exploratory associations; science_measurement stays false",
        ],
        "input_sha256": {
            str(propensity_path): sha256(propensity_path),
            str(ipw_summary_path): sha256(ipw_summary_path),
            str(covariance_path): sha256(covariance_path),
            str(covariance_summary_path): sha256(covariance_summary_path),
        },
        "covariance_scope_sha256": sha256(scope_path),
    }
    summary_path = output / "summary.json"
    summary_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True, allow_nan=False) + "\n",
                            encoding="utf-8")
    print(json.dumps({
        "status": summary["status"], "science_measurement": False,
        "results": {model: {band: {
            "corrected_r": val["pixel_noise_corrected_partial_r"],
            "corrected_ci95": val["pixel_noise_corrected_ipw_partial_r_ci95"],
            "paired_change_ci95": val.get("paired_change_ci95"),
        } for band, val in item["association_by_filter"].items()}
                    for model, item in results.items()},
        "output": str(output),
    }, indent=2))
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ipw-variance-dir", type=Path, default=IPW_VARIANCE_OUTPUT)
    parser.add_argument("--covariance-dir", type=Path, default=DEFAULT_COVARIANCE)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--replicas", type=int, default=2000)
    parser.add_argument("--seed", type=int, default=20261107)
    args = parser.parse_args()
    if args.replicas < 1000:
        parser.error("at least 1,000 connected-block replicates are required")
    run(args.output, args.ipw_variance_dir, args.covariance_dir, args.replicas, args.seed)


if __name__ == "__main__":
    main()
