#!/usr/bin/env python3
"""Fit cluster-robust multivariable models for annulus support selection.

Fits unpenalized logistic regressions to saved candidate-level diagnostics,
with a cluster-ID sandwich covariance. This is a descriptive selection model,
not a completeness correction or science measurement.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.optimize import minimize
from scipy.special import expit
from scipy.stats import t as student_t


DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
INPUT_DIR = DATA / "icl_satellite_mask_selection_2026-11-05_01"
INPUT = INPUT_DIR / "candidate_selection_inputs.tsv"
DEFAULT_OUTPUT = DATA / "icl_satellite_mask_selection_adjusted_2026-11-05_01"
FILTERS = ("i", "i2")
PREDICTORS = ("cluster_radius_kpc", "log_detection_snr", "redshift", "log_richness")
OUTCOMES = ("coverage_pass", "profile_retained")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    import csv
    with Path(path).open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def fit_cluster_logit(x: np.ndarray, y: np.ndarray, clusters: np.ndarray) -> dict:
    if x.ndim != 2 or y.shape != (x.shape[0],) or clusters.shape != y.shape:
        raise ValueError("invalid logistic input shapes")
    n, p = x.shape
    unique_clusters = np.unique(clusters)
    g = len(unique_clusters)
    if n <= p or g < 10 or set(np.unique(y)) != {0., 1.}:
        raise ValueError("insufficient observations, clusters, or outcome variation")

    def objective(beta: np.ndarray) -> tuple[float, np.ndarray]:
        eta = x @ beta
        # Stable Bernoulli negative log-likelihood and score.
        loss = np.logaddexp(0., eta).sum() - np.dot(y, eta)
        probability = expit(eta)
        gradient = x.T @ (probability - y)
        return float(loss), gradient

    optimized = minimize(lambda beta: objective(beta), np.zeros(p), jac=True,
                         method="BFGS", options={"gtol": 1e-6, "maxiter": 2000})
    beta = optimized.x
    if not optimized.success or np.linalg.norm(optimized.jac, ord=np.inf) > 1e-5:
        raise RuntimeError(
            f"cluster-logit optimization did not converge: {optimized.message}; "
            f"gradient_inf_norm={np.linalg.norm(optimized.jac, ord=np.inf):.3g}"
        )
    probability = expit(x @ beta)
    weights = probability * (1. - probability)
    information = x.T @ (weights[:, None] * x)
    bread = np.linalg.pinv(information, rcond=1e-12)
    meat = np.zeros((p, p), dtype=float)
    for group in unique_clusters:
        selected = clusters == group
        score = x[selected].T @ (y[selected] - probability[selected])
        meat += np.outer(score, score)
    correction = (g / (g - 1.)) * ((n - 1.) / (n - p))
    covariance = correction * bread @ meat @ bread
    errors = np.sqrt(np.maximum(np.diag(covariance), 0.))
    critical = float(student_t.ppf(.975, g - 1))
    intervals = np.column_stack((beta - critical * errors, beta + critical * errors))
    names = ("intercept", *PREDICTORS)
    coefficients = {}
    for index, name in enumerate(names):
        lo, hi = intervals[index]
        coefficients[name] = {
            "log_odds_per_predictor_sd": float(beta[index]),
            "cluster_robust_se": float(errors[index]),
            "ci95_log_odds": [float(lo), float(hi)],
            "odds_ratio": float(np.exp(beta[index])),
            "ci95_odds_ratio": [float(np.exp(lo)), float(np.exp(hi))],
        }
    predicted = probability >= .5
    return {
        "n_candidates": n, "n_cluster_blocks": g,
        "outcome_fraction": float(y.mean()),
        "optimizer_success": bool(optimized.success),
        "optimizer_message": str(optimized.message),
        "gradient_inf_norm": float(np.linalg.norm(optimized.jac, ord=np.inf)),
        "information_condition_number": float(np.linalg.cond(information)),
        "brier_score": float(np.mean((probability - y) ** 2)),
        "threshold_0p5_accuracy_descriptive": float(np.mean(predicted == y)),
        "coefficients": coefficients,
        "covariance_method": "cluster-ID sandwich, finite-sample correction, t critical value with G-1 degrees of freedom",
    }


def run(output: Path) -> dict:
    if output.exists():
        raise FileExistsError(f"refusing to overwrite {output}")
    if not INPUT.is_file():
        raise FileNotFoundError(INPUT)
    rows = read_tsv(INPUT)
    derived = []
    for row in rows:
        snr = max(float(row["detection_snr"]), 0.)
        richness = float(row["richness"])
        if richness <= 0:
            raise ValueError("richness must be positive for log transform")
        derived.append({
            **row,
            "log_detection_snr": float(np.log1p(snr)),
            "log_richness": float(np.log(richness)),
        })

    results = {}
    standardization = {}
    for band in FILTERS:
        selected = [row for row in derived if row["filter"] == band]
        raw = np.asarray([[float(row[name]) for name in PREDICTORS] for row in selected])
        means = raw.mean(axis=0)
        scales = raw.std(axis=0, ddof=0)
        if np.any(~np.isfinite(raw)) or np.any(scales <= 0):
            raise ValueError(f"non-finite or invariant predictor in {band}")
        standardized = (raw - means) / scales
        design = np.column_stack((np.ones(len(selected)), standardized))
        clusters = np.asarray([row["cluster_id"] for row in selected], dtype=str)
        results[band] = {}
        for outcome in OUTCOMES:
            y = np.asarray([float(row[outcome]) for row in selected])
            fit = fit_cluster_logit(design, y, clusters)
            results[band][outcome] = fit
        standardization[band] = {
            name: {"mean": float(mean), "sd": float(scale)}
            for name, mean, scale in zip(PREDICTORS, means, scales)
        }

    output.mkdir(parents=True, exist_ok=False)
    summary = {
        "schema": "astrafek-icl-satellite-mask-selection-adjusted-v1",
        "status": "complete_offline_diagnostic", "science_measurement": False,
        "hsc_network_used": False, "science_image_pixels_read": False,
        "input_sha256": {str(INPUT): sha256(INPUT)},
        "script_sha256": sha256(Path(__file__).resolve()),
        "predictors_per_sd": list(PREDICTORS),
        "transformations": {"detection_snr": "log1p(max(SNR_NOMINAL,0))", "richness": "natural_log"},
        "predictor_standardization": standardization,
        "results": results,
        "limitations": [
            "Coefficients are adjusted descriptive associations under a linear logit specification.",
            "Cluster-robust intervals account for candidate dependence within cluster ID, but not PSF, sky, calibration, or mask-model uncertainty.",
            "The model is not validated against injection/recovery and must not be used as an inverse-probability science weight yet.",
            "The two filters are modeled separately; coefficient intervals are not adjusted for multiple comparisons.",
        ],
    }
    summary_path = output / "summary.json"
    summary_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True)+"\n",
                            encoding="utf-8")
    compact = {}
    for band, fits in results.items():
        compact[band] = {
            outcome: {name: {
                "odds_ratio": fit["coefficients"][name]["odds_ratio"],
                "ci95": fit["coefficients"][name]["ci95_odds_ratio"],
            } for name in PREDICTORS}
            for outcome, fit in fits.items()
        }
    print(json.dumps({"status": summary["status"], "output": str(output),
                      "adjusted_odds_ratios": compact}, ensure_ascii=False, indent=2))
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    run(args.output.resolve())


if __name__ == "__main__":
    main()
