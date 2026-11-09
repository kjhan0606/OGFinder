#!/usr/bin/env python3
"""Assess sensitivity of the ICL--satellite association to delivered variance."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from analyze_icl_satellite_supported_94 import read_tsv, sha256
from evaluate_prestack_pair_support import pair_block_labels
from audit_icl_target_control_outer_sky import partial_pearson

ROOT = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
SUPPORT_DIR = ROOT / "icl_satellite_supported_94_association_2026-11-05_01"
SUPPORT = SUPPORT_DIR / "supported_pairs.tsv"
SUPPORT_SUMMARY = SUPPORT_DIR / "summary.json"
VARIANCE_DIR = ROOT / "icl_variance_radial_balance_full783_2026-11-06_01"
VARIANCE = VARIANCE_DIR / "target_control_radial_comparison.tsv"
VARIANCE_SUMMARY = VARIANCE_DIR / "summary.json"
DEFAULT_OUTPUT = ROOT / "icl_variance_adjusted_satellite_assoc_2026-11-06_01"
EXPECTED_SUPPORT_SHA = "79f7c4c6b18b29ad837e91da4a4359433b09a564cad62efaf60b4cce655685da"
EXPECTED_VARIANCE_SHA = "db77c9a32797aafa79efaac7f7e18d7e14a8f79c96b18e8fe296aceb6f35830d"
MODEL_COVARIATES = {
    "baseline_z_richness": ("redshift", "richness"),
    "plus_log_variance_ratio": ("redshift", "richness", "log_variance_ratio"),
    "plus_variance_ratio_and_mean_depth": (
        "redshift", "richness", "log_variance_ratio", "mean_log_variance",
    ),
}


def connected_labels(rows: list[dict[str, str]]) -> np.ndarray:
    target = [{"ra_deg": float(row["target_ra_deg"]), "dec_deg": float(row["target_dec_deg"])}
              for row in rows]
    control = [{"ra_deg": float(row["control_ra_deg"]), "dec_deg": float(row["control_dec_deg"])}
               for row in rows]
    _, labels, *_ = pair_block_labels(target, control, np.ones(len(rows), dtype=bool))
    return np.asarray(labels, dtype=str)


def estimate(rows: list[dict[str, str]], model: str, indices: np.ndarray | None = None) -> float:
    chosen = rows if indices is None else [rows[int(i)] for i in indices]
    covariates = [np.asarray([float(row[field]) for row in chosen])
                  for field in MODEL_COVARIATES[model]]
    x = np.asarray([float(row["icl_delta_200_300"]) for row in chosen])
    y = np.asarray([float(row["satellite_delta_50_75"]) for row in chosen])
    return partial_pearson(x, y, covariates)


def bootstrap_models(rows: list[dict[str, str]], labels: np.ndarray,
                     replicas: int, seed: int) -> dict:
    blocks = np.unique(labels)
    if len(blocks) < 6:
        raise ValueError(f"at least six connected spatial blocks required; got {len(blocks)}")
    by_block = {block: np.flatnonzero(labels == block) for block in blocks}
    rng = np.random.default_rng(seed)
    draws = {model: [] for model in MODEL_COVARIATES}
    deltas = {model: [] for model in tuple(MODEL_COVARIATES)[1:]}
    for _ in range(replicas):
        selected = rng.choice(blocks, size=len(blocks), replace=True)
        indices = np.concatenate([by_block[block] for block in selected])
        estimates = {model: estimate(rows, model, indices) for model in MODEL_COVARIATES}
        if not all(np.isfinite(value) for value in estimates.values()):
            continue
        for model, value in estimates.items():
            draws[model].append(float(value))
        baseline = estimates["baseline_z_richness"]
        for model in deltas:
            deltas[model].append(float(estimates[model] - baseline))
    if len(draws["baseline_z_richness"]) < replicas * 0.8:
        raise ValueError(f"fewer than 80% valid joint block-bootstrap draws: {len(draws['baseline_z_richness'])}/{replicas}")
    result = {}
    for model in MODEL_COVARIATES:
        values = np.asarray(draws[model])
        result[model] = {
            "partial_r": float(estimate(rows, model)),
            "paired_block_bootstrap_ci95": np.percentile(values, [2.5, 97.5]).tolist(),
            "replicas_valid": len(values),
        }
    for model, values in deltas.items():
        draws_array = np.asarray(values)
        result[model]["change_from_baseline"] = float(
            result[model]["partial_r"] - result["baseline_z_richness"]["partial_r"]
        )
        result[model]["paired_change_ci95"] = np.percentile(
            draws_array, [2.5, 97.5]
        ).tolist()
    result["resampling"] = {
        "n_connected_5deg_blocks": int(len(blocks)),
        "largest_block_pairs": int(max(len(indices) for indices in by_block.values())),
        "replicas_requested": replicas,
        "replicas_valid_joint": len(draws["baseline_z_richness"]),
        "seed": seed,
        "rule": "resample connected components sharing target or control equal-area 5-degree sky cells; same draws for all models",
    }
    return result


def run(output: Path, replicas: int, seed: int) -> dict:
    output = output.resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite existing result directory: {output}")
    for path in (SUPPORT, SUPPORT_SUMMARY, VARIANCE, VARIANCE_SUMMARY):
        if not path.is_file():
            raise FileNotFoundError(path)
    if sha256(SUPPORT) != EXPECTED_SUPPORT_SHA or sha256(VARIANCE) != EXPECTED_VARIANCE_SHA:
        raise ValueError("pinned supported-pair or full-sample variance input has changed")
    support_summary = json.loads(SUPPORT_SUMMARY.read_text(encoding="utf-8"))
    variance_summary = json.loads(VARIANCE_SUMMARY.read_text(encoding="utf-8"))
    if support_summary.get("status") != "complete_diagnostic" or support_summary.get("science_measurement") is not False:
        raise ValueError("94-pair support parent is not a verified diagnostic")
    if variance_summary.get("science_measurement") is not False:
        raise ValueError("variance parent must remain diagnostic-only")

    support = read_tsv(SUPPORT)
    comparison = read_tsv(VARIANCE)
    if len(support) != 94 or any(row.get("science_measurement") != "False" for row in support):
        raise ValueError("expected 94 supported diagnostic pairs")
    selected = [row for row in comparison
                if float(row["r_min_kpc"]) == 200.0 and float(row["r_max_kpc"]) == 300.0]
    var_by_key = {}
    for row in selected:
        key = (row["cluster_id"], row["filter"], row["control_id"])
        if key in var_by_key:
            raise ValueError(f"duplicate variance row in 200-300 kpc bin: {key}")
        var_by_key[key] = row
    joined = []
    for row in support:
        key = (row["cluster_id"], row["filter"], row["control_id"])
        var = var_by_key.get(key)
        if var is None:
            raise ValueError(f"missing 200-300 kpc variance comparison for supported pair: {key}")
        for field in ("icl_delta_200_300", "satellite_delta_50_75", "redshift"):
            if abs(float(row[field]) - float(var[field])) > 1e-8:
                raise ValueError(f"support/variance parent values disagree for {key}: {field}")
        ratio = float(var["log_target_control_variance_ratio"])
        target_variance = float(var["target_variance"])
        control_variance = float(var["control_variance"])
        if not np.isfinite([ratio, target_variance, control_variance]).all() or min(target_variance, control_variance) <= 0:
            raise ValueError(f"invalid variance metrics for {key}")
        joined.append({
            **row,
            "log_variance_ratio": str(ratio),
            "mean_log_variance": str(0.5 * (np.log(target_variance) + np.log(control_variance))),
            "target_variance_200_300": str(target_variance),
            "control_variance_200_300": str(control_variance),
            "science_measurement": "False",
        })
    joined.sort(key=lambda row: (row["filter"], int(row["cluster_id"])))
    if len(joined) != 94 or len(var_by_key) < 94:
        raise ValueError(f"expected exact 94-pair join, got {len(joined)}")
    filter_counts = {band: sum(row["filter"] == band for row in joined) for band in ("i", "i2")}
    if filter_counts != {"i": 26, "i2": 68}:
        raise ValueError(f"unexpected filter counts: {filter_counts}")

    results = {}
    for band_index, band in enumerate(("i", "i2")):
        rows = [row for row in joined if row["filter"] == band]
        labels = connected_labels(rows)
        results[band] = {
            "n_pairs": len(rows),
            "variance_ratio_summary": {
                "median_target_control_ratio": float(np.median(np.exp([float(row["log_variance_ratio"]) for row in rows]))),
                "median_log_ratio": float(np.median([float(row["log_variance_ratio"]) for row in rows])),
            },
            "models": bootstrap_models(rows, labels, replicas, seed + band_index * 100),
        }

    output.mkdir(parents=True, exist_ok=False)
    joined_path = output / "joined_supported_pairs.tsv"
    with joined_path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(joined[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(joined)
    summary = {
        "schema": "astrafek-icl-variance-adjusted-satellite-association-v1",
        "status": "complete_diagnostic",
        "science_measurement": False,
        "hsc_query_submitted": False,
        "image_pixels_read": False,
        "n_joined_pairs": len(joined),
        "filter_counts": filter_counts,
        "variance_bin_kpc": [200, 300],
        "models": {
            "baseline_z_richness": "partial Pearson r of ICL 200-300 kpc versus satellite 50-75 kpc contrasts, controlling redshift and richness",
            "plus_log_variance_ratio": "baseline plus log(target/control delivered VARIANCE-plane median ratio) at 200-300 kpc",
            "plus_variance_ratio_and_mean_depth": "previous model plus mean(log(target variance), log(control variance)) as a coarse pair-depth proxy",
        },
        "results_by_filter": results,
        "inputs_sha256": {str(path): sha256(path) for path in (SUPPORT, SUPPORT_SUMMARY, VARIANCE, VARIANCE_SUMMARY)},
        "outputs_sha256": {str(joined_path): sha256(joined_path)},
        "limitations": [
            "the selected support sample contains only 26 i and 68 i2 pairs; model-adjusted correlations are low-power sensitivity diagnostics",
            "the delivered variance plane omits the unknown effective-gain/source-Poisson term and is not a calibrated total-noise estimate",
            "adjusting for variance balance does not prove confounding, mediation, or causation",
            "mean pair variance is only a coarse depth proxy; it does not measure the ICL estimator covariance, coadd weights, sky transfer, PSF wings, or completeness",
            "connected 5-degree cells are approximate spatial blocks, not validated contiguous HSC fields",
        ],
    }
    summary_path = output / "summary.json"
    summary_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    summary["summary_sha256"] = sha256(summary_path)
    print(json.dumps({
        "status": summary["status"], "science_measurement": False,
        "results": {band: {model: value["partial_r"] for model, value in result["models"].items()
                            if model in MODEL_COVARIATES}
                     for band, result in results.items()},
        "output": str(output), "summary_sha256": summary["summary_sha256"],
    }, indent=2))
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--replicas", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=20261106)
    args = parser.parse_args()
    if args.replicas < 1000:
        parser.error("at least 1,000 block-bootstrap replicates are required")
    run(args.output, args.replicas, args.seed)


if __name__ == "__main__":
    main()
