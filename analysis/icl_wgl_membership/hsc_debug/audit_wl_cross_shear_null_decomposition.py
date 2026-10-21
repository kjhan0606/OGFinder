#!/usr/bin/env python3
"""Decompose the existing HSC cluster-only cross-shear null diagnostic.

This audit consumes only the saved tomographic summary. It does not read
source rows, contact HSC CAS, or turn the cluster-only stack into a science
measurement.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.stats import chi2, norm


DEFAULT_INPUT = Path(
    "/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/"
    "wl_tomographic_consistency_2026-10-20_01/summary.json"
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def holm(p_values):
    values = np.asarray(p_values, dtype=float)
    if values.ndim != 1 or not np.all(np.isfinite(values)) or np.any((values < 0) | (values > 1)):
        raise ValueError("Holm adjustment requires finite p-values in [0, 1]")
    order = np.argsort(values, kind="stable")
    adjusted = np.empty_like(values)
    running = 0.0
    count = len(values)
    for rank, index in enumerate(order):
        running = max(running, (count - rank) * values[index])
        adjusted[index] = min(1.0, running)
    return adjusted


def _joint_test(values, covariance):
    values = np.asarray(values, dtype=float)
    covariance = np.asarray(covariance, dtype=float)
    if covariance.shape != (len(values), len(values)):
        raise ValueError("joint-test covariance shape does not match the vector")
    rank = int(np.linalg.matrix_rank(covariance))
    if rank == 0:
        raise ValueError("joint-test covariance has zero rank")
    statistic = float(values @ np.linalg.pinv(covariance, hermitian=True) @ values)
    return {"chi2": statistic, "degrees_of_freedom": rank,
            "p_value": float(chi2.sf(statistic, rank))}


def audit(summary: dict, source_sha256: str) -> dict:
    if summary.get("schema") != "astrafek-hsc-y3-tomographic-consistency-jackknife-v1":
        raise ValueError("unexpected tomographic summary schema")
    if summary.get("science_measurement") is not False:
        raise ValueError("input summary must remain explicitly non-science")
    profiles = summary.get("profile")
    covariance = np.asarray(summary.get("covariance"), dtype=float)
    if not isinstance(profiles, list) or not profiles:
        raise ValueError("tomographic profile is missing")
    if covariance.ndim != 2 or covariance.shape[0] != covariance.shape[1]:
        raise ValueError("profile covariance must be square")
    if not np.all(np.isfinite(covariance)) or not np.allclose(covariance, covariance.T, rtol=1e-10, atol=1e-12):
        raise ValueError("profile covariance must be finite and symmetric")
    bins = sorted({int(row["source_zbin"]) for row in profiles})
    radial_keys = sorted({(float(row["r_min_mpc"]), float(row["r_max_mpc"])) for row in profiles})
    if bins != [1, 2, 3, 4] or len(radial_keys) < 1:
        raise ValueError("expected HSC Y3 source bins 1..4 and radial profile bins")
    n_bin, n_radial = len(bins), len(radial_keys)
    expected_dimension = n_bin * n_radial * 2
    if covariance.shape != (expected_dimension, expected_dimension):
        raise ValueError("covariance dimensions do not match bin/radius/tangential-cross profile order")
    by_key = {(int(row["source_zbin"]), float(row["r_min_mpc"]), float(row["r_max_mpc"])): row
              for row in profiles}
    if len(by_key) != expected_dimension // 2:
        raise ValueError("duplicate or missing tomographic profile rows")

    cross_indices = []
    cross_values = []
    cross_errors = []
    for source_bin in bins:
        for radial_index, (r_min, r_max) in enumerate(radial_keys):
            row = by_key[(source_bin, r_min, r_max)]
            index = (source_bin - 1) * n_radial * 2 + radial_index * 2 + 1
            value = float(row["cross_delta_sigma_weak_msun_mpc2"])
            error = float(row["cross_jackknife_error_msun_mpc2"])
            if not np.isfinite(value) or not np.isfinite(error) or error <= 0:
                raise ValueError("cross-shear profile and error must be finite with positive error")
            if not np.isclose(error**2, covariance[index, index], rtol=1e-8, atol=1e-8):
                raise ValueError("reported cross-shear error disagrees with covariance diagonal")
            cross_indices.append(index)
            cross_values.append(value)
            cross_errors.append(error)
    cross_indices = np.asarray(cross_indices, dtype=int)
    cross_values = np.asarray(cross_values, dtype=float)
    cross_errors = np.asarray(cross_errors, dtype=float)
    cross_covariance = covariance[np.ix_(cross_indices, cross_indices)]

    global_test = _joint_test(cross_values, cross_covariance)
    prior_test = summary.get("joint_cross_null", {})
    for field in ("chi2", "degrees_of_freedom", "p_value"):
        if field not in prior_test or not np.isclose(global_test[field], prior_test[field], rtol=1e-9, atol=1e-12):
            raise ValueError(f"recomputed global cross null disagrees with input: {field}")

    radial_tests = []
    for radial_index, (r_min, r_max) in enumerate(radial_keys):
        local = np.asarray([source_bin * n_radial + radial_index for source_bin in range(n_bin)])
        test = _joint_test(cross_values[local], cross_covariance[np.ix_(local, local)])
        radial_tests.append({"r_min_mpc": r_min, "r_max_mpc": r_max, **test})
    radial_adjusted = holm([row["p_value"] for row in radial_tests])
    for row, adjusted in zip(radial_tests, radial_adjusted):
        row["p_holm_radial_bins"] = float(adjusted)

    source_tests = []
    for bin_index, source_bin in enumerate(bins):
        local = np.arange(bin_index * n_radial, (bin_index + 1) * n_radial)
        test = _joint_test(cross_values[local], cross_covariance[np.ix_(local, local)])
        source_tests.append({"source_zbin": source_bin, **test})
    source_adjusted = holm([row["p_value"] for row in source_tests])
    for row, adjusted in zip(source_tests, source_adjusted):
        row["p_holm_four_source_bins"] = float(adjusted)

    marginal_p = 2.0 * norm.sf(np.abs(cross_values / cross_errors))
    marginal_adjusted = holm(marginal_p)
    marginals = []
    for position, (value, error, p_value, adjusted) in enumerate(
            zip(cross_values, cross_errors, marginal_p, marginal_adjusted)):
        source_bin = bins[position // n_radial]
        radial_index = position % n_radial
        r_min, r_max = radial_keys[radial_index]
        marginals.append({
            "source_zbin": source_bin, "r_min_mpc": r_min, "r_max_mpc": r_max,
            "cross_delta_sigma_msun_mpc2": float(value),
            "jackknife_error_msun_mpc2": float(error),
            "standardized_residual": float(value / error),
            "p_value_normal_marginal": float(p_value),
            "p_holm_28_marginals": float(adjusted),
        })
    return {
        "schema": "astrafek-wl-cross-shear-null-decomposition-v1",
        "status": "validated_cluster_only_null_diagnostic",
        "science_measurement": False,
        "random_lens_subtraction_applied": False,
        "input_summary_sha256": source_sha256,
        "n_source_bins": n_bin,
        "n_radial_bins": n_radial,
        "n_cross_components": len(cross_values),
        "n_occupied_cells": int(summary["n_occupied_cells"]),
        "global_cross_null": global_test,
        "radial_bin_joint_tests": radial_tests,
        "source_bin_joint_tests": source_tests,
        "marginal_tests": marginals,
        "method_notes": [
            "Radial and source-bin tests use the corresponding marginal covariance submatrix.",
            "Marginal component p-values use a Gaussian approximation and are Holm-adjusted across 28 components.",
            "These are exploratory localization diagnostics; they do not identify a physical cause or replace random-lens subtraction.",
            "The saved jackknife covariance uses 5-degree equal-area cells and is not an observed-footprint covariance.",
        ],
        "limitations": [
            "cluster-only stack; no random-lens subtraction",
            "fiducial n(z), shape-calibration, and source-selection systematics remain",
            "equal-area cells can combine disconnected survey regions",
            "science_measurement remains false",
        ],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    input_path = args.input.resolve()
    output_dir = args.output_dir.resolve()
    if output_dir.exists():
        raise FileExistsError(f"refusing to overwrite output directory: {output_dir}")
    summary = json.loads(input_path.read_text(encoding="utf-8"))
    result = audit(summary, sha256(input_path))
    output_dir.mkdir(parents=True, exist_ok=False)
    result["analysis_script_sha256"] = sha256(Path(__file__).resolve())
    output_path = output_dir / "cross_shear_null_decomposition.json"
    with output_path.open("x", encoding="utf-8") as stream:
        json.dump(result, stream, ensure_ascii=False, indent=2, allow_nan=False)
        stream.write("\n")
    result["output_sha256"] = sha256(output_path)
    print(json.dumps({"status": result["status"], "science_measurement": False,
                      "global_cross_null": result["global_cross_null"],
                      "radial_bin_joint_tests": result["radial_bin_joint_tests"],
                      "source_bin_joint_tests": result["source_bin_joint_tests"],
                      "output": str(output_path), "output_sha256": result["output_sha256"]},
                     ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
