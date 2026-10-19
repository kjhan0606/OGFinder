#!/usr/bin/env python3
"""Test whether source-bin offsets differ between lens-redshift halves.

This post-processes the saved Astrafek lens-redshift tomography summary. At
each radius it tests three difference-in-differences contrasts, allowing each
lens-redshift cohort to have its own common lensing amplitude. It uses the
saved joint spatial-jackknife covariance and requires no new catalog data.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.stats import chi2

EXPECTED_SCHEMA = "astrafek-hsc-y3-lens-redshift-tomographic-jackknife-v1"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def holm_adjust(p_values):
    values = np.asarray(p_values, dtype=float)
    if values.ndim != 1 or not np.all(np.isfinite(values)) or np.any((values < 0) | (values > 1)):
        raise ValueError("p-values must be a finite vector in [0, 1]")
    order = np.argsort(values)
    adjusted = np.empty_like(values)
    running = 0.0
    for rank, index in enumerate(order):
        running = max(running, (len(values) - rank) * values[index])
        adjusted[index] = min(1.0, running)
    return adjusted


def analyze_interaction(summary):
    if (summary.get("schema") != EXPECTED_SCHEMA or summary.get("status") != "exploratory"
            or summary.get("science_measurement") is not False):
        raise ValueError("input is not a completed exploratory lens-redshift tomography summary")
    edges = np.asarray(summary.get("radial_edges_mpc", []), dtype=float)
    if len(edges) < 2 or not np.all(np.isfinite(edges)) or np.any(np.diff(edges) <= 0):
        raise ValueError("invalid radial edges")
    n_radial = len(edges) - 1
    n_features = 2 * 4 * 2 * n_radial
    if summary.get("profile_order") != (
        "lens_redshift_cohort lower_z then upper_z; source_zbin 1..4; radial bin; tangential and cross"
    ):
        raise ValueError("unexpected profile ordering")
    covariance = np.asarray(summary.get("covariance"), dtype=float)
    if covariance.shape != (n_features, n_features) or not np.all(np.isfinite(covariance)):
        raise ValueError("covariance shape or values do not match profile ordering")
    covariance = 0.5 * (covariance + covariance.T)
    profile = summary.get("profile")
    if not isinstance(profile, list) or len(profile) != 2 * 4 * n_radial:
        raise ValueError("profile does not contain both cohorts, four source bins, and every radius")
    # The full covariance includes both tangential and cross components; this
    # interaction uses only tangential values, so cross slots are zero-filled
    # and receive no contrast weight below.
    values = np.zeros(n_features, dtype=float)
    cohort_index = {"lower_z": 0, "upper_z": 1}
    seen = set()
    for row in profile:
        cohort = row.get("lens_redshift_cohort")
        source_bin = int(row.get("source_zbin", 0))
        r_min, r_max = float(row["r_min_mpc"]), float(row["r_max_mpc"])
        radial_matches = np.flatnonzero(np.isclose(edges[:-1], r_min) & np.isclose(edges[1:], r_max))
        radial = int(radial_matches[0]) if len(radial_matches) == 1 else -1
        key = (cohort, source_bin, radial)
        if cohort not in cohort_index or source_bin not in (1, 2, 3, 4) or not (0 <= radial < n_radial):
            raise ValueError(f"invalid profile key: {key}")
        if key in seen:
            raise ValueError(f"duplicate profile key: {key}")
        seen.add(key)
        index = cohort_index[cohort] * 4 * 2 * n_radial + (source_bin - 1) * 2 * n_radial + 2 * radial
        values[index] = float(row["delta_sigma_weak_msun_mpc2"])
    if len(seen) != 2 * 4 * n_radial or not np.all(np.isfinite(values)):
        raise ValueError("profile keys are incomplete or profile values are non-finite")

    tests = []
    all_contrast_rows = []
    for radial in range(n_radial):
        matrix = np.zeros((3, n_features), dtype=float)
        for k, source_bin in enumerate((2, 3, 4)):
            low_bin1 = (0 * 4 + 0) * 2 * n_radial + 2 * radial
            low_bin_k = (0 * 4 + source_bin - 1) * 2 * n_radial + 2 * radial
            high_bin1 = (1 * 4 + 0) * 2 * n_radial + 2 * radial
            high_bin_k = (1 * 4 + source_bin - 1) * 2 * n_radial + 2 * radial
            matrix[k, low_bin_k] = 1.0
            matrix[k, low_bin1] = -1.0
            matrix[k, high_bin_k] = -1.0
            matrix[k, high_bin1] = 1.0
            all_contrast_rows.append(matrix[k].copy())
        contrast = matrix @ values
        contrast_cov = matrix @ covariance @ matrix.T
        contrast_cov = 0.5 * (contrast_cov + contrast_cov.T)
        rank = int(np.linalg.matrix_rank(contrast_cov))
        if rank == 0:
            statistic = 0.0 if np.allclose(contrast, 0.0) else float("inf")
            p_value = 1.0 if statistic == 0.0 else 0.0
        else:
            statistic = float(contrast @ np.linalg.pinv(contrast_cov, hermitian=True) @ contrast)
            p_value = float(chi2.sf(statistic, rank))
        tests.append({
            "radial_bin": radial,
            "r_min_mpc": float(edges[radial]),
            "r_max_mpc": float(edges[radial + 1]),
            "contrast_order": ["(bin 2 - bin 1)_lower - (bin 2 - bin 1)_upper",
                               "(bin 3 - bin 1)_lower - (bin 3 - bin 1)_upper",
                               "(bin 4 - bin 1)_lower - (bin 4 - bin 1)_upper"],
            "difference_in_differences_msun_mpc2": contrast.tolist(),
            "chi2": statistic,
            "degrees_of_freedom": rank,
            "p_value": p_value,
            "covariance_rank": rank,
        })
    adjusted = holm_adjust([test["p_value"] for test in tests])
    for test, p_adjusted in zip(tests, adjusted):
        test["holm_p_value_7_radial_tests"] = float(p_adjusted)

    global_matrix = np.asarray(all_contrast_rows, dtype=float)
    global_contrast = global_matrix @ values
    global_covariance = global_matrix @ covariance @ global_matrix.T
    global_covariance = 0.5 * (global_covariance + global_covariance.T)
    global_rank = int(np.linalg.matrix_rank(global_covariance))
    if global_rank == 0:
        global_statistic = 0.0 if np.allclose(global_contrast, 0.0) else float("inf")
        global_p = 1.0 if global_statistic == 0.0 else 0.0
    else:
        global_statistic = float(global_contrast @ np.linalg.pinv(global_covariance, hermitian=True) @ global_contrast)
        global_p = float(chi2.sf(global_statistic, global_rank))

    return {
        "schema": "astrafek-hsc-y3-lens-redshift-tomographic-interaction-v1",
        "status": "exploratory",
        "science_measurement": False,
        "null_hypothesis": "relative source-bin DeltaSigma offsets (bins 2-4 versus bin 1) are the same in both lens-redshift cohorts; each cohort has its own common amplitude",
        "method": "difference-in-differences Wald tests using the saved full joint sky-cell jackknife covariance; one three-contrast test per radius and one global test across radii",
        "radial_edges_mpc": edges.tolist(),
        "n_input_lenses": summary["lens_redshift_split"]["n_lenses_total"],
        "lens_redshift_split": summary["lens_redshift_split"],
        "per_radial_bin_interaction": tests,
        "global_interaction": {
            "chi2": global_statistic,
            "degrees_of_freedom": global_rank,
            "p_value": global_p,
            "covariance_rank": global_rank,
        },
        "limitations": [
            "This tests changes in relative source-bin offsets, not equality of the physical lensing amplitude between redshift cohorts.",
            "The result inherits the cluster-only stack, fiducial n(z), equal-area jackknife, and unpropagated redshift-calibration uncertainty.",
            "The Wald tests use a Gaussian covariance approximation and are diagnostics, not lensing detections.",
        ],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        parser.error(f"refusing to overwrite existing output: {args.output}")
    source = json.loads(args.summary.read_text(encoding="utf-8"))
    result = analyze_interaction(source)
    result["input"] = {
        "summary_path": str(args.summary.resolve()),
        "summary_sha256": sha256(args.summary),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("x", encoding="utf-8") as stream:
        json.dump(result, stream, indent=2, ensure_ascii=False, allow_nan=False)
        stream.write("\n")
    print(json.dumps({
        "status": result["status"],
        "global_interaction": result["global_interaction"],
        "output": str(args.output.resolve()),
    }, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
