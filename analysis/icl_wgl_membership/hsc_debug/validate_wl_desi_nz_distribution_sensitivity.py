#!/usr/bin/env python3
"""Independently validate the saved DESI/HSC n(z) jackknife comparison."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np


def sha(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def validate(summary_path: Path, baseline_path: Path, desi_directory: Path) -> dict:
    summary_path = Path(summary_path).resolve()
    baseline_path = Path(baseline_path).resolve()
    desi_directory = Path(desi_directory).resolve()
    result = json.loads(summary_path.read_text(encoding="utf-8"))
    baseline = json.loads(baseline_path.read_text(encoding="utf-8"))
    manifest_path = desi_directory / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if (result.get("schema") != "astrafek-hsc-pdr3-desi-nz-spatial-jackknife-sensitivity-v1"
            or result.get("status") != "complete_cluster_only_nz_sensitivity"
            or result.get("science_measurement") is not False
            or result.get("random_lens_subtraction_applied") is not False):
        raise ValueError("DESI result has an unexpected schema or scientific status")
    if (baseline.get("schema") != "astrafek-hsc-y3-tomographic-consistency-jackknife-v1"
            or baseline.get("science_measurement") is not False):
        raise ValueError("saved SACC baseline has an unexpected schema or status")
    if manifest.get("commit") != "abc70af29910e49fce17a63ef3083160d515fb6c":
        raise ValueError("calibration manifest is not pinned to the approved revision")

    inp, base_inp = result["input"], baseline["input"]
    same_fields = (
        "row_plan_sha256", "row_batch_sha256", "n_pair_rows", "n_lenses",
        "masked_value", "n_z_sacc_sha256", "radial_edges_mpc",
    )
    for key in same_fields:
        if inp.get(key) != base_inp.get(key):
            raise ValueError(f"DESI and SACC baseline differ in input provenance: {key}")
    desi_cosmo = inp.get("cosmology", {})
    sacc_cosmo = base_inp.get("cosmology", {})
    if (desi_cosmo.get("H0_km_s_Mpc") != sacc_cosmo.get("H0_km_s_mpc")
            or desi_cosmo.get("Omega_m") != sacc_cosmo.get("Omega_m")):
        raise ValueError("DESI and SACC baseline cosmologies differ")
    if (inp.get("baseline_summary_sha256") != sha(baseline_path)
            or inp.get("baseline_summary_path") != str(baseline_path)):
        raise ValueError("stored SACC baseline path/hash does not verify")
    if (inp.get("n_pair_rows") != 3069544 or inp.get("n_lenses") != 800
            or len(inp.get("row_batch_sha256", {})) != 32
            or result.get("photoz_posterior_samples_available") is not False):
        raise ValueError("shape-row provenance or DESI posterior-sample state is unexpected")
    if inp.get("desi_calibration_directory") != str(desi_directory):
        raise ValueError("result points to a different DESI calibration directory")

    manifest_hashes = {}
    for row in manifest.get("files", []):
        path = desi_directory / row["path"]
        if (not path.is_file() or path.stat().st_size != row["bytes"]
                or sha(path) != row["sha256"]):
            raise ValueError(f"approved DESI artifact hash/size mismatch: {row['path']}")
        manifest_hashes[str(path)] = row["sha256"]
    manifest_hashes[str(manifest_path)] = sha(manifest_path)
    if manifest_hashes != inp.get("desi_artifact_sha256") or len(manifest.get("files", [])) != 4:
        raise ValueError("DESI output artifact hashes do not match the four-file manifest")

    expected_names = {
        f"{variant}__{kind}"
        for variant in ("0.3_3_v_1p1", "1_5_v_1p1")
        for kind in ("mean", "mean_mu", "lower", "upper", "lower_mu", "upper_mu")
    }
    scenarios = result.get("scenarios", {})
    if set(scenarios) != expected_names:
        raise ValueError("saved result does not contain the complete 12-scenario grid")
    if baseline.get("n_occupied_cells") != 52 or baseline.get("n_input_pairs_in_tomographic_bins") != 2254646:
        raise ValueError("SACC baseline does not have the expected 52 cells and assigned-pair count")
    base_profile = {
        (int(row["source_zbin"]), float(row["r_min_mpc"]), float(row["r_max_mpc"])): row
        for row in baseline["profile"]
    }
    if len(base_profile) != 28:
        raise ValueError("SACC baseline profile does not contain the 4x7 expected bins")

    central_checks, envelope_checks = {}, {}
    for name, scenario in scenarios.items():
        jk = scenario["spatial_jackknife"]
        rows = jk.get("profile", [])
        if len(rows) != 28:
            raise ValueError(f"{name} has the wrong profile length")
        keys = {(int(r["source_zbin"]), float(r["r_min_mpc"]), float(r["r_max_mpc"])) for r in rows}
        if keys != set(base_profile):
            raise ValueError(f"{name} profile keys differ from SACC baseline")
        comparisons = {
            (int(r["source_zbin"]), float(r["r_min_mpc"]), float(r["r_max_mpc"])): r
            for r in scenario["vs_sacc_baseline"]
        }
        if set(comparisons) != set(base_profile):
            raise ValueError(f"{name} comparison profile keys are incomplete")
        for key, row in ((k, next(v for v in rows if (int(v["source_zbin"]), float(v["r_min_mpc"]), float(v["r_max_mpc"])) == k)) for k in keys):
            comp, base = comparisons[key], base_profile[key]
            observed = float(row["delta_sigma_weak_msun_mpc2"])
            expected = float(base["delta_sigma_weak_msun_mpc2"])
            if (not np.isclose(comp["desi_delta_sigma_msun_mpc2"], observed, rtol=1e-12, atol=1e-6)
                    or not np.isclose(comp["sacc_delta_sigma_msun_mpc2"], expected, rtol=1e-12, atol=1e-6)
                    or not np.isclose(comp["change_msun_mpc2"], observed - expected, rtol=1e-10, atol=1e-6)):
                raise ValueError(f"{name} saved comparison disagrees with the profile")
        fractions = np.array([r["fractional_change"] for r in scenario["vs_sacc_baseline"] if r["fractional_change"] is not None], dtype=float)
        if len(fractions) != 28 or not np.all(np.isfinite(fractions)):
            raise ValueError(f"{name} fractional profile changes are incomplete")

        if scenario["is_central_estimate"]:
            covariance = np.asarray(jk.get("covariance"), dtype=float)
            if (jk.get("n_occupied_cells") != 52 or jk.get("n_input_pairs_in_tomographic_bins") != 2254646
                    or covariance.shape != (56, 56) or not np.all(np.isfinite(covariance))):
                raise ValueError(f"{name} central covariance dimensions/counts are invalid")
            scale = max(float(np.max(np.abs(covariance))), 1.0)
            relative_asymmetry = float(np.max(np.abs(covariance - covariance.T)) / scale)
            rank = int(np.linalg.matrix_rank(covariance))
            if relative_asymmetry > 1e-12 or rank != 51:
                raise ValueError(f"{name} covariance symmetry or rank check failed")
            if jk.get("source_bin_counts") != baseline.get("source_bin_counts"):
                raise ValueError(f"{name} source-bin row counts differ from the SACC baseline")
            central_checks[name] = {
                "n_occupied_cells": 52,
                "n_assigned_pairs": 2254646,
                "covariance_shape": list(covariance.shape),
                "covariance_rank": rank,
                "relative_covariance_asymmetry": relative_asymmetry,
                "fractional_change_median_by_source_bin": {
                    str(bin_index): float(np.median([
                        comp["fractional_change"] for key, comp in comparisons.items() if key[0] == bin_index
                    ])) for bin_index in range(1, 5)
                },
                "consistency_p_0p1_0p7_mpc": [
                    {"r_min_mpc": row["r_min_mpc"], "r_max_mpc": row["r_max_mpc"], "p_value": row["p_value"]}
                    for row in jk["tomographic_delta_sigma_consistency"] if row["r_max_mpc"] <= 0.7
                ],
            }
        else:
            if scenario.get("uncertainty_method", "").startswith("none") is False:
                raise ValueError(f"{name} pointwise envelope unexpectedly claims a covariance")
            if any(row.get("desi_spatial_jackknife_error_msun_mpc2") is not None for row in scenario["vs_sacc_baseline"]):
                raise ValueError(f"{name} pointwise envelope unexpectedly has jackknife errors")
            envelope_checks[name] = {
                "n_profile_bins": len(rows),
                "fractional_change_min": float(np.min(fractions)),
                "fractional_change_max": float(np.max(fractions)),
            }

    audit = {
        "schema": "astrafek-hsc-pdr3-desi-nz-sensitivity-validation-v1",
        "status": "validated_cluster_only_sensitivity",
        "science_measurement": False,
        "source_summary": str(summary_path),
        "source_summary_sha256": sha(summary_path),
        "baseline_summary_sha256": sha(baseline_path),
        "approved_desi_commit": manifest["commit"],
        "approved_desi_artifacts": 4,
        "shape_pair_rows": 3069544,
        "tomographic_pair_rows": 2254646,
        "spatial_cells": 52,
        "n_scenarios": len(scenarios),
        "central_scenario_checks": central_checks,
        "pointwise_envelope_checks": envelope_checks,
        "limitations": [
            "Validation confirms hashes, provenance, profile arithmetic, covariance symmetry/rank, and row/block counts only.",
            "The profiles remain cluster-only; no random-lens subtraction is included.",
            "Pointwise lower/upper curves are sensitivity envelopes and not joint posterior samples or covariances.",
        ],
    }
    return audit


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--baseline-summary", type=Path, required=True)
    parser.add_argument("--desi-directory", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(f"refusing to overwrite {args.output}")
    audit = validate(args.summary, args.baseline_summary, args.desi_directory)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(audit, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": audit["status"], "n_scenarios": audit["n_scenarios"],
        "spatial_cells": audit["spatial_cells"], "output": str(args.output.resolve()),
        "sha256": sha(args.output),
    }, ensure_ascii=False))


if __name__ == "__main__":
    main()
