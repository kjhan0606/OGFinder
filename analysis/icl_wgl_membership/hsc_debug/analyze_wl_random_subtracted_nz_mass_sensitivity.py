#!/usr/bin/env python3
"""Re-stack and refit the random-subtracted WL signal under HSC Y3 n(z) shifts.

This offline sensitivity reuses only the previously approved source shape rows,
matched random-center rows, and local SACC n(z). Each scenario re-estimates the
cluster and random DeltaSigma profiles using the same 435-lens common-support
sample, then refits the reduced-shear NFW model. The saved baseline spatial
jackknife covariance is held fixed to isolate the effect of the n(z) shift;
scenario-specific covariance and calibrated spline posterior uncertainty are
not propagated. The science gate remains false.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys
import time

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
ROOT = HERE.parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from audit_wl_published_nz_shift_sensitivity import (
    PUBLISHED_SHIFT_SCENARIOS,
    translate_nz_curve,
)
from fit_wl_random_stack_nfw_mass import (
    DATA, FIT_PROFILE_BINS, LENS_PLAN, NZ_SACC, RANDOM_SUMMARY,
    TOMO_SUMMARY, load_pairs_and_provenance,
)
from fit_wl_random_stack_reduced_shear_sensitivity import (
    fit_profile as fit_reduced_shear_profile,
    prepare_rows as prepare_reduced_shear_rows,
)
from fetch_ancillary_catalogues import sha
from ogfkit.wlcluster import Cosmology
from run_wl_random_shape_stack import verify_and_load_pairs as verify_random_pairs
from s19a_shape_stack import (
    load_hsc_y3_nz_sacc,
    stack_s19a_tomographic_pairs,
    subtract_random_lens_profile,
)

RANDOM_ROW_PLAN = DATA / "wl_random_shape_rows_lens_reuse_2026-11-03_01"
RANDOM_ROWS = RANDOM_ROW_PLAN / "results"
RANDOM_COUNT_RESULTS = DATA / "wl_random_shape_count_results_lens_reuse_2026-11-03_01"
DEFAULT_OUTPUT = DATA / "wl_random_subtracted_nz_mass_sensitivity_2026-11-06_01"
RANDOM_EXPORT_AUDIT = DATA / "wl_random_shape_export_integrity_2026-11-05_01/summary.json"
BASE_NFW = DATA / "wl_random_stack_nfw_mass_2026-11-05_01/summary.json"
EXPECTED_RANDOM_STACK_SHA = "83c5484454033b2652adcba5bd0699f28082fa6983eefa960613da75491bafc8"
EXPECTED_TOMO_SHA = "47419ac229ffec9aab2a001f629252627126a456812b01cf11908ad44122fe14"
EXPECTED_BASE_NFW_SHA = "d3ca78a24fc030f64694247fd45d33f607d89c18fb102b6b5a6879d53c52cd4b"
EXPECTED_RANDOM_EXPORT_AUDIT_SHA = "b7a371a6ff5e4c8d546cbc88ec81a2289f57c09816fd5874eef135e2cc7da27d"
RADIAL_EDGES_MPC = np.asarray([0.1, 0.2, 0.4, 0.7, 1.0, 1.5, 2.0, 3.0], dtype=float)
TANGENTIAL_INDEX = np.arange(0, 14, 2)
CONCENTRATIONS = (3.0, 4.0, 5.0)
MASKED_B_MODE_VALUE = 0
FIDUCIAL_SHIFT = {1: 0.0, 2: 0.0, 3: 0.0, 4: 0.0}
ALL_SCENARIOS = {"fiducial": FIDUCIAL_SHIFT, **PUBLISHED_SHIFT_SCENARIOS}


def selected_bins(pairs: dict) -> dict:
    bins = np.asarray(pairs["source_zbin"], dtype=int)
    keep = (bins >= 2) & (bins <= 4)
    if not np.any(keep):
        raise ValueError("source-bin 1 exclusion removed every pair")
    return {name: np.asarray(values)[keep] for name, values in pairs.items()}


def profile_values(result: dict) -> np.ndarray:
    rows = result.get("profile", [])
    if len(rows) != len(RADIAL_EDGES_MPC) - 1:
        raise ValueError("stack result does not contain all seven radial bins")
    values = np.asarray([
        [row["delta_sigma_weak_msun_mpc2"],
         row["cross_delta_sigma_weak_msun_mpc2"]]
        for row in rows
    ], dtype=float)
    if values.shape != (7, 2) or not np.all(np.isfinite(values)):
        raise ValueError("stack result has missing or non-finite profile bins")
    return values


def translate_all(z_grid: np.ndarray, base_nz: np.ndarray,
                  shifts: dict[int, float]) -> tuple[np.ndarray, dict]:
    shifted = np.asarray(base_nz, dtype=float).copy()
    details = {}
    if shifted.shape != (4, len(z_grid)):
        raise ValueError("expected the four official HSC Y3 source-bin curves")
    for source_bin in range(1, 5):
        shift = float(shifts[source_bin])
        if shift:
            curve, detail = translate_nz_curve(z_grid, shifted[source_bin - 1], shift)
            shifted[source_bin - 1] = curve
        else:
            curve = shifted[source_bin - 1]
            norm = float(np.trapezoid(curve, z_grid))
            mean = float(np.trapezoid(z_grid * curve, z_grid) / norm)
            detail = {
                "requested_mean_shift": 0.0,
                "mean_redshift_before": mean,
                "mean_redshift_after": mean,
                "realized_mean_shift_after_grid_truncation": 0.0,
                "renormalized_integral": 1.0,
            }
        details[str(source_bin)] = detail
    return shifted, details


def run(output: Path, *, verbose: bool = True) -> dict:
    output = Path(output).resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite output: {output}")
    for path, expected in (
        (RANDOM_SUMMARY, EXPECTED_RANDOM_STACK_SHA),
        (TOMO_SUMMARY, EXPECTED_TOMO_SHA),
        (BASE_NFW, EXPECTED_BASE_NFW_SHA),
        (RANDOM_EXPORT_AUDIT, EXPECTED_RANDOM_EXPORT_AUDIT_SHA),
    ):
        if sha(path) != expected:
            raise ValueError(f"pinned input changed: {path}")

    start = time.monotonic()
    target_pairs, random_doc, _tomo = load_pairs_and_provenance()
    random_pairs, random_provenance = verify_random_pairs(
        RANDOM_ROWS, RANDOM_ROW_PLAN, MASKED_B_MODE_VALUE,
    )
    common_lenses = {str(value) for value in np.unique(target_pairs["lens_id"])}
    center_to_lens = random_provenance["random_center_matched_lens_ids"]
    selected_center_ids = {
        center_id for center_id, lens_id in center_to_lens.items()
        if str(lens_id) in common_lenses
    }
    if len(common_lenses) != 435 or len(selected_center_ids) != 2175:
        raise ValueError("verified common sample is not 435 lenses with five random centers each")
    random_keep = np.isin(random_pairs["lens_id"].astype(str), list(selected_center_ids))
    random_pairs = {name: np.asarray(values)[random_keep]
                    for name, values in random_pairs.items()}
    if not set(np.unique(random_pairs["lens_id"].astype(str))).issubset(selected_center_ids):
        raise ValueError("random-pair filtering retained an unrequested center")
    target_pairs = selected_bins(target_pairs)
    random_pairs = selected_bins(random_pairs)

    z_grid, baseline_nz = load_hsc_y3_nz_sacc(NZ_SACC)
    cosmology = Cosmology(H0=70.0, Om=0.3)
    tomo = json.loads(TOMO_SUMMARY.read_text(encoding="utf-8"))
    base_tomo = tomo["block_width_results"]["5.0"]["source_zbin_1_exclusion_sensitivity"]
    baseline_covariance = np.asarray(base_tomo["profile_covariance_msun2_mpc4"], dtype=float)
    if baseline_covariance.shape != (14, 14) or not np.all(np.isfinite(baseline_covariance)):
        raise ValueError("saved baseline source-bin-excluded covariance is invalid")
    tangential_covariance = baseline_covariance[np.ix_(TANGENTIAL_INDEX, TANGENTIAL_INDEX)]
    if np.min(np.linalg.eigvalsh(tangential_covariance)) <= 0:
        raise ValueError("baseline tangential covariance is not positive definite")
    baseline_expected = np.asarray(base_tomo["profile_msun_mpc2"], dtype=float)[:, 0]

    all_profiles, all_fits, all_shift_details = {}, {}, {}
    fiducial_comparison = None
    diagnostics = []
    for scenario_name, scenario_shifts in ALL_SCENARIOS.items():
        shifted_nz, details = translate_all(z_grid, baseline_nz, scenario_shifts)
        lens_stack = stack_s19a_tomographic_pairs(
            target_pairs, z_grid, shifted_nz, RADIAL_EDGES_MPC, cosmology,
        )
        random_stack = stack_s19a_tomographic_pairs(
            random_pairs, z_grid, shifted_nz, RADIAL_EDGES_MPC, cosmology,
        )
        subtracted = subtract_random_lens_profile(lens_stack, random_stack)
        profile = np.asarray([
            row["delta_sigma_random_subtracted_msun_mpc2"]
            for row in subtracted["profile"]
        ], dtype=float)
        cross = np.asarray([
            row["cross_delta_sigma_random_subtracted_msun_mpc2"]
            for row in subtracted["profile"]
        ], dtype=float)
        if profile.shape != (7,) or not np.all(np.isfinite(profile)):
            raise ValueError(f"invalid random-subtracted profile in {scenario_name}")
        if scenario_name == "fiducial":
            relative = np.max(np.abs(profile - baseline_expected) /
                              np.maximum(np.abs(baseline_expected), 1.0))
            fiducial_comparison = {
                "max_abs_difference_msun_mpc2": float(np.max(np.abs(profile - baseline_expected))),
                "max_fractional_difference": float(relative),
                "recomputed_pair_level_profile_msun_mpc2": profile.tolist(),
                "saved_tomography_block_tensor_profile_msun_mpc2": baseline_expected.tolist(),
                "interpretation": "The tomography jackknife uses block-tensor aggregation while this runner recomputes the pair-level HSC estimator; retain any difference as estimator sensitivity.",
            }

        model_rows = prepare_reduced_shear_rows(
            target_pairs, z_grid, shifted_nz, RADIAL_EDGES_MPC, cosmology,
        )
        fits = {
            str(concentration): fit_reduced_shear_profile(
                model_rows, cosmology, profile, tangential_covariance,
                concentration,
            )
            for concentration in CONCENTRATIONS
        }
        all_profiles[scenario_name] = {
            "source_bin_1_excluded_delta_sigma_msun_mpc2": profile.tolist(),
            "source_bin_1_excluded_cross_msun_mpc2": cross.tolist(),
            "change_from_central_delta_sigma_msun_mpc2": None,
            "target_lens_count": len(common_lenses),
            "random_center_count": len(selected_center_ids),
        }
        all_fits[scenario_name] = fits
        all_shift_details[scenario_name] = details
        diagnostics.append({
            "scenario": scenario_name,
            "max_abs_cross_msun_mpc2": float(np.max(np.abs(cross))),
            "min_realized_shift_bin3": details["3"]["realized_mean_shift_after_grid_truncation"],
            "min_realized_shift_bin4": details["4"]["realized_mean_shift_after_grid_truncation"],
        })
        if verbose:
            masses = {c: fits[c]["best_fit_M200c_Msun"] for c in fits}
            print(json.dumps({"scenario": scenario_name,
                              "profile_first_bin": float(profile[0]),
                              "mass_Msun_by_concentration": masses}, ensure_ascii=False),
                  flush=True)

    central = np.asarray(all_profiles["central"]["source_bin_1_excluded_delta_sigma_msun_mpc2"])
    for profile in all_profiles.values():
        shifted_profile = np.asarray(profile["source_bin_1_excluded_delta_sigma_msun_mpc2"])
        profile["change_from_central_delta_sigma_msun_mpc2"] = (shifted_profile - central).tolist()

    base_mass = json.loads(BASE_NFW.read_text(encoding="utf-8"))
    result = {
        "schema": "astrafek-wl-random-subtracted-nz-mass-sensitivity-v1",
        "status": "complete_exploratory_sensitivity",
        "science_measurement": False,
        "hsc_network_used": False,
        "images_read": False,
        "new_data_requested": False,
        "sample": {
            "common_cluster_lenses": len(common_lenses),
            "matched_random_centers": len(selected_center_ids),
            "source_bins": [2, 3, 4],
            "radial_edges_mpc": RADIAL_EDGES_MPC.tolist(),
        },
        "method": "recompute HSC Y3 tomographic cluster and matched-random DeltaSigma stacks for each constant n(z) translation scenario, subtract the random profile, then refit exact-pair reduced-shear NFW models at fixed concentration",
        "scenarios": {
            name: {
                "requested_mean_redshift_shifts": {
                    str(index): float(value) for index, value in shifts.items()
                },
                "realized_shift_details": all_shift_details[name],
                "profile": all_profiles[name],
                "reduced_shear_nfw_fits_by_concentration": all_fits[name],
            }
            for name, shifts in ALL_SCENARIOS.items()
        },
        "mass_shift_vs_fiducial_fraction": {
            concentration: {
                name: float(
                    all_fits[name][concentration]["best_fit_M200c_Msun"]
                    / all_fits["central"][concentration]["best_fit_M200c_Msun"] - 1.0
                )
                for name in ALL_SCENARIOS
            }
            for concentration in map(str, CONCENTRATIONS)
        },
        "baseline_mass_screen_by_concentration": base_mass["by_jackknife_block_width_deg"]["5.0"],
        "uncertainty_treatment": "saved 5-degree source-bin-1-excluded tangential spatial-jackknife covariance held fixed across n(z) scenarios; scenario-specific covariance not recomputed",
        "fiducial_estimator_comparison": fiducial_comparison,
        "limitations": [
            "Published n(z) mean-shift endpoints are applied as constant translations of only source bins 3 and 4, not the authors' calibrated spline curves or joint posterior samples.",
            "The five scenarios form a rectangular grid of marginal 68% endpoints and do not encode covariance between bin shifts.",
            "The saved baseline covariance is held fixed; changed estimator weights can alter the covariance and that uncertainty is not propagated.",
            "NFW masses are effective single-halo stack scales with fixed concentration, not mean cluster masses or a mass calibration.",
            "Miscentering, magnification, two-halo structure, calibrated random-aperture matching, and field-covariance validation remain incomplete.",
        ],
        "input_sha256": {
            str(RANDOM_SUMMARY): sha(RANDOM_SUMMARY),
            str(TOMO_SUMMARY): sha(TOMO_SUMMARY),
            str(BASE_NFW): sha(BASE_NFW),
            str(RANDOM_EXPORT_AUDIT): sha(RANDOM_EXPORT_AUDIT),
            str(NZ_SACC): sha(NZ_SACC),
            str(LENS_PLAN / "row_export_plan.json"): sha(LENS_PLAN / "row_export_plan.json"),
            str(RANDOM_ROW_PLAN / "random_shape_row_plan.json"): sha(RANDOM_ROW_PLAN / "random_shape_row_plan.json"),
        },
        "elapsed_seconds": time.monotonic() - start,
        "diagnostics": diagnostics,
    }
    output.mkdir(parents=True, exist_ok=False)
    summary_path = output / "summary.json"
    summary_path.write_text(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False) + "\n",
                             encoding="utf-8")
    result["summary_sha256"] = sha(summary_path)
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--summary-only", type=Path,
                        help="print compact results from an existing machine summary")
    args = parser.parse_args()
    if args.summary_only:
        result = json.loads(args.summary_only.read_text(encoding="utf-8"))
    else:
        result = run(args.output)
    compact = {
        name: {
            concentration: {
                "mass_Msun": scenario["reduced_shear_nfw_fits_by_concentration"][concentration]["best_fit_M200c_Msun"],
                "conditional_log10_interval": scenario["reduced_shear_nfw_fits_by_concentration"][concentration]["conditional_profile_interval_log10_M200c"],
                "chi2_dof": [scenario["reduced_shear_nfw_fits_by_concentration"][concentration]["chi2"], scenario["reduced_shear_nfw_fits_by_concentration"][concentration]["dof"]],
            }
            for concentration in map(str, CONCENTRATIONS)
        }
        for name, scenario in result["scenarios"].items()
    }
    summary_hash = result.get("summary_sha256")
    if summary_hash is None and args.summary_only:
        summary_hash = sha(args.summary_only)
    print(json.dumps({"status": result["status"], "science_measurement": False,
                      "mass_fits": compact, "mass_shift_vs_fiducial_fraction": result["mass_shift_vs_fiducial_fraction"],
                      "summary_sha256": summary_hash}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
