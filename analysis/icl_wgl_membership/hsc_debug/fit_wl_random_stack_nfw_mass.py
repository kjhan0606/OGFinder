#!/usr/bin/env python3
"""Fit an exploratory single-halo NFW scale to the random-subtracted HSC stack.

The forward model averages NFW DeltaSigma over the exact matched lens-source
pairs and the same weak-shear estimator weights. It uses source bins 2--4,
excludes the innermost and outermost radial bins, and fits fixed-concentration
models against the saved paired spatial-jackknife covariance. This is a
conditional diagnostic, not a calibrated science mass measurement.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import sys

import numpy as np
from scipy.optimize import minimize_scalar
from scipy.optimize import brentq
from scipy.stats import chi2

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
ROOT = HERE.parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from ogfkit.wlcluster import Cosmology, nfw_sigma
from run_wl_random_shape_stack import (
    validate_cluster_stack_compatibility,
)
from run_wl_shape_stack import verify_and_load_pairs as verify_and_load_lens_pairs
from s19a_shape_stack import (
    C_KMS, G_MPC, hsc_selection_bias, load_hsc_y3_nz_sacc,
    pair_geometry, stack_s19a_tomographic_pairs,
    tomographic_nz_lensing_efficiency,
)
from fetch_ancillary_catalogues import sha


DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
RANDOM_SUMMARY = DATA / "wl_random_shape_stack_lens_reuse_2026-11-04_01/summary.json"
RANDOM_EXPORT_AUDIT = DATA / "wl_random_shape_export_integrity_2026-11-05_01/summary.json"
TOMO_SUMMARY = DATA / "wl_tomographic_weight_boost_2026-11-05_04/summary.json"
ROWS = DATA / "wl_shape_rows_export_2026-10-18_03"
PREFLIGHT = DATA / "wl_shape_preflight_2026-10-18_01"
LENS_PLAN = DATA / "wl_shape_rows_plan_2026-10-18_04"
NZ_SACC = DATA / "hsc_y3_nz_sacc_2026-10-19_01/hsc_y3_fourier_space_data_vector.sacc"
DEFAULT_OUTPUT = DATA / "wl_random_stack_nfw_mass_2026-11-05_01"
EXPECTED_RANDOM_ROW_PLAN_SHA = "fc883444d2b23663d5f6034d773b76db948add29f1cda0c23384d17cfed2996a"
EXPECTED_RAND_EXPORT_AUDIT_SHA = "b7a371a6ff5e4c8d546cbc88ec81a2289f57c09816fd5874eef135e2cc7da27d"
EXPECTED_TOMO_SUMMARY_SHA = "47419ac229ffec9aab2a001f629252627126a456812b01cf11908ad44122fe14"
EXPECTED_FULL_LENS_SUMMARY = DATA / "wl_cluster_stack_bmode0_2026-10-19_01/summary.json"
MASKED_B_MODE_VALUE = 0
FIT_PROFILE_BINS = (1, 2, 3, 4)  # 0.2--1.5 proper Mpc; omit inner/outer bins.
FIXED_CONCENTRATIONS = (3.0, 4.0, 5.0)
def load_pairs_and_provenance() -> tuple[dict, dict, dict]:
    random_doc = json.loads(RANDOM_SUMMARY.read_text(encoding="utf-8"))
    export_audit = json.loads(RANDOM_EXPORT_AUDIT.read_text(encoding="utf-8"))
    tomo = json.loads(TOMO_SUMMARY.read_text(encoding="utf-8"))
    full_stack = json.loads(EXPECTED_FULL_LENS_SUMMARY.read_text(encoding="utf-8"))
    if sha(RANDOM_SUMMARY) != "83c5484454033b2652adcba5bd0699f28082fa6983eefa960613da75491bafc8":
        raise ValueError("random-subtracted stack summary hash changed")
    if (export_audit.get("status") != "complete"
            or export_audit.get("science_measurement") is not False
            or export_audit.get("plan", {}).get("sha256") != EXPECTED_RANDOM_ROW_PLAN_SHA
            or len(export_audit.get("batches", [])) != 124
            or export_audit.get("summary", {}).get("n_source_pairs") != 9010324):
        raise ValueError("random source-row export is not the audited complete 124-batch product")
    if sha(RANDOM_EXPORT_AUDIT) != EXPECTED_RAND_EXPORT_AUDIT_SHA:
        raise ValueError("random row export integrity audit hash changed")
    if (tomo.get("status") != "complete_diagnostic"
            or tomo.get("science_measurement") is not False
            or tomo.get("source_bin_zero_excluded") is not True
            or sha(TOMO_SUMMARY) != EXPECTED_TOMO_SUMMARY_SHA):
        raise ValueError("tomographic source-bin exclusion diagnostic is incomplete or changed")
    if tomo.get("input_sha256", {}).get(str(RANDOM_SUMMARY)) != sha(RANDOM_SUMMARY):
        raise ValueError("tomographic diagnostic does not reference the saved random-subtracted stack")
    if (random_doc.get("status") != "exploratory"
            or random_doc.get("science_measurement") is not False
            or random_doc.get("random_lens_subtraction_applied") is not True
            or random_doc.get("spatial_block_covariance_applied") is not True
            or random_doc.get("input", {}).get("random_row_plan_sha256") != EXPECTED_RANDOM_ROW_PLAN_SHA
            or random_doc.get("input", {}).get("n_matched_lenses_in_cluster_and_random_stacks") != 435):
        raise ValueError("random-subtracted stack provenance does not match the analyzed sample")
    validate_cluster_stack_compatibility(
        full_stack, sha(NZ_SACC), MASKED_B_MODE_VALUE,
        [0.1, 0.2, 0.4, 0.7, 1.0, 1.5, 2.0, 3.0], 70.0, 0.3,
    )
    lens_pairs, lens_provenance = verify_and_load_lens_pairs(
        ROWS, PREFLIGHT, LENS_PLAN, MASKED_B_MODE_VALUE,
    )
    if lens_provenance.get("row_plan_sha256") != full_stack.get("input", {}).get("row_plan_sha256"):
        raise ValueError("verified cluster source rows differ from the saved cluster-stack plan")

    target_ids = {str(v) for v in np.unique(lens_pairs["lens_id"])}
    random_input = random_doc["input"]
    center_lens = random_input["random_center_matched_lens_ids"]
    center_counts = random_input["random_center_expected_source_counts"]
    positive_random_lenses = {
        str(center_lens[center]) for center, count in center_counts.items() if int(count) > 0
    }
    common_ids = target_ids & positive_random_lenses & {str(x) for x in random_input["matched_lens_ids"]}
    if (len(common_ids) != 435
            or hashlib.sha256("\n".join(sorted(common_ids)).encode()).hexdigest()
            != random_input["matched_sample"]["common_lens_ids_sha256"]):
        raise ValueError("reconstructed common lens sample does not match the saved random-subtracted stack")
    selected = np.isin(lens_pairs["lens_id"].astype(str), list(common_ids))
    pairs = {name: values[selected] for name, values in lens_pairs.items()}
    if {str(x) for x in np.unique(pairs["lens_id"])} != common_ids:
        raise ValueError("filtered lens source rows do not contain the exact 435-lens common sample")
    return pairs, random_doc, tomo


def prepare_estimator_weights(pairs: dict, z_grid: np.ndarray, nz_by_bin: np.ndarray,
                              edges: np.ndarray, cosmology: Cosmology) -> list[dict]:
    """Precompute exact pair geometry and estimator weights for NFW evaluation."""
    source_bin = np.asarray(pairs["source_zbin"], dtype=int)
    keep_bin = (source_bin >= 2) & (source_bin <= 4)
    if not np.any(keep_bin):
        raise ValueError("no source bins 2--4 are available for NFW modeling")
    sub = {key: np.asarray(value)[keep_bin] for key, value in pairs.items()}
    beta = tomographic_nz_lensing_efficiency(
        sub["lens_z"], sub["source_zbin"], z_grid, nz_by_bin, cosmology,
    )
    radius, _phi = pair_geometry(sub, cosmology)
    shape_amp = np.hypot(sub["e1"], sub["e2"])
    shape_ok = ((sub["sigma_e"] >= 0.0) & (sub["sigma_e"] <= 0.4)
                & (shape_amp < 2.0))
    bmode = np.asarray(sub["b_mode_pass"], dtype=bool)
    source_ok = shape_ok & (sub["mag_a10"] <= 25.5) & (sub["resolution"] >= 0.3)
    usable = bmode & source_ok & (beta > 0) & (radius >= edges[0]) & (radius <= edges[-1])
    chi_factor = 4.0 * math.pi * G_MPC / (C_KMS * C_KMS)
    d_l = np.asarray(cosmology.angular_diameter(sub["lens_z"]), dtype=float)
    sigma_crit_inv = chi_factor * d_l * beta
    usable &= sigma_crit_inv > 0
    prepared = []
    for index, (lo, hi) in enumerate(zip(edges[:-1], edges[1:])):
        in_bin = usable & (radius >= lo) & ((radius < hi) if index < len(edges) - 2 else (radius <= hi))
        if not np.any(in_bin):
            continue
        w = np.asarray(sub["weight"][in_bin], dtype=float)
        beta_k = sigma_crit_inv[in_bin]
        lensing_weight = w * beta_k * beta_k
        rms = np.asarray(sub["rms_e"][in_bin], dtype=float)
        responsivity = float(1.0 - np.sum(lensing_weight * rms * rms) / np.sum(lensing_weight))
        if not np.isfinite(responsivity) or responsivity <= 0:
            raise ValueError(f"non-positive responsivity in radial bin {index}")
        m = np.asarray(sub["m"][in_bin], dtype=float)
        selection = hsc_selection_bias(
            lensing_weight, sub["mag_a10"][in_bin], sub["resolution"][in_bin],
            sub["psf_shape11"][in_bin], sub["psf_shape22"][in_bin],
            sub["psf_shape12"][in_bin],
        )
        msel = float(selection["m_sel"])
        w_shear = w * (2.0 * responsivity * (1.0 + m) * (1.0 + msel)) ** 2
        effective_weight = w_shear * beta_k * beta_k

        prepared.append({
            "bin": index,
            "radius_mpc": np.asarray(radius[in_bin], dtype=float),
            "lens_z": np.asarray(sub["lens_z"][in_bin], dtype=float),
            "effective_weight": np.asarray(effective_weight, dtype=float),
            "weight_sum": float(np.sum(effective_weight)),
        })
    if len(prepared) != len(edges) - 1:
        raise ValueError("failed to prepare every radial-bin estimator")
    return prepared


def estimator_model_profiles(prepared: list[dict], cosmology: Cosmology,
                             mass_msun: float, concentration: float) -> np.ndarray:
    """Forward-model DeltaSigma under the exact saved weak-shear estimator weights."""
    model = np.full(len(prepared), np.nan)
    for item in prepared:
        z = item["lens_z"]
        hz = cosmology.H0 * np.sqrt(cosmology.Om * (1.0 + z) ** 3 + (1.0 - cosmology.Om))
        rho_crit = 3.0 * hz * hz / (8.0 * math.pi * G_MPC)
        r200 = (3.0 * mass_msun / (4.0 * math.pi * 200.0 * rho_crit)) ** (1.0 / 3.0)
        rs = r200 / concentration
        delta_c = (200.0 / 3.0) * concentration ** 3 / (
            math.log1p(concentration) - concentration / (1.0 + concentration))
        rho_s = delta_c * rho_crit
        _sigma, delta_sigma = nfw_sigma(item["radius_mpc"], rho_s, rs)
        model[item["bin"]] = float(
            np.sum(item["effective_weight"] * delta_sigma) / item["weight_sum"])
    return model


def profile_fit(prepared, edges, cosmology, observed, covariance,
                concentration, fit_bins=FIT_PROFILE_BINS) -> dict:
    indices = np.asarray(fit_bins, dtype=int)
    y = np.asarray(observed, dtype=float)[indices]
    cov = np.asarray(covariance, dtype=float)[np.ix_(indices, indices)]
    if cov.shape != (len(indices), len(indices)) or not np.all(np.isfinite(cov)):
        raise ValueError("invalid tangential jackknife covariance submatrix")
    if np.min(np.linalg.eigvalsh(cov)) <= 0:
        raise ValueError("tangential jackknife covariance is not positive definite")
    inv_cov = np.linalg.inv(cov)

    cache = {}
    def model(log_mass):
        key = float(log_mass)
        if key not in cache:
            cache[key] = estimator_model_profiles(
                prepared, cosmology, 10.0 ** key, concentration,
            )[indices]
        return cache[key]

    def objective(log_mass):
        delta = y - model(log_mass)
        return float(delta @ inv_cov @ delta)

    fit = minimize_scalar(objective, bounds=(12.5, 16.0), method="bounded",
                          options={"xatol": 1e-5, "maxiter": 250})
    if not fit.success or not np.isfinite(fit.fun):
        raise RuntimeError("NFW profile fit failed")
    log_best, chi2_best = float(fit.x), float(fit.fun)
    target = chi2_best + 1.0
    def crossing(direction):
        inner = log_best
        step = 0.02
        bound = 12.5 if direction < 0 else 16.0
        while True:
            candidate = log_best + direction * step
            bounded = candidate <= bound if direction < 0 else candidate >= bound
            outer = float(bound) if bounded else float(candidate)
            if objective(outer) >= target:
                left, right = sorted((inner, outer))
                root = brentq(lambda value: objective(value) - target, left, right, xtol=1e-6)
                return float(root), False
            if bounded:
                return float(bound), True
            inner = outer
            step *= 2.0
    low, low_bounded = crossing(-1.0)
    high, high_bounded = crossing(+1.0)
    ci_status = ("conditional_delta_chi2_1_profile_interval_bound_limited"
                 if low_bounded or high_bounded else "conditional_delta_chi2_1_profile_interval")
    dof = len(indices) - 1
    return {
        "concentration_fixed": float(concentration),
        "fit_profile_bins_zero_based": indices.tolist(),
        "fit_radial_range_mpc": [float(edges[indices[0]]), float(edges[indices[-1] + 1])],
        "mass_definition": "M200c in Msun",
        "best_fit_log10_M200c_Msun": log_best,
        "best_fit_M200c_Msun": float(10.0 ** log_best),
        "conditional_profile_interval_status": ci_status,
        "conditional_profile_interval_log10_M200c": [low, high] if low is not None and high is not None else None,
        "chi2": chi2_best,
        "dof": dof,
        "nominal_chi2_p": float(chi2.sf(chi2_best, dof)) if dof > 0 else None,
        "profile_at_best_fit_msun_mpc2": model(log_best).tolist(),
        "observed_profile_fit_bins_msun_mpc2": y.tolist(),
        "covariance_min_eigenvalue": float(np.linalg.eigvalsh(cov).min()),
        "covariance_condition_number": float(np.linalg.cond(cov)),
        "fit_statistic": "full-covariance chi-square using the saved paired endpoint jackknife covariance; covariance uncertainty not propagated",
    }


def run(output: Path) -> dict:
    output = Path(output).resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite NFW mass output: {output}")
    pairs, random_doc, tomo = load_pairs_and_provenance()
    z_grid, nz_by_bin = load_hsc_y3_nz_sacc(NZ_SACC)
    cosmology = Cosmology(H0=70.0, Om=0.3)
    edges = np.asarray(random_doc["input"]["radial_edges_mpc"], dtype=float)
    if len(edges) != 8 or len(FIT_PROFILE_BINS) != 4:
        raise ValueError("unexpected saved radial-bin contract")
    prepared = prepare_estimator_weights(pairs, z_grid, nz_by_bin, edges, cosmology)
    result_by_width = {}
    for block_width in ("5.0", "2.5"):
        block = tomo["block_width_results"][block_width]
        excluded = block["source_zbin_1_exclusion_sensitivity"]
        observed_full = np.asarray(excluded["profile_msun_mpc2"], dtype=float)
        covariance_full = np.asarray(excluded["profile_covariance_msun2_mpc4"], dtype=float)
        if observed_full.shape != (7, 2) or covariance_full.shape != (14, 14):
            raise ValueError(f"unexpected saved source-bin-excluded profile shape for {block_width}°")
        obs_t = observed_full[:, 0]
        covariance_t = covariance_full[np.ix_(np.arange(0, 14, 2), np.arange(0, 14, 2))]
        fits = {}
        for concentration in FIXED_CONCENTRATIONS:
            fits[str(concentration)] = profile_fit(
                prepared, edges, cosmology,
                obs_t, covariance_t, concentration,
            )
        result_by_width[block_width] = {
            "n_lenses": 435,
            "n_random_centers": 2175,
            "n_endpoint_blocks": int(block["n_endpoint_blocks"]),
            "source_redshift_bins_used": [2, 3, 4],
            "source_zbin_1_excluded": True,
            "fits_by_fixed_concentration": fits,
            "omitted_profile_bins": {"inner": 0, "outer": [5, 6]},
        }

    report = {
        "schema": "astrafek-wl-random-subtracted-nfw-mass-screen-v1",
        "status": "complete_exploratory_model_screen",
        "science_measurement": False,
        "hsc_network_used": False,
        "images_read": False,
        "new_data_requested": False,
        "estimator": "single spherical NFW weak-shear DeltaSigma model averaged over the exact 435-lens source pairs, pair-wise HSC Y3 n(z) efficiencies, actual projected radii, and the same per-bin HSC shear/selection weights",
        "data_profile": "saved random-subtracted stack with source bin 1 excluded",
        "mass_interpretation": "effective stack-weighted single-halo M200c scale; not the arithmetic mean cluster mass",
        "cosmology": {"H0_km_s_Mpc": 70.0, "Omega_m": 0.3},
        "radial_fit_range_mpc": [0.2, 1.5],
        "fit_bins": list(FIT_PROFILE_BINS),
        "concentrations_fixed": list(FIXED_CONCENTRATIONS),
        "uncertainty": "conditional Delta-chi-square profile interval from paired endpoint-jackknife covariance; does not include concentration, source n(z), reduced-shear, miscentering, two-halo, random-aperture matching, or covariance uncertainty",
        "by_jackknife_block_width_deg": result_by_width,
        "input_sha256": {
            str(RANDOM_SUMMARY): sha(RANDOM_SUMMARY),
            str(RANDOM_EXPORT_AUDIT): sha(RANDOM_EXPORT_AUDIT),
            str(TOMO_SUMMARY): sha(TOMO_SUMMARY),
            str(EXPECTED_FULL_LENS_SUMMARY): sha(EXPECTED_FULL_LENS_SUMMARY),
            str(NZ_SACC): sha(NZ_SACC),
            str(LENS_PLAN / "row_export_plan.json"): sha(LENS_PLAN / "row_export_plan.json"),
        },
        "limitations": [
            "This fit uses a single-mass NFW forward model; the real stack spans a richness-selected mass distribution and is nonlinear in mass.",
            "The 0.2--1.5 Mpc one-halo range and fixed concentrations 3, 4, and 5 are diagnostic choices, not a validated mass-calibration model.",
            "The estimator uses weak shear; reduced shear and magnification are not modeled.",
            "The random-center depth/mask/PSF/sky match and photo-z systematic uncertainty remain unresolved.",
            "The nominal chi-square probabilities and profile intervals condition on the saved spatial covariance and do not constitute detection or science-grade mass calibration.",
        ],
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.mkdir(parents=True, exist_ok=False)
    path = output / "summary.json"
    path.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    report["summary_sha256"] = sha(path)
    return report


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    result = run(args.output)
    print(json.dumps({
        "status": result["status"], "science_measurement": False,
        "by_jackknife_block_width_deg": {
            width: {
                c: {
                    "M200c_Msun": fit["best_fit_M200c_Msun"],
                    "conditional_profile_interval_log10_M200c": fit["conditional_profile_interval_log10_M200c"],
                    "chi2_dof": [fit["chi2"], fit["dof"]],
                }
                for c, fit in value["fits_by_fixed_concentration"].items()
            }
            for width, value in result["by_jackknife_block_width_deg"].items()
        },
        "summary_sha256": result["summary_sha256"],
    }, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
