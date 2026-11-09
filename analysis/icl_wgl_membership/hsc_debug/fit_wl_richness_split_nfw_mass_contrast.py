#!/usr/bin/env python3
"""Fit two richness-group NFW masses to their random-subtracted contrast.

The observable is the redshift-matched high-minus-low DeltaSigma contrast
with its saved shared spatial-jackknife covariance. Model profiles are averaged
over the exact target lens-source rows and estimator weights in each richness
group; the mean random-aperture one-halo signal is zero. This is an exploratory
single-mass-per-group diagnostic, not a richness-mass calibration.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import sys

import numpy as np
from scipy.optimize import brentq, least_squares, minimize_scalar
from scipy.stats import chi2

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
ROOT = HERE.parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from fetch_ancillary_catalogues import sha
from fit_wl_random_stack_nfw_mass import DATA, load_pairs_and_provenance
from ogfkit.wlcluster import Cosmology, nfw_sigma
from s19a_shape_stack import (
    C_KMS, G_MPC, hsc_selection_bias, load_hsc_y3_nz_sacc,
    pair_geometry, tomographic_nz_lensing_efficiency,
)


SPLIT_SUMMARY = DATA / "wl_random_subtracted_richness_split_2026-11-05_02/summary.json"
OUTPUT = DATA / "wl_richness_split_nfw_mass_contrast_2026-11-05_03"
EXPECTED_SPLIT_SHA = "b38a78f6005dc17d8a12ddd6117debd63a005d0b21ff8dbda05961bd31eb65b7"
SOURCEBIN1_SUMMARY = DATA / "wl_random_subtracted_richness_split_zbin1_excluded_2026-11-05_01/summary.json"
EXPECTED_SOURCEBIN1_SHA = "814e9162aad25cc6c80c65a19ddfc52b88f1acc6f7d324723253395cea456df6"
NZ_SACC = DATA / "hsc_y3_nz_sacc_2026-10-19_01/hsc_y3_fourier_space_data_vector.sacc"
EDGES_MPC = np.asarray([0.1, 0.2, 0.4, 0.7, 1.0, 1.5, 2.0, 3.0])
FIT_BINS = np.asarray([1, 2, 3, 4], dtype=int)  # 0.2--1.5 proper Mpc
CONCENTRATIONS = (3.0, 4.0, 5.0)


def prepare_estimator_rows(pairs, selected_lenses, z_grid, nz_by_bin, cosmology,
                           source_bin_min=1):
    """Mirror the HSC Y3 stack estimator for source bins 1--4 or 2--4."""
    lens_id = np.asarray(pairs["lens_id"]).astype(str)
    source_bin = np.asarray(pairs["source_zbin"], dtype=int)
    select = (np.isin(lens_id, list(selected_lenses))
              & (source_bin >= source_bin_min) & (source_bin <= 4))
    sub = {key: np.asarray(value)[select] for key, value in pairs.items()}
    bins = source_bin[select]
    if len(bins) == 0:
        raise ValueError("no source-bin 1--4 rows in one richness cohort")
    beta = tomographic_nz_lensing_efficiency(
        sub["lens_z"], bins, z_grid, nz_by_bin, cosmology,
    )
    radius, _phi = pair_geometry(sub, cosmology)
    shape_amp = np.hypot(sub["e1"], sub["e2"])
    shape_ok = ((sub["sigma_e"] >= 0.0) & (sub["sigma_e"] <= 0.4)
                & (shape_amp < 2.0))
    source_ok = shape_ok & (sub["mag_a10"] <= 25.5) & (sub["resolution"] >= 0.3)
    usable = (np.asarray(sub["b_mode_pass"], dtype=bool) & source_ok & (beta > 0)
              & (radius >= EDGES_MPC[0]) & (radius <= EDGES_MPC[-1]))
    chi_factor = 4.0 * math.pi * G_MPC / (C_KMS * C_KMS)
    d_l = np.asarray(cosmology.angular_diameter(sub["lens_z"]), dtype=float)
    sigma_crit_inv = chi_factor * d_l * beta
    usable &= sigma_crit_inv > 0
    prepared = []
    used_lenses = set()
    for index, (lo, hi) in enumerate(zip(EDGES_MPC[:-1], EDGES_MPC[1:])):
        in_bin = usable & (radius >= lo) & (
            (radius < hi) if index < len(EDGES_MPC) - 2 else (radius <= hi)
        )
        if not np.any(in_bin):
            raise ValueError(f"empty usable lens-source bin {index}")
        weight = np.asarray(sub["weight"][in_bin], dtype=float)
        q = sigma_crit_inv[in_bin]
        lensing_weight = weight * q * q
        rms = np.asarray(sub["rms_e"][in_bin], dtype=float)
        responsivity = float(1.0 - np.sum(lensing_weight * rms * rms) / np.sum(lensing_weight))
        if not np.isfinite(responsivity) or responsivity <= 0:
            raise ValueError(f"invalid HSC responsivity in radial bin {index}")
        selection = hsc_selection_bias(
            lensing_weight, sub["mag_a10"][in_bin], sub["resolution"][in_bin],
            sub["psf_shape11"][in_bin], sub["psf_shape22"][in_bin],
            sub["psf_shape12"][in_bin],
        )
        m = np.asarray(sub["m"][in_bin], dtype=float)
        w_shear = weight * (
            2.0 * responsivity * (1.0 + m) * (1.0 + float(selection["m_sel"]))
        ) ** 2
        prepared.append({
            "bin": index,
            "radius_mpc": np.asarray(radius[in_bin], dtype=float),
            "lens_z": np.asarray(sub["lens_z"][in_bin], dtype=float),
            "effective_weight": np.asarray(w_shear * q * q, dtype=float),
            "weight_sum": float(np.sum(w_shear * q * q)),
        })
        used_lenses.update(str(x) for x in np.unique(np.asarray(sub["lens_id"])[in_bin]))
    return prepared, used_lenses


def model_profile(prepared, cosmology, mass_msun, concentration):
    values = np.full(len(prepared), np.nan)
    for item in prepared:
        z = item["lens_z"]
        hz = cosmology.H0 * np.sqrt(cosmology.Om * (1.0 + z) ** 3 + (1.0 - cosmology.Om))
        rho_crit = 3.0 * hz * hz / (8.0 * math.pi * G_MPC)
        r200 = (3.0 * mass_msun / (4.0 * math.pi * 200.0 * rho_crit)) ** (1.0 / 3.0)
        rs = r200 / concentration
        delta_c = (200.0 / 3.0) * concentration**3 / (
            np.log1p(concentration) - concentration / (1.0 + concentration)
        )
        _sigma, delta_sigma = nfw_sigma(item["radius_mpc"], delta_c * rho_crit, rs)
        values[item["bin"]] = float(
            np.sum(item["effective_weight"] * delta_sigma) / item["weight_sum"]
        )
    return values


def fit_two_group_masses(prepared_low, prepared_high, observed, covariance,
                         cosmology, concentration):
    indices = FIT_BINS
    y = np.asarray(observed, dtype=float)[indices]
    cov = np.asarray(covariance, dtype=float)[np.ix_(indices, indices)]
    if np.min(np.linalg.eigvalsh(cov)) <= 0:
        raise ValueError("saved contrast covariance is not positive definite")
    chol = np.linalg.cholesky(cov)
    cache = {"low": {}, "high": {}}

    def one(group, prepared, log_mass):
        key = float(log_mass)
        if key not in cache[group]:
            cache[group][key] = model_profile(
                prepared, cosmology, 10.0**key, concentration,
            )[indices]
        return cache[group][key]

    def residual(log_masses):
        pred = (one("high", prepared_high, log_masses[1])
                - one("low", prepared_low, log_masses[0]))
        return np.linalg.solve(chol, y - pred)

    fit = least_squares(
        residual, x0=np.asarray([14.0, 14.3]), bounds=([12.5, 12.5], [16.0, 16.0]),
        method="trf", x_scale="jac", ftol=2e-6, xtol=2e-6, gtol=2e-6,
        max_nfev=120,
    )
    if not fit.success or not np.all(np.isfinite(fit.x)):
        raise RuntimeError(f"two-mass NFW fit failed: {fit.message}")
    best = fit.x.copy()
    best_chi2 = float(residual(best) @ residual(best))
    target = best_chi2 + 1.0

    def profiled_chi(log_mass, fixed_index):
        free_index = 1 - fixed_index

        def objective(free_mass):
            trial = best.copy()
            trial[fixed_index] = log_mass
            trial[free_index] = float(free_mass)
            value = residual(trial)
            return float(value @ value)

        prof = minimize_scalar(objective, bounds=(12.5, 16.0), method="bounded",
                               options={"xatol": 2e-4, "maxiter": 100})
        if not prof.success:
            raise RuntimeError("profile likelihood minimization failed")
        return float(prof.fun)

    def interval(fixed_index):
        center = float(best[fixed_index])
        ends = []
        for direction in (-1.0, 1.0):
            inner, step = center, 0.02
            bound = 12.5 if direction < 0 else 16.0
            while True:
                trial = center + direction * step
                bounded = trial <= bound if direction < 0 else trial >= bound
                outer = float(bound) if bounded else float(trial)
                if profiled_chi(outer, fixed_index) >= target:
                    root = brentq(
                        lambda value: profiled_chi(value, fixed_index) - target,
                        min(inner, outer), max(inner, outer), xtol=2e-4,
                    )
                    ends.append(float(root))
                    break
                if bounded:
                    ends.append(bound)
                    break
                inner, step = outer, step * 2.0
        return ends

    def ratio_profile(log_ratio):
        # Profile over the common mass scale while preserving both mass bounds.
        low_min = max(12.5, 12.5 - log_ratio)
        low_max = min(16.0, 16.0 - log_ratio)
        if low_min >= low_max:
            return float("inf")

        def objective(log_low):
            trial = np.asarray([float(log_low), float(log_low + log_ratio)])
            value = residual(trial)
            return float(value @ value)

        prof = minimize_scalar(objective, bounds=(low_min, low_max), method="bounded",
                               options={"xatol": 2e-4, "maxiter": 120})
        if not prof.success:
            raise RuntimeError("mass-ratio profile minimization failed")
        return float(prof.fun)

    best_log_ratio = float(best[1] - best[0])

    def ratio_interval():
        ends = []
        limits = (-3.5, 3.5)
        for direction in (-1.0, 1.0):
            inner, step = best_log_ratio, 0.02
            bound = limits[0] if direction < 0 else limits[1]
            while True:
                trial = best_log_ratio + direction * step
                bounded = trial <= bound if direction < 0 else trial >= bound
                outer = float(bound) if bounded else float(trial)
                if ratio_profile(outer) >= target:
                    root = brentq(
                        lambda value: ratio_profile(value) - target,
                        min(inner, outer), max(inner, outer), xtol=2e-4,
                    )
                    ends.append(float(root))
                    break
                if bounded:
                    ends.append(bound)
                    break
                inner, step = outer, step * 2.0
        return ends

    low_profile = model_profile(prepared_low, cosmology, 10.0**best[0], concentration)
    high_profile = model_profile(prepared_high, cosmology, 10.0**best[1], concentration)
    return {
        "concentration_fixed": float(concentration),
        "lower_richness_log10_M200c_Msun": float(best[0]),
        "lower_richness_M200c_Msun": float(10.0**best[0]),
        "lower_richness_conditional_delta_chi2_1_interval_log10_M200c": interval(0),
        "upper_richness_log10_M200c_Msun": float(best[1]),
        "upper_richness_M200c_Msun": float(10.0**best[1]),
        "upper_richness_conditional_delta_chi2_1_interval_log10_M200c": interval(1),
        "mass_ratio_upper_over_lower": float(10.0 ** (best[1] - best[0])),
        "mass_ratio_conditional_delta_chi2_1_interval_log10": ratio_interval(),
        "delta_chi2_at_equal_group_masses": float(ratio_profile(0.0) - best_chi2),
        "conditional_equal_mass_profile_p_1d": float(
            chi2.sf(ratio_profile(0.0) - best_chi2, 1)
        ),
        "chi2": best_chi2,
        "dof": int(len(indices) - 2),
        "nominal_conditional_fit_p": float(chi2.sf(best_chi2, len(indices) - 2)),
        "observed_high_minus_low_fit_bins_msun_mpc2": y.tolist(),
        "model_high_minus_low_fit_bins_msun_mpc2": (high_profile - low_profile)[indices].tolist(),
        "fit_profile_bins_zero_based": indices.tolist(),
        "fit_radius_mpc_proper": [0.2, 1.5],
    }


def run(output: Path, exclude_source_bin_1=False):
    output = Path(output).resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite mass-contrast output: {output}")
    base_split = json.loads(SPLIT_SUMMARY.read_text(encoding="utf-8"))
    if (sha(SPLIT_SUMMARY) != EXPECTED_SPLIT_SHA
            or base_split.get("status") != "complete_exploratory"
            or base_split.get("science_measurement") is not False
            or base_split.get("n_matched_lens_pairs") != 204):
        raise ValueError("richness-split input is changed or outside its audited scope")
    if exclude_source_bin_1:
        split = json.loads(SOURCEBIN1_SUMMARY.read_text(encoding="utf-8"))
        if (sha(SOURCEBIN1_SUMMARY) != EXPECTED_SOURCEBIN1_SHA
                or split.get("status") != "complete_exploratory_sensitivity"
                or split.get("science_measurement") is not False
                or split.get("source_bin_1_excluded_from_target_and_random") is not True
                or split.get("provenance", {}).get("all_source_bins_split_summary_sha256") != EXPECTED_SPLIT_SHA):
            raise ValueError("source-bin-1-excluded profile is changed or incompatible")
    else:
        split = base_split
    target_pairs, _random_doc, _tomo = load_pairs_and_provenance()
    low_ids = {str(row["lower_lens_id"]) for row in base_split["matched_lens_pairs"]}
    high_ids = {str(row["upper_lens_id"]) for row in base_split["matched_lens_pairs"]}
    all_ids = {str(x) for x in np.unique(target_pairs["lens_id"])}
    if not low_ids <= all_ids or not high_ids <= all_ids or low_ids & high_ids:
        raise ValueError("matched richness lens IDs do not match the verified target rows")
    nz_z, nz_bins = load_hsc_y3_nz_sacc(NZ_SACC)
    cosmology = Cosmology(H0=70.0, Om=0.3)
    prepared_low, used_low = prepare_estimator_rows(
        target_pairs, low_ids, nz_z, nz_bins, cosmology,
        source_bin_min=2 if exclude_source_bin_1 else 1,
    )
    prepared_high, used_high = prepare_estimator_rows(
        target_pairs, high_ids, nz_z, nz_bins, cosmology,
        source_bin_min=2 if exclude_source_bin_1 else 1,
    )
    if len(used_low) != 204 or len(used_high) != 204:
        raise ValueError("target model rows do not include all 204 usable redshift-matched lenses")

    observed = np.asarray([
        row["upper_minus_lower_delta_sigma_msun_mpc2"] for row in split["profile"]
    ], dtype=float)
    full_cov = np.asarray(split["upper_minus_lower_covariance"], dtype=float)
    if full_cov.shape != (14, 14):
        raise ValueError("richness contrast covariance does not match seven tangential/cross bins")
    tangent_indices = np.arange(0, 14, 2)
    covariance_t = full_cov[np.ix_(tangent_indices, tangent_indices)]
    fits = {
        str(concentration): fit_two_group_masses(
            prepared_low, prepared_high, observed, covariance_t,
            cosmology, concentration,
        )
        for concentration in CONCENTRATIONS
    }
    result = {
        "schema": "astrafek-wl-redshift-matched-richness-two-mass-nfw-contrast-v1",
        "status": "complete_exploratory_model_screen",
        "science_measurement": False,
        "hsc_network_used": False,
        "images_read": False,
        "new_data_requested": False,
        "analysis_script": str(Path(__file__).resolve()),
        "analysis_script_sha256": sha(Path(__file__).resolve()),
        "method": "fit M200c_low and M200c_high simultaneously to the random-subtracted high-minus-low DeltaSigma contrast using the saved shared spatial-jackknife covariance and exact group-specific weak-shear estimator weights",
        "data_source_bins": [2, 3, 4] if exclude_source_bin_1 else [1, 2, 3, 4],
        "source_bin_1_excluded": bool(exclude_source_bin_1),
        "mass_definition": "effective stack-weighted single-halo M200c in Msun",
        "radial_fit_range_mpc_proper": [0.2, 1.5],
        "n_lenses_per_group": 204,
        "concentrations_fixed": list(CONCENTRATIONS),
        "fits_by_fixed_concentration": fits,
        "input_sha256": {
            str(SPLIT_SUMMARY): sha(SPLIT_SUMMARY),
            **({str(SOURCEBIN1_SUMMARY): sha(SOURCEBIN1_SUMMARY)} if exclude_source_bin_1 else {}),
            str(NZ_SACC): sha(NZ_SACC),
        },
        "limitations": [
            "single spherical NFW mass per richness group does not describe the richness-selected mass distribution",
            "source bin 1 is retained to match the measured split profile; the earlier full-sample NFW screen excluded it, so this fit requires a source-bin-1-excluded sensitivity before interpretation",
            "the 0.2--1.5 Mpc range may include two-halo signal at its outer end",
            "reduced shear, magnification, miscentering, photo-z calibration uncertainty, and covariance uncertainty are not propagated in this contrast fit",
            "matched-aperture depth/mask/PSF/sky and validated contiguous survey-region covariance remain unresolved",
            "conditional profile intervals and nominal chi-square probabilities are diagnostic only; this is not a calibrated richness-mass relation or detection",
        ],
    }
    output.mkdir(parents=True, exist_ok=False)
    path = output / "summary.json"
    path.write_text(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False) + "\n")
    result["summary_sha256"] = sha(path)
    print(json.dumps({
        "status": result["status"],
        "science_measurement": False,
        "fits": {
            c: {
                "low_M200c_Msun": fit["lower_richness_M200c_Msun"],
                "high_M200c_Msun": fit["upper_richness_M200c_Msun"],
                "mass_ratio": fit["mass_ratio_upper_over_lower"],
                "chi2_dof": [fit["chi2"], fit["dof"]],
            }
            for c, fit in fits.items()
        },
        "summary_sha256": result["summary_sha256"],
    }, indent=2))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument("--exclude-source-bin-1", action="store_true")
    args = parser.parse_args()
    run(args.output, args.exclude_source_bin_1)


if __name__ == "__main__":
    main()
