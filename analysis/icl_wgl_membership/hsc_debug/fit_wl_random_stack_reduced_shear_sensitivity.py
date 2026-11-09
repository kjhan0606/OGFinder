#!/usr/bin/env python3
"""Reduced-shear sensitivity of the saved random-subtracted stack.

This is a local model sensitivity, not an HSC query or science measurement.
It retains the exact source-pair radii, lens redshifts, weak-shear estimator
weights, calibrated tomographic n(z), fitted radial bins, and saved spatial
jackknife covariance used by the adjacent single-halo NFW screen. The source
redshift integral uses moments through ninth order, expanding the reduced-
shear denominator through eighth order; the geometric-series tail is bounded.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.optimize import brentq, minimize_scalar
from scipy.stats import chi2

from fit_wl_random_stack_nfw_mass import (
    DATA, EXPECTED_FULL_LENS_SUMMARY, FIT_PROFILE_BINS, LENS_PLAN, NZ_SACC,
    RANDOM_EXPORT_AUDIT, RANDOM_SUMMARY, TOMO_SUMMARY,
    load_pairs_and_provenance,
)
from fetch_ancillary_catalogues import sha
from ogfkit.wlcluster import Cosmology, nfw_sigma
from s19a_shape_stack import (
    C_KMS, G_MPC, hsc_selection_bias, load_hsc_y3_nz_sacc,
    pair_geometry, tomographic_nz_lensing_efficiency,
)


DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
BASE_NFW = DATA / "wl_random_stack_nfw_mass_2026-11-05_01/summary.json"
DEFAULT_OUTPUT = DATA / "wl_random_stack_reduced_shear_sensitivity_2026-11-05_03"
SERIES_ORDER = 8
EXPECTED_BASE_NFW_SHA = "d3ca78a24fc030f64694247fd45d33f607d89c18fb102b6b5a6879d53c52cd4b"


def beta_moments_by_pair(lens_z, source_bin, z_grid, nz_by_bin, cosmology,
                         max_power=SERIES_ORDER + 1):
    """Return calibrated n(z)-weighted <(D_ls/D_s)^p>, p=1..max_power."""
    lens_z = np.asarray(lens_z, dtype=float)
    source_bin = np.asarray(source_bin, dtype=int)
    z = np.asarray(z_grid, dtype=float)
    nz = np.asarray(nz_by_bin, dtype=float)
    norm = np.trapezoid(nz, z, axis=1)
    pdf = nz / norm[:, None]
    chi_s = np.asarray(cosmology.comoving(z), dtype=float)
    unique_lens_z, inverse = np.unique(lens_z, return_inverse=True)
    chi_l = np.asarray(cosmology.comoving(unique_lens_z), dtype=float)
    result_by_lens_bin = np.zeros(
        (len(unique_lens_z), nz.shape[0], max_power), dtype=float,
    )
    for bin_index in range(1, nz.shape[0] + 1):
        efficiency = np.divide(
            chi_s[None, :] - chi_l[:, None], chi_s[None, :],
            out=np.zeros((len(unique_lens_z), len(z)), dtype=float),
            where=(z[None, :] > unique_lens_z[:, None]) & (chi_s[None, :] > 0),
        )
        efficiency = np.clip(efficiency, 0.0, 1.0)
        for power in range(1, max_power + 1):
            result_by_lens_bin[:, bin_index - 1, power - 1] = np.trapezoid(
                pdf[bin_index - 1][None, :] * efficiency**power, z, axis=1,
            )
    out = np.zeros((len(lens_z), max_power), dtype=float)
    used = (source_bin >= 1) & (source_bin <= nz.shape[0])
    out[used] = result_by_lens_bin[inverse[used], source_bin[used] - 1, :]
    return out


def prepare_rows(pairs, z_grid, nz_by_bin, edges, cosmology):
    """Mirror the saved stack cuts/calibration and add <beta^2>/<beta>."""
    source_bin = np.asarray(pairs["source_zbin"], dtype=int)
    keep = (source_bin >= 2) & (source_bin <= 4)
    sub = {key: np.asarray(value)[keep] for key, value in pairs.items()}
    bins = source_bin[keep]
    beta1 = tomographic_nz_lensing_efficiency(
        sub["lens_z"], bins, z_grid, nz_by_bin, cosmology,
    )
    beta_moments = beta_moments_by_pair(
        sub["lens_z"], bins, z_grid, nz_by_bin, cosmology,
    )
    radius, phi = pair_geometry(sub, cosmology)
    amplitude = np.hypot(sub["e1"], sub["e2"])
    shape_ok = ((sub["sigma_e"] >= 0.0) & (sub["sigma_e"] <= 0.4)
                & (amplitude < 2.0))
    bmode = np.asarray(sub["b_mode_pass"], dtype=bool)
    source_ok = shape_ok & (sub["mag_a10"] <= 25.5) & (sub["resolution"] >= 0.3)
    chi_factor = 4.0 * math.pi * G_MPC / (C_KMS * C_KMS)
    d_l = np.asarray(cosmology.angular_diameter(sub["lens_z"]), dtype=float)
    geometry = chi_factor * d_l
    q_eff = geometry * beta1
    usable = (bmode & source_ok & (beta1 > 0) & (beta_moments[:, 1] > 0)
              & (radius >= edges[0]) & (radius <= edges[-1]) & (q_eff > 0))
    rows = []
    for index, (lo, hi) in enumerate(zip(edges[:-1], edges[1:])):
        in_bin = usable & (radius >= lo) & ((radius < hi) if index < len(edges) - 2 else (radius <= hi))
        if not np.any(in_bin):
            raise ValueError(f"no usable source rows in radial bin {index}")
        weight = np.asarray(sub["weight"][in_bin], dtype=float)
        q = q_eff[in_bin]
        lensing_weight = weight * q * q
        rms = np.asarray(sub["rms_e"][in_bin], dtype=float)
        responsivity = float(1.0 - np.sum(lensing_weight * rms * rms) / np.sum(lensing_weight))
        if not np.isfinite(responsivity) or responsivity <= 0:
            raise ValueError(f"invalid responsivity in radial bin {index}")
        selection = hsc_selection_bias(
            lensing_weight, sub["mag_a10"][in_bin], sub["resolution"][in_bin],
            sub["psf_shape11"][in_bin], sub["psf_shape22"][in_bin],
            sub["psf_shape12"][in_bin],
        )
        m = np.asarray(sub["m"][in_bin], dtype=float)
        w_shear = weight * (
            2.0 * responsivity * (1.0 + m) * (1.0 + float(selection["m_sel"]))
        ) ** 2
        rows.append({
            "bin": index,
            "radius_mpc": np.asarray(radius[in_bin], dtype=float),
            "lens_z": np.asarray(sub["lens_z"][in_bin], dtype=float),
            "geometry": np.asarray(geometry[in_bin], dtype=float),
            "q_eff": np.asarray(q, dtype=float),
            "beta_moment_ratios": np.asarray(
                beta_moments[in_bin] / beta1[in_bin, None], dtype=float,
            ),
            "weight": np.asarray(w_shear, dtype=float),
            "weight_sum_q2": float(np.sum(w_shear * q * q)),
        })
    return rows


def model_profile(rows, cosmology, mass_msun, concentration, series_order=SERIES_ORDER):
    values = np.full(len(rows), np.nan)
    max_kappa_bound = 0.0
    for item in rows:
        z = item["lens_z"]
        hz = cosmology.H0 * np.sqrt(cosmology.Om * (1.0 + z) ** 3 + (1.0 - cosmology.Om))
        rho_crit = 3.0 * hz * hz / (8.0 * math.pi * G_MPC)
        r200 = (3.0 * mass_msun / (4.0 * math.pi * 200.0 * rho_crit)) ** (1.0 / 3.0)
        rs = r200 / concentration
        delta_c = (200.0 / 3.0) * concentration**3 / (
            math.log1p(concentration) - concentration / (1.0 + concentration))
        sigma, delta_sigma = nfw_sigma(item["radius_mpc"], delta_c * rho_crit, rs)
        kappa_max_beta1 = np.asarray(sigma) * item["geometry"]
        max_kappa_bound = max(max_kappa_bound, float(np.max(kappa_max_beta1)))
        # g = gamma/(1-kappa); integrate beta over calibrated n(z) with
        # sum_n kappa_inf**n <beta**(n+1)> / <beta>.
        kappa_inf = np.asarray(sigma) * item["geometry"]
        correction = np.zeros_like(kappa_inf)
        for order in range(series_order + 1):
            correction += (kappa_inf**order
                           * item["beta_moment_ratios"][:, order])
        ds_reduced = np.asarray(delta_sigma) * correction
        values[item["bin"]] = float(
            np.sum(item["weight"] * item["q_eff"]**2 * ds_reduced)
            / item["weight_sum_q2"]
        )
    return values, max_kappa_bound


def fit_profile(rows, cosmology, observed, covariance, concentration,
                series_order=SERIES_ORDER):
    indices = np.asarray(FIT_PROFILE_BINS, dtype=int)
    y = np.asarray(observed, dtype=float)[indices]
    cov = np.asarray(covariance, dtype=float)[np.ix_(indices, indices)]
    if np.min(np.linalg.eigvalsh(cov)) <= 0:
        raise ValueError("saved tangential covariance is not positive definite")
    inv_cov = np.linalg.inv(cov)
    cache = {}

    def evaluate(log_mass):
        key = float(log_mass)
        if key not in cache:
            cache[key] = model_profile(
                rows, cosmology, 10.0**key, concentration, series_order,
            )
        return cache[key]

    def objective(log_mass):
        pred, _ = evaluate(log_mass)
        delta = y - pred[indices]
        return float(delta @ inv_cov @ delta)

    result = minimize_scalar(objective, bounds=(12.5, 16.0), method="bounded",
                             options={"xatol": 1e-5, "maxiter": 250})
    if not result.success or not np.isfinite(result.fun):
        raise RuntimeError("reduced-shear sensitivity fit failed")
    log_best, best_chi2 = float(result.x), float(result.fun)
    target = best_chi2 + 1.0

    def crossing(direction):
        inner, step = log_best, 0.02
        bound = 12.5 if direction < 0 else 16.0
        while True:
            candidate = log_best + direction * step
            bounded = candidate <= bound if direction < 0 else candidate >= bound
            outer = float(bound) if bounded else float(candidate)
            if objective(outer) >= target:
                return float(brentq(lambda value: objective(value) - target,
                                    min(inner, outer), max(inner, outer), xtol=1e-6))
            if bounded:
                return float(bound)
            inner, step = outer, step * 2.0

    lower, upper = crossing(-1.0), crossing(+1.0)
    pred, max_kappa_bound = evaluate(log_best)
    return {
        "concentration_fixed": float(concentration),
        "best_fit_M200c_Msun": float(10.0**log_best),
        "conditional_profile_interval_log10_M200c": [lower, upper],
        "chi2": best_chi2,
        "dof": len(indices) - 1,
        "nominal_chi2_p_conditional": float(chi2.sf(best_chi2, len(indices) - 1)),
        "model_fit_bins_msun_mpc2": pred[indices].tolist(),
        "max_kappa_bound_beta_le_1_at_best_fit": max_kappa_bound,
        "relative_truncation_bound": (
            max_kappa_bound**(series_order + 1) / (1.0 - max_kappa_bound)
            if max_kappa_bound < 1.0 else None
        ),
        "reduced_shear_series_order": int(series_order),
    }


def run(output: Path):
    output = Path(output).resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite reduced-shear output: {output}")
    base = json.loads(BASE_NFW.read_text(encoding="utf-8"))
    if sha(BASE_NFW) != EXPECTED_BASE_NFW_SHA or base.get("science_measurement") is not False:
        raise ValueError("base NFW screen changed or has an unexpected science gate")
    pairs, _random_doc, tomo = load_pairs_and_provenance()
    z_grid, nz_by_bin = load_hsc_y3_nz_sacc(NZ_SACC)
    cosmology = Cosmology(H0=70.0, Om=0.3)
    edges = np.asarray([0.1, 0.2, 0.4, 0.7, 1.0, 1.5, 2.0, 3.0], dtype=float)
    rows = prepare_rows(pairs, z_grid, nz_by_bin, edges, cosmology)
    result_by_width = {}
    for width in ("5.0", "2.5"):
        block = tomo["block_width_results"][width]
        excluded = block["source_zbin_1_exclusion_sensitivity"]
        observed = np.asarray(excluded["profile_msun_mpc2"], dtype=float)[:, 0]
        full_cov = np.asarray(excluded["profile_covariance_msun2_mpc4"], dtype=float)
        cov_t = full_cov[np.ix_(np.arange(0, 14, 2), np.arange(0, 14, 2))]
        result_by_width[width] = {
            str(c): fit_profile(rows, cosmology, observed, cov_t, c, SERIES_ORDER)
            for c in (3.0, 4.0, 5.0)
        }
    out = {
        "schema": "astrafek-wl-random-subtracted-reduced-shear-sensitivity-v1",
        "status": "complete_exploratory_sensitivity",
        "science_measurement": False,
        "hsc_network_used": False,
        "images_read": False,
        "new_data_requested": False,
        "model": "single spherical NFW; source-redshift-integrated reduced shear using calibrated n(z) beta moments through power 9 (geometric series order 8); exact estimator-pair radii and weak-shear calibration weights retained",
        "data_profile": "saved random-subtracted 435-lens stack, source bin 1 excluded",
        "fit_radial_range_mpc_proper": [0.2, 1.5],
        "source_bins": [2, 3, 4],
        "jackknife_block_width_deg": result_by_width,
        "limitations": [
            "reduced-shear geometric series is truncated at order 8; per-fit worst-case beta<=1 tail bounds are reported",
            "single-mass NFW does not represent the richness-selected mass distribution or its nonlinear stack average",
            "concentration, magnification, miscentering, two-halo structure, photo-z calibration uncertainty, random-aperture matching, and covariance uncertainty remain unresolved",
            "conditional profile intervals and nominal chi-square probabilities are exploratory and do not constitute a detection or mass calibration",
        ],
        "input_sha256": {
            str(BASE_NFW): sha(BASE_NFW),
            str(RANDOM_SUMMARY): sha(RANDOM_SUMMARY),
            str(RANDOM_EXPORT_AUDIT): sha(RANDOM_EXPORT_AUDIT),
            str(TOMO_SUMMARY): sha(TOMO_SUMMARY),
            str(EXPECTED_FULL_LENS_SUMMARY): sha(EXPECTED_FULL_LENS_SUMMARY),
            str(NZ_SACC): sha(NZ_SACC),
            str(LENS_PLAN / "row_export_plan.json"): sha(LENS_PLAN / "row_export_plan.json"),
        },
        "science_measurement_gate_reason": "exploratory reduced-shear sensitivity only; independent depth, photo-z, random-aperture, concentration-distribution, and covariance gates remain open",
        "analysis_script": str(Path(__file__).resolve()),
        "analysis_script_sha256": sha(Path(__file__).resolve()),
    }
    output.mkdir(parents=True, exist_ok=False)
    path = output / "summary.json"
    path.write_text(json.dumps(out, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    out["summary_sha256"] = sha(path)
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    result = run(args.output)
    compact = {
        width: {
            concentration: {
                "mass_Msun": fit["best_fit_M200c_Msun"],
                "chi2_dof": [fit["chi2"], fit["dof"]],
                "max_kappa": fit["max_kappa_bound_beta_le_1_at_best_fit"],
                "series_tail_bound": fit["relative_truncation_bound"],
            }
            for concentration, fit in models.items()
        }
        for width, models in result["jackknife_block_width_deg"].items()
    }
    print(json.dumps({"science_measurement": False, "fits": compact,
                      "summary_sha256": result["summary_sha256"]}, indent=2))


if __name__ == "__main__":
    main()
