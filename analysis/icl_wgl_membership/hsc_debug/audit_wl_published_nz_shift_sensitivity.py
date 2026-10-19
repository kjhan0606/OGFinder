#!/usr/bin/env python3
"""Approximate HSC Y3 tomography sensitivity to published DESI mean shifts.

This applies simple constant translations to only the existing local SACC
curves for source bins 3 and 4, using the central and marginal 68% endpoint
higher-redshift shifts reported by Choppin de Janvry et al. (2026). It is not
a substitute for the authors' calibrated spline distributions; those files
require separate user approval and are not downloaded here.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
ROOT = HERE.parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from audit_wl_tomographic_consistency import tomographic_spatial_jackknife
from audit_wl_tomographic_lens_redshift import holm_adjust
from fetch_ancillary_catalogues import sha
from ogfkit.wlcluster import Cosmology
from run_wl_shape_stack import verify_and_load_pairs
from s19a_shape_stack import load_hsc_y3_nz_sacc


PUBLISHED_MEAN_SHIFTS = {1: 0.0, 2: 0.0, 3: 0.043, 4: 0.050}
# The paper reports Delta-z in the opposite sign convention: negative values
# mean the calibrated n(z) shifts to higher redshift. Convert its asymmetric
# intervals to positive higher-z translation amplitudes.
PUBLISHED_HIGH_Z_SHIFT_INTERVALS = {3: (0.021, 0.063), 4: (0.038, 0.062)}
PUBLISHED_SHIFT_SCENARIOS = {
    "central": {1: 0.0, 2: 0.0, 3: 0.043, 4: 0.050},
    "low3_low4": {1: 0.0, 2: 0.0, 3: 0.021, 4: 0.038},
    "low3_high4": {1: 0.0, 2: 0.0, 3: 0.021, 4: 0.062},
    "high3_low4": {1: 0.0, 2: 0.0, 3: 0.063, 4: 0.038},
    "high3_high4": {1: 0.0, 2: 0.0, 3: 0.063, 4: 0.062},
}
_integral = getattr(np, "trapezoid", np.trapz)


def translate_nz_curve(redshift, nz, mean_shift):
    """Shift an n(z) curve by ``mean_shift`` in redshift and renormalize."""
    redshift = np.asarray(redshift, dtype=float)
    nz = np.asarray(nz, dtype=float)
    shift = float(mean_shift)
    if (redshift.ndim != 1 or nz.shape != redshift.shape or len(redshift) < 2
            or not np.all(np.isfinite(redshift)) or not np.all(np.diff(redshift) > 0)
            or not np.all(np.isfinite(nz)) or np.any(nz < 0)
            or not np.isfinite(shift)):
        raise ValueError("redshift grid and n(z) must be finite, non-negative, and matching")
    shifted = np.interp(redshift - shift, redshift, nz, left=0.0, right=0.0)
    norm_before = float(_integral(nz, redshift))
    norm_after = float(_integral(shifted, redshift))
    if norm_before <= 0 or norm_after <= 0:
        raise ValueError("n(z) has no positive finite integral after translation")
    shifted /= norm_after
    mean_before = float(_integral(redshift * nz, redshift) / norm_before)
    mean_after = float(_integral(redshift * shifted, redshift))
    return shifted, {
        "requested_mean_shift": shift,
        "mean_redshift_before": mean_before,
        "mean_redshift_after": mean_after,
        "realized_mean_shift_after_grid_truncation": mean_after - mean_before,
        "integral_before_normalization": norm_before,
        "integral_after_truncation_before_renormalization": norm_after,
        "renormalized_integral": float(np.trapezoid(shifted, redshift)),
    }


def compare_profiles(baseline, shifted):
    def by_key(rows):
        return {
            (int(row["source_zbin"]), float(row["r_min_mpc"]), float(row["r_max_mpc"])): row
            for row in rows
        }
    base_rows, shifted_rows = by_key(baseline), by_key(shifted)
    if base_rows.keys() != shifted_rows.keys():
        raise ValueError("baseline and shifted profiles have different source-bin/radius keys")
    output = []
    for key in sorted(base_rows):
        left, right = base_rows[key], shifted_rows[key]
        base_value = float(left["delta_sigma_weak_msun_mpc2"])
        shifted_value = float(right["delta_sigma_weak_msun_mpc2"])
        output.append({
            "source_zbin": key[0],
            "r_min_mpc": key[1],
            "r_max_mpc": key[2],
            "baseline_delta_sigma_msun_mpc2": base_value,
            "shifted_delta_sigma_msun_mpc2": shifted_value,
            "delta_sigma_change_msun_mpc2": shifted_value - base_value,
            "baseline_jackknife_error_msun_mpc2": float(left["delta_sigma_jackknife_error_msun_mpc2"]),
            "shifted_jackknife_error_msun_mpc2": float(right["delta_sigma_jackknife_error_msun_mpc2"]),
        })
    return output


def run_sensitivity(pairs, provenance, nz_redshift, nz_by_bin, baseline,
                    radial_edges_mpc, cosmology, *, cell_width_deg,
                    mean_shifts=PUBLISHED_MEAN_SHIFTS, scenario_label="central"):
    if baseline.get("schema") != "astrafek-hsc-y3-tomographic-consistency-jackknife-v1":
        raise ValueError("baseline summary is not the expected full-sample tomography audit")
    if baseline.get("science_measurement") is not False:
        raise ValueError("baseline has an unexpected science gate")
    baseline_input = baseline.get("input", {})
    for key in ("row_plan_sha256", "row_batch_sha256", "n_pair_rows", "n_lenses", "masked_value"):
        observed_key = "b_mode_masked_raw_value" if key == "masked_value" else key
        if baseline_input.get(observed_key) != provenance.get(key):
            raise ValueError(f"baseline and current verified input differ at {key}")

    shifts = {}
    shifted = np.asarray(nz_by_bin, dtype=float).copy()
    if shifted.shape != (4, len(nz_redshift)):
        raise ValueError("expected the four official HSC Y3 n(z) curves")
    for source_bin in (1, 2, 3, 4):
        delta = float(mean_shifts[source_bin])
        if delta:
            curve, detail = translate_nz_curve(nz_redshift, shifted[source_bin - 1], delta)
            shifted[source_bin - 1] = curve
        else:
            curve = shifted[source_bin - 1]
            norm = float(_integral(curve, nz_redshift))
            mean = float(_integral(nz_redshift * curve, nz_redshift) / norm)
            detail = {
                "requested_mean_shift": 0.0,
                "mean_redshift_before": mean,
                "mean_redshift_after": mean,
                "realized_mean_shift_after_grid_truncation": 0.0,
                "renormalized_integral": 1.0,
            }
        shifts[str(source_bin)] = detail

    result = tomographic_spatial_jackknife(
        pairs, nz_redshift, shifted, radial_edges_mpc, cosmology,
        cell_width_deg=cell_width_deg,
    )
    raw_p = [float(row["p_value"]) for row in result["tomographic_delta_sigma_consistency"]]
    adjusted = holm_adjust(raw_p)
    for row, value in zip(result["tomographic_delta_sigma_consistency"], adjusted):
        row["holm_p_value_7_radial_tests"] = float(value)

    baseline_test_rows = baseline["tomographic_delta_sigma_consistency"]
    baseline_adjusted = holm_adjust([float(row["p_value"]) for row in baseline_test_rows])
    baseline_tests = {
        (float(row["r_min_mpc"]), float(row["r_max_mpc"])): row
        for row in baseline_test_rows
    }
    consistency_comparison = []
    for row in result["tomographic_delta_sigma_consistency"]:
        base = baseline_tests[(float(row["r_min_mpc"]), float(row["r_max_mpc"]))]
        base_index = baseline_test_rows.index(base)
        consistency_comparison.append({
            "r_min_mpc": row["r_min_mpc"],
            "r_max_mpc": row["r_max_mpc"],
            "baseline_chi2": base["chi2"],
            "baseline_p_value": base["p_value"],
            "baseline_holm_p_value_7_radial_tests": float(baseline_adjusted[base_index]),
            "shifted_chi2": row["chi2"],
            "shifted_p_value": row["p_value"],
            "shifted_holm_p_value_7_radial_tests": row["holm_p_value_7_radial_tests"],
        })
    return {
        "schema": "astrafek-hsc-y3-published-nz-shift-sensitivity-v1",
        "status": "approximate_sensitivity",
        "science_measurement": False,
        "scenario": scenario_label,
        "method": "translate only local SACC HSC Y3 source-bin 3 and 4 curves by a published central or marginal-interval endpoint shift toward higher redshift, renormalize on the SACC grid, and rerun the same four-bin joint spatial jackknife",
        "published_shift_source": "Choppin de Janvry et al. 2026, arXiv:2511.18133v3; paper-convention Delta z_3=-0.043(+0.022/-0.020), Delta z_4=-0.050(+0.012/-0.012); negative Delta z is described as a shift toward higher redshift",
        "applied_mean_redshift_shifts": shifts,
        "bins_1_and_2_treatment": "unchanged; bin 1's small lower-z shift is not quantified in this approximation, while bin 2 is reported consistent with photometric calibration",
        "shifted_tomography": result,
        "baseline_vs_shifted_consistency": consistency_comparison,
        "baseline_vs_shifted_profiles": compare_profiles(baseline["profile"], result["profile"]),
        "limitations": [
            "A constant translation of two fiducial SACC curves is an approximation, not the authors' calibrated spline distributions or posterior samples.",
            "The five scenarios span the reported marginal 68% endpoints as a rectangular sensitivity grid; no joint posterior or covariance between the two shifts is available here.",
            "The bin-1 low-redshift correction is not included.",
            "This remains a cluster-only stack with no random-lens subtraction or connected-footprint covariance.",
            "No DESI distribution files, new catalog rows, or images are downloaded by this analysis.",
        ],
        "input": {
            **provenance,
            "baseline_row_plan_sha256": baseline_input.get("row_plan_sha256"),
            "n_z_sacc_sha256": sha(Path(baseline_input["n_z_sacc_path"])),
            "n_z_sacc_path": baseline_input["n_z_sacc_path"],
        },
    }


def run_published_shift_grid(pairs, provenance, nz_redshift, nz_by_bin, baseline,
                             radial_edges_mpc, cosmology, *, cell_width_deg):
    scenarios = {}
    for label, shifts in PUBLISHED_SHIFT_SCENARIOS.items():
        scenarios[label] = run_sensitivity(
            pairs, provenance, nz_redshift, nz_by_bin, baseline,
            radial_edges_mpc, cosmology, cell_width_deg=cell_width_deg,
            mean_shifts=shifts, scenario_label=label,
        )
    radii = [(float(row["r_min_mpc"]), float(row["r_max_mpc"]))
             for row in next(iter(scenarios.values()))["baseline_vs_shifted_consistency"]]
    consistency_envelope = []
    for r_min, r_max in radii:
        rows = [next(row for row in scenario["baseline_vs_shifted_consistency"]
                     if (float(row["r_min_mpc"]), float(row["r_max_mpc"])) == (r_min, r_max))
                for scenario in scenarios.values()]
        consistency_envelope.append({
            "r_min_mpc": r_min,
            "r_max_mpc": r_max,
            "min_shifted_chi2": min(float(row["shifted_chi2"]) for row in rows),
            "max_shifted_chi2": max(float(row["shifted_chi2"]) for row in rows),
            "min_shifted_holm_p_value": min(float(row["shifted_holm_p_value_7_radial_tests"]) for row in rows),
            "max_shifted_holm_p_value": max(float(row["shifted_holm_p_value_7_radial_tests"]) for row in rows),
        })
    return {
        "schema": "astrafek-hsc-y3-published-nz-shift-uncertainty-grid-v1",
        "status": "approximate_sensitivity_grid",
        "science_measurement": False,
        "shift_intervals_higher_redshift": {
            str(k): list(v) for k, v in PUBLISHED_HIGH_Z_SHIFT_INTERVALS.items()
        },
        "grid_definition": "central point plus all four corners of the rectangular marginal 68% shift intervals for bins 3 and 4",
        "scenarios": scenarios,
        "tomographic_consistency_envelope": consistency_envelope,
        "input": scenarios["central"]["input"],
        "limitations": [
            "Each curve is a constant translation of a local SACC fiducial n(z), not the authors' calibrated spline distribution or posterior samples.",
            "The marginal endpoint rectangle is not a joint credible region; covariance between bin shifts is not propagated.",
            "This remains a cluster-only stack with no random-lens subtraction or connected-footprint covariance.",
            "No DESI distribution files, new catalog rows, or images are downloaded by this analysis.",
        ],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rows", type=Path, required=True)
    parser.add_argument("--preflight", type=Path, required=True)
    parser.add_argument("--plan-dir", type=Path, required=True)
    parser.add_argument("--nz-sacc", type=Path, required=True)
    parser.add_argument("--baseline-summary", type=Path, required=True)
    parser.add_argument("--masked-value", type=int, choices=(0, 1), required=True)
    parser.add_argument("--cell-width-deg", type=float, default=5.0)
    parser.add_argument("--radial-edges-mpc", type=float, nargs="+", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--h0", type=float, default=70.0)
    parser.add_argument("--omega-m", type=float, default=0.3)
    args = parser.parse_args()
    if args.output.exists():
        parser.error(f"refusing to overwrite existing output directory: {args.output}")
    pairs, provenance = verify_and_load_pairs(
        args.rows, args.preflight, args.plan_dir, args.masked_value,
    )
    baseline = json.loads(args.baseline_summary.read_text(encoding="utf-8"))
    baseline_input = baseline.get("input", {})
    if baseline_input.get("n_z_sacc_sha256") != sha(args.nz_sacc):
        parser.error("baseline summary and local SACC n(z) file have different checksums")
    if baseline_input.get("radial_edges_mpc") != [float(value) for value in args.radial_edges_mpc]:
        parser.error("baseline summary and requested radial edges differ")
    if baseline_input.get("cosmology") != {"H0_km_s_mpc": args.h0, "Omega_m": args.omega_m}:
        parser.error("baseline summary and requested cosmology differ")
    nz_redshift, nz_by_bin = load_hsc_y3_nz_sacc(args.nz_sacc)
    result = run_published_shift_grid(
        pairs, provenance, nz_redshift, nz_by_bin, baseline,
        args.radial_edges_mpc, Cosmology(H0=args.h0, Om=args.omega_m),
        cell_width_deg=args.cell_width_deg,
    )
    result["input"].update({
        "runner_sha256": sha(Path(__file__).resolve()),
        "baseline_summary_path": str(args.baseline_summary.resolve()),
        "baseline_summary_sha256": sha(args.baseline_summary),
        "cosmology": {"H0_km_s_mpc": args.h0, "Omega_m": args.omega_m},
        "radial_edges_mpc": [float(value) for value in args.radial_edges_mpc],
        "cell_width_deg": args.cell_width_deg,
    })
    args.output.mkdir(parents=True, exist_ok=False)
    path = args.output / "summary.json"
    with path.open("x", encoding="utf-8") as stream:
        json.dump(result, stream, indent=2, ensure_ascii=False, allow_nan=False)
        stream.write("\n")
    print(json.dumps({
        "status": result["status"],
        "science_measurement": result["science_measurement"],
        "scenarios": {key: value["applied_mean_redshift_shifts"]
                      for key, value in result["scenarios"].items()},
        "consistency_envelope": result["tomographic_consistency_envelope"],
        "output": str(path),
    }, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
