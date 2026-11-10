#!/usr/bin/env python3
"""Propagate HSC S19A edge-selection coefficient errors one at a time.

The runner preserves the pinned 52-lens/260-control sample and recomputes the
random-subtracted profile after perturbing each released selection-bias
coefficient by its quoted 1-sigma error. It does not combine multiplicative
and additive calibration uncertainties into a joint covariance.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np

import analyze_s23b_wl_fixed_common_support_spatial_covariance as base
from fetch_ancillary_catalogues import sha
from ogfkit.wlcluster import Cosmology
from s19a_shape_stack import (
    HSC_SELECTION_BIAS_COEFFICIENTS,
    HSC_SELECTION_BIAS_COEFFICIENT_ERRORS,
    load_hsc_y3_nz_sacc,
)
from submit_s23b_wl_common_shape_rows import DEFAULT_PLAN_DIR, DEFAULT_RESULTS


DATA = base.DATA
DEFAULT_QA_DIR = base.QA_DIR
DEFAULT_SUPPORT_DIR = base.SUPPORT_DIR
DEFAULT_PROFILE_DIR = base.FIXED_PROFILE_DIR
DEFAULT_NZ = base.NZ_PATH
DEFAULT_OUTPUT = DATA / "wl_s23b_fixed52_selection_coefficient_sensitivity_2026-11-11_01"
COEFFICIENTS = dict(HSC_SELECTION_BIAS_COEFFICIENTS)
COEFFICIENT_ERRORS = dict(HSC_SELECTION_BIAS_COEFFICIENT_ERRORS)


def read_json(path: Path) -> dict:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def write_json(path: Path, payload: dict) -> None:
    with Path(path).open("x", encoding="utf-8") as stream:
        json.dump(payload, stream, indent=2, ensure_ascii=False, allow_nan=False)
        stream.write("\n")


def write_tsv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty table: {path.name}")
    with Path(path).open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(
            stream, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def evaluate_profile(lens_pairs: dict, control_pairs: dict, lens_ids: set[str],
                     controls_by_lens: dict[str, list[str]], lens_cache: dict,
                     control_cache: dict, lens_rows_by_set: dict[str, np.ndarray],
                     control_rows_by_center: dict[str, np.ndarray],
                     cosmology: Cosmology, coefficients: dict[str, float]) -> np.ndarray:
    per_set = []
    for lens_id in sorted(lens_ids):
        target_indices = lens_rows_by_set.get(lens_id, np.asarray([], dtype=np.int64))
        center_indices = [
            control_rows_by_center[control_id]
            for control_id in sorted(controls_by_lens[lens_id])
            if control_id in control_rows_by_center
        ]
        if len(center_indices) < base.MIN_CONTROLS_PER_ANNULUS:
            raise ValueError(f"fixed lens set {lens_id} has too few control centers with source rows")
        random_indices = np.concatenate(center_indices)
        target_result = base.cached_stack(
            lens_pairs, lens_cache, target_indices, cosmology,
            selection_bias_coefficients=coefficients,
        )
        random_result = base.cached_stack(
            control_pairs, control_cache, random_indices, cosmology,
            selection_bias_coefficients=coefficients,
        )
        delta, n_controls = base.profile_values(target_result, random_result)
        if (delta is None or np.any(~np.isfinite(delta))
                or np.any(n_controls < base.MIN_CONTROLS_PER_ANNULUS)):
            raise ValueError(f"selection-coefficient case changed fixed support for {lens_id}")
        per_set.append(delta)
    return np.mean(np.stack(per_set), axis=0)


def run(args) -> dict:
    output = args.output.resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite existing sensitivity output: {output}")
    for label, path, expected in (
        ("QA summary", args.qa_dir / "summary.json", base.EXPECTED_QA_SUMMARY_SHA256),
        ("fixed-profile summary", args.fixed_profile_dir / "summary.json",
         base.EXPECTED_FIXED_PROFILE_SUMMARY_SHA256),
        ("HSC Y3 n(z)", args.nz_sacc, base.EXPECTED_NZ_SHA256),
    ):
        observed = sha(path)
        if observed != expected:
            raise ValueError(f"{label} hash mismatch: {observed}")
    qa_summary = read_json(args.qa_dir / "summary.json")
    fixed_summary = read_json(args.fixed_profile_dir / "summary.json")
    if (qa_summary.get("science_measurement") is not False
            or fixed_summary.get("science_measurement") is not False
            or fixed_summary.get("covariance_estimated") is not False):
        raise ValueError("a pinned input artifact unexpectedly promotes science or contains covariance")

    lens_ids, fixed_selection = base.fixed_sample_ids(args.support_dir.resolve())
    lens_pairs, all_control_pairs, center_to_lens, row_provenance = base.load_approved_pairs(
        args.plan_dir.resolve(), args.rows.resolve(),
    )
    controls_by_lens: dict[str, list[str]] = {lens_id: [] for lens_id in lens_ids}
    for center_id, lens_id in center_to_lens.items():
        if lens_id in lens_ids:
            controls_by_lens[lens_id].append(center_id)
    if any(len(centers) != 5 for centers in controls_by_lens.values()):
        raise ValueError("a fixed lens set no longer has exactly five matched controls")
    selected_control_ids = {
        center_id for centers in controls_by_lens.values() for center_id in centers
    }
    lens_pairs = base.select_rows(lens_pairs, lens_ids)
    control_pairs = base.select_rows(all_control_pairs, selected_control_ids)
    del all_control_pairs

    nz_z, nz_bins = load_hsc_y3_nz_sacc(args.nz_sacc.resolve())
    cosmology = Cosmology(H0=70.0, Om=0.3)
    lens_cache = base.cache_geometry_and_efficiency(lens_pairs, nz_z, nz_bins, cosmology)
    control_cache = base.cache_geometry_and_efficiency(control_pairs, nz_z, nz_bins, cosmology)
    lens_rows_by_set = base.group_indices(lens_pairs["lens_id"])
    control_rows_by_center = base.group_indices(control_pairs["lens_id"])
    baseline = evaluate_profile(
        lens_pairs, control_pairs, lens_ids, controls_by_lens, lens_cache,
        control_cache, lens_rows_by_set, control_rows_by_center, cosmology, COEFFICIENTS,
    )
    reference_rows = fixed_summary.get("radial_profile_msun_mpc2", [])
    if len(reference_rows) != len(base.RADIAL_EDGES_MPC) - 1:
        raise ValueError("pinned fixed-profile summary does not have seven radial bins")
    reference = np.asarray([
        [row["delta_sigma_mean_msun_mpc2"], row["cross_delta_sigma_mean_msun_mpc2"]]
        for row in reference_rows
    ], dtype=float)
    if not np.allclose(baseline, reference, rtol=2e-9, atol=1e-3):
        raise ValueError("recomputed baseline differs from the pinned fixed 52-set profile")

    responses: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
    for name, sigma in COEFFICIENT_ERRORS.items():
        plus = dict(COEFFICIENTS)
        minus = dict(COEFFICIENTS)
        plus[name] += sigma
        minus[name] -= sigma
        plus_profile = evaluate_profile(
            lens_pairs, control_pairs, lens_ids, controls_by_lens, lens_cache,
            control_cache, lens_rows_by_set, control_rows_by_center, cosmology, plus,
        )
        minus_profile = evaluate_profile(
            lens_pairs, control_pairs, lens_ids, controls_by_lens, lens_cache,
            control_cache, lens_rows_by_set, control_rows_by_center, cosmology, minus,
        )
        responses[name] = (plus_profile, minus_profile, (plus_profile - minus_profile) / 2.0)

    rows = []
    for name, (plus_profile, minus_profile, symmetric_shift) in responses.items():
        for index, (lo, hi) in enumerate(zip(
                base.RADIAL_EDGES_MPC[:-1], base.RADIAL_EDGES_MPC[1:])):
            rows.append({
                "coefficient": name,
                "nominal_value": COEFFICIENTS[name],
                "quoted_1sigma": COEFFICIENT_ERRORS[name],
                "radial_bin": index,
                "r_min_mpc": float(lo),
                "r_max_mpc": float(hi),
                "baseline_delta_sigma_msun_mpc2": float(baseline[index, 0]),
                "plus_1sigma_delta_sigma_msun_mpc2": float(plus_profile[index, 0]),
                "minus_1sigma_delta_sigma_msun_mpc2": float(minus_profile[index, 0]),
                "symmetric_1sigma_delta_sigma_response_msun_mpc2": float(symmetric_shift[index, 0]),
                "baseline_cross_msun_mpc2": float(baseline[index, 1]),
                "plus_1sigma_cross_msun_mpc2": float(plus_profile[index, 1]),
                "minus_1sigma_cross_msun_mpc2": float(minus_profile[index, 1]),
                "symmetric_1sigma_cross_response_msun_mpc2": float(symmetric_shift[index, 1]),
            })

    output.mkdir(parents=True, exist_ok=False)
    write_tsv(output / "coefficient_sensitivity.tsv", rows)
    summary = {
        "schema": "astrafek-s23b-wl-selection-coefficient-sensitivity-v1",
        "status": "complete_one_at_a_time_selection_coefficient_sensitivity",
        "science_measurement": False,
        "network_used": False,
        "cas_query_submitted": False,
        "new_data_downloaded": False,
        "images_read": False,
        "halo_mass_fit_performed": False,
        "sample": {
            **fixed_selection,
            "n_control_centers": len(selected_control_ids),
            "lens_source_pair_rows": int(len(lens_pairs["lens_id"])),
            "control_source_pair_rows": int(len(control_pairs["lens_id"])),
            "source_reuse_audit": base.verify_repeated_source_rows(lens_pairs, control_pairs),
            "baseline_matches_pinned_profile": True,
            "baseline_max_absolute_difference_msun_mpc2": float(np.max(np.abs(baseline - reference))),
        },
        "method": {
            "profile": "equal mean over fixed 52 lens sets of target minus pooled matched-control DeltaSigma; recalibrated for each coefficient perturbation",
            "coefficients_perturbed_one_at_a_time": COEFFICIENTS,
            "quoted_coefficient_1sigma_errors": COEFFICIENT_ERRORS,
            "response": "symmetric finite difference of the complete estimator between plus and minus one quoted coefficient sigma",
            "joint_coefficient_covariance": "not estimated; multiplicative/additive cross-covariance is not supplied by the released utility",
            "spatial_covariance": "not estimated by this sensitivity runner",
        },
        "provenance": {
            "runner_path": str(Path(__file__).resolve()),
            "runner_sha256": sha(Path(__file__).resolve()),
            "approved_row_plan_sha256": row_provenance["plan_sha256"],
            "approved_row_status_sha256": row_provenance["row_status_sha256"],
            "approved_row_summary_sha256": row_provenance["row_summary_sha256"],
            "approved_source_row_batch_sha256": row_provenance["batch_sha256"],
            "support_summary_sha256": sha(args.support_dir / "summary.json"),
            "support_rows_sha256": sha(args.support_dir / "per_lens_radial_support.tsv"),
            "qa_summary_sha256": sha(args.qa_dir / "summary.json"),
            "fixed_profile_summary_sha256": sha(args.fixed_profile_dir / "summary.json"),
            "n_z_sacc_sha256": sha(args.nz_sacc),
            "selection_calibration_source": "https://github.com/mr-superonion/hsc-y3-calib/blob/main/utils.py#get_sel_bias",
        },
        "limitations": [
            "reported coefficients are varied one at a time; the output is not a joint systematic covariance or total error budget",
            "the released utility assumes the two edge-cut contributions are independent within the quoted marginal m_sel and a_sel errors, but does not provide a cross-covariance between multiplicative and additive terms",
            "photo-z posterior and dilution uncertainty, reduced shear, magnification, boost/member dilution, random-aperture image matching, PSF uncertainty, and validated contiguous-footprint covariance remain unpropagated",
            "no detection significance or halo mass is inferred",
        ],
        "output_directory": str(output),
        "outputs": {"coefficient_sensitivity_sha256": sha(output / "coefficient_sensitivity.tsv")},
    }
    assessment = [
        "# HSC S19A selection-coefficient sensitivity — fixed S23B 52-lens sample",
        "",
        "This offline sensitivity analysis reuses only the previously approved shape-row export and pinned Astrafek profile/support products. The baseline was required to reproduce the fixed 52-lens profile before any perturbation. No CAS request, data download, or image read was made. `science_measurement=false`.",
        "",
        "Each HSC selection-bias coefficient was shifted independently by its released 1-sigma error in both directions. Every case recomputes the calibrated target and pooled matched-control profiles and averages target-minus-control profiles over the same 52 lens sets.",
        "",
        "The released calibration utility combines the magnitude-edge and resolution-edge contributions in quadrature for the marginal multiplicative and additive errors. It does not publish their cross-covariance. The table therefore reports coefficient-by-coefficient estimator responses; these responses are not a joint systematic covariance and are not added to the spatial uncertainty.",
        "",
        "| Coefficient | Radial bin (proper Mpc) | Baseline ΔΣ | Symmetric 1σ ΔΣ response | Symmetric 1σ cross response |",
        "|---|---:|---:|---:|---:|",
    ]
    assessment.extend(
        f"| {row['coefficient']} | {row['r_min_mpc']:.1f}–{row['r_max_mpc']:.1f} | "
        f"{row['baseline_delta_sigma_msun_mpc2'] / 1e13:+.4f} | "
        f"{row['symmetric_1sigma_delta_sigma_response_msun_mpc2'] / 1e13:+.4f} | "
        f"{row['symmetric_1sigma_cross_response_msun_mpc2'] / 1e13:+.4f} |"
        for row in rows
    )
    assessment.extend([
        "",
        "Keep `science_measurement=false`. Before interpreting the weak-lensing stack, resolve the independent photo-z/dilution, image-matched random, reduced-shear, member-boost, PSF, and spatial-covariance gates.",
        "",
        f"- Runner SHA-256: `{summary['provenance']['runner_sha256']}`",
        f"- Approved row-plan SHA-256: `{row_provenance['plan_sha256']}`",
        f"- Fixed lens-set membership SHA-256: `{fixed_selection['fixed_lens_ids_sha256']}`",
        f"- Output directory: `{output}`",
        "",
    ])
    (output / "assessment.md").write_text("\n".join(assessment), encoding="utf-8")
    summary["outputs"]["assessment_sha256"] = sha(output / "assessment.md")
    write_json(output / "summary.json", summary)
    print(json.dumps({
        "status": summary["status"],
        "science_measurement": False,
        "baseline_matches_pinned_profile": summary["sample"]["baseline_matches_pinned_profile"],
        "coefficient_sensitivities": {
            name: np.max(np.abs(value[2]), axis=0).tolist()
            for name, value in responses.items()
        },
        "output_directory": str(output),
        "summary_sha256": sha(output / "summary.json"),
    }, ensure_ascii=False, indent=2), flush=True)
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan-dir", type=Path, default=DEFAULT_PLAN_DIR)
    parser.add_argument("--rows", type=Path, default=DEFAULT_RESULTS)
    parser.add_argument("--qa-dir", type=Path, default=DEFAULT_QA_DIR)
    parser.add_argument("--support-dir", type=Path, default=DEFAULT_SUPPORT_DIR)
    parser.add_argument("--fixed-profile-dir", type=Path, default=DEFAULT_PROFILE_DIR)
    parser.add_argument("--nz-sacc", type=Path, default=DEFAULT_NZ)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    run(parser.parse_args())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
