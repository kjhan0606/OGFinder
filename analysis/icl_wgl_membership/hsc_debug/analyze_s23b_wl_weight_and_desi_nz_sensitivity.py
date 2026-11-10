#!/usr/bin/env python3
"""Compare equal-lens random weighting and existing S19A DESI n(z) curves.

This offline sensitivity runner reuses only the approved 80-lens/400-control
HSC shape-row export, the pinned SACC baseline, and the four previously
approved DESI/HSC n(z) artifacts. It issues no query, downloads no data, reads
no image pixels, and records science_measurement=false.
"""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import analyze_s23b_wl_random_subtracted_profile as base
from audit_wl_desi_nz_distribution_sensitivity import load_desi_curves
from fetch_ancillary_catalogues import sha
from ogfkit.wlcluster import Cosmology
from s19a_shape_stack import load_hsc_y3_nz_sacc, stack_s19a_tomographic_pairs


DATA = base.DATA
DEFAULT_DESI = DATA / "desi_y3_hsc_nz_abc70af2_2026-10-20"
DEFAULT_BASELINE = DATA / "wl_s23b_random_subtracted_profile_2026-11-11_01/summary.json"
DEFAULT_OUTPUT = DATA / "wl_s23b_weight_and_desi_nz_sensitivity_2026-11-11_03"
EXPECTED_BASELINE_SHA256 = "2c899b0090dd94a34f1fd6b0cb111e3d1ada696da2178f98e9fd5b659e699edf"
EXPECTED_BASELINE_RUNNER_SHA256 = "7e7ca33190f59fb9daeee91994231127c58facfa9d452554758a43fccf046c65"
SOURCE_BINS = {2, 3, 4}
PRIMARY_CASE = "controls_ge4_bins_2_3_4_5deg"
CURVE_KINDS = ("mean", "mean_mu", "lower", "upper", "lower_mu", "upper_mu")


def read_json(path: Path) -> dict:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def write_json(path: Path, payload: dict) -> None:
    Path(path).write_text(
        json.dumps(payload, ensure_ascii=False, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with Path(path).open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def stack_difference(lens_pairs: dict, control_pairs: dict, nz_z: np.ndarray,
                     nz_by_bin: np.ndarray, cosmology: Cosmology) -> dict:
    lens = stack_s19a_tomographic_pairs(
        lens_pairs, nz_z, nz_by_bin, base.RADIAL_EDGES_MPC, cosmology,
        apply_b_mode_mask=True,
    )
    control = stack_s19a_tomographic_pairs(
        control_pairs, nz_z, nz_by_bin, base.RADIAL_EDGES_MPC, cosmology,
        apply_b_mode_mask=True,
    )
    rows = []
    for index, (lens_row, control_row) in enumerate(zip(lens["profile"], control["profile"])):
        if lens_row["bin"] != control_row["bin"]:
            raise ValueError("lens/control stack radial-bin mismatch")
        rows.append({
            "bin": index,
            "r_min_mpc": lens_row["r_min_mpc"],
            "r_max_mpc": lens_row["r_max_mpc"],
            "delta_sigma_random_subtracted_msun_mpc2": (
                float(lens_row["delta_sigma_weak_msun_mpc2"])
                - float(control_row["delta_sigma_weak_msun_mpc2"])
                if lens_row["delta_sigma_weak_msun_mpc2"] is not None
                and control_row["delta_sigma_weak_msun_mpc2"] is not None else None
            ),
            "cross_delta_sigma_random_subtracted_msun_mpc2": (
                float(lens_row["cross_delta_sigma_weak_msun_mpc2"])
                - float(control_row["cross_delta_sigma_weak_msun_mpc2"])
                if lens_row["cross_delta_sigma_weak_msun_mpc2"] is not None
                and control_row["cross_delta_sigma_weak_msun_mpc2"] is not None else None
            ),
            "delta_sigma_jackknife_error_msun_mpc2": None,
            "cross_jackknife_error_msun_mpc2": None,
            "n_lens_pairs": int(lens_row["n_pairs"]),
            "n_random_pairs": int(control_row["n_pairs"]),
        })
    return {"profile": rows, "lens_stack": lens, "control_stack": control}


def equal_lens_weight_sensitivity(lens_pairs: dict, control_pairs: dict,
                                  lens_ids: set[str], control_ids: set[str],
                                  center_to_lens: dict[str, str], nz_z: np.ndarray,
                                  nz_by_bin: np.ndarray, cosmology: Cosmology,
                                  baseline_profile: list[dict]) -> tuple[dict, list[dict]]:
    """Compare pooled and equal-paired estimates on identical radial support."""
    selected_lens = base.subset_pairs(lens_pairs, lens_ids, SOURCE_BINS)
    selected_control = base.subset_pairs(control_pairs, control_ids, SOURCE_BINS)
    full_pooled = stack_difference(selected_lens, selected_control, nz_z, nz_by_bin, cosmology)
    for observed, expected in zip(full_pooled["profile"], baseline_profile):
        delta = observed["delta_sigma_random_subtracted_msun_mpc2"]
        target = float(expected["delta_sigma_random_subtracted_msun_mpc2"])
        if delta is None or not np.isclose(delta, target, rtol=1e-10, atol=1e-5):
            raise ValueError("recomputed SACC pooled profile differs from the pinned baseline")

    lens_center_ids = np.asarray(selected_lens["lens_id"]).astype(str)
    center_ids = np.asarray(selected_control["lens_id"]).astype(str)
    matched_lens = np.asarray([center_to_lens[center_id] for center_id in center_ids], dtype=object)
    group_profiles = []
    for lens_id in sorted(lens_ids):
        lens_indices = np.flatnonzero(lens_center_ids == lens_id)
        control_indices = np.flatnonzero(matched_lens == lens_id)
        lens_result = None
        control_result = None
        if len(lens_indices):
            lens_subset = {key: np.asarray(value)[lens_indices] for key, value in selected_lens.items()}
            lens_result = stack_s19a_tomographic_pairs(
                lens_subset, nz_z, nz_by_bin, base.RADIAL_EDGES_MPC, cosmology,
                apply_b_mode_mask=True,
            )
        if len(control_indices):
            control_subset = {key: np.asarray(value)[control_indices] for key, value in selected_control.items()}
            control_result = stack_s19a_tomographic_pairs(
                control_subset, nz_z, nz_by_bin, base.RADIAL_EDGES_MPC, cosmology,
                apply_b_mode_mask=True,
            )
        n_controls = (len(set(np.asarray(control_subset["lens_id"]).astype(str)))
                      if control_result is not None else 0)
        for index, (lens_row, control_row) in enumerate(zip(
                lens_result["profile"] if lens_result is not None else [None] * (len(base.RADIAL_EDGES_MPC) - 1),
                control_result["profile"] if control_result is not None else [None] * (len(base.RADIAL_EDGES_MPC) - 1))):
            group_profiles.append({
                "matched_lens_id": lens_id,
                "n_control_centers_with_selected_rows": n_controls,
                "radial_bin": index,
                "r_min_mpc": float(base.RADIAL_EDGES_MPC[index]),
                "r_max_mpc": float(base.RADIAL_EDGES_MPC[index + 1]),
                "lens_delta_sigma_msun_mpc2": lens_row["delta_sigma_weak_msun_mpc2"] if lens_row else None,
                "lens_cross_msun_mpc2": lens_row["cross_delta_sigma_weak_msun_mpc2"] if lens_row else None,
                "random_delta_sigma_msun_mpc2": control_row["delta_sigma_weak_msun_mpc2"] if control_row else None,
                "random_cross_msun_mpc2": control_row["cross_delta_sigma_weak_msun_mpc2"] if control_row else None,
                "n_lens_pairs": int(lens_row["n_pairs"]) if lens_row else 0,
                "n_control_pairs": int(control_row["n_pairs"]) if control_row else 0,
            })

    rows = []
    for index, base_row in enumerate(baseline_profile):
        group_rows = [row for row in group_profiles if row["radial_bin"] == index]
        common = [row for row in group_rows
                  if row["lens_delta_sigma_msun_mpc2"] is not None
                  and row["random_delta_sigma_msun_mpc2"] is not None
                  and row["lens_cross_msun_mpc2"] is not None
                  and row["random_cross_msun_mpc2"] is not None]
        if not common:
            raise ValueError(f"no matched lens sets have common target/control support in radial bin {index}")
        common_ids = {row["matched_lens_id"] for row in common}
        common_control_ids = {
            center_id for center_id in control_ids
            if center_to_lens[center_id] in common_ids
        }
        common_lens_pairs = base.subset_pairs(lens_pairs, common_ids, SOURCE_BINS)
        common_control_pairs = base.subset_pairs(control_pairs, common_control_ids, SOURCE_BINS)
        common_pooled = stack_difference(
            common_lens_pairs, common_control_pairs, nz_z, nz_by_bin, cosmology,
        )
        pooled_row = common_pooled["profile"][index]
        pooled_lens_row = common_pooled["lens_stack"]["profile"][index]
        pooled_control_row = common_pooled["control_stack"]["profile"][index]
        if any(row[key] is None for row, key in (
                (pooled_row, "delta_sigma_random_subtracted_msun_mpc2"),
                (pooled_row, "cross_delta_sigma_random_subtracted_msun_mpc2"),
                (pooled_lens_row, "delta_sigma_weak_msun_mpc2"),
                (pooled_lens_row, "cross_delta_sigma_weak_msun_mpc2"),
                (pooled_control_row, "delta_sigma_weak_msun_mpc2"),
                (pooled_control_row, "cross_delta_sigma_weak_msun_mpc2"))):
            raise ValueError(f"common-support pooled stack is empty in radial bin {index}")
        random_t = np.asarray([float(row["random_delta_sigma_msun_mpc2"]) for row in common])
        random_x = np.asarray([float(row["random_cross_msun_mpc2"]) for row in common])
        paired_t = np.asarray([
            float(row["lens_delta_sigma_msun_mpc2"]) - float(row["random_delta_sigma_msun_mpc2"])
            for row in common
        ])
        paired_x = np.asarray([
            float(row["lens_cross_msun_mpc2"]) - float(row["random_cross_msun_mpc2"])
            for row in common
        ])
        pooled_t = float(pooled_row["delta_sigma_random_subtracted_msun_mpc2"])
        pooled_x = float(pooled_row["cross_delta_sigma_random_subtracted_msun_mpc2"])
        equal_random_t = float(pooled_lens_row["delta_sigma_weak_msun_mpc2"] - np.mean(random_t))
        equal_random_x = float(pooled_lens_row["cross_delta_sigma_weak_msun_mpc2"] - np.mean(random_x))
        equal_pair_t, equal_pair_x = float(np.mean(paired_t)), float(np.mean(paired_x))
        rows.append({
            "radial_bin": index,
            "r_min_mpc": float(pooled_row["r_min_mpc"]),
            "r_max_mpc": float(pooled_row["r_max_mpc"]),
            "n_matched_lens_sets": len(common_ids),
            "common_lens_ids_sha256": hashlib.sha256("\n".join(sorted(common_ids)).encode()).hexdigest(),
            "pooled_delta_sigma_msun_mpc2": pooled_t,
            "equal_control_weight_delta_sigma_msun_mpc2": equal_random_t,
            "equal_control_weight_delta_sigma_change_msun_mpc2": equal_random_t - pooled_t,
            "equal_control_weight_delta_sigma_fractional_change": (equal_random_t / pooled_t - 1.0) if pooled_t else None,
            "equal_paired_lens_set_delta_sigma_msun_mpc2": equal_pair_t,
            "equal_paired_lens_set_delta_sigma_change_msun_mpc2": equal_pair_t - pooled_t,
            "equal_paired_lens_set_delta_sigma_fractional_change": (equal_pair_t / pooled_t - 1.0) if pooled_t else None,
            "pooled_cross_msun_mpc2": pooled_x,
            "equal_control_weight_cross_msun_mpc2": equal_random_x,
            "equal_paired_lens_set_cross_msun_mpc2": equal_pair_x,
            "equal_weight_random_delta_sigma_q16_msun_mpc2": float(np.percentile(random_t, 16)),
            "equal_weight_random_delta_sigma_median_msun_mpc2": float(np.median(random_t)),
            "equal_weight_random_delta_sigma_q84_msun_mpc2": float(np.percentile(random_t, 84)),
            "equal_weight_paired_delta_sigma_q16_msun_mpc2": float(np.percentile(paired_t, 16)),
            "equal_weight_paired_delta_sigma_median_msun_mpc2": float(np.median(paired_t)),
            "equal_weight_paired_delta_sigma_q84_msun_mpc2": float(np.percentile(paired_t, 84)),
            "n_lens_pairs": int(pooled_row["n_lens_pairs"]),
            "n_control_pairs_pooled": int(pooled_row["n_random_pairs"]),
            "uncertainty": "common-support per-set spread is descriptive only; no jackknife covariance for weighting variants",
        })
    return {"profile": rows, "n_matched_lens_sets": len(lens_ids), "group_profiles": group_profiles}, rows


def compare_to_sacc(candidate: dict, baseline_profile: list[dict],
                    case_label: str, curve_kind: str, uncertainty_method: str) -> list[dict]:
    rows = []
    for current, base_row in zip(candidate["profile"], baseline_profile):
        value = current["delta_sigma_random_subtracted_msun_mpc2"]
        cross = current["cross_delta_sigma_random_subtracted_msun_mpc2"]
        baseline = float(base_row["delta_sigma_random_subtracted_msun_mpc2"])
        baseline_cross = float(base_row["cross_delta_sigma_random_subtracted_msun_mpc2"])
        rows.append({
            "case_label": case_label,
            "curve_kind": curve_kind,
            "radial_bin": int(current["bin"]),
            "r_min_mpc": float(current["r_min_mpc"]),
            "r_max_mpc": float(current["r_max_mpc"]),
            "sacc_delta_sigma_msun_mpc2": baseline,
            "candidate_delta_sigma_msun_mpc2": value,
            "delta_sigma_difference_msun_mpc2": (float(value) - baseline) if value is not None else None,
            "delta_sigma_fractional_change": (float(value) / baseline - 1.0) if value is not None and baseline else None,
            "candidate_delta_sigma_error_msun_mpc2": current.get("delta_sigma_jackknife_error_msun_mpc2"),
            "sacc_cross_msun_mpc2": baseline_cross,
            "candidate_cross_msun_mpc2": cross,
            "cross_difference_msun_mpc2": (float(cross) - baseline_cross) if cross is not None else None,
            "n_lens_pairs": int(current["n_lens_pairs"]),
            "n_control_pairs": int(current["n_random_pairs"]),
            "uncertainty_method": uncertainty_method,
        })
    return rows


def make_plot(path: Path, baseline_profile: list[dict], equal_rows: list[dict],
              cases: dict) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FuncFormatter

    radius = np.sqrt(base.RADIAL_EDGES_MPC[:-1] * base.RADIAL_EDGES_MPC[1:])
    baseline_t = np.asarray([row["delta_sigma_random_subtracted_msun_mpc2"] for row in baseline_profile]) / 1e13
    baseline_err = np.asarray([row["delta_sigma_jackknife_error_msun_mpc2"] for row in baseline_profile]) / 1e13
    figure, axes = plt.subplots(1, 2, figsize=(13.2, 5.6), sharex=True)
    axes[0].errorbar(radius, baseline_t, yerr=baseline_err, marker="o", color="#1769aa",
                     capsize=2, lw=1.4, label="SACC pooled; full sample")
    axes[0].plot(radius,
                 np.asarray([row["pooled_delta_sigma_msun_mpc2"] for row in equal_rows]) / 1e13,
                 marker="s", color="#6c757d", lw=1.1, label="Pooled; radius-matched support")
    axes[0].plot(radius,
                 np.asarray([row["equal_control_weight_delta_sigma_msun_mpc2"] for row in equal_rows]) / 1e13,
                 marker="^", color="#e07a1f", lw=1.3, label="Equal-control; pooled lens")
    axes[0].plot(radius,
                 np.asarray([row["equal_paired_lens_set_delta_sigma_msun_mpc2"] for row in equal_rows]) / 1e13,
                 marker="v", color="#7950a1", lw=1.3, label="Equal lens-set target−control")
    for variant, label, color, marker in (
        ("0.3_3_v_1p1", "DESI 0.3–3 mean", "#2a9d8f", "^"),
        ("1_5_v_1p1", "DESI 1–5 mean", "#a23b72", "D"),
    ):
        case = cases[f"{variant}__mean"]
        profile = case["profile"]
        values = np.asarray([row["delta_sigma_random_subtracted_msun_mpc2"] for row in profile]) / 1e13
        errors = np.asarray([row["delta_sigma_jackknife_error_msun_mpc2"] for row in profile]) / 1e13
        axes[0].errorbar(radius, values, yerr=errors, marker=marker, color=color,
                         capsize=2, lw=1.2, label=label)
    axes[0].axhline(0.0, color="0.3", lw=0.8)
    axes[0].set_ylabel(r"$\Delta\Sigma_t$ [$10^{13}\,M_\odot\,\mathrm{Mpc}^{-2}$]")
    axes[0].text(0.03, 0.96, "Random-weighting and n(z) profiles", transform=axes[0].transAxes,
                 va="top", ha="left", fontsize=11)
    axes[0].legend(loc="upper right", frameon=False, fontsize=8)

    for variant, kind, label, color, marker in (
        ("0.3_3_v_1p1", "mean", "0.3–3 mean", "#2a9d8f", "o"),
        ("0.3_3_v_1p1", "mean_mu", "0.3–3 mean μ-corrected", "#2a9d8f", "s"),
        ("1_5_v_1p1", "mean", "1–5 mean", "#a23b72", "o"),
        ("1_5_v_1p1", "mean_mu", "1–5 mean μ-corrected", "#a23b72", "s"),
    ):
        rows = cases[f"{variant}__{kind}"]["vs_sacc_baseline"]
        pct = np.asarray([row["delta_sigma_fractional_change"] for row in rows]) * 100.0
        axes[1].plot(radius, pct, marker=marker, color=color, lw=1.2, label=label)
    envelope_values = np.asarray([
        [row["delta_sigma_fractional_change"] * 100.0 for row in cases[key]["vs_sacc_baseline"]]
        for key in cases if cases[key]["is_pointwise_marginal_envelope"]
    ])
    axes[1].fill_between(radius, np.min(envelope_values, axis=0), np.max(envelope_values, axis=0),
                         color="#777777", alpha=0.18,
                         label="DESI pointwise marginal envelope range")
    axes[1].axhline(0.0, color="0.3", lw=0.8)
    axes[1].set_ylabel(r"Tangential $\Delta\Sigma_t$ change vs SACC [%]")
    axes[1].text(0.03, 0.96, "Calibration sensitivity", transform=axes[1].transAxes,
                 va="top", ha="left", fontsize=11)
    axes[1].legend(loc="upper right", frameon=False, fontsize=7.5)

    for axis in axes:
        axis.set_xscale("log")
        axis.set_xticks([0.1, 0.2, 0.4, 0.7, 1.0, 1.5, 2.0, 3.0])
        axis.get_xaxis().set_major_formatter(FuncFormatter(lambda value, _pos: f"{value:g}"))
        axis.set_xlabel("Projected radius R [proper Mpc]")
        axis.grid(True, which="major", alpha=0.22)
    figure.text(0.5, 0.015,
                "Exploratory S19A shape stack; marginal envelope is not a joint confidence band.",
                ha="center", fontsize=9)
    figure.tight_layout(rect=(0, 0.05, 1, 1))
    figure.savefig(path, dpi=180)
    plt.close(figure)


def write_assessment(path: Path, summary: dict) -> None:
    equal_rows = summary["equal_lens_weight"]["profile"]
    n_common_min = min(row["n_matched_lens_sets"] for row in equal_rows)
    n_common_max = max(row["n_matched_lens_sets"] for row in equal_rows)
    case_rows = []
    for label, case in summary["desi_cases"].items():
        pct = np.asarray([
            row["delta_sigma_fractional_change"] * 100.0
            for row in case["vs_sacc_baseline"]
            if row["delta_sigma_fractional_change"] is not None
        ])
        case_rows.append((label, case, pct))
    rows = [
        "# S23B random-weight and S19A DESI n(z) sensitivity (2026-11-11)",
        "",
        "This offline Astrafek analysis reused only the approved 80-lens/400-control HSC S19A shape-row export, the existing HSC Y3 SACC baseline, and the four already approved DESI/HSC n(z) products. It made no CAS query, network request, data download, or image read. All results are exploratory and `science_measurement=false`.",
        "",
        "## Equal matched-lens-set weighting",
        "",
        f"The reference pools lens/control rows over supported matched sets. Comparisons at each radius use the same lens sets with both target and control support; this subset varies from {n_common_min} to {n_common_max} sets by radius. The equal-control estimate changes only the control side by averaging per-set random profiles; the equal-paired estimate averages each lens-set lens-minus-control profile. The common-support pooled estimate is the reference. Per-set quantiles describe heterogeneity, not uncertainty. Values are in 10^13 Msun Mpc^-2.",
        "",
        "| R proper (Mpc) | Common pooled ΔΣt | Equal-control ΔΣt | Equal-paired ΔΣt | Change (%) | Sets |",
        "|---:|---:|---:|---:|---:|---:|",
    ]
    for row in equal_rows:
        rows.append(
            f"| {row['r_min_mpc']:.1f}–{row['r_max_mpc']:.1f} | "
            f"{row['pooled_delta_sigma_msun_mpc2']/1e13:.4f} | "
            f"{row['equal_control_weight_delta_sigma_msun_mpc2']/1e13:.4f} | "
            f"{row['equal_paired_lens_set_delta_sigma_msun_mpc2']/1e13:.4f} | "
            f"{100*row['equal_paired_lens_set_delta_sigma_fractional_change']:+.2f} | "
            f"{row['n_matched_lens_sets']} |"
        )
    rows += [
        "",
        "## DESI/HSC n(z) comparison with the SACC baseline",
        "",
        "Each case uses the same 74-lens/367-control primary sample and source bins 2–4. Central mean and magnification-corrected curves are rerun with the same 5-degree paired endpoint jackknife. Published lower/upper curves are marginal envelopes only and have no paired covariance; their range is not a joint confidence interval.",
        "",
        "| DESI scale cut and curve | Median ΔΣt change across radii (%) | Minimum–maximum across radii (%) | Covariance |",
        "|---|---:|---:|---|",
    ]
    for label, case, pct in case_rows:
        median = float(np.median(pct)) if len(pct) else float("nan")
        minmax = f"{np.min(pct):+.2f} to {np.max(pct):+.2f}" if len(pct) else "n/a"
        covariance = "5° endpoint jackknife" if case["is_central_estimate"] else "not estimated; marginal envelope"
        rows.append(f"| {label} | {median:+.2f} | {minmax} | {covariance} |")
    rows += [
        "",
        "The central-curve changes are calibration sensitivities, not a posterior uncertainty measurement. The published products contain pointwise lower/upper curves but no joint posterior draws. The DESI products are S19A-specific and improve the source-selection match relative to treating the HSC Y3 SACC curve as exact, but they do not remove photo-z dilution uncertainty or the cross-shear issue seen in the prior profile.",
        "",
        "## Interpretation limits",
        "",
        "The equal-set profiles are weighting sensitivities and have no survey-region covariance of their own. The 5-degree cells remain approximate and can combine disconnected HSC fields. Neither this sensitivity nor the previous profile propagates a full n(z) posterior, reduced shear, magnification uncertainty, cluster-member dilution/boost, random-aperture depth/mask/PSF/sky mismatch, or selection-response uncertainty. No halo mass was fit; keep `science_measurement=false`.",
        "",
        "## Reproducibility",
        "",
        f"- Runner SHA-256: `{summary['provenance']['runner_sha256']}`",
        f"- Baseline profile summary SHA-256: `{summary['provenance']['baseline_summary_sha256']}`",
        f"- Approved source-row plan SHA-256: `{summary['provenance']['row_plan_sha256']}`",
        f"- DESI release commit: `{summary['provenance']['desi_commit']}`",
        f"- Output directory: `{summary['output_directory']}`",
        "- Profile and sensitivity table: `profile_sensitivity.tsv` and `equal_lens_control_profiles.tsv`",
        "",
    ]
    path.write_text("\n".join(rows), encoding="utf-8")


def run(args) -> dict:
    output = args.output.resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite sensitivity output: {output}")
    if sha(args.baseline_summary) != EXPECTED_BASELINE_SHA256:
        raise ValueError("pinned S23B baseline profile summary changed")
    baseline = read_json(args.baseline_summary)
    baseline_provenance = baseline.get("provenance", {})
    if (baseline.get("science_measurement") is not False
            or baseline.get("network_used") is not False
            or baseline_provenance.get("runner_sha256") != EXPECTED_BASELINE_RUNNER_SHA256):
        raise ValueError("baseline profile has unexpected provenance or science state")
    baseline_profile = baseline["profile_cases"][PRIMARY_CASE]["profile"]

    lens_pairs, control_pairs, center_to_lens, row_provenance = base.load_approved_pairs(
        args.plan_dir.resolve(), args.rows.resolve(),
    )
    support, support_summary = base.support_policies(args.qa_dir.resolve(), center_to_lens)
    lens_ids, control_ids = support[4]
    if (len(lens_ids) != 74 or len(control_ids) != 367
            or len(baseline_profile) != len(base.RADIAL_EDGES_MPC) - 1):
        raise ValueError("primary matched sample differs from the pinned baseline")

    nz_z, sacc_nz = load_hsc_y3_nz_sacc(args.nz_sacc.resolve())
    desi_z, desi_curves, desi_hashes = load_desi_curves(args.desi_directory.resolve())
    cosmology = Cosmology(H0=70.0, Om=0.3)

    equal_lens, equal_set_rows = equal_lens_weight_sensitivity(
        lens_pairs, control_pairs, lens_ids, control_ids, center_to_lens,
        nz_z, sacc_nz, cosmology, baseline_profile,
    )
    selected_lens = base.subset_pairs(lens_pairs, lens_ids, SOURCE_BINS)
    selected_control = base.subset_pairs(control_pairs, control_ids, SOURCE_BINS)
    desi_cases, comparison_rows = {}, []
    for variant in sorted(desi_curves):
        for kind in CURVE_KINDS:
            label = f"{variant}__{kind}"
            curve = desi_curves[variant][kind]
            central = kind in ("mean", "mean_mu")
            if central:
                result = base.run_profile_case(
                    lens_pairs, control_pairs, lens_ids, control_ids, SOURCE_BINS,
                    desi_z, curve, cosmology, 5.0, label,
                )
                uncertainty = "paired 5-degree endpoint jackknife"
            else:
                result = stack_difference(selected_lens, selected_control, desi_z, curve, cosmology)
                uncertainty = "not estimated; published pointwise marginal envelope"
                result.update({
                    "science_measurement": False,
                    "case_label": label,
                    "source_bins": sorted(SOURCE_BINS),
                    "n_selected_lenses": len(lens_ids),
                    "n_selected_control_centers": len(control_ids),
                })
            comparisons = compare_to_sacc(result, baseline_profile, label, kind, uncertainty)
            case_record = {
                "case_label": label,
                "scale_cut_mpc_over_h": [0.3, 3.0] if variant == "0.3_3_v_1p1" else [1.0, 5.0],
                "curve_kind": kind,
                "is_central_estimate": central,
                "is_magnification_corrected": kind.endswith("_mu"),
                "is_pointwise_marginal_envelope": kind.startswith(("lower", "upper")),
                "uncertainty_method": uncertainty,
                **result,
                "vs_sacc_baseline": comparisons,
            }
            desi_cases[label] = case_record
            comparison_rows.extend(comparisons)
            print(f"completed n(z) case {label}", flush=True)

    output.mkdir(parents=True, exist_ok=False)
    equal_tsv_rows = []
    for group_row in equal_lens["group_profiles"]:
        equal_tsv_rows.append(group_row)
    write_tsv(
        output / "equal_lens_control_profiles.tsv", equal_tsv_rows,
        ["matched_lens_id", "n_control_centers_with_selected_rows", "radial_bin",
         "r_min_mpc", "r_max_mpc", "lens_delta_sigma_msun_mpc2",
         "lens_cross_msun_mpc2", "random_delta_sigma_msun_mpc2",
         "random_cross_msun_mpc2", "n_lens_pairs", "n_control_pairs"],
    )
    write_tsv(
        output / "profile_sensitivity.tsv", comparison_rows,
        ["case_label", "curve_kind", "radial_bin", "r_min_mpc", "r_max_mpc",
         "sacc_delta_sigma_msun_mpc2", "candidate_delta_sigma_msun_mpc2",
         "delta_sigma_difference_msun_mpc2", "delta_sigma_fractional_change",
         "candidate_delta_sigma_error_msun_mpc2", "sacc_cross_msun_mpc2",
         "candidate_cross_msun_mpc2", "cross_difference_msun_mpc2",
         "n_lens_pairs", "n_control_pairs", "uncertainty_method"],
    )
    summary = {
        "schema": "astrafek-s23b-wl-weight-and-desi-nz-sensitivity-v1",
        "status": "complete_exploratory_sensitivity",
        "science_measurement": False,
        "network_used": False,
        "cas_query_submitted": False,
        "new_data_downloaded": False,
        "images_read": False,
        "halo_mass_fit_performed": False,
        "scope": "equal matched-lens random weighting and existing S19A DESI/HSC n(z) curves for the 80-lens pilot",
        "sample": {
            "n_input_lenses": 80,
            "n_input_control_centers": 400,
            "n_selected_lens_sets": len(lens_ids),
            "n_selected_control_centers": len(control_ids),
            "source_bins": sorted(SOURCE_BINS),
            "radial_edges_proper_mpc": base.RADIAL_EDGES_MPC.tolist(),
            "support_policy": support_summary["policies"]["4"],
        },
        "equal_lens_weight": equal_lens,
        "desi_cases": desi_cases,
        "limitations": [
            "published DESI lower/upper n(z) arrays are pointwise marginal envelopes, not joint posterior draws",
            "equal-lens-set weighting sensitivity has no independent survey-region covariance; per-set scatter is descriptive",
            "5-degree endpoint cells are approximate and may combine disconnected HSC fields",
            "reduced shear, magnification uncertainty, cluster-member dilution/boost, random-aperture depth/mask/PSF/sky matching, and selection-response uncertainty are not fully propagated",
            "the analysis is exploratory and science_measurement remains false",
        ],
        "provenance": {
            "runner_path": str(Path(__file__).resolve()),
            "runner_sha256": sha(Path(__file__).resolve()),
            "baseline_runner_sha256": baseline_provenance["runner_sha256"],
            "baseline_summary_path": str(args.baseline_summary.resolve()),
            "baseline_summary_sha256": sha(args.baseline_summary),
            "row_plan_sha256": row_provenance["plan_sha256"],
            "row_status_sha256": row_provenance["row_status_sha256"],
            "row_summary_sha256": row_provenance["row_summary_sha256"],
            "row_batch_sha256": row_provenance["batch_sha256"],
            "qa_summary_sha256": sha(args.qa_dir.resolve() / "summary.json"),
            "per_center_support_sha256": sha(args.qa_dir.resolve() / "per_center.tsv"),
            "sacc_nz_sha256": sha(args.nz_sacc.resolve()),
            "desi_commit": "abc70af29910e49fce17a63ef3083160d515fb6c",
            "desi_artifact_sha256": desi_hashes,
            "science_measurement": False,
        },
        "output_directory": str(output),
    }
    write_assessment(output / "assessment.md", summary)
    make_plot(output / "sensitivity.png", baseline_profile, equal_lens["profile"], desi_cases)
    summary["output_sha256"] = {
        name: sha(output / name)
        for name in ("assessment.md", "equal_lens_control_profiles.tsv", "profile_sensitivity.tsv", "sensitivity.png")
    }
    write_json(output / "summary.json", summary)
    print(json.dumps({
        "status": summary["status"],
        "science_measurement": False,
        "output_directory": str(output),
        "summary_sha256": sha(output / "summary.json"),
        "equal_lens_weight_profile": equal_lens["profile"],
        "desi_median_percent_change": {
            label: float(np.median([
                row["delta_sigma_fractional_change"] * 100.0
                for row in case["vs_sacc_baseline"]
                if row["delta_sigma_fractional_change"] is not None
            ]))
            for label, case in desi_cases.items()
        },
    }, ensure_ascii=False, indent=2), flush=True)
    return summary


def main() -> int:
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan-dir", type=Path, default=base.DEFAULT_PLAN_DIR)
    parser.add_argument("--rows", type=Path, default=base.DEFAULT_RESULTS)
    parser.add_argument("--qa-dir", type=Path, default=base.DEFAULT_QA)
    parser.add_argument("--nz-sacc", type=Path, default=base.DEFAULT_NZ)
    parser.add_argument("--desi-directory", type=Path, default=DEFAULT_DESI)
    parser.add_argument("--baseline-summary", type=Path, default=DEFAULT_BASELINE)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    run(parser.parse_args())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
