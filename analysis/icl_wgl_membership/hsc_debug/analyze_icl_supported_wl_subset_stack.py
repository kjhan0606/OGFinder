#!/usr/bin/env python3
"""Assess the exploratory ICL-supported WL stacks using their full covariance."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import chi2


DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
DEFAULT_STACK = DATA / "icl_supported_wl_subset_stack_2026-11-05_01/summary.json"
DEFAULT_OUTPUT = DATA / "icl_supported_wl_subset_assessment_2026-11-05_01"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def null_test(profile: list[dict], covariance: np.ndarray, offset: int) -> dict:
    indices = np.arange(offset, 14, 2)
    values = np.asarray([
        row["delta_sigma_random_subtracted_msun_mpc2"] if offset == 0
        else row["cross_delta_sigma_random_subtracted_msun_mpc2"]
        for row in profile
    ], dtype=float)
    sub_cov = covariance[np.ix_(indices, indices)]
    eigenvalues = np.linalg.eigvalsh(sub_cov)
    scale = max(float(np.max(np.abs(eigenvalues))), 1.0)
    rank = int(np.count_nonzero(eigenvalues > scale * 1e-10))
    if rank == 0 or not np.all(np.isfinite(values)):
        raise ValueError("invalid profile/covariance for null test")
    statistic = float(values @ np.linalg.pinv(sub_cov, rcond=1e-10) @ values)
    p_value = float(chi2.sf(statistic, rank))
    errors = np.sqrt(np.maximum(np.diag(sub_cov), 0.0))
    return {
        "null_hypothesis": "all seven radial-bin values equal zero",
        "chi_square": statistic,
        "covariance_rank": rank,
        "degrees_of_freedom": rank,
        "p_value_conditional_on_saved_covariance": p_value,
        "profile_to_marginal_error": (values / errors).tolist(),
        "covariance_min_eigenvalue": float(np.min(eigenvalues)),
        "covariance_max_eigenvalue": float(np.max(eigenvalues)),
    }


def run(stack_path: Path, output: Path) -> dict:
    if output.exists():
        raise FileExistsError(f"refusing to overwrite subset assessment: {output}")
    stack_path = stack_path.resolve()
    summary = json.loads(stack_path.read_text(encoding="utf-8"))
    if (summary.get("status") != "complete_exploratory"
            or summary.get("science_measurement") is not False
            or summary.get("photoz_uncertainty_propagated", False)):
        raise ValueError("unexpected stack status or unsupported science gate")
    by_filter = {}
    for filt, item in summary["by_filter"].items():
        jackknife = item["paired_endpoint_jackknife"]
        profile = jackknife["profile"]
        covariance = np.asarray(jackknife["covariance"], dtype=float)
        if covariance.shape != (14, 14) or len(profile) != 7:
            raise ValueError(f"unexpected seven-bin covariance/profile for {filt}")
        by_filter[filt] = {
            "n_lenses": item["n_icl_profile_supported_lenses"],
            "n_occupied_blocks": jackknife["n_occupied_blocks"],
            "tangential_null": null_test(profile, covariance, 0),
            "cross_null": null_test(profile, covariance, 1),
            "profile": profile,
            "profile_covariance_sha256": hashlib.sha256(covariance.tobytes()).hexdigest(),
        }

    output.mkdir(parents=True)
    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.8), sharex=True, sharey=True,
                             constrained_layout=True)
    colors = {"tangential": "#1769aa", "cross": "#c65d17"}
    for ax, filt in zip(axes, ("i", "i2")):
        result = by_filter[filt]
        profile = result["profile"]
        radius = np.asarray([np.sqrt(r["r_min_mpc"] * r["r_max_mpc"]) for r in profile])
        tangential = np.asarray([r["delta_sigma_random_subtracted_msun_mpc2"] for r in profile]) / 1e12
        tangential_err = np.asarray([r["delta_sigma_jackknife_error_msun_mpc2"] for r in profile]) / 1e12
        cross = np.asarray([r["cross_delta_sigma_random_subtracted_msun_mpc2"] for r in profile]) / 1e12
        cross_err = np.asarray([r["cross_jackknife_error_msun_mpc2"] for r in profile]) / 1e12
        ax.errorbar(radius, tangential, yerr=tangential_err, color=colors["tangential"],
                    marker="o", linestyle="-", capsize=2, label=r"Tangential $\Delta\Sigma$")
        ax.errorbar(radius, cross, yerr=cross_err, color=colors["cross"],
                    marker="s", linestyle="--", capsize=2, label=r"Cross $\Delta\Sigma_\times$")
        ax.axhline(0.0, color="0.35", linewidth=0.8)
        ax.set_xscale("log")
        ax.set_xlim(0.09, 3.2)
        ax.set_title(f"{filt}: {result['n_lenses']} lenses", loc="left", fontsize=10)
        ax.legend(loc="upper right", fontsize=8, frameon=False)
        ax.set_xlabel(r"$R$ [proper Mpc]")
        ax.grid(axis="y", alpha=0.18)
    axes[0].set_ylabel(r"$\Delta\Sigma$ [$10^{12}\,M_\odot\,\mathrm{Mpc}^{-2}$]")
    fig.suptitle("ICL-supported cluster subset: random-subtracted weak lensing", fontsize=11)
    figure = output / "icl_supported_wl_profiles.png"
    fig.savefig(figure, dpi=180)
    plt.close(fig)

    report = {
        "schema": "astrafek-icl-supported-wl-subset-assessment-v1",
        "status": "complete_exploratory",
        "science_measurement": False,
        "network_used": False,
        "stack_summary": str(stack_path),
        "stack_summary_sha256": sha256(stack_path),
        "by_filter": by_filter,
        "figure": str(figure),
        "figure_sha256": sha256(figure),
        "interpretation": {
            "i": "The cross component rejects zero under the saved 5-degree jackknife covariance; this is a serious systematic or covariance warning, so the tangential amplitude is not interpretable as a lensing detection.",
            "i2": "The cross component is consistent with zero under this covariance and the tangential component is nonzero, but this remains exploratory pending block-width, footprint, random-aperture, photo-z, reduced-shear, and boost validation.",
        },
        "limitations": summary["limitations"],
    }
    report_path = output / "assessment.json"
    report_path.write_text(json.dumps(report, ensure_ascii=False, indent=2, allow_nan=False) + "\n",
                           encoding="utf-8")
    print(json.dumps({
        "status": report["status"],
        "science_measurement": False,
        "by_filter": {
            filt: {
                "n_lenses": item["n_lenses"],
                "tangential_chi2": item["tangential_null"]["chi_square"],
                "tangential_p": item["tangential_null"]["p_value_conditional_on_saved_covariance"],
                "cross_chi2": item["cross_null"]["chi_square"],
                "cross_p": item["cross_null"]["p_value_conditional_on_saved_covariance"],
            } for filt, item in by_filter.items()
        },
        "figure": str(figure),
        "assessment": str(report_path),
    }, ensure_ascii=False))
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stack", type=Path, default=DEFAULT_STACK)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    run(args.stack, args.output)
