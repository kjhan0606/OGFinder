#!/usr/bin/env python3
"""Visualize the saved joint source-density-boost and DeltaSigma covariance.

The input is the corrected v2 diagnostic summary. This figure is descriptive:
it does not apply a boost correction or change the science-measurement gate.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


SCHEMA = "astrafek-wl-source-density-boost-joint-covariance-v2"
METRICS = {
    "count_boost": "Count",
    "shape_weight_boost": "Shape weight",
    "sigma_crit_weight_boost": r"$\Sigma_{\rm crit}^{-2}$ weight",
    "estimator_weight_boost": "S19A estimator weight",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def plot(summary_path: Path, output_path: Path, block_width: str = "2.5") -> dict:
    summary_path = summary_path.resolve()
    output_path = output_path.resolve()
    data = json.loads(summary_path.read_text(encoding="utf-8"))
    if (data.get("schema") != SCHEMA or data.get("status") != "complete_diagnostic"
            or data.get("science_measurement") is not False
            or data.get("boost_correction_applied") is not False
            or data.get("hsc_network_used") is not False):
        raise ValueError("input is not the saved offline, uncorrected v2 diagnostic")
    result = data["results_by_block_width_deg"][block_width]
    joint = result["joint_boost_delta_sigma_jackknife_covariance"]
    if (joint.get("saved_stack_point_reproduced") is not True
            or joint.get("saved_stack_covariance_reproduced_on_pair_endpoint_blocks") is not True
            or result.get("n_lenses") != 435 or result.get("n_random_centers") != 2175):
        raise ValueError("joint covariance does not pass the saved-stack reproduction gates")

    order = joint["parameter_order"]
    corr = np.asarray(joint["correlation"], dtype=float)
    cov = np.asarray(joint["covariance"], dtype=float)
    if corr.shape != (42, 42) or cov.shape != (42, 42) or len(order) != 42:
        raise ValueError("expected the 42-parameter joint covariance")
    if (not np.isfinite(corr).all() or not np.isfinite(cov).all()
            or not np.allclose(corr, corr.T, atol=1e-10)
            or not np.allclose(cov, cov.T, atol=1e-10)
            or not np.allclose(np.diag(corr), 1.0, atol=1e-8)):
        raise ValueError("joint covariance is non-finite, asymmetric, or has invalid correlation diagonal")

    rows = data["profile_rows"]
    rmin = np.asarray([row["r_min_mpc"] for row in rows], dtype=float)
    rmax = np.asarray([row["r_max_mpc"] for row in rows], dtype=float)
    radius = np.sqrt(rmin * rmax)
    xerr = np.vstack((radius - rmin, rmax - radius))
    metric_cols = {
        name: np.asarray([row[f"{name}_{block_width.replace('.', 'p')}deg"] for row in rows], dtype=float)
        for name in METRICS
    }
    metric_err_cols = {
        name: np.asarray([row[f"{name}_{block_width.replace('.', 'p')}deg_jackknife_sigma"] for row in rows], dtype=float)
        for name in METRICS
    }
    delta_indices = [i for i, item in enumerate(order)
                     if item["metric"] in ("delta_sigma_random_subtracted_msun_mpc2",
                                           "cross_delta_sigma_random_subtracted_msun_mpc2")]
    delta_cov = cov[np.ix_(delta_indices, delta_indices)]
    delta_sigma = np.asarray([pair[0] for pair in result["delta_sigma_random_subtracted_full_sample"]], dtype=float) / 1e13
    cross_sigma = np.asarray([pair[1] for pair in result["delta_sigma_random_subtracted_full_sample"]], dtype=float) / 1e13
    delta_err = np.sqrt(np.diag(delta_cov)[0::2]) / 1e13
    cross_err = np.sqrt(np.diag(delta_cov)[1::2]) / 1e13
    if not all(np.isfinite(x).all() for x in (radius, xerr, delta_sigma, cross_sigma, delta_err, cross_err)):
        raise ValueError("profile contains a non-finite value")

    fig = plt.figure(figsize=(11.2, 9.4), constrained_layout=True)
    grid = fig.add_gridspec(3, 1, height_ratios=(1.0, 1.05, 2.15))
    ax_boost = fig.add_subplot(grid[0])
    ax_ds = fig.add_subplot(grid[1], sharex=ax_boost)
    ax_corr = fig.add_subplot(grid[2])

    colors = {"count_boost": "#777777", "shape_weight_boost": "#d08b00",
              "sigma_crit_weight_boost": "#3a923a", "estimator_weight_boost": "#1768ac"}
    for name, label in METRICS.items():
        ax_boost.errorbar(radius, metric_cols[name], xerr=xerr,
                          yerr=metric_err_cols[name], marker="o", markersize=3.5,
                          linewidth=1.1, capsize=2, color=colors[name], label=label)
    ax_boost.axhline(1.0, color="0.25", linestyle="--", linewidth=0.8)
    ax_boost.set_ylabel("Density ratio")
    ax_boost.set_ylim(0.8, max(2.2, float(np.max(metric_cols["shape_weight_boost"] + metric_err_cols["shape_weight_boost"]) * 1.08)))
    ax_boost.legend(loc="upper right", ncol=2, frameon=True, fontsize=8)
    ax_boost.text(0.015, 0.95, f"{block_width}° delete-one blocks; 435 lenses / 2,175 random centers",
                  transform=ax_boost.transAxes, va="top", fontsize=8,
                  bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.82, "pad": 2})

    ax_ds.errorbar(radius, delta_sigma, xerr=xerr, yerr=delta_err, color="#1768ac",
                   fmt="o-", markersize=3.7, capsize=2, linewidth=1.1, label=r"Tangential $\Delta\Sigma$")
    ax_ds.errorbar(radius, cross_sigma, xerr=xerr, yerr=cross_err, color="#b54a24",
                   fmt="s--", markersize=3.4, capsize=2, linewidth=1.0, label="Cross component")
    ax_ds.axhline(0, color="0.3", linestyle=":", linewidth=0.8)
    ax_ds.set_ylabel(r"$\Delta\Sigma$ [$10^{13}\,M_\odot\,\mathrm{Mpc}^{-2}$]")
    ax_ds.legend(loc="upper right", ncol=2, frameon=True, fontsize=8)
    ax_ds.text(0.015, 0.95, "Random-subtracted profile; no boost correction applied",
               transform=ax_ds.transAxes, va="top", fontsize=8,
               bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.82, "pad": 2})

    image = ax_corr.imshow(corr, origin="upper", vmin=-1, vmax=1, cmap="coolwarm", interpolation="nearest")
    labels = []
    for item in order:
        metric = item["metric"]
        radial = int(item["radial_bin"]) + 1
        short = {"count_boost": "N", "shape_weight_boost": "W", "sigma_crit_weight_boost": "S",
                 "estimator_weight_boost": "E",
                 "delta_sigma_random_subtracted_msun_mpc2": "T",
                 "cross_delta_sigma_random_subtracted_msun_mpc2": "X"}[metric]
        labels.append(f"{short}{radial}")
    ticks = [0, 6, 13, 20, 27, 34, 41]
    ax_corr.set_xticks(ticks, [labels[i] for i in ticks], fontsize=8)
    ax_corr.set_yticks(ticks, [labels[i] for i in ticks], fontsize=8)
    ax_corr.set_xlabel("Parameters: N=count, W=shape weight, S=$\\Sigma_{\\rm crit}^{-2}$, E=S19A estimator, T=tangential, X=cross")
    ax_corr.set_ylabel("Joint covariance parameter")
    for boundary in (6.5, 13.5, 20.5, 27.5):
        ax_corr.axhline(boundary, color="black", linewidth=0.6, alpha=0.65)
        ax_corr.axvline(boundary, color="black", linewidth=0.6, alpha=0.65)
    colorbar = fig.colorbar(image, ax=ax_corr, fraction=0.025, pad=0.012)
    colorbar.set_label("Jackknife correlation")
    ax_corr.text(0.99, 0.01, "Exploratory covariance; science gate closed", transform=ax_corr.transAxes,
                 ha="right", va="bottom", color="black", fontsize=8,
                 bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.8, "pad": 2})

    for ax in (ax_boost, ax_ds):
        ax.set_xscale("log")
        ax.set_xlim(0.09, 3.3)
        ax.set_xticks([0.1, 0.3, 1.0, 3.0], labels=["0.1", "0.3", "1", "3"])
        ax.grid(color="0.88", linewidth=0.65)
        ax.set_axisbelow(True)
        ax.tick_params(direction="in", top=True, right=True)
    ax_boost.tick_params(labelbottom=False)
    ax_ds.set_xlabel("Proper projected radius [Mpc]")
    ax_corr.set_title("Full 42-parameter joint delete-one correlation", fontsize=10, loc="left")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    if output_path.exists():
        raise FileExistsError(f"refusing to overwrite existing figure: {output_path}")
    fig.savefig(output_path, dpi=180, facecolor="white")
    plt.close(fig)
    return {"summary": str(summary_path), "summary_sha256": sha256(summary_path),
            "figure": str(output_path), "figure_sha256": sha256(output_path),
            "block_width_deg": float(block_width), "n_lenses": result["n_lenses"],
            "n_random_centers": result["n_random_centers"],
            "n_joint_parameters": len(order), "science_measurement": False}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--block-width-deg", choices=("2.5", "5.0"), default="2.5")
    args = parser.parse_args()
    print(json.dumps(plot(args.summary, args.output, args.block_width_deg), indent=2))


if __name__ == "__main__":
    main()
