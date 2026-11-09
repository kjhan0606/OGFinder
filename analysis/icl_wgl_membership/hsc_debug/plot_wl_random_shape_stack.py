#!/usr/bin/env python3
"""Plot the matched random-subtracted HSC weak-lensing profile."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def plot(summary_path: Path, output_path: Path) -> None:
    summary_path = summary_path.resolve()
    output_path = output_path.resolve()
    data = json.loads(summary_path.read_text(encoding="utf-8"))
    if (data.get("schema") != "hsc-s19a-random-center-stack-run-v1"
            or data.get("science_measurement") is not False
            or data.get("random_lens_subtraction_applied") is not True
            or data.get("spatial_block_covariance_applied") is not True):
        raise ValueError("expected the completed exploratory random-subtracted stack")
    profile = data.get("paired_spatial_jackknife", {}).get("profile")
    sample = data.get("input", {}).get("matched_sample", {})
    jackknife = data.get("paired_spatial_jackknife", {})
    n_blocks = jackknife.get("n_occupied_blocks")
    if (not profile or len(profile) != 7
            or sample.get("n_common_shape_support_lenses") != 435
            or not isinstance(n_blocks, int) or n_blocks < 3):
        raise ValueError("unexpected radial profile or common-support lens sample")

    r_min = np.asarray([row["r_min_mpc"] for row in profile], dtype=float)
    r_max = np.asarray([row["r_max_mpc"] for row in profile], dtype=float)
    radius = np.sqrt(r_min * r_max)
    xerr = np.vstack((radius - r_min, r_max - radius))
    ds = np.asarray([row["delta_sigma_random_subtracted_msun_mpc2"] for row in profile]) / 1e13
    ds_err = np.asarray([row["delta_sigma_jackknife_error_msun_mpc2"] for row in profile]) / 1e13
    cross = np.asarray([row["cross_delta_sigma_random_subtracted_msun_mpc2"] for row in profile]) / 1e13
    cross_err = np.asarray([row["cross_jackknife_error_msun_mpc2"] for row in profile]) / 1e13
    if not all(np.isfinite(values).all() for values in (radius, xerr, ds, ds_err, cross, cross_err)):
        raise ValueError("profile contains non-finite values")

    fig, (ax_ds, ax_cross) = plt.subplots(
        1, 2, figsize=(9.6, 4.4), sharex=True,
        gridspec_kw={"wspace": 0.02}, constrained_layout=True,
    )
    for ax in (ax_ds, ax_cross):
        ax.set_xscale("log")
        ax.set_xlim(0.09, 3.3)
        ax.set_xticks([0.1, 0.3, 1.0, 3.0], labels=["0.1", "0.3", "1", "3"])
        ax.grid(axis="both", color="0.87", linewidth=0.7)
        ax.set_axisbelow(True)
        ax.axhline(0.0, color="0.28", linewidth=0.9, linestyle="--")
        ax.tick_params(direction="in", top=True, right=True)
    ax_ds.errorbar(radius, ds, xerr=xerr, yerr=ds_err, fmt="o", color="#1768ac",
                   markersize=4.5, capsize=2.5, linewidth=1.2)
    ax_cross.errorbar(radius, cross, xerr=xerr, yerr=cross_err, fmt="s", color="#b54a24",
                      markersize=4.2, capsize=2.5, linewidth=1.2)
    ax_ds.set_ylim(0, 19)
    cross_lim = max(2.0, float(np.max(np.abs(cross) + cross_err) * 1.15))
    ax_cross.set_ylim(-cross_lim, cross_lim)
    ax_ds.set_ylabel(r"$\Delta\Sigma$ [$10^{13}\,M_\odot\,\mathrm{Mpc}^{-2}$]")
    ax_cross.set_ylabel(r"Cross component [$10^{13}\,M_\odot\,\mathrm{Mpc}^{-2}$]")
    ax_ds.set_xlabel("Proper projected radius [Mpc]")
    ax_cross.set_xlabel("Proper projected radius [Mpc]")
    ax_ds.text(0.97, 0.96, r"Random-subtracted $\Delta\Sigma$", transform=ax_ds.transAxes,
               ha="right", va="top", fontsize=10, fontweight="semibold",
               bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.85, "pad": 2})
    ax_cross.text(0.97, 0.96, "Cross component", transform=ax_cross.transAxes,
                  ha="right", va="top", fontsize=10, fontweight="semibold",
                  bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.85, "pad": 2})
    ax_ds.text(0.04, 0.04,
               f"435 common-support lenses\n{n_blocks} equal-area sky blocks",
               transform=ax_ds.transAxes, ha="left", va="bottom", fontsize=8,
               bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.85, "pad": 2})
    ax_cross.text(0.96, 0.04, "Exploratory; science gate closed", transform=ax_cross.transAxes,
                  ha="right", va="bottom", fontsize=8,
                  bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.85, "pad": 2})
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=180, facecolor="white")
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    plot(args.summary, args.output)
    print(json.dumps({"summary": str(args.summary.resolve()),
                      "figure": str(args.output.resolve())}, ensure_ascii=False))


if __name__ == "__main__":
    main()
