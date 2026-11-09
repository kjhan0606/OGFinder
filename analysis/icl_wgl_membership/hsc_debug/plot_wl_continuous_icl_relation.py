#!/usr/bin/env python3
"""Plot the saved continuous-ICL regression for the i2 weak-lensing sample."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

ROOT = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
RESULT = ROOT / "wl_continuous_icl_relation_i2_2026-11-06_02/summary.json"
EXPECTED_RESULT_SHA256 = "12c9456799e2671e3a0452c09a52465bdc6434c5b7546659c4bf21ac340f5fea"
DEFAULT_PNG = RESULT.parent / "wl_continuous_icl_relation_i2_v2.png"
DEFAULT_PDF = RESULT.parent / "wl_continuous_icl_relation_i2_v2.pdf"


def sha256(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def plot(png_path: Path, pdf_path: Path) -> None:
    png_path, pdf_path = Path(png_path).resolve(), Path(pdf_path).resolve()
    if png_path.exists() or pdf_path.exists():
        raise FileExistsError("refusing to overwrite an existing continuous ICL/WL figure")
    if sha256(RESULT) != EXPECTED_RESULT_SHA256:
        raise ValueError("the pinned continuous ICL/WL result changed")
    summary = json.loads(RESULT.read_text(encoding="utf-8"))
    if (summary.get("status") != "complete_exploratory"
            or summary.get("science_measurement") is not False
            or summary.get("filter") != "i2"
            or summary.get("n_common_support_lenses") != 70):
        raise ValueError("unexpected result scope or science gate")

    rows = summary["regression_rows"]
    if len(rows) != 14:
        raise ValueError("expected seven radial bins for each shear component")
    fig, axes = plt.subplots(2, 1, figsize=(7.2, 6.9), sharex=True, sharey=True)
    fig.subplots_adjust(left=0.15, right=0.98, top=0.97, bottom=0.15, hspace=0.035)
    colors = {"HC3": "#2369a1", "Block deletion": "#d56b27"}
    offsets = {"HC3": 0.975, "Block deletion": 1.025}
    components = (("tangential", "Tangential"), ("cross", "Cross"))
    all_centers = []
    for ax, (component, title) in zip(axes, components):
        selected = [row for row in rows if row["component"] == component]
        selected.sort(key=lambda row: row["bin"])
        centers = np.asarray([
            np.sqrt(float(row["r_min_mpc"]) * float(row["r_max_mpc"]))
            for row in selected
        ])
        all_centers.extend(centers.tolist())
        estimate = np.asarray([row["icl_slope_per_1sd_msun_mpc2"] / 1e12 for row in selected])
        err_hc3 = np.asarray([row["hc3_standard_error_msun_mpc2"] / 1e12 for row in selected])
        err_block = np.asarray([row["spatial_jackknife_error_msun_mpc2"] / 1e12 for row in selected])
        ax.axhline(0.0, color="#555555", linewidth=0.9, linestyle="--", zorder=0)
        ax.errorbar(centers * offsets["HC3"], estimate, yerr=err_hc3,
                    fmt="o", markersize=5, capsize=2.3, linewidth=1.1,
                    color=colors["HC3"], label="Lens-level HC3 SE", zorder=3)
        ax.errorbar(centers * offsets["Block deletion"], estimate, yerr=err_block,
                    fmt="D", markersize=4.2, capsize=2.3, linewidth=1.1,
                    color=colors["Block deletion"], label="27-cell delete-one SE", zorder=4)
        ax.set_xscale("log")
        ax.set_xlim(0.09, 3.25)
        ax.set_ylim(-220, 160)
        ax.set_ylabel("$d\\Delta\\Sigma / d(\\mathrm{ICL}/1\\,\\sigma)$\n($10^{12}\\,M_\\odot\\,\\mathrm{Mpc}^{-2}$)")
        ax.grid(axis="y", color="#dddddd", linewidth=0.6)
        ax.text(0.035, 0.94, title, transform=ax.transAxes, ha="left", va="top",
                fontsize=11, fontweight="semibold")
    axes[0].text(0.98, 0.94,
                 "2–3 Mpc: HC3 p=0.021; 7-bin BH q=0.149",
                 transform=axes[0].transAxes, ha="right", va="top", fontsize=8.5,
                 bbox={"boxstyle": "round,pad=0.25", "facecolor": "white", "edgecolor": "#cccccc", "alpha": 0.92})
    axes[1].text(0.98, 0.94, "Cross slopes: nominal p ≥ 0.20",
                 transform=axes[1].transAxes, ha="right", va="top", fontsize=8.5,
                 bbox={"boxstyle": "round,pad=0.25", "facecolor": "white", "edgecolor": "#cccccc", "alpha": 0.92})
    axes[0].legend(loc="lower right", fontsize=8, frameon=True, framealpha=0.94)
    axes[1].set_xlabel("Projected radius (proper Mpc)")
    axes[1].set_xticks([0.1, 0.2, 0.4, 0.7, 1.0, 1.5, 2.0, 3.0])
    axes[1].set_xticklabels(["0.1", "0.2", "0.4", "0.7", "1", "1.5", "2", "3"])
    fig.text(0.5, 0.045,
             "N=70; radial N=58–70; five random centers/lens; z and ln richness controlled. Approximate blocks; exploratory.",
             ha="center", va="bottom", fontsize=7.7, color="#444444")
    png_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(png_path, dpi=300, facecolor="white")
    fig.savefig(pdf_path, facecolor="white")
    plt.close(fig)
    print(json.dumps({
        "result_sha256": sha256(RESULT),
        "png": str(png_path), "png_sha256": sha256(png_path),
        "pdf": str(pdf_path), "pdf_sha256": sha256(pdf_path),
    }, indent=2))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--png", type=Path, default=DEFAULT_PNG)
    parser.add_argument("--pdf", type=Path, default=DEFAULT_PDF)
    args = parser.parse_args()
    plot(args.png, args.pdf)


if __name__ == "__main__":
    main()
