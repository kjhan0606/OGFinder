#!/usr/bin/env python3
"""Plot S23B target unmasked-area support by filter, radius, and mask size."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


DEFAULT_SUMMARY = Path(
    "/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/"
    "s23b_icl_mask_support_assessment_2026-11-07_01/summary.json"
)
DEFAULT_OUTPUT = DEFAULT_SUMMARY.parent / "radial_unmasked_area_by_filter_v4.png"
FACTORS = ("2.5", "3.5", "5.0")
COLORS = {"2.5": "#0072B2", "3.5": "#D55E00", "5.0": "#009E73"}


def plot(summary_path: Path, output: Path) -> None:
    summary_path, output = summary_path.resolve(), output.resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite {output}")
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    if (summary.get("schema") != "astrafek-s23b-icl-mask-support-assessment-v1"
            or summary.get("science_measurement") is not False
            or summary.get("images_ready_for_stacking") is not False):
        raise ValueError("input is not the pinned descriptive support assessment")
    bins = ["50-100", "100-150", "150-200", "200-300", "300-400", "400-500"]
    x = np.array([75, 125, 175, 250, 350, 450], dtype=float)
    plt.rcParams.update({"font.size": 10, "axes.labelsize": 10, "legend.fontsize": 9})
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 4.6), sharex=True, sharey=True,
                             gridspec_kw={"wspace": 0.03})
    for ax, band in zip(axes, ("i", "i2")):
        for factor in FACTORS:
            stats = summary["radial_unmasked_area_fraction"][factor][band]
            median = np.array([stats[b]["median"] for b in bins])
            low = np.array([stats[b]["p10"] for b in bins])
            high = np.array([stats[b]["p90"] for b in bins])
            color = COLORS[factor]
            ax.fill_between(x, low, high, color=color, alpha=0.13, linewidth=0)
            ax.plot(x, median, color=color, marker="o", linewidth=1.8,
                    markersize=4.5, label=f"Kron {factor}")
        x_label = 0.05
        align = "left"
        ax.text(x_label, 0.94, f"HSC filter {band}", transform=ax.transAxes,
                ha=align, va="top", fontsize=10)
        ax.set_xlim(50, 500)
        ax.set_ylim(0, 1)
        ax.set_xticks([50, 100, 200, 300, 400, 500])
        ax.set_yticks(np.linspace(0, 1, 6))
        ax.grid(True, color="#b8b8b8", alpha=0.4, linewidth=0.6)
        ax.legend(loc="upper right", frameon=True, facecolor="white", framealpha=0.92)
    # At the shared panel boundary, the left-panel 500-kpc and right-panel
    # 50-kpc tick labels occupy the same position. Keep both tick marks and
    # suppress only the duplicate boundary label on the right panel.
    axes[1].set_xticklabels(["", "100", "200", "300", "400", "500"])
    axes[0].set_ylabel("Unmasked annular area fraction")
    fig.text(0.5, 0.095, "Cluster-centric radius (proper kpc)",
             ha="center", va="center", fontsize=10)
    fig.text(0.5, 0.035,
             "Lines show medians; shaded bands show 10th–90th percentiles. "
             "Automated-mask support only; no selection cut adopted.",
             ha="center", va="bottom", fontsize=9)
    fig.subplots_adjust(left=0.09, right=0.99, top=0.98, bottom=0.22, wspace=0.03)
    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, dpi=180, facecolor="white")
    plt.close(fig)
    print(json.dumps({"output": str(output), "n_targets": summary["n_mask_support_success"],
                      "science_measurement": False}, indent=2))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", type=Path, default=DEFAULT_SUMMARY)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    plot(args.summary, args.output)


if __name__ == "__main__":
    main()
