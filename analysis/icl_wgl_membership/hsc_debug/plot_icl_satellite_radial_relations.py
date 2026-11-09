#!/usr/bin/env python3
"""Plot continuous and median-split ICL associations with satellite profiles."""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
DEFAULT_INPUT = ROOT / "icl_satellite_radial_relations_joint_2026-11-06_01"
RANGES = ((5, 10), (10, 15), (15, 20), (20, 30), (30, 40), (40, 50), (50, 75))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args()
    source = args.input.resolve()
    output = (args.output or source / "icl_satellite_radial_relations.png").resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite plot: {output}")
    with (source / "icl_satellite_radial_relation_profiles.tsv").open(
            encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))

    plt.rcParams.update({"font.size": 10, "axes.titlesize": 10,
                         "axes.labelsize": 10, "xtick.labelsize": 9})
    fig, axes = plt.subplots(2, 2, figsize=(10.8, 7.2), sharex="col")
    centers = np.asarray([np.sqrt(lo * hi) for lo, hi in RANGES])
    methods = (("continuous_icl", "Continuous ICL slope per 1σ"),
               ("median_split", "Median-split high − low"))
    for col, band in enumerate(("i", "i2")):
        for row_idx, (relation, title) in enumerate(methods):
            ax = axes[row_idx, col]
            subset = [r for r in rows if r["filter"] == band and r["relation"] == relation]
            subset_by_bin = {r["radial_bin_kpc"]: r for r in subset}
            vals, lows, highs = [], [], []
            for lo, hi in RANGES:
                item = subset_by_bin[f"{lo}-{hi}"]
                vals.append(float(item["adjusted_effect_njy_arcsec2"]))
                lows.append(float(item["cr1_t95_low"]))
                highs.append(float(item["cr1_t95_high"]))
            vals, lows, highs = map(np.asarray, (vals, lows, highs))
            ax.errorbar(centers, vals, yerr=np.vstack((vals - lows, highs - vals)),
                        fmt="o-", color="#426c9e", capsize=3, markersize=4,
                        linewidth=1.2)
            ax.axhline(0.0, color="#777777", linewidth=1, linestyle="--")
            ax.set_xscale("log")
            ax.set_title(f"{band} band: {title}")
            ax.set_ylabel("Adjusted satellite excess\n[nJy arcsec⁻²]")
            ax.grid(alpha=0.2, which="both")
            if row_idx == 1:
                ax.set_xlabel("Projected radius [kpc]")
                ax.set_xticks(centers, [f"{lo}–{hi}" for lo, hi in RANGES])
                ax.tick_params(axis="x", labelrotation=35)
    fig.text(0.5, 0.012,
             "Whiskers are unadjusted CR1 95% intervals; p-values correct jointly across 28 tests.",
             ha="center", fontsize=9)
    fig.tight_layout(rect=(0, 0.035, 1, 1))
    fig.savefig(output, dpi=180, bbox_inches="tight")
    print(output)


if __name__ == "__main__":
    main()
