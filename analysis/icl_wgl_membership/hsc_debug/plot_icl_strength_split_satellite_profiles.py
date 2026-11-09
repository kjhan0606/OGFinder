#!/usr/bin/env python3
"""Plot the saved ICL-strength split satellite profile contrasts."""
from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


ROOT = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
DEFAULT_INPUT = ROOT / "icl_satellite_icl_strength_stacked_profiles_2026-11-06_01"
RANGES = ((5, 10), (10, 15), (15, 20), (20, 30), (30, 40), (40, 50), (50, 75))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args()
    source = args.input.resolve()
    output = (args.output or source / "icl_strength_split_satellite_profiles.png").resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite plot: {output}")
    with (source / "icl_strength_split_satellite_profile.tsv").open(
            encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))

    plt.rcParams.update({"font.size": 10, "axes.titlesize": 10,
                         "axes.labelsize": 10, "legend.fontsize": 9})
    fig, axes = plt.subplots(2, 2, figsize=(10.8, 7.0), sharex="col",
                             gridspec_kw={"height_ratios": [1.0, 0.92]})
    colors = {"low": "#3569a8", "high": "#d05a36"}
    for col, band in enumerate(("i", "i2")):
        subset = [row for row in rows if row["filter"] == band]
        centers = np.asarray([np.sqrt(lo * hi) for lo, hi in RANGES])
        low = np.asarray([float(next(r for r in subset if r["radial_bin_kpc"] == f"{lo}-{hi}")[
            "low_icl_mean_satellite_excess_njy_arcsec2"]) for lo, hi in RANGES])
        high = np.asarray([float(next(r for r in subset if r["radial_bin_kpc"] == f"{lo}-{hi}")[
            "high_icl_mean_satellite_excess_njy_arcsec2"]) for lo, hi in RANGES])
        contrast = np.asarray([float(next(r for r in subset if r["radial_bin_kpc"] == f"{lo}-{hi}")[
            "adjusted_high_minus_low_njy_arcsec2"]) for lo, hi in RANGES])
        lower = np.asarray([float(next(r for r in subset if r["radial_bin_kpc"] == f"{lo}-{hi}")[
            "cr1_t95_low"]) for lo, hi in RANGES])
        upper = np.asarray([float(next(r for r in subset if r["radial_bin_kpc"] == f"{lo}-{hi}")[
            "cr1_t95_high"]) for lo, hi in RANGES])

        top, bottom = axes[0, col], axes[1, col]
        top.plot(centers, low, "o-", color=colors["low"], label="Low ICL")
        top.plot(centers, high, "s-", color=colors["high"], label="High ICL")
        top.set_xscale("log")
        top.set_yscale("log")
        top.set_title(f"{band} band: cluster-mean satellite excess")
        top.set_ylabel(r"Excess surface brightness [nJy arcsec$^{-2}$]")
        top.legend(loc="upper right", frameon=False)
        top.grid(alpha=0.2, which="both")

        bottom.errorbar(centers, contrast,
                        yerr=np.vstack((contrast - lower, upper - contrast)),
                        fmt="o", color="#3d3d3d", capsize=3, markersize=4)
        bottom.axhline(0.0, color="#777777", linewidth=1, linestyle="--")
        bottom.set_xscale("log")
        bottom.set_xlabel("Projected radius [kpc]")
        bottom.set_ylabel("Adjusted high − low\n[nJy arcsec⁻²]")
        bottom.set_title("Redshift- and richness-adjusted contrast")
        bottom.grid(alpha=0.2, which="both")
        bottom.set_xticks(centers, [f"{lo}–{hi}" for lo, hi in RANGES])
        bottom.tick_params(axis="x", labelrotation=35)

    fig.tight_layout()
    fig.savefig(output, dpi=180, bbox_inches="tight")
    print(output)


if __name__ == "__main__":
    main()
