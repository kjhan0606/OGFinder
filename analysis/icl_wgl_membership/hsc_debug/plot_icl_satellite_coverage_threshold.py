#!/usr/bin/env python3
"""Plot radial support and 50--75 kpc association sensitivity to coverage floor."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def load_radial(path: Path) -> dict[tuple[str, str], dict]:
    import csv
    with path.open(encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    return {(r["satellite_radial_bin_kpc"], r["filter"]): r for r in rows}


def plot(summary_path: Path, support60_path: Path, support40_path: Path,
         output_path: Path) -> dict:
    summary_path, support60_path = summary_path.resolve(), support60_path.resolve()
    support40_path, output_path = support40_path.resolve(), output_path.resolve()
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    if (summary.get("schema") != "astrafek-icl-satellite-coverage-threshold-sensitivity-v1"
            or summary.get("status") != "complete_exploratory_sensitivity"
            or summary.get("science_measurement") is not False):
        raise ValueError("unexpected coverage-threshold comparison summary")
    s60, s40 = load_radial(support60_path), load_radial(support40_path)
    bins = ["5-10", "10-15", "15-20", "20-30", "30-40", "40-50", "50-75"]
    if len(s60) != 14 or len(s40) != 14:
        raise ValueError("expected seven annuli by two filters for both floors")

    fig, (ax_support, ax_assoc) = plt.subplots(
        1, 2, figsize=(10.3, 4.5), gridspec_kw={"width_ratios": (1.4, 1.0)},
        constrained_layout=True,
    )
    x = np.arange(len(bins), dtype=float)
    width = 0.36
    for offset, floor, table, color, hatch in (
            (-width/2, "60% floor", s60, "#1768ac", ""),
            (width/2, "40% floor", s40, "#d08b00", "//")):
        total = np.asarray([int(table[(rb, "i")]["common_cluster_filter_pairs"])
                            + int(table[(rb, "i2")]["common_cluster_filter_pairs"])
                            for rb in bins])
        bars = ax_support.bar(x + offset, total, width=width, label=floor,
                              color=color, edgecolor="white", linewidth=0.5,
                              hatch=hatch)
        for idx in (0, len(bins)-1):
            ax_support.text(bars[idx].get_x()+bars[idx].get_width()/2,
                            bars[idx].get_height()+3, str(int(total[idx])),
                            ha="center", va="bottom", fontsize=8)
    ax_support.set_xticks(x, [b.replace("-", "–") for b in bins], rotation=35, ha="right")
    ax_support.set_ylabel("ICL/satellite common pairs")
    ax_support.set_xlabel("Satellite projected radius [kpc]")
    ax_support.set_ylim(0, 250)
    ax_support.grid(axis="y", color="0.88", linewidth=0.65)
    ax_support.set_axisbelow(True)
    ax_support.legend(loc="upper right", frameon=True, fontsize=8)
    ax_support.text(0.03, 0.96, "800 parent images; ICL annulus 200–300 kpc",
                    transform=ax_support.transAxes, va="top", fontsize=8,
                    bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.85, "pad": 2})

    colors = {"i": "#1768ac", "i2": "#b54a24"}
    markers = {"60%": "o", "40%": "s"}
    offsets = {"i": -0.10, "i2": 0.10}
    for band, xpos in (("i", 0.), ("i2", 1.)):
        result = summary["results_by_filter"][band]
        for floor, key, ci_key in (
                ("60%", "partial_r_60pct", "partial_r_60pct_block_bootstrap_ci95"),
                ("40%", "partial_r_40pct", "partial_r_40pct_block_bootstrap_ci95")):
            estimate = float(result[key])
            low, high = map(float, result[ci_key])
            yerr = np.asarray([[estimate-low], [high-estimate]])
            ax_assoc.errorbar(xpos+offsets[band], estimate, yerr=yerr,
                              fmt=markers[floor], color=colors[band], markersize=6,
                              capsize=3, linewidth=1.2,
                              label=f"{band}, {floor} floor")
            ax_assoc.annotate(f"n={result['n_'+floor.replace('%','pct')+'_support']}",
                              (xpos+offsets[band], high), xytext=(0, 5),
                              textcoords="offset points", ha="center", fontsize=7,
                              color=colors[band])
    ax_assoc.axhline(0, color="0.3", linestyle="--", linewidth=0.9)
    ax_assoc.set_xticks([0, 1], ["i", "i2"])
    ax_assoc.set_xlim(-0.35, 1.35)
    ax_assoc.set_ylim(-0.85, 0.8)
    ax_assoc.set_ylabel("Partial correlation r")
    ax_assoc.set_xlabel("Filter; ICL 200–300 / satellite 50–75 kpc")
    ax_assoc.grid(axis="y", color="0.88", linewidth=0.65)
    ax_assoc.set_axisbelow(True)
    ax_assoc.legend(loc="lower left", frameon=True, fontsize=7, ncol=2)
    ax_assoc.text(0.97, 0.97, "5° block bootstrap; 2,000 draws",
                  transform=ax_assoc.transAxes, ha="right", va="top", fontsize=8,
                  bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.85, "pad": 2})
    ax_assoc.text(0.03, 0.84, "Diagnostic; neither floor adopted",
                  transform=ax_assoc.transAxes, ha="left", va="top", fontsize=7,
                  bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.85, "pad": 2})
    fig.suptitle("Support recovery trades against estimator sensitivity", fontsize=11)
    if output_path.exists():
        raise FileExistsError(f"refusing to overwrite {output_path}")
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=180, facecolor="white")
    plt.close(fig)
    return {"figure": str(output_path), "figure_sha256": sha256(output_path),
            "summary_sha256": sha256(summary_path),
            "support60_sha256": sha256(support60_path),
            "support40_sha256": sha256(support40_path), "science_measurement": False}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--support60", type=Path, required=True)
    parser.add_argument("--support40", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(plot(args.summary, args.support60, args.support40, args.output), indent=2))


if __name__ == "__main__":
    main()
