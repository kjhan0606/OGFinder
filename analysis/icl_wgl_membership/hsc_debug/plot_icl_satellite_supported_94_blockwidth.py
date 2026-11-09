#!/usr/bin/env python3
"""Plot the saved ICL--satellite block-width bootstrap comparison."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import NullLocator


DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
DEFAULT_INPUT = DATA / "icl_satellite_supported_94_blockwidth_2026-11-05_01/summary.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def plot(input_path: Path, output_path: Path) -> Path:
    input_path, output_path = Path(input_path).resolve(), Path(output_path).resolve()
    if output_path.exists():
        raise FileExistsError(f"refusing to overwrite figure: {output_path}")
    summary = json.loads(input_path.read_text(encoding="utf-8"))
    if (summary.get("schema") != "astrafek-icl-satellite-supported-94-blockwidth-v1"
            or summary.get("status") != "complete_diagnostic"
            or summary.get("science_measurement") is not False):
        raise ValueError("unexpected input block-width diagnostic")
    widths = np.asarray(summary["block_widths_deg"], dtype=float)
    fig, axes = plt.subplots(
        1, 2, figsize=(9.2, 4.2), sharex=True, sharey=True,
        gridspec_kw={"wspace": 0.0},
    )
    colors = {"previous74": "#3572A5", "support94": "#C44E52"}
    names = {"previous74": "Prior 74-pair subset", "support94": "Full 94-pair support"}
    offsets = {"previous74": 0.97, "support94": 1.03}
    for axis, band, title in zip(axes, ("i", "i2"), ("(a) HSC i", "(b) HSC i2")):
        for sample in ("previous74", "support94"):
            estimates, lower, upper = [], [], []
            for width in widths:
                item = summary["results_by_width_and_filter"][str(float(width))][band][sample]
                estimate = float(item["partial_r_z_richness"])
                ci = item["block_bootstrap_ci95"]
                estimates.append(estimate)
                lower.append(estimate - float(ci[0]))
                upper.append(float(ci[1]) - estimate)
            axis.errorbar(
                widths * offsets[sample], estimates, yerr=[lower, upper],
                fmt="o-", color=colors[sample], lw=1.4, ms=4, capsize=2.5,
                label=names[sample],
            )
        axis.axhline(0.0, color="0.35", lw=1.0, ls="--")
        axis.set_xscale("log")
        axis.xaxis.set_minor_locator(NullLocator())
        axis.set_xticks(widths, labels=[f"{width:g}" for width in widths])
        axis.set_xlim(2.1, 12.0)
        axis.set_ylim(-0.65, 0.70)
        axis.grid(alpha=0.20, which="both")
        axis.text(0.04, 0.96, title, transform=axis.transAxes,
                  va="top", ha="left", fontsize=10)
        axis.set_xlabel("Equal-area cell width [deg]")
    axes[0].set_ylabel("Partial correlation r\n(controls: redshift, richness)")
    axes[0].legend(loc="lower left", frameon=False, fontsize=9)
    axes[1].yaxis.tick_right()
    axes[1].yaxis.set_label_position("right")
    axes[1].set_ylabel("Partial correlation r\n(controls: redshift, richness)")
    fig.suptitle("ICL–satellite association: block-width sensitivity", y=1.02)
    fig.subplots_adjust(left=0.10, right=0.90, bottom=0.16, top=0.88)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(json.dumps({
        "status": "complete_diagnostic",
        "science_measurement": False,
        "input_summary_sha256": sha256(input_path),
        "figure": str(output_path),
        "figure_sha256": sha256(output_path),
    }, indent=2))
    return output_path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args()
    output = args.output or args.input.with_name("icl_satellite_94_blockwidth.png")
    plot(args.input, output)


if __name__ == "__main__":
    main()
