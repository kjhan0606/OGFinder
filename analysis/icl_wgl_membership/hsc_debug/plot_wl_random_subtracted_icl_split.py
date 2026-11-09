#!/usr/bin/env python3
"""Plot an existing exploratory random-subtracted ICL/WL split summary."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
DEFAULT_INPUT = DATA / "wl_random_subtracted_icl_split_2026-11-05_02/summary.json"


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
    if (summary.get("schema") != "astrafek-hsc-random-subtracted-icl-split-v1"
            or summary.get("status") != "complete_exploratory"
            or summary.get("science_measurement") is not False):
        raise ValueError("input is not the expected exploratory ICL/WL split")
    rows = summary["profile"]
    rmin = np.asarray([row["r_min_mpc"] for row in rows], dtype=float)
    rmax = np.asarray([row["r_max_mpc"] for row in rows], dtype=float)
    radius = np.sqrt(rmin * rmax)
    low = np.asarray([row["low_icl_delta_sigma_msun_mpc2"] for row in rows]) / 1e12
    high = np.asarray([row["high_icl_delta_sigma_msun_mpc2"] for row in rows]) / 1e12
    contrast = np.asarray([row["high_minus_low_delta_sigma_msun_mpc2"] for row in rows]) / 1e12
    error = np.asarray([row["high_minus_low_delta_sigma_error_msun_mpc2"] for row in rows]) / 1e12

    fig, (ax_profile, ax_contrast) = plt.subplots(
        1, 2, figsize=(10.0, 4.4), sharex=True,
        gridspec_kw={"width_ratios": [1.05, 1.0], "wspace": 0.0},
    )
    color_low, color_high, color_delta = "#3572A5", "#C44E52", "#333333"
    ax_profile.plot(radius, low, "o-", color=color_low, label="Low ICL", lw=1.5, ms=4)
    ax_profile.plot(radius, high, "s-", color=color_high, label="High ICL", lw=1.5, ms=4)
    ax_profile.set_xscale("log")
    ax_profile.set_xlim(float(rmin.min()), float(rmax.max()))
    ax_profile.set_ylabel(r"$\Delta\Sigma\ [10^{12}\ M_\odot\,\mathrm{Mpc}^{-2}]$")
    ax_profile.set_title("", pad=4)
    ax_profile.text(0.04, 0.06, "(a) Random-subtracted profiles",
                    transform=ax_profile.transAxes, va="bottom", ha="left", fontsize=10)
    ax_profile.legend(loc="upper right", frameon=False, fontsize=9)
    ax_profile.grid(alpha=0.18, which="both")

    ax_contrast.axhline(0.0, color="0.35", lw=1.0, ls="--")
    ax_contrast.errorbar(radius, contrast, yerr=error, fmt="o", color=color_delta,
                         ecolor="#777777", elinewidth=1.2, capsize=2.5, ms=4)
    ax_contrast.set_xscale("log")
    ax_contrast.set_xlim(float(rmin.min()), float(rmax.max()))
    ax_contrast.set_ylabel(r"High $-$ low $\Delta\Sigma\ [10^{12}\ M_\odot\,\mathrm{Mpc}^{-2}]$", labelpad=7)
    ax_contrast.yaxis.set_label_position("right")
    ax_contrast.yaxis.tick_right()
    ax_contrast.text(0.04, 0.96, "(b) Contrast with paired block error",
                     transform=ax_contrast.transAxes, va="top", ha="left", fontsize=10)
    ax_contrast.grid(alpha=0.18, which="both")
    ax_profile.set_xlabel("Proper projected radius [Mpc]")
    ax_contrast.set_xlabel("Proper projected radius [Mpc]")
    fig.suptitle(
        f"HSC {summary['filter']} ICL split; source bins 2–4; "
        f"{summary['n_matched_lens_pairs']} matched pairs; "
        f"{summary['block_width_deg']:g}° blocks",
        y=1.01,
    )
    fig.subplots_adjust(left=0.10, right=0.91, bottom=0.16, top=0.88)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    print(json.dumps({
        "status": "complete_exploratory",
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
    output = args.output or args.input.with_name("wl_random_subtracted_icl_split.png")
    plot(args.input, output)


if __name__ == "__main__":
    main()
