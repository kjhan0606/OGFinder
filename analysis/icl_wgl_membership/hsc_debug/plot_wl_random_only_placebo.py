#!/usr/bin/env python3
"""Plot the saved random-subtracted lensing profile against its random null."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
STACK = DATA / "wl_random_shape_stack_lens_reuse_2026-11-04_01/summary.json"
PLACEBO_DIR = DATA / "wl_random_only_placebo_2026-11-06_01"
PLACEBO = PLACEBO_DIR / "summary.json"
PNG = PLACEBO_DIR / "random_subtracted_vs_random_only.png"
PDF = PLACEBO_DIR / "random_subtracted_vs_random_only.pdf"
MANIFEST = PLACEBO_DIR / "figure_manifest.json"
SCALE = 1e13


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def load_and_validate() -> tuple[dict, dict, list[dict], list[dict]]:
    for path in (STACK, PLACEBO):
        if not path.is_file():
            raise FileNotFoundError(path)
    stack, placebo = (json.loads(path.read_text(encoding="utf-8"))
                      for path in (STACK, PLACEBO))
    if (stack.get("science_measurement") is not False
            or stack.get("random_lens_subtraction_applied") is not True
            or placebo.get("science_measurement") is not False
            or placebo.get("status") != "complete_exploratory_diagnostic"):
        raise ValueError("plot inputs have unexpected completion or science-gate status")
    subtracted = stack["paired_spatial_jackknife"]["profile"]
    random = placebo["profile"]
    if len(subtracted) != 7 or len(random) != 7:
        raise ValueError("expected the seven saved HSC lensing radial bins")
    for index, (a, b) in enumerate(zip(subtracted, random, strict=True)):
        for field_a, field_b in (
            ("random_delta_sigma_msun_mpc2", "random_only_delta_sigma_msun_mpc2"),
            ("random_cross_delta_sigma_msun_mpc2", "random_only_cross_delta_sigma_msun_mpc2"),
        ):
            if not np.isclose(a[field_a], b[field_b], rtol=1e-12, atol=1e-8):
                raise ValueError(f"random-only source does not match paired stack in bin {index}")
        if (a["r_min_mpc"] != b["r_min_mpc"] or a["r_max_mpc"] != b["r_max_mpc"]):
            raise ValueError(f"radial-bin mismatch in bin {index}")
    return stack, placebo, subtracted, random


def plot(output_png: Path = PNG, output_pdf: Path = PDF) -> dict:
    output_png, output_pdf = Path(output_png).resolve(), Path(output_pdf).resolve()
    if output_png.exists() or output_pdf.exists():
        raise FileExistsError("refusing to overwrite an existing figure")
    stack, placebo, subtracted, random = load_and_validate()
    radius = np.asarray([np.sqrt(row["r_min_mpc"] * row["r_max_mpc"])
                         for row in subtracted])
    plt.rcParams.update({"font.size": 11, "axes.titlesize": 11,
                         "axes.labelsize": 11, "legend.fontsize": 9})
    fig, axes = plt.subplots(2, 1, figsize=(7.4, 6.8), sharex=True,
                             gridspec_kw={"hspace": 0.0})
    series = (
        ("Tangential $\\Delta\\Sigma$", "delta_sigma_random_subtracted_msun_mpc2",
         "random_only_delta_sigma_msun_mpc2", "delta_sigma_jackknife_error_msun_mpc2",
         "jackknife_error_delta_sigma_msun_mpc2"),
        ("Cross component", "cross_delta_sigma_random_subtracted_msun_mpc2",
         "random_only_cross_delta_sigma_msun_mpc2", "cross_jackknife_error_msun_mpc2",
         "jackknife_error_cross_msun_mpc2"),
    )
    for ax, (title, stack_value, placebo_value, stack_err, placebo_err) in zip(axes, series, strict=True):
        y_stack = np.asarray([row[stack_value] for row in subtracted]) / SCALE
        e_stack = np.asarray([row[stack_err] for row in subtracted]) / SCALE
        y_random = np.asarray([row[placebo_value] for row in random]) / SCALE
        e_random = np.asarray([row[placebo_err] for row in random]) / SCALE
        ax.axhline(0.0, color="#555555", linewidth=0.9, linestyle="--", zorder=1)
        ax.errorbar(radius, y_stack, yerr=e_stack, marker="o", color="#215a9a",
                    markersize=4.5, linewidth=1.7, capsize=2.4,
                    label="Cluster $-$ matched random", zorder=3)
        ax.errorbar(radius, y_random, yerr=e_random, marker="s", color="#c75027",
                    markersize=4.2, linewidth=1.4, capsize=2.4,
                    label="Random centers only", zorder=2)
        ax.set_xscale("log")
        ax.set_xlim(0.09, 3.25)
        ax.set_ylabel(r"$\Delta\Sigma$ ($10^{13}$ M$_\odot$ Mpc$^{-2}$)")
        ax.text(0.025, 0.94, title, transform=ax.transAxes, ha="left", va="top",
                fontsize=11, fontweight="semibold")
        ax.grid(True, which="major", color="#d9dde3", linewidth=0.7)
        ax.grid(True, which="minor", axis="x", color="#eceef1", linewidth=0.5)
        ax.tick_params(which="both", direction="in", top=True, right=True)
        ax.legend(loc="upper right", frameon=True, framealpha=0.9,
                  edgecolor="#b7bec8", borderpad=0.35, handlelength=1.6)
    axes[-1].set_xlabel("Projected radius $R$ (proper Mpc)")
    fig.text(0.99, 0.012,
             "Exploratory; 5° equal-area block jackknife. Random null does not validate cluster signal.",
             ha="right", va="bottom", fontsize=8, color="#4b5563")
    fig.subplots_adjust(left=0.13, right=0.98, top=0.98, bottom=0.09)
    fig.savefig(output_png, dpi=220, facecolor="white")
    fig.savefig(output_pdf, facecolor="white")
    plt.close(fig)
    manifest = {
        "schema": "astrafek-wl-random-placebo-figure-v1",
        "status": "complete_exploratory_figure",
        "science_measurement": False,
        "input_summaries": {str(STACK): sha256(STACK), str(PLACEBO): sha256(PLACEBO)},
        "output_files": {str(output_png): sha256(output_png), str(output_pdf): sha256(output_pdf)},
        "random_centers": placebo["n_random_centers"],
        "cluster_lenses": stack["input"]["matched_sample"]["n_common_shape_support_lenses"],
        "caption": "Cluster-minus-matched-random DeltaSigma and the random-only placebo profile, both with their saved 5-degree spatial-block jackknife errors.",
        "limitations": [
            "Random-only null consistency does not establish validity of the cluster-lens signal.",
            "The equal-area blocks are not validated contiguous HSC survey regions.",
            "All profiles remain exploratory and science_measurement is false.",
        ],
    }
    MANIFEST.write_text(json.dumps(manifest, ensure_ascii=False, indent=2) + "\n",
                        encoding="utf-8")
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--png", type=Path, default=PNG)
    parser.add_argument("--pdf", type=Path, default=PDF)
    args = parser.parse_args()
    print(json.dumps(plot(args.png, args.pdf), ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
