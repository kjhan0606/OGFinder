#!/usr/bin/env python3
"""Plot the audited HSC Y3 tomographic DeltaSigma profiles."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


EXPECTED_SCHEMA = "astrafek-hsc-y3-tomographic-consistency-jackknife-v1"


def load_profile(summary_path: Path):
    result = json.loads(Path(summary_path).read_text(encoding="utf-8"))
    if result.get("schema") != EXPECTED_SCHEMA or result.get("status") != "exploratory":
        raise ValueError("input is not a completed Astrafek tomographic consistency audit")
    if result.get("science_measurement") is not False:
        raise ValueError("tomographic audit has an unexpected science gate")
    profiles = result.get("profile")
    if not isinstance(profiles, list) or len(profiles) != 28:
        raise ValueError("expected 4 source bins by 7 radial bins")
    arrays = {}
    for source_bin in range(1, 5):
        rows = [row for row in profiles if row.get("source_zbin") == source_bin]
        if len(rows) != 7:
            raise ValueError(f"source bin {source_bin} does not have 7 radial bins")
        rows.sort(key=lambda row: float(row["r_min_mpc"]))
        radius = np.array([
            float(row["r_min_mpc"]) for row in rows
        ]), np.array([
            float(row["r_max_mpc"]) for row in rows
        ])
        tangent = np.array([float(row["delta_sigma_weak_msun_mpc2"]) for row in rows])
        tangent_error = np.array([float(row["delta_sigma_jackknife_error_msun_mpc2"]) for row in rows])
        cross = np.array([float(row["cross_delta_sigma_weak_msun_mpc2"]) for row in rows])
        cross_error = np.array([float(row["cross_jackknife_error_msun_mpc2"]) for row in rows])
        if any(not np.all(np.isfinite(values)) for values in (*radius, tangent, tangent_error, cross, cross_error)):
            raise ValueError(f"source bin {source_bin} contains non-finite profile values")
        if np.any(radius[0] <= 0) or np.any(radius[1] <= radius[0]) or np.any(tangent_error <= 0) or np.any(cross_error <= 0):
            raise ValueError(f"source bin {source_bin} contains invalid radii or errors")
        arrays[source_bin] = {
            "r_min": radius[0], "r_max": radius[1],
            "tangent": tangent, "tangent_error": tangent_error,
            "cross": cross, "cross_error": cross_error,
        }
    return arrays


def plot(summary_path: Path, output_path: Path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    if output_path.exists():
        raise FileExistsError(f"refusing to overwrite figure: {output_path}")
    profiles = load_profile(summary_path)
    colors = ("#4c78a8", "#f58518", "#54a24b", "#e45756")
    markers = ("o", "s", "^", "D")
    fig, axes = plt.subplots(2, 1, figsize=(8.4, 8.2), sharex=True, constrained_layout=True)
    for source_bin, color, marker in zip(range(1, 5), colors, markers):
        row = profiles[source_bin]
        x = np.sqrt(row["r_min"] * row["r_max"])
        xerr = np.vstack((x - row["r_min"], row["r_max"] - x))
        axes[0].errorbar(
            x, row["tangent"] / 1e12, xerr=xerr,
            yerr=row["tangent_error"] / 1e12, color=color, marker=marker,
            linestyle="-", linewidth=1.2, markersize=4, capsize=2,
            label=f"HSC Y3 source bin {source_bin}",
        )
        axes[1].errorbar(
            x, row["cross"] / 1e12, xerr=xerr,
            yerr=row["cross_error"] / 1e12, color=color, marker=marker,
            linestyle="-", linewidth=1.0, markersize=4, capsize=2,
            label=f"source bin {source_bin}",
        )
    axes[0].set_xscale("log")
    axes[0].set_ylabel(r"Tangential $\Delta\Sigma$ [$10^{12}\,M_\odot\,\mathrm{Mpc}^{-2}$]")
    axes[1].set_xscale("log")
    axes[1].axhline(0.0, color="black", linewidth=0.8, alpha=0.7)
    axes[1].set_ylabel(r"Cross $\Delta\Sigma$ [$10^{12}\,M_\odot\,\mathrm{Mpc}^{-2}$]")
    axes[1].set_xlabel("Proper projected radius [Mpc]")
    for ax in axes:
        ax.grid(True, which="both", alpha=0.22)
        ax.legend(frameon=False, ncol=2, fontsize=8)
    fig.suptitle("HSC Y3 source-bin lensing profiles\nCluster-only; random-lens subtraction not applied", fontsize=12)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=200, metadata={"Title": "HSC Y3 tomographic weak-lensing profiles"})
    plt.close(fig)
    return output_path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    output = plot(args.summary.resolve(), args.output.resolve())
    print(f"wrote={output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
