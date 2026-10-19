#!/usr/bin/env python3
"""Plot HSC Y3 source-bin profiles split by lens redshift."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


EXPECTED_SCHEMA = "astrafek-hsc-y3-lens-redshift-tomographic-jackknife-v1"


def load_profiles(summary_path: Path):
    result = json.loads(Path(summary_path).read_text(encoding="utf-8"))
    if result.get("schema") != EXPECTED_SCHEMA or result.get("status") != "exploratory":
        raise ValueError("input is not a completed Astrafek lens-redshift tomography audit")
    if result.get("science_measurement") is not False:
        raise ValueError("lens-redshift audit has an unexpected science gate")
    split = result.get("lens_redshift_split")
    if not isinstance(split, dict) or split.get("n_lenses_lower_z", 0) < 3 or split.get("n_lenses_upper_z", 0) < 3:
        raise ValueError("lens-redshift split metadata is missing or incomplete")
    records = result.get("profile")
    if not isinstance(records, list) or len(records) != 56:
        raise ValueError("expected 2 lens-redshift cohorts by 4 source bins by 7 radial bins")
    profiles = {}
    for cohort in ("lower_z", "upper_z"):
        profiles[cohort] = {}
        for source_bin in range(1, 5):
            rows = [row for row in records
                    if row.get("lens_redshift_cohort") == cohort
                    and row.get("source_zbin") == source_bin]
            if len(rows) != 7:
                raise ValueError(f"{cohort} source bin {source_bin} does not have 7 radial bins")
            rows.sort(key=lambda row: float(row["r_min_mpc"]))
            r_min = np.asarray([float(row["r_min_mpc"]) for row in rows])
            r_max = np.asarray([float(row["r_max_mpc"]) for row in rows])
            tangent = np.asarray([float(row["delta_sigma_weak_msun_mpc2"]) for row in rows])
            tangent_error = np.asarray([float(row["delta_sigma_jackknife_error_msun_mpc2"]) for row in rows])
            cross = np.asarray([float(row["cross_delta_sigma_weak_msun_mpc2"]) for row in rows])
            cross_error = np.asarray([float(row["cross_jackknife_error_msun_mpc2"]) for row in rows])
            arrays = (r_min, r_max, tangent, tangent_error, cross, cross_error)
            if any(not np.all(np.isfinite(values)) for values in arrays):
                raise ValueError(f"{cohort} source bin {source_bin} has non-finite values")
            if np.any(r_min <= 0) or np.any(r_max <= r_min) or np.any(tangent_error <= 0) or np.any(cross_error <= 0):
                raise ValueError(f"{cohort} source bin {source_bin} has invalid radii or errors")
            profiles[cohort][source_bin] = {
                "r_min": r_min, "r_max": r_max,
                "tangent": tangent, "tangent_error": tangent_error,
                "cross": cross, "cross_error": cross_error,
            }
    return result, profiles


def plot(summary_path: Path, output_path: Path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    if output_path.exists():
        raise FileExistsError(f"refusing to overwrite figure: {output_path}")
    result, profiles = load_profiles(summary_path)
    split = result["lens_redshift_split"]
    colors = ("#4c78a8", "#f58518", "#54a24b", "#e45756")
    markers = ("o", "s", "^", "D")
    fig, axes = plt.subplots(2, 2, figsize=(11, 8), sharex="col", constrained_layout=True)
    for cohort_index, cohort in enumerate(("lower_z", "upper_z")):
        for source_bin, color, marker in zip(range(1, 5), colors, markers):
            row = profiles[cohort][source_bin]
            x = np.sqrt(row["r_min"] * row["r_max"])
            xerr = np.vstack((x - row["r_min"], row["r_max"] - x))
            for axis, key, error_key in (
                (axes[cohort_index, 0], "tangent", "tangent_error"),
                (axes[cohort_index, 1], "cross", "cross_error"),
            ):
                axis.errorbar(
                    x, row[key] / 1e12, xerr=xerr, yerr=row[error_key] / 1e12,
                    color=color, marker=marker, linestyle="-", linewidth=1.1,
                    markersize=4, capsize=2, label=f"source bin {source_bin}",
                )
    for column, ylabel in enumerate((
            r"Tangential $\Delta\Sigma$ [$10^{12}\,M_\odot\,\mathrm{Mpc}^{-2}$]",
            r"Cross $\Delta\Sigma$ [$10^{12}\,M_\odot\,\mathrm{Mpc}^{-2}$]")):
        for row in range(2):
            axes[row, column].set_xscale("log")
            axes[row, column].set_ylabel(ylabel)
            axes[row, column].grid(True, which="both", alpha=0.22)
            axes[row, column].legend(frameon=False, ncol=2, fontsize=8)
    for column in range(2):
        axes[1, column].set_xlabel("Proper projected radius [Mpc]")
        axes[1, column].axhline(0.0, color="black", linewidth=0.8, alpha=0.7)
    axes[0, 0].set_title(
        f"Lower-z lenses (N={split['n_lenses_lower_z']}; z≤{split['lower_z_max']:.4f})"
    )
    axes[0, 1].set_title("Lower-z lenses")
    axes[1, 0].set_title(
        f"Higher-z lenses (N={split['n_lenses_upper_z']}; z≥{split['upper_z_min']:.4f})"
    )
    axes[1, 1].set_title("Higher-z lenses")
    fig.suptitle(
        "HSC Y3 source-bin consistency by lens redshift\n"
        "Cluster-only; random-lens subtraction not applied; diagnostic, not a detection",
        fontsize=12,
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=200, metadata={"Title": "HSC Y3 lens-redshift tomography profiles"})
    plt.close(fig)
    return output_path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    print(f"wrote={plot(args.summary.resolve(), args.output.resolve())}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
