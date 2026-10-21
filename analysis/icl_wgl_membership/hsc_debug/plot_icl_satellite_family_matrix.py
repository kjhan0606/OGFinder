#!/usr/bin/env python3
"""Plot the full exploratory ICL–satellite radial hypothesis family."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


EXPECTED_ESTIMATORS = ("equal_satellite", "camira_weighted")
EXPECTED_FILTERS = ("i", "i2")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _read_tsv(path: Path):
    with path.open(encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        if not reader.fieldnames:
            raise ValueError(f"missing TSV header: {path}")
        return reader.fieldnames, list(reader)


def load_family(summary_path: Path, association_path: Path, support_path: Path):
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    if summary.get("schema") != "astrafek-icl-satellite-weighted-radial-matrix-familywise-v1":
        raise ValueError("unexpected family summary schema")
    if summary.get("science_measurement") is not False:
        raise ValueError("source family must remain exploratory")
    for path in (association_path, support_path):
        expected = summary.get("output_sha256", {}).get(str(path))
        if not expected or sha256(path) != expected:
            raise ValueError(f"input checksum mismatch: {path}")

    association_fields, association = _read_tsv(association_path)
    support_fields, support = _read_tsv(support_path)
    required_assoc = {
        "satellite_estimator", "filter", "icl_radial_bin_kpc", "satellite_radial_bin_kpc",
        "status", "n_clusters", "n_spatial_blocks", "partial_pearson_r_z_richness",
        "p_holm_all_eligible_cells", "p_maxT_single_step_complete_null",
    }
    required_support = {
        "estimator", "filter", "icl_radial_bin_kpc", "satellite_radial_bin_kpc",
        "n_clusters", "status",
    }
    if not required_assoc.issubset(association_fields) or not required_support.issubset(support_fields):
        raise ValueError("family TSV schema is incomplete")
    expected_n = int(summary.get("n_expected_cells", 0))
    if expected_n != 168 or len(association) != expected_n or len(support) != expected_n:
        raise ValueError("expected exactly 168 association and support rows")

    def key(row, estimator_key):
        return (row[estimator_key], row["filter"], row["icl_radial_bin_kpc"],
                row["satellite_radial_bin_kpc"])

    associations = {key(row, "satellite_estimator"): row for row in association}
    supports = {key(row, "estimator"): row for row in support}
    if len(associations) != expected_n or len(supports) != expected_n:
        raise ValueError("duplicate family cell")
    if set(associations) != set(supports):
        raise ValueError("association and support cells do not match")
    estimators = set(EXPECTED_ESTIMATORS)
    if {item[0] for item in associations} != estimators:
        raise ValueError("unexpected estimator family")
    if {item[1] for item in associations} != set(EXPECTED_FILTERS):
        raise ValueError("unexpected filter family")
    if any(row["status"] != "exploratory_inference" for row in association):
        raise ValueError("an association cell is not marked exploratory")
    for ident, row in supports.items():
        assoc = associations[ident]
        if int(row["n_clusters"]) != int(assoc["n_clusters"]):
            raise ValueError("association/support cluster count mismatch")
    icl_bins = sorted({item[2] for item in associations}, key=lambda x: int(x.split("-")[0]))
    satellite_bins = sorted({item[3] for item in associations}, key=lambda x: int(x.split("-")[0]))
    if len(icl_bins) != 6 or len(satellite_bins) != 7:
        raise ValueError("expected the complete 6 by 7 radial grid")
    return summary, associations, supports, icl_bins, satellite_bins


def plot_family(summary_path: Path, association_path: Path, support_path: Path,
                output_dir: Path):
    if output_dir.exists():
        raise FileExistsError(f"refusing to overwrite output directory: {output_dir}")
    summary, associations, supports, icl_bins, satellite_bins = load_family(
        summary_path, association_path, support_path)
    output_dir.mkdir(parents=True, exist_ok=False)

    fig, axes = plt.subplots(2, 2, figsize=(16, 11), sharex=True, sharey=True,
                             constrained_layout=True)
    cmap = plt.get_cmap("RdBu_r").copy()
    cmap.set_bad("#e6e6e6")
    image = None
    significant = 0
    panel_counts = {}
    for row_index, estimator in enumerate(EXPECTED_ESTIMATORS):
        for col_index, filter_name in enumerate(EXPECTED_FILTERS):
            ax = axes[row_index, col_index]
            matrix = np.full((len(satellite_bins), len(icl_bins)), np.nan)
            panel_significant = 0
            for iy, sat_bin in enumerate(satellite_bins):
                for ix, icl_bin in enumerate(icl_bins):
                    item = associations[(estimator, filter_name, icl_bin, sat_bin)]
                    r = float(item["partial_pearson_r_z_richness"])
                    n = int(item["n_clusters"])
                    ph = float(item["p_holm_all_eligible_cells"])
                    pm = float(item["p_maxT_single_step_complete_null"])
                    matrix[iy, ix] = r
                    supported = supports[(estimator, filter_name, icl_bin, sat_bin)]["status"] == "eligible"
                    if not supported:
                        matrix[iy, ix] = np.nan
                    if ph < 0.05 or pm < 0.05:
                        significant += 1
                        panel_significant += 1
                    color = "white" if abs(r) >= 0.34 else "#17212b"
                    marker = "*" if ph < 0.05 or pm < 0.05 else ""
                    ax.text(ix, iy, f"{r:+.2f}{marker}\nn={n}", ha="center", va="center",
                            fontsize=8, color=color)
            image = ax.imshow(matrix, origin="upper", aspect="auto", cmap=cmap,
                              vmin=-0.5, vmax=0.5, interpolation="nearest")
            estimator_label = "Equal satellite weights" if estimator == "equal_satellite" else "CAMIRA member weights"
            ax.set_title(f"{estimator_label} · {filter_name} band\nfamily-adjusted cells: {panel_significant}/42",
                         fontsize=11, pad=9)
            ax.set_xticks(range(len(icl_bins)), icl_bins, rotation=35, ha="right")
            ax.set_yticks(range(len(satellite_bins)), satellite_bins)
            ax.set_xlabel("ICL annulus (kpc)")
            ax.set_ylabel("Satellite annulus (kpc)")
            ax.set_xticks(np.arange(-.5, len(icl_bins), 1), minor=True)
            ax.set_yticks(np.arange(-.5, len(satellite_bins), 1), minor=True)
            ax.grid(which="minor", color="white", linewidth=1.3)
            ax.tick_params(which="minor", bottom=False, left=False)

    fig.suptitle("Exploratory ICL–satellite partial correlations across the full radial family\n"
                 "adjusted for redshift and richness; annotation gives r and cluster count",
                 fontsize=15, fontweight="bold")
    colorbar = fig.colorbar(image, ax=axes, shrink=0.83, pad=0.02)
    colorbar.set_label("Partial Pearson correlation r")
    fig.text(0.5, -0.01,
             f"168-cell family: {significant} cells pass Holm or max-T p < 0.05. "
             "All panels are exploratory; PSF/sky and selection systematics remain.",
             ha="center", fontsize=10)
    png = output_dir / "icl_satellite_radial_family.png"
    pdf = output_dir / "icl_satellite_radial_family.pdf"
    fig.savefig(png, dpi=180, bbox_inches="tight")
    fig.savefig(pdf, bbox_inches="tight")
    plt.close(fig)

    result = {
        "schema": "astrafek-icl-satellite-family-plot-v1",
        "status": "exploratory_visualization",
        "science_measurement": False,
        "n_hypotheses": len(associations),
        "n_cells_with_holm_or_maxT_p_below_0_05": significant,
        "inputs": {
            str(summary_path): sha256(summary_path),
            str(association_path): sha256(association_path),
            str(support_path): sha256(support_path),
        },
        "outputs": {str(png): sha256(png), str(pdf): sha256(pdf)},
        "family_summary_status": summary.get("status"),
    }
    (output_dir / "plot_manifest.json").write_text(
        json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--summary", type=Path, required=True)
    parser.add_argument("--association", type=Path, required=True)
    parser.add_argument("--support", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args(argv)
    result = plot_family(args.summary, args.association, args.support, args.output_dir)
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
