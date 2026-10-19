#!/usr/bin/env python3
"""Create broad-annulus summaries and diagnostic plots for a verified ICL stack."""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
from matplotlib import pyplot as plt
from astropy.io import fits
import numpy as np


REGIONS_KPC = ((0, 50), (50, 100), (100, 200), (200, 300), (300, 400), (400, 500))


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def broad_pair_profiles(rows: list[dict[str, str]]) -> list[dict]:
    output = []
    for row in rows:
        bins = json.loads(row["profile_json"])
        band = row["filter"]
        for lo, hi in REGIONS_KPC:
            selected = [b for b in bins if b["r_min_kpc"] >= lo and b["r_max_kpc"] <= hi]
            n_pixels = sum(int(b["n_valid_pixels"]) for b in selected)
            geom = sum(int(b["geometric_pixels"]) for b in selected)
            record = {"target_id": row["target_id"], "control_id": row["control_id"],
                      "filter": band, "r_min_kpc": lo, "r_max_kpc": hi,
                      "n_valid_pixels": n_pixels, "geometric_pixels": geom,
                      "common_coverage_fraction": n_pixels / geom if geom else 0.0}
            for source_key, output_key in (("target_mean_njy_arcsec2", "target"),
                                           ("control_mean_njy_arcsec2", "control"),
                                           ("delta_mean_njy_arcsec2", "delta")):
                weights = np.asarray([int(b["n_valid_pixels"]) for b in selected], dtype=float)
                values = np.asarray([np.nan if b[source_key] is None else float(b[source_key])
                                     for b in selected], dtype=float)
                usable = (weights > 0) & np.isfinite(values)
                record[f"{output_key}_njy_arcsec2"] = (
                    float(np.average(values[usable], weights=weights[usable])) if usable.any() else None)
            output.append(record)
    return output


def aggregate_broad_profiles(pair_rows: list[dict]) -> list[dict]:
    output = []
    for band in ("i", "i2"):
        for lo, hi in REGIONS_KPC:
            selected = [r for r in pair_rows if r["filter"] == band
                        and r["r_min_kpc"] == lo and r["r_max_kpc"] == hi
                        and r["delta_njy_arcsec2"] is not None]
            row = {"filter": band, "r_min_kpc": lo, "r_max_kpc": hi,
                   "n_pairs": len(selected),
                   "median_pair_common_coverage": float(np.median(
                       [r["common_coverage_fraction"] for r in selected])) if selected else None}
            for key in ("target", "control", "delta"):
                values = np.asarray([r[f"{key}_njy_arcsec2"] for r in selected], dtype=float)
                row[f"mean_{key}_njy_arcsec2"] = float(values.mean()) if values.size else None
                row[f"pair_sd_{key}_njy_arcsec2"] = float(values.std(ddof=1)) if values.size > 1 else None
                row[f"pair_sem_{key}_njy_arcsec2"] = float(values.std(ddof=1) / math.sqrt(values.size)) if values.size > 1 else None
            sem = row["pair_sem_delta_njy_arcsec2"]
            row["descriptive_delta_over_pair_sem"] = (
                row["mean_delta_njy_arcsec2"] / sem if sem and np.isfinite(sem) else None)
            output.append(row)
    return output


def plot_band(work: Path, band: str, profile_rows: list[dict], out_path: Path) -> dict:
    product = work / f"icl_target_control_stack_{band}.fits"
    with fits.open(product, memmap=False, checksum=True) as hdus:
        target = np.asarray(hdus["TARGET"].data, dtype=float)
        control = np.asarray(hdus["CONTROL"].data, dtype=float)
        delta = np.asarray(hdus["DELTA"].data, dtype=float)
        n_pair = np.asarray(hdus["NPAIR"].data, dtype=int)
        kpc_max = float(hdus[0].header["KPCMAX"])
        kpc = np.linspace(-kpc_max, kpc_max, target.shape[0], dtype=float)
    selected = [row for row in profile_rows if row["filter"] == band]
    radius = np.asarray([(float(r["r_min_kpc"]) + float(r["r_max_kpc"])) / 2 for r in selected])
    means = np.asarray([np.nan if r["mean_delta_njy_arcsec2"] is None
                        else float(r["mean_delta_njy_arcsec2"]) for r in selected])
    sem = np.asarray([np.nan if r["pair_sem_delta_njy_arcsec2"] is None
                      else float(r["pair_sem_delta_njy_arcsec2"]) for r in selected])
    coverage = np.asarray([np.nan if r["median_pair_common_coverage"] is None
                           else float(r["median_pair_common_coverage"]) for r in selected])
    extent = [float(kpc.min()), float(kpc.max()), float(kpc.min()), float(kpc.max())]
    p995 = float(np.nanpercentile(np.abs(delta), 99.5))
    if not np.isfinite(p995) or p995 <= 0:
        p995 = 1.0
    fig, axes = plt.subplots(2, 2, figsize=(11, 9), constrained_layout=True)
    for ax, data, title in ((axes[0, 0], target, "Target mean"),
                            (axes[0, 1], control, "Matched control mean")):
        finite = data[np.isfinite(data)]
        lim = float(np.nanpercentile(finite, 99.5)) if finite.size else 1.0
        im = ax.imshow(data, origin="lower", extent=extent, cmap="magma", vmin=0, vmax=lim)
        ax.set_title(f"{title} (nJy/arcsec²)")
        ax.set_xlabel("proper kpc"); ax.set_ylabel("proper kpc")
        fig.colorbar(im, ax=ax, shrink=.8)
    im = axes[1, 0].imshow(delta, origin="lower", extent=extent, cmap="coolwarm",
                           vmin=-p995, vmax=p995)
    axes[1, 0].set_title("Target − matched control (nJy/arcsec²)")
    axes[1, 0].set_xlabel("proper kpc"); axes[1, 0].set_ylabel("proper kpc")
    fig.colorbar(im, ax=axes[1, 0], shrink=.8)
    imn = axes[1, 1].imshow(n_pair, origin="lower", extent=extent, cmap="viridis")
    axes[1, 1].set_title("Pair count per grid pixel")
    axes[1, 1].set_xlabel("proper kpc"); axes[1, 1].set_ylabel("proper kpc")
    fig.colorbar(imn, ax=axes[1, 1], shrink=.8)
    for ax in axes.ravel():
        ax.set_aspect("equal")
    fig.suptitle(f"Exploratory HSC-{band} ICL stack; common target/control support")
    fig.savefig(out_path, dpi=160)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 5), constrained_layout=True)
    ax.axhline(0.0, color="black", linewidth=.8)
    ax.plot(radius, means, marker="o", color="#215f9a", label="equal-pair mean")
    ax.fill_between(radius, means - sem, means + sem, color="#4e8ec5", alpha=.25,
                    label="between-pair SEM")
    ax.set_xlabel("projected radius (proper kpc)")
    ax.set_ylabel("target − matched control (nJy/arcsec²)")
    ax.set_title(f"HSC-{band}: paired radial profile; n and coverage vary with radius")
    ax.legend(loc="best")
    twin = ax.twinx()
    twin.plot(radius, coverage, color="#b56b22", linestyle="--", alpha=.8,
              label="median common-area coverage")
    twin.set_ylabel("median pair common-area coverage")
    fig.savefig(out_path.with_name(out_path.stem + "_radial.png"), dpi=160)
    plt.close(fig)
    central = float(delta[delta.shape[0] // 2, delta.shape[1] // 2])
    return {"map_figure": str(out_path), "radial_figure": str(out_path.with_name(out_path.stem + "_radial.png")),
            "central_pixel_delta_njy_arcsec2": central,
            "delta_percentile_99_5_abs_njy_arcsec2": p995,
            "pair_count_center": int(n_pair[n_pair.shape[0] // 2, n_pair.shape[1] // 2]),
            "pair_count_range": [int(n_pair.min()), int(n_pair.max())]}


def evaluate(work: Path) -> dict:
    verification = read_json(work / "independent_verification.json")
    if verification.get("status") != "verified":
        raise ValueError("independent stack verification must pass before evaluation")
    pairs = read_tsv(work / "paired_radial_profiles.tsv")
    pair_profiles = broad_pair_profiles(pairs)
    aggregates = aggregate_broad_profiles(pair_profiles)
    pair_output = work / "broad_radial_pair_profiles.tsv"
    aggregate_output = work / "broad_radial_profiles.tsv"
    for path, rows in ((pair_output, pair_profiles), (aggregate_output, aggregates)):
        if path.exists():
            raise FileExistsError(f"refusing to overwrite {path}")
        with path.open("x", encoding="utf-8", newline="") as stream:
            fields = list(rows[0])
            writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", extrasaction="raise")
            writer.writeheader(); writer.writerows(rows)
    figures = {band: plot_band(work, band, aggregates,
                               work / f"icl_stack_diagnostic_{band}.png")
               for band in ("i", "i2")}
    report = {"schema": "astrafek-hsc-icl-2d-stack-evaluation-v1",
              "status": "completed_exploratory", "verification_sha256": verification["summary_sha256"],
              "n_pairs": len(pairs), "n_by_filter": {band: sum(row["filter"] == band for row in pairs)
                                                       for band in ("i", "i2")},
              "annuli_kpc": [list(x) for x in REGIONS_KPC],
              "broad_radial_profiles": str(aggregate_output),
              "broad_radial_pair_profiles": str(pair_output), "figures": figures,
              "uncertainty_warning": "pair SEM describes pair-to-pair scatter only; it is not a full covariance-aware confidence interval or detection significance.",
              "science_measurement": False}
    report_path = work / "evaluation.json"
    if report_path.exists():
        raise FileExistsError(f"refusing to overwrite {report_path}")
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))
    return report


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--work", type=Path, required=True)
    args = parser.parse_args(argv)
    evaluate(args.work.resolve())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
