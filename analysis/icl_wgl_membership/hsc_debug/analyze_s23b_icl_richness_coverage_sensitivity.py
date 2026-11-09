#!/usr/bin/env python3
"""Sensitivity of the S23B richness relation to common-mask area support.

This diagnostic compares the pinned primary richness model with a model that
also adjusts for the number of common target/control pixels in the inner
annulus. It uses paired block-bootstrap draws for the coefficient difference.
Coverage adjustment is a sensitivity, not a preferred causal estimand.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
RICHNESS_DIR = DATA / "s23b_icl_richness_relation_2026-11-07_01"
RICHNESS_SUMMARY = RICHNESS_DIR / "summary.json"
RICHNESS_PAIRS = RICHNESS_DIR / "pair_richness_contrasts.tsv"
RELATIVE_DIR = DATA / "s23b_icl_relative_radial_contrasts_2026-11-07_02"
RELATIVE_SUMMARY = RELATIVE_DIR / "summary.json"
RELATIVE_PAIRS = RELATIVE_DIR / "paired_relative_contrasts.tsv"
SPATIAL_SUMMARY = DATA / "s23b_icl_spatial_blocks_2026-11-07_02/summary.json"
DEFAULT_OUTPUT = DATA / "s23b_icl_richness_coverage_sensitivity_2026-11-07_01"

EXPECTED_HASHES = {
    str(RICHNESS_SUMMARY): "144da86644e8fadfff95d18662af4744d1b4eef0c558e5750d7f7f5200e45f16",
    str(RICHNESS_PAIRS): "69ffbd8f5a39658ef1e5fd3f0ee27872a9c1cb8eca9e5f8a0106f5b63dd3c9fc",
    str(RELATIVE_SUMMARY): "a47c6679ad318abcbee5ee5ab46e1a984976bd01c6f1dcbbe3d6de99b89306a8",
    str(RELATIVE_PAIRS): "c5de5f32b51c6899c8d17fe2b3e4b0a3966d810afd8db5ee27ba649fc0642800",
    str(SPATIAL_SUMMARY): "423bb5e3ad52f3b349dedc74b00a2c05ebb4589e20402abb68b4e99d17199afc",
}
BANDS = ("i", "i2")
CONTRASTS = (
    "50-300_minus_300-800",
    "100-300_minus_300-800",
    "200-300_minus_300-800",
)
WIDTHS_DEG = (2.5, 5.0, 10.0)
PHASES = (0.0, 0.5)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def read_tsv(path: Path) -> list[dict]:
    with Path(path).open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f"refusing to write an empty table: {path}")
    with Path(path).open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t",
                                lineterminator="\n", extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def solve_bootstrap(X_base: np.ndarray, X_cov: np.ndarray, y: np.ndarray,
                    labels: np.ndarray, replicas: int, seed: int):
    groups = np.unique(labels)
    group_idx = {label: index for index, label in enumerate(groups)}
    row_group = np.asarray([group_idx[label] for label in labels], dtype=int)
    if len(groups) < X_cov.shape[1] + 2:
        raise ValueError("too few supported sky blocks for coverage sensitivity")

    beta_base, _, rank_base, _ = np.linalg.lstsq(X_base, y, rcond=None)
    beta_cov, _, rank_cov, _ = np.linalg.lstsq(X_cov, y, rcond=None)
    if rank_base != X_base.shape[1] or rank_cov != X_cov.shape[1]:
        raise ValueError("full-sample model is rank deficient")

    rng = np.random.default_rng(seed)
    block_draws = rng.multinomial(
        len(groups), np.full(len(groups), 1.0 / len(groups)), size=replicas)
    weights = block_draws[:, row_group].astype(float, copy=False)

    def fit_distribution(X: np.ndarray) -> tuple[np.ndarray, int]:
        n_rep, n_par = replicas, X.shape[1]
        gram = np.empty((n_rep, n_par, n_par), dtype=float)
        for left in range(n_par):
            for right in range(left, n_par):
                value = weights @ (X[:, left] * X[:, right])
                gram[:, left, right] = value
                gram[:, right, left] = value
        rhs = np.column_stack([weights @ (X[:, index] * y)
                               for index in range(n_par)])
        condition = np.linalg.cond(gram)
        good = np.isfinite(condition) & (condition < 1e10)
        beta = np.full((n_rep, n_par), np.nan, dtype=float)
        if np.any(good):
            beta[good] = np.linalg.solve(gram[good], rhs[good, :, None])[..., 0]
        return beta, int(good.sum())

    draw_base, _ = fit_distribution(X_base)
    draw_cov, _ = fit_distribution(X_cov)
    good = (np.all(np.isfinite(draw_base), axis=1)
            & np.all(np.isfinite(draw_cov), axis=1))
    if int(good.sum()) < replicas * 0.9:
        raise ValueError(f"only {int(good.sum())}/{replicas} paired bootstrap fits")
    slope_base = draw_base[good, 1]
    slope_cov = draw_cov[good, 1]
    delta = slope_cov - slope_base
    ci_base = np.percentile(slope_base, [2.5, 97.5])
    ci_cov = np.percentile(slope_cov, [2.5, 97.5])
    ci_delta = np.percentile(delta, [2.5, 97.5])
    return {
        "base_slope": float(beta_base[1]),
        "base_ci_low": float(ci_base[0]),
        "base_ci_high": float(ci_base[1]),
        "coverage_adjusted_slope": float(beta_cov[1]),
        "coverage_adjusted_ci_low": float(ci_cov[0]),
        "coverage_adjusted_ci_high": float(ci_cov[1]),
        "slope_change": float(beta_cov[1] - beta_base[1]),
        "paired_slope_change_ci_low": float(ci_delta[0]),
        "paired_slope_change_ci_high": float(ci_delta[1]),
        "n_connected_blocks": int(len(groups)),
        "bootstrap_replicas_valid": int(good.sum()),
    }


def run(output: Path, replicas: int, seed: int) -> dict:
    output = Path(output).resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite coverage-sensitivity output: {output}")
    if replicas < 1000:
        raise ValueError("at least 1,000 paired block-bootstrap replicates are required")
    paths = (RICHNESS_SUMMARY, RICHNESS_PAIRS, RELATIVE_SUMMARY,
             RELATIVE_PAIRS, SPATIAL_SUMMARY)
    for path in paths:
        if not path.is_file():
            raise FileNotFoundError(path)
    inputs = {str(path): sha256(path) for path in paths}
    if inputs != EXPECTED_HASHES:
        raise ValueError("one or more pinned richness/relative-profile inputs changed")

    richness_summary = read_json(RICHNESS_SUMMARY)
    relative_summary = read_json(RELATIVE_SUMMARY)
    spatial = read_json(SPATIAL_SUMMARY)
    if (richness_summary.get("status") != "complete_exploratory_richness_relation"
            or richness_summary.get("science_measurement") is not False
            or relative_summary.get("science_measurement") is not False
            or spatial.get("science_measurement") is not False):
        raise ValueError("parent products are not the pinned exploratory analyses")

    richness_rows = read_tsv(RICHNESS_PAIRS)
    support_rows = read_tsv(RELATIVE_PAIRS)
    support = {}
    for row in support_rows:
        key = (row["target_id"], row["filter"], row["contrast"])
        if key in support:
            raise ValueError(f"duplicate support row: {key}")
        support[key] = int(row["inner_common_pixels"])
    if len(support) != 3066:
        raise ValueError("unexpected relative-profile support row count")

    results = []
    pair_sensitivity_rows = []
    layout_counter = 0
    for width in WIDTHS_DEG:
        for phase in PHASES:
            for band in BANDS:
                layout = [row for row in richness_rows
                          if row["filter"] == band
                          and float(row["block_width_deg"]) == width
                          and float(row["grid_phase_fraction"]) == phase]
                if not layout:
                    raise ValueError(f"missing pair-layout rows for {band}/{width}/{phase}")
                if len({row["target_id"] for row in layout}) != len(layout):
                    raise ValueError("duplicate target within a richness layout")
                expected = next(item for item in spatial["layout_summary"]
                                if item["filter"] == band
                                and float(item["block_width_deg"]) == width
                                and float(item["grid_phase_fraction"]) == phase)
                labels = np.asarray([row["connected_block_id"] for row in layout], dtype=str)
                group_sizes = [int(np.count_nonzero(labels == label))
                               for label in np.unique(labels)]
                if (len(group_sizes) != int(expected["n_connected_blocks"])
                        or max(group_sizes) != int(expected["largest_block_pairs"])):
                    raise ValueError("pair block labels do not reproduce the pinned spatial audit")
                for contrast_index, contrast in enumerate(CONTRASTS):
                    selected = []
                    for row, label in zip(layout, labels):
                        value = row.get(contrast, "")
                        if value in ("", "None"):
                            continue
                        count = support.get((row["target_id"], band, contrast))
                        if count is None or count <= 0:
                            raise ValueError("measured contrast lacks positive common-pixel support")
                        selected.append((row, label, float(value), count))
                    if not selected:
                        raise ValueError("no supported pairs in richness contrast")
                    sublabels = np.asarray([item[1] for item in selected], dtype=str)
                    y = np.asarray([item[2] for item in selected], dtype=float)
                    log_richness = np.asarray([float(item[0]["log10_richness"])
                                               for item in selected])
                    z = np.asarray([float(item[0]["z"]) for item in selected])
                    delta_ebv = np.asarray([float(item[0]["delta_ebv_sfd"])
                                            for item in selected])
                    log_support = np.log10(np.asarray([item[3] for item in selected], dtype=float))
                    base = np.column_stack((
                        np.ones(len(selected)),
                        log_richness - log_richness.mean(),
                        (z - z.mean()) / 0.1,
                        (delta_ebv - delta_ebv.mean()) / 0.01,
                    ))
                    cov_adjusted = np.column_stack((
                        base,
                        log_support - log_support.mean(),
                    ))
                    fitted = solve_bootstrap(
                        base, cov_adjusted, y, sublabels, replicas,
                        seed + layout_counter * 1009 + contrast_index)
                    fitted.update({
                        "filter": band,
                        "contrast": contrast,
                        "block_width_deg": width,
                        "grid_phase_fraction": phase,
                        "n_pairs": len(selected),
                        "support_richness_pearson_r": float(
                            np.corrcoef(log_richness, log_support)[0, 1]),
                        "support_geometric_mean_pixels": float(10 ** log_support.mean()),
                        "bootstrap_replicas_requested": replicas,
                        "adjustment": "add log10 common target-control pixels in inner annulus to base model",
                        "interpretation": "sensitivity only; common-mask area may be a selection variable",
                        "science_measurement": False,
                    })
                    results.append(fitted)
                    for row, label, value, count in selected:
                        pair_sensitivity_rows.append({
                            "target_id": row["target_id"], "control_id": row["control_id"],
                            "filter": band, "contrast": contrast,
                            "block_width_deg": width, "grid_phase_fraction": phase,
                            "connected_block_id": label,
                            "contrast_njy_arcsec2": value,
                            "inner_common_pixels": count,
                            "log10_richness": float(row["log10_richness"]),
                            "science_measurement": False,
                        })
                layout_counter += 1

    output.mkdir(parents=True, exist_ok=False)
    result_path = output / "coverage_adjusted_richness_slopes.tsv"
    pair_path = output / "pair_coverage_sensitivity.tsv"
    write_tsv(result_path, results)
    write_tsv(pair_path, pair_sensitivity_rows)
    primary = [row for row in results if row["block_width_deg"] == 5.0
               and row["grid_phase_fraction"] == 0.0]
    fig_path = output / "coverage_adjusted_richness_slopes.png"
    render(primary, fig_path)
    result = {
        "schema": "astrafek-s23b-icl-richness-coverage-sensitivity-v1",
        "status": "complete_exploratory_coverage_sensitivity",
        "science_measurement": False,
        "network_used": False,
        "image_pixels_read": False,
        "primary_5deg_phase0": primary,
        "limitations": [
            "Adjusting for common-mask area is a sensitivity, not a preferred causal estimator or a completeness correction.",
            "The support pixels come from automatic masks and are dependent on local source crowding and target/control image quality.",
            "This analysis still uses CAMIRA richness, not member coordinates or satellite profiles.",
            "Residual sky, extinction, PSF wings, coadd covariance, centering, and sample-selection uncertainties remain.",
            "Block-bootstrap intervals are exploratory; no multiplicity-adjusted claim is made.",
        ],
        "inputs_sha256": inputs,
        "outputs": {
            result_path.name: {"path": str(result_path), "sha256": sha256(result_path)},
            pair_path.name: {"path": str(pair_path), "sha256": sha256(pair_path)},
            fig_path.name: {"path": str(fig_path), "sha256": sha256(fig_path)},
        },
        "runner": {"path": str(Path(__file__).resolve()),
                   "sha256": sha256(Path(__file__).resolve())},
    }
    (output / "summary.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False) + "\n",
        encoding="utf-8")
    assessment = [
        "# S23B ICL–richness support-area sensitivity — 2026-11-07",
        "",
        "**Status:** exploratory sensitivity; `science_measurement=false`. No network or science-image pixels were used.",
        "",
        "The pinned primary richness model is compared with the same model plus log10 common target/control pixels in the inner annulus. Estimates use paired connected-block bootstrap draws, preserving the same block resample for both models. The extra support term is a diagnostic adjustment only; common area may be a selection variable and this model is not promoted as the preferred causal estimand.",
        "",
        "| Filter | Relative annulus | Pairs / blocks | Primary slope [95% interval] | Coverage-adjusted slope [95% interval] | Paired slope change [95% interval] | corr(log richness, log support) |",
        "|---|---|---:|---:|---:|---:|---:|",
    ]
    for row in primary:
        name = row["contrast"].split("_minus_")[0].replace("-", "–")
        assessment.append(
            f"| {row['filter']} | {name} | {row['n_pairs']} / {row['n_connected_blocks']} | "
            f"{row['base_slope']:.4g} [{row['base_ci_low']:.4g}, {row['base_ci_high']:.4g}] | "
            f"{row['coverage_adjusted_slope']:.4g} [{row['coverage_adjusted_ci_low']:.4g}, {row['coverage_adjusted_ci_high']:.4g}] | "
            f"{row['slope_change']:.4g} [{row['paired_slope_change_ci_low']:.4g}, {row['paired_slope_change_ci_high']:.4g}] | "
            f"{row['support_richness_pearson_r']:.3f} |")
    assessment.extend([
        "",
        "The paired slope-change interval is the 95% block-bootstrap interval for (coverage-adjusted richness slope − primary richness slope), using identical block draws. It is a sensitivity comparison, not a causal effect.",
        "",
        "## Limits",
        "",
        "- Common-mask pixel count is an imperfect support proxy; it does not model which masked sources or background structures remove area.",
        "- Conditioning on coverage can introduce selection bias; compare with the primary model rather than replacing it.",
        "- CAMIRA richness does not reveal satellite positions, luminosity, or radial structure. The member-file download remains approval-gated.",
        "- Keep all results exploratory and `science_measurement=false`.",
        "",
        "## Provenance",
        "",
        *[f"- `{path}` — SHA-256 `{digest}`" for path, digest in inputs.items()],
        f"- Runner: `{Path(__file__).resolve()}` — SHA-256 `{sha256(Path(__file__).resolve())}`",
        "",
    ])
    (output / "assessment.md").write_text("\n".join(assessment), encoding="utf-8")
    return result


def render(rows: list[dict], path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.6), sharey=True)
    x = np.arange(len(CONTRASTS))
    labels = ("50–300", "100–300", "200–300")
    for ax, band in zip(axes, BANDS):
        selected = [next(row for row in rows if row["filter"] == band
                         and row["contrast"] == contrast) for contrast in CONTRASTS]
        for offset, key, low_key, high_key, color, name in (
            (-0.08, "base_slope", "base_ci_low", "base_ci_high", "#1f5a85", "primary"),
            (0.08, "coverage_adjusted_slope", "coverage_adjusted_ci_low",
             "coverage_adjusted_ci_high", "#d17c00", "coverage adjusted"),
        ):
            y = np.asarray([row[key] for row in selected])
            lo = np.asarray([row[low_key] for row in selected])
            hi = np.asarray([row[high_key] for row in selected])
            ax.errorbar(x + offset, y, yerr=np.vstack((y - lo, hi - y)),
                        fmt="o", color=color, capsize=3, linewidth=1.2,
                        label=name)
        ax.axhline(0.0, color="0.35", linewidth=0.9)
        ax.set_xticks(x, labels)
        ax.grid(axis="y", alpha=0.22)
        ax.text(0.04, 0.96, f"HSC-{band}", transform=ax.transAxes,
                ha="left", va="top", fontsize=11)
        ax.set_xlabel("Inner radius range [proper kpc], vs 300–800 kpc")
    axes[0].set_ylabel("Richness slope [nJy arcsec⁻² dex⁻¹]")
    axes[1].legend(loc="upper right", fontsize=9, frameon=True)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--replicas", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=20261110)
    args = parser.parse_args()
    result = run(args.output, args.replicas, args.seed)
    compact = [{k: row[k] for k in (
        "filter", "contrast", "n_pairs", "n_connected_blocks", "base_slope",
        "base_ci_low", "base_ci_high", "coverage_adjusted_slope",
        "coverage_adjusted_ci_low", "coverage_adjusted_ci_high", "slope_change",
        "paired_slope_change_ci_low", "paired_slope_change_ci_high")}
        for row in result["primary_5deg_phase0"]]
    print(json.dumps({"status": result["status"],
                      "science_measurement": result["science_measurement"],
                      "primary_5deg_phase0": compact,
                      "output": str(Path(args.output).resolve())},
                     ensure_ascii=False, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
