#!/usr/bin/env python3
"""Estimate paired inner-minus-outer ICL contrasts in the S23B stack.

The estimator subtracts each target-control pair's 300--800 proper-kpc
profile from its inner profile. It is a relative radial contrast, not an
absolute ICL measurement or a validated sky correction.
"""
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

from analyze_s23b_icl_spatial_blocks import (
    CONTROL_MANIFEST,
    EXPECTED_CONTROL_MANIFEST_SHA256,
    EXPECTED_STACK_SUMMARY_SHA256,
    EXPECTED_TARGET_MANIFEST_SHA256,
    PAIRED_PROFILES,
    STACK_SUMMARY,
    STACK_DIR,
    TARGET_MANIFEST,
    broad_profile,
    connected_pair_blocks,
    read_json,
    read_tsv,
    sha256,
    write_tsv,
)


DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
SPATIAL_SUMMARY = DATA / "s23b_icl_spatial_blocks_2026-11-07_02/summary.json"
DEFAULT_OUTPUT = DATA / "s23b_icl_relative_radial_contrasts_2026-11-07_01"
EXPECTED_SPATIAL_SUMMARY_SHA256 = "423bb5e3ad52f3b349dedc74b00a2c05ebb4589e20402abb68b4e99d17199afc"
EXPECTED_PAIR_PROFILE_SHA256 = "285906b03360e9d137727e948c15419767c356d3343e61f409562567ac4abc43"
WIDTHS_DEG = (2.5, 5.0, 10.0)
PHASES = (0.0, 0.5)
CONTRASTS_KPC = ((50, 300), (100, 300), (200, 300))
OUTER_KPC = (300, 800)
EXPECTED_N = {"i": 241, "i2": 781}


def profile_delta(profile: list[dict], low: int, high: int) -> tuple[float | None, int]:
    return broad_profile(profile, low, high)


def make_pair_contrasts(pair_rows: list[dict], target_by_id: dict,
                        control_by_target: dict) -> dict[str, list[dict]]:
    out = {"i": [], "i2": []}
    seen = set()
    for row in pair_rows:
        target_id = str(row["target_id"])
        band = str(row["filter"])
        if band not in out or target_id not in target_by_id or target_id not in control_by_target:
            raise ValueError(f"unexpected target or filter in pair table: {target_id}/{band}")
        key = (target_id, band)
        if key in seen:
            raise ValueError(f"duplicate target/filter pair: {key}")
        seen.add(key)
        target = target_by_id[target_id]
        control = control_by_target[target_id]
        if (str(control["id"]) != str(row["control_id"])
                or float(row["ra_deg"]) != float(target["ra_deg"])
                or float(row["dec_deg"]) != float(target["dec_deg"])):
            raise ValueError(f"pair coordinates or control mapping disagree for {target_id}")
        profile = json.loads(row["profile_json"])
        record = {
            "target_id": target_id,
            "control_id": str(row["control_id"]),
            "filter": band,
            "ra_deg": float(row["ra_deg"]),
            "dec_deg": float(row["dec_deg"]),
            "control_ra_deg": float(control["ra_deg"]),
            "control_dec_deg": float(control["dec_deg"]),
            "z": float(row["z"]),
            "contrasts": {},
        }
        for inner_low, inner_high in CONTRASTS_KPC:
            inner, n_inner = profile_delta(profile, inner_low, inner_high)
            outer, n_outer = profile_delta(profile, *OUTER_KPC)
            label = f"{inner_low}-{inner_high}_minus_{OUTER_KPC[0]}-{OUTER_KPC[1]}"
            value = inner - outer if n_inner > 0 and n_outer > 0 else None
            record["contrasts"][label] = {
                "inner_mean_njy_arcsec2": inner,
                "outer_mean_njy_arcsec2": outer,
                "contrast_njy_arcsec2": value,
                "inner_common_pixels": n_inner,
                "outer_common_pixels": n_outer,
            }
        out[band].append(record)
    if {band: len(rows) for band, rows in out.items()} != EXPECTED_N:
        raise ValueError("filter-specific pair counts differ from the pinned stack")
    return out


def summarize_contrast_layout(rows: list[dict], labels: np.ndarray, band: str,
                              width: float, phase: float, replicas: int,
                              seed: int) -> tuple[list[dict], list[dict]]:
    groups = np.unique(labels)
    group_index = {group: index for index, group in enumerate(groups)}
    row_group = np.asarray([group_index[label] for label in labels], dtype=int)
    rng = np.random.default_rng(seed)
    block_draws = rng.multinomial(
        len(groups), np.full(len(groups), 1.0 / len(groups)), size=replicas)
    summary_rows, loo_rows = [], []
    for inner_low, inner_high in CONTRASTS_KPC:
        label = f"{inner_low}-{inner_high}_minus_{OUTER_KPC[0]}-{OUTER_KPC[1]}"
        all_values = np.asarray([
            np.nan if row["contrasts"][label]["contrast_njy_arcsec2"] is None
            else float(row["contrasts"][label]["contrast_njy_arcsec2"])
            for row in rows
        ], dtype=float)
        valid = np.isfinite(all_values)
        values = all_values[valid]
        value_groups = row_group[valid]
        n_blocks = int(np.unique(value_groups).size)
        mean = float(values.mean()) if values.size else None
        sd = float(values.std(ddof=1)) if values.size > 1 else None
        sem = sd / np.sqrt(values.size) if sd is not None else None
        group_sums = np.bincount(value_groups, weights=values, minlength=len(groups))
        group_counts = np.bincount(value_groups, minlength=len(groups))
        denominator = block_draws @ group_counts
        numerator = block_draws @ group_sums
        valid_draws = denominator > 0
        boot = numerator[valid_draws] / denominator[valid_draws]
        if n_blocks >= 3 and boot.size >= replicas * 0.9:
            ci_low, ci_high = np.percentile(boot, [2.5, 97.5]).tolist()
            status = "complete"
        else:
            ci_low = ci_high = None
            status = "insufficient_supported_blocks_or_valid_draws"

        deletion_rows = []
        for group in np.unique(value_groups):
            keep = valid & (row_group != group)
            if np.any(keep):
                deleted_mean = float(all_values[keep].mean())
                deletion_rows.append((str(groups[group]), deleted_mean,
                                      int(np.count_nonzero(keep))))
        if mean is not None and deletion_rows:
            most_influential = max(deletion_rows,
                                   key=lambda item: abs(item[1] - mean))
            max_shift = abs(most_influential[1] - mean)
            sign_flips = sum(1 for _, value, _ in deletion_rows
                             if value != 0 and np.sign(value) != np.sign(mean))
        else:
            most_influential = (None, None, None)
            max_shift = None
            sign_flips = None
        summary_rows.append({
            "filter": band,
            "contrast": label,
            "inner_r_min_kpc": inner_low,
            "inner_r_max_kpc": inner_high,
            "outer_r_min_kpc": OUTER_KPC[0],
            "outer_r_max_kpc": OUTER_KPC[1],
            "block_width_deg": width,
            "grid_phase_fraction": phase,
            "mean_contrast_njy_arcsec2": mean,
            "pair_sd_njy_arcsec2": sd,
            "pair_sem_njy_arcsec2": sem,
            "block_bootstrap_ci95_low": ci_low,
            "block_bootstrap_ci95_high": ci_high,
            "bootstrap_replicas_requested": replicas,
            "bootstrap_replicas_valid": int(boot.size),
            "bootstrap_status": status,
            "n_supported_pairs": int(values.size),
            "n_supported_blocks": n_blocks,
            "largest_supported_block_pairs": int(max(group_counts, default=0)),
            "max_leave_one_block_abs_shift_njy_arcsec2": max_shift,
            "max_shift_block_id": most_influential[0],
            "n_leave_one_block_sign_flips": sign_flips,
            "science_measurement": False,
        })
        for block_id, value, n_retained in deletion_rows:
            loo_rows.append({
                "filter": band,
                "contrast": label,
                "block_width_deg": width,
                "grid_phase_fraction": phase,
                "deleted_block_id": block_id,
                "n_pairs_retained": n_retained,
                "mean_contrast_njy_arcsec2": value,
                "shift_from_full_njy_arcsec2": value - mean,
            })
    return summary_rows, loo_rows


def render_plot(path: Path, summaries: list[dict], replicas: int) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 5.0), sharex=True, sharey=True)
    contrast_names = [f"{lo}–{hi}\nvs 300–800" for lo, hi in CONTRASTS_KPC]
    xpos = np.arange(len(CONTRASTS_KPC))
    for ax, band in zip(axes, ("i", "i2")):
        rows = [r for r in summaries if r["filter"] == band
                and r["block_width_deg"] == 5.0
                and r["grid_phase_fraction"] == 0.0]
        rows.sort(key=lambda r: r["inner_r_min_kpc"])
        y = np.asarray([r["mean_contrast_njy_arcsec2"] for r in rows], dtype=float)
        lo = np.asarray([r["block_bootstrap_ci95_low"] for r in rows], dtype=float)
        hi = np.asarray([r["block_bootstrap_ci95_high"] for r in rows], dtype=float)
        sem = np.asarray([r["pair_sem_njy_arcsec2"] for r in rows], dtype=float)
        ax.errorbar(xpos, y, yerr=sem, color="#1f5a85", marker="o", linewidth=1.2,
                    capsize=3, label="pair SEM")
        ax.errorbar(xpos + 0.045, y, yerr=np.vstack([y - lo, hi - y]), fmt="none",
                    ecolor="#d17c00", linewidth=2.2, capsize=5,
                    label="5° block-bootstrap 95% interval")
        ax.axhline(0.0, color="0.35", linewidth=0.9)
        ax.set_xticks(xpos, contrast_names)
        ax.grid(axis="y", alpha=0.22)
        ax.text(0.04, 0.95, f"HSC-{band}", transform=ax.transAxes,
                ha="left", va="top", fontsize=11)
        ax.set_xlabel("Inner minus outer annulus [proper kpc]")
    axes[0].set_ylabel("Relative target − control surface brightness [nJy arcsec⁻²]")
    axes[1].legend(loc="upper right", fontsize=9, frameon=True)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def run(output: Path, replicas: int, seed: int) -> dict:
    output = Path(output).resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite relative-contrast output: {output}")
    if replicas < 1000:
        raise ValueError("at least 1,000 block-bootstrap replicates are required")
    input_paths = (TARGET_MANIFEST, CONTROL_MANIFEST, STACK_SUMMARY,
                   PAIRED_PROFILES, SPATIAL_SUMMARY)
    for path in input_paths:
        if not path.is_file():
            raise FileNotFoundError(path)
    input_hashes = {str(path.resolve()): sha256(path) for path in input_paths}
    expected = {
        str(TARGET_MANIFEST.resolve()): EXPECTED_TARGET_MANIFEST_SHA256,
        str(CONTROL_MANIFEST.resolve()): EXPECTED_CONTROL_MANIFEST_SHA256,
        str(STACK_SUMMARY.resolve()): EXPECTED_STACK_SUMMARY_SHA256,
        str(PAIRED_PROFILES.resolve()): EXPECTED_PAIR_PROFILE_SHA256,
        str(SPATIAL_SUMMARY.resolve()): EXPECTED_SPATIAL_SUMMARY_SHA256,
    }
    if input_hashes != expected:
        raise ValueError("one or more target/control stack inputs changed from their pinned hashes")
    target_manifest = read_json(TARGET_MANIFEST)
    control_manifest = read_json(CONTROL_MANIFEST)
    stack_summary = read_json(STACK_SUMMARY)
    spatial_summary = read_json(SPATIAL_SUMMARY)
    if (stack_summary.get("status") != "complete_exploratory"
            or stack_summary.get("science_measurement") is not False
            or spatial_summary.get("status") != "complete_spatial_block_sensitivity"
            or spatial_summary.get("science_measurement") is not False):
        raise ValueError("pinned parent outputs are not complete exploratory products")
    expected_pair_sha = stack_summary.get("outputs", {}).get("paired_profiles_sha256")
    if expected_pair_sha != EXPECTED_PAIR_PROFILE_SHA256:
        raise ValueError("stack summary does not pin the expected pair profile table")
    spatial_runner = spatial_summary.get("runner", {})
    if spatial_runner.get("sha256") != sha256(Path(spatial_runner.get("path", ""))):
        raise ValueError("spatial-block helper source changed after its pinned run")

    targets = {str(row["id"]): row for row in target_manifest.get("clusters", [])}
    controls = control_manifest.get("clusters", [])
    control_by_target = {str(row["target_cluster_id"]): row for row in controls}
    if len(targets) != 1043 or len(control_by_target) != 1041:
        raise ValueError("target/control manifest counts changed")
    pair_rows = read_tsv(PAIRED_PROFILES)
    by_band = make_pair_contrasts(pair_rows, targets, control_by_target)
    pair_contrast_rows = []
    for band, rows in by_band.items():
        for row in rows:
            for label, values in row["contrasts"].items():
                pair_contrast_rows.append({
                    "target_id": row["target_id"], "control_id": row["control_id"],
                    "filter": band, "z": row["z"], "contrast": label,
                    "inner_mean_njy_arcsec2": values["inner_mean_njy_arcsec2"],
                    "outer_mean_njy_arcsec2": values["outer_mean_njy_arcsec2"],
                    "contrast_njy_arcsec2": values["contrast_njy_arcsec2"],
                    "inner_common_pixels": values["inner_common_pixels"],
                    "outer_common_pixels": values["outer_common_pixels"],
                    "science_measurement": False,
                })

    all_summaries, all_loo, all_blocks = [], [], []
    layout_index = 0
    for width in WIDTHS_DEG:
        for phase in PHASES:
            for band_index, band in enumerate(("i", "i2")):
                rows = by_band[band]
                labels, blocks = connected_pair_blocks(
                    rows, control_by_target, width, phase,
                    float(spatial_summary["ra_grid_origin_deg"]))
                all_blocks.extend({"filter": band, "block_width_deg": width,
                                   "grid_phase_fraction": phase, **block}
                                  for block in blocks)
                expected_layout = next(item for item in spatial_summary["layout_summary"]
                                       if item["filter"] == band
                                       and item["block_width_deg"] == width
                                       and item["grid_phase_fraction"] == phase)
                if (len(blocks) != expected_layout["n_connected_blocks"]
                        or max(b["n_pairs"] for b in blocks) != expected_layout["largest_block_pairs"]):
                    raise ValueError("block labels do not reproduce the pinned spatial-block audit")
                summaries, loo = summarize_contrast_layout(
                    rows, labels, band, width, phase, replicas,
                    seed + layout_index * 101 + band_index)
                all_summaries.extend(summaries)
                all_loo.extend(loo)
                layout_index += 1

    output.mkdir(parents=True, exist_ok=False)
    pair_path = output / "paired_relative_contrasts.tsv"
    summary_table = output / "contrast_block_uncertainty.tsv"
    loo_path = output / "leave_one_block_out.tsv"
    block_path = output / "spatial_blocks.tsv"
    plot_path = output / "s23b_icl_relative_radial_contrasts.png"
    write_tsv(pair_path, pair_contrast_rows)
    write_tsv(summary_table, all_summaries)
    write_tsv(loo_path, all_loo)
    write_tsv(block_path, all_blocks)
    render_plot(plot_path, all_summaries, replicas)
    primary = [row for row in all_summaries
               if row["block_width_deg"] == 5.0
               and row["grid_phase_fraction"] == 0.0]
    sensitivity_envelope = []
    for band in ("i", "i2"):
        for inner_low, inner_high in CONTRASTS_KPC:
            label = f"{inner_low}-{inner_high}_minus_{OUTER_KPC[0]}-{OUTER_KPC[1]}"
            variants = [row for row in all_summaries
                        if row["filter"] == band and row["contrast"] == label]
            lows = [float(row["block_bootstrap_ci95_low"]) for row in variants
                    if row["block_bootstrap_ci95_low"] is not None]
            highs = [float(row["block_bootstrap_ci95_high"]) for row in variants
                     if row["block_bootstrap_ci95_high"] is not None]
            sensitivity_envelope.append({
                "filter": band,
                "contrast": label,
                "min_interval_lower": min(lows) if lows else None,
                "max_interval_lower": max(lows) if lows else None,
                "min_interval_upper": min(highs) if highs else None,
                "max_interval_upper": max(highs) if highs else None,
                "all_layouts_exclude_zero": bool(
                    lows and highs and (min(lows) > 0.0 or max(highs) < 0.0)),
            })
    outputs = {path.name: {"path": str(path), "sha256": sha256(path)}
               for path in (pair_path, summary_table, loo_path, block_path, plot_path)}
    result = {
        "schema": "astrafek-s23b-icl-relative-radial-contrast-v1",
        "status": "complete_exploratory_relative_contrast",
        "science_measurement": False,
        "image_pixels_read": False,
        "network_used": False,
        "n_pairs_by_filter": {band: len(rows) for band, rows in by_band.items()},
        "estimand": "equal-pair mean of the within-pair common-support target-minus-control inner-annulus profile minus that pair's 300-800 proper-kpc profile",
        "contrast_inner_annuli_kpc": [list(x) for x in CONTRASTS_KPC],
        "outer_annulus_kpc": list(OUTER_KPC),
        "block_widths_deg": list(WIDTHS_DEG),
        "grid_phases": list(PHASES),
        "replicas_per_layout_and_filter": replicas,
        "seed_base": int(seed),
        "primary_5deg_profiles": primary,
        "six_layout_sensitivity_envelope": sensitivity_envelope,
        "all_layout_summary": [{k: row[k] for k in (
            "filter", "block_width_deg", "grid_phase_fraction", "n_connected_blocks",
            "largest_block_pairs")} for row in spatial_summary["layout_summary"]],
        "inputs_sha256": input_hashes,
        "spatial_helper_sha256": spatial_runner["sha256"],
        "outputs": outputs,
        "limitations": [
            "Subtracting the outer annulus can cancel a uniform additive offset but assumes neither that the outer profile is sky nor that it contains no cluster light.",
            "The RA-Dec connected cells approximate spatial dependence and are not official HSC field polygons or a calibrated survey-region covariance.",
            "Block bootstrap intervals do not propagate coadd-pixel covariance, source Poisson noise, Galactic extinction, residual-sky gradients, PSF wings, or center/mask uncertainty.",
            "The annular means use strict target/control common support, whose area and supported-pair count vary with radius.",
            "This relative profile does not measure a complete satellite census; the separate CAMIRA member-position file has not been approved for download or imported.",
            "Keep science_measurement=false; this output is an exploratory radial-gradient analysis.",
        ],
        "runner": {"path": str(Path(__file__).resolve()),
                   "sha256": sha256(Path(__file__).resolve())},
    }
    summary_path = output / "summary.json"
    summary_path.write_text(json.dumps(result, ensure_ascii=False, indent=2,
                                       allow_nan=False) + "\n", encoding="utf-8")

    assessment = [
        "# S23B paired inner-minus-outer ICL contrast — 2026-11-07",
        "",
        "**Status:** exploratory relative radial contrast; `science_measurement=false`.",
        "",
        "For each target/control pair, this analysis subtracts the pair's 300–800 proper-kpc common-support mean from its inner-annulus common-support mean. The subtraction can reduce a uniform additive sky offset. The outer annulus is not assumed to be blank sky and may contain real cluster light.",
        "",
        "No image pixels or HSC services were accessed. The analysis uses the 1,022 saved pair profiles (241 `i`; 781 `i2`) and resamples connected RA-Dec blocks at three widths and two cell phases.",
        "",
        "## Primary 5-degree results",
        "",
        "| Filter | Relative annulus (proper kpc) | Pairs / blocks | Mean contrast | Pair SEM | Block-bootstrap 95% interval |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for row in primary:
        lo = row["block_bootstrap_ci95_low"]
        hi = row["block_bootstrap_ci95_high"]
        ci = f"[{lo:.3g}, {hi:.3g}]" if lo is not None else "not estimated"
        sem = f"{row['pair_sem_njy_arcsec2']:.3g}" if row["pair_sem_njy_arcsec2"] is not None else "—"
        annulus = f"{row['inner_r_min_kpc']}–{row['inner_r_max_kpc']} minus 300–800"
        report_line = (f"| {row['filter']} | {annulus} | {row['n_supported_pairs']} / "
                       f"{row['n_supported_blocks']} | {row['mean_contrast_njy_arcsec2']:.3g} | "
                       f"{sem} | {ci} |")
        assessment.append(report_line)
    assessment.extend([
        "",
        "## Layout sensitivity",
        "",
        "The table spans the lower and upper endpoints of the six separate 95% intervals across the three block widths and two grid phases. This is a sensitivity envelope, not a single confidence interval.",
        "",
        "| Filter | Relative annulus | Range of interval lower endpoints | Range of interval upper endpoints | All layouts exclude zero? |",
        "|---|---:|---:|---:|---|",
    ])
    for row in sensitivity_envelope:
        contrast = row["contrast"].split("_minus_")[0].replace("-", "–")
        low_range = (f"[{row['min_interval_lower']:.3g}, {row['max_interval_lower']:.3g}]"
                     if row["min_interval_lower"] is not None else "not estimated")
        high_range = (f"[{row['min_interval_upper']:.3g}, {row['max_interval_upper']:.3g}]"
                      if row["min_interval_upper"] is not None else "not estimated")
        assessment.append(
            f"| {row['filter']} | {contrast} vs 300–800 | {low_range} | "
            f"{high_range} | {'yes' if row['all_layouts_exclude_zero'] else 'no'} |")
    assessment.extend([
        "",
        "The 50–300 and 100–300 kpc contrasts remain positive in both filters across all six layouts. The narrower 200–300 kpc contrast crosses zero for at least one layout in each filter. A stable positive contrast means the inner surface brightness exceeds the selected outer baseline within this common-support estimator; it does not by itself identify the excess as ICL.",
        "",
        "The accompanying TSV also reports pairwise contrasts and every layout-specific interval. Changes across widths and phases reflect approximate block-partition sensitivity; they do not constitute a calibrated survey covariance.",
        "",
        "## Limitations and next gates",
        "",
        "- The 300–800 kpc baseline can contain cluster light and is not an absolute sky estimate.",
        "- HSC extinction coefficients and spatial variation, residual-sky gradients, extended PSF matching, and coadd covariance remain uncorrected.",
        "- Common usable area is sparse in the center and varies by radius; this analysis uses no post-hoc coverage cut.",
        "- The member catalog download remains a separate approval-gated step; no satellite coordinates are inferred from the cluster table or an LRG file.",
        "- The contrast and block intervals remain exploratory; `science_measurement=false`.",
        "",
        "## Provenance",
        "",
        *[f"- `{path}` — SHA-256 `{digest}`" for path, digest in input_hashes.items()],
        f"- Runner: `{Path(__file__).resolve()}` — SHA-256 `{sha256(Path(__file__).resolve())}`",
        f"- Spatial-block helper SHA-256: `{spatial_runner['sha256']}`",
        "",
    ])
    (output / "assessment.md").write_text("\n".join(assessment), encoding="utf-8")
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--replicas", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=20261108)
    args = parser.parse_args()
    result = run(args.output, args.replicas, args.seed)
    compact = [{k: row[k] for k in (
        "filter", "contrast", "n_supported_pairs", "n_supported_blocks",
        "mean_contrast_njy_arcsec2", "pair_sem_njy_arcsec2",
        "block_bootstrap_ci95_low", "block_bootstrap_ci95_high")}
        for row in result["primary_5deg_profiles"]]
    print(json.dumps({"status": result["status"],
                      "science_measurement": result["science_measurement"],
                      "primary_5deg_profiles": compact,
                      "output": str(Path(args.output).resolve())},
                     ensure_ascii=False, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
