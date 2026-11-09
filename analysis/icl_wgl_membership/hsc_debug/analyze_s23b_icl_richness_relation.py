#!/usr/bin/env python3
"""Estimate the S23B ICL radial-contrast relation with CAMIRA richness.

This is an exploratory cluster-level characterization, not a satellite-member
profile. It models the saved paired ICL contrast against log CAMIRA richness,
adjusting for redshift and the residual SFD E(B-V) mismatch between each target
and its matched random control. Uncertainty is resampled by connected sky
blocks over the same six RA-Dec grid layouts used by the pinned spatial audit.
No HSC requests or science-image pixels are used.
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

from analyze_s23b_icl_spatial_blocks import (
    connected_pair_blocks,
    read_json,
    sha256,
    write_tsv,
)


DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
TARGET_MANIFEST = DATA / "s23b_icl_target_image_plan_2026-11-07_01/target_manifest.json"
CONTROL_MANIFEST = DATA / "s23b_icl_control_image_plan_2026-11-07_02/control_manifest.json"
STACK_SUMMARY = DATA / "s23b_icl_target_control_stack_2026-11-07_01/summary.json"
PAIR_PROFILES = DATA / "s23b_icl_target_control_stack_2026-11-07_01/paired_radial_profiles.tsv"
SPATIAL_SUMMARY = DATA / "s23b_icl_spatial_blocks_2026-11-07_02/summary.json"
RELATIVE_SUMMARY = DATA / "s23b_icl_relative_radial_contrasts_2026-11-07_02/summary.json"
DEFAULT_OUTPUT = DATA / "s23b_icl_richness_relation_2026-11-07_01"

EXPECTED_HASHES = {
    str(TARGET_MANIFEST): "e487582ba40cd3e58da67704242db4d7219e381fb08f4b33fed47ce4ff268b44",
    str(CONTROL_MANIFEST): "ae9eb406cfc179d9980031bb67f85f8ee71c4580ba59e97713278823db5c06df",
    str(STACK_SUMMARY): "c8bac8dbdd0a886604c216b2410725eb97883b0cdc1a91fd47e7a7cfcaaaa20b",
    str(PAIR_PROFILES): "285906b03360e9d137727e948c15419767c356d3343e61f409562567ac4abc43",
    str(SPATIAL_SUMMARY): "423bb5e3ad52f3b349dedc74b00a2c05ebb4589e20402abb68b4e99d17199afc",
    str(RELATIVE_SUMMARY): "a47c6679ad318abcbee5ee5ab46e1a984976bd01c6f1dcbbe3d6de99b89306a8",
}
EXPECTED_RELATIVE_TABLE_SHA256 = "c5de5f32b51c6899c8d17fe2b3e4b0a3966d810afd8db5ee27ba649fc0642800"
EXPECTED_BLOCK_HELPER_SHA256 = "d7a21f68988fb782e5bcda147ff9a9040ded1b8c3f6a12ac883598beb2cd0db5"
BLOCK_WIDTHS_DEG = (2.5, 5.0, 10.0)
GRID_PHASES = (0.0, 0.5)
BANDS = ("i", "i2")
CONTRASTS = (
    "50-300_minus_300-800",
    "100-300_minus_300-800",
    "200-300_minus_300-800",
)
EXPECTED_N = {"i": 241, "i2": 781}


def read_tsv(path: Path) -> list[dict]:
    with Path(path).open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def as_float(value, label: str) -> float:
    number = float(value)
    if not math.isfinite(number):
        raise ValueError(f"non-finite {label}")
    return number


def verify_inputs() -> tuple[dict, dict, dict, dict, dict, Path, dict]:
    paths = (TARGET_MANIFEST, CONTROL_MANIFEST, STACK_SUMMARY, PAIR_PROFILES,
             SPATIAL_SUMMARY, RELATIVE_SUMMARY)
    for path in paths:
        if not path.is_file():
            raise FileNotFoundError(path)
    actual = {str(path): sha256(path) for path in paths}
    if actual != EXPECTED_HASHES:
        raise ValueError("one or more pinned S23B inputs changed")

    target = read_json(TARGET_MANIFEST)
    control = read_json(CONTROL_MANIFEST)
    stack = read_json(STACK_SUMMARY)
    spatial = read_json(SPATIAL_SUMMARY)
    relative = read_json(RELATIVE_SUMMARY)
    rel_entry = relative.get("outputs", {}).get("paired_relative_contrasts.tsv", {})
    relative_table = Path(rel_entry.get("path", ""))
    if (not relative_table.is_file()
            or sha256(relative_table) != EXPECTED_RELATIVE_TABLE_SHA256
            or rel_entry.get("sha256") != EXPECTED_RELATIVE_TABLE_SHA256):
        raise ValueError("paired relative-contrast table changed from the pinned run")
    helper = Path(__file__).resolve().with_name("analyze_s23b_icl_spatial_blocks.py")
    if sha256(helper) != EXPECTED_BLOCK_HELPER_SHA256:
        raise ValueError("connected-block helper source changed")
    if (stack.get("status") != "complete_exploratory"
            or stack.get("science_measurement") is not False
            or relative.get("status") != "complete_exploratory_relative_contrast"
            or relative.get("science_measurement") is not False
            or spatial.get("science_measurement") is not False):
        raise ValueError("parent outputs are not the pinned exploratory products")
    return target, control, stack, spatial, relative, relative_table, actual


def make_pair_records(target_manifest: dict, control_manifest: dict,
                      relative_table: Path) -> dict[str, list[dict]]:
    targets = {str(row["id"]): row for row in target_manifest.get("clusters", [])}
    controls = {str(row["target_cluster_id"]): row
                for row in control_manifest.get("clusters", [])}
    profiles = {(str(row["target_id"]), str(row["filter"])): row
                for row in read_tsv(PAIR_PROFILES)}
    if (len(targets) != 1043 or len(controls) != 1041
            or len(profiles) != 1022):
        raise ValueError("pinned target/control/profile counts changed")

    rows_by_key: dict[tuple[str, str], dict] = {}
    for pair in read_tsv(relative_table):
        target_id, band = str(pair["target_id"]), str(pair["filter"])
        if band not in BANDS or target_id not in targets or target_id not in controls:
            raise ValueError(f"unexpected target/filter key {target_id}/{band}")
        key = (target_id, band)
        if key not in rows_by_key:
            target, control = targets[target_id], controls[target_id]
            profile = profiles.get(key)
            if profile is None or str(profile["control_id"]) != str(pair["control_id"]):
                raise ValueError(f"pair-profile link mismatch for {target_id}/{band}")
            if (str(control["id"]) != str(pair["control_id"])
                    or str(control["target_cluster_id"]) != target_id):
                raise ValueError(f"control manifest mapping mismatch for {target_id}")
            if str(pair["science_measurement"]).lower() != "false":
                raise ValueError("input contrast unexpectedly claims science measurement")
            richness = as_float(target.get("richness"), f"richness/{target_id}")
            if richness <= 0:
                raise ValueError(f"non-positive CAMIRA richness for {target_id}")
            target_ebv = as_float(control.get("target_ebv_sfd"), f"target E(B-V)/{target_id}")
            control_ebv = as_float(control.get("candidate_ebv_sfd"), f"control E(B-V)/{target_id}")
            rows_by_key[key] = {
                "target_id": target_id,
                "control_id": str(pair["control_id"]),
                "filter": band,
                "ra_deg": as_float(profile["ra_deg"], f"RA/{target_id}"),
                "dec_deg": as_float(profile["dec_deg"], f"Dec/{target_id}"),
                "control_ra_deg": as_float(control["ra_deg"], f"control RA/{target_id}"),
                "control_dec_deg": as_float(control["dec_deg"], f"control Dec/{target_id}"),
                "z": as_float(target.get("z"), f"redshift/{target_id}"),
                "richness": richness,
                "log10_richness": math.log10(richness),
                "target_ebv_sfd": target_ebv,
                "control_ebv_sfd": control_ebv,
                "delta_ebv_sfd": target_ebv - control_ebv,
                "contrasts": {},
            }
        if pair["contrast"] not in CONTRASTS:
            raise ValueError(f"unexpected contrast label {pair['contrast']}")
        raw_value = pair.get("contrast_njy_arcsec2", "")
        value = (None if raw_value in (None, "", "None") else
                 as_float(raw_value, f"contrast/{target_id}/{pair['contrast']}"))
        rows_by_key[key]["contrasts"][pair["contrast"]] = value

    if len(rows_by_key) != 1022:
        raise ValueError("relative-contrast rows do not map to 1,022 unique pairs")
    out = {band: [] for band in BANDS}
    for (target_id, band), row in rows_by_key.items():
        if set(row["contrasts"]) != set(CONTRASTS):
            raise ValueError(f"missing radial contrast for {target_id}/{band}")
        out[band].append(row)
    for band in BANDS:
        out[band].sort(key=lambda row: row["target_id"])
        if len(out[band]) != EXPECTED_N[band]:
            raise ValueError(f"unexpected {band} pair count")
    return out


def design_matrix(rows: list[dict]) -> tuple[np.ndarray, np.ndarray, dict]:
    y_rich = np.asarray([row["log10_richness"] for row in rows], dtype=float)
    z = np.asarray([row["z"] for row in rows], dtype=float)
    delta_ebv = np.asarray([row["delta_ebv_sfd"] for row in rows], dtype=float)
    X = np.column_stack((
        np.ones(len(rows)),
        y_rich - y_rich.mean(),
        (z - z.mean()) / 0.1,
        (delta_ebv - delta_ebv.mean()) / 0.01,
    ))
    scaling = {
        "log10_richness_center": float(y_rich.mean()),
        "redshift_center": float(z.mean()),
        "redshift_unit": 0.1,
        "delta_ebv_center": float(delta_ebv.mean()),
        "delta_ebv_unit": 0.01,
    }
    return X, y_rich, scaling


def block_bootstrap_slopes(X: np.ndarray, y: np.ndarray, labels: np.ndarray,
                           replicas: int, seed: int) -> tuple[float, float, float, int, int]:
    groups = np.unique(labels)
    group_index = {label: index for index, label in enumerate(groups)}
    row_group = np.asarray([group_index[label] for label in labels], dtype=int)
    if len(groups) < X.shape[1] + 2:
        raise ValueError("too few connected sky blocks for adjusted regression")
    full, _, rank, _ = np.linalg.lstsq(X, y, rcond=None)
    if rank != X.shape[1]:
        raise ValueError("full-sample richness model is rank deficient")

    rng = np.random.default_rng(seed)
    draws = rng.multinomial(len(groups), np.full(len(groups), 1.0 / len(groups)),
                            size=replicas)
    weights = draws[:, row_group].astype(float, copy=False)
    n_boot, n_param = replicas, X.shape[1]
    gram = np.empty((n_boot, n_param, n_param), dtype=float)
    for left in range(n_param):
        for right in range(left, n_param):
            values = weights @ (X[:, left] * X[:, right])
            gram[:, left, right] = values
            gram[:, right, left] = values
    rhs = np.column_stack([weights @ (X[:, index] * y)
                           for index in range(n_param)])
    condition = np.linalg.cond(gram)
    good = np.isfinite(condition) & (condition < 1e10)
    coeffs = np.full((n_boot, n_param), np.nan, dtype=float)
    if np.any(good):
        coeffs[good] = np.linalg.solve(gram[good], rhs[good, :, None])[..., 0]
    slopes = coeffs[:, 1]
    slopes = slopes[np.isfinite(slopes)]
    if slopes.size < replicas * 0.9:
        raise ValueError(f"only {slopes.size}/{replicas} valid block-bootstrap fits")
    low, high = np.percentile(slopes, [2.5, 97.5]).tolist()
    return float(full[1]), float(low), float(high), int(len(groups)), int(slopes.size)


def calculate(rows_by_band: dict[str, list[dict]], spatial: dict,
              replicas: int, seed_base: int) -> tuple[list[dict], list[dict]]:
    summary_rows, pair_rows = [], []
    layout_index = 0
    ra_origin = float(spatial["ra_grid_origin_deg"])
    for width in BLOCK_WIDTHS_DEG:
        for phase in GRID_PHASES:
            for band_index, band in enumerate(BANDS):
                rows = rows_by_band[band]
                control_rows = {
                    row["target_id"]: {
                        "ra_deg": row["control_ra_deg"],
                        "dec_deg": row["control_dec_deg"],
                    }
                    for row in rows
                }
                labels, blocks = connected_pair_blocks(
                    rows, control_rows, width, phase, ra_origin)
                expected_layout = next(
                    item for item in spatial["layout_summary"]
                    if item["filter"] == band
                    and float(item["block_width_deg"]) == width
                    and float(item["grid_phase_fraction"]) == phase)
                if (len(blocks) != int(expected_layout["n_connected_blocks"])
                        or max(int(block["n_pairs"]) for block in blocks)
                        != int(expected_layout["largest_block_pairs"])):
                    raise ValueError("connected blocks do not reproduce the spatial audit")
                X_all, _, scaling = design_matrix(rows)
                for contrast_index, contrast in enumerate(CONTRASTS):
                    keep = np.asarray([row["contrasts"][contrast] is not None
                                       for row in rows], dtype=bool)
                    fit_rows = [row for row, valid in zip(rows, keep) if valid]
                    X = X_all[keep]
                    y = np.asarray([row["contrasts"][contrast] for row in fit_rows], dtype=float)
                    fit_labels = labels[keep]
                    slope, low, high, n_blocks, valid = block_bootstrap_slopes(
                        X, y, fit_labels, replicas,
                        seed_base + layout_index * 1009 + band_index * 101 + contrast_index)
                    summary_rows.append({
                        "filter": band,
                        "contrast": contrast,
                        "block_width_deg": width,
                        "grid_phase_fraction": phase,
                        "n_pairs": len(fit_rows),
                        "n_connected_blocks": n_blocks,
                        "richness_slope_njy_arcsec2_per_dex": slope,
                        "block_bootstrap_ci95_low": low,
                        "block_bootstrap_ci95_high": high,
                        "bootstrap_replicas_requested": replicas,
                        "bootstrap_replicas_valid": valid,
                        "model": "equal-pair OLS: relative contrast ~ log10(CAMIRA richness) + redshift + target-control delta E(B-V)",
                        "delta_ebv_used_as_covariate_only": True,
                        "science_measurement": False,
                    })
                for row, label in zip(rows, labels):
                    pair_rows.append({
                        "target_id": row["target_id"], "control_id": row["control_id"],
                        "filter": band, "ra_deg": row["ra_deg"], "dec_deg": row["dec_deg"],
                        "z": row["z"], "richness": row["richness"],
                        "log10_richness": row["log10_richness"],
                        "target_ebv_sfd": row["target_ebv_sfd"],
                        "control_ebv_sfd": row["control_ebv_sfd"],
                        "delta_ebv_sfd": row["delta_ebv_sfd"],
                        "connected_block_id": label,
                        "block_width_deg": width,
                        "grid_phase_fraction": phase,
                        "scaling_json": json.dumps(scaling, sort_keys=True),
                        **{contrast: (row["contrasts"][contrast]
                                      if row["contrasts"][contrast] is not None else "")
                           for contrast in CONTRASTS},
                        "science_measurement": False,
                    })
                layout_index += 1
    return summary_rows, pair_rows


def sensitivity(rows: list[dict]) -> list[dict]:
    result = []
    for band in BANDS:
        for contrast in CONTRASTS:
            subset = [row for row in rows if row["filter"] == band
                      and row["contrast"] == contrast]
            lowers = [float(row["block_bootstrap_ci95_low"]) for row in subset]
            uppers = [float(row["block_bootstrap_ci95_high"]) for row in subset]
            result.append({
                "filter": band,
                "contrast": contrast,
                "min_ci_lower_across_layouts": min(lowers),
                "max_ci_lower_across_layouts": max(lowers),
                "min_ci_upper_across_layouts": min(uppers),
                "max_ci_upper_across_layouts": max(uppers),
                "all_layouts_ci_above_zero": min(lowers) > 0.0,
                "all_layouts_ci_below_zero": max(uppers) < 0.0,
            })
    return result


def render_plot(path: Path, rows: list[dict]) -> None:
    labels = ("50–300", "100–300", "200–300")
    fig, axes = plt.subplots(1, 2, figsize=(10.6, 4.6), sharey=True)
    x = np.arange(len(CONTRASTS))
    for ax, band in zip(axes, BANDS):
        selected = [next(row for row in rows if row["filter"] == band
                         and row["contrast"] == contrast
                         and row["block_width_deg"] == 5.0
                         and row["grid_phase_fraction"] == 0.0)
                    for contrast in CONTRASTS]
        estimates = np.asarray([r["richness_slope_njy_arcsec2_per_dex"]
                                for r in selected])
        lower = np.asarray([r["block_bootstrap_ci95_low"] for r in selected])
        upper = np.asarray([r["block_bootstrap_ci95_high"] for r in selected])
        ax.errorbar(x, estimates, yerr=np.vstack((estimates - lower, upper - estimates)),
                    fmt="o", color="#1f5a85", capsize=4, linewidth=1.5)
        ax.axhline(0.0, color="0.35", linewidth=0.9)
        ax.set_xticks(x, labels)
        ax.grid(axis="y", alpha=0.22)
        ax.text(0.04, 0.95, f"HSC-{band}", transform=ax.transAxes,
                ha="left", va="top", fontsize=11)
        ax.set_xlabel("Inner radius range [proper kpc], vs 300–800 kpc")
    axes[0].set_ylabel("Adjusted richness slope [nJy arcsec⁻² dex⁻¹]")
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def run(output: Path, replicas: int, seed: int) -> dict:
    output = Path(output).resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite richness-relation output: {output}")
    if replicas < 1000:
        raise ValueError("at least 1,000 block-bootstrap replicates are required")
    target, control, stack, spatial, relative, relative_table, input_hashes = verify_inputs()
    rows_by_band = make_pair_records(target, control, relative_table)
    summaries, pairs = calculate(rows_by_band, spatial, replicas, seed)
    envelopes = sensitivity(summaries)

    output.mkdir(parents=True, exist_ok=False)
    pair_path = output / "pair_richness_contrasts.tsv"
    slope_path = output / "block_bootstrap_richness_slopes.tsv"
    envelope_path = output / "layout_sensitivity.tsv"
    figure_path = output / "s23b_icl_richness_slopes.png"
    write_tsv(pair_path, pairs)
    write_tsv(slope_path, summaries)
    write_tsv(envelope_path, envelopes)
    render_plot(figure_path, summaries)
    outputs = {p.name: {"path": str(p), "sha256": sha256(p)}
               for p in (pair_path, slope_path, envelope_path, figure_path)}
    main = [row for row in summaries if row["block_width_deg"] == 5.0
            and row["grid_phase_fraction"] == 0.0]
    result = {
        "schema": "astrafek-s23b-icl-richness-relation-v1",
        "status": "complete_exploratory_richness_relation",
        "science_measurement": False,
        "image_pixels_read": False,
        "network_used": False,
        "n_pairs_by_filter": {band: len(rows) for band, rows in rows_by_band.items()},
        "estimand": "equal-pair OLS slope of within-pair inner-minus-300-800-kpc target-control surface-brightness contrast versus log10 CAMIRA richness, adjusted for redshift and target-control SFD E(B-V) mismatch",
        "richness_is_satellite_profile": False,
        "covariates": ["redshift", "target-minus-control SFD E(B-V) as a matching diagnostic; no extinction correction applied"],
        "block_widths_deg": list(BLOCK_WIDTHS_DEG),
        "grid_phases": list(GRID_PHASES),
        "replicas_per_layout_filter_contrast": replicas,
        "seed_base": seed,
        "primary_5deg_phase0_slopes": main,
        "layout_sensitivity_envelope": envelopes,
        "inputs_sha256": input_hashes,
        "relative_contrast_table_sha256": EXPECTED_RELATIVE_TABLE_SHA256,
        "connected_block_helper_sha256": EXPECTED_BLOCK_HELPER_SHA256,
        "outputs": outputs,
        "limitations": [
            "CAMIRA richness is used only as a cluster-scale proxy; no member positions, satellite luminosities, or satellite radial profiles are included.",
            "The 300-800 kpc reference annulus may contain real cluster light; these outcomes are relative radial contrasts, not absolute ICL.",
            "The residual E(B-V) mismatch is a covariate only; neither target nor control flux receives an extinction correction.",
            "Connected rectangular RA-Dec blocks approximate spatial dependence, not official HSC survey-region covariance.",
            "The regression does not propagate residual sky, PSF wings, source and sky Poisson noise, coadd covariance, centering/mask uncertainty, measurement error in richness, or selection effects.",
            "No p-values or multiplicity-adjusted claims are reported; this is exploratory and science_measurement remains false.",
        ],
        "runner": {"path": str(Path(__file__).resolve()),
                   "sha256": sha256(Path(__file__).resolve())},
    }
    (output / "summary.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False) + "\n",
        encoding="utf-8")
    assessment = [
        "# S23B ICL radial contrast versus CAMIRA richness — 2026-11-07",
        "",
        "**Status:** exploratory, `science_measurement=false`. No HSC requests or science-image pixels were used.",
        "",
        "The outcome is the saved within-pair target-minus-control surface-brightness contrast in an inner annulus minus the same pair's 300–800 proper-kpc contrast. Equal-pair OLS estimates the slope against log10 CAMIRA richness, adjusting for redshift and the target-minus-control SFD E(B−V) mismatch. The E(B−V) mismatch is only a covariate; no flux correction is applied.",
        "",
        "CAMIRA richness is a cluster-scale proxy and is not a satellite-member measurement. The member catalog remains separately approval-gated. This analysis cannot establish an ICL–satellite-profile association.",
        "",
        "## Primary 5-degree, phase-zero results",
        "",
        "| Filter | Inner annulus vs 300–800 kpc | Pairs / blocks | Slope per dex richness | 95% block-bootstrap interval |",
        "|---|---|---:|---:|---:|",
    ]
    for row in main:
        name = row["contrast"].split("_minus_")[0].replace("-", "–")
        assessment.append(
            f"| {row['filter']} | {name} | {row['n_pairs']} / {row['n_connected_blocks']} | "
            f"{row['richness_slope_njy_arcsec2_per_dex']:.4g} | "
            f"[{row['block_bootstrap_ci95_low']:.4g}, {row['block_bootstrap_ci95_high']:.4g}] |")
    assessment.extend([
        "",
        "Intervals resample connected target/control sky blocks, not individual clusters. They are exploratory and do not include a multiple-testing adjustment.",
        "",
        "## Interpretation limits",
        "",
        "- The outer 300–800 kpc annulus can contain cluster-envelope light, so the dependent variable is not an absolute ICL surface brightness.",
        "- Fluxes are not corrected for Galactic extinction. A residual E(B−V) difference is included as a nuisance covariate only.",
        "- Sky transfer, extended PSF wings, coadd covariance, mask/centering uncertainty, richness measurement uncertainty, and sample selection remain unmodeled.",
        "- A richness trend cannot substitute for satellite positions, member completeness, or a member-centered ICL/satellite cross-profile.",
        "- Keep `science_measurement=false` and do not treat an interval crossing zero as evidence of no relation.",
        "",
        "## Provenance",
        "",
        *[f"- `{path}` — SHA-256 `{digest}`" for path, digest in input_hashes.items()],
        f"- Relative-pair contrast table SHA-256: `{EXPECTED_RELATIVE_TABLE_SHA256}`",
        f"- Connected-block helper SHA-256: `{EXPECTED_BLOCK_HELPER_SHA256}`",
        f"- Runner: `{Path(__file__).resolve()}` — SHA-256 `{sha256(Path(__file__).resolve())}`",
        "",
    ])
    (output / "assessment.md").write_text("\n".join(assessment), encoding="utf-8")
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--replicas", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=20261109)
    args = parser.parse_args()
    result = run(args.output, args.replicas, args.seed)
    print(json.dumps({
        "status": result["status"],
        "science_measurement": result["science_measurement"],
        "n_pairs_by_filter": result["n_pairs_by_filter"],
        "primary_5deg_phase0_slopes": [
            {key: row[key] for key in (
                "filter", "contrast", "n_pairs", "n_connected_blocks",
                "richness_slope_njy_arcsec2_per_dex",
                "block_bootstrap_ci95_low", "block_bootstrap_ci95_high")}
            for row in result["primary_5deg_phase0_slopes"]],
        "output": str(Path(args.output).resolve()),
    }, ensure_ascii=False, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
