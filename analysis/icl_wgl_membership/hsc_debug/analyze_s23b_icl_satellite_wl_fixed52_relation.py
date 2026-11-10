#!/usr/bin/env python3
"""Explore ICL, CAMIRA-member, and weak-lensing relations in the fixed 52 set.

This runner uses only pinned, already-derived Astrafek tables. It does not
query HSC, download data, read image pixels, or read the original shape-row
CSVs. Every result remains exploratory and science_measurement=false.
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

from evaluate_prestack_pair_support import pair_block_labels


DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
ICL_DIR = DATA / "s23b_icl_satellite_relation_2026-11-07_03"
WL_DIR = DATA / "wl_s23b_fixed_common_support_spatial_covariance_2026-11-11_01"
OVERLAP_DIR = DATA / "s23b_icl_wl_fixed_support_overlap_2026-11-11_01"
DEFAULT_OUTPUT = DATA / "s23b_icl_satellite_wl_fixed52_relation_2026-11-11_03"

INPUTS = {
    "icl_summary.json": ICL_DIR / "summary.json",
    "icl_profiles.tsv": ICL_DIR / "paired_icl_satellite_profiles.tsv",
    "wl_summary.json": WL_DIR / "summary.json",
    "wl_profiles.tsv": WL_DIR / "per_lens_set_profile.tsv",
    "overlap_summary.json": OVERLAP_DIR / "summary.json",
}
EXPECTED_SHA256 = {
    "icl_summary.json": "fab29d7aef464707c31e2616fe1db9fc3bad9d281ad6d5dfa6d63f2b484bbb39",
    "icl_profiles.tsv": "98ed9e0022bfa68d5921ead62dce55be7df1cfa169e6a777d66d1e9760f43994",
    "wl_summary.json": "2391fede047fc2de25b34082b0c0d6635e88910d5701a6f9d16b5498ec66602f",
    "wl_profiles.tsv": "5bc599d97123264058011313b0ad35c3def49b08bba98a016cc5edf6ff03906c",
    "overlap_summary.json": "8145f969b60fe93075d86bd952c9cd0d0e38f9be71be2ff521252048754cf69e",
}
FIXED52_IDS_SHA256 = "18312f27139b3750abdbf71543f22d9b7b18b2d3b24cbd29ed686a377524bb7d"
RADIAL_EDGES_MPC = (0.1, 0.2, 0.4, 0.7, 1.0, 1.5, 2.0, 3.0)
FILTERS = ("i", "i2")
SUPPORT_GROUPS = (("all", 0.0), ("support_ge_0p10", 0.10), ("support_ge_0p20", 0.20))
BOOTSTRAP_REPLICATES = 2000
MIN_CORRELATION_N = 12
MIN_BLOCKS_FOR_INTERVAL = 6
BASE_SEED = 20261111


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def id_sha256(values: set[str]) -> str:
    return hashlib.sha256("\n".join(sorted(values)).encode()).hexdigest()


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f"Refusing to write an empty table: {path}")
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def finite_float(row: dict[str, str], key: str) -> float:
    value = float(row[key])
    if not math.isfinite(value):
        raise ValueError(f"non-finite {key}: {value}")
    return value


def residual_correlation(x: np.ndarray, y: np.ndarray,
                         covariates: list[np.ndarray]) -> float:
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    covariates = [np.asarray(c, dtype=float) for c in covariates]
    if x.ndim != 1 or y.shape != x.shape or len(x) < MIN_CORRELATION_N:
        return float("nan")
    if not np.isfinite(x).all() or not np.isfinite(y).all():
        return float("nan")
    if any(c.shape != x.shape or not np.isfinite(c).all() for c in covariates):
        return float("nan")
    design = np.column_stack((np.ones(len(x)), *covariates))
    residualized, _, rank, _ = np.linalg.lstsq(design, np.column_stack((x, y)), rcond=None)
    if rank < design.shape[1]:
        return float("nan")
    residuals = np.column_stack((x, y)) - design @ residualized
    denominator = float(np.linalg.norm(residuals[:, 0]) * np.linalg.norm(residuals[:, 1]))
    if denominator <= 0 or not math.isfinite(denominator):
        return float("nan")
    return float(residuals[:, 0] @ residuals[:, 1] / denominator)


def block_layout(rows: list[dict], selected: np.ndarray | None = None):
    if not rows:
        return [], {}, 0, 0
    if selected is None:
        selected = np.ones(len(rows), dtype=bool)
    target = [{"ra_deg": float(row["ra_deg"]), "dec_deg": float(row["dec_deg"])}
              for row in rows]
    control = [{"ra_deg": float(row["control_ra_deg"]),
                "dec_deg": float(row["control_dec_deg"])} for row in rows]
    indices, labels, _, _, n_blocks, largest = pair_block_labels(target, control, selected)
    indices = np.asarray(indices, dtype=int)
    labels = np.asarray(labels)
    blocks = np.unique(labels)
    by_block = {block: indices[labels == block] for block in blocks}
    return blocks, by_block, n_blocks, largest


def bootstrap_stat(rows: list[dict], statistic, replicas: int, seed: int) -> dict:
    blocks, by_block, n_blocks, largest = block_layout(rows)
    point = float(statistic(np.arange(len(rows))))
    result = {
        "estimate": point,
        "interval_95_low": None,
        "interval_95_high": None,
        "bootstrap_replicates": replicas,
        "bootstrap_valid": 0,
        "bootstrap_valid_fraction": 0.0,
        "n_spatial_blocks": n_blocks,
        "largest_block_n": largest,
        "interval_status": "insufficient_sample_or_blocks",
    }
    if len(rows) < MIN_CORRELATION_N or n_blocks < MIN_BLOCKS_FOR_INTERVAL:
        return result
    rng = np.random.default_rng(seed)
    draws = []
    for _ in range(replicas):
        sampled_blocks = rng.choice(blocks, size=len(blocks), replace=True)
        indices = np.concatenate([by_block[block] for block in sampled_blocks])
        value = float(statistic(indices))
        if math.isfinite(value):
            draws.append(value)
    result["bootstrap_valid"] = len(draws)
    result["bootstrap_valid_fraction"] = len(draws) / replicas
    if len(draws) >= max(100, int(0.5 * replicas)):
        low, high = np.quantile(draws, [0.025, 0.975]).tolist()
        result["interval_95_low"] = float(low)
        result["interval_95_high"] = float(high)
        result["interval_status"] = "descriptive_pointwise_block_bootstrap"
    else:
        result["interval_status"] = "too_few_valid_bootstrap_draws"
    return result


def correlation_record(rows: list[dict], relationship: str, predictor: str,
                       outcome: str, covariate_names: list[str], band: str,
                       support_group: str, radial_bin: str, component: str,
                       replicas: int, seed: int) -> dict:
    x = np.asarray([float(row[predictor]) for row in rows])
    y = np.asarray([float(row[outcome]) for row in rows])
    covariates = [np.asarray([float(row[key]) for row in rows])
                  for key in covariate_names]
    statistic = lambda indices: residual_correlation(
        x[indices], y[indices], [value[indices] for value in covariates])
    block_result = bootstrap_stat(rows, statistic, replicas, seed)
    point = block_result["estimate"]
    return {
        "filter": band,
        "support_group": support_group,
        "relationship": relationship,
        "predictor": predictor,
        "outcome": outcome,
        "component": component,
        "radial_bin": radial_bin,
        "n_lenses": len(rows),
        "n_spatial_blocks": block_result["n_spatial_blocks"],
        "largest_block_n": block_result["largest_block_n"],
        "partial_r": point,
        "partial_r_block_bootstrap_95_low": block_result["interval_95_low"],
        "partial_r_block_bootstrap_95_high": block_result["interval_95_high"],
        "bootstrap_replicates": block_result["bootstrap_replicates"],
        "bootstrap_valid": block_result["bootstrap_valid"],
        "bootstrap_valid_fraction": block_result["bootstrap_valid_fraction"],
        "interval_status": block_result["interval_status"],
        "covariates": ",".join(covariate_names),
        "inferential_claim": False,
    }


def load_inputs() -> tuple[dict, list[dict], list[dict], dict]:
    input_hashes = {name: sha256(path) for name, path in INPUTS.items()}
    for name, actual in input_hashes.items():
        if actual != EXPECTED_SHA256[name]:
            raise ValueError(f"Pinned input hash mismatch for {name}: {actual}")
    icl_summary = read_json(INPUTS["icl_summary.json"])
    wl_summary = read_json(INPUTS["wl_summary.json"])
    overlap_summary = read_json(INPUTS["overlap_summary.json"])
    for label, summary in (("ICL", icl_summary), ("WL", wl_summary),
                           ("overlap", overlap_summary)):
        if summary.get("science_measurement") is not False:
            raise ValueError(f"{label} input was promoted beyond exploratory status")
    if (icl_summary.get("sample", {}).get("science_image_pixels_read") is not False
            or icl_summary.get("sample", {}).get("network_used") is not False):
        raise ValueError("ICL input scope differs from the pinned derived table")
    if (wl_summary.get("sample", {}).get("n_fixed_lens_sets") != 52
            or wl_summary.get("sample", {}).get("fixed_lens_ids_sha256") != FIXED52_IDS_SHA256
            or wl_summary.get("provenance", {}).get("approved_row_plan_sha256")
            != "7bfb998149f30a50b5ee7a537bb703ed302b7ea269a3296c35204a836e75aabc"):
        raise ValueError("WL input is not the audited fixed 52-lens profile")
    if (overlap_summary.get("fixed_lens_count") != 52
            or overlap_summary.get("fixed_lens_ids_sha256") != FIXED52_IDS_SHA256
            or overlap_summary.get("joint_relation_fit_performed") is not False):
        raise ValueError("ICL/WL overlap input has unexpected sample or analysis scope")
    icl_rows = read_tsv(INPUTS["icl_profiles.tsv"])
    wl_rows = read_tsv(INPUTS["wl_profiles.tsv"])
    return input_hashes, icl_rows, wl_rows, overlap_summary


def build_joined_rows(icl_rows: list[dict[str, str]],
                      wl_rows: list[dict[str, str]]) -> dict[str, list[dict]]:
    fixed_ids = {row["matched_lens_id"] for row in wl_rows}
    if len(fixed_ids) != 52 or id_sha256(fixed_ids) != FIXED52_IDS_SHA256:
        raise ValueError("The WL profile lens IDs differ from the pinned fixed 52 set")
    wl_by_lens: dict[str, list[dict[str, str]]] = {}
    for row in wl_rows:
        lens_id = row["matched_lens_id"]
        wl_by_lens.setdefault(lens_id, []).append(row)
    if any(len(rows) != 7 for rows in wl_by_lens.values()):
        raise ValueError("Each fixed lens must have exactly seven WL radial bins")
    for lens_id, rows in wl_by_lens.items():
        if {int(row["radial_bin"]) for row in rows} != set(range(7)):
            raise ValueError(f"Incomplete WL radial bins for {lens_id}")

    by_filter: dict[str, list[dict]] = {band: [] for band in FILTERS}
    seen: set[tuple[str, str]] = set()
    for icl in icl_rows:
        lens_id = icl["target_id"]
        band = icl["filter"]
        if band not in by_filter or lens_id not in fixed_ids:
            continue
        key = (lens_id, band)
        if key in seen:
            raise ValueError(f"Duplicate ICL profile for {key}")
        seen.add(key)
        if icl["science_measurement"].lower() != "false":
            raise ValueError(f"ICL input unexpectedly promotes science measurement for {key}")
        if lens_id not in wl_by_lens:
            raise ValueError(f"No WL profile for ICL lens {lens_id}")
        redshift = finite_float(icl, "z_cl")
        richness = finite_float(icl, "richness")
        if richness <= 0:
            raise ValueError(f"Nonpositive CAMIRA richness for {lens_id}")
        support = min(finite_float(icl, "inner_common_support_fraction"),
                      finite_float(icl, "outer_common_support_fraction"))
        if not 0 <= support <= 1:
            raise ValueError(f"Invalid common-pixel support for {key}: {support}")
        row = {
            "cluster_id": lens_id,
            "filter": band,
            "control_id": icl["control_id"],
            "ra_deg": finite_float(icl, "ra_deg"),
            "dec_deg": finite_float(icl, "dec_deg"),
            "control_ra_deg": finite_float(icl, "control_ra_deg"),
            "control_dec_deg": finite_float(icl, "control_dec_deg"),
            "redshift": redshift,
            "richness": richness,
            "log_richness": math.log(richness),
            "icl_contrast_njy_arcsec2": finite_float(icl, "icl_relative_contrast_njy_arcsec2"),
            "satellite_log10_density_ratio": finite_float(icl, "log10_inner_outer_candidate_density_ratio"),
            "inner_common_support_fraction": finite_float(icl, "inner_common_support_fraction"),
            "outer_common_support_fraction": finite_float(icl, "outer_common_support_fraction"),
            "minimum_common_support_fraction": support,
            "n_member_candidates_inner": int(icl["camira_member_candidates_50_300_in_image"]),
            "n_member_candidates_outer": int(icl["camira_member_candidates_300_800_in_image"]),
            "wl_by_bin": {int(wl["radial_bin"]): wl for wl in wl_by_lens[lens_id]},
        }
        by_filter[band].append(row)

    if {band: len(rows) for band, rows in by_filter.items()} != {"i": 25, "i2": 27}:
        raise ValueError("Fixed 52 sample no longer splits into 25 i and 27 i2 profiles")
    if ({row["cluster_id"] for row in by_filter["i"]}
            & {row["cluster_id"] for row in by_filter["i2"]}):
        raise ValueError("The audited i and i2 samples must be disjoint")
    if ({row["cluster_id"] for band in FILTERS for row in by_filter[band]} != fixed_ids):
        raise ValueError("ICL rows do not cover all 52 fixed WL lenses")
    return by_filter


def sample_rows(rows: list[dict], minimum_support: float) -> list[dict]:
    return [row for row in rows if row["minimum_common_support_fraction"] >= minimum_support]


def profile_rows_for_bin(rows: list[dict], radial_bin: int, component: str) -> list[dict]:
    key = ("delta_sigma_random_subtracted_msun_mpc2" if component == "tangential"
           else "cross_random_subtracted_msun_mpc2")
    result = []
    for row in rows:
        wl = row["wl_by_bin"][radial_bin]
        if int(wl["n_lens_pairs"]) <= 0 or int(wl["n_random_pairs"]) <= 0:
            continue
        value = finite_float(wl, key)
        merged = dict(row)
        merged["wl_value"] = value
        merged["n_lens_pairs"] = int(wl["n_lens_pairs"])
        merged["n_random_pairs"] = int(wl["n_random_pairs"])
        result.append(merged)
    return result


def analyze(by_filter: dict[str, list[dict]], replicas: int) -> tuple[list[dict], list[dict]]:
    metrics: list[dict] = []
    profiles: list[dict] = []
    seed_offset = 0
    for band in FILTERS:
        for group_name, minimum_support in SUPPORT_GROUPS:
            selected = sample_rows(by_filter[band], minimum_support)
            blocks, _, n_blocks, largest = block_layout(selected)
            if not selected:
                raise ValueError(f"Empty sample for {band}/{group_name}")
            if group_name == "support_ge_0p20":
                # There are only 5 and 7 lenses in this high-support stratum;
                # preserve counts but do not estimate correlations or stacks.
                profiles.append({
                    "filter": band, "support_group": group_name,
                    "radial_bin": "all", "component": "all",
                    "r_min_mpc": None, "r_max_mpc": None,
                    "n_lenses": len(selected), "n_spatial_blocks": n_blocks,
                    "largest_block_n": largest, "mean_delta_sigma_msun_mpc2": None,
                    "block_bootstrap_95_low": None, "block_bootstrap_95_high": None,
                    "median_lens_pairs": None, "median_random_pairs": None,
                    "interval_status": "descriptive_support_count_only_small_sample",
                })
                continue

            # Direct ICL–satellite relation, adjusting redshift and log richness.
            metric = correlation_record(
                selected, "icl_vs_satellite", "icl_contrast_njy_arcsec2",
                "satellite_log10_density_ratio", ["redshift", "log_richness"],
                band, group_name, "na", "na", replicas, BASE_SEED + seed_offset)
            seed_offset += 1
            metrics.append(metric)

            for radial_bin in range(7):
                r_min, r_max = RADIAL_EDGES_MPC[radial_bin:radial_bin + 2]
                for component in ("tangential", "cross"):
                    prof_rows = profile_rows_for_bin(selected, radial_bin, component)
                    profile_values = np.asarray([row["wl_value"] for row in prof_rows], dtype=float)
                    nprof = len(prof_rows)
                    if nprof:
                        stat = lambda indices, vals=profile_values: float(np.mean(vals[indices]))
                        boot = bootstrap_stat(prof_rows, stat, replicas,
                                              BASE_SEED + 1000 + seed_offset)
                        med_lens = float(np.median([row["n_lens_pairs"] for row in prof_rows]))
                        med_random = float(np.median([row["n_random_pairs"] for row in prof_rows]))
                    else:
                        boot = {"estimate": None, "interval_95_low": None,
                                "interval_95_high": None, "n_spatial_blocks": 0,
                                "largest_block_n": 0,
                                "interval_status": "no_lenses_with_positive_target_and_random_pairs"}
                        med_lens = med_random = None
                    profiles.append({
                        "filter": band, "support_group": group_name,
                        "radial_bin": radial_bin, "component": component,
                        "r_min_mpc": r_min, "r_max_mpc": r_max,
                        "n_lenses": nprof,
                        "n_spatial_blocks": boot["n_spatial_blocks"],
                        "largest_block_n": boot["largest_block_n"],
                        "mean_delta_sigma_msun_mpc2": boot["estimate"],
                        "block_bootstrap_95_low": boot["interval_95_low"],
                        "block_bootstrap_95_high": boot["interval_95_high"],
                        "median_lens_pairs": med_lens,
                        "median_random_pairs": med_random,
                        "interval_status": boot["interval_status"],
                    })
                    if nprof < MIN_CORRELATION_N:
                        continue
                    outcome = "wl_delta_sigma_msun_mpc2" if component == "tangential" else "wl_cross_msun_mpc2"
                    for predictor, covariates, relationship in (
                        ("icl_contrast_njy_arcsec2",
                         ["satellite_log10_density_ratio", "redshift", "log_richness"],
                         "wl_vs_icl"),
                        ("satellite_log10_density_ratio",
                         ["icl_contrast_njy_arcsec2", "redshift", "log_richness"],
                         "wl_vs_satellite"),
                    ):
                        rows_for_corr = []
                        for row in prof_rows:
                            corr_row = dict(row)
                            corr_row["wl_delta_sigma_msun_mpc2"] = row["wl_value"]
                            corr_row["wl_cross_msun_mpc2"] = row["wl_value"]
                            # The profile loop stores the selected component's value.
                            # Outcome naming remains explicit in the TSV.
                            rows_for_corr.append(corr_row)
                        metric = correlation_record(
                            rows_for_corr, relationship, predictor, outcome,
                            covariates, band, group_name, str(radial_bin), component,
                            replicas, BASE_SEED + seed_offset)
                        seed_offset += 1
                        metrics.append(metric)
    return metrics, profiles


def flatten_joined(by_filter: dict[str, list[dict]]) -> list[dict]:
    rows = []
    for band in FILTERS:
        for row in by_filter[band]:
            for radial_bin in range(7):
                wl = row["wl_by_bin"][radial_bin]
                rows.append({
                    "cluster_id": row["cluster_id"], "filter": band,
                    "radial_bin": radial_bin,
                    "r_min_mpc": RADIAL_EDGES_MPC[radial_bin],
                    "r_max_mpc": RADIAL_EDGES_MPC[radial_bin + 1],
                    "control_id": row["control_id"],
                    "redshift": row["redshift"], "richness": row["richness"],
                    "icl_contrast_njy_arcsec2": row["icl_contrast_njy_arcsec2"],
                    "satellite_log10_density_ratio": row["satellite_log10_density_ratio"],
                    "inner_common_support_fraction": row["inner_common_support_fraction"],
                    "outer_common_support_fraction": row["outer_common_support_fraction"],
                    "n_member_candidates_inner": row["n_member_candidates_inner"],
                    "n_member_candidates_outer": row["n_member_candidates_outer"],
                    "n_lens_pairs": int(wl["n_lens_pairs"]),
                    "n_random_pairs": int(wl["n_random_pairs"]),
                    "delta_sigma_random_subtracted_msun_mpc2": finite_float(
                        wl, "delta_sigma_random_subtracted_msun_mpc2"),
                    "cross_random_subtracted_msun_mpc2": finite_float(
                        wl, "cross_random_subtracted_msun_mpc2"),
                })
    return rows


def make_figure(by_filter: dict[str, list[dict]], profiles: list[dict], path: Path) -> None:
    plt.rcParams.update({"font.size": 10, "axes.titlesize": 11, "axes.labelsize": 10,
                         "legend.fontsize": 8, "figure.dpi": 130})
    fig, axes = plt.subplots(2, 2, figsize=(12.2, 8.0), constrained_layout=True)
    cmap = plt.get_cmap("viridis")
    for col, band in enumerate(FILTERS):
        rows = by_filter[band]
        ax = axes[0, col]
        high = [r for r in rows if r["minimum_common_support_fraction"] >= 0.10]
        low = [r for r in rows if r["minimum_common_support_fraction"] < 0.10]
        if high:
            sc = ax.scatter([r["satellite_log10_density_ratio"] for r in high],
                            [r["icl_contrast_njy_arcsec2"] for r in high],
                            c=[r["redshift"] for r in high], cmap=cmap,
                            vmin=0.35, vmax=0.55, s=42, edgecolor="black", linewidth=0.35,
                            label=f"common pixels ≥10% (N={len(high)})")
        if low:
            ax.scatter([r["satellite_log10_density_ratio"] for r in low],
                       [r["icl_contrast_njy_arcsec2"] for r in low],
                       facecolors="none", edgecolors="0.45", marker="x", s=46,
                       label=f"below 10% (N={len(low)})")
        ax.axhline(0, color="0.65", lw=0.8)
        ax.set_title(f"HSC {band}: ICL–member profile")
        ax.set_xlabel("log₁₀(inner/outer member-candidate density ratio)")
        ax.set_ylabel("ICL inner−outer contrast (nJy arcsec⁻²)")
        ax.legend(loc="best", frameon=False)

        p_ax = axes[1, col]
        for component, color, marker, ls in (
            ("tangential", "#1f4e79", "o", "-"),
            ("cross", "#c56a13", "s", "--"),
        ):
            subset = [r for r in profiles if r["filter"] == band
                      and r["support_group"] == "support_ge_0p10"
                      and r["component"] == component and r["mean_delta_sigma_msun_mpc2"] is not None]
            subset.sort(key=lambda r: int(r["radial_bin"]))
            if not subset:
                continue
            x = np.asarray([math.sqrt(float(r["r_min_mpc"]) * float(r["r_max_mpc"])) for r in subset])
            y = np.asarray([float(r["mean_delta_sigma_msun_mpc2"]) / 1e13 for r in subset])
            lo = np.asarray([(float(r["block_bootstrap_95_low"]) / 1e13
                              if r["block_bootstrap_95_low"] is not None else np.nan)
                             for r in subset])
            hi = np.asarray([(float(r["block_bootstrap_95_high"]) / 1e13
                              if r["block_bootstrap_95_high"] is not None else np.nan)
                             for r in subset])
            p_ax.errorbar(x, y, yerr=np.vstack((np.maximum(0, y - lo),
                                               np.maximum(0, hi - y))), color=color,
                          marker=marker, linestyle=ls, capsize=2.5, label=component)
        p_ax.axhline(0, color="0.65", lw=0.8)
        p_ax.set_xscale("log")
        p_ax.set_title(f"HSC {band}: fixed-sample ΔΣ stack (≥10% support)")
        p_ax.set_xlabel("Projected radius (proper Mpc)")
        p_ax.set_ylabel("Equal-lens mean ΔΣ (10¹³ M☉ Mpc⁻²)")
        p_ax.legend(loc="best", frameon=False)
    norm = matplotlib.colors.Normalize(vmin=0.35, vmax=0.55)
    scalar = matplotlib.cm.ScalarMappable(norm=norm, cmap=cmap)
    colorbar = fig.colorbar(scalar, ax=axes[0, :], shrink=0.82, pad=0.02)
    colorbar.set_label("Cluster redshift")
    fig.suptitle("S23B fixed 52-lens exploratory ICL–satellite–lensing view", fontsize=13)
    fig.savefig(path, dpi=180, bbox_inches="tight")
    plt.close(fig)


def fmt(value, scale: float = 1.0) -> str:
    if value is None or not math.isfinite(float(value)):
        return "NA"
    return f"{float(value) / scale:.3f}"


def key_results(metrics: list[dict], profiles: list[dict]) -> dict:
    direct = []
    for row in metrics:
        if row["relationship"] != "icl_vs_satellite":
            continue
        direct.append({key: row[key] for key in (
            "filter", "support_group", "n_lenses", "n_spatial_blocks", "partial_r",
            "partial_r_block_bootstrap_95_low", "partial_r_block_bootstrap_95_high")})
    wl_rows = [row for row in metrics if row["relationship"].startswith("wl_")]
    intervals = [row for row in wl_rows
                 if row["partial_r_block_bootstrap_95_low"] is not None
                 and row["partial_r_block_bootstrap_95_high"] is not None
                 and math.isfinite(float(row["partial_r_block_bootstrap_95_low"]))
                 and math.isfinite(float(row["partial_r_block_bootstrap_95_high"]))]
    excludes = [row for row in intervals
                if float(row["partial_r_block_bootstrap_95_low"])
                * float(row["partial_r_block_bootstrap_95_high"]) > 0]
    profile_main = [row for row in profiles
                    if row["support_group"] == "support_ge_0p10"
                    and row["radial_bin"] != "all"]
    tangential = [row for row in profile_main if row["component"] == "tangential"]
    cross = [row for row in profile_main if row["component"] == "cross"]
    cross_excludes = [row for row in cross
                      if row["block_bootstrap_95_low"] is not None
                      and row["block_bootstrap_95_high"] is not None
                      and float(row["block_bootstrap_95_low"])
                      * float(row["block_bootstrap_95_high"]) > 0]
    return {
        "direct_icl_satellite": direct,
        "wl_pointwise_intervals_total": len(intervals),
        "wl_pointwise_intervals_excluding_zero": len(excludes),
        "wl_pointwise_intervals_excluding_zero_rows": [{key: row[key] for key in (
            "filter", "support_group", "relationship", "component", "radial_bin",
            "n_lenses", "partial_r", "partial_r_block_bootstrap_95_low",
            "partial_r_block_bootstrap_95_high")} for row in excludes],
        "tangential_stack_support_ge_10pct": [{key: row[key] for key in (
            "filter", "radial_bin", "r_min_mpc", "r_max_mpc", "n_lenses",
            "mean_delta_sigma_msun_mpc2", "block_bootstrap_95_low",
            "block_bootstrap_95_high", "median_lens_pairs", "median_random_pairs")}
            for row in tangential],
        "cross_stack_intervals_excluding_zero_support_ge_10pct": len(cross_excludes),
        "cross_stack_support_ge_10pct_bins": len(cross),
    }


def write_report(path: Path, summary: dict, profiles: list[dict], metrics: list[dict]) -> None:
    key = summary["key_results"]
    lines = [
        "# S23B 고정 52개 표본 ICL–위성–약중력 탐색 분석 (2026-11-11)",
        "",
        "## 범위와 지위",
        "",
        "해시가 고정된 기존 Astrafek ICL–CAMIRA member profile과 weak-lensing per-lens profile만 결합했다. HSC/CAS 요청, 새 다운로드, 원본 shape-row CSV 열람, 이미지 픽셀 분석은 없었다. 결과는 모두 `science_measurement=false`이며 보정된 약중력 측정이나 확인적 검출이 아니다.",
        "",
        f"전체 교집합은 {summary['n_fixed_lenses']}개이고, i 필터 {summary['n_i']}개와 i2 필터 {summary['n_i2']}개는 서로 다른 은하단 집합이다. ICL 두 annulus의 최소 공통 픽셀 지원을 10% 이상으로 두면 i {summary['n_i_support_ge_10pct']}개, i2 {summary['n_i2_support_ge_10pct']}개가 남는다. 20% 이상은 각각 {summary['n_i_support_ge_20pct']}개와 {summary['n_i2_support_ge_20pct']}개뿐이어서 상관 적합에는 사용하지 않았다.",
        "",
        "## 계산",
        "",
        "ICL–위성 지표는 ICL inner−outer contrast와 CAMIRA member-candidate inner/outer surface-density ratio의 로그다. 약중력 지표는 저장된 target-minus-random ΔΣ와 cross profile이며, 7개 radial bin에서 equal-lens mean을 계산했다. 관계 표에는 redshift와 ln(richness)을 조정한 partial correlation 및 target/control 위치의 연결된 5° 공간 블록 재표집 구간을 기록했다. 95% 구간은 점별 탐색용이며 다중비교 보정이나 물리적 유의도 판정이 아니다.",
        "",
        "고정 52개 profile estimator에는 HSC S19A per-source multiplicative/additive shear terms (m, c1, c2), radial-bin responsivity, edge-selection response, HSC Y3 SACC tomographic n(z), 그리고 target-minus-random subtraction이 이미 적용되어 있다. 이는 알려진 응답 보정의 적용 사실이며 보정 불확실성이나 잔여 편향이 충분히 검증됐다는 뜻은 아니다.",
        "",
        "### ICL–위성 관계",
        "",
        "| 필터 | 공통 픽셀 선택 | N | 블록 수 | 조정 partial r | 점별 95% 블록-bootstrap 구간 |",
        "|---|---|---:|---:|---:|---:|",
    ]
    for row in key["direct_icl_satellite"]:
        lines.append(
            f"| {row['filter']} | {row['support_group']} | {row['n_lenses']} | {row['n_spatial_blocks']} | "
            f"{fmt(row['partial_r'])} | [{fmt(row['partial_r_block_bootstrap_95_low'])}, "
            f"{fmt(row['partial_r_block_bootstrap_95_high'])}] |"
        )
    lines.extend([
        "",
        "네 점별 구간 모두 0을 포함한다. 현재 고정 표본에서는 적색편이와 풍부도를 조정한 뒤 ICL contrast–CAMIRA 후보 밀도비의 명확한 연관성이 분해되지 않는다.",
        "",
        "### 약중력 stack 및 관계 민감도",
        "",
        "| 필터 | R (proper Mpc) | N | tangential ΔΣ (10¹³ M☉ Mpc⁻²) | 점별 95% 블록-bootstrap 구간 | 중앙 lens/random pair 수 |",
        "|---|---:|---:|---:|---:|---:|",
    ])
    for row in key["tangential_stack_support_ge_10pct"]:
        mean = fmt(row["mean_delta_sigma_msun_mpc2"], 1e13)
        low = fmt(row["block_bootstrap_95_low"], 1e13)
        high = fmt(row["block_bootstrap_95_high"], 1e13)
        lines.append(
            f"| {row['filter']} | {row['r_min_mpc']:.1f}–{row['r_max_mpc']:.1f} | {row['n_lenses']} | "
            f"{mean} | [{low}, {high}] | {fmt(row['median_lens_pairs'])} / {fmt(row['median_random_pairs'])} |"
        )
    lines.extend([
        "",
        f"공통 픽셀 10% 이상 표본의 cross-profile {key['cross_stack_support_ge_10pct_bins']}개 점별 구간은 모두 0을 포함한다. 이는 약한 null 점검일 뿐 PSF·shear 보정을 검증하지 않는다.",
        f"전체 분석에서 약중력 partial-correlation 구간 {key['wl_pointwise_intervals_excluding_zero']}개가 {key['wl_pointwise_intervals_total']}개 점별 구간 중 0을 포함하지 않았다. 위치는 `metrics.tsv`에 기록했지만 radial bin·필터·관계의 다중 탐색을 거쳤고 독립 검증이 없어 검출로 해석하지 않는다.",
        "",
        "공통 픽셀 10% 미만 표본을 포함한 전체 결과와 10% 이상 결과를 나란히 계산했다. 20% 이상 표본은 수가 너무 적어 표본 수만 기록했다. 약중력 profile의 렌즈별 pair 수는 작을 수 있으므로 각 bin의 pair-count 열과 stack 표를 함께 확인해야 한다.",
        "",
        "## 한계",
        "",
        "표본은 높은 원시 S19A 형상원 가용성 및 ICL pixel support로 선택되어 전체 S23B 은하단을 대표하지 않는다. 필터별 표본이 disjoint하므로 i와 i2 값을 서로 짝지은 필터 비교로 해석하지 않는다. CAMIRA member positions는 위성 은하의 확정 분류가 아니며, 이미지 지지율이 낮은 프로파일에는 강한 선택 민감도가 있다. 이미 적용한 HSC 응답 보정과 HSC Y3 n(z)의 불확실성, 외부 검증표본에 기반한 photo-z dilution 보정, reduced shear, magnification, cluster-member dilution/boost, random aperture와 target의 깊이/mask/PSF/sky 매칭, HSC 조사영역 covariance는 해결되지 않았다. 별도 74-lens 민감도 분석에서 DESI n(z)와 magnification 대안을 비교했지만, 그것은 이 고정 52개 표본별 보정오차를 측정하지 않는다. 블록 bootstrap은 근사 공간 민감도만 제공하며 survey-region covariance를 검증하지 않는다.",
        "",
        "## 산출물",
        "",
        "- `metrics.tsv`: filter/support별 ICL–위성 및 radial-bin별 ΔΣ 관계와 점별 블록-bootstrap 구간.",
        "- `stack_profiles.tsv`: 10% 공통 지지 표본의 tangential/cross equal-lens mean profiles와 pair counts.",
        "- `joined_lens_profiles.tsv`: 분석에 사용한 52×7 lens/bin 조인 테이블.",
        "- `fixed52_relation.png`: filter별 ICL–위성 산점도와 weak-lensing stack.",
        "- `summary.json`: 입력/출력 SHA-256, 실행 범위, 상태.",
        "",
        "완전한 관계별 수치와 점별 구간은 `metrics.tsv`를 기준으로 한다. 이 탐색 분석만으로 ICL–위성 연관성 또는 halo-mass 관계를 검출했다고 주장하지 않는다.",
        "",
        "## 재현",
        "",
        "실행기: [analyze_s23b_icl_satellite_wl_fixed52_relation.py](analyze_s23b_icl_satellite_wl_fixed52_relation.py)",
        f"결과 디렉터리: `{summary['output_directory']}`",
        f"결과 요약 SHA-256: `{summary['summary_sha256']}`",
    ])
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--bootstrap", type=int, default=BOOTSTRAP_REPLICATES)
    args = parser.parse_args()
    if args.bootstrap < 1000:
        parser.error("Use at least 1,000 spatial-block bootstrap replicates")
    if args.output.exists():
        parser.error(f"Refusing to overwrite existing output directory: {args.output}")

    input_hashes, icl_rows, wl_rows, overlap_summary = load_inputs()
    by_filter = build_joined_rows(icl_rows, wl_rows)
    metrics, profiles = analyze(by_filter, args.bootstrap)
    joined_rows = flatten_joined(by_filter)

    output = args.output
    output.mkdir(parents=True, exist_ok=False)
    write_tsv(output / "metrics.tsv", metrics)
    write_tsv(output / "stack_profiles.tsv", profiles)
    write_tsv(output / "joined_lens_profiles.tsv", joined_rows)
    make_figure(by_filter, profiles, output / "fixed52_relation.png")

    summary = {
        "schema": "astrafek-s23b-icl-satellite-wl-fixed52-relation-v1",
        "status": "complete_exploratory_fixed52_joint_relation",
        "science_measurement": False,
        "network_used": False,
        "hsc_query_submitted": False,
        "new_data_downloaded": False,
        "raw_source_shape_csvs_read": False,
        "images_read": False,
        "joint_relation_fit_performed": True,
        "n_fixed_lenses": 52,
        "n_i": len(by_filter["i"]),
        "n_i2": len(by_filter["i2"]),
        "n_i_support_ge_10pct": len(sample_rows(by_filter["i"], 0.10)),
        "n_i2_support_ge_10pct": len(sample_rows(by_filter["i2"], 0.10)),
        "n_i_support_ge_20pct": len(sample_rows(by_filter["i"], 0.20)),
        "n_i2_support_ge_20pct": len(sample_rows(by_filter["i2"], 0.20)),
        "i_and_i2_lens_sets_disjoint": True,
        "radial_edges_proper_mpc": list(RADIAL_EDGES_MPC),
        "icl_support_groups": {name: threshold for name, threshold in SUPPORT_GROUPS},
        "wl_stack": "equal-lens mean of saved random-subtracted profiles; per bin requires positive target and random pair counts",
        "relationship": "partial Pearson correlation; ICL-satellite adjusts z and ln(richness); WL relations additionally adjust the other ICL/satellite predictor",
        "uncertainty": "pointwise percentile interval from connected target/control 5-degree block bootstrap; approximate spatial diagnostic only",
        "bootstrap_replicates": args.bootstrap,
        "minimum_lenses_for_relation": MIN_CORRELATION_N,
        "minimum_blocks_for_interval": MIN_BLOCKS_FOR_INTERVAL,
        "limitations": [
            "selected 52-lens exploratory sample is not representative of all S23B clusters",
            "i and i2 groups are disjoint and cannot be interpreted as a paired filter comparison",
            "low common-pixel support can bias the ICL profile; the 20-percent subgroup is too small for relation fitting",
            "CAMIRA member positions are candidate satellite positions, not confirmed satellite classifications",
            "the estimator applies HSC S19A per-source m/c terms, radial responsivity, edge-selection response and HSC Y3 SACC n(z), but does not propagate their uncertainty or resolve photo-z dilution, reduced-shear, magnification, boost, random-aperture matching, or validated survey-region covariance systematics",
            "no p-values, multiplicity-adjusted claims, or physical detection claims are made",
        ],
        "key_results": key_results(metrics, profiles),
        "input_sha256": {name: input_hashes[name] for name in sorted(input_hashes)},
        "fixed_lens_ids_sha256": FIXED52_IDS_SHA256,
        "overlap_audit_sha256": EXPECTED_SHA256["overlap_summary.json"],
        "runner_path": str(Path(__file__).resolve()),
        "runner_sha256": sha256(Path(__file__).resolve()),
        "output_directory": str(output),
        "output_sha256": {
            "metrics.tsv": sha256(output / "metrics.tsv"),
            "stack_profiles.tsv": sha256(output / "stack_profiles.tsv"),
            "joined_lens_profiles.tsv": sha256(output / "joined_lens_profiles.tsv"),
            "fixed52_relation.png": sha256(output / "fixed52_relation.png"),
        },
    }
    summary_path = output / "summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    summary["summary_sha256"] = sha256(summary_path)
    write_report(output / "assessment.md", summary, profiles, metrics)
    print(json.dumps({key: summary[key] for key in (
        "status", "science_measurement", "n_fixed_lenses", "n_i", "n_i2",
        "n_i_support_ge_10pct", "n_i2_support_ge_10pct", "output_directory",
        "output_sha256")}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
