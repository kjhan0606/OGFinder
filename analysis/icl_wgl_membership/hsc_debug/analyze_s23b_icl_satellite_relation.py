#!/usr/bin/env python3
"""Relate S23B relative ICL excess to CAMIRA member-candidate concentration.

The analysis uses only the approved S23B member file, its catalog/WCS-joined
position table, and saved target/control radial profiles. It does not request
or download data and reads no science-image pixels. CAMIRA positions are
unweighted candidate-member counts; the two undocumented trailing catalog
fields are not interpreted. The fitted ICL quantity is a relative
50--300 minus 300--800 proper-kpc target/control contrast, not an absolute
ICL measurement.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from astropy.cosmology import FlatLambdaCDM

from analyze_s23b_icl_spatial_blocks import (
    connected_pair_blocks,
    read_json,
    read_tsv,
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
RELATIVE_TABLE = DATA / "s23b_icl_relative_radial_contrasts_2026-11-07_02/paired_relative_contrasts.tsv"
MEMBER_IMPORT = DATA / "s23b_camira_member_catalog_2026-11-07_01/import_2026-11-07_02"
MEMBER_SUMMARY = MEMBER_IMPORT / "summary.json"
MEMBER_POSITIONS = MEMBER_IMPORT / "target_member_positions.tsv"
MEMBER_COUNTS = MEMBER_IMPORT / "target_member_footprint_counts.tsv"
DEFAULT_OUTPUT = DATA / "s23b_icl_satellite_relation_2026-11-07_03"
BLOCK_HELPER = Path(__file__).resolve().with_name("analyze_s23b_icl_spatial_blocks.py")

EXPECTED_HASHES = {
    str(TARGET_MANIFEST): "e487582ba40cd3e58da67704242db4d7219e381fb08f4b33fed47ce4ff268b44",
    str(CONTROL_MANIFEST): "ae9eb406cfc179d9980031bb67f85f8ee71c4580ba59e97713278823db5c06df",
    str(STACK_SUMMARY): "c8bac8dbdd0a886604c216b2410725eb97883b0cdc1a91fd47e7a7cfcaaaa20b",
    str(PAIR_PROFILES): "285906b03360e9d137727e948c15419767c356d3343e61f409562567ac4abc43",
    str(SPATIAL_SUMMARY): "423bb5e3ad52f3b349dedc74b00a2c05ebb4589e20402abb68b4e99d17199afc",
    str(RELATIVE_SUMMARY): "a47c6679ad318abcbee5ee5ab46e1a984976bd01c6f1dcbbe3d6de99b89306a8",
    str(RELATIVE_TABLE): "c5de5f32b51c6899c8d17fe2b3e4b0a3966d810afd8db5ee27ba649fc0642800",
    str(MEMBER_SUMMARY): "3d1442420aeae051f5a970c6b717a93c69fdaa57cebb096063c46d8afe5f4d4e",
    str(MEMBER_POSITIONS): "d1eb604922146c4d953ba4cad2f92d4c515fb4d4f10cd553f62b683535bfa6ae",
    str(MEMBER_COUNTS): "fcd1eaf30ab419d45b3d4c5438f74dffc098f49d96ae53a7b076a0d28a43a4ee",
}
EXPECTED_BLOCK_HELPER_SHA256 = "d7a21f68988fb782e5bcda147ff9a9040ded1b8c3f6a12ac883598beb2cd0db5"
EXPECTED_ROWS = {"i": 241, "i2": 781}
BLOCK_WIDTHS = (2.5, 5.0, 10.0)
GRID_PHASES = (0.0, 0.5)
CONTRAST_LABEL = "50-300_minus_300-800"
INNER_AREA = math.pi * (300.0**2 - 50.0**2)
OUTER_AREA = math.pi * (800.0**2 - 300.0**2)
N_BOOTSTRAP = 5_000
SEED = 20261107


def verify_inputs() -> tuple[dict, dict, dict, dict, list[dict], list[dict], dict]:
    paths = [Path(path) for path in EXPECTED_HASHES]
    actual = {str(path): sha256(path) for path in paths}
    if actual != EXPECTED_HASHES:
        changed = [path for path in EXPECTED_HASHES if actual.get(path) != EXPECTED_HASHES[path]]
        raise ValueError(f"pinned inputs changed: {changed}")
    if sha256(BLOCK_HELPER) != EXPECTED_BLOCK_HELPER_SHA256:
        raise ValueError("connected sky-block helper changed")
    target = read_json(TARGET_MANIFEST)
    control = read_json(CONTROL_MANIFEST)
    stack = read_json(STACK_SUMMARY)
    spatial = read_json(SPATIAL_SUMMARY)
    relative = read_json(RELATIVE_SUMMARY)
    member = read_json(MEMBER_SUMMARY)
    profiles = read_tsv(PAIR_PROFILES)
    relative_rows = read_tsv(RELATIVE_TABLE)
    member_rows = read_tsv(MEMBER_POSITIONS)
    if (target.get("science_measurement") is not False
            or control.get("science_measurement") is not False
            or stack.get("science_measurement") is not False
            or relative.get("science_measurement") is not False
            or member.get("science_measurement") is not False):
        raise ValueError("an exploratory input unexpectedly claims a science measurement")
    if (member.get("ready_for_exploratory_satellite_analysis") is not True
            or member.get("catalog_join", {}).get("target_clusters_with_member_rows") != 1043):
        raise ValueError("member import is not a complete 1,043-target catalog join")
    if len(target.get("clusters", [])) != 1043 or len(control.get("clusters", [])) != 1041:
        raise ValueError("target/control manifest counts changed")
    if len(profiles) != 1022 or len(member_rows) != 45634:
        raise ValueError("paired profile or target member-position count changed")
    return target, control, stack, spatial, profiles, relative_rows, member


def annular_profile(profile: list[dict], low: int, high: int) -> tuple[float | None, int, int, float]:
    selected = [
        row for row in profile
        if float(row["r_min_kpc"]) >= low and float(row["r_max_kpc"]) <= high
    ]
    if not selected:
        raise ValueError(f"saved profile has no bins within {low}-{high} kpc")
    valid = sum(int(row["n_valid_pixels"]) for row in selected)
    geometric = sum(int(row["geometric_pixels"]) for row in selected)
    if geometric <= 0:
        raise ValueError(f"non-positive geometric pixel count in {low}-{high} kpc")
    if valid == 0:
        return None, 0, geometric, 0.0
    total = sum(
        float(row["delta_mean_njy_arcsec2"]) * int(row["n_valid_pixels"])
        for row in selected if int(row["n_valid_pixels"]) > 0
    )
    return total / valid, valid, geometric, valid / geometric


def read_member_counts(member_rows: list[dict], target_manifest: dict) -> dict[str, list[dict]]:
    target_by_id = {str(row["id"]): row for row in target_manifest["clusters"]}
    cosmology = FlatLambdaCDM(
        H0=float(target_manifest["h0"]),
        Om0=float(target_manifest["omega_m"]),
    )
    grouped = {target_id: [] for target_id in target_by_id}
    for row in member_rows:
        target_id = str(row["target_cluster_id"])
        if target_id not in grouped:
            raise ValueError(f"member position refers to unknown target {target_id}")
        target = target_by_id[target_id]
        z = float(target["z_cl"])
        scale = float(cosmology.kpc_proper_per_arcmin(z).value) / 60.0
        radius_kpc = float(row["cluster_offset_arcsec"]) * scale
        grouped[target_id].append({
            "radius_kpc": radius_kpc,
            "inside_image": str(row["inside_image_pixel_centers"]).lower() == "true",
        })
    if any(not rows for rows in grouped.values()):
        missing = [target_id for target_id, rows in grouped.items() if not rows]
        raise ValueError(f"target member-position groups are empty: {missing[:10]}")
    return grouped


def counts_in_annulus(rows: list[dict], low: float, high: float) -> tuple[int, int, float]:
    selected = [row for row in rows if low <= row["radius_kpc"] < high]
    in_image = sum(bool(row["inside_image"]) for row in selected)
    total = len(selected)
    fraction = in_image / total if total else 1.0
    return total, in_image, fraction


def make_pair_rows(target_manifest: dict, control_manifest: dict, profiles: list[dict],
                   relative_rows: list[dict], member_positions: list[dict]) -> list[dict]:
    targets = {str(row["id"]): row for row in target_manifest["clusters"]}
    controls = {str(row["target_cluster_id"]): row for row in control_manifest["clusters"]}
    members_by_target = read_member_counts(member_positions, target_manifest)
    rel_lookup = {}
    for row in relative_rows:
        key = (str(row["target_id"]), str(row["filter"]), str(row["contrast"]))
        if key in rel_lookup:
            raise ValueError(f"duplicate saved relative contrast {key}")
        rel_lookup[key] = row

    output = []
    seen = set()
    for pair in profiles:
        target_id = str(pair["target_id"])
        band = str(pair["filter"])
        key = (target_id, band)
        if key in seen:
            raise ValueError(f"duplicate target/filter pair {key}")
        seen.add(key)
        if target_id not in targets or target_id not in controls:
            raise ValueError(f"pair has no pinned target/control row: {target_id}")
        target, control = targets[target_id], controls[target_id]
        if (str(control["id"]) != str(pair["control_id"])
                or str(control["target_cluster_id"]) != target_id):
            raise ValueError(f"target/control pair provenance mismatch for {target_id}")
        member_rows = members_by_target[target_id]
        n_inner, n_inner_image, f_inner = counts_in_annulus(member_rows, 50.0, 300.0)
        n_outer, n_outer_image, f_outer = counts_in_annulus(member_rows, 300.0, 800.0)
        density_ratio = math.log10(
            ((n_inner_image + 0.5) / INNER_AREA)
            / ((n_outer_image + 0.5) / OUTER_AREA)
        )
        profile = json.loads(pair["profile_json"])
        inner, n_inner_pix, geom_inner, support_inner = annular_profile(profile, 50, 300)
        outer, n_outer_pix, geom_outer, support_outer = annular_profile(profile, 300, 800)
        if inner is None or outer is None:
            contrast = None
        else:
            contrast = inner - outer
        if contrast is not None:
            expected = rel_lookup.get((target_id, band, CONTRAST_LABEL))
            if expected is None or str(expected.get("science_measurement")).lower() != "false":
                raise ValueError(f"saved contrast table missing {target_id}/{band}")
            saved = float(expected["contrast_njy_arcsec2"])
            if not math.isclose(contrast, saved, rel_tol=0.0, abs_tol=1e-9):
                raise ValueError(f"radial contrast does not reproduce for {target_id}/{band}")
        richness = float(target["richness"])
        z = float(target["z_cl"])
        target_ebv = float(control["target_ebv_sfd"])
        control_ebv = float(control["candidate_ebv_sfd"])
        output.append({
            "target_id": target_id,
            "control_id": str(pair["control_id"]),
            "filter": band,
            "ra_deg": float(target["ra_deg"]),
            "dec_deg": float(target["dec_deg"]),
            "z_cl": z,
            "richness": richness,
            "log10_richness": math.log10(richness),
            "delta_ebv_sfd_target_minus_control": target_ebv - control_ebv,
            "icl_inner_50_300_njy_arcsec2": inner,
            "icl_outer_300_800_njy_arcsec2": outer,
            "icl_relative_contrast_njy_arcsec2": contrast,
            "inner_common_pixels": n_inner_pix,
            "inner_geometric_pixels": geom_inner,
            "inner_common_support_fraction": support_inner,
            "outer_common_pixels": n_outer_pix,
            "outer_geometric_pixels": geom_outer,
            "outer_common_support_fraction": support_outer,
            "camira_member_candidates_50_300_catalog": n_inner,
            "camira_member_candidates_50_300_in_image": n_inner_image,
            "member_position_footprint_fraction_50_300": f_inner,
            "camira_member_candidates_300_800_catalog": n_outer,
            "camira_member_candidates_300_800_in_image": n_outer_image,
            "member_position_footprint_fraction_300_800": f_outer,
            "log10_inner_outer_candidate_density_ratio": density_ratio,
            "science_measurement": False,
        })
    if len(output) != 1022:
        raise ValueError(f"expected 1,022 paired profiles, found {len(output)}")
    valid = [row for row in output if row["icl_relative_contrast_njy_arcsec2"] is not None]
    counts = {band: sum(row["filter"] == band for row in valid) for band in ("i", "i2")}
    if counts != {"i": 232, "i2": 756}:
        raise ValueError(f"relative-contrast support changed: {counts}")
    return output


def model_matrix(rows: list[dict]) -> tuple[np.ndarray, np.ndarray, dict]:
    y = np.asarray([row["icl_relative_contrast_njy_arcsec2"] for row in rows], dtype=float)
    x = np.asarray([row["log10_inner_outer_candidate_density_ratio"] for row in rows], dtype=float)
    richness = np.asarray([row["log10_richness"] for row in rows], dtype=float)
    redshift = np.asarray([row["z_cl"] for row in rows], dtype=float)
    delta_ebv = np.asarray([row["delta_ebv_sfd_target_minus_control"] for row in rows], dtype=float)
    common_support = np.minimum(
        np.asarray([row["inner_common_support_fraction"] for row in rows], dtype=float),
        np.asarray([row["outer_common_support_fraction"] for row in rows], dtype=float),
    )
    member_footprint = np.minimum(
        np.asarray([row["member_position_footprint_fraction_50_300"] for row in rows], dtype=float),
        np.asarray([row["member_position_footprint_fraction_300_800"] for row in rows], dtype=float),
    )
    columns = [np.ones(len(rows)), x - x.mean()]
    names = ["intercept", "satellite_density_contrast"]
    metadata = {
        "satellite_density_contrast_center": float(x.mean()),
        "log10_richness_center": float(richness.mean()),
        "redshift_center": float(redshift.mean()),
        "redshift_unit": 0.1,
        "delta_ebv_center": float(delta_ebv.mean()),
        "delta_ebv_unit": 0.01,
        "minimum_common_support_center": float(common_support.mean()),
        "minimum_member_footprint_fraction_center": float(member_footprint.mean()),
    }
    candidates = [
        ("log10_richness", richness - richness.mean()),
        ("redshift_per_0p1", (redshift - redshift.mean()) / 0.1),
        ("delta_ebv_per_0p01", (delta_ebv - delta_ebv.mean()) / 0.01),
        ("minimum_common_support_fraction", common_support - common_support.mean()),
        ("minimum_member_position_footprint_fraction", member_footprint - member_footprint.mean()),
    ]
    for name, values in candidates:
        if float(np.ptp(values)) > 1e-10:
            columns.append(values)
            names.append(name)
    X = np.column_stack(columns)
    if np.linalg.matrix_rank(X) != X.shape[1]:
        raise ValueError(f"rank-deficient adjusted model: {names}")
    metadata["model_columns"] = names
    return X, y, metadata


def block_bootstrap_slope(X: np.ndarray, y: np.ndarray, labels: np.ndarray,
                          replicas: int, seed: int) -> tuple[float, float, float, int, int]:
    groups = np.unique(labels)
    if len(groups) < X.shape[1] + 2:
        raise ValueError(f"too few connected blocks ({len(groups)}) for {X.shape[1]} model columns")
    full, _, rank, _ = np.linalg.lstsq(X, y, rcond=None)
    if rank != X.shape[1]:
        raise ValueError("full-sample adjusted model is rank deficient")
    group_index = {label: i for i, label in enumerate(groups)}
    row_group = np.asarray([group_index[label] for label in labels], dtype=int)
    rng = np.random.default_rng(seed)
    draws = rng.multinomial(len(groups), np.full(len(groups), 1.0 / len(groups)), size=replicas)
    weights = draws[:, row_group].astype(float, copy=False)
    n_param = X.shape[1]
    gram = np.empty((replicas, n_param, n_param), dtype=float)
    for left in range(n_param):
        for right in range(left, n_param):
            values = weights @ (X[:, left] * X[:, right])
            gram[:, left, right] = values
            gram[:, right, left] = values
    rhs = np.column_stack([weights @ (X[:, index] * y) for index in range(n_param)])
    condition = np.linalg.cond(gram)
    good = np.isfinite(condition) & (condition < 1e10)
    beta = np.full((replicas, n_param), np.nan, dtype=float)
    if np.any(good):
        beta[good] = np.linalg.solve(gram[good], rhs[good, :, None])[..., 0]
    slopes = beta[:, 1]
    slopes = slopes[np.isfinite(slopes)]
    if slopes.size < replicas * 0.9:
        raise ValueError(f"only {slopes.size}/{replicas} valid connected-block fits")
    low, high = np.percentile(slopes, [2.5, 97.5])
    return float(full[1]), float(low), float(high), int(len(groups)), int(slopes.size)


def estimate(rows_by_band: dict[str, list[dict]], spatial: dict,
             replicas: int) -> tuple[list[dict], dict]:
    results = []
    fit_cache = {}
    origin = float(spatial["ra_grid_origin_deg"])
    if not math.isclose(origin, 84.056247280625, rel_tol=0.0, abs_tol=1e-9):
        raise ValueError("S23B spatial-grid origin differs from pinned audit")
    for width_index, width in enumerate(BLOCK_WIDTHS):
        for phase_index, phase in enumerate(GRID_PHASES):
            for band_index, band in enumerate(("i", "i2")):
                rows = rows_by_band[band]
                controls = {
                    row["target_id"]: {
                        "ra_deg": row["control_ra_deg"],
                        "dec_deg": row["control_dec_deg"],
                    } for row in rows
                }
                labels, blocks = connected_pair_blocks(
                    rows, controls, width, phase, origin
                )
                X, y, metadata = model_matrix(rows)
                slope, low, high, n_blocks, valid = block_bootstrap_slope(
                    X, y, labels, replicas,
                    SEED + 1000 * width_index + 100 * phase_index + band_index,
                )
                row = {
                    "filter": band,
                    "block_width_deg": width,
                    "grid_phase_fraction": phase,
                    "n_pairs": len(rows),
                    "n_connected_blocks": n_blocks,
                    "n_bootstrap_requested": replicas,
                    "n_bootstrap_valid": valid,
                    "density_ratio_slope_njy_arcsec2_per_dex": slope,
                    "block_bootstrap_ci95_low": low,
                    "block_bootstrap_ci95_high": high,
                    "ci_excludes_zero": bool(low > 0 or high < 0),
                    "model_columns": ",".join(metadata["model_columns"]),
                    "science_measurement": False,
                }
                results.append(row)
                if width == 5.0 and phase == 0.0:
                    fit_cache[band] = (rows, X, y, metadata, slope, low, high)
    return results, fit_cache


def render_plot(path: Path, fit_cache: dict) -> None:
    figure, axes = plt.subplots(1, 2, figsize=(11.5, 5.0), sharey=True)
    for axis, band in zip(axes, ("i", "i2")):
        rows, X, y, metadata, slope, low, high = fit_cache[band]
        x = X[:, 1] + metadata["satellite_density_contrast_center"]
        nuisance = X[:, [0] + list(range(2, X.shape[1]))]
        x_resid = x - nuisance @ np.linalg.lstsq(nuisance, x, rcond=None)[0]
        y_resid = y - nuisance @ np.linalg.lstsq(nuisance, y, rcond=None)[0]
        axis.scatter(x_resid, y_resid, s=14, alpha=0.35, color="#1f5a85",
                     linewidths=0)
        xline = np.linspace(float(np.min(x_resid)), float(np.max(x_resid)), 100)
        axis.plot(xline, slope * xline, color="#b04a24", linewidth=1.8)
        axis.axhline(0.0, color="0.4", linewidth=0.8)
        axis.axvline(0.0, color="0.4", linewidth=0.8)
        axis.text(
            0.04, 0.96,
            f"HSC-{band}   N={len(rows)}\n"
            f"adjusted slope={slope:.2f} [{low:.2f}, {high:.2f}]",
            transform=axis.transAxes, ha="left", va="top", fontsize=10,
        )
        axis.grid(alpha=0.18)
        axis.set_xlabel("Residual log satellite-density ratio")
    axes[0].set_ylabel("Residual relative ICL contrast [nJy arcsec$^{-2}$]")
    figure.suptitle(
        "CAMIRA member-candidate concentration vs relative ICL",
        y=0.995, fontsize=11,
    )
    figure.tight_layout(rect=(0, 0, 1, 0.97))
    figure.savefig(path, dpi=180)
    plt.close(figure)


def leave_one_block_out(rows_by_band: dict[str, list[dict]], spatial: dict,
                        replicas: int) -> tuple[list[dict], dict]:
    output = []
    summaries = {}
    origin = float(spatial["ra_grid_origin_deg"])
    for band in ("i", "i2"):
        rows = rows_by_band[band]
        controls = {
            row["target_id"]: {
                "ra_deg": row["control_ra_deg"],
                "dec_deg": row["control_dec_deg"],
            } for row in rows
        }
        labels, _ = connected_pair_blocks(rows, controls, 5.0, 0.0, origin)
        X, y, _ = model_matrix(rows)
        full = np.linalg.lstsq(X, y, rcond=None)[0][1]
        deletions = []
        for block in np.unique(labels):
            keep = labels != block
            beta, _, rank, _ = np.linalg.lstsq(X[keep], y[keep], rcond=None)
            if rank != X.shape[1]:
                raise ValueError(f"leave-one-block model is rank deficient: {band}/{block}")
            row = {
                "filter": band,
                "block_width_deg": 5.0,
                "grid_phase_fraction": 0.0,
                "deleted_block_id": str(block),
                "n_pairs_retained": int(np.count_nonzero(keep)),
                "density_ratio_slope_njy_arcsec2_per_dex": float(beta[1]),
                "shift_from_full_njy_arcsec2_per_dex": float(beta[1] - full),
                "sign_changed": bool(np.sign(beta[1]) != np.sign(full)),
            }
            output.append(row)
            deletions.append(row)
        summaries[band] = {
            "n_blocks": len(deletions),
            "n_deleted_block_slopes_below_zero": sum(
                row["density_ratio_slope_njy_arcsec2_per_dex"] < 0 for row in deletions
            ),
            "n_sign_changes": sum(row["sign_changed"] for row in deletions),
            "minimum_deleted_block_slope": min(
                row["density_ratio_slope_njy_arcsec2_per_dex"] for row in deletions
            ),
            "maximum_deleted_block_slope": max(
                row["density_ratio_slope_njy_arcsec2_per_dex"] for row in deletions
            ),
            "largest_absolute_slope_shift": max(
                abs(row["shift_from_full_njy_arcsec2_per_dex"]) for row in deletions
            ),
        }
        residual = y - X @ np.linalg.lstsq(X, y, rcond=None)[0]
        top = np.argsort(np.abs(residual))[::-1][:5]
        summaries[band]["largest_absolute_model_residuals"] = [
            {
                "target_id": rows[index]["target_id"],
                "residual_njy_arcsec2": float(residual[index]),
                "relative_icl_contrast_njy_arcsec2":
                    float(rows[index]["icl_relative_contrast_njy_arcsec2"]),
            }
            for index in top
        ]
    return output, summaries


def support_floor_sensitivity(rows_by_band: dict[str, list[dict]], spatial: dict,
                              replicas: int) -> list[dict]:
    output = []
    origin = float(spatial["ra_grid_origin_deg"])
    floors = (0.05, 0.10, 0.20)
    for floor_index, floor in enumerate(floors):
        for width_index, width in enumerate(BLOCK_WIDTHS):
            for phase_index, phase in enumerate(GRID_PHASES):
                for band_index, band in enumerate(("i", "i2")):
                    rows = [
                        row for row in rows_by_band[band]
                        if min(row["inner_common_support_fraction"],
                               row["outer_common_support_fraction"]) >= floor
                    ]
                    row = {
                        "filter": band,
                        "minimum_common_support_fraction_both_annuli": floor,
                        "block_width_deg": width,
                        "grid_phase_fraction": phase,
                        "n_pairs": len(rows),
                        "n_connected_blocks": None,
                        "n_bootstrap_requested": replicas,
                        "n_bootstrap_valid": 0,
                        "density_ratio_slope_njy_arcsec2_per_dex": None,
                        "block_bootstrap_ci95_low": None,
                        "block_bootstrap_ci95_high": None,
                        "reason": "",
                        "status": "skipped",
                        "science_measurement": False,
                    }
                    try:
                        if not rows:
                            raise ValueError("no pairs meet the support floor")
                        controls = {
                            item["target_id"]: {
                                "ra_deg": item["control_ra_deg"],
                                "dec_deg": item["control_dec_deg"],
                            } for item in rows
                        }
                        labels, _ = connected_pair_blocks(
                            rows, controls, width, phase, origin
                        )
                        row["n_connected_blocks"] = int(np.unique(labels).size)
                        X, y, metadata = model_matrix(rows)
                        slope, low, high, n_blocks, valid = block_bootstrap_slope(
                            X, y, labels, replicas,
                            SEED + 10_000 + 1000 * floor_index
                            + 100 * width_index + 10 * phase_index + band_index,
                        )
                        row.update({
                            "n_connected_blocks": n_blocks,
                            "n_bootstrap_valid": valid,
                            "density_ratio_slope_njy_arcsec2_per_dex": slope,
                            "block_bootstrap_ci95_low": low,
                            "block_bootstrap_ci95_high": high,
                            "ci_excludes_zero": bool(low > 0 or high < 0),
                            "model_columns": ",".join(metadata["model_columns"]),
                            "status": "complete",
                        })
                    except ValueError as error:
                        row["status"] = "not_estimable"
                        row["reason"] = str(error)
                    output.append(row)
    return output


def run(output: Path, replicas: int) -> dict:
    if output.exists():
        raise FileExistsError(f"refusing to overwrite {output}")
    target, control, stack, spatial, profiles, relative_rows, member_summary = verify_inputs()
    member_positions = read_tsv(MEMBER_POSITIONS)
    pair_rows = make_pair_rows(
        target, control, profiles, relative_rows, member_positions
    )
    valid_rows = [
        row for row in pair_rows
        if row["icl_relative_contrast_njy_arcsec2"] is not None
    ]
    rows_by_band = {
        band: [row for row in valid_rows if row["filter"] == band]
        for band in ("i", "i2")
    }
    for band, expected in (("i", 232), ("i2", 756)):
        if len(rows_by_band[band]) != expected:
            raise ValueError(f"unexpected supported pair count for {band}")
    for row in pair_rows:
        target_row = next(x for x in target["clusters"] if str(x["id"]) == row["target_id"])
        control_row = next(x for x in control["clusters"]
                           if str(x["target_cluster_id"]) == row["target_id"])
        row["control_ra_deg"] = float(control_row["ra_deg"])
        row["control_dec_deg"] = float(control_row["dec_deg"])
    results, fit_cache = estimate(rows_by_band, spatial, replicas)
    deletion_rows, deletion_summary = leave_one_block_out(rows_by_band, spatial, replicas)
    support_rows = support_floor_sensitivity(rows_by_band, spatial, replicas)

    output.mkdir(parents=True, exist_ok=False)
    pair_table = output / "paired_icl_satellite_profiles.tsv"
    slope_table = output / "connected_block_slope_sensitivity.tsv"
    deletion_table = output / "primary_5deg_leave_one_block_out.tsv"
    support_table = output / "common_support_floor_sensitivity.tsv"
    figure_path = output / "adjusted_icl_satellite_relation.png"
    write_tsv(pair_table, valid_rows)
    write_tsv(slope_table, results)
    write_tsv(deletion_table, deletion_rows)
    write_tsv(support_table, support_rows)
    render_plot(figure_path, fit_cache)

    primary = {
        band: next(row for row in results
                   if row["filter"] == band
                   and row["block_width_deg"] == 5.0
                   and row["grid_phase_fraction"] == 0.0)
        for band in ("i", "i2")
    }
    valid_support = {}
    for band in ("i", "i2"):
        rows = rows_by_band[band]
        valid_support[band] = {
            "n_pairs": len(rows),
            "median_inner_common_support_fraction": float(np.median([
                row["inner_common_support_fraction"] for row in rows
            ])),
            "median_outer_common_support_fraction": float(np.median([
                row["outer_common_support_fraction"] for row in rows
            ])),
            "median_member_position_footprint_fraction_50_300": float(np.median([
                row["member_position_footprint_fraction_50_300"] for row in rows
            ])),
            "median_member_position_footprint_fraction_300_800": float(np.median([
                row["member_position_footprint_fraction_300_800"] for row in rows
            ])),
            "median_candidates_50_300_in_image": float(np.median([
                row["camira_member_candidates_50_300_in_image"] for row in rows
            ])),
            "median_candidates_300_800_in_image": float(np.median([
                row["camira_member_candidates_300_800_in_image"] for row in rows
            ])),
        }
    summary = {
        "schema": "astrafek-s23b-icl-satellite-relation-v1",
        "status": "complete_exploratory_direct_member_position_relation",
        "science_measurement": False,
        "estimator": {
            "icl_outcome": "common-support target-minus-control mean at 50-300 kpc minus that at 300-800 kpc",
            "satellite_predictor": "log10 of (50-300 kpc CAMIRA candidate-member surface density / 300-800 kpc density), with 0.5-count continuity correction",
            "member_weights": "unweighted candidate counts; headerless trailing fields are not interpreted",
            "regression_adjustments": [
                "log10 CAMIRA richness",
                "cluster photometric redshift",
                "target-minus-control SFD E(B-V) mismatch",
                "minimum of inner/outer common-pixel support fractions",
                "minimum of inner/outer member-position image-footprint fractions when nonconstant",
            ],
            "uncertainty": "5,000 connected RA-Dec sky-block bootstrap resamples per filter and each 2.5/5/10 degree layout with two grid phases",
            "primary_layout": "5 degree, phase 0",
            "h0": float(target["h0"]),
            "omega_m": float(target["omega_m"]),
        },
        "sample": {
            "all_filter_matched_pairs": len(pair_rows),
            "supported_pairs_by_filter": {band: len(rows_by_band[band]) for band in ("i", "i2")},
            "total_target_member_positions_imported": member_summary["catalog_join"]["target_member_rows"],
            "science_image_pixels_read": False,
            "network_used": False,
        },
        "supported_sample_diagnostics": valid_support,
        "primary_layout_results": primary,
        "spatial_partition_sensitivity": {
            band: {
                "n_layouts": sum(row["filter"] == band for row in results),
                "n_layouts_ci_excluding_zero": sum(
                    row["filter"] == band and row["ci_excludes_zero"] for row in results
                ),
                "ci_excludes_zero_in_all_layouts": all(
                    row["ci_excludes_zero"] for row in results if row["filter"] == band
                ),
                "slope_point_estimate_is_layout_invariant": True,
            }
            for band in ("i", "i2")
        },
        "primary_layout_leave_one_block_out": deletion_summary,
        "common_support_floor_sensitivity_primary_layout": {
            band: [
                row for row in support_rows
                if row["filter"] == band
                and row["block_width_deg"] == 5.0
                and row["grid_phase_fraction"] == 0.0
            ]
            for band in ("i", "i2")
        },
        "limitations": [
            "CAMIRA member candidates are photometric catalog positions, not spectroscopically confirmed satellites.",
            "Candidate counts are unweighted; their completeness and membership-probability calibration are not established here.",
            "The ICL outcome is a relative radial contrast and the 300-800 kpc reference may contain cluster light.",
            "The source catalog positions inside a FITS footprint are not necessarily on usable, unmasked pixels.",
            "The ICL profile has sparse common-pixel support; the linear support adjustment does not correct selection bias.",
            "Residual sky, extinction, extended PSF wings, full coadd covariance, centering, and mask sensitivity remain unresolved.",
        ],
        "inputs_sha256": {
            **EXPECTED_HASHES,
            str(BLOCK_HELPER): EXPECTED_BLOCK_HELPER_SHA256,
        },
        "output_sha256": {
            pair_table.name: sha256(pair_table),
            slope_table.name: sha256(slope_table),
            deletion_table.name: sha256(deletion_table),
            support_table.name: sha256(support_table),
            figure_path.name: sha256(figure_path),
        },
        "analysis_script_sha256": sha256(Path(__file__).resolve()),
        "bootstrap_replicas_per_layout": replicas,
    }
    summary_path = output / "summary.json"
    summary_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2, allow_nan=False) + "\n",
                             encoding="utf-8")
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--bootstrap-replicas", type=int, default=N_BOOTSTRAP)
    args = parser.parse_args()
    if args.bootstrap_replicas < 1000:
        parser.error("use at least 1,000 connected-block resamples")
    summary = run(args.output, args.bootstrap_replicas)
    print(json.dumps(summary, ensure_ascii=False, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as error:
        print(f"analyze_s23b_icl_satellite_relation: {type(error).__name__}: {error}",
              file=sys.stderr)
        sys.exit(1)
