#!/usr/bin/env python3
"""Measure a raw CAMIRA member-field-8 weighting sensitivity for S23B.

This offline companion repeats the saved direct ICL--candidate-density model,
replacing unweighted in-image candidate counts with sums of the approved
member file's eighth raw field. Public CAMIRA literature makes this field
consistent with the membership weight, but the exact headerless file column
order is not documented; the result is therefore a field-8 sensitivity, not
a calibrated membership-weighted science measurement. No network or image
pixels are used.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from astropy.cosmology import FlatLambdaCDM

import analyze_s23b_icl_satellite_relation as base


DEFAULT_OUTPUT = base.DATA / "s23b_icl_satellite_field8_weighted_density_2026-11-07_02"
WEIGHT_COLUMN = "member_field8_raw"
INNER_AREA = base.INNER_AREA
OUTER_AREA = base.OUTER_AREA


def weighted_annulus_sums(member_positions: list[dict], target_manifest: dict) -> dict[str, dict]:
    """Sum raw field 8 for in-footprint candidates in the two fixed annuli."""
    cosmology = FlatLambdaCDM(
        H0=float(target_manifest["h0"]), Om0=float(target_manifest["omega_m"])
    )
    target_rows = {str(row["id"]): row for row in target_manifest["clusters"]}
    output = {
        target_id: {
            "n_inner_catalog": 0, "n_inner_in_image": 0, "sum_field8_inner": 0.0,
            "n_outer_catalog": 0, "n_outer_in_image": 0, "sum_field8_outer": 0.0,
        }
        for target_id in target_rows
    }
    scales = {
        target_id: float(cosmology.kpc_proper_per_arcmin(float(row["z_cl"])).value) / 60.0
        for target_id, row in target_rows.items()
    }
    for row in member_positions:
        target_id = str(row["target_cluster_id"])
        if target_id not in output:
            raise ValueError(f"unknown target in member table: {target_id}")
        radius = float(row["cluster_offset_arcsec"]) * scales[target_id]
        annulus = "inner" if 50.0 <= radius < 300.0 else "outer" if 300.0 <= radius < 800.0 else None
        if annulus is None:
            continue
        record = output[target_id]
        record[f"n_{annulus}_catalog"] += 1
        if str(row["inside_image_pixel_centers"]).lower() != "true":
            continue
        weight = float(row[WEIGHT_COLUMN])
        if not math.isfinite(weight) or weight < 0.001 or weight > 1.0:
            raise ValueError(f"raw field 8 outside published-threshold range: {weight}")
        record[f"n_{annulus}_in_image"] += 1
        record[f"sum_field8_{annulus}"] += weight
    return output


def weighted_predictor(inner_sum: float, outer_sum: float) -> float:
    """Return the same annular density ratio with 0.5 pseudo-counts."""
    return math.log10(((inner_sum + 0.5) / INNER_AREA) /
                      ((outer_sum + 0.5) / OUTER_AREA))


def add_weighted_predictor(pair_rows: list[dict], member_positions: list[dict],
                           target_manifest: dict) -> list[dict]:
    sums = weighted_annulus_sums(member_positions, target_manifest)
    output = []
    for source in pair_rows:
        row = dict(source)
        target_id = str(row["target_id"])
        member = sums[target_id]
        row["unweighted_log10_inner_outer_density_ratio"] = row[
            "log10_inner_outer_candidate_density_ratio"
        ]
        row["weighted_log10_inner_outer_density_ratio_field8"] = weighted_predictor(
            member["sum_field8_inner"], member["sum_field8_outer"]
        )
        row["field8_weight_sum_50_300_in_image"] = member["sum_field8_inner"]
        row["field8_weight_sum_300_800_in_image"] = member["sum_field8_outer"]
        row["field8_catalog_candidate_count_50_300"] = member["n_inner_catalog"]
        row["field8_in_image_candidate_count_50_300"] = member["n_inner_in_image"]
        row["field8_catalog_candidate_count_300_800"] = member["n_outer_catalog"]
        row["field8_in_image_candidate_count_300_800"] = member["n_outer_in_image"]
        output.append(row)
    return output


def fit_block_resamples(X: np.ndarray, y: np.ndarray, labels: np.ndarray,
                        block_draws: np.ndarray) -> tuple[float, np.ndarray]:
    groups = np.unique(labels)
    if len(groups) < X.shape[1] + 2:
        raise ValueError(f"too few connected blocks ({len(groups)}) for {X.shape[1]} columns")
    full, _, rank, _ = np.linalg.lstsq(X, y, rcond=None)
    if rank != X.shape[1]:
        raise ValueError("full-sample model is rank deficient")
    if block_draws.ndim != 2 or block_draws.shape[1] != len(groups):
        raise ValueError("block-resample matrix does not match the connected-block labels")
    group_index = {label: index for index, label in enumerate(groups)}
    row_group = np.asarray([group_index[label] for label in labels], dtype=int)
    replicas = block_draws.shape[0]
    weights = block_draws[:, row_group].astype(float, copy=False)
    n_param = X.shape[1]
    gram = np.empty((replicas, n_param, n_param), dtype=float)
    for left in range(n_param):
        for right in range(left, n_param):
            value = weights @ (X[:, left] * X[:, right])
            gram[:, left, right] = value
            gram[:, right, left] = value
    rhs = np.column_stack([weights @ (X[:, index] * y) for index in range(n_param)])
    condition = np.linalg.cond(gram)
    good = np.isfinite(condition) & (condition < 1e10)
    beta = np.full((replicas, n_param), np.nan, dtype=float)
    if np.any(good):
        beta[good] = np.linalg.solve(gram[good], rhs[good, :, None])[..., 0]
    slopes = beta[:, 1]
    slopes = slopes[np.isfinite(slopes)]
    if slopes.size < replicas * 0.9:
        raise ValueError(f"only {slopes.size}/{replicas} valid block-bootstrap fits")
    return float(full[1]), slopes


def paired_layout_fits(rows_by_filter: dict[str, list[dict]], spatial: dict,
                       replicas: int, support_floor: float | None = None) -> list[dict]:
    output = []
    origin = float(spatial["ra_grid_origin_deg"])
    if not math.isclose(origin, 84.056247280625, rel_tol=0.0, abs_tol=1e-9):
        raise ValueError("spatial-grid origin differs from the pinned direct-relation analysis")
    for width_index, width in enumerate(base.BLOCK_WIDTHS):
        for phase_index, phase in enumerate(base.GRID_PHASES):
            for band_index, band in enumerate(("i", "i2")):
                rows = rows_by_filter[band]
                try:
                    controls = {row["target_id"]: {"ra_deg": row["control_ra_deg"],
                                                   "dec_deg": row["control_dec_deg"]}
                                for row in rows}
                    labels, _ = base.connected_pair_blocks(rows, controls, width, phase, origin)
                    weighted_rows = []
                    unweighted_rows = []
                    for row in rows:
                        weighted = dict(row)
                        weighted["log10_inner_outer_candidate_density_ratio"] = row[
                            "weighted_log10_inner_outer_density_ratio_field8"
                        ]
                        unweighted = dict(row)
                        unweighted["log10_inner_outer_candidate_density_ratio"] = row[
                            "unweighted_log10_inner_outer_density_ratio"
                        ]
                        weighted_rows.append(weighted)
                        unweighted_rows.append(unweighted)
                    X_weighted, y, meta = base.model_matrix(weighted_rows)
                    X_unweighted, y_unweighted, meta_unweighted = base.model_matrix(unweighted_rows)
                    if (not np.array_equal(y, y_unweighted)
                            or meta["model_columns"] != meta_unweighted["model_columns"]):
                        raise ValueError("weighted/unweighted fits do not share outcomes and controls")
                    seed = base.SEED + 1000 * width_index + 100 * phase_index + band_index
                    rng = np.random.default_rng(seed)
                    n_blocks = len(np.unique(labels))
                    block_draws = rng.multinomial(
                        n_blocks, np.full(n_blocks, 1.0 / n_blocks), size=replicas
                    )
                    slope_weighted, draws_weighted = fit_block_resamples(
                        X_weighted, y, labels, block_draws
                    )
                    slope_unweighted, draws_unweighted = fit_block_resamples(
                        X_unweighted, y, labels, block_draws
                    )
                    delta = draws_weighted - draws_unweighted
                    low_w, high_w = np.percentile(draws_weighted, [2.5, 97.5])
                    low_u, high_u = np.percentile(draws_unweighted, [2.5, 97.5])
                    low_d, high_d = np.percentile(delta, [2.5, 97.5])
                    output.append({
                        "filter": band,
                        "minimum_common_support_fraction_both_annuli": support_floor,
                        "block_width_deg": width,
                        "grid_phase_fraction": phase,
                        "n_pairs": len(rows),
                        "n_connected_blocks": int(np.unique(labels).size),
                        "bootstrap_replicas": replicas,
                        "weighted_slope_field8_njy_arcsec2_per_dex": slope_weighted,
                        "weighted_ci95_low": float(low_w),
                        "weighted_ci95_high": float(high_w),
                        "weighted_ci_excludes_zero": bool(low_w > 0 or high_w < 0),
                        "unweighted_slope_njy_arcsec2_per_dex": slope_unweighted,
                        "unweighted_ci95_low": float(low_u),
                        "unweighted_ci95_high": float(high_u),
                        "paired_weighted_minus_unweighted_slope": slope_weighted - slope_unweighted,
                        "paired_delta_ci95_low": float(low_d),
                        "paired_delta_ci95_high": float(high_d),
                        "model_columns": ",".join(meta["model_columns"]),
                        "status": "complete",
                        "reason": "",
                        "science_measurement": False,
                    })
                except (ValueError, np.linalg.LinAlgError) as error:
                    output.append({
                        "filter": band,
                        "minimum_common_support_fraction_both_annuli": support_floor,
                        "block_width_deg": width,
                        "grid_phase_fraction": phase,
                        "n_pairs": len(rows),
                        "n_connected_blocks": None,
                        "bootstrap_replicas": replicas,
                        "weighted_slope_field8_njy_arcsec2_per_dex": None,
                        "weighted_ci95_low": None,
                        "weighted_ci95_high": None,
                        "weighted_ci_excludes_zero": None,
                        "unweighted_slope_njy_arcsec2_per_dex": None,
                        "unweighted_ci95_low": None,
                        "unweighted_ci95_high": None,
                        "paired_weighted_minus_unweighted_slope": None,
                        "paired_delta_ci95_low": None,
                        "paired_delta_ci95_high": None,
                        "model_columns": "",
                        "status": "not_estimable",
                        "reason": str(error),
                        "science_measurement": False,
                    })
    return output


def support_floor_fits(rows_by_filter: dict[str, list[dict]], spatial: dict,
                       replicas: int) -> list[dict]:
    output = []
    for floor in (0.05, 0.10, 0.20):
        selected = {
            band: [row for row in rows if min(
                float(row["inner_common_support_fraction"]),
                float(row["outer_common_support_fraction"]),
            ) >= floor]
            for band, rows in rows_by_filter.items()
        }
        fit_rows = paired_layout_fits(selected, spatial, replicas, support_floor=floor)
        for row in fit_rows:
            output.append(row)
    return output


def render_plot(path: Path, rows_by_filter: dict[str, list[dict]], layout_rows: list[dict]) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 5.0), sharey=True)
    for ax, band in zip(axes, ("i", "i2")):
        rows = rows_by_filter[band]
        weighted = [dict(row, log10_inner_outer_candidate_density_ratio=
                         row["weighted_log10_inner_outer_density_ratio_field8"]) for row in rows]
        X, y, meta = base.model_matrix(weighted)
        x = X[:, 1] + meta["satellite_density_contrast_center"]
        nuisance = X[:, [0] + list(range(2, X.shape[1]))]
        x_resid = x - nuisance @ np.linalg.lstsq(nuisance, x, rcond=None)[0]
        y_resid = y - nuisance @ np.linalg.lstsq(nuisance, y, rcond=None)[0]
        ax.scatter(x_resid, y_resid, s=14, alpha=0.35, color="#1f5a85", linewidths=0)
        primary = next(row for row in layout_rows if row["filter"] == band
                       and row["block_width_deg"] == 5.0
                       and row["grid_phase_fraction"] == 0.0)
        xline = np.linspace(float(np.min(x_resid)), float(np.max(x_resid)), 100)
        ax.plot(xline, primary["weighted_slope_field8_njy_arcsec2_per_dex"] * xline,
                color="#b04a24", linewidth=1.8)
        ax.axhline(0.0, color="0.4", linewidth=0.8)
        ax.axvline(0.0, color="0.4", linewidth=0.8)
        ax.text(0.04, 0.96,
                f"HSC-{band}   N={len(rows)}\n"
                f"field-8 slope={primary['weighted_slope_field8_njy_arcsec2_per_dex']:.2f} "
                f"[{primary['weighted_ci95_low']:.2f}, {primary['weighted_ci95_high']:.2f}]",
                transform=ax.transAxes, ha="left", va="top", fontsize=10)
        ax.grid(alpha=0.18)
        ax.set_xlabel("Residual log field-8-weighted density ratio")
    axes[0].set_ylabel("Residual relative ICL contrast [nJy arcsec$^{-2}$]")
    fig.suptitle("Raw CAMIRA member-field-8 density sensitivity", y=0.995, fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    fig.savefig(path, dpi=180)
    plt.close(fig)


def run(output: Path, replicas: int) -> dict:
    if output.exists():
        raise FileExistsError(f"refusing to overwrite {output}")
    target, control, _stack, spatial, profiles, relative_rows, member_summary = base.verify_inputs()
    member_positions = base.read_tsv(base.MEMBER_POSITIONS)
    pair_rows = base.make_pair_rows(target, control, profiles, relative_rows, member_positions)
    pair_rows = add_weighted_predictor(pair_rows, member_positions, target)
    controls = {str(row["target_cluster_id"]): row for row in control["clusters"]}
    for row in pair_rows:
        match = controls[row["target_id"]]
        row["control_ra_deg"] = float(match["ra_deg"])
        row["control_dec_deg"] = float(match["dec_deg"])
    valid_rows = [row for row in pair_rows if row["icl_relative_contrast_njy_arcsec2"] is not None]
    rows_by_filter = {band: [row for row in valid_rows if row["filter"] == band]
                      for band in ("i", "i2")}
    if {band: len(rows) for band, rows in rows_by_filter.items()} != {"i": 232, "i2": 756}:
        raise ValueError("weighted sensitivity must retain the exact unweighted supported sample")

    layout_rows = paired_layout_fits(rows_by_filter, spatial, replicas)
    support_rows = support_floor_fits(rows_by_filter, spatial, replicas)
    output.mkdir(parents=True, exist_ok=False)
    pair_path = output / "paired_field8_weighted_density.tsv"
    layout_path = output / "paired_layout_sensitivity.tsv"
    support_path = output / "common_support_floor_sensitivity.tsv"
    plot_path = output / "field8_weighted_density_relation.png"
    base.write_tsv(pair_path, valid_rows)
    base.write_tsv(layout_path, layout_rows)
    base.write_tsv(support_path, support_rows)
    render_plot(plot_path, rows_by_filter, layout_rows)
    primary = {band: next(row for row in layout_rows if row["filter"] == band
                          and row["block_width_deg"] == 5.0
                          and row["grid_phase_fraction"] == 0.0)
               for band in ("i", "i2")}
    source_file = base.DATA / "s23b_camira_member_catalog_2026-11-07_01/camira_s23b_wide_sm_member.dat.gz"
    summary = {
        "schema": "astrafek-s23b-icl-satellite-field8-weighted-density-v1",
        "status": "complete_exploratory_raw_field8_weight_sensitivity",
        "science_measurement": False,
        "field_interpretation": {
            "raw_field_used": 8,
            "likely_semantics": "CAMIRA membership weight w_mem, based on range and published threshold consistency",
            "file_specific_column_order_documented": False,
            "stored_name": "member_field8_raw",
            "published_threshold_check": "all imported rows are in [0.001000069, 1.0], consistent with w_mem >= 10^-3",
        },
        "estimator": {
            "icl_outcome": "same common-support target-minus-control contrast: 50-300 minus 300-800 proper kpc",
            "predictor": "log10 ratio of inner/outer sums of raw field 8 for in-image CAMIRA candidates, divided by geometric annular areas, with 0.5 pseudo-counts",
            "unweighted_comparator": "the same sample, count-density ratio, 0.5-count correction, and adjustment set as analyze_s23b_icl_satellite_relation.py",
            "adjustments": ["log10 richness", "cluster redshift", "target-control SFD E(B-V) mismatch", "minimum ICL common-pixel support", "minimum member image-footprint fraction"],
            "uncertainty": f"paired connected RA-Dec block bootstrap, {replicas} replicates per 2.5/5/10 degree layout and two phases",
        },
        "sample": {
            "supported_pairs_by_filter": {band: len(rows) for band, rows in rows_by_filter.items()},
            "image_pixels_read": False,
            "network_used": False,
            "new_data_downloaded": False,
            "member_rows": member_summary["catalog_join"]["target_member_rows"],
        },
        "primary_layout_results": primary,
        "common_support_floor_results": {
            str(floor): [row for row in support_rows
                         if row.get("minimum_common_support_fraction_both_annuli") == floor]
            for floor in (0.05, 0.10, 0.20)
        },
        "limitations": [
            "The member file is headerless and no official file-specific column order was located; field 8 remains a sensitivity variable.",
            "The published member threshold and field range support, but do not alone prove, the field-8 mapping to w_mem.",
            "CAMIRA candidates are photometric and the weights do not calibrate purity or image-detection completeness.",
            "The ICL quantity is relative and the outer reference annulus may contain cluster light.",
            "Residual sky, extinction, extended PSF wings, source Poisson noise, coadd covariance, centering, and mask selection remain unresolved.",
        ],
        "inputs_sha256": {
            **base.EXPECTED_HASHES,
            str(base.BLOCK_HELPER): base.EXPECTED_BLOCK_HELPER_SHA256,
            str(source_file): base.sha256(source_file),
            str(Path(__file__).resolve()): base.sha256(Path(__file__).resolve()),
        },
        "output_sha256": {
            pair_path.name: base.sha256(pair_path),
            layout_path.name: base.sha256(layout_path),
            support_path.name: base.sha256(support_path),
            plot_path.name: base.sha256(plot_path),
        },
    }
    summary_path = output / "summary.json"
    summary_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2,
                                       sort_keys=True, allow_nan=False) + "\n",
                            encoding="utf-8")
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--bootstrap-replicas", type=int, default=base.N_BOOTSTRAP)
    args = parser.parse_args()
    if args.bootstrap_replicas < 1000:
        parser.error("use at least 1,000 block-bootstrap replicates")
    print(json.dumps(run(args.output, args.bootstrap_replicas), ensure_ascii=False,
                     indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as error:
        print(f"analyze_s23b_icl_satellite_field8_weighted_density: {type(error).__name__}: {error}",
              file=sys.stderr)
        sys.exit(1)
