#!/usr/bin/env python3
"""Sensitivity of the S23B ICL/member-position relation to CAMIRA field 8.

This runner repeats the pinned direct association after replacing unit candidate
counts with sums of the imported member-file field 8. The public S23B paper
describes a CAMIRA membership-weight field w_mem and uses a 1e-3 threshold for
its member-galaxy count, but the downloaded headerless file has no local column
dictionary. Therefore field 8 is treated only as an unverified weight-like
sensitivity, never as a confirmed w_mem mapping or a science measurement.
"""
from __future__ import annotations

import argparse
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

import analyze_s23b_icl_satellite_relation as base


DATA = base.DATA
DEFAULT_OUTPUT = DATA / "s23b_icl_satellite_field8_weight_sensitivity_2026-11-07_01"
INNER_AREA = base.INNER_AREA
OUTER_AREA = base.OUTER_AREA
SEED = 20261108


def field8_annular_sums(member_rows: list[dict], target_manifest: dict) -> tuple[dict, dict]:
    targets = {str(row["id"]): row for row in target_manifest["clusters"]}
    cosmology = FlatLambdaCDM(
        H0=float(target_manifest["h0"]),
        Om0=float(target_manifest["omega_m"]),
    )
    by_target: dict[str, dict] = {
        ident: {"inner_sum": 0.0, "outer_sum": 0.0,
                "inner_n": 0, "outer_n": 0, "all_sum": 0.0,
                "richness": float(row["richness"])}
        for ident, row in targets.items()
    }
    raw_values = []
    source_rows = set()
    for row in member_rows:
        ident = str(row["target_cluster_id"])
        if ident not in by_target:
            raise ValueError(f"member row refers to unpinned target {ident}")
        source_row = str(row["source_row"])
        key = (ident, source_row)
        if key in source_rows:
            raise ValueError(f"duplicate source row in target group: {key}")
        source_rows.add(key)
        weight_like = float(row["member_field8_raw"])
        if not math.isfinite(weight_like) or not 0.0 <= weight_like <= 1.0:
            raise ValueError(f"field 8 outside finite [0,1] range: {ident}/{source_row}")
        raw_values.append(weight_like)
        group = by_target[ident]
        group["all_sum"] += weight_like
        z = float(row["z_cl"])
        kpc_per_arcsec = float(cosmology.kpc_proper_per_arcmin(z).value) / 60.0
        radius_kpc = float(row["cluster_offset_arcsec"]) * kpc_per_arcsec
        in_image = str(row["inside_image_pixel_centers"]).lower() == "true"
        if not in_image:
            continue
        if 50.0 <= radius_kpc < 300.0:
            group["inner_sum"] += weight_like
            group["inner_n"] += 1
        elif 300.0 <= radius_kpc < 800.0:
            group["outer_sum"] += weight_like
            group["outer_n"] += 1

    predictors = {}
    ratios = []
    for ident, group in by_target.items():
        density_ratio = math.log10(
            ((group["inner_sum"] + 0.5) / INNER_AREA)
            / ((group["outer_sum"] + 0.5) / OUTER_AREA)
        )
        predictors[ident] = {
            "inner_field8_sum_in_image": group["inner_sum"],
            "outer_field8_sum_in_image": group["outer_sum"],
            "inner_candidate_count_in_image": group["inner_n"],
            "outer_candidate_count_in_image": group["outer_n"],
            "log10_inner_outer_field8_weighted_density_ratio": density_ratio,
            "field8_sum_all_target_members": group["all_sum"],
            "field8_sum_over_catalog_richness": (
                group["all_sum"] / group["richness"] if group["richness"] else None
            ),
        }
        ratios.append(predictors[ident]["field8_sum_over_catalog_richness"])

    audit = {
        "field_name_in_import": "member_field8_raw",
        "semantic_status": "unverified_weight_like_field; sensitivity only",
        "n_member_rows": len(raw_values),
        "n_unique_source_rows": len(source_rows),
        "value_min_median_max": [
            float(np.min(raw_values)), float(np.median(raw_values)), float(np.max(raw_values))
        ],
        "fraction_in_closed_unit_interval": float(np.mean(
            (np.asarray(raw_values) >= 0.0) & (np.asarray(raw_values) <= 1.0)
        )),
        "all_member_field8_sum_over_catalog_richness_q10_median_q90": [
            float(np.percentile(ratios, 10)),
            float(np.median(ratios)),
            float(np.percentile(ratios, 90)),
        ],
        "interpretation": (
            "The field is numerically consistent with a probability-like weight and "
            "the S23B paper discusses CAMIRA membership weights, but no field-name "
            "mapping is asserted for this headerless local file."
        ),
        "science_measurement": False,
    }
    return predictors, audit


def replace_predictor(rows_by_band: dict[str, list[dict]]) -> dict[str, list[dict]]:
    output = {}
    for band, rows in rows_by_band.items():
        output[band] = []
        for row in rows:
            copy = dict(row)
            copy["unweighted_log10_inner_outer_candidate_density_ratio"] = (
                copy["log10_inner_outer_candidate_density_ratio"]
            )
            copy["log10_inner_outer_candidate_density_ratio"] = (
                copy["log10_inner_outer_field8_weighted_density_ratio"]
            )
            output[band].append(copy)
    return output


def render_plot(path: Path, fit_cache: dict) -> None:
    figure, axes = plt.subplots(1, 2, figsize=(11.5, 5.0), sharey=True)
    for axis, band in zip(axes, ("i", "i2")):
        rows, X, y, metadata, slope, low, high = fit_cache[band]
        nuisance = X[:, [0] + list(range(2, X.shape[1]))]
        x_resid = X[:, 1] - nuisance @ np.linalg.lstsq(nuisance, X[:, 1], rcond=None)[0]
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
        axis.set_xlabel("Residual log field-8 weighted density ratio")
    axes[0].set_ylabel("Residual relative ICL contrast [nJy arcsec$^{-2}$]")
    figure.suptitle(
        "Field-8 weighted member-candidate density vs relative ICL",
        y=0.995, fontsize=11,
    )
    figure.tight_layout(rect=(0, 0, 1, 0.97))
    figure.savefig(path, dpi=180)
    plt.close(figure)


def run(output: Path, replicas: int) -> dict:
    if output.exists():
        raise FileExistsError(f"refusing to overwrite {output}")
    target, control, stack, spatial, profiles, relative_rows, member_summary = base.verify_inputs()
    member_rows = base.read_tsv(base.MEMBER_POSITIONS)
    predictors, field8_audit = field8_annular_sums(member_rows, target)
    pair_rows = base.make_pair_rows(
        target, control, profiles, relative_rows, member_rows
    )
    for row in pair_rows:
        row.update(predictors[row["target_id"]])
        control_row = next(
            item for item in control["clusters"]
            if str(item["target_cluster_id"]) == row["target_id"]
        )
        row["control_ra_deg"] = float(control_row["ra_deg"])
        row["control_dec_deg"] = float(control_row["dec_deg"])
    valid_rows = [
        row for row in pair_rows
        if row["icl_relative_contrast_njy_arcsec2"] is not None
    ]
    raw_rows_by_band = {
        band: [row for row in valid_rows if row["filter"] == band]
        for band in ("i", "i2")
    }
    for band, expected in (("i", 232), ("i2", 756)):
        if len(raw_rows_by_band[band]) != expected:
            raise ValueError(f"unexpected supported pair count for {band}")
    rows_by_band = replace_predictor(raw_rows_by_band)

    results, fit_cache = base.estimate(rows_by_band, spatial, replicas)
    deletion_rows, deletion_summary = base.leave_one_block_out(
        rows_by_band, spatial, replicas
    )
    support_rows = base.support_floor_sensitivity(
        rows_by_band, spatial, replicas
    )
    output.mkdir(parents=True, exist_ok=False)
    pair_table = output / "paired_field8_weighted_relation.tsv"
    slope_table = output / "connected_block_field8_weighted_sensitivity.tsv"
    deletion_table = output / "primary_field8_weighted_leave_one_block_out.tsv"
    support_table = output / "field8_weighted_support_floor_sensitivity.tsv"
    figure_path = output / "adjusted_field8_weighted_relation.png"

    base.write_tsv(pair_table, valid_rows)
    base.write_tsv(slope_table, results)
    base.write_tsv(deletion_table, deletion_rows)
    base.write_tsv(support_table, support_rows)
    render_plot(figure_path, fit_cache)

    primary = {
        band: next(row for row in results
                   if row["filter"] == band
                   and row["block_width_deg"] == 5.0
                   and row["grid_phase_fraction"] == 0.0)
        for band in ("i", "i2")
    }
    summary = {
        "schema": "astrafek-s23b-icl-satellite-field8-weight-sensitivity-v1",
        "status": "complete_exploratory_field8_weight_sensitivity",
        "science_measurement": False,
        "estimator": {
            "icl_outcome": "common-support target-minus-control relative ICL contrast at 50-300 minus 300-800 proper kpc",
            "candidate_predictor": "log10 ratio of inner/outer surface density from sums of raw member file field 8; +0.5 field-unit continuity correction",
            "field8_semantics": "unverified; this is a sensitivity only and does not assign a column name",
            "regression_adjustments": [
                "log10 CAMIRA richness", "cluster photometric redshift",
                "target-minus-control SFD E(B-V) mismatch",
                "minimum common-pixel support fraction",
            ],
            "uncertainty": "5,000 connected RA-Dec block bootstrap resamples per filter over six 2.5/5/10 degree layouts and two phases",
            "primary_layout": "5 degree, phase 0",
        },
        "field8_raw_audit": field8_audit,
        "sample": {
            "all_filter_matched_pairs": len(pair_rows),
            "supported_pairs_by_filter": {
                band: len(rows_by_band[band]) for band in ("i", "i2")
            },
            "total_member_positions": member_summary["catalog_join"]["target_member_rows"],
            "network_used": False,
            "science_image_pixels_read": False,
        },
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
            }
            for band in ("i", "i2")
        },
        "leave_one_block_out": deletion_summary,
        "support_floor_sensitivity_primary_layout": {
            band: [
                row for row in support_rows if row["filter"] == band
                and row["block_width_deg"] == 5.0
                and row["grid_phase_fraction"] == 0.0
            ] for band in ("i", "i2")
        },
        "limitations": [
            "The headerless member-file field 8 is not formally mapped to w_mem here; results are a raw-field sensitivity.",
            "The S23B outer 300-800 kpc annulus may contain cluster light.",
            "Common-pixel support is sparse and selected; the linear support covariate does not correct selection bias.",
            "Residual sky, extinction, extended PSF wings, complete coadd covariance, centering, and member completeness remain open.",
        ],
        "source_context": {
            "S23B_paper": "https://arxiv.org/html/2512.13954",
            "repository_README": "https://raw.githubusercontent.com/oguri/cluster_catalogs/main/hsc_s23b_camira/README.md",
        },
        "inputs_sha256": {
            **base.EXPECTED_HASHES,
            str(base.BLOCK_HELPER): base.EXPECTED_BLOCK_HELPER_SHA256,
            str(Path(__file__).resolve()): base.sha256(Path(__file__).resolve()),
        },
        "output_sha256": {
            pair_table.name: base.sha256(pair_table),
            slope_table.name: base.sha256(slope_table),
            deletion_table.name: base.sha256(deletion_table),
            support_table.name: base.sha256(support_table),
            figure_path.name: base.sha256(figure_path),
        },
        "analysis_script_sha256": base.sha256(Path(__file__).resolve()),
        "bootstrap_replicas_per_layout": replicas,
    }
    (output / "summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--bootstrap-replicas", type=int, default=base.N_BOOTSTRAP)
    args = parser.parse_args()
    if args.bootstrap_replicas < 1000:
        parser.error("use at least 1,000 connected-block resamples")
    summary = run(args.output, args.bootstrap_replicas)
    print(json.dumps({
        "status": summary["status"],
        "output": str(args.output),
        "science_measurement": summary["science_measurement"],
        "field8_raw_audit": summary["field8_raw_audit"],
        "primary_layout_results": summary["primary_layout_results"],
        "spatial_partition_sensitivity": summary["spatial_partition_sensitivity"],
        "support_floor_sensitivity_primary_layout":
            summary["support_floor_sensitivity_primary_layout"],
    }, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as error:
        print(
            f"analyze_s23b_icl_satellite_field8_weight_sensitivity: "
            f"{type(error).__name__}: {error}",
            file=sys.stderr,
        )
        sys.exit(1)
