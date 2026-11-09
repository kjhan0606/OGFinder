#!/usr/bin/env python3
"""Run a common-block restricted wild-bootstrap audit of the S23B relation.

The calculation reuses saved pair-level measurements. It performs no HSC
request, data download, or image-pixel read and leaves the science flag false.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from analyze_s23b_icl_spatial_blocks import connected_pair_blocks


DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
SOURCE_DIR = DATA / "s23b_icl_satellite_relation_2026-11-07_03"
SOURCE_SUMMARY = SOURCE_DIR / "summary.json"
PAIR_TABLE = SOURCE_DIR / "paired_icl_satellite_profiles.tsv"
SPATIAL_SUMMARY = DATA / "s23b_icl_spatial_blocks_2026-11-07_02/summary.json"
DEFAULT_OUTPUT = DATA / "s23b_icl_satellite_relation_wild_2026-11-07_01"
EXPECTED_PAIR_TABLE_SHA256 = "98ed9e0022bfa68d5921ead62dce55be7df1cfa169e6a777d66d1e9760f43994"
EXPECTED_SPATIAL_SUMMARY_SHA256 = "423bb5e3ad52f3b349dedc74b00a2c05ebb4589e20402abb68b4e99d17199afc"
WEBB_WEIGHTS = np.asarray([
    -math.sqrt(1.5), -1.0, -math.sqrt(0.5),
    math.sqrt(0.5), 1.0, math.sqrt(1.5),
], dtype=float)
FILTERS = ("i", "i2")
SUPPORT_FLOORS = (None, 0.05, 0.10, 0.20)
REPLICAS = 10_000
SEED = 20261107
MIN_BLOCKS = 6


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t",
                                extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def verify_inputs() -> tuple[list[dict[str, str]], dict, dict, dict[str, str]]:
    for path in (SOURCE_SUMMARY, PAIR_TABLE, SPATIAL_SUMMARY):
        if not path.is_file():
            raise FileNotFoundError(path)
    summary = json.loads(SOURCE_SUMMARY.read_text(encoding="utf-8"))
    spatial = json.loads(SPATIAL_SUMMARY.read_text(encoding="utf-8"))
    if (summary.get("schema") != "astrafek-s23b-icl-satellite-relation-v1"
            or summary.get("status") != "complete_exploratory_direct_member_position_relation"
            or summary.get("science_measurement") is not False
            or summary.get("sample", {}).get("network_used") is not False
            or summary.get("sample", {}).get("science_image_pixels_read") is not False
            or summary.get("output_sha256", {}).get(PAIR_TABLE.name) != EXPECTED_PAIR_TABLE_SHA256
            or sha256(PAIR_TABLE) != EXPECTED_PAIR_TABLE_SHA256):
        raise ValueError("paired relation input does not match its saved exploratory summary")
    if (sha256(SPATIAL_SUMMARY) != EXPECTED_SPATIAL_SUMMARY_SHA256
            or spatial.get("science_measurement") is not False
            or not math.isclose(float(spatial.get("ra_grid_origin_deg", -1)),
                                84.056247280625, abs_tol=1e-9)):
        raise ValueError("connected-block spatial grid does not match the audited input")
    rows = read_tsv(PAIR_TABLE)
    if len(rows) != 988 or any(str(row.get("science_measurement", "")).lower() != "false"
                               for row in rows):
        raise ValueError("expected 988 exploratory pair/filter records")
    by_pair = {(row["target_id"], row["filter"]) for row in rows}
    if len(by_pair) != len(rows) or {row["filter"] for row in rows} != set(FILTERS):
        raise ValueError("duplicate pair/filter rows or unexpected filter set")
    inputs = {
        str(SOURCE_SUMMARY): sha256(SOURCE_SUMMARY),
        str(PAIR_TABLE): sha256(PAIR_TABLE),
        str(SPATIAL_SUMMARY): sha256(SPATIAL_SUMMARY),
        str(Path(__file__).resolve().with_name("analyze_s23b_icl_spatial_blocks.py")):
            sha256(Path(__file__).resolve().with_name("analyze_s23b_icl_spatial_blocks.py")),
    }
    return rows, summary, spatial, inputs


def design_matrices(rows: list[dict[str, str]]) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[str]]:
    """Return full and restricted designs matching the published S23B fit."""
    n = len(rows)
    x = np.asarray([float(row["log10_inner_outer_candidate_density_ratio"])
                    for row in rows])
    y = np.asarray([float(row["icl_relative_contrast_njy_arcsec2"])
                    for row in rows])
    richness = np.asarray([float(row["log10_richness"]) for row in rows])
    redshift = np.asarray([float(row["z_cl"]) for row in rows])
    delta_ebv = np.asarray([float(row["delta_ebv_sfd_target_minus_control"])
                            for row in rows])
    common_support = np.minimum(
        np.asarray([float(row["inner_common_support_fraction"]) for row in rows]),
        np.asarray([float(row["outer_common_support_fraction"]) for row in rows]),
    )
    member_footprint = np.minimum(
        np.asarray([float(row["member_position_footprint_fraction_50_300"]) for row in rows]),
        np.asarray([float(row["member_position_footprint_fraction_300_800"]) for row in rows]),
    )

    nuisance_cols = [np.ones(n)]
    names = ["intercept"]
    covariates = [
        ("log10_richness", richness - richness.mean()),
        ("redshift_per_0p1", (redshift - redshift.mean()) / 0.1),
        ("delta_ebv_per_0p01", (delta_ebv - delta_ebv.mean()) / 0.01),
        ("minimum_common_support_fraction", common_support - common_support.mean()),
        ("minimum_member_position_footprint_fraction", member_footprint - member_footprint.mean()),
    ]
    for name, values in covariates:
        if float(np.ptp(values)) > 1e-10:
            nuisance_cols.append(values)
            names.append(name)
    nuisance = np.column_stack(nuisance_cols)
    x_centered = x - x.mean()
    full = np.column_stack((nuisance, x_centered))
    names.append("satellite_density_contrast")
    if (np.linalg.matrix_rank(nuisance) != nuisance.shape[1]
            or np.linalg.matrix_rank(full) != full.shape[1]):
        raise ValueError(f"rank-deficient model: {names}")
    return full, nuisance, y, names


def cluster_t(beta: np.ndarray, residual: np.ndarray, design: np.ndarray,
              inverse: np.ndarray, labels: np.ndarray) -> float:
    groups = np.unique(labels)
    n, k = design.shape
    g = len(groups)
    if n <= k or g < 2:
        return math.nan
    scores = np.vstack([
        design[labels == group].T @ residual[labels == group]
        for group in groups
    ])
    meat = scores.T @ scores
    correction = (g / (g - 1)) * ((n - 1) / (n - k))
    variance = correction * (inverse @ meat @ inverse)[-1, -1]
    if not math.isfinite(float(variance)) or variance <= 0:
        return math.nan
    return float(beta[-1] / math.sqrt(float(variance)))


def prepare_group(rows: list[dict[str, str]], labels: np.ndarray,
                  global_block_labels: np.ndarray, minimum_blocks: int) -> tuple[dict | None, dict]:
    full, nuisance, y, names = design_matrices(rows)
    if len(labels) != len(rows):
        raise ValueError("group labels do not align with model rows")
    unique = np.unique(labels)
    n, k = full.shape
    base = {"n_pairs": n, "n_spatial_blocks": len(unique),
            "n_model_columns": k, "model_columns": ",".join(names)}
    if n <= k + 2:
        return None, {**base, "status": "insufficient_pairs"}
    if len(unique) < minimum_blocks:
        return None, {**base, "status": "insufficient_spatial_blocks"}
    full_bread = np.linalg.inv(full.T @ full)
    nuisance_bread = np.linalg.inv(nuisance.T @ nuisance)
    null_beta = nuisance_bread @ nuisance.T @ y
    null_fit = nuisance @ null_beta
    residual_null = y - null_fit
    leverage = np.sum((nuisance @ nuisance_bread) * nuisance, axis=1)
    residual_hc2 = residual_null / np.sqrt(np.maximum(1.0 - leverage, 0.05))
    beta = full_bread @ full.T @ y
    residual = y - full @ beta
    observed_t = cluster_t(beta, residual, full, full_bread, labels)
    if not math.isfinite(observed_t):
        return None, {**base, "status": "invalid_observed_cluster_variance"}
    x = full[:, -1]
    x_resid = x - nuisance @ (nuisance_bread @ nuisance.T @ x)
    denom = float(np.linalg.norm(x_resid) * np.linalg.norm(residual_null))
    partial_r = float(x_resid @ residual_null / denom) if denom > 0 else math.nan
    return {
        "rows": rows, "labels": labels, "full": full, "nuisance": nuisance,
        "y": y, "full_bread": full_bread, "null_fit": null_fit,
        "residual_hc2": residual_hc2, "beta": float(beta[-1]),
        "observed_t": observed_t, "partial_r": partial_r,
        "global_block_labels": global_block_labels,
        "model_column_names": names,
    }, {**base, "status": "eligible", "slope": float(beta[-1]),
        "t_cluster_cr1": observed_t, "partial_r": partial_r}


def bootstrap_t(group: dict, weights_by_block: dict[str, float]) -> float:
    labels = group["labels"]
    weights = np.asarray([weights_by_block[str(label)] for label in labels])
    y_star = group["null_fit"] + group["residual_hc2"] * weights
    design, bread = group["full"], group["full_bread"]
    beta = bread @ design.T @ y_star
    residual = y_star - design @ beta
    return cluster_t(beta, residual, design, bread, labels)


def holm_adjust(p_values: np.ndarray) -> np.ndarray:
    order = np.argsort(p_values)
    adjusted = np.empty(len(p_values), dtype=float)
    running = 0.0
    total = len(p_values)
    for rank, index in enumerate(order):
        running = max(running, (total - rank) * float(p_values[index]))
        adjusted[index] = min(1.0, running)
    return adjusted


def run(output: Path, replicas: int, seed: int, minimum_blocks: int) -> dict:
    if output.exists():
        raise FileExistsError(f"refusing to overwrite {output}")
    all_rows, source_summary, spatial, input_hashes = verify_inputs()
    origin = float(spatial["ra_grid_origin_deg"])
    control_by_target = {
        row["target_id"]: {"ra_deg": float(row["control_ra_deg"]),
                            "dec_deg": float(row["control_dec_deg"])}
        for row in all_rows
    }
    # Construct one partition from the union of both filters. Every test and
    # support-floor subset then uses the same block multiplier for each sky area.
    global_labels, global_blocks = connected_pair_blocks(
        all_rows, control_by_target, 5.0, 0.0, origin
    )
    global_blocks = np.asarray(global_labels, dtype=str)
    groups: dict[tuple[str, float | None], dict] = {}
    record_by_key: dict[tuple[str, float | None], dict] = {}
    for band in FILTERS:
        for floor in SUPPORT_FLOORS:
            selected_indices = [
                index for index, row in enumerate(all_rows)
                if row["filter"] == band
                and (floor is None or min(float(row["inner_common_support_fraction"]),
                                          float(row["outer_common_support_fraction"])) >= floor)
            ]
            selected_rows = [all_rows[index] for index in selected_indices]
            selected_labels = global_blocks[selected_indices]
            key = (band, floor)
            group, record = prepare_group(selected_rows, selected_labels,
                                          global_blocks, minimum_blocks)
            record_by_key[key] = {
                "filter": band,
                "common_support_floor_both_annuli": floor,
                "support_definition": "any positive common support" if floor is None else f">={floor:.0%} in both annuli",
                "status": record["status"],
                "n_pairs": record["n_pairs"],
                "n_spatial_blocks": record["n_spatial_blocks"],
                "n_model_columns": record["n_model_columns"],
                "model_columns": record["model_columns"],
                "slope_njy_arcsec2_per_dex": record.get("slope"),
                "t_cluster_cr1": record.get("t_cluster_cr1"),
                "partial_r": record.get("partial_r"),
                "p_wild_cluster_raw": None,
                "p_holm_across_eligible_groups": None,
                "p_maxT_single_step_complete_null": None,
                "replicas": replicas,
                "science_measurement": False,
            }
            if group is not None:
                groups[key] = group
    if len(groups) < 2:
        raise ValueError("fewer than two filter/support groups meet the inference floor")

    global_unique = np.unique(global_blocks)
    rng = np.random.default_rng(seed)
    t_draws = {key: np.empty(replicas, dtype=float) for key in groups}
    max_t = np.empty(replicas, dtype=float)
    for replica in range(replicas):
        draw = rng.choice(WEBB_WEIGHTS, size=len(global_unique), replace=True)
        weights = dict(zip((str(value) for value in global_unique), draw.tolist()))
        maximum = 0.0
        for key, group in groups.items():
            value = bootstrap_t(group, weights)
            if not math.isfinite(value):
                raise ValueError(f"non-finite wild-bootstrap draw for {key} at {replica}")
            t_draws[key][replica] = value
            maximum = max(maximum, abs(value))
        max_t[replica] = maximum

    ordered_keys = list(groups)
    raw_p = np.asarray([
        (1 + int(np.count_nonzero(np.abs(t_draws[key]) >= abs(groups[key]["observed_t"]))))
        / (replicas + 1) for key in ordered_keys
    ])
    max_t_p = np.asarray([
        (1 + int(np.count_nonzero(max_t >= abs(groups[key]["observed_t"]))))
        / (replicas + 1) for key in ordered_keys
    ])
    holm_p = holm_adjust(raw_p)
    eligible_index = {key: index for index, key in enumerate(ordered_keys)}
    for key, group in groups.items():
        index = eligible_index[key]
        record = record_by_key[key]
        record.update({
            "status": "exploratory_wild_bootstrap_inference",
            "n_pairs": len(group["rows"]),
            "n_spatial_blocks": int(np.unique(group["labels"]).size),
            "n_model_columns": int(group["full"].shape[1]),
            "model_columns": ",".join(group["model_column_names"]),
            "slope_njy_arcsec2_per_dex": group["beta"],
            "t_cluster_cr1": group["observed_t"],
            "partial_r": group["partial_r"],
            "p_wild_cluster_raw": float(raw_p[index]),
            "p_holm_across_eligible_groups": float(holm_p[index]),
            "p_maxT_single_step_complete_null": float(max_t_p[index]),
        })

    records = [record_by_key[(band, floor)]
               for band in FILTERS for floor in SUPPORT_FLOORS]

    output.mkdir(parents=True, exist_ok=False)
    result_path = output / "wild_association.tsv"
    fields = ["filter", "common_support_floor_both_annuli", "support_definition",
              "status", "n_pairs", "n_spatial_blocks", "n_model_columns",
              "model_columns", "slope_njy_arcsec2_per_dex", "t_cluster_cr1",
              "partial_r", "p_wild_cluster_raw", "p_holm_across_eligible_groups",
              "p_maxT_single_step_complete_null", "replicas", "science_measurement"]
    write_tsv(result_path, records, fields)
    null_path = output / "maxT_null.tsv"
    write_tsv(null_path, [{"replica": i + 1, "max_abs_t": float(value)}
                          for i, value in enumerate(max_t)], ["replica", "max_abs_t"])
    summary = {
        "schema": "astrafek-s23b-icl-satellite-wild-association-v1",
        "status": "complete_exploratory",
        "science_measurement": False,
        "network_used": False,
        "image_pixels_read": False,
        "new_data_downloaded": False,
        "method": {
            "model": "relative ICL contrast ~ log inner/outer CAMIRA candidate-density ratio + log richness + redshift + target-control delta E(B-V) + common support + member-footprint term when nonconstant",
            "null": "restricted nuisance model without satellite-density term; HC2-adjusted residuals",
            "bootstrap": "Webb six-point restricted wild-cluster bootstrap-t; one shared set of 5-degree connected target/control block weights across filters and support floors",
            "multiplicity": "Holm-adjusted marginal bootstrap p-values and single-step max absolute t under the complete null",
            "family": "2 filters x 4 prespecified support definitions; only groups with enough pairs and at least the configured number of blocks enter the p-value family",
            "replicas": replicas,
            "seed": seed,
            "minimum_blocks": minimum_blocks,
            "n_common_spatial_blocks": int(len(global_unique)),
        },
        "n_expected_groups": len(FILTERS) * len(SUPPORT_FLOORS),
        "n_eligible_groups": len(groups),
        "n_raw_p_below_0p05": int(np.count_nonzero(raw_p < 0.05)),
        "n_holm_p_below_0p05": int(np.count_nonzero(holm_p < 0.05)),
        "n_maxT_p_below_0p05": int(np.count_nonzero(max_t_p < 0.05)),
        "inputs_sha256": input_hashes,
        "output_sha256": {result_path.name: sha256(result_path),
                          null_path.name: sha256(null_path)},
        "analysis_script_sha256": sha256(Path(__file__).resolve()),
        "limitations": [
            "The single-step max-T result is calibrated under the complete null; subset pivotality is unverified, so strong familywise control is not claimed.",
            "The connected 5-degree blocks are an approximate spatial partition; they do not represent official HSC field polygons or a calibrated survey covariance.",
            "This calculation reuses a coverage-selected, photometric CAMIRA member sample and relative ICL outcome; it does not correct membership completeness or residual sky, extinction, PSF wings, masks, or coadd covariance.",
            "The raw field-8 weighting is not used; candidate counts are unweighted, matching the pinned primary model.",
            "science_measurement remains false.",
        ],
    }
    summary_path = output / "summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True,
                                       ensure_ascii=False, allow_nan=False) + "\n",
                            encoding="utf-8")
    (output / "assessment.md").write_text(render_assessment(summary, records, output),
                                          encoding="utf-8")
    return summary


def render_assessment(summary: dict, records: list[dict], output: Path) -> str:
    lines = [
        "# S23B ICL–satellite relation: restricted wild-cluster inference",
        "",
        "## Scope and method",
        "",
        "This local CPU analysis reuses the saved 988 supported target/filter records. It made no network request, downloaded no data, and read no image pixels. It repeats the pinned adjusted model with unweighted CAMIRA candidate counts, then tests the satellite-density coefficient under a restricted Webb six-point wild-cluster bootstrap. The same 5-degree connected target/control block weights are shared across both filters and all support floors.",
        "",
        "The eight planned groups are HSC `i` and `i2` at any positive common support, and at 5%, 10%, and 20% common support in both annuli. The marginal bootstrap p-values are Holm-adjusted; a single-step max-T p-value is also reported under the complete null. The result retains `science_measurement=false`.",
        "",
        "## Results",
        "",
        "| Filter | Common support | Pairs / blocks | Slope | CR1 t | Raw wild p | Holm p | max-T p |",
        "|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    by_key = {(row.get("filter"), row.get("common_support_floor_both_annuli")): row
              for row in records}
    floor_list = (None, 0.05, 0.10, 0.20)
    for band in FILTERS:
        for floor in floor_list:
            row = by_key[(band, floor)]
            if row.get("status") == "exploratory_wild_bootstrap_inference":
                text = "any positive" if floor is None else f">={floor:.0%}"
                lines.append(
                    f"| {band} | {text} | {row['n_pairs']} / {row['n_spatial_blocks']} | "
                    f"{row['slope_njy_arcsec2_per_dex']:.2f} | {row['t_cluster_cr1']:.2f} | "
                    f"{row['p_wild_cluster_raw']:.4f} | {row['p_holm_across_eligible_groups']:.4f} | "
                    f"{row['p_maxT_single_step_complete_null']:.4f} |"
                )
            else:
                lines.append(f"| {band} | {row['support_definition']} | {row['n_pairs']} / {row['n_spatial_blocks']} | — | — | — | — | — |")
    lines += [
        "",
        f"Among {summary['n_eligible_groups']} eligible groups, {summary['n_raw_p_below_0p05']} raw p-values are below 0.05, {summary['n_holm_p_below_0p05']} remain below 0.05 after Holm correction, and {summary['n_maxT_p_below_0p05']} are below 0.05 under the complete-null max-T diagnostic.",
        "",
        "## Interpretation and limits",
        "",
        "This is a stricter null-model check than the percentile block-bootstrap intervals. A small raw p-value is not treated as a detection unless it also survives the stated family adjustment. The max-T calculation does not establish strong familywise control because subset pivotality across support-selected samples has not been shown. The support cuts are sensitivity conditions on the same images, not independent replication.",
        "",
        "The current result still uses a coverage-selected, photometric candidate sample and a relative ICL contrast. Sky transfer, extinction, PSF wings, masks, membership completeness, and full coadd covariance remain outside this model. Keep `science_measurement=false`.",
        "",
        "## Provenance",
        "",
        f"- Input pair table SHA-256: `{summary['inputs_sha256'][str(PAIR_TABLE)]}`",
        f"- Spatial-grid summary SHA-256: `{summary['inputs_sha256'][str(SPATIAL_SUMMARY)]}`",
        f"- Analysis script SHA-256: `{summary['analysis_script_sha256']}`",
        f"- Replicas: {summary['method']['replicas']:,}; seed: {summary['method']['seed']}",
        f"- Output summary SHA-256: `{sha256(output / 'summary.json')}`",
        "",
        "The paired association table is in `wild_association.tsv`; the 10,000-draw complete-null max-T reference distribution is in `maxT_null.tsv`.",
    ]
    return "\n".join(lines) + "\n"


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--replicas", type=int, default=REPLICAS)
    parser.add_argument("--seed", type=int, default=SEED)
    parser.add_argument("--min-blocks", type=int, default=MIN_BLOCKS)
    args = parser.parse_args()
    if args.replicas < 1000:
        parser.error("at least 1,000 restricted wild-bootstrap draws are required")
    if args.min_blocks < 2:
        parser.error("minimum block count must be at least 2")
    summary = run(args.output.resolve(), args.replicas, args.seed, args.min_blocks)
    print(json.dumps({key: summary[key] for key in (
        "status", "science_measurement", "n_expected_groups", "n_eligible_groups",
        "n_raw_p_below_0p05", "n_holm_p_below_0p05", "n_maxT_p_below_0p05",
        "output_sha256")}, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
