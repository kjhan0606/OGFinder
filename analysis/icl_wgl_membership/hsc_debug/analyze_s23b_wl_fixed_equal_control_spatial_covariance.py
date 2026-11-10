#!/usr/bin/env python3
"""Compare matched-random weighting rules with the same fixed-sample covariance.

This offline analysis reuses the approved S23B HSC S19A shape-row export and
the pinned Astrafek/HSC Y3 products. It fixes the 52 lens sets supported in all
seven radial bins, compares source-pair-pooled controls with an equal mean of
positive individual control-center profiles, and evaluates both with paired
delete-one endpoint-cell covariance. Equal-area cells are exploratory spatial
blocks, not validated contiguous HSC survey regions. The science flag remains
false.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
import sys

import numpy as np
from scipy.stats import chi2 as chi2_distribution

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import analyze_s23b_wl_fixed_common_support_spatial_covariance as base
from fetch_ancillary_catalogues import sha
from ogfkit.wlcluster import Cosmology
from s19a_shape_stack import load_hsc_y3_nz_sacc
from submit_s23b_wl_common_shape_rows import DEFAULT_PLAN_DIR, DEFAULT_RESULTS


DATA = base.DATA
DEFAULT_OUTPUT = DATA / "wl_s23b_equal_random_center_spatial_covariance_2026-11-11_01"
DEFAULT_POOLED_SUMMARY = DATA / "wl_s23b_fixed_common_support_spatial_covariance_2026-11-11_02/summary.json"
EXPECTED_POOLED_SUMMARY_SHA256 = "acdd4da69e4d7fe3dbb52d8584b3996bd836a41196da59aecfe3fcdf3e1fc6e9"
VARIANTS = ("pooled_controls", "equal_control_centers", "equal_minus_pooled")


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with Path(path).open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def finite_profile(result) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    nbin = len(base.RADIAL_EDGES_MPC) - 1
    values = np.full((nbin, 2), np.nan, dtype=float)
    pairs = np.zeros(nbin, dtype=np.int64)
    centers = np.zeros(nbin, dtype=np.int64)
    if result is None:
        return values, pairs, centers
    for b, row in enumerate(result["profile"]):
        pairs[b] = int(row["n_pairs"])
        centers[b] = int(row["n_lenses"])
        tangential = row["delta_sigma_weak_msun_mpc2"]
        cross = row["cross_delta_sigma_weak_msun_mpc2"]
        if tangential is not None and cross is not None:
            values[b] = (float(tangential), float(cross))
    return values, pairs, centers


def equal_center_controls(control_pairs, control_cache, indices_by_center,
                          control_ids: list[str], cosmology: Cosmology,
                          keep_mask: np.ndarray | None = None) -> dict:
    """Average individually calibrated control centers within each radial bin."""
    nbin = len(base.RADIAL_EDGES_MPC) - 1
    values: list[list[tuple[float, float]]] = [[] for _ in range(nbin)]
    pair_counts = np.zeros(nbin, dtype=np.int64)
    for center_id in control_ids:
        indices = indices_by_center.get(center_id)
        if indices is None or not len(indices):
            continue
        if keep_mask is not None:
            indices = indices[keep_mask[indices]]
        if not len(indices):
            continue
        result = base.cached_stack(control_pairs, control_cache, indices, cosmology)
        center_values, center_pairs, center_counts = finite_profile(result)
        for b in range(nbin):
            if center_counts[b] > 0 and np.all(np.isfinite(center_values[b])):
                values[b].append((float(center_values[b, 0]), float(center_values[b, 1])))
                pair_counts[b] += center_pairs[b]
    means = np.full((nbin, 2), np.nan, dtype=float)
    n_positive = np.zeros(nbin, dtype=np.int64)
    for b, rows in enumerate(values):
        n_positive[b] = len(rows)
        if rows:
            means[b] = np.mean(np.asarray(rows, dtype=float), axis=0)
    return {"values": means, "n_centers": n_positive, "n_pairs": pair_counts}


def subtracted_profile(target_result, pooled_result, equal_result) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    target, target_pairs, _ = finite_profile(target_result)
    pooled, pooled_pairs, pooled_centers = finite_profile(pooled_result)
    equal, equal_centers = equal_result["values"], equal_result["n_centers"]
    pooled_delta = target - pooled
    equal_delta = target - equal
    nbin = len(base.RADIAL_EDGES_MPC) - 1
    for b in range(nbin):
        if target_pairs[b] == 0 or not np.all(np.isfinite(target[b])):
            pooled_delta[b] = np.nan
            equal_delta[b] = np.nan
        if pooled_pairs[b] == 0 or not np.all(np.isfinite(pooled[b])):
            pooled_delta[b] = np.nan
        if not np.all(np.isfinite(equal[b])):
            equal_delta[b] = np.nan
    if np.any(pooled_centers != equal_centers):
        raise ValueError("pooled and equal-center estimators disagree on positive-control support")
    return pooled_delta, equal_delta, equal_centers


def flatten(profile: np.ndarray) -> np.ndarray:
    return np.asarray(profile, dtype=float).reshape(-1)


def covariance_or_none(replicates: list[np.ndarray], minimum_sets: int,
                       nblocks: int) -> tuple[np.ndarray | None, dict]:
    return base.covariance_result(replicates, minimum_sets, nblocks)


def zero_vector_diagnostic(values: np.ndarray, covariance: np.ndarray | None,
                           component_indices: np.ndarray, label: str) -> dict:
    """Return a conditional quadratic-form check for one covariance subvector."""
    if covariance is None:
        return {"status": "not_estimable_covariance_missing", "label": label}
    vector = np.asarray(values, dtype=float)
    indices = np.asarray(component_indices, dtype=np.int64)
    subcovariance = np.asarray(covariance, dtype=float)[np.ix_(indices, indices)]
    if (vector.shape != (len(indices),) or not np.all(np.isfinite(vector))
            or not np.all(np.isfinite(subcovariance))):
        return {"status": "not_estimable_nonfinite_vector_or_covariance", "label": label}
    eigenvalues, eigenvectors = np.linalg.eigh((subcovariance + subcovariance.T) * 0.5)
    tolerance = max(float(np.max(np.abs(eigenvalues))), 1.0) * 1e-10
    keep = eigenvalues > tolerance
    rank = int(np.count_nonzero(keep))
    if rank == 0 or float(np.min(eigenvalues)) < -tolerance:
        return {"status": "not_estimable_rank_zero_or_non_psd", "label": label,
                "covariance_rank": rank}
    projected = eigenvectors[:, keep].T @ vector
    statistic = float(np.sum(projected * projected / eigenvalues[keep]))
    return {
        "status": "conditional_nominal_zero_vector_diagnostic",
        "label": label,
        "chi_square": statistic,
        "degrees_of_freedom": rank,
        "nominal_p_value": float(chi2_distribution.sf(statistic, rank)),
        "covariance_rank": rank,
        "warning": "conditional on the approximate equal-area endpoint jackknife; not a confirmatory test",
    }


def run(args) -> dict:
    output = args.output.resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite existing output: {output}")
    if sha(args.pooled_summary) != EXPECTED_POOLED_SUMMARY_SHA256:
        raise ValueError("pinned pooled spatial-covariance summary hash differs")
    pooled_summary = base.read_json(args.pooled_summary)
    if (pooled_summary.get("science_measurement") is not False
            or pooled_summary.get("sample", {}).get("n_fixed_lens_sets") != base.EXPECTED_FIXED_LENSES):
        raise ValueError("pinned pooled reference has unexpected science state or sample size")

    lens_ids, fixed_selection = base.fixed_sample_ids(args.support_dir)
    lens_pairs, all_control_pairs, center_to_lens, row_provenance = base.load_approved_pairs(
        args.plan_dir.resolve(), args.rows.resolve(),
    )
    controls_by_lens: dict[str, list[str]] = {lens_id: [] for lens_id in lens_ids}
    for control_id, matched_lens in center_to_lens.items():
        if matched_lens in lens_ids:
            controls_by_lens[matched_lens].append(control_id)
    for lens_id, centers in controls_by_lens.items():
        if len(centers) != 5:
            raise ValueError(f"fixed lens set {lens_id} does not have exactly five matched controls")
        controls_by_lens[lens_id] = sorted(centers)
    selected_control_ids = {center for centers in controls_by_lens.values() for center in centers}

    lens_pairs = base.select_rows(lens_pairs, lens_ids)
    control_pairs = base.select_rows(all_control_pairs, selected_control_ids)
    del all_control_pairs
    nz_z, nz_bins = load_hsc_y3_nz_sacc(args.nz_sacc.resolve())
    cosmology = Cosmology(H0=70.0, Om=0.3)
    lens_cache = base.cache_geometry_and_efficiency(lens_pairs, nz_z, nz_bins, cosmology)
    control_cache = base.cache_geometry_and_efficiency(control_pairs, nz_z, nz_bins, cosmology)
    source_reuse = base.verify_repeated_source_rows(lens_pairs, control_pairs)
    lens_rows_by_set = base.group_indices(lens_pairs["lens_id"])
    control_rows_by_center = base.group_indices(control_pairs["lens_id"])

    full_by_lens = {}
    per_lens_rows = []
    pooled_set_values, equal_set_values = {}, {}
    for lens_id in sorted(lens_ids):
        target_indices = lens_rows_by_set.get(lens_id, np.asarray([], dtype=np.int64))
        center_ids = controls_by_lens[lens_id]
        center_indices = {
            center_id: control_rows_by_center.get(center_id, np.asarray([], dtype=np.int64))
            for center_id in center_ids
        }
        random_indices = np.concatenate([center_indices[center_id] for center_id in center_ids])
        target_result = base.cached_stack(lens_pairs, lens_cache, target_indices, cosmology)
        pooled_result = base.cached_stack(control_pairs, control_cache, random_indices, cosmology)
        equal_result = equal_center_controls(
            control_pairs, control_cache, center_indices, center_ids, cosmology,
        )
        pooled_delta, equal_delta, n_centers = subtracted_profile(
            target_result, pooled_result, equal_result,
        )
        if (np.any(~np.isfinite(pooled_delta)) or np.any(~np.isfinite(equal_delta))
                or np.any(n_centers < base.MIN_CONTROLS_PER_ANNULUS)):
            raise ValueError(f"full fixed sample lacks four positive centers in all bins for {lens_id}")
        pooled_set_values[lens_id] = pooled_delta
        equal_set_values[lens_id] = equal_delta
        full_by_lens[lens_id] = {
            "target": target_result, "pooled": pooled_result, "equal": equal_result,
            "target_indices": target_indices, "random_indices": random_indices,
            "center_indices": center_indices, "center_ids": center_ids,
        }
        target_values, target_pairs, _ = finite_profile(target_result)
        pooled_values, pooled_pairs, _ = finite_profile(pooled_result)
        for b, (lo, hi) in enumerate(zip(base.RADIAL_EDGES_MPC[:-1], base.RADIAL_EDGES_MPC[1:])):
            per_lens_rows.append({
                "matched_lens_id": lens_id,
                "radial_bin": b,
                "r_min_mpc": float(lo), "r_max_mpc": float(hi),
                "n_lens_pairs": int(target_pairs[b]),
                "n_random_pairs": int(pooled_pairs[b]),
                "n_positive_random_centers": int(n_centers[b]),
                "lens_delta_sigma_msun_mpc2": float(target_values[b, 0]),
                "pooled_random_delta_sigma_msun_mpc2": float(pooled_values[b, 0]),
                "equal_center_random_delta_sigma_msun_mpc2": float(equal_result["values"][b, 0]),
                "pooled_subtracted_delta_sigma_msun_mpc2": float(pooled_delta[b, 0]),
                "equal_center_subtracted_delta_sigma_msun_mpc2": float(equal_delta[b, 0]),
                "paired_estimator_difference_msun_mpc2": float(equal_delta[b, 0] - pooled_delta[b, 0]),
                "lens_cross_msun_mpc2": float(target_values[b, 1]),
                "pooled_random_cross_msun_mpc2": float(pooled_values[b, 1]),
                "equal_center_random_cross_msun_mpc2": float(equal_result["values"][b, 1]),
                "pooled_subtracted_cross_msun_mpc2": float(pooled_delta[b, 1]),
                "equal_center_subtracted_cross_msun_mpc2": float(equal_delta[b, 1]),
                "paired_cross_difference_msun_mpc2": float(equal_delta[b, 1] - pooled_delta[b, 1]),
            })

    full_profiles = {
        "pooled_controls": np.mean(np.stack(list(pooled_set_values.values())), axis=0),
        "equal_control_centers": np.mean(np.stack(list(equal_set_values.values())), axis=0),
    }
    full_profiles["equal_minus_pooled"] = (
        full_profiles["equal_control_centers"] - full_profiles["pooled_controls"]
    )
    # The pooled reference is independent of block width; use its recorded 5-degree profile.
    pooled_ref_width = next(
        item for item in pooled_summary["block_width_results"]
        if float(item["block_width_deg"]) == 5.0
    )
    prior_profile = np.asarray([
        [row["delta_sigma_msun_mpc2"], row["cross_delta_sigma_msun_mpc2"]]
        for row in pooled_ref_width["full_profile"]
    ], dtype=float)
    if not np.allclose(full_profiles["pooled_controls"], prior_profile, rtol=2e-9, atol=1e-3):
        raise ValueError("pooled estimator does not reproduce the pinned 52-set reference profile")

    per_width = []
    deletion_rows, support_rows, replicate_rows, covariance_rows = [], [], [], []
    cross_diagnostic_rows = []
    for width in args.block_widths:
        lens_center_blocks = base.equal_area_block_ids(
            lens_pairs["lens_ra_deg"], lens_pairs["lens_dec_deg"], cell_width_deg=width,
        )
        lens_source_blocks = base.equal_area_block_ids(
            lens_pairs["source_ra_deg"], lens_pairs["source_dec_deg"], cell_width_deg=width,
        )
        control_center_blocks = base.equal_area_block_ids(
            control_pairs["lens_ra_deg"], control_pairs["lens_dec_deg"], cell_width_deg=width,
        )
        control_source_blocks = base.equal_area_block_ids(
            control_pairs["source_ra_deg"], control_pairs["source_dec_deg"], cell_width_deg=width,
        )
        base.verify_repeated_cells(
            np.concatenate((lens_pairs["object_id"], control_pairs["object_id"])),
            np.concatenate((lens_source_blocks, control_source_blocks)),
            f"source at {width:g} degrees",
        )
        lens_active, control_active = lens_cache["active"], control_cache["active"]
        block_ids = np.unique(np.concatenate((
            lens_center_blocks[lens_active], lens_source_blocks[lens_active],
            control_center_blocks[control_active], control_source_blocks[control_active],
        )))
        if len(block_ids) < 3:
            raise ValueError(f"only {len(block_ids)} occupied {width:g}-degree endpoint cells")
        affected_by_lens = {}
        for lens_id in sorted(lens_ids):
            entry = full_by_lens[lens_id]
            target_active = entry["target_indices"][lens_active[entry["target_indices"]]]
            random_active = entry["random_indices"][control_active[entry["random_indices"]]]
            cells = []
            if target_active.size:
                cells.extend((lens_center_blocks[target_active], lens_source_blocks[target_active]))
            if random_active.size:
                cells.extend((control_center_blocks[random_active], control_source_blocks[random_active]))
            affected_by_lens[lens_id] = (
                set(np.unique(np.concatenate(cells)).astype(int).tolist()) if cells else set()
            )

        replicate_sets = {variant: [] for variant in VARIANTS}
        counts_by_block = []
        min_controls_by_block = []
        print(f"width {width:g} deg: {len(block_ids)} occupied endpoint cells", flush=True)
        for block_index, block_id in enumerate(block_ids, start=1):
            keep_lens = (lens_center_blocks != block_id) & (lens_source_blocks != block_id)
            keep_control = (control_center_blocks != block_id) & (control_source_blocks != block_id)
            accepted = {variant: {} for variant in VARIANTS}
            finite_sets = np.zeros(len(base.RADIAL_EDGES_MPC) - 1, dtype=int)
            ge4_sets = np.zeros(len(base.RADIAL_EDGES_MPC) - 1, dtype=int)
            min_controls_seen = []
            deleted_lens = int(np.count_nonzero(~keep_lens))
            deleted_control = int(np.count_nonzero(~keep_control))
            affected_sets = 0

            for lens_id in sorted(lens_ids):
                entry = full_by_lens[lens_id]
                if int(block_id) in affected_by_lens[lens_id]:
                    affected_sets += 1
                    target_indices = entry["target_indices"][keep_lens[entry["target_indices"]]]
                    random_indices = entry["random_indices"][keep_control[entry["random_indices"]]]
                    target_result = base.cached_stack(lens_pairs, lens_cache, target_indices, cosmology)
                    pooled_result = base.cached_stack(control_pairs, control_cache, random_indices, cosmology)
                    equal_result = equal_center_controls(
                        control_pairs, control_cache, entry["center_indices"], entry["center_ids"],
                        cosmology, keep_mask=keep_control,
                    )
                else:
                    target_result, pooled_result, equal_result = (
                        entry["target"], entry["pooled"], entry["equal"],
                    )
                pooled_delta, equal_delta, n_centers = subtracted_profile(
                    target_result, pooled_result, equal_result,
                )
                finite = np.all(np.isfinite(pooled_delta), axis=1) & np.all(np.isfinite(equal_delta), axis=1)
                enough = n_centers >= base.MIN_CONTROLS_PER_ANNULUS
                finite_sets += finite.astype(int)
                ge4_sets += (finite & enough).astype(int)
                min_controls_seen.extend(n_centers[finite & enough].tolist())
                if np.all(finite & enough):
                    accepted["pooled_controls"][lens_id] = pooled_delta
                    accepted["equal_control_centers"][lens_id] = equal_delta
                    accepted["equal_minus_pooled"][lens_id] = equal_delta - pooled_delta

            common_ids = set(accepted["pooled_controls"])
            if any(set(accepted[variant]) != common_ids for variant in VARIANTS):
                raise ValueError("estimator variants do not share the same delete-one lens-set intersection")
            n_common = len(common_ids)
            if n_common >= 3:
                for variant in VARIANTS:
                    values = np.mean(np.stack(list(accepted[variant].values())), axis=0)
                    replicate_sets[variant].append(flatten(values))
                counts_by_block.append(n_common)
                min_controls = min(min_controls_seen) if min_controls_seen else 0
                min_controls_by_block.append(min_controls)
                status = "complete"
            else:
                values = None
                min_controls = min(min_controls_seen) if min_controls_seen else 0
                status = "fewer_than_three_common_lens_sets"

            for b in range(len(base.RADIAL_EDGES_MPC) - 1):
                support_rows.append({
                    "block_width_deg": float(width), "block_id": int(block_id), "radial_bin": b,
                    "r_min_mpc": float(base.RADIAL_EDGES_MPC[b]),
                    "r_max_mpc": float(base.RADIAL_EDGES_MPC[b + 1]),
                    "finite_lens_and_both_random_estimators": int(finite_sets[b]),
                    "at_least_four_positive_random_centers": int(ge4_sets[b]),
                    "common_sets_all_seven_bins": n_common,
                    "minimum_positive_centers_in_common_sets": min_controls,
                })
                for variant in VARIANTS:
                    val = values[b] if values is not None else np.asarray([np.nan, np.nan])
                    replicate_rows.append({
                        "block_width_deg": float(width), "block_id": int(block_id),
                        "variant": variant, "radial_bin": b,
                        "r_min_mpc": float(base.RADIAL_EDGES_MPC[b]),
                        "r_max_mpc": float(base.RADIAL_EDGES_MPC[b + 1]),
                        "status": status, "n_common_lens_sets": n_common,
                        "delta_sigma_msun_mpc2": float(val[0]) if values is not None else None,
                        "cross_delta_sigma_msun_mpc2": float(val[1]) if values is not None else None,
                    })
            deletion_rows.append({
                "block_width_deg": float(width), "block_id": int(block_id),
                "affected_lens_sets": affected_sets,
                "deleted_lens_pair_rows": deleted_lens,
                "deleted_control_pair_rows": deleted_control,
                "common_lens_sets_all_bins": n_common,
                "minimum_positive_centers_in_common_sets": min_controls,
                "status": status,
            })
            if block_index == 1 or block_index % 5 == 0 or block_index == len(block_ids):
                print(f"  completed {block_index}/{len(block_ids)} cells", flush=True)

        variant_results = {}
        full_for_width = []
        for variant in VARIANTS:
            covariance, diagnostics = covariance_or_none(
                replicate_sets[variant],
                min(counts_by_block) if counts_by_block else 0,
                len(block_ids),
            )
            errors = np.sqrt(np.diag(covariance)) if covariance is not None else np.full(base.EXPECTED_DIMENSION, np.nan)
            values = full_profiles[variant]
            diagnostic_label = (
                "paired equal-center minus pooled control cross component"
                if variant == "equal_minus_pooled"
                else f"{variant} cross-shear profile"
            )
            cross_diagnostic = zero_vector_diagnostic(
                values[:, 1], covariance, np.arange(1, base.EXPECTED_DIMENSION, 2),
                diagnostic_label,
            )
            cross_diagnostic_rows.append({
                "block_width_deg": float(width), "variant": variant,
                **cross_diagnostic,
            })
            radial = []
            for b, (lo, hi) in enumerate(zip(base.RADIAL_EDGES_MPC[:-1], base.RADIAL_EDGES_MPC[1:])):
                radial.append({
                    "radial_bin": b, "r_min_mpc": float(lo), "r_max_mpc": float(hi),
                    "n_full_lens_sets": base.EXPECTED_FIXED_LENSES,
                    "delta_sigma_msun_mpc2": float(values[b, 0]),
                    "delta_sigma_error_msun_mpc2": float(errors[2 * b]) if covariance is not None else None,
                    "cross_delta_sigma_msun_mpc2": float(values[b, 1]),
                    "cross_delta_sigma_error_msun_mpc2": float(errors[2 * b + 1]) if covariance is not None else None,
                    "n_delete_blocks": len(block_ids),
                    "minimum_common_sets_after_delete": diagnostics.get("minimum_common_lens_sets_across_deletions"),
                    "covariance_status": diagnostics["status"],
                })
                full_for_width.append({"variant": variant, **radial[-1]})
            variant_results[variant] = {
                "covariance": diagnostics,
                "cross_component_zero_diagnostic": cross_diagnostic,
                "full_profile": radial,
                "covariance_matrix": covariance.tolist() if covariance is not None else None,
            }
            if covariance is not None:
                for i in range(base.EXPECTED_DIMENSION):
                    for j in range(base.EXPECTED_DIMENSION):
                        covariance_rows.append({
                            "block_width_deg": float(width), "variant": variant,
                            "component_i": "tangential" if i % 2 == 0 else "cross",
                            "radial_bin_i": i // 2,
                            "component_j": "tangential" if j % 2 == 0 else "cross",
                            "radial_bin_j": j // 2,
                            "covariance_msun2_mpc4": float(covariance[i, j]),
                        })
        per_width.append({
            "block_width_deg": float(width),
            "n_occupied_endpoint_cells": int(len(block_ids)),
            "active_lens_pair_rows": int(np.count_nonzero(lens_active)),
            "active_control_pair_rows": int(np.count_nonzero(control_active)),
            "minimum_common_sets_after_delete": min(counts_by_block) if counts_by_block else 0,
            "maximum_common_sets_after_delete": max(counts_by_block) if counts_by_block else 0,
            "minimum_positive_centers_after_delete": min(min_controls_by_block) if min_controls_by_block else 0,
            "estimators": variant_results,
        })

    output.mkdir(parents=True, exist_ok=False)
    write_tsv(output / "per_lens_set_profiles.tsv", per_lens_rows, list(per_lens_rows[0]))
    write_tsv(output / "block_deletion_support.tsv", deletion_rows, list(deletion_rows[0]))
    write_tsv(output / "block_support_by_radius.tsv", support_rows, list(support_rows[0]))
    write_tsv(output / "block_replicates.tsv", replicate_rows, list(replicate_rows[0]))
    write_tsv(output / "estimator_profiles.tsv", [
        {"block_width_deg": item["block_width_deg"], **row}
        for item in per_width for variant in VARIANTS
        for row in item["estimators"][variant]["full_profile"]
    ], [
        "block_width_deg", "radial_bin", "r_min_mpc", "r_max_mpc", "n_full_lens_sets",
        "delta_sigma_msun_mpc2", "delta_sigma_error_msun_mpc2",
        "cross_delta_sigma_msun_mpc2", "cross_delta_sigma_error_msun_mpc2",
        "n_delete_blocks", "minimum_common_sets_after_delete", "covariance_status",
    ])
    write_tsv(output / "cross_component_diagnostics.tsv", cross_diagnostic_rows,
              list(cross_diagnostic_rows[0]))
    if covariance_rows:
        write_tsv(output / "covariance_matrices.tsv", covariance_rows, list(covariance_rows[0]))

    summary = {
        "schema": "astrafek-s23b-wl-fixed-equal-control-spatial-covariance-v2",
        "status": "complete_exploratory_fixed_52_random_weighting_covariance",
        "science_measurement": False,
        "network_used": False,
        "cas_query_submitted": False,
        "new_data_downloaded": False,
        "images_read": False,
        "full_photoz_pdfs_included": False,
        "halo_mass_fit_performed": False,
        "scope": "pooled versus equal-control-center random subtraction on fixed 52 lens sets with paired endpoint-cell covariance",
        "sample": {
            **fixed_selection,
            "n_control_centers_selected": len(selected_control_ids),
            "n_lens_source_pair_rows": int(len(lens_pairs["lens_id"])),
            "n_control_source_pair_rows": int(len(control_pairs["lens_id"])),
            "source_reuse_audit": source_reuse,
            "pooled_profile_reproduces_pinned_reference": True,
            "pooled_reference_max_absolute_difference_msun_mpc2": float(np.max(np.abs(full_profiles["pooled_controls"] - prior_profile))),
        },
        "estimator": {
            "source_bins": sorted(base.SOURCE_BINS),
            "radial_edges_proper_mpc": base.RADIAL_EDGES_MPC.tolist(),
            "outer_lens_set_weighting": "equal arithmetic mean over the same fixed 52 lens sets",
            "pooled_controls": "calibrate the union of available matched-control source pairs by source-pair lensing weights, then subtract from the separately calibrated lens profile",
            "equal_control_centers": "calibrate each matched-control center separately, average the positive center profiles equally within each lens set and radial bin, then subtract from the separately calibrated lens profile",
            "support_rule": "at least four positive matched random centers per lens set in every radial bin; the same seven-bin lens-set intersection is used for all estimator variants and deletions",
            "endpoint_rule": "delete a center-source pair when either its center endpoint or source endpoint falls in the deleted equal-area RA/sin(Dec) cell; recalibrate both roles after each deletion",
            "block_widths_deg": args.block_widths,
        },
        "block_width_results": per_width,
        "limitations": [
            "equal-area cells are not validated contiguous HSC regions and may combine disconnected fields",
            "the compared random weighting rules do not address whether either random-aperture population matches target depth, mask, PSF, sky, boost, or member contamination",
            "fiducial n(z) uncertainty, photo-z dilution, reduced shear, magnification, and selection-response uncertainty are not propagated",
            "covariance remains an exploratory block-jackknife sensitivity and cannot establish a detection or halo mass",
            "no direct pixel image read, PSF model convolution, or completeness correction was performed",
        ],
        "provenance": {
            "runner_path": str(Path(__file__).resolve()),
            "runner_sha256": sha(Path(__file__).resolve()),
            "pooled_covariance_summary_sha256": sha(args.pooled_summary),
            "approved_row_plan_sha256": row_provenance["plan_sha256"],
            "approved_row_status_sha256": row_provenance["row_status_sha256"],
            "approved_row_summary_sha256": row_provenance["row_summary_sha256"],
            "approved_source_row_batch_sha256": row_provenance["batch_sha256"],
            "radial_support_summary_sha256": sha(args.support_dir / "summary.json"),
            "per_lens_radial_support_sha256": sha(args.support_dir / "per_lens_radial_support.tsv"),
            "qa_summary_sha256": sha(args.qa_dir / "summary.json"),
            "n_z_sacc_sha256": sha(args.nz_sacc),
            "fixed_lens_ids_sha256": fixed_selection["fixed_lens_ids_sha256"],
        },
        "output_directory": str(output),
    }
    assessment = [
        "# Fixed-52 matched-random weighting with shared spatial covariance (2026-11-11)", "",
        "This offline Astrafek analysis uses only the already approved HSC S19A shape-row CSVs, pinned support/QA artifacts, the existing HSC Y3 SACC n(z), and the prior pooled-profile covariance summary. It performs no CAS query, network request, download, image read, or halo-mass fit. `science_measurement=false`.", "",
        "For each of the same 52 lens sets, the target profile is calibrated separately. The pooled random estimator combines all supported control source pairs by their lensing weights. The equal-center estimator calibrates each random center separately and gives each positive center equal weight in that annulus. Both subtract from the same target profile and use an equal outer mean over the fixed lens sets.", "",
        "Every delete-one replicate removes a pair if either endpoint is in the deleted equal-area RA/sin(Dec) cell, recalibrates target and controls, applies the same at-least-four-control floor in every annulus, and uses the common seven-bin lens-set intersection for all three reported variants. The estimator-difference error is paired across the same spatial deletions; it is a sensitivity diagnostic, not a significance test.", "",
        "| R proper (Mpc) | Pooled ΔΣt | Pooled error 2.5° / 5° | Equal-center ΔΣt | Equal-center error 2.5° / 5° | Equal−pooled | Paired-difference error 2.5° / 5° |", "|---:|---:|---:|---:|---:|---:|---:|",
    ]
    by_width = {float(item["block_width_deg"]): item["estimators"] for item in per_width}
    for b, (lo, hi) in enumerate(zip(base.RADIAL_EDGES_MPC[:-1], base.RADIAL_EDGES_MPC[1:])):
        p25 = by_width[2.5]["pooled_controls"]["full_profile"][b]
        p5 = by_width[5.0]["pooled_controls"]["full_profile"][b]
        e25 = by_width[2.5]["equal_control_centers"]["full_profile"][b]
        e5 = by_width[5.0]["equal_control_centers"]["full_profile"][b]
        d25 = by_width[2.5]["equal_minus_pooled"]["full_profile"][b]
        d5 = by_width[5.0]["equal_minus_pooled"]["full_profile"][b]
        scale = 1e13
        assessment.append(
            f"| {lo:.1f}–{hi:.1f} | {p5['delta_sigma_msun_mpc2']/scale:.4f} | "
            f"{p25['delta_sigma_error_msun_mpc2']/scale:.4f} / {p5['delta_sigma_error_msun_mpc2']/scale:.4f} | "
            f"{e5['delta_sigma_msun_mpc2']/scale:.4f} | "
            f"{e25['delta_sigma_error_msun_mpc2']/scale:.4f} / {e5['delta_sigma_error_msun_mpc2']/scale:.4f} | "
            f"{d5['delta_sigma_msun_mpc2']/scale:+.4f} | "
            f"{d25['delta_sigma_error_msun_mpc2']/scale:.4f} / {d5['delta_sigma_error_msun_mpc2']/scale:.4f} |"
        )
    assessment.extend([
        "",
        "## Cross-component conditional zero-profile diagnostics",
        "",
        "The nominal quadratic-form values below use the cross-component submatrix of each paired endpoint covariance. They condition on the approximate equal-area cells; the nominal p-values are not detection or calibration claims.",
        "",
        "| Cell width | Pooled cross profile χ²/dof (p) | Equal-center cross profile χ²/dof (p) | Equal−pooled cross contrast χ²/dof (p) |",
        "|---:|---:|---:|---:|",
    ])
    for item in per_width:
        estimates = item["estimators"]
        def diagnostic_text(variant):
            diag = estimates[variant]["cross_component_zero_diagnostic"]
            if diag.get("status") != "conditional_nominal_zero_vector_diagnostic":
                return diag.get("status", "not available")
            return f"{diag['chi_square']:.3f}/{diag['degrees_of_freedom']} ({diag['nominal_p_value']:.4f})"
        assessment.append(
            f"| {item['block_width_deg']:g}° | {diagnostic_text('pooled_controls')} | "
            f"{diagnostic_text('equal_control_centers')} | {diagnostic_text('equal_minus_pooled')} |"
        )
    assessment.extend([
        "",
        "The paired-difference covariance reuses each block deletion and common set intersection, so it directly quantifies how the weighting choice moves the profile under this approximate spatial partition. It does not validate the partition as an HSC survey footprint. A numerical PSD/full-rank gate is reported per estimator and cell width; a passed gate remains exploratory.",
        "",
        "The data do not include full photo-z PDFs or image-level target/control matching. No estimator is selected as final from this comparison. Keep `science_measurement=false`; mass fitting and detection claims remain deferred until footprint-aware covariance, random-field matching, and calibration/injection gates pass.",
        "",
        f"- Fixed lens-set membership SHA-256: `{fixed_selection['fixed_lens_ids_sha256']}`",
        f"- Approved source-row plan SHA-256: `{row_provenance['plan_sha256']}`",
        f"- Runner SHA-256: `{summary['provenance']['runner_sha256']}`",
        f"- Output directory: `{output}`", "",
    ])
    (output / "assessment.md").write_text("\n".join(assessment), encoding="utf-8")
    names = ["assessment.md", "per_lens_set_profiles.tsv", "block_deletion_support.tsv",
             "block_support_by_radius.tsv", "block_replicates.tsv", "estimator_profiles.tsv",
             "cross_component_diagnostics.tsv"]
    if covariance_rows:
        names.append("covariance_matrices.tsv")
    summary["output_sha256"] = {name: sha(output / name) for name in names}
    base.write_json(output / "summary.json", summary)
    print(json.dumps({
        "status": summary["status"], "science_measurement": False,
        "fixed_lens_sets": base.EXPECTED_FIXED_LENSES,
        "pooled_reproduces_reference": summary["sample"]["pooled_profile_reproduces_pinned_reference"],
        "block_widths": [
            {"width_deg": item["block_width_deg"], "occupied_cells": item["n_occupied_endpoint_cells"],
             "minimum_sets": item["minimum_common_sets_after_delete"],
             "estimators": {name: {"rank": result["covariance"].get("covariance_rank"),
                                    "PSD": result["covariance"].get("positive_semidefinite_with_relative_tolerance_1e-10"),
                                    "gate": result["covariance"].get("exploratory_covariance_quality_gate_passed")}
                             for name, result in item["estimators"].items()}}
            for item in per_width
        ],
        "output_directory": str(output),
        "summary_sha256": sha(output / "summary.json"),
    }, ensure_ascii=False, indent=2), flush=True)
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan-dir", type=Path, default=DEFAULT_PLAN_DIR)
    parser.add_argument("--rows", type=Path, default=DEFAULT_RESULTS)
    parser.add_argument("--support-dir", type=Path, default=base.SUPPORT_DIR)
    parser.add_argument("--qa-dir", type=Path, default=base.QA_DIR)
    parser.add_argument("--nz-sacc", type=Path, default=base.NZ_PATH)
    parser.add_argument("--pooled-summary", type=Path, default=DEFAULT_POOLED_SUMMARY)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--block-widths", type=float, nargs="+", default=[2.5, 5.0])
    args = parser.parse_args()
    if sorted(set(args.block_widths)) != [2.5, 5.0]:
        parser.error("this analysis requires exactly the 2.5- and 5-degree block widths")
    args.block_widths = sorted(args.block_widths)
    run(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
