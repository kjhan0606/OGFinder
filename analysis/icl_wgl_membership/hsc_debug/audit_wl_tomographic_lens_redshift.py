#!/usr/bin/env python3
"""Test HSC Y3 source-bin consistency in low/high lens-redshift halves.

The lens split is defined only from the unique lens redshifts and IDs, before
looking at shear profiles. All source-bin and lens-half profiles share the
same delete-one equal-area sky-cell jackknife. This is a diagnostic, not a
weak-lensing detection or calibrated mass measurement.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
ROOT = HERE.parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from audit_wl_tomographic_consistency import fit_common_value
from fetch_ancillary_catalogues import sha
from ogfkit.wlcluster import Cosmology
from run_wl_cluster_spatial_jackknife import (
    delete_block_covariance,
    equal_area_block_ids,
    _profile_vector,
)
from run_wl_shape_stack import verify_and_load_pairs
from s19a_shape_stack import load_hsc_y3_nz_sacc, stack_s19a_tomographic_pairs


def _lens_half_labels(lens_ids: np.ndarray, lens_z: np.ndarray):
    """Return a deterministic equal-count split ordered by (lens_z, lens_id)."""
    lens_ids = np.asarray(lens_ids).astype(str)
    lens_z = np.asarray(lens_z, dtype=float)
    if lens_ids.ndim != 1 or lens_z.shape != lens_ids.shape or len(lens_ids) < 6:
        raise ValueError("lens IDs/redshifts must be matching vectors for at least six lenses")
    if len(set(lens_ids.tolist())) != len(lens_ids):
        raise ValueError("lens IDs must be unique in the split input")
    if not np.all(np.isfinite(lens_z)) or np.any(lens_z < 0):
        raise ValueError("lens redshifts must be finite and nonnegative")
    order = np.lexsort((lens_ids, lens_z))
    split = len(order) // 2
    if split < 3 or len(order) - split < 3:
        raise ValueError("both lens-redshift halves must contain at least three lenses")
    labels = np.empty(len(order), dtype=np.int8)
    labels[order[:split]] = 0
    labels[order[split:]] = 1
    boundary = float((lens_z[order[split - 1]] + lens_z[order[split]]) / 2.0)
    return labels, {
        "split_rule": "sort unique lenses by (lens_z, lens_id); first floor(N/2) are lower-z",
        "n_lenses_total": int(len(order)),
        "n_lenses_lower_z": int(split),
        "n_lenses_upper_z": int(len(order) - split),
        "lower_z_max": float(lens_z[order[split - 1]]),
        "upper_z_min": float(lens_z[order[split]]),
        "midpoint_between_boundary_lenses": boundary,
        "ties_at_split_boundary": bool(lens_z[order[split - 1]] == lens_z[order[split]]),
    }


def holm_adjust(p_values):
    """Holm family-wise adjustment, returned in input order."""
    values = np.asarray(p_values, dtype=float)
    if values.ndim != 1 or not np.all(np.isfinite(values)) or np.any((values < 0) | (values > 1)):
        raise ValueError("p-values must be a finite vector in [0, 1]")
    order = np.argsort(values)
    adjusted = np.empty_like(values)
    running = 0.0
    count = len(values)
    for rank, index in enumerate(order):
        running = max(running, (count - rank) * values[index])
        adjusted[index] = min(1.0, running)
    return adjusted


def lens_redshift_tomographic_jackknife(
    pairs, nz_redshift, nz_by_bin, radial_edges_mpc, cosmology, *, cell_width_deg=5.0,
):
    """Stack four source bins in two lens-z halves using shared sky deletions."""
    required = (
        "lens_id", "lens_z", "source_zbin", "lens_ra_deg", "lens_dec_deg",
        "source_ra_deg", "source_dec_deg",
    )
    missing = [key for key in required if key not in pairs]
    if missing:
        raise ValueError(f"pair table lacks fields: {missing}")

    lens_ids_all = np.asarray(pairs["lens_id"]).astype(str)
    lens_z_all = np.asarray(pairs["lens_z"], dtype=float)
    if lens_ids_all.ndim != 1 or lens_z_all.shape != lens_ids_all.shape:
        raise ValueError("lens_id and lens_z must be matching one-dimensional arrays")
    unique_lens_ids, inverse = np.unique(lens_ids_all, return_inverse=True)
    unique_lens_z = np.empty(len(unique_lens_ids), dtype=float)
    for index in range(len(unique_lens_ids)):
        values = lens_z_all[inverse == index]
        if not np.all(np.isfinite(values)) or not np.all(values == values[0]):
            raise ValueError(f"lens redshift is invalid or inconsistent for {unique_lens_ids[index]}")
        unique_lens_z[index] = values[0]
    lens_labels, split = _lens_half_labels(unique_lens_ids, unique_lens_z)
    lens_label_by_id = dict(zip(unique_lens_ids.tolist(), lens_labels.tolist()))
    labels = np.fromiter((lens_label_by_id[key] for key in lens_ids_all), dtype=np.int8,
                         count=len(lens_ids_all))

    source_bins_raw = np.asarray(pairs["source_zbin"])
    try:
        source_bins_numeric = source_bins_raw.astype(float)
    except (TypeError, ValueError, OverflowError):
        raise ValueError("source_zbin must be integer values") from None
    if (not np.all(np.isfinite(source_bins_numeric))
            or not np.array_equal(source_bins_numeric, source_bins_numeric.astype(int))
            or np.any((source_bins_numeric < 0) | (source_bins_numeric > 4))):
        raise ValueError("source_zbin must use HSC Y3 integer values 0..4")
    source_bins = source_bins_numeric.astype(int)
    n_bin_zero = int(np.count_nonzero(source_bins == 0))
    assigned = source_bins > 0
    if not np.any(assigned):
        raise ValueError("no source rows belong to HSC Y3 bins 1..4")
    pairs = {key: np.asarray(value)[assigned] for key, value in pairs.items()}
    source_bins = source_bins[assigned]
    labels = labels[assigned]

    lens_cells = equal_area_block_ids(
        pairs["lens_ra_deg"], pairs["lens_dec_deg"], cell_width_deg=cell_width_deg,
    )
    source_cells = equal_area_block_ids(
        pairs["source_ra_deg"], pairs["source_dec_deg"], cell_width_deg=cell_width_deg,
    )
    block_ids = np.unique(np.concatenate((lens_cells, source_cells)))
    if len(block_ids) < 3:
        raise ValueError("spatial jackknife requires at least three occupied sky cells")

    all_keys = tuple(pairs)
    full_results = []
    cohort_bin_masks = []
    for cohort in range(2):
        masks = []
        results = []
        for source_bin in range(1, 5):
            mask = (labels == cohort) & (source_bins == source_bin)
            if not np.any(mask):
                raise ValueError(f"lens cohort {cohort} source bin {source_bin} has no input rows")
            masks.append(mask)
            subset = {key: np.asarray(pairs[key])[mask] for key in all_keys}
            results.append(stack_s19a_tomographic_pairs(
                subset, nz_redshift, nz_by_bin, radial_edges_mpc, cosmology,
            ))
        cohort_bin_masks.append(masks)
        full_results.append(results)

    full_vector = np.concatenate([
        _profile_vector(result) for results in full_results for result in results
    ])
    replicates = []
    deletion_rows = []
    for block in block_ids:
        keep = (lens_cells != block) & (source_cells != block)
        replicate_vectors = []
        remaining_by_cohort_bin = []
        deleted_rows = 0
        for masks in cohort_bin_masks:
            for mask in masks:
                keep_bin = keep[mask]
                subset = {
                    key: np.asarray(pairs[key])[mask][keep_bin]
                    for key in all_keys
                }
                replicate = stack_s19a_tomographic_pairs(
                    subset, nz_redshift, nz_by_bin, radial_edges_mpc, cosmology,
                )
                replicate_vectors.append(_profile_vector(replicate))
                remaining_by_cohort_bin.append(int(replicate["counts"]["used_pairs"]))
                deleted_rows += int(np.count_nonzero(~keep_bin))
        replicates.append(np.concatenate(replicate_vectors))
        deletion_rows.append({
            "block_id": int(block),
            "deleted_pair_rows": deleted_rows,
            "used_pairs_remaining_lower_z_source_bins_1_to_4": remaining_by_cohort_bin[:4],
            "used_pairs_remaining_upper_z_source_bins_1_to_4": remaining_by_cohort_bin[4:],
        })
    covariance = delete_block_covariance(np.asarray(replicates, dtype=float))
    n_radial = len(radial_edges_mpc) - 1
    profile = []
    consistency = []
    raw_p_values = []
    for cohort, cohort_name in enumerate(("lower_z", "upper_z")):
        cohort_base = cohort * 4 * 2 * n_radial
        for bin_index, result in enumerate(full_results[cohort]):
            for radial_index, row in enumerate(result["profile"]):
                vector_index = cohort_base + bin_index * 2 * n_radial + 2 * radial_index
                profile.append({
                    "lens_redshift_cohort": cohort_name,
                    "source_zbin": bin_index + 1,
                    **row,
                    "delta_sigma_jackknife_error_msun_mpc2": float(np.sqrt(covariance[vector_index, vector_index])),
                    "cross_jackknife_error_msun_mpc2": float(np.sqrt(covariance[vector_index + 1, vector_index + 1])),
                })
        for radial_index in range(n_radial):
            indices = [
                cohort_base + bin_index * 2 * n_radial + 2 * radial_index
                for bin_index in range(4)
            ]
            values = [float(full_vector[index]) for index in indices]
            fit = fit_common_value(values, covariance[np.ix_(indices, indices)])
            raw_p_values.append(fit["p_value"])
            consistency.append({
                "lens_redshift_cohort": cohort_name,
                "radial_bin": radial_index,
                "r_min_mpc": float(radial_edges_mpc[radial_index]),
                "r_max_mpc": float(radial_edges_mpc[radial_index + 1]),
                "source_bin_delta_sigma_msun_mpc2": values,
                **fit,
            })
    # Keep all 14 predeclared tests in the family. A non-computable p-value
    # contributes the conservative bound 1 to the adjustment.
    adjusted = holm_adjust([1.0 if value is None else value for value in raw_p_values])
    adjusted_by_index = {
        index: (None if raw_p_values[index] is None else float(adjusted[index]))
        for index in range(len(raw_p_values))
    }
    for index, row in enumerate(consistency):
        row["holm_p_value_14_tests"] = adjusted_by_index.get(index)

    return {
        "schema": "astrafek-hsc-y3-lens-redshift-tomographic-jackknife-v1",
        "status": "exploratory",
        "science_measurement": False,
        "random_lens_subtraction_applied": False,
        "n_z_uncertainty_propagated": False,
        "method": "four HSC Y3 source-bin DeltaSigma stacks in two equal-count lens-redshift cohorts; one shared delete-one equal-area sky-cell jackknife; each replicate recomputes shape and selection calibration",
        "lens_redshift_split": split,
        "cell_width_equatorial_deg": float(cell_width_deg),
        "n_occupied_cells": int(len(block_ids)),
        "n_input_unassigned_source_bin_zero_excluded": n_bin_zero,
        "n_input_pairs_in_tomographic_bins": int(len(source_bins)),
        "radial_edges_mpc": [float(value) for value in radial_edges_mpc],
        "source_bin_counts_by_lens_redshift_cohort": [
            {"lens_redshift_cohort": name,
             "counts_by_source_bin": [result["counts"] for result in results]}
            for name, results in zip(("lower_z", "upper_z"), full_results)
        ],
        "profile_order": "lens_redshift_cohort lower_z then upper_z; source_zbin 1..4; radial bin; tangential and cross",
        "profile": profile,
        "tomographic_delta_sigma_consistency": consistency,
        "covariance_shape": list(covariance.shape),
        "covariance_rank": int(np.linalg.matrix_rank(covariance)),
        "block_deletions": deletion_rows,
        "limitations": [
            "This is the cluster-only stack and has no random-lens subtraction.",
            "The HSC Y3 n(z) curves are fiducial; their posterior and systematic uncertainty is not propagated.",
            "Equal-area cells are an exploratory fallback, not connected HSC survey regions.",
            "The equal-count lens-redshift halves are a diagnostic split, not externally defined physical populations.",
            "Holm adjustment covers the 14 radial-bin consistency tests; p-values are diagnostics, not detections.",
        ],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rows", type=Path, required=True)
    parser.add_argument("--preflight", type=Path, required=True)
    parser.add_argument("--plan-dir", type=Path, required=True)
    parser.add_argument("--nz-sacc", type=Path, required=True)
    parser.add_argument("--masked-value", type=int, choices=(0, 1), required=True)
    parser.add_argument("--cell-width-deg", type=float, default=5.0)
    parser.add_argument("--radial-edges-mpc", type=float, nargs="+", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--h0", type=float, default=70.0)
    parser.add_argument("--omega-m", type=float, default=0.3)
    args = parser.parse_args()
    if args.output.exists():
        parser.error(f"refusing to overwrite existing output directory: {args.output}")
    if len(args.radial_edges_mpc) < 2 or any(
            right <= left for left, right in zip(args.radial_edges_mpc, args.radial_edges_mpc[1:])):
        parser.error("radial edges must be strictly increasing and contain at least two values")
    pairs, provenance = verify_and_load_pairs(
        args.rows, args.preflight, args.plan_dir, args.masked_value,
    )
    nz_redshift, nz_by_bin = load_hsc_y3_nz_sacc(args.nz_sacc)
    cosmology = Cosmology(H0=args.h0, Om=args.omega_m)
    result = lens_redshift_tomographic_jackknife(
        pairs, nz_redshift, nz_by_bin, args.radial_edges_mpc, cosmology,
        cell_width_deg=args.cell_width_deg,
    )
    args.output.mkdir(parents=True, exist_ok=False)
    payload = {
        **result,
        "input": {
            **provenance,
            "n_z_sacc_path": str(args.nz_sacc.resolve()),
            "n_z_sacc_sha256": sha(args.nz_sacc),
            "runner_sha256": sha(Path(__file__).resolve()),
            "tomographic_runner_sha256": sha(Path(__file__).with_name("audit_wl_tomographic_consistency.py")),
            "shape_adapter_sha256": sha(Path(__file__).with_name("s19a_shape_stack.py")),
            "spatial_jackknife_sha256": sha(Path(__file__).with_name("run_wl_cluster_spatial_jackknife.py")),
            "radial_edges_mpc": [float(value) for value in args.radial_edges_mpc],
            "cosmology": {"H0_km_s_mpc": args.h0, "Omega_m": args.omega_m},
            "b_mode_masked_raw_value": args.masked_value,
        },
    }
    summary = args.output / "summary.json"
    with summary.open("x", encoding="utf-8") as stream:
        json.dump(payload, stream, indent=2, ensure_ascii=False, allow_nan=False)
        stream.write("\n")
    print(json.dumps({
        "status": result["status"],
        "science_measurement": result["science_measurement"],
        "lens_redshift_split": result["lens_redshift_split"],
        "occupied_cells": result["n_occupied_cells"],
        "source_bin_used_pairs": result["source_bin_counts_by_lens_redshift_cohort"],
        "consistency": result["tomographic_delta_sigma_consistency"],
        "output": str(summary),
    }, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
