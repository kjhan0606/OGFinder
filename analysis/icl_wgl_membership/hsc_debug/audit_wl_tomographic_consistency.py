#!/usr/bin/env python3
"""Test cluster-only DeltaSigma consistency across HSC Y3 source bins.

Uses only the previously approved lens-source rows and the local HSC Y3 SACC
file. A joint sky-cell jackknife keeps covariance between source bins and
radial bins. This remains exploratory until matched random subtraction and
redshift-calibration uncertainties are included.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import sys

import numpy as np
from scipy.stats import chi2

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
ROOT = HERE.parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from fetch_ancillary_catalogues import sha
from ogfkit.wlcluster import Cosmology
from run_wl_cluster_spatial_jackknife import (
    delete_block_covariance,
    equal_area_block_ids,
    _profile_vector,
)
from run_wl_shape_stack import verify_and_load_pairs
from s19a_shape_stack import load_hsc_y3_nz_sacc, stack_s19a_tomographic_pairs


def fit_common_value(values, covariance):
    """GLS consistency test for several estimates of one common quantity."""
    values = np.asarray(values, dtype=float)
    covariance = np.asarray(covariance, dtype=float)
    if values.ndim != 1 or len(values) < 2 or covariance.shape != (len(values), len(values)):
        raise ValueError("values and covariance must describe at least two estimates")
    if not np.all(np.isfinite(values)) or not np.all(np.isfinite(covariance)):
        raise ValueError("values and covariance must be finite")
    covariance = 0.5 * (covariance + covariance.T)
    inverse = np.linalg.pinv(covariance, hermitian=True)
    ones = np.ones(len(values))
    denominator = float(ones @ inverse @ ones)
    if denominator <= 0:
        raise ValueError("covariance has no positive information on the common mode")
    common = float((ones @ inverse @ values) / denominator)
    residual = values - common
    statistic = float(residual @ inverse @ residual)
    rank = int(np.linalg.matrix_rank(covariance))
    dof = rank - 1
    return {
        "common_value": common,
        "chi2": statistic,
        "degrees_of_freedom": dof,
        "p_value": float(chi2.sf(statistic, dof)) if dof > 0 else None,
        "covariance_rank": rank,
    }


def tomographic_spatial_jackknife(pairs, nz_redshift, nz_by_bin, radial_edges_mpc,
                                  cosmology, *, cell_width_deg=5.0):
    """Stack each HSC source bin and jointly jackknife all profile components."""
    required = (
        "source_zbin", "lens_ra_deg", "lens_dec_deg", "source_ra_deg", "source_dec_deg",
    )
    missing = [key for key in required if key not in pairs]
    if missing:
        raise ValueError(f"pair table lacks fields: {missing}")
    bins_raw = np.asarray(pairs["source_zbin"])
    try:
        bins_numeric = bins_raw.astype(float)
    except (TypeError, ValueError, OverflowError):
        raise ValueError("source_zbin must be integer values") from None
    if (not np.all(np.isfinite(bins_numeric))
            or not np.array_equal(bins_numeric, bins_numeric.astype(int))
            or np.any((bins_numeric < 0) | (bins_numeric > 4))):
        raise ValueError("source_zbin must use HSC Y3 integer values 0..4")
    bins = bins_numeric.astype(int)
    unassigned_count = int(np.count_nonzero(bins == 0))
    assigned = bins > 0
    if not np.any(assigned):
        raise ValueError("no source rows belong to HSC Y3 bins 1..4")
    pairs = {key: np.asarray(value)[assigned] for key, value in pairs.items()}
    bins = bins[assigned]
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
    bin_masks = []
    for source_bin in range(1, 5):
        mask = bins == source_bin
        if not np.any(mask):
            raise ValueError(f"source tomographic bin {source_bin} has no input rows")
        bin_masks.append(mask)
        subset = {key: np.asarray(pairs[key])[mask] for key in all_keys}
        full_results.append(stack_s19a_tomographic_pairs(
            subset, nz_redshift, nz_by_bin, radial_edges_mpc, cosmology,
        ))

    full_vector = np.concatenate([_profile_vector(result) for result in full_results])
    replicates = []
    deletion_rows = []
    for block in block_ids:
        keep = (lens_cells != block) & (source_cells != block)
        replicate_vectors = []
        remaining_by_bin = []
        deleted_rows = 0
        for mask in bin_masks:
            keep_bin = keep[mask]
            subset = {
                key: np.asarray(pairs[key])[mask][keep_bin]
                for key in all_keys
            }
            replicate = stack_s19a_tomographic_pairs(
                subset, nz_redshift, nz_by_bin, radial_edges_mpc, cosmology,
            )
            replicate_vectors.append(_profile_vector(replicate))
            remaining_by_bin.append(int(replicate["counts"]["used_pairs"]))
            deleted_rows += int(np.count_nonzero(~keep_bin))
        replicates.append(np.concatenate(replicate_vectors))
        deletion_rows.append({
            "block_id": int(block),
            "deleted_pair_rows": deleted_rows,
            "used_pairs_remaining_by_source_bin": remaining_by_bin,
        })
    covariance = delete_block_covariance(np.asarray(replicates, dtype=float))
    n_radial = len(radial_edges_mpc) - 1
    profile = []
    tangent_indices_by_radial = []
    for bin_index, result in enumerate(full_results):
        for radial_index, row in enumerate(result["profile"]):
            vector_index = bin_index * 2 * n_radial + 2 * radial_index
            profile.append({
                "source_zbin": bin_index + 1,
                **row,
                "delta_sigma_jackknife_error_msun_mpc2": float(np.sqrt(covariance[vector_index, vector_index])),
                "cross_jackknife_error_msun_mpc2": float(np.sqrt(covariance[vector_index + 1, vector_index + 1])),
            })
    for radial_index in range(n_radial):
        indices = [bin_index * 2 * n_radial + 2 * radial_index for bin_index in range(4)]
        tangent_indices_by_radial.append(indices)
    consistency = []
    for radial_index, indices in enumerate(tangent_indices_by_radial):
        values = [float(full_vector[index]) for index in indices]
        block = covariance[np.ix_(indices, indices)]
        consistency.append({
            "radial_bin": radial_index,
            "r_min_mpc": float(radial_edges_mpc[radial_index]),
            "r_max_mpc": float(radial_edges_mpc[radial_index + 1]),
            "source_bin_delta_sigma_msun_mpc2": values,
            **fit_common_value(values, block),
        })
    cross_indices = [
        bin_index * 2 * n_radial + 2 * radial_index + 1
        for bin_index in range(4) for radial_index in range(n_radial)
    ]
    cross = full_vector[cross_indices]
    cross_cov = covariance[np.ix_(cross_indices, cross_indices)]
    cross_inverse = np.linalg.pinv(cross_cov, hermitian=True)
    cross_dof = int(np.linalg.matrix_rank(cross_cov))
    cross_chi2 = float(cross @ cross_inverse @ cross)
    return {
        "schema": "astrafek-hsc-y3-tomographic-consistency-jackknife-v1",
        "status": "exploratory",
        "science_measurement": False,
        "random_lens_subtraction_applied": False,
        "n_z_uncertainty_propagated": False,
        "method": "four separate HSC Y3 source-bin DeltaSigma stacks with one shared delete-one equal-area sky-cell jackknife; every replicate recomputes shape and selection calibration",
        "cell_width_equatorial_deg": float(cell_width_deg),
        "n_occupied_cells": int(len(block_ids)),
        "n_input_unassigned_source_bin_zero_excluded": unassigned_count,
        "n_input_pairs_in_tomographic_bins": int(len(bins)),
        "radial_edges_mpc": [float(value) for value in radial_edges_mpc],
        "source_bin_counts": [result["counts"] for result in full_results],
        "profile_order": "source_zbin 1..4, then radial bin, then tangential and cross",
        "profile": profile,
        "tomographic_delta_sigma_consistency": consistency,
        "joint_cross_null": {
            "chi2": cross_chi2,
            "degrees_of_freedom": cross_dof,
            "p_value": float(chi2.sf(cross_chi2, cross_dof)) if cross_dof > 0 else None,
            "covariance_rank": cross_dof,
        },
        "covariance": covariance.tolist(),
        "block_deletions": deletion_rows,
        "limitations": [
            "This is the cluster-only stack and has no random-lens subtraction.",
            "The HSC Y3 n(z) curves are fiducial and their posterior/systematic uncertainty is not propagated.",
            "Equal-area cells are an exploratory fallback, not connected HSC survey regions.",
            "The GLS tomography-consistency p-values use a Gaussian covariance approximation and are diagnostics, not detections.",
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
    result = tomographic_spatial_jackknife(
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
        "occupied_cells": result["n_occupied_cells"],
        "source_bin_used_pairs": [row["used_pairs"] for row in result["source_bin_counts"]],
        "tomography_consistency_chi2_by_radial_bin": [row["chi2"] for row in result["tomographic_delta_sigma_consistency"]],
        "output": str(summary),
    }, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
