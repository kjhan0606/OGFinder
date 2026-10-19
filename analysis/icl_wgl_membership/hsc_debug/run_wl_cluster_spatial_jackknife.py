#!/usr/bin/env python3
"""Estimate cluster-only HSC Y3 stack covariance with sky-block deletion.

This is an intermediate diagnostic while the matched random-lens stack is
being assembled. It jointly removes lens-source pairs whenever either end is
in a deleted equal-area sky cell, so repeated sources and lenses in the same
cell are not treated as independent rows. The result remains exploratory.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.stats import chi2

from run_wl_shape_stack import verify_and_load_pairs
from fetch_ancillary_catalogues import sha
from s19a_shape_stack import load_hsc_y3_nz_sacc, stack_s19a_tomographic_pairs
from ogfkit.wlcluster import Cosmology


def equal_area_block_ids(ra_deg, dec_deg, *, cell_width_deg=5.0):
    """Map finite sky coordinates to reproducible approximately square cells.

    RA is divided into ``cell_width_deg`` intervals and sin(Dec) into intervals
    with the same equatorial angular width. This gives equal-area cells and is
    an exploratory spatial jackknife, not a survey-region partition.
    """
    ra = np.asarray(ra_deg, dtype=float)
    dec = np.asarray(dec_deg, dtype=float)
    width = float(cell_width_deg)
    if (ra.ndim != 1 or dec.shape != ra.shape or not np.all(np.isfinite(ra))
            or not np.all(np.isfinite(dec)) or np.any((ra < 0) | (ra >= 360))
            or np.any(np.abs(dec) > 90)):
        raise ValueError("sky coordinates must be finite RA/Dec vectors in degrees")
    if not np.isfinite(width) or not 0 < width <= 180:
        raise ValueError("cell_width_deg must be finite and in (0, 180]")
    n_ra = int(round(360.0 / width))
    n_sin_dec = max(1, int(round(2.0 / np.sin(np.deg2rad(width)))))
    ra_width = 360.0 / n_ra
    ra_index = np.floor(((ra + 0.5 * ra_width) % 360.0) / ra_width).astype(int)
    sin_dec = np.sin(np.deg2rad(dec))
    dec_index = np.minimum(((sin_dec + 1.0) * 0.5 * n_sin_dec).astype(int), n_sin_dec - 1)
    return dec_index * n_ra + ra_index


def delete_block_covariance(replicates):
    """Return standard delete-one-block covariance for [block, observable]."""
    values = np.asarray(replicates, dtype=float)
    if values.ndim != 2 or values.shape[0] < 3 or values.shape[1] < 1:
        raise ValueError("replicates must have shape [at least 3 blocks, observables]")
    if not np.all(np.isfinite(values)):
        raise ValueError("jackknife replicates must be finite")
    mean = values.mean(axis=0)
    centered = values - mean
    covariance = (values.shape[0] - 1.0) / values.shape[0] * centered.T @ centered
    return covariance


def _profile_vector(result):
    vector = []
    for row in result["profile"]:
        tangential = row["delta_sigma_weak_msun_mpc2"]
        cross = row["cross_delta_sigma_weak_msun_mpc2"]
        if tangential is None or cross is None or not np.isfinite(tangential) or not np.isfinite(cross):
            raise ValueError("every radial bin must remain populated in each leave-one-block run")
        vector.extend((float(tangential), float(cross)))
    return np.asarray(vector, dtype=float)


def cluster_spatial_jackknife(pairs, nz_redshift, nz_by_bin, radial_edges_mpc,
                              cosmology, *, cell_width_deg=5.0):
    """Recompute the full calibrated cluster stack after each spatial deletion."""
    required = ("lens_ra_deg", "lens_dec_deg", "source_ra_deg", "source_dec_deg")
    missing = [key for key in required if key not in pairs]
    if missing:
        raise ValueError(f"pair table lacks sky coordinates: {missing}")
    lens_cells = equal_area_block_ids(
        pairs["lens_ra_deg"], pairs["lens_dec_deg"], cell_width_deg=cell_width_deg,
    )
    source_cells = equal_area_block_ids(
        pairs["source_ra_deg"], pairs["source_dec_deg"], cell_width_deg=cell_width_deg,
    )
    block_ids = np.unique(np.concatenate((lens_cells, source_cells)))
    if len(block_ids) < 3:
        raise ValueError("spatial jackknife requires at least three occupied sky cells")

    full = stack_s19a_tomographic_pairs(
        pairs, nz_redshift, nz_by_bin, radial_edges_mpc, cosmology,
    )
    full_vector = _profile_vector(full)
    replicates = []
    deletion_rows = []
    for block in block_ids:
        keep = (lens_cells != block) & (source_cells != block)
        replicate = stack_s19a_tomographic_pairs(
            {key: np.asarray(value)[keep] for key, value in pairs.items()},
            nz_redshift, nz_by_bin, radial_edges_mpc, cosmology,
        )
        replicates.append(_profile_vector(replicate))
        deletion_rows.append({
            "block_id": int(block),
            "deleted_pair_rows": int(np.count_nonzero(~keep)),
            "deleted_lens_ids": int(len(set(np.asarray(pairs["lens_id"])[lens_cells == block]))),
            "deleted_source_ids": int(len(set(np.asarray(pairs["object_id"])[source_cells == block]))),
            "used_pairs_remaining": int(replicate["counts"]["used_pairs"]),
        })
    replicate_array = np.asarray(replicates, dtype=float)
    covariance = delete_block_covariance(replicate_array)
    errors = np.sqrt(np.diag(covariance))
    profile = []
    cross_indices = []
    for index, row in enumerate(full["profile"]):
        tangential_index, cross_index = 2 * index, 2 * index + 1
        cross_indices.append(cross_index)
        profile.append({
            **row,
            "delta_sigma_jackknife_error_msun_mpc2": float(errors[tangential_index]),
            "cross_delta_sigma_jackknife_error_msun_mpc2": float(errors[cross_index]),
            "cross_jackknife_z": float(row["cross_delta_sigma_weak_msun_mpc2"] / errors[cross_index])
            if errors[cross_index] > 0 else None,
        })
    cross = full_vector[cross_indices]
    cross_cov = covariance[np.ix_(cross_indices, cross_indices)]
    inverse = np.linalg.pinv(cross_cov, hermitian=True)
    cross_chi2 = float(cross @ inverse @ cross)
    cross_dof = int(np.linalg.matrix_rank(cross_cov))
    return {
        "schema": "astrafek-hsc-y3-cluster-spatial-jackknife-v1",
        "status": "exploratory",
        "science_measurement": False,
        "random_lens_subtraction_applied": False,
        "photoz_uncertainty_propagated": False,
        "reduced_shear_applied": False,
        "method": "delete one equal-area sky cell; remove every pair whose lens or source lies in the deleted cell; recompute all per-bin calibration",
        "cell_width_equatorial_deg": float(cell_width_deg),
        "n_occupied_cells": int(len(block_ids)),
        "full_stack_counts": full["counts"],
        "profile": profile,
        "profile_order": "[DeltaSigma bin 0, cross bin 0, DeltaSigma bin 1, cross bin 1, ...]",
        "covariance": covariance.tolist(),
        "cross_null": {
            "chi2": cross_chi2,
            "degrees_of_freedom": cross_dof,
            "p_value": float(chi2.sf(cross_chi2, cross_dof)) if cross_dof > 0 else None,
            "covariance_rank": cross_dof,
            "reference_distribution": "chi-square under a Gaussian covariance approximation; diagnostic only",
        },
        "block_deletions": deletion_rows,
        "limitations": [
            "equal-area grid cells are a fallback and may combine disconnected survey regions",
            "cluster-only covariance does not include random-lens subtraction or matched field variance",
            "fiducial n(z), shear calibration, source selection, and reduced-shear uncertainties are not propagated",
            "cross-shear chi-square is a diagnostic and does not establish a science detection",
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
    if not args.nz_sacc.is_file():
        parser.error(f"SACC file does not exist: {args.nz_sacc}")
    if len(args.radial_edges_mpc) < 2 or any(
            right <= left for left, right in zip(args.radial_edges_mpc, args.radial_edges_mpc[1:])):
        parser.error("radial edges must be strictly increasing and contain at least two values")
    pairs, provenance = verify_and_load_pairs(
        args.rows, args.preflight, args.plan_dir, args.masked_value,
    )
    nz_redshift, nz_by_bin = load_hsc_y3_nz_sacc(args.nz_sacc)
    cosmology = Cosmology(H0=args.h0, Om=args.omega_m)
    result = cluster_spatial_jackknife(
        pairs, nz_redshift, nz_by_bin, args.radial_edges_mpc, cosmology,
        cell_width_deg=args.cell_width_deg,
    )
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    payload = {
        **result,
        "input": {
            **provenance,
            "n_z_sacc_path": str(args.nz_sacc.resolve()),
            "n_z_sacc_sha256": sha(args.nz_sacc),
            "radial_edges_mpc": [float(value) for value in args.radial_edges_mpc],
            "cosmology": {"H0_km_s_mpc": args.h0, "Omega_m": args.omega_m},
            "b_mode_masked_raw_value": args.masked_value,
        },
    }
    summary = output / "summary.json"
    with summary.open("x", encoding="utf-8") as stream:
        json.dump(payload, stream, indent=2, ensure_ascii=False, allow_nan=False)
        stream.write("\n")
    print(json.dumps({
        "status": result["status"], "science_measurement": False,
        "occupied_cells": result["n_occupied_cells"],
        "used_pairs": result["full_stack_counts"]["used_pairs"],
        "cross_null": result["cross_null"], "output": str(summary),
    }, ensure_ascii=False))


if __name__ == "__main__":
    main()
