#!/usr/bin/env python3
"""Measure how source-object reuse affects HSC cluster-stack jackknife errors.

Using only the already-approved lens--source pair export, compare three
spatial delete-one covariance constructions: lens-endpoint cells, source-
endpoint cells, and joint cells deleting a pair when either endpoint belongs
to the removed cell. Reused source objects therefore leave together in the
source and joint variants. This is a covariance sensitivity, not a science
measurement or a substitute for matched random-lens subtraction.
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

from audit_wl_tomographic_consistency import delete_block_covariance, equal_area_block_ids, fit_common_value
from audit_wl_tomographic_lens_redshift import holm_adjust
from fetch_ancillary_catalogues import sha
from ogfkit.wlcluster import Cosmology
from run_wl_shape_stack import verify_and_load_pairs
from s19a_shape_stack import (
    load_hsc_y3_nz_sacc,
    pair_geometry,
    stack_s19a_tomographic_pairs,
    tomographic_nz_lensing_efficiency,
)


def selected_pair_mask(pairs, nz_redshift, nz_by_bin, radial_edges_mpc, cosmology):
    """Reproduce the published shape, B-mode, n(z), and radius selections."""
    bins = np.asarray(pairs["source_zbin"], dtype=int)
    e1, e2 = np.asarray(pairs["e1"], float), np.asarray(pairs["e2"], float)
    sigma = np.asarray(pairs["sigma_e"], float)
    magnitude = np.asarray(pairs["mag_a10"], float)
    resolution = np.asarray(pairs["resolution"], float)
    shape = ((sigma >= 0.0) & (sigma <= 0.4) & (np.hypot(e1, e2) < 2.0)
             & (magnitude <= 25.5) & (resolution >= 0.3))
    bmode = np.asarray(pairs["b_mode_pass"], dtype=bool)
    radius, _ = pair_geometry(pairs, cosmology)
    beta = tomographic_nz_lensing_efficiency(
        pairs["lens_z"], bins, nz_redshift, nz_by_bin, cosmology,
    )
    return ((bins >= 1) & (bins <= 4) & shape & bmode & (beta > 0.0)
            & (radius >= radial_edges_mpc[0]) & (radius <= radial_edges_mpc[-1]))


def source_reuse_summary(pairs, used, cell_width_deg):
    ids = np.asarray(pairs["object_id"], dtype=str)[used]
    if not len(ids):
        raise ValueError("no selected pairs remain for source-reuse audit")
    _, inverse, multiplicity = np.unique(ids, return_inverse=True, return_counts=True)
    repeated = multiplicity > 1
    lens_cells = equal_area_block_ids(
        np.asarray(pairs["lens_ra_deg"], float)[used],
        np.asarray(pairs["lens_dec_deg"], float)[used],
        cell_width_deg=cell_width_deg,
    )
    source_cells = equal_area_block_ids(
        np.asarray(pairs["source_ra_deg"], float)[used],
        np.asarray(pairs["source_dec_deg"], float)[used],
        cell_width_deg=cell_width_deg,
    )
    cell_values, lens_cell_inverse = np.unique(lens_cells, return_inverse=True)
    _, source_cell_inverse = np.unique(source_cells, return_inverse=True)
    n_cells = int(max(len(cell_values), int(source_cell_inverse.max()) + 1))
    source_cell_codes = np.unique(inverse * n_cells + source_cell_inverse)
    source_cell_counts = np.bincount(source_cell_codes // n_cells, minlength=len(multiplicity))
    lens_cell_codes = np.unique(inverse * n_cells + lens_cell_inverse)
    lens_cell_counts = np.bincount(lens_cell_codes // n_cells, minlength=len(multiplicity))
    return {
        "selected_pair_rows": int(len(ids)),
        "unique_source_objects": int(len(multiplicity)),
        "source_objects_reused_across_lens_pairs": int(np.count_nonzero(repeated)),
        "fraction_unique_sources_reused": float(np.mean(repeated)),
        "pair_rows_from_reused_sources": int(multiplicity[repeated].sum()),
        "fraction_pair_rows_from_reused_sources": float(multiplicity[repeated].sum() / len(ids)),
        "maximum_lens_pair_multiplicity_per_source": int(multiplicity.max()),
        "reused_sources_spanning_multiple_lens_endpoint_cells": int(
            np.count_nonzero(repeated & (lens_cell_counts > 1))),
        "reused_sources_spanning_multiple_source_endpoint_cells": int(
            np.count_nonzero(repeated & (source_cell_counts > 1))),
        "unique_lens_endpoint_cells": int(len(cell_values)),
        "unique_source_endpoint_cells": int(len(np.unique(source_cells))),
        "source_cell_width_deg": float(cell_width_deg),
        "source_endpoint_cell_assignment_is_unique_per_object_id": bool(
            np.all(source_cell_counts == 1)),
    }


def jackknife_mode(pairs, nz_redshift, nz_by_bin, radial_edges_mpc, cosmology,
                  *, mode, cell_width_deg):
    bins = np.asarray(pairs["source_zbin"], dtype=int)
    assigned = bins > 0
    work = {key: np.asarray(value)[assigned] for key, value in pairs.items()}
    bins = bins[assigned]
    lens_cells = equal_area_block_ids(
        work["lens_ra_deg"], work["lens_dec_deg"], cell_width_deg=cell_width_deg,
    )
    source_cells = equal_area_block_ids(
        work["source_ra_deg"], work["source_dec_deg"], cell_width_deg=cell_width_deg,
    )
    if mode == "lens_only":
        blocks = np.unique(lens_cells)
    elif mode == "source_only":
        blocks = np.unique(source_cells)
    elif mode == "joint_endpoints":
        blocks = np.unique(np.concatenate((lens_cells, source_cells)))
    else:
        raise ValueError(f"unknown jackknife mode: {mode}")
    if len(blocks) < 3:
        raise ValueError(f"{mode} has fewer than three occupied deletion blocks")

    bin_masks = [bins == source_bin for source_bin in range(1, 5)]
    full_results = [stack_s19a_tomographic_pairs(
        {key: value[mask] for key, value in work.items()},
        nz_redshift, nz_by_bin, radial_edges_mpc, cosmology,
    ) for mask in bin_masks]
    full_vector = np.concatenate([
        np.asarray([
            component
            for row in result["profile"]
            for component in (row["delta_sigma_weak_msun_mpc2"],
                              row["cross_delta_sigma_weak_msun_mpc2"])
        ], dtype=float)
        for result in full_results
    ])

    replicates = []
    deletion_rows = []
    for block in blocks:
        if mode == "lens_only":
            keep = lens_cells != block
        elif mode == "source_only":
            keep = source_cells != block
        else:
            keep = (lens_cells != block) & (source_cells != block)
        vector_parts = []
        remaining = []
        for mask in bin_masks:
            selected = mask & keep
            result = stack_s19a_tomographic_pairs(
                {key: value[selected] for key, value in work.items()},
                nz_redshift, nz_by_bin, radial_edges_mpc, cosmology,
            )
            vector_parts.extend(
                component
                for row in result["profile"]
                for component in (row["delta_sigma_weak_msun_mpc2"],
                                  row["cross_delta_sigma_weak_msun_mpc2"])
            )
            remaining.append(int(result["counts"]["used_pairs"]))
        replicates.append(vector_parts)
        deletion_rows.append({
            "block_id": int(block),
            "deleted_pair_rows": int(np.count_nonzero(~keep)),
            "used_pairs_remaining_by_source_bin": remaining,
        })

    covariance = delete_block_covariance(np.asarray(replicates, dtype=float))
    n_radial = len(radial_edges_mpc) - 1
    consistency = []
    profile = []
    for source_index, result in enumerate(full_results):
        for radial_index, row in enumerate(result["profile"]):
            position = source_index * 2 * n_radial + 2 * radial_index
            profile.append({
                "source_zbin": source_index + 1,
                **row,
                "delta_sigma_jackknife_error_msun_mpc2": float(np.sqrt(covariance[position, position])),
                "cross_jackknife_error_msun_mpc2": float(np.sqrt(covariance[position + 1, position + 1])),
            })
    for radial_index in range(n_radial):
        indices = [source_index * 2 * n_radial + 2 * radial_index
                   for source_index in range(4)]
        values = [float(full_vector[index]) for index in indices]
        consistency.append({
            "radial_bin": radial_index,
            "r_min_mpc": float(radial_edges_mpc[radial_index]),
            "r_max_mpc": float(radial_edges_mpc[radial_index + 1]),
            "source_bin_delta_sigma_msun_mpc2": values,
            **fit_common_value(values, covariance[np.ix_(indices, indices)]),
        })
    adjusted = holm_adjust([row["p_value"] for row in consistency])
    for row, p_value in zip(consistency, adjusted):
        row["holm_p_value_7_radial_tests"] = float(p_value)
    return {
        "mode": mode,
        "n_deletion_blocks": int(len(blocks)),
        "covariance_dimension": int(covariance.shape[0]),
        "covariance_rank": int(np.linalg.matrix_rank(covariance)),
        "relative_covariance_symmetry_error": float(
            np.max(np.abs(covariance - covariance.T))
            / max(float(np.max(np.abs(covariance))), 1e-300)),
        "profile": profile,
        "tomographic_consistency": consistency,
        "covariance": covariance.tolist(),
        "block_deletions": deletion_rows,
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
        parser.error("radial edges must be strictly increasing")

    pairs, provenance = verify_and_load_pairs(
        args.rows, args.preflight, args.plan_dir, args.masked_value,
    )
    nz_redshift, nz_by_bin = load_hsc_y3_nz_sacc(args.nz_sacc)
    cosmology = Cosmology(H0=args.h0, Om=args.omega_m)
    used = selected_pair_mask(pairs, nz_redshift, nz_by_bin,
                              args.radial_edges_mpc, cosmology)
    reuse = source_reuse_summary(pairs, used, args.cell_width_deg)
    modes = {}
    for mode in ("lens_only", "source_only", "joint_endpoints"):
        print(f"Running {mode} delete-one covariance", flush=True)
        modes[mode] = jackknife_mode(
            pairs, nz_redshift, nz_by_bin, args.radial_edges_mpc, cosmology,
            mode=mode, cell_width_deg=args.cell_width_deg,
        )
    payload = {
        "schema": "astrafek-hsc-y3-source-reuse-covariance-audit-v1",
        "status": "verified_covariance_sensitivity",
        "science_measurement": False,
        "random_lens_subtraction_applied": False,
        "method": "compare lens-only, source-only, and joint endpoint equal-area delete-one jackknife covariances; report exact source-object reuse among selected lens-source pairs",
        "cell_width_deg": args.cell_width_deg,
        "source_reuse": reuse,
        "covariance_modes": modes,
        "input": {
            **provenance,
            "n_z_sacc_path": str(args.nz_sacc.resolve()),
            "n_z_sacc_sha256": sha(args.nz_sacc),
            "runner_sha256": sha(Path(__file__).resolve()),
            "tomographic_consistency_sha256": sha(HERE / "audit_wl_tomographic_consistency.py"),
            "spatial_jackknife_sha256": sha(HERE / "run_wl_cluster_spatial_jackknife.py"),
            "shape_adapter_sha256": sha(HERE / "s19a_shape_stack.py"),
            "radial_edges_mpc": [float(value) for value in args.radial_edges_mpc],
            "cosmology": {"H0_km_s_mpc": args.h0, "Omega_m": args.omega_m},
        },
        "limitations": [
            "Spatial delete-one covariance measures large-scale block sensitivity, not a complete source-level covariance model.",
            "The source-cell construction is an exploratory equal-area partition, not connected HSC survey regions.",
            "No random-lens subtraction, calibrated n(z) posterior propagation, or reduced-shear correction is applied.",
        ],
    }
    args.output.mkdir(parents=True, exist_ok=False)
    path = args.output / "summary.json"
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
                    encoding="utf-8")
    print(json.dumps({
        "status": payload["status"],
        "science_measurement": payload["science_measurement"],
        "source_reuse": reuse,
        "mode_p_values": {mode: [row["holm_p_value_7_radial_tests"]
                                  for row in result["tomographic_consistency"]]
                          for mode, result in modes.items()},
        "output": str(path),
    }, ensure_ascii=False), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
