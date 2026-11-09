#!/usr/bin/env python3
"""Measure matched lens/random weak-lensing sensitivity to HSC tract regions.

This is a leave-one-connected-tract-center-component-out influence analysis.
Patch centers define an approximate footprint graph; they are not polygon
boundaries or official field labels, so the result is not a field covariance.
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
PROJECT_ANALYSIS = HERE.parents[2]
if str(PROJECT_ANALYSIS) not in sys.path:
    sys.path.insert(0, str(PROJECT_ANALYSIS))

from fetch_ancillary_catalogues import sha
from ogfkit.wlcluster import Cosmology
from run_wl_hsc_tract_region_influence import (
    nearest_region_labels,
    tract_centers_from_patch_qa,
    tract_region_labels,
)
from run_wl_random_shape_stack import (
    verify_and_load_pairs as verify_and_load_random_pairs,
    verify_saved_cluster_profile,
    validate_cluster_stack_compatibility,
)
from run_wl_shape_stack import verify_and_load_pairs as verify_and_load_lens_pairs
from s19a_shape_stack import (
    load_hsc_y3_nz_sacc,
    stack_s19a_tomographic_pairs,
    subtract_random_lens_profile,
)


ROOT = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
DEFAULT_LENS_ROWS = ROOT / "wl_shape_rows_export_2026-10-18_03"
DEFAULT_LENS_PREFLIGHT = ROOT / "wl_shape_preflight_2026-10-18_01"
DEFAULT_LENS_PLAN = ROOT / "wl_shape_rows_plan_2026-10-18_04"
DEFAULT_RANDOM_ROWS = ROOT / "wl_random_shape_rows_lens_reuse_2026-11-03_01/results"
DEFAULT_RANDOM_PLAN = ROOT / "wl_random_shape_rows_lens_reuse_2026-11-03_01"
DEFAULT_LENS_STACK = ROOT / "wl_cluster_stack_bmode0_2026-10-19_01/summary.json"
DEFAULT_NZ = ROOT / "hsc_y3_nz_sacc_2026-10-19_01/hsc_y3_fourier_space_data_vector.sacc"
DEFAULT_PATCH_QA = ROOT / "catalogues_2026-10-13_01/patch_qa.csv"
DEFAULT_OUTPUT = ROOT / "wl_random_tract_region_influence_2026-11-06_01"
RADIAL_EDGES_MPC = [0.1, 0.2, 0.4, 0.7, 1.0, 1.5, 2.0, 3.0]


def run(args) -> dict:
    output = args.output.resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite output directory: {output}")
    output.mkdir(parents=True)

    random_pairs, random_provenance = verify_and_load_random_pairs(
        args.random_rows, args.random_plan, args.masked_value,
    )
    lens_pairs, lens_provenance = verify_and_load_lens_pairs(
        args.lens_rows, args.lens_preflight, args.lens_plan, args.masked_value,
    )
    lens_payload = json.loads(args.lens_stack.read_text(encoding="utf-8"))
    nz_redshift, nz_by_bin = load_hsc_y3_nz_sacc(args.nz_sacc)
    radial_edges = RADIAL_EDGES_MPC
    cosmology = Cosmology(H0=args.h0, Om=args.omega_m)
    validate_cluster_stack_compatibility(
        lens_payload, sha(args.nz_sacc), args.masked_value, radial_edges,
        args.h0, args.omega_m,
    )

    # Confirm that the current verified shape export still reproduces the
    # saved 800-lens cluster stack before restricting to random common support.
    full_lens_stack = stack_s19a_tomographic_pairs(
        lens_pairs, nz_redshift, nz_by_bin, radial_edges, cosmology,
    )
    verify_saved_cluster_profile(lens_payload["stack"], full_lens_stack["profile"])

    centers_to_lens = random_provenance["random_center_matched_lens_ids"]
    center_counts = random_provenance["random_center_expected_source_counts"]
    selected_lenses = set(random_provenance["matched_lens_ids"])
    cluster_lenses = set(map(str, np.unique(lens_pairs["lens_id"])))
    positive_random_lenses = {
        centers_to_lens[center_id]
        for center_id, count in center_counts.items() if int(count) > 0
    }
    common_lenses = selected_lenses & cluster_lenses & positive_random_lenses
    if len(common_lenses) != 435:
        raise ValueError(f"expected 435 common lens IDs, found {len(common_lenses)}")
    lens_keep = np.isin(lens_pairs["lens_id"].astype(str), list(common_lenses))
    lens = {key: value[lens_keep] for key, value in lens_pairs.items()}
    center_to_lens = random_provenance["random_center_matched_lens_ids"]
    random_centers = {
        center_id for center_id, lens_id in center_to_lens.items()
        if lens_id in common_lenses
    }
    random_keep = np.isin(random_pairs["lens_id"].astype(str), list(random_centers))
    random = {key: value[random_keep] for key, value in random_pairs.items()}
    if len(random_centers) != 5 * len(common_lenses):
        raise ValueError("common-support lenses do not have exactly five random centers each")

    # Build the same approximate connected HSC tract-center components used by
    # the existing Astrafek influence audit, then assign both stack endpoints.
    tract_keys, tract_xyz = tract_centers_from_patch_qa(args.patch_qa)
    tract_regions = tract_region_labels(tract_xyz, args.adjacency_deg)
    lens_ids, first_rows = np.unique(lens["lens_id"].astype(str), return_index=True)
    lens_region_lookup = dict(zip(
        lens_ids,
        nearest_region_labels(
            lens["lens_ra_deg"][first_rows], lens["lens_dec_deg"][first_rows],
            tract_xyz, tract_regions,
        ).tolist(),
    ))
    lens_region = np.fromiter(
        (lens_region_lookup[str(value)] for value in lens["lens_id"]),
        dtype=np.int32, count=len(lens["lens_id"]),
    )
    random_center_region = nearest_region_labels(
        random["lens_ra_deg"], random["lens_dec_deg"], tract_xyz, tract_regions,
    )
    lens_source_region = nearest_region_labels(
        lens["source_ra_deg"], lens["source_dec_deg"], tract_xyz, tract_regions,
    )
    random_source_region = nearest_region_labels(
        random["source_ra_deg"], random["source_dec_deg"], tract_xyz, tract_regions,
    )
    occupied = np.unique(np.concatenate((lens_region, random_center_region,
                                         lens_source_region, random_source_region)))
    if len(occupied) < 2:
        raise ValueError("fewer than two connected tract-center components overlap the stack")

    lens_full = stack_s19a_tomographic_pairs(
        lens, nz_redshift, nz_by_bin, radial_edges, cosmology,
    )
    random_full = stack_s19a_tomographic_pairs(
        random, nz_redshift, nz_by_bin, radial_edges, cosmology,
    )
    full = subtract_random_lens_profile(lens_full, random_full)["profile"]
    if len(full) != 7 or any(row["delta_sigma_random_subtracted_msun_mpc2"] is None
                             for row in full):
        raise ValueError("full-support random-subtracted profile is incomplete")

    deletions = []
    for region in occupied:
        keep_lens = (lens_region != region) & (lens_source_region != region)
        keep_random = (random_center_region != region) & (random_source_region != region)
        lens_replica = stack_s19a_tomographic_pairs(
            {key: value[keep_lens] for key, value in lens.items()},
            nz_redshift, nz_by_bin, radial_edges, cosmology,
        )
        random_replica = stack_s19a_tomographic_pairs(
            {key: value[keep_random] for key, value in random.items()},
            nz_redshift, nz_by_bin, radial_edges, cosmology,
        )
        profile = subtract_random_lens_profile(lens_replica, random_replica)["profile"]
        rows = []
        for base, replica in zip(full, profile):
            t_key = "delta_sigma_random_subtracted_msun_mpc2"
            x_key = "cross_delta_sigma_random_subtracted_msun_mpc2"
            t_full, t_loo = base[t_key], replica[t_key]
            x_full, x_loo = base[x_key], replica[x_key]
            rows.append({
                "bin": base["bin"], "r_min_mpc": base["r_min_mpc"],
                "r_max_mpc": base["r_max_mpc"],
                "delta_sigma_full": t_full, "delta_sigma_leave_region_out": t_loo,
                "delta_sigma_shift": t_loo - t_full,
                "delta_sigma_fractional_shift": (t_loo - t_full) / t_full if t_full else None,
                "cross_full": x_full, "cross_leave_region_out": x_loo,
                "cross_shift": x_loo - x_full,
                "lens_pairs_remaining": replica["n_lens_pairs"],
                "random_pairs_remaining": replica["n_random_pairs"],
            })
        deleted_lenses = sorted(
            {str(value) for value in lens["lens_id"][lens_region == region]}
        )
        deleted_centers = sorted(
            {str(value) for value in random["lens_id"][random_center_region == region]}
        )
        deleted_sources = (
            set(map(str, lens["object_id"][lens_source_region == region]))
            | set(map(str, random["object_id"][random_source_region == region]))
        )
        deletions.append({
            "region_id": int(region),
            "deleted_lens_count": len(deleted_lenses),
            "deleted_lens_ids": deleted_lenses,
            "deleted_random_center_count": len(deleted_centers),
            "deleted_source_ids_count": len(deleted_sources),
            "deleted_lens_pair_rows": int(np.count_nonzero(~keep_lens)),
            "deleted_random_pair_rows": int(np.count_nonzero(~keep_random)),
            "profile": rows,
        })

    result = {
        "schema": "astrafek-hsc-s19a-random-subtracted-tract-region-influence-v1",
        "status": "exploratory_field_influence",
        "science_measurement": False,
        "random_lens_subtraction_applied": True,
        "covariance_estimated": False,
        "method": "delete each occupied connected HSC tract-center component; on both lens and random stacks remove every pair whose lens/random center or source is assigned to that component; recompute complete calibrated estimators",
        "tract_assignment": "nearest spherical-mean patch-QA tract center, then connected component of tract centers",
        "tract_adjacency_deg": float(args.adjacency_deg),
        "n_tracts": len(tract_keys),
        "n_tract_components": int(len(np.unique(tract_regions))),
        "n_occupied_stack_components": int(len(occupied)),
        "common_support_lenses": len(common_lenses),
        "common_support_random_centers": len(random_centers),
        "lens_pair_rows": int(len(lens["lens_id"])),
        "random_pair_rows": int(len(random["lens_id"])),
        "full_profile": full,
        "region_deletions": deletions,
        "inputs": {
            "lens_row_plan_sha256": lens_provenance["row_plan_sha256"],
            "lens_pair_rows": lens_provenance["n_pair_rows"],
            "random_row_plan_sha256": random_provenance["random_row_plan_sha256"],
            "random_count_plan_sha256": random_provenance["count_plan_sha256"],
            "n_random_source_pairs_full_export": random_provenance["n_source_pairs"],
            "lens_stack_sha256": sha(args.lens_stack),
            "n_z_sacc_sha256": sha(args.nz_sacc),
            "patch_qa_sha256": sha(args.patch_qa),
            "masked_value": args.masked_value,
            "radial_edges_mpc": radial_edges,
            "cosmology": {"H0_km_s_Mpc": args.h0, "Omega_m": args.omega_m},
        },
        "limitations": [
            "patch-center proximity is an approximate footprint graph, not polygon overlap or official survey field IDs",
            "few occupied components permit influence analysis only, not a stable covariance estimate",
            "random-center validity is inherited from the upstream random selection and matching audit",
            "photo-z calibration, reduced shear, magnification, shear calibration, PSF and selection uncertainties are not propagated",
        ],
    }
    with (output / "summary.json").open("x", encoding="utf-8") as stream:
        json.dump(result, stream, ensure_ascii=False, indent=2, allow_nan=False)
        stream.write("\n")
    print(json.dumps({
        "status": result["status"], "science_measurement": False,
        "common_support_lenses": len(common_lenses),
        "tract_components": result["n_tract_components"],
        "occupied_components": len(occupied),
        "lens_pair_rows": len(lens["lens_id"]),
        "random_pair_rows": len(random["lens_id"]),
        "output": str(output / "summary.json"),
    }, ensure_ascii=False))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lens-rows", type=Path, default=DEFAULT_LENS_ROWS)
    parser.add_argument("--lens-preflight", type=Path, default=DEFAULT_LENS_PREFLIGHT)
    parser.add_argument("--lens-plan", type=Path, default=DEFAULT_LENS_PLAN)
    parser.add_argument("--random-rows", type=Path, default=DEFAULT_RANDOM_ROWS)
    parser.add_argument("--random-plan", type=Path, default=DEFAULT_RANDOM_PLAN)
    parser.add_argument("--lens-stack", type=Path, default=DEFAULT_LENS_STACK)
    parser.add_argument("--nz-sacc", type=Path, default=DEFAULT_NZ)
    parser.add_argument("--patch-qa", type=Path, default=DEFAULT_PATCH_QA)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--masked-value", type=int, choices=(0, 1), default=0)
    parser.add_argument("--adjacency-deg", type=float, default=2.0)
    parser.add_argument("--h0", type=float, default=70.0)
    parser.add_argument("--omega-m", type=float, default=0.3)
    args = parser.parse_args()
    run(args)


if __name__ == "__main__":
    main()
