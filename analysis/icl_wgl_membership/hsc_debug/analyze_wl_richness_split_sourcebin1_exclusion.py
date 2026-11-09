#!/usr/bin/env python3
"""Repeat the richness-split random-subtracted stack excluding HSC source bin 1.

This sensitivity keeps the exact 204 redshift-matched target pairs and their
five preassigned random centers fixed, removes source_zbin=1 from both target
and random rows, and recomputes all estimator calibrations and the shared
spatial jackknife. It reads verified local exports only.
"""
from __future__ import annotations

import argparse
import hashlib
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
from fit_wl_random_stack_nfw_mass import DATA, load_pairs_and_provenance
from ogfkit.wlcluster import Cosmology
from run_wl_cluster_spatial_jackknife import (
    delete_block_covariance, equal_area_block_ids, _profile_vector,
)
from run_wl_random_shape_stack import verify_and_load_pairs as verify_random_pairs
from s19a_shape_stack import load_hsc_y3_nz_sacc, stack_s19a_tomographic_pairs


RICHNESS_SUMMARY = DATA / "wl_random_subtracted_richness_split_2026-11-05_02/summary.json"
EXPECTED_RICHNESS_SHA = "b38a78f6005dc17d8a12ddd6117debd63a005d0b21ff8dbda05961bd31eb65b7"
RANDOM_PLAN = DATA / "wl_random_shape_rows_lens_reuse_2026-11-03_01"
RANDOM_ROWS = RANDOM_PLAN / "results"
NZ_SACC = DATA / "hsc_y3_nz_sacc_2026-10-19_01/hsc_y3_fourier_space_data_vector.sacc"
DEFAULT_OUTPUT = DATA / "wl_random_subtracted_richness_split_zbin1_excluded_2026-11-05_01"
RADIAL_EDGES = np.asarray([0.1, 0.2, 0.4, 0.7, 1.0, 1.5, 2.0, 3.0])


def subset(pairs, mask):
    return {key: np.asarray(value)[mask] for key, value in pairs.items()}


def stack(pairs, nz_z, nz_bins, cosmology):
    result = stack_s19a_tomographic_pairs(
        pairs, nz_z, nz_bins, RADIAL_EDGES, cosmology,
    )
    return result, _profile_vector(result)


def run(output: Path, block_width_deg=5.0):
    output = Path(output).resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite source-bin sensitivity: {output}")
    richness = json.loads(RICHNESS_SUMMARY.read_text(encoding="utf-8"))
    if (sha(RICHNESS_SUMMARY) != EXPECTED_RICHNESS_SHA
            or richness.get("status") != "complete_exploratory"
            or richness.get("science_measurement") is not False
            or richness.get("n_matched_lens_pairs") != 204):
        raise ValueError("fixed richness-matched sample summary is changed or invalid")
    target_pairs, random_doc, _tomo = load_pairs_and_provenance()
    random_pairs, random_provenance = verify_random_pairs(
        RANDOM_ROWS, RANDOM_PLAN, masked_value=0,
    )
    low_ids = {str(row["lower_lens_id"]) for row in richness["matched_lens_pairs"]}
    high_ids = {str(row["upper_lens_id"]) for row in richness["matched_lens_pairs"]}
    selected_lenses = low_ids | high_ids
    centers = {
        center for center, lens_id in random_doc["input"]["random_center_matched_lens_ids"].items()
        if lens_id in selected_lenses
    }
    center_to_lens = random_doc["input"]["random_center_matched_lens_ids"]
    per_lens = {}
    for center in centers:
        lens_id = center_to_lens[center]
        per_lens[lens_id] = per_lens.get(lens_id, 0) + 1
    if set(per_lens.values()) != {5} or len(centers) != 2040:
        raise ValueError("the fixed matched sample must retain exactly five random centers per lens")
    target_id = np.asarray(target_pairs["lens_id"]).astype(str)
    control_id = np.asarray(random_pairs["lens_id"]).astype(str)
    target_groups = [subset(target_pairs, np.isin(target_id, list(ids)))
                     for ids in (low_ids, high_ids)]
    control_groups = [
        subset(random_pairs, np.fromiter(
            (center_to_lens[value] in ids for value in control_id),
            dtype=bool, count=len(control_id),
        ))
        for ids in (low_ids, high_ids)
    ]
    # Keep the lens sample and aperture assignments fixed; exclude source bin 1
    # from both halves of the random-subtracted estimator.
    target_groups = [subset(group, np.asarray(group["source_zbin"], dtype=int) >= 2)
                     for group in target_groups]
    control_groups = [subset(group, np.asarray(group["source_zbin"], dtype=int) >= 2)
                      for group in control_groups]
    nz_z, nz_bins = load_hsc_y3_nz_sacc(NZ_SACC)
    cosmology = Cosmology(H0=70.0, Om=0.3)
    target_full = [stack(group, nz_z, nz_bins, cosmology) for group in target_groups]
    control_full = [stack(group, nz_z, nz_bins, cosmology) for group in control_groups]
    group_delta = [target_full[i][1] - control_full[i][1] for i in range(2)]
    contrast = group_delta[1] - group_delta[0]
    endpoints = [*target_groups, *control_groups]
    block_data = []
    for group in endpoints:
        block_data.append((
            equal_area_block_ids(group["lens_ra_deg"], group["lens_dec_deg"],
                                 cell_width_deg=block_width_deg),
            equal_area_block_ids(group["source_ra_deg"], group["source_dec_deg"],
                                 cell_width_deg=block_width_deg),
        ))
    block_ids = np.unique(np.concatenate([x for pair in block_data for x in pair]))
    if len(block_ids) < 3:
        raise ValueError("at least three shared delete-one cells are required")
    replicates, deletions = [], []
    for block_id in block_ids:
        vectors, deleted = [], []
        for endpoint, (lens_cells, source_cells) in zip(endpoints, block_data):
            keep = (lens_cells != block_id) & (source_cells != block_id)
            result, vector = stack(subset(endpoint, keep), nz_z, nz_bins, cosmology)
            vectors.append(vector)
            deleted.append(int(np.count_nonzero(~keep)))
        low = vectors[0] - vectors[2]
        high = vectors[1] - vectors[3]
        replicates.append(high - low)
        deletions.append({
            "block_id": int(block_id),
            "target_lower_rows_deleted": deleted[0],
            "target_upper_rows_deleted": deleted[1],
            "random_lower_rows_deleted": deleted[2],
            "random_upper_rows_deleted": deleted[3],
        })
    covariance = delete_block_covariance(np.asarray(replicates, dtype=float))
    errors = np.sqrt(np.maximum(np.diag(covariance), 0.0))
    profile = []
    for index in range(len(RADIAL_EDGES) - 1):
        profile.append({
            "bin": index,
            "r_min_mpc": float(RADIAL_EDGES[index]),
            "r_max_mpc": float(RADIAL_EDGES[index + 1]),
            "lower_richness_random_subtracted_delta_sigma_msun_mpc2": float(group_delta[0][2 * index]),
            "upper_richness_random_subtracted_delta_sigma_msun_mpc2": float(group_delta[1][2 * index]),
            "upper_minus_lower_delta_sigma_msun_mpc2": float(contrast[2 * index]),
            "upper_minus_lower_delta_sigma_error_msun_mpc2": float(errors[2 * index]),
            "upper_minus_lower_cross_msun_mpc2": float(contrast[2 * index + 1]),
            "upper_minus_lower_cross_error_msun_mpc2": float(errors[2 * index + 1]),
            "lower_target_used_lenses": int(target_full[0][0]["counts"]["unique_used_lenses"]),
            "upper_target_used_lenses": int(target_full[1][0]["counts"]["unique_used_lenses"]),
            "lower_target_used_pairs": int(target_full[0][0]["counts"]["used_pairs"]),
            "upper_target_used_pairs": int(target_full[1][0]["counts"]["used_pairs"]),
            "lower_random_used_centers": int(control_full[0][0]["counts"]["unique_used_lenses"]),
            "upper_random_used_centers": int(control_full[1][0]["counts"]["unique_used_lenses"]),
        })
    tangent = np.arange(0, len(contrast), 2)
    cov_t = covariance[np.ix_(tangent, tangent)]
    eig = np.linalg.eigvalsh(cov_t)
    if eig.min() <= 0:
        raise ValueError("source-bin-excluded tangential contrast covariance is not positive definite")
    y = contrast[tangent]
    wald = float(y @ np.linalg.solve(cov_t, y))
    result = {
        "schema": "astrafek-wl-richness-split-sourcebin1-exclusion-v1",
        "status": "complete_exploratory_sensitivity",
        "science_measurement": False,
        "hsc_network_used": False,
        "images_read": False,
        "new_data_requested": False,
        "analysis_script": str(Path(__file__).resolve()),
        "analysis_script_sha256": sha(Path(__file__).resolve()),
        "fixed_sample": "same 204 redshift-matched richness pairs and five preassigned random centers per lens as wl_random_subtracted_richness_split_2026-11-05_02",
        "source_bins_used": [2, 3, 4],
        "source_bin_1_excluded_from_target_and_random": True,
        "block_width_deg": float(block_width_deg),
        "n_shared_delete_one_blocks": int(len(block_ids)),
        "n_lenses_per_richness_group": 204,
        "n_random_centers": 2040,
        "n_random_centers_with_rows": int(len(set(np.unique(random_pairs["lens_id"]).astype(str)) & centers)),
        "profile": profile,
        "contrast_covariance_order": "[DeltaSigma bin 0, cross bin 0, DeltaSigma bin 1, cross bin 1, ...]",
        "upper_minus_lower_covariance": covariance.tolist(),
        "conditional_tangential_wald": {
            "chi2": wald,
            "rank": int(len(y)),
            "nominal_p": float(chi2.sf(wald, len(y))),
            "interpretation": "conditional diagnostic only; equal-area cells are approximate, not detection significance",
        },
        "block_deletions": deletions,
        "provenance": {
            "all_source_bins_split_summary_sha256": sha(RICHNESS_SUMMARY),
            "random_row_plan_sha256": random_provenance["random_row_plan_sha256"],
            "random_count_plan_sha256": random_provenance["count_plan_sha256"],
            "random_batch_sha256": random_provenance["random_row_batch_sha256"],
            "random_stack_summary_sha256": sha(DATA / "wl_random_shape_stack_lens_reuse_2026-11-04_01/summary.json"),
            "hsc_y3_nz_sacc_sha256": sha(NZ_SACC),
        },
        "limitations": [
            "keeps the same selected 204-pair lens sample; it does not correct source-support selection",
            "matched-aperture depth/mask/PSF/sky equivalence and contiguous HSC covariance are unresolved",
            "fiducial n(z), shear calibration, reduced shear, magnification, richness-mass distribution, and covariance uncertainty are not propagated",
            "not a detection or science-grade mass calibration; science_measurement remains false",
        ],
    }
    output.mkdir(parents=True, exist_ok=False)
    path = output / "summary.json"
    path.write_text(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False) + "\n")
    print(json.dumps({
        "status": result["status"],
        "science_measurement": False,
        "source_bins_used": result["source_bins_used"],
        "shared_blocks": len(block_ids),
        "wald": result["conditional_tangential_wald"],
        "summary_sha256": sha(path),
    }, indent=2))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--block-width-deg", type=float, default=5.0)
    args = parser.parse_args()
    run(args.output, args.block_width_deg)


if __name__ == "__main__":
    main()
