#!/usr/bin/env python3
"""Compare redshift-matched richness groups in the random-subtracted HSC stack.

All HSC shape rows were previously approved and downloaded. This runner only
reads the verified local exports; it does not contact CAS or request data.
The parent 800-lens richness split is fixed first, then common random-support
lenses are paired one-to-one in redshift without consulting shear profiles.
Cluster and five matched random centers per lens are stacked separately with a
shared spatial delete-one jackknife.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
from scipy.optimize import linear_sum_assignment
from scipy.stats import chi2

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
ROOT = HERE.parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from audit_wl_richness_split import richness_labels
from fetch_ancillary_catalogues import sha
from fit_wl_random_stack_nfw_mass import DATA, load_pairs_and_provenance
from ogfkit.wlcluster import Cosmology
from run_wl_cluster_spatial_jackknife import (
    delete_block_covariance, equal_area_block_ids, _profile_vector,
)
from run_wl_random_shape_stack import verify_and_load_pairs as verify_random_pairs
from s19a_shape_stack import (
    load_hsc_y3_nz_sacc, pair_geometry,
    stack_s19a_tomographic_pairs, tomographic_nz_lensing_efficiency,
)


RANDOM_PLAN = DATA / "wl_random_shape_rows_lens_reuse_2026-11-03_01"
RANDOM_ROWS = RANDOM_PLAN / "results"
RANDOM_SUMMARY = DATA / "wl_random_shape_stack_lens_reuse_2026-11-04_01/summary.json"
SAMPLE_MANIFEST = DATA / "sample_manifest.json"
NZ_SACC = DATA / "hsc_y3_nz_sacc_2026-10-19_01/hsc_y3_fourier_space_data_vector.sacc"
DEFAULT_OUTPUT = DATA / "wl_random_subtracted_richness_split_2026-11-05_02"
RADIAL_EDGES_MPC = np.asarray([0.1, 0.2, 0.4, 0.7, 1.0, 1.5, 2.0, 3.0])
MAX_DELTA_Z = 0.02


def make_support_matched_pairs(target_pairs, parent_labels, sample_manifest,
                               eligible_lens_ids, max_delta_z=MAX_DELTA_Z):
    manifest = json.loads(sample_manifest.read_text(encoding="utf-8"))
    meta = {str(row["id"]): row for row in manifest["clusters"]}
    supported = sorted({str(value) for value in eligible_lens_ids})
    groups = {
        0: [lens_id for lens_id in supported if parent_labels[lens_id] == 0],
        1: [lens_id for lens_id in supported if parent_labels[lens_id] == 1],
    }
    low, high = groups[0], groups[1]
    z_low = np.asarray([float(meta[lens_id]["z"]) for lens_id in low])
    z_high = np.asarray([float(meta[lens_id]["z"]) for lens_id in high])
    distance = np.abs(z_low[:, None] - z_high[None, :])
    # A dummy assignment costs 1; allowed matches cost <0.02 and forbidden
    # matches cost 2. Thus assignment maximizes pair count before minimizing dz.
    costs = np.ones((len(low), len(high) + len(low)), dtype=float)
    costs[:, :len(high)] = np.where(distance <= max_delta_z, distance, 2.0)
    low_idx, columns = linear_sum_assignment(costs)
    pairs = []
    for i, j in zip(low_idx, columns):
        if j < len(high) and distance[i, j] <= max_delta_z:
            low_id, high_id = low[i], high[j]
            pairs.append({
                "pair_id": f"support_zmatch_{len(pairs):04d}",
                "lower_lens_id": low_id,
                "upper_lens_id": high_id,
                "lower_z": float(z_low[i]),
                "upper_z": float(z_high[j]),
                "abs_delta_z": float(distance[i, j]),
                "lower_richness": float(meta[low_id]["richness"]),
                "upper_richness": float(meta[high_id]["richness"]),
            })
    if len(pairs) < 10:
        raise ValueError(f"too few redshift-matched common-support pairs: {len(pairs)}")
    return pairs, groups, {
        "split_rule": "retain the predeclared lower/upper richness labels from the full parent-800 split",
        "n_usable_common_support_lenses_before_match": len(supported),
        "n_lower_support_lenses": len(low),
        "n_upper_support_lenses": len(high),
        "max_abs_delta_z": float(max_delta_z),
        "n_redshift_matched_pairs": len(pairs),
        "n_unmatched_lower_support": len(low) - len(pairs),
        "n_unmatched_upper_support": len(high) - len(pairs),
        "mean_abs_delta_z": float(np.mean([x["abs_delta_z"] for x in pairs])),
        "median_abs_delta_z": float(np.median([x["abs_delta_z"] for x in pairs])),
        "max_observed_abs_delta_z": float(max(x["abs_delta_z"] for x in pairs)),
        "lower_richness_median": float(np.median([x["lower_richness"] for x in pairs])),
        "upper_richness_median": float(np.median([x["upper_richness"] for x in pairs])),
        "matching": "maximum-cardinality 1:1 assignment within the common-support sample, minimizing total |delta-z|",
    }


def subset(pairs, mask):
    return {key: np.asarray(value)[mask] for key, value in pairs.items()}


def stack_vector(pairs, nz_z, nz_bins, cosmology):
    result = stack_s19a_tomographic_pairs(
        pairs, nz_z, nz_bins, RADIAL_EDGES_MPC, cosmology,
    )
    return result, _profile_vector(result)


def active_lens_ids(pairs, nz_z, nz_bins, cosmology):
    """Mirror the stack's source, B-mode, tomography and radial support cuts."""
    source_bin = np.asarray(pairs["source_zbin"], dtype=int)
    beta = tomographic_nz_lensing_efficiency(
        pairs["lens_z"], source_bin, nz_z, nz_bins, cosmology,
    )
    radius, _phi = pair_geometry(pairs, cosmology)
    shape = np.hypot(pairs["e1"], pairs["e2"])
    usable = (
        np.asarray(pairs["b_mode_pass"], dtype=bool)
        & (pairs["sigma_e"] >= 0.0) & (pairs["sigma_e"] <= 0.4)
        & (shape < 2.0) & (pairs["mag_a10"] <= 25.5)
        & (pairs["resolution"] >= 0.3) & (beta > 0.0)
        & (radius >= RADIAL_EDGES_MPC[0])
        & (radius <= RADIAL_EDGES_MPC[-1])
    )
    return {str(value) for value in np.unique(np.asarray(pairs["lens_id"])[usable])}


def run(output: Path, block_width_deg=5.0):
    output = Path(output).resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite richness-split result: {output}")
    parent_labels, split_doc = richness_labels(SAMPLE_MANIFEST)
    target_pairs, random_doc, _tomo = load_pairs_and_provenance()
    if sha(RANDOM_SUMMARY) != "83c5484454033b2652adcba5bd0699f28082fa6983eefa960613da75491bafc8":
        raise ValueError("random-subtracted stack summary hash changed")
    support_ids = {str(value) for value in np.unique(target_pairs["lens_id"])}
    center_to_lens = random_doc["input"]["random_center_matched_lens_ids"]
    center_counts = random_doc["input"]["random_center_expected_source_counts"]
    positive_random_lenses = {
        center_to_lens[center_id]
        for center_id, count in center_counts.items() if int(count) > 0
    }
    if support_ids != support_ids & positive_random_lenses:
        raise ValueError("target source support and random-center positive support differ")
    common_hash = hashlib.sha256("\n".join(sorted(support_ids)).encode()).hexdigest()
    if common_hash != random_doc["input"]["matched_sample"]["common_lens_ids_sha256"]:
        raise ValueError("reconstructed common-support lens IDs do not match the saved stack")

    random_pairs, random_provenance = verify_random_pairs(
        RANDOM_ROWS, RANDOM_PLAN, masked_value=0,
    )
    nz_z, nz_bins = load_hsc_y3_nz_sacc(NZ_SACC)
    cosmology = Cosmology(H0=70.0, Om=0.3)
    active_target_ids = active_lens_ids(target_pairs, nz_z, nz_bins, cosmology)
    active_random_centers = active_lens_ids(random_pairs, nz_z, nz_bins, cosmology)
    active_random_lenses = {
        center_to_lens[center_id] for center_id in active_random_centers
        if center_id in center_to_lens
    }
    eligible_ids = support_ids & active_target_ids & active_random_lenses
    if len(eligible_ids) < 20:
        raise ValueError(f"too few target/control pairs with usable tomography: {len(eligible_ids)}")
    pairs, groups, matching_doc = make_support_matched_pairs(
        target_pairs, parent_labels, SAMPLE_MANIFEST, eligible_ids,
    )
    matching_doc.update({
        "n_raw_common_support_lenses": len(support_ids),
        "n_active_target_lenses": len(active_target_ids),
        "n_lenses_with_at_least_one_active_random_center": len(active_random_lenses),
        "n_lenses_with_active_target_and_random_support": len(eligible_ids),
        "target_usable_lens_ids_sha256": hashlib.sha256(
            "\n".join(sorted(active_target_ids)).encode()
        ).hexdigest(),
        "random_usable_center_ids_sha256": hashlib.sha256(
            "\n".join(sorted(active_random_centers)).encode()
        ).hexdigest(),
    })
    lower_ids = {pair["lower_lens_id"] for pair in pairs}
    upper_ids = {pair["upper_lens_id"] for pair in pairs}
    if lower_ids & upper_ids:
        raise ValueError("a lens appears in both matched richness cohorts")
    chosen_ids = lower_ids | upper_ids
    target_id = np.asarray(target_pairs["lens_id"]).astype(str)
    target_lower = np.isin(target_id, list(lower_ids))
    target_upper = np.isin(target_id, list(upper_ids))
    target_cohorts = [subset(target_pairs, target_lower), subset(target_pairs, target_upper)]

    center_ids = {
        center for center, lens_id in center_to_lens.items() if lens_id in chosen_ids
    }
    centers_per_lens = {}
    for center in center_ids:
        lens_id = center_to_lens[center]
        centers_per_lens[lens_id] = centers_per_lens.get(lens_id, 0) + 1
    if set(centers_per_lens.values()) != {5}:
        raise ValueError("every matched lens must retain five approved random centers")
    random_id = np.asarray(random_pairs["lens_id"]).astype(str)
    random_keep = np.isin(random_id, list(center_ids))
    random_selected = subset(random_pairs, random_keep)
    random_selected_ids = set(np.unique(random_selected["lens_id"]).astype(str))
    expected_nonempty_centers = {
        center for center in center_ids if int(center_counts[center]) > 0
    }
    if random_selected_ids != expected_nonempty_centers:
        raise ValueError("verified random rows do not match expected positive-support centers")
    random_center_labels = {
        center: parent_labels[center_to_lens[center]] for center in center_ids
    }
    random_cohorts = [
        subset(random_selected, np.fromiter(
            (random_center_labels[value] == group
             for value in np.asarray(random_selected["lens_id"]).astype(str)),
            dtype=bool, count=len(random_selected["lens_id"]),
        ))
        for group in (0, 1)
    ]

    full_target = [stack_vector(cohort, nz_z, nz_bins, cosmology) for cohort in target_cohorts]
    full_random = [stack_vector(cohort, nz_z, nz_bins, cosmology) for cohort in random_cohorts]
    group_delta = [full_target[i][1] - full_random[i][1] for i in range(2)]
    full_contrast = group_delta[1] - group_delta[0]

    endpoint_pairs = [*target_cohorts, *random_cohorts]
    endpoint_blocks = []
    for endpoint in endpoint_pairs:
        lens_cells = equal_area_block_ids(
            endpoint["lens_ra_deg"], endpoint["lens_dec_deg"], cell_width_deg=block_width_deg,
        )
        source_cells = equal_area_block_ids(
            endpoint["source_ra_deg"], endpoint["source_dec_deg"], cell_width_deg=block_width_deg,
        )
        endpoint_blocks.append((lens_cells, source_cells))
    block_ids = np.unique(np.concatenate([
        value for pair in endpoint_blocks for value in pair
    ]))
    if len(block_ids) < 3:
        raise ValueError("at least three occupied shared sky blocks are required")

    replicates, deletions = [], []
    for block_id in block_ids:
        vectors = []
        counts_deleted = []
        for endpoint, (lens_cells, source_cells) in zip(endpoint_pairs, endpoint_blocks):
            keep = (lens_cells != block_id) & (source_cells != block_id)
            result, vector = stack_vector(subset(endpoint, keep), nz_z, nz_bins, cosmology)
            vectors.append(vector)
            counts_deleted.append(int(np.count_nonzero(~keep)))
        low = vectors[0] - vectors[2]
        high = vectors[1] - vectors[3]
        replicates.append(high - low)
        deletions.append({
            "block_id": int(block_id),
            "target_lower_rows_deleted": counts_deleted[0],
            "target_upper_rows_deleted": counts_deleted[1],
            "random_lower_rows_deleted": counts_deleted[2],
            "random_upper_rows_deleted": counts_deleted[3],
        })
    covariance = delete_block_covariance(np.asarray(replicates, dtype=float))
    errors = np.sqrt(np.maximum(np.diag(covariance), 0.0))
    profile = []
    for index in range(len(RADIAL_EDGES_MPC) - 1):
        lo, hi = full_target[0][0]["profile"][index], full_target[1][0]["profile"][index]
        profile.append({
            "bin": index,
            "r_min_mpc": float(lo["r_min_mpc"]),
            "r_max_mpc": float(lo["r_max_mpc"]),
            "lower_richness_delta_sigma_msun_mpc2": float(group_delta[0][2 * index]),
            "upper_richness_delta_sigma_msun_mpc2": float(group_delta[1][2 * index]),
            "upper_minus_lower_delta_sigma_msun_mpc2": float(full_contrast[2 * index]),
            "upper_minus_lower_delta_sigma_error_msun_mpc2": float(errors[2 * index]),
            "upper_minus_lower_cross_msun_mpc2": float(full_contrast[2 * index + 1]),
            "upper_minus_lower_cross_error_msun_mpc2": float(errors[2 * index + 1]),
            "lower_target_lenses": int(full_target[0][0]["counts"]["unique_used_lenses"]),
            "upper_target_lenses": int(full_target[1][0]["counts"]["unique_used_lenses"]),
            "lower_target_source_pairs": int(full_target[0][0]["counts"]["used_pairs"]),
            "upper_target_source_pairs": int(full_target[1][0]["counts"]["used_pairs"]),
            "lower_random_source_pairs": int(full_random[0][0]["counts"]["used_pairs"]),
            "upper_random_source_pairs": int(full_random[1][0]["counts"]["used_pairs"]),
        })
    tangential = np.arange(0, len(full_contrast), 2)
    cov_t = covariance[np.ix_(tangential, tangential)]
    y_t = full_contrast[tangential]
    eigenvalues = np.linalg.eigvalsh(cov_t)
    if np.min(eigenvalues) <= 0:
        raise ValueError("richness-contrast tangential covariance is not positive definite")
    statistic = float(y_t @ np.linalg.solve(cov_t, y_t))
    rank = int(len(y_t))
    summary = {
        "schema": "astrafek-hsc-random-subtracted-richness-split-v1",
        "status": "complete_exploratory",
        "science_measurement": False,
        "hsc_network_used": False,
        "images_read": False,
        "new_data_requested": False,
        "analysis_script": str(Path(__file__).resolve()),
        "analysis_script_sha256": sha(Path(__file__).resolve()),
        "method": "predeclared parent richness split; maximum-cardinality 1:1 redshift match within exact target/rand common support; separate target and five-per-lens random stacks; shared delete-one sky-block covariance",
        "parent_richness_split": split_doc,
        "support_matching": matching_doc,
        "n_matched_lens_pairs": len(pairs),
        "matched_lens_pairs": pairs,
        "n_random_centers_per_matched_lens": 5,
        "n_random_centers_in_selected_pair_sample": 10 * len(pairs),
        "n_random_centers_with_rows": len(random_selected_ids),
        "block_width_deg": float(block_width_deg),
        "n_shared_delete_one_blocks": int(len(block_ids)),
        "profile": profile,
        "contrast_covariance_order": "[DeltaSigma bin 0, cross bin 0, DeltaSigma bin 1, cross bin 1, ...]",
        "upper_minus_lower_covariance": covariance.tolist(),
        "conditional_tangential_wald": {
            "chi2": statistic,
            "rank": rank,
            "nominal_p": float(chi2.sf(statistic, rank)),
            "interpretation": "conditional diagnostic only; equal-area cells are approximate, and this is not a detection significance",
        },
        "provenance": {
            "random_stack_summary_sha256": sha(RANDOM_SUMMARY),
            "common_lens_ids_sha256": common_hash,
            "sample_manifest_sha256": sha(SAMPLE_MANIFEST),
            "hsc_y3_nz_sacc_sha256": sha(NZ_SACC),
            "random_row_plan_sha256": random_provenance["random_row_plan_sha256"],
            "random_count_plan_sha256": random_provenance["count_plan_sha256"],
            "random_batch_sha256": random_provenance["random_row_batch_sha256"],
            "target_stack_sha256": random_doc["input"]["target_stack_sha256"],
            "target_row_plan_sha256": random_doc["input"]["lens_row_plan_sha256"],
        },
        "limitations": [
            "only lenses with at least one estimator-usable target pair and one estimator-usable matched-random center enter the common sample; additional redshift matching further conditions this subset",
            "five random apertures per lens share approximate equal-area cells; the covariance is not based on validated contiguous HSC regions",
            "fiducial n(z), shear calibration, reduced shear, magnification, concentration/mass distribution, and random-aperture depth/mask/PSF/sky uncertainty are not propagated",
            "this random-subtracted richness contrast is exploratory and is not a richness-mass calibration or ICL-satellite detection",
        ],
    }
    output.mkdir(parents=True, exist_ok=False)
    path = output / "summary.json"
    path.write_text(json.dumps(summary, ensure_ascii=False, indent=2, allow_nan=False) + "\n")
    print(json.dumps({
        "status": summary["status"],
        "science_measurement": False,
        "matched_lens_pairs": len(pairs),
        "shared_blocks": len(block_ids),
        "wald": summary["conditional_tangential_wald"],
        "summary_sha256": sha(path),
    }, indent=2))
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--block-width-deg", type=float, default=5.0)
    args = parser.parse_args()
    run(args.output, args.block_width_deg)


if __name__ == "__main__":
    main()
