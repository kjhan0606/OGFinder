#!/usr/bin/env python3
"""Compare ICL-bright and ICL-faint supported samples in the random-subtracted HSC stack.

All HSC shape rows were previously approved and downloaded. This runner only
reads the verified local exports; it does not contact CAS or request data.
Within each filter, ICL-supported lenses in the exact random-subtracted WL
common sample are divided at the median target-minus-control 200--300 kpc
surface-brightness contrast. High and low ICL groups are matched one-to-one
in redshift and log richness without consulting shear profiles. Cluster and
five matched random centers per lens are stacked separately with a shared
spatial delete-one jackknife. This is an exploratory conditional association,
not an ICL detection or halo-mass calibration.
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
ICL_SUPPORT = DATA / "icl_annulus_support_2026-10-18_01/pair_metrics_k2p5_c40_r50_300.tsv"
ICL_WL_OVERLAP = DATA / "icl_wl_sample_overlap_audit_2026-11-05_03/summary.json"
NZ_SACC = DATA / "hsc_y3_nz_sacc_2026-10-19_01/hsc_y3_fourier_space_data_vector.sacc"
DEFAULT_OUTPUT = DATA / "wl_random_subtracted_icl_split_2026-11-05_01"
RADIAL_EDGES_MPC = np.asarray([0.1, 0.2, 0.4, 0.7, 1.0, 1.5, 2.0, 3.0])
SOURCE_BINS = (2, 3, 4)
MAX_DELTA_Z = 0.03
MAX_ABS_DELTA_LOG_RICHNESS = 0.25


def make_support_matched_pairs(target_pairs, sample_manifest, eligible_lens_ids,
                               icl_values, max_delta_z=MAX_DELTA_Z,
                               max_abs_delta_log_richness=MAX_ABS_DELTA_LOG_RICHNESS):
    manifest = json.loads(sample_manifest.read_text(encoding="utf-8"))
    meta = {str(row["id"]): row for row in manifest["clusters"]}
    supported = sorted({str(value) for value in eligible_lens_ids})
    cutoff = float(np.median([icl_values[lens_id] for lens_id in supported]))
    groups = {
        0: [lens_id for lens_id in supported if icl_values[lens_id] <= cutoff],
        1: [lens_id for lens_id in supported if icl_values[lens_id] > cutoff],
    }
    low, high = groups[0], groups[1]
    z_low = np.asarray([float(meta[lens_id]["z"]) for lens_id in low])
    z_high = np.asarray([float(meta[lens_id]["z"]) for lens_id in high])
    richness_low = np.asarray([float(meta[lens_id]["richness"]) for lens_id in low])
    richness_high = np.asarray([float(meta[lens_id]["richness"]) for lens_id in high])
    distance_z = np.abs(z_low[:, None] - z_high[None, :])
    distance_log_richness = np.abs(
        np.log(richness_low[:, None]) - np.log(richness_high[None, :])
    )
    allowed = ((distance_z <= max_delta_z)
               & (distance_log_richness <= max_abs_delta_log_richness))
    # Allowed costs are <0.2 and a dummy costs 1, so assignment maximizes
    # cardinality first, then minimizes normalized redshift/richness distance.
    costs = np.ones((len(low), len(high) + len(low)), dtype=float)
    match_cost = 0.1 * (
        distance_z / max_delta_z
        + distance_log_richness / max_abs_delta_log_richness
    )
    costs[:, :len(high)] = np.where(allowed, match_cost, 2.0)
    low_idx, columns = linear_sum_assignment(costs)
    pairs = []
    for i, j in zip(low_idx, columns):
        if j < len(high) and allowed[i, j]:
            low_id, high_id = low[i], high[j]
            pairs.append({
                "pair_id": f"icl_wl_match_{len(pairs):04d}",
                "lower_lens_id": low_id,
                "upper_lens_id": high_id,
                "lower_z": float(z_low[i]),
                "upper_z": float(z_high[j]),
                "abs_delta_z": float(distance_z[i, j]),
                "abs_delta_log_richness": float(distance_log_richness[i, j]),
                "lower_richness": float(meta[low_id]["richness"]),
                "upper_richness": float(meta[high_id]["richness"]),
                "lower_icl_delta_sb_200_300": float(icl_values[low_id]),
                "upper_icl_delta_sb_200_300": float(icl_values[high_id]),
            })
    if len(pairs) < 10:
        raise ValueError(f"too few redshift/richness-matched common-support pairs: {len(pairs)}")
    return pairs, groups, {
        "split_rule": "filter-specific median split of target-minus-control ICL 200-300 kpc surface brightness among common-support lenses",
        "icl_median_cut": cutoff,
        "n_low_icl_support_lenses": len(low),
        "n_high_icl_support_lenses": len(high),
        "n_usable_common_support_lenses_before_match": len(supported),
        "n_lower_support_lenses": len(low),
        "n_upper_support_lenses": len(high),
        "max_abs_delta_z": float(max_delta_z),
        "max_abs_delta_log_richness": float(max_abs_delta_log_richness),
        "n_redshift_matched_pairs": len(pairs),
        "n_unmatched_lower_support": len(low) - len(pairs),
        "n_unmatched_upper_support": len(high) - len(pairs),
        "mean_abs_delta_z": float(np.mean([x["abs_delta_z"] for x in pairs])),
        "median_abs_delta_z": float(np.median([x["abs_delta_z"] for x in pairs])),
        "max_observed_abs_delta_z": float(max(x["abs_delta_z"] for x in pairs)),
        "lower_icl_delta_median": float(np.median([x["lower_icl_delta_sb_200_300"] for x in pairs])),
        "upper_icl_delta_median": float(np.median([x["upper_icl_delta_sb_200_300"] for x in pairs])),
        "lower_richness_median": float(np.median([x["lower_richness"] for x in pairs])),
        "upper_richness_median": float(np.median([x["upper_richness"] for x in pairs])),
        "matching": "maximum-cardinality 1:1 assignment within z and log-richness calipers; normalized distances minimize total match cost",
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


def run(output: Path, filter_name="i2", block_width_deg=5.0,
        max_delta_z=MAX_DELTA_Z,
        max_abs_delta_log_richness=MAX_ABS_DELTA_LOG_RICHNESS):
    output = Path(output).resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite ICL-split result: {output}")
    if filter_name not in {"i", "i2"}:
        raise ValueError("filter must be i or i2")
    if sha(ICL_SUPPORT) != "86e6d453df847ba4c9316c4fdc188565b1bc292bb41b0c5ba8ad310864f02aa5":
        raise ValueError("ICL support table hash changed")
    if sha(ICL_WL_OVERLAP) != "98d94834b895f1275337f494a801aab7e021048ef502d4f7074635bb26b25e12":
        raise ValueError("ICL/WL overlap summary hash changed")
    overlap_doc = json.loads(ICL_WL_OVERLAP.read_text(encoding="utf-8"))
    if overlap_doc.get("status") != "complete_exploratory" or overlap_doc.get("science_measurement") is not False:
        raise ValueError("ICL/WL overlap summary is not a closed exploratory audit")
    filter_doc = overlap_doc["icl_profile_support_by_filter"][filter_name]
    if filter_doc.get("n_in_verified_random_subtracted_common_support", 0) < 10:
        raise ValueError("too few ICL/WL common-support lenses in this filter")
    icl_values = {}
    import csv
    with ICL_SUPPORT.open(encoding="utf-8-sig", newline="") as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            if row["filter"] != filter_name or not row["delta_sb_200_300"]:
                continue
            if float(row["coverage_200_300"]) < 0.40:
                continue
            lens_id = str(row["id"])
            if lens_id in icl_values:
                raise ValueError(f"duplicate ICL support row: {filter_name}/{lens_id}")
            icl_values[lens_id] = float(row["delta_sb_200_300"])
    expected_parent = int(overlap_doc["icl_profile_support_by_filter"][filter_name]["n_icl_profile_supported_clusters"])
    if len(icl_values) != expected_parent:
        raise ValueError(f"ICL support row count mismatch: expected {expected_parent}, found {len(icl_values)}")
    approved_common_ids = set(map(str, filter_doc["icl_supported_ids_in_common_wl_support"]))
    if not approved_common_ids <= set(icl_values):
        raise ValueError("ICL/WL overlap IDs do not map into supported ICL measurements")
    target_pairs, random_doc, _tomo = load_pairs_and_provenance()
    support_ids = {str(value) for value in np.unique(target_pairs["lens_id"])}
    if sha(RANDOM_SUMMARY) != "83c5484454033b2652adcba5bd0699f28082fa6983eefa960613da75491bafc8":
        raise ValueError("random-subtracted stack summary hash changed")
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
    target_pairs = subset(target_pairs, np.isin(target_pairs["source_zbin"], SOURCE_BINS))
    random_pairs = subset(random_pairs, np.isin(random_pairs["source_zbin"], SOURCE_BINS))
    nz_z, nz_bins = load_hsc_y3_nz_sacc(NZ_SACC)
    cosmology = Cosmology(H0=70.0, Om=0.3)
    active_target_ids = active_lens_ids(target_pairs, nz_z, nz_bins, cosmology)
    active_random_centers = active_lens_ids(random_pairs, nz_z, nz_bins, cosmology)
    active_random_lenses = {
        center_to_lens[center_id] for center_id in active_random_centers
        if center_id in center_to_lens
    }
    eligible_ids = (support_ids & active_target_ids & active_random_lenses
                    & approved_common_ids)
    if len(eligible_ids) < 20:
        raise ValueError(f"too few ICL/WL pairs with usable tomography: {len(eligible_ids)}")
    pairs, groups, matching_doc = make_support_matched_pairs(
        target_pairs, SAMPLE_MANIFEST, eligible_ids, icl_values,
        max_delta_z=max_delta_z,
        max_abs_delta_log_richness=max_abs_delta_log_richness,
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
        raise ValueError("a lens appears in both matched ICL cohorts")
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
    source_bin_centers = set(np.unique(random_pairs["lens_id"]).astype(str))
    expected_nonempty_centers = center_ids & source_bin_centers
    if random_selected_ids != expected_nonempty_centers:
        raise ValueError("verified random rows do not match expected positive-support centers")
    lens_labels = {lens_id: 0 for lens_id in lower_ids}
    lens_labels.update({lens_id: 1 for lens_id in upper_ids})
    random_center_labels = {center: lens_labels[center_to_lens[center]] for center in center_ids}
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
            "low_icl_delta_sigma_msun_mpc2": float(group_delta[0][2 * index]),
            "high_icl_delta_sigma_msun_mpc2": float(group_delta[1][2 * index]),
            "high_minus_low_delta_sigma_msun_mpc2": float(full_contrast[2 * index]),
            "high_minus_low_delta_sigma_error_msun_mpc2": float(errors[2 * index]),
            "high_minus_low_cross_msun_mpc2": float(full_contrast[2 * index + 1]),
            "high_minus_low_cross_error_msun_mpc2": float(errors[2 * index + 1]),
            "low_icl_target_lenses": int(full_target[0][0]["counts"]["unique_used_lenses"]),
            "high_icl_target_lenses": int(full_target[1][0]["counts"]["unique_used_lenses"]),
            "low_icl_target_source_pairs": int(full_target[0][0]["counts"]["used_pairs"]),
            "high_icl_target_source_pairs": int(full_target[1][0]["counts"]["used_pairs"]),
            "low_icl_random_source_pairs": int(full_random[0][0]["counts"]["used_pairs"]),
            "high_icl_random_source_pairs": int(full_random[1][0]["counts"]["used_pairs"]),
        })
    tangential = np.arange(0, len(full_contrast), 2)
    cov_t = covariance[np.ix_(tangential, tangential)]
    y_t = full_contrast[tangential]
    eigenvalues = np.linalg.eigvalsh(cov_t)
    if np.min(eigenvalues) <= 0:
        raise ValueError("ICL-contrast tangential covariance is not positive definite")
    statistic = float(y_t @ np.linalg.solve(cov_t, y_t))
    rank = int(len(y_t))
    cross_indices = np.arange(1, len(full_contrast), 2)
    cov_x = covariance[np.ix_(cross_indices, cross_indices)]
    y_x = full_contrast[cross_indices]
    if np.min(np.linalg.eigvalsh(cov_x)) <= 0:
        raise ValueError("ICL-contrast cross-component covariance is not positive definite")
    cross_statistic = float(y_x @ np.linalg.solve(cov_x, y_x))
    summary = {
        "schema": "astrafek-hsc-random-subtracted-icl-split-v1",
        "status": "complete_exploratory",
        "science_measurement": False,
        "hsc_network_used": False,
        "images_read": False,
        "new_data_requested": False,
        "analysis_script": str(Path(__file__).resolve()),
        "analysis_script_sha256": sha(Path(__file__).resolve()),
        "method": "filter-specific median split of measured target-minus-control 200-300 kpc ICL surface-brightness contrast; maximum-cardinality 1:1 matching within redshift and log-richness calipers; separate target and five-per-lens random stacks; shared delete-one sky-block covariance",
        "filter": filter_name,
        "source_tomographic_bins": list(SOURCE_BINS),
        "icl_metric": "delta_sb_200_300 from the K2.5/40% annulus-support table (target minus matched control; nJy arcsec^-2)",
        "icl_wl_overlap_summary_sha256": sha(ICL_WL_OVERLAP),
        "icl_support_table_sha256": sha(ICL_SUPPORT),
        "n_icl_wl_common_lenses_before_matching": len(eligible_ids),
        "icl_split_and_matching": matching_doc,
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
        "high_minus_low_covariance": covariance.tolist(),
        "conditional_tangential_wald": {
            "chi2": statistic,
            "rank": rank,
            "nominal_p": float(chi2.sf(statistic, rank)),
            "interpretation": "conditional diagnostic only; equal-area cells are approximate, and this is not a detection significance",
        },
        "conditional_cross_null": {
            "chi2": cross_statistic,
            "rank": int(len(y_x)),
            "nominal_p": float(chi2.sf(cross_statistic, len(y_x))),
            "interpretation": "systematics diagnostic under the same approximate block covariance; not a calibrated probability",
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
            "this is a selected, estimator-supported ICL subset and does not represent the full 800-cluster parent sample",
            "matching on redshift and log richness does not control all mass, environment, centering, or image-selection differences",
            "the ICL aperture signal is not corrected for physical injected-source Poisson variance, extended PSF wings, sky transfer, or completeness",
            "five random apertures per lens share approximate equal-area cells; covariance is not based on validated contiguous HSC regions",
            "this conditional ICL/WL contrast is exploratory and is not a confirmed ICL-lensing relation",
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
        "cross_null": summary["conditional_cross_null"],
        "summary_sha256": sha(path),
    }, indent=2))
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--filter", choices=("i", "i2"), default="i2")
    parser.add_argument("--block-width-deg", type=float, default=5.0)
    parser.add_argument("--max-delta-z", type=float, default=MAX_DELTA_Z)
    parser.add_argument("--max-abs-delta-log-richness", type=float,
                        default=MAX_ABS_DELTA_LOG_RICHNESS)
    args = parser.parse_args()
    run(args.output, args.filter, args.block_width_deg, args.max_delta_z,
        args.max_abs_delta_log_richness)


if __name__ == "__main__":
    main()
