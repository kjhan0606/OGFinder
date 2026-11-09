#!/usr/bin/env python3
"""Redshift-match a parent-sample richness split and compare HSC Y3 stacks.

Uses only the already-approved HSC source-row export. Matching is performed on
the parent lens manifest before loading or examining shear profiles. The result
is a cluster-only sensitivity diagnostic, not a mass measurement.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.optimize import linear_sum_assignment

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
import sys
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from fetch_ancillary_catalogues import sha
from ogfkit.wlcluster import Cosmology
from run_wl_cluster_spatial_jackknife import equal_area_block_ids, delete_block_covariance, _profile_vector
from run_wl_shape_stack import verify_and_load_pairs
from s19a_shape_stack import load_hsc_y3_nz_sacc, stack_s19a_tomographic_pairs


def make_matched_split(sample_manifest: Path, *, caliper: float):
    manifest = json.loads(sample_manifest.read_text(encoding="utf-8"))
    clusters = manifest.get("clusters")
    if not isinstance(clusters, list) or len(clusters) != 800:
        raise ValueError("the parent manifest must contain exactly 800 clusters")
    rows = []
    for item in clusters:
        lens_id = str(item["id"])
        z = float(item["z"])
        richness = float(item["richness"])
        if not np.isfinite(z) or not np.isfinite(richness):
            raise ValueError(f"non-finite parent metadata for lens {lens_id}")
        rows.append({"id": lens_id, "z": z, "richness": richness})
    if len({row["id"] for row in rows}) != len(rows):
        raise ValueError("parent lens IDs are not unique")
    order = sorted(range(800), key=lambda i: (rows[i]["richness"], rows[i]["id"]))
    low = [rows[i] for i in order[:400]]
    high = [rows[i] for i in order[400:]]
    z_low = np.asarray([row["z"] for row in low])
    z_high = np.asarray([row["z"] for row in high])
    distance = np.abs(z_low[:, None] - z_high[None, :])
    # One dummy column per low-richness lens permits unmatched objects. Since
    # every allowed match costs < caliper < 1 and a dummy costs 1, assignment
    # first maximizes the number of pairs and then minimizes total |delta z|.
    costs = np.ones((400, 800), dtype=float)
    costs[:, :400] = np.where(distance <= caliper, distance, 2.0)
    low_idx, columns = linear_sum_assignment(costs)
    matched = []
    for i, j in zip(low_idx, columns):
        if j < 400 and distance[i, j] <= caliper:
            matched.append({
                "pair_id": f"zmatch_{len(matched):04d}",
                "lower_lens_id": low[i]["id"],
                "upper_lens_id": high[j]["id"],
                "lower_z": low[i]["z"],
                "upper_z": high[j]["z"],
                "abs_delta_z": float(distance[i, j]),
                "lower_richness": low[i]["richness"],
                "upper_richness": high[j]["richness"],
            })
    if len(matched) < 3:
        raise ValueError("fewer than three redshift-matched pairs")
    return matched, {
        "method": "balanced parent richness split; maximum-cardinality 1:1 minimum-|delta-z| assignment",
        "split_rule": "sort all 800 parent lenses by (CAMIRA richness, lens_id); first 400 are lower-richness",
        "caliper_abs_delta_z": float(caliper),
        "n_parent_per_richness_group": 400,
        "n_matched_pairs": len(matched),
        "n_unmatched_per_group": 400 - len(matched),
        "mean_abs_delta_z": float(np.mean([x["abs_delta_z"] for x in matched])),
        "median_abs_delta_z": float(np.median([x["abs_delta_z"] for x in matched])),
        "max_abs_delta_z": float(np.max([x["abs_delta_z"] for x in matched])),
        "lower_richness_range": [float(min(x["lower_richness"] for x in matched)), float(max(x["lower_richness"] for x in matched))],
        "upper_richness_range": [float(min(x["upper_richness"] for x in matched)), float(max(x["upper_richness"] for x in matched))],
        "richness_difference_mean_upper_minus_lower": float(np.mean([x["upper_richness"] - x["lower_richness"] for x in matched])),
    }


def run(args):
    if args.output.exists():
        raise FileExistsError(f"refusing to overwrite {args.output}")
    matched, matching = make_matched_split(args.sample_manifest, caliper=args.caliper)
    pairs, provenance = verify_and_load_pairs(args.rows, args.preflight, args.plan_dir, args.masked_value)
    lens_ids = np.asarray(pairs["lens_id"]).astype(str)
    present = set(lens_ids)
    # A parent lens can have no exported rows after the approved shape-source
    # cone/count selection. Keep a matched pair only when both lenses have
    # rows, so the final stack remains balanced in lens count.
    supported = [row for row in matched
                 if row["lower_lens_id"] in present and row["upper_lens_id"] in present]
    if len(supported) < 3:
        raise ValueError("fewer than three matched pairs have source-row support on both lenses")
    pre_support_matching = dict(matching)
    matching["n_pairs_before_source_support_filter"] = len(matched)
    matching["n_pairs_without_complete_source_row_support"] = len(matched) - len(supported)
    matched = supported
    z_separations = [row["abs_delta_z"] for row in matched]
    matching.update({
        "n_matched_pairs": len(matched),
        "n_unmatched_per_group_after_support_filter": 400 - len(matched),
        "mean_abs_delta_z_after_support_filter": float(np.mean(z_separations)),
        "median_abs_delta_z_after_support_filter": float(np.median(z_separations)),
        "max_abs_delta_z_after_support_filter": float(np.max(z_separations)),
        "mean_richness_difference_upper_minus_lower_after_support_filter": float(np.mean([row["upper_richness"] - row["lower_richness"] for row in matched])),
        "pre_support_matching_summary": pre_support_matching,
    })
    low_ids = {row["lower_lens_id"] for row in matched}
    high_ids = {row["upper_lens_id"] for row in matched}
    low_mask = np.isin(lens_ids, list(low_ids))
    high_mask = np.isin(lens_ids, list(high_ids))
    if np.any(low_mask & high_mask):
        raise ValueError("a lens appears in both matched groups")
    nz_z, nz_bins = load_hsc_y3_nz_sacc(args.nz_sacc)
    cosmology = Cosmology(H0=args.h0, Om=args.omega_m)
    lens_cells = equal_area_block_ids(pairs["lens_ra_deg"], pairs["lens_dec_deg"], cell_width_deg=args.cell_width_deg)
    source_cells = equal_area_block_ids(pairs["source_ra_deg"], pairs["source_dec_deg"], cell_width_deg=args.cell_width_deg)
    block_ids = np.unique(np.concatenate((lens_cells[low_mask | high_mask], source_cells[low_mask | high_mask])))
    if len(block_ids) < 3:
        raise ValueError("at least three occupied spatial blocks are required")

    def subset(mask):
        return {key: np.asarray(value)[mask] for key, value in pairs.items()}

    cohort_masks = (low_mask, high_mask)
    names = ("lower_richness", "upper_richness")
    full_results = []
    full_vectors = []
    cohort_counts = []
    for name, mask in zip(names, cohort_masks):
        result = stack_s19a_tomographic_pairs(subset(mask), nz_z, nz_bins, args.radial_edges_mpc, cosmology)
        full_results.append(result)
        full_vectors.append(_profile_vector(result))
        used = np.asarray(pairs["lens_id"])[mask & (np.asarray(pairs["source_zbin"]).astype(int) > 0)]
        cohort_counts.append({
            "cohort": name,
            "matched_parent_lenses": len(matched),
            "input_pairs": int(np.count_nonzero(mask)),
            "unique_lenses_with_assigned_source_zbin": int(len(np.unique(used))),
            "stack_counts": result["counts"],
        })

    replicates = []
    for block_id in block_ids:
        vectors = []
        for mask in cohort_masks:
            keep = mask & (lens_cells != block_id) & (source_cells != block_id)
            result = stack_s19a_tomographic_pairs(subset(keep), nz_z, nz_bins, args.radial_edges_mpc, cosmology)
            vectors.append(_profile_vector(result))
        replicates.append(np.concatenate(vectors))
    covariance = delete_block_covariance(np.asarray(replicates))
    nbin = len(args.radial_edges_mpc) - 1
    transform = np.zeros((2 * nbin, 4 * nbin), dtype=float)
    for i in range(2 * nbin):
        transform[i, 2 * nbin + i] = 1.0
        transform[i, i] = -1.0
    contrast = full_vectors[1] - full_vectors[0]
    contrast_cov = transform @ covariance @ transform.T
    errors = np.sqrt(np.maximum(np.diag(contrast_cov), 0.0))
    profile = []
    for i, (lo, hi) in enumerate(zip(full_results[0]["profile"], full_results[1]["profile"])):
        profile.append({
            "bin": i, "r_min_mpc": lo["r_min_mpc"], "r_max_mpc": lo["r_max_mpc"],
            "lower_delta_sigma_msun_mpc2": lo["delta_sigma_weak_msun_mpc2"],
            "upper_delta_sigma_msun_mpc2": hi["delta_sigma_weak_msun_mpc2"],
            "upper_minus_lower_delta_sigma_msun_mpc2": float(contrast[2*i]),
            "upper_minus_lower_jackknife_error_msun_mpc2": float(errors[2*i]),
            "lower_cross_msun_mpc2": lo["cross_delta_sigma_weak_msun_mpc2"],
            "upper_cross_msun_mpc2": hi["cross_delta_sigma_weak_msun_mpc2"],
            "upper_minus_lower_cross_msun_mpc2": float(contrast[2*i+1]),
            "upper_minus_lower_cross_jackknife_error_msun_mpc2": float(errors[2*i+1]),
            "n_lower_pairs": lo["n_pairs"], "n_upper_pairs": hi["n_pairs"],
        })
    payload = {
        "schema": "astrafek-hsc-y3-redshift-matched-richness-diagnostic-v1",
        "status": "exploratory", "science_measurement": False,
        "analysis_script_sha256": sha(Path(__file__).resolve()),
        "random_lens_subtraction_applied": False,
        "photoz_uncertainty_propagated": False,
        "reduced_shear_applied": False,
        "matching": matching,
        "matched_lens_pairs": matched,
        "method": "redshift-matched cluster-only HSC Y3 tomographic stacks; shared delete-one equal-area sky-cell jackknife; all per-bin calibrations recomputed",
        "cell_width_equatorial_deg": float(args.cell_width_deg),
        "n_occupied_cells": int(len(block_ids)),
        "cohorts": cohort_counts,
        "profile": profile,
        "contrast_covariance_order": "[DeltaSigma bin 0, cross bin 0, DeltaSigma bin 1, cross bin 1, ...]",
        "upper_minus_lower_covariance": contrast_cov.tolist(),
        "input": {
            **provenance,
            "sample_manifest": str(args.sample_manifest.resolve()), "sample_manifest_sha256": sha(args.sample_manifest),
            "n_z_sacc": str(args.nz_sacc.resolve()), "n_z_sacc_sha256": sha(args.nz_sacc),
            "radial_edges_mpc": [float(value) for value in args.radial_edges_mpc],
            "cosmology": {"H0_km_s_mpc": args.h0, "Omega_m": args.omega_m},
            "b_mode_masked_raw_value": args.masked_value,
        },
        "limitations": [
            "cluster-only profile; matched random-lens subtraction and matched-field variance are absent",
            "equal-area blocks are not validated contiguous HSC survey regions",
            "matching controls lens redshift only; remaining selection, richness, and source-support differences remain",
            "fiducial n(z), shear calibration, source-selection, and reduced-shear uncertainties are not propagated",
            "the contrast is a robustness diagnostic, not a richness-mass calibration or detection significance",
        ],
    }
    args.output.mkdir(parents=True, exist_ok=False)
    target = args.output / "summary.json"
    with target.open("x", encoding="utf-8") as stream:
        json.dump(payload, stream, indent=2, ensure_ascii=False, allow_nan=False)
        stream.write("\n")
    print(json.dumps({"status": "exploratory", "science_measurement": False,
                      "matched_pairs": len(matched), "occupied_cells": int(len(block_ids)),
                      "low_used_pairs": cohort_counts[0]["stack_counts"]["used_pairs"],
                      "high_used_pairs": cohort_counts[1]["stack_counts"]["used_pairs"],
                      "output": str(target)}, ensure_ascii=False))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rows", type=Path, required=True)
    parser.add_argument("--preflight", type=Path, required=True)
    parser.add_argument("--plan-dir", type=Path, required=True)
    parser.add_argument("--sample-manifest", type=Path, required=True)
    parser.add_argument("--nz-sacc", type=Path, required=True)
    parser.add_argument("--masked-value", type=int, choices=(0, 1), required=True)
    parser.add_argument("--caliper", type=float, default=0.02)
    parser.add_argument("--cell-width-deg", type=float, default=5.0)
    parser.add_argument("--radial-edges-mpc", type=float, nargs="+", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--h0", type=float, default=70.0)
    parser.add_argument("--omega-m", type=float, default=0.3)
    run(parser.parse_args())


if __name__ == "__main__":
    main()
