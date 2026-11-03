#!/usr/bin/env python3
"""Compare cluster-only weak-lensing profiles across a predeclared richness split.

The split is fixed from the complete 800-lens parent manifest before inspecting
the shape profiles. Both samples use the same equal-area delete-one jackknife.
This is a population-heterogeneity diagnostic, not a mass measurement.
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

from fetch_ancillary_catalogues import sha
from ogfkit.wlcluster import Cosmology
from run_wl_cluster_spatial_jackknife import equal_area_block_ids, delete_block_covariance, _profile_vector
from run_wl_shape_stack import verify_and_load_pairs
from s19a_shape_stack import load_hsc_y3_nz_sacc, stack_s19a_tomographic_pairs


def richness_labels(sample_manifest: Path):
    manifest = json.loads(sample_manifest.read_text(encoding="utf-8"))
    clusters = manifest.get("clusters")
    if not isinstance(clusters, list) or len(clusters) != 800:
        raise ValueError("the parent manifest must contain exactly 800 clusters")
    rows = []
    for item in clusters:
        lens_id = str(item["id"])
        richness = float(item["richness"])
        if not np.isfinite(richness):
            raise ValueError(f"non-finite richness for lens {lens_id}")
        rows.append((lens_id, richness))
    if len({row[0] for row in rows}) != 800:
        raise ValueError("parent lens IDs are not unique")
    # Exactly balanced deterministic split; ID breaks ties at the boundary.
    order = sorted(range(len(rows)), key=lambda index: (rows[index][1], rows[index][0]))
    lower_ids = {rows[index][0] for index in order[:400]}
    label_by_id = {lens_id: (0 if lens_id in lower_ids else 1) for lens_id, _ in rows}
    lower_max = rows[order[399]][1]
    upper_min = rows[order[400]][1]
    return label_by_id, {
        "split_rule": "sort all 800 parent lenses by (CAMIRA richness, lens_id); first 400 are lower-richness",
        "n_parent_lenses": 800,
        "n_parent_lower": 400,
        "n_parent_upper": 400,
        "lower_richness_max_inclusive": lower_max,
        "upper_richness_min_inclusive": upper_min,
        "ties_at_boundary_split_by_lens_id": bool(lower_max == upper_min),
        "richness_definition": "CAMIRA catalog richness from the existing parent sample manifest",
    }


def run(args):
    if args.output.exists():
        raise FileExistsError(f"refusing to overwrite {args.output}")
    pairs, provenance = verify_and_load_pairs(
        args.rows, args.preflight, args.plan_dir, args.masked_value,
    )
    labels_by_id, split = richness_labels(args.sample_manifest)
    lens_id = np.asarray(pairs["lens_id"]).astype(str)
    try:
        labels = np.fromiter((labels_by_id[value] for value in lens_id), dtype=np.int8, count=len(lens_id))
    except KeyError as exc:
        raise ValueError(f"shape export contains a lens absent from the parent manifest: {exc}") from None
    nz_z, nz_bins = load_hsc_y3_nz_sacc(args.nz_sacc)
    cosmology = Cosmology(H0=args.h0, Om=args.omega_m)
    lens_cells = equal_area_block_ids(
        pairs["lens_ra_deg"], pairs["lens_dec_deg"], cell_width_deg=args.cell_width_deg,
    )
    source_cells = equal_area_block_ids(
        pairs["source_ra_deg"], pairs["source_dec_deg"], cell_width_deg=args.cell_width_deg,
    )
    block_ids = np.unique(np.concatenate((lens_cells, source_cells)))
    if len(block_ids) < 3:
        raise ValueError("at least three occupied spatial blocks are required")

    def subset(mask):
        return {key: np.asarray(value)[mask] for key, value in pairs.items()}

    cohort_masks = [labels == 0, labels == 1]
    names = ("lower_richness", "upper_richness")
    full_results = []
    full_vectors = []
    cohort_counts = []
    for name, mask in zip(names, cohort_masks):
        result = stack_s19a_tomographic_pairs(
            subset(mask), nz_z, nz_bins, args.radial_edges_mpc, cosmology,
        )
        full_results.append(result)
        full_vectors.append(_profile_vector(result))
        used_lenses = np.unique(np.asarray(pairs["lens_id"])[mask & (np.asarray(pairs["source_zbin"]).astype(int) > 0)])
        cohort_counts.append({
            "cohort": name,
            "parent_lenses": 400,
            "input_pairs": int(np.count_nonzero(mask)),
            "unique_lenses_with_assigned_source_zbin": int(len(used_lenses)),
            "stack_counts": result["counts"],
        })

    replicates = []
    deletions = []
    for block_id in block_ids:
        keep = (lens_cells != block_id) & (source_cells != block_id)
        vectors = []
        deleted_counts = []
        for mask in cohort_masks:
            keep_cohort = mask & keep
            result = stack_s19a_tomographic_pairs(
                subset(keep_cohort), nz_z, nz_bins, args.radial_edges_mpc, cosmology,
            )
            vectors.append(_profile_vector(result))
            deleted_counts.append(int(np.count_nonzero(mask & ~keep)))
        replicates.append(np.concatenate(vectors))
        deletions.append({
            "block_id": int(block_id),
            "lower_richness_pairs_deleted": deleted_counts[0],
            "upper_richness_pairs_deleted": deleted_counts[1],
        })

    covariance = delete_block_covariance(np.asarray(replicates))
    nbin = len(args.radial_edges_mpc) - 1
    transform = np.zeros((2 * nbin, 4 * nbin), dtype=float)
    for index in range(2 * nbin):
        transform[index, 2 * nbin + index] = 1.0
        transform[index, index] = -1.0
    contrast = full_vectors[1] - full_vectors[0]
    contrast_cov = transform @ covariance @ transform.T
    errors = np.sqrt(np.maximum(np.diag(contrast_cov), 0.0))
    profile = []
    for radial in range(nbin):
        lower_row = full_results[0]["profile"][radial]
        upper_row = full_results[1]["profile"][radial]
        profile.append({
            "bin": radial,
            "r_min_mpc": lower_row["r_min_mpc"],
            "r_max_mpc": lower_row["r_max_mpc"],
            "lower_richness_delta_sigma_msun_mpc2": lower_row["delta_sigma_weak_msun_mpc2"],
            "upper_richness_delta_sigma_msun_mpc2": upper_row["delta_sigma_weak_msun_mpc2"],
            "upper_minus_lower_delta_sigma_msun_mpc2": float(contrast[2 * radial]),
            "upper_minus_lower_delta_sigma_jackknife_error_msun_mpc2": float(errors[2 * radial]),
            "lower_richness_cross_msun_mpc2": lower_row["cross_delta_sigma_weak_msun_mpc2"],
            "upper_richness_cross_msun_mpc2": upper_row["cross_delta_sigma_weak_msun_mpc2"],
            "upper_minus_lower_cross_msun_mpc2": float(contrast[2 * radial + 1]),
            "upper_minus_lower_cross_jackknife_error_msun_mpc2": float(errors[2 * radial + 1]),
            "n_lower_pairs": lower_row["n_pairs"],
            "n_upper_pairs": upper_row["n_pairs"],
        })
    payload = {
        "schema": "astrafek-hsc-y3-richness-split-diagnostic-v1",
        "status": "exploratory",
        "science_measurement": False,
        "random_lens_subtraction_applied": False,
        "photoz_uncertainty_propagated": False,
        "reduced_shear_applied": False,
        "method": "balanced parent-sample richness split; HSC Y3 tomographic stack; shared delete-one equal-area sky-cell jackknife",
        "split": split,
        "cohorts": cohort_counts,
        "cell_width_equatorial_deg": float(args.cell_width_deg),
        "n_occupied_cells": int(len(block_ids)),
        "profile": profile,
        "contrast_covariance_order": "[DeltaSigma bin 0, cross bin 0, DeltaSigma bin 1, cross bin 1, ...]",
        "upper_minus_lower_covariance": contrast_cov.tolist(),
        "block_deletions": deletions,
        "input": {
            **provenance,
            "sample_manifest": str(args.sample_manifest.resolve()),
            "sample_manifest_sha256": sha(args.sample_manifest),
            "n_z_sacc": str(args.nz_sacc.resolve()),
            "n_z_sacc_sha256": sha(args.nz_sacc),
            "radial_edges_mpc": [float(value) for value in args.radial_edges_mpc],
            "cosmology": {"H0_km_s_mpc": args.h0, "Omega_m": args.omega_m},
            "b_mode_masked_raw_value": args.masked_value,
        },
        "limitations": [
            "cluster-only profile; matched random-lens subtraction and matched-field variance are absent",
            "equal-area blocks are not validated contiguous HSC survey regions",
            "split uses parent-sample CAMIRA richness, but lens redshift and selection effects can differ between cohorts",
            "fiducial n(z), shear calibration, source-selection, and reduced-shear uncertainties are not propagated",
            "the contrast is a heterogeneity diagnostic, not a richness-mass calibration or detection significance",
        ],
    }
    args.output.mkdir(parents=True, exist_ok=False)
    target = args.output / "summary.json"
    with target.open("x", encoding="utf-8") as stream:
        json.dump(payload, stream, indent=2, ensure_ascii=False, allow_nan=False)
        stream.write("\n")
    print(json.dumps({
        "status": "exploratory",
        "science_measurement": False,
        "n_occupied_cells": int(len(block_ids)),
        "low_used_pairs": cohort_counts[0]["stack_counts"]["used_pairs"],
        "high_used_pairs": cohort_counts[1]["stack_counts"]["used_pairs"],
        "output": str(target),
    }, ensure_ascii=False))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rows", type=Path, required=True)
    parser.add_argument("--preflight", type=Path, required=True)
    parser.add_argument("--plan-dir", type=Path, required=True)
    parser.add_argument("--sample-manifest", type=Path, required=True)
    parser.add_argument("--nz-sacc", type=Path, required=True)
    parser.add_argument("--masked-value", type=int, choices=(0, 1), required=True)
    parser.add_argument("--cell-width-deg", type=float, default=5.0)
    parser.add_argument("--radial-edges-mpc", type=float, nargs="+", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--h0", type=float, default=70.0)
    parser.add_argument("--omega-m", type=float, default=0.3)
    run(parser.parse_args())


if __name__ == "__main__":
    main()
