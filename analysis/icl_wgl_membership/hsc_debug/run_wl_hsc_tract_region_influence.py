#!/usr/bin/env python3
"""Measure field-scale influence on the exploratory cluster-only HSC lensing stack.

HSC patch-QA tract centers define a reproducible approximate footprint graph.
The graph is used only for leave-one-connected-region-out influence checks;
with few regions it is not a defensible full covariance estimator.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree
from scipy.spatial.distance import cdist

from fetch_ancillary_catalogues import sha
from run_wl_shape_stack import verify_and_load_pairs
from s19a_shape_stack import load_hsc_y3_nz_sacc, stack_s19a_tomographic_pairs
from ogfkit.wlcluster import Cosmology


def unit_vectors(ra_deg, dec_deg):
    ra = np.deg2rad(np.asarray(ra_deg, dtype=float))
    dec = np.deg2rad(np.asarray(dec_deg, dtype=float))
    if (ra.ndim != 1 or dec.shape != ra.shape or not np.isfinite(ra).all()
            or not np.isfinite(dec).all() or np.any((dec < -np.pi / 2) | (dec > np.pi / 2))):
        raise ValueError("RA/Dec must be finite one-dimensional vectors with valid declinations")
    return np.column_stack((np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)))


def tract_centers_from_patch_qa(path: Path):
    """Return spherical mean center vectors keyed by unique HSC tract ID."""
    points = {}
    with path.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(
            (line[2:] if line.startswith("# ") else line for line in stream)
        )
        required = {"skymap_id", "tract", "ra", "dec"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError("patch-QA table lacks skymap_id, tract, ra, or dec")
        for row in reader:
            if not row.get("skymap_id") or not row.get("tract") or not row.get("ra") or not row.get("dec"):
                continue
            # skymap_id is the per-patch identifier here; grouping by it would
            # incorrectly create one "tract" for every patch.
            key = str(row["tract"])
            points.setdefault(key, []).append((float(row["ra"]), float(row["dec"])))
    if not points:
        raise ValueError("patch-QA table contains no tract centers")
    keys = sorted(points)
    vectors = []
    for key in keys:
        coords = np.asarray(points[key], dtype=float)
        vector = unit_vectors(coords[:, 0], coords[:, 1]).mean(axis=0)
        norm = np.linalg.norm(vector)
        if not np.isfinite(norm) or norm == 0:
            raise ValueError(f"invalid spherical center for tract {key}")
        vectors.append(vector / norm)
    return keys, np.asarray(vectors)


def tract_region_labels(centers_xyz, adjacency_deg=2.0):
    """Connected-component IDs for tracts whose centers are within a threshold."""
    centers = np.asarray(centers_xyz, dtype=float)
    if (centers.ndim != 2 or centers.shape[1] != 3 or len(centers) < 3
            or not np.isfinite(centers).all()):
        raise ValueError("tract centers must be a finite [N,3] array")
    if not np.isfinite(adjacency_deg) or not 0 < adjacency_deg <= 180:
        raise ValueError("adjacency threshold must be in (0, 180] degrees")
    chord = 2.0 * np.sin(np.deg2rad(adjacency_deg) / 2.0)
    pairs = cKDTree(centers).query_pairs(chord, output_type="ndarray")
    parent = np.arange(len(centers), dtype=int)

    def root(index):
        while parent[index] != index:
            parent[index] = parent[parent[index]]
            index = parent[index]
        return index

    for first, second in pairs:
        a, b = root(int(first)), root(int(second))
        if a != b:
            parent[max(a, b)] = min(a, b)
    roots = np.asarray([root(index) for index in range(len(centers))])
    _, labels = np.unique(roots, return_inverse=True)
    return labels


def nearest_region_labels(ra_deg, dec_deg, centers_xyz, tract_regions):
    """Assign sky coordinates to the region of their nearest HSC tract center."""
    vectors = unit_vectors(ra_deg, dec_deg)
    indices = cKDTree(np.asarray(centers_xyz, dtype=float)).query(vectors, k=1)[1]
    return np.asarray(tract_regions, dtype=int)[indices]


def _profile_vector(result):
    values = []
    for row in result["profile"]:
        tangential = row["delta_sigma_weak_msun_mpc2"]
        cross = row["cross_delta_sigma_weak_msun_mpc2"]
        if tangential is None or cross is None or not np.isfinite(tangential) or not np.isfinite(cross):
            raise ValueError("a radial bin became empty in a leave-one-region replicate")
        values.extend((float(tangential), float(cross)))
    return np.asarray(values)


def field_region_influence(pairs, nz_redshift, nz_by_bin, radial_edges_mpc,
                           cosmology, lens_manifest, patch_qa, *, adjacency_deg=2.0):
    tract_keys, centers = tract_centers_from_patch_qa(patch_qa)
    tract_regions = tract_region_labels(centers, adjacency_deg)
    manifest_rows = lens_manifest.get("clusters", [])
    if len(manifest_rows) != 800:
        raise ValueError(f"expected the fixed 800-lens manifest, got {len(manifest_rows)}")
    lens_ids = [str(row["id"]) for row in manifest_rows]
    if len(set(lens_ids)) != 800:
        raise ValueError("lens manifest contains duplicate IDs")
    lens_region_values = nearest_region_labels(
        [row["ra_deg"] for row in manifest_rows],
        [row["dec_deg"] for row in manifest_rows], centers, tract_regions,
    )
    lens_region = dict(zip(lens_ids, lens_region_values.tolist()))
    unknown = sorted(set(map(str, pairs["lens_id"])) - set(lens_region))
    if unknown:
        raise ValueError(f"pair table contains lens IDs absent from manifest: {unknown[:5]}")

    pair_lens_region = np.fromiter(
        (lens_region[str(value)] for value in pairs["lens_id"]), dtype=np.int32,
        count=len(pairs["lens_id"]),
    )
    pair_source_region = nearest_region_labels(
        pairs["source_ra_deg"], pairs["source_dec_deg"], centers, tract_regions,
    )
    occupied = np.unique(np.concatenate((pair_lens_region, pair_source_region)))
    if len(occupied) < 3:
        raise ValueError("at least three occupied HSC tract regions are required")

    full = stack_s19a_tomographic_pairs(
        pairs, nz_redshift, nz_by_bin, radial_edges_mpc, cosmology,
    )
    full_vector = _profile_vector(full)
    deletions = []
    for region in occupied:
        keep = (pair_lens_region != region) & (pair_source_region != region)
        replicate = stack_s19a_tomographic_pairs(
            {key: np.asarray(value)[keep] for key, value in pairs.items()},
            nz_redshift, nz_by_bin, radial_edges_mpc, cosmology,
        )
        vector = _profile_vector(replicate)
        profile = []
        for index, row in enumerate(full["profile"]):
            t_index = 2 * index
            x_index = t_index + 1
            t = float(row["delta_sigma_weak_msun_mpc2"])
            x = float(row["cross_delta_sigma_weak_msun_mpc2"])
            profile.append({
                "r_min_mpc": row["r_min_mpc"], "r_max_mpc": row["r_max_mpc"],
                "full_delta_sigma_msun_mpc2": t,
                "loo_delta_sigma_msun_mpc2": float(vector[t_index]),
                "delta_sigma_shift_msun_mpc2": float(vector[t_index] - t),
                "delta_sigma_fractional_shift": float((vector[t_index] - t) / t) if t else None,
                "full_cross_delta_sigma_msun_mpc2": x,
                "loo_cross_delta_sigma_msun_mpc2": float(vector[x_index]),
                "cross_delta_sigma_shift_msun_mpc2": float(vector[x_index] - x),
            })
        deleted_lens_ids = sorted({str(value) for value in pairs["lens_id"][pair_lens_region == region]})
        deleted_source_ids = set(map(str, pairs["object_id"][pair_source_region == region]))
        deletions.append({
            "region_id": int(region),
            "deleted_lenses": len(deleted_lens_ids),
            "deleted_lens_ids": deleted_lens_ids,
            "deleted_sources": len(deleted_source_ids),
            "deleted_pair_rows": int(np.count_nonzero(~keep)),
            "used_pairs_remaining": int(replicate["counts"]["used_pairs"]),
            "profile": profile,
        })

    occupied_by_region = {}
    for lens_id, region in lens_region.items():
        occupied_by_region.setdefault(int(region), []).append(lens_id)
    return {
        "schema": "astrafek-hsc-y3-tract-region-influence-v1",
        "status": "exploratory_field_influence",
        "science_measurement": False,
        "random_lens_subtraction_applied": False,
        "method": "leave one approximate connected HSC tract-center component out; delete rows when either lens or source is assigned to the component's nearest tract center",
        "tract_adjacency_threshold_deg": float(adjacency_deg),
        "n_tracts": len(tract_keys),
        "n_tract_components": int(len(np.unique(tract_regions))),
        "n_occupied_pair_components": int(len(occupied)),
        "tract_center_assignment": "nearest spherical-mean patch-QA center, then connected component of tract centers",
        "region_lens_counts": {str(region): len(ids) for region, ids in sorted(occupied_by_region.items())},
        "full_stack_counts": full["counts"],
        "full_profile": full["profile"],
        "region_deletions": deletions,
        "limitations": [
            "tract-center proximity is an approximate footprint graph, not a polygon-overlap map or official field label",
            "the small number of connected regions is insufficient for a stable full covariance estimate",
            "source region is assigned to the nearest tract center; boundary assignment uncertainty is not propagated",
            "cluster-only stack omits random-lens subtraction, matched-field variance, n(z) uncertainty, shear-calibration uncertainty, and reduced-shear correction",
            "leave-one-region shifts are influence diagnostics, not confidence intervals or detection significance",
        ],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rows", type=Path, required=True)
    parser.add_argument("--preflight", type=Path, required=True)
    parser.add_argument("--plan-dir", type=Path, required=True)
    parser.add_argument("--lens-manifest", type=Path, required=True)
    parser.add_argument("--patch-qa", type=Path, required=True)
    parser.add_argument("--nz-sacc", type=Path, required=True)
    parser.add_argument("--masked-value", type=int, choices=(0, 1), required=True)
    parser.add_argument("--adjacency-deg", type=float, default=2.0)
    parser.add_argument("--radial-edges-mpc", type=float, nargs="+", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--h0", type=float, default=70.0)
    parser.add_argument("--omega-m", type=float, default=0.3)
    args = parser.parse_args()
    if args.output.exists():
        parser.error(f"refusing to overwrite existing output directory: {args.output}")
    for path in (args.lens_manifest, args.patch_qa, args.nz_sacc):
        if not path.is_file():
            parser.error(f"required input does not exist: {path}")
    if len(args.radial_edges_mpc) < 2 or any(
            right <= left for left, right in zip(args.radial_edges_mpc, args.radial_edges_mpc[1:])):
        parser.error("radial edges must be strictly increasing and contain at least two values")
    pairs, provenance = verify_and_load_pairs(
        args.rows, args.preflight, args.plan_dir, args.masked_value,
    )
    nz_redshift, nz_by_bin = load_hsc_y3_nz_sacc(args.nz_sacc)
    cosmology = Cosmology(H0=args.h0, Om=args.omega_m)
    manifest = json.loads(args.lens_manifest.read_text(encoding="utf-8"))
    result = field_region_influence(
        pairs, nz_redshift, nz_by_bin, args.radial_edges_mpc, cosmology,
        manifest, args.patch_qa, adjacency_deg=args.adjacency_deg,
    )
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    payload = {
        **result,
        "input": {
            **provenance,
            "lens_manifest_path": str(args.lens_manifest.resolve()),
            "lens_manifest_sha256": sha(args.lens_manifest),
            "patch_qa_path": str(args.patch_qa.resolve()),
            "patch_qa_sha256": sha(args.patch_qa),
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
        "tracts": result["n_tracts"], "tract_components": result["n_tract_components"],
        "occupied_components": result["n_occupied_pair_components"],
        "region_lens_counts": result["region_lens_counts"],
        "output": str(summary),
    }, ensure_ascii=False))


if __name__ == "__main__":
    main()
