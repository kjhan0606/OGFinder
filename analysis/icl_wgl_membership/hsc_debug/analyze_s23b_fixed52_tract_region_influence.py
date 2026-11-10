#!/usr/bin/env python3
"""Assess regional influence for the fixed 52-lens S23B weak-lensing sample.

Uses only the previously approved shape-row export and pinned n(z). By default,
saved HSC patch-QA centers define approximate connected regions. With
``--geometry-dir``, official HSC PDR3 Wide patch polygons assign endpoints to
their parent tract IDs for single-tract influence blocks. Both modes are
influence audits, not survey-region covariance estimates. The science flag
remains false.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np

import analyze_s23b_wl_fixed_common_support_spatial_covariance as base
import analyze_s23b_wl_fixed_equal_control_spatial_covariance as equal
from fetch_ancillary_catalogues import sha
from ogfkit.wlcluster import Cosmology
from run_wl_hsc_tract_region_influence import (
    nearest_region_labels,
    tract_centers_from_patch_qa,
    tract_region_labels,
)
from hsc_pdr3_tract_patch_geometry import load_index
from s19a_shape_stack import load_hsc_y3_nz_sacc
from submit_s23b_wl_common_shape_rows import DEFAULT_PLAN_DIR, DEFAULT_RESULTS


DATA = base.DATA
DEFAULT_SUPPORT = base.SUPPORT_DIR
DEFAULT_NZ = base.NZ_PATH
DEFAULT_PATCH_QA = DATA / "catalogues_2026-10-13_01/patch_qa.csv"
DEFAULT_REFERENCE = DATA / "wl_s23b_equal_random_center_spatial_covariance_2026-11-11_02/summary.json"
DEFAULT_OUTPUT = DATA / "wl_s23b_fixed52_tract_region_influence_2026-11-11_01"
DEFAULT_EXACT_OUTPUT = DATA / "wl_s23b_fixed52_exact_tract_influence_2026-11-11_01"
EXPECTED_REFERENCE_SHA256 = "00c8427ef6810b78f0900da1d515e3fca2ed77233bd012ebbeddbf5bc5e5c033"
UNASSIGNED_GEOMETRY = "__unassigned__"
CROSS_FIELD_GEOMETRY = "__cross_field_ambiguous__"


def write_json(path: Path, payload: dict) -> None:
    with path.open("x", encoding="utf-8") as stream:
        json.dump(payload, stream, ensure_ascii=False, indent=2, allow_nan=False)
        stream.write("\n")


def write_tsv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty table: {path.name}")
    with path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(
            stream, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def read_reference(path: Path) -> dict[str, np.ndarray]:
    if sha(path) != EXPECTED_REFERENCE_SHA256:
        raise ValueError("pinned fixed-52 pooled/equal estimator reference hash differs")
    payload = base.read_json(path)
    if (payload.get("science_measurement") is not False
            or payload.get("sample", {}).get("n_fixed_lens_sets") != base.EXPECTED_FIXED_LENSES
            or payload.get("provenance", {}).get("approved_row_plan_sha256")
            != base.APPROVED_PLAN_SHA256):
        raise ValueError("pinned reference has unexpected sample or science state")
    width = next(
        item for item in payload["block_width_results"]
        if float(item["block_width_deg"]) == 5.0
    )
    return {
        name: np.asarray([
            [row["delta_sigma_msun_mpc2"], row["cross_delta_sigma_msun_mpc2"]]
            for row in width["estimators"][name]["full_profile"]
        ], dtype=float)
        for name in ("pooled_controls", "equal_control_centers")
    }


def center_region_map(pairs: dict, centers_xyz: np.ndarray,
                      component_ids: np.ndarray) -> tuple[dict[str, int], np.ndarray]:
    ids, first = np.unique(np.asarray(pairs["lens_id"]).astype(str), return_index=True)
    values = nearest_region_labels(
        np.asarray(pairs["lens_ra_deg"], dtype=float)[first],
        np.asarray(pairs["lens_dec_deg"], dtype=float)[first],
        centers_xyz, component_ids,
    )
    lookup = dict(zip(ids.tolist(), values.astype(int).tolist()))
    row_values = np.fromiter(
        (lookup[str(value)] for value in np.asarray(pairs["lens_id"]).astype(str)),
        dtype=np.int32, count=len(pairs["lens_id"]),
    )
    return lookup, row_values


def assign_tract_keys(index, ra_deg, dec_deg) -> tuple[np.ndarray, dict]:
    """Assign coordinates through official patch polygons to their parent tracts."""
    assigned = index.assign(ra_deg, dec_deg, level="tract")
    keys = np.full(len(ra_deg), UNASSIGNED_GEOMETRY, dtype=object)
    for row, (field, tract, cross_field) in enumerate(zip(
            assigned.field, assigned.tract, assigned.cross_field_ambiguous)):
        if cross_field:
            keys[row] = CROSS_FIELD_GEOMETRY
        elif field is not None and int(tract) >= 0:
            keys[row] = f"{field}:tract-{int(tract)}"
    return keys, {
        "n_points": int(len(keys)),
        "unassigned": int(np.count_nonzero(keys == UNASSIGNED_GEOMETRY)),
        "cross_field_ambiguous": int(np.count_nonzero(keys == CROSS_FIELD_GEOMETRY)),
        "multiple_patch_polygon_matches": assigned.overlap_count,
    }


def exact_region_maps(index, lens_pairs: dict, control_pairs: dict) -> tuple[dict, dict]:
    """Return exact tract IDs for aperture centers and unique source positions.

    Source coordinates are assigned once per unique object_id, then expanded
    back to center/source pair rows. This avoids repeating point-in-polygon
    work for sources reused by multiple apertures.
    """
    def center_key_map(pairs):
        ids, first = np.unique(np.asarray(pairs["lens_id"]).astype(str), return_index=True)
        keys, audit = assign_tract_keys(
            index,
            np.asarray(pairs["lens_ra_deg"], dtype=float)[first],
            np.asarray(pairs["lens_dec_deg"], dtype=float)[first],
        )
        return dict(zip(ids.tolist(), keys.tolist())), audit

    lens_center_keys, lens_center_audit = center_key_map(lens_pairs)
    control_center_keys, control_center_audit = center_key_map(control_pairs)
    lens_n = len(lens_pairs["object_id"])
    source_ids = np.concatenate((
        np.asarray(lens_pairs["object_id"]).astype(str),
        np.asarray(control_pairs["object_id"]).astype(str),
    ))
    source_ra = np.concatenate((
        np.asarray(lens_pairs["source_ra_deg"], dtype=float),
        np.asarray(control_pairs["source_ra_deg"], dtype=float),
    ))
    source_dec = np.concatenate((
        np.asarray(lens_pairs["source_dec_deg"], dtype=float),
        np.asarray(control_pairs["source_dec_deg"], dtype=float),
    ))
    unique_ids, first, inverse = np.unique(source_ids, return_index=True, return_inverse=True)
    unique_source_keys, source_audit = assign_tract_keys(
        index, source_ra[first], source_dec[first],
    )
    all_keys = set(unique_source_keys.tolist())
    all_keys.update(lens_center_keys.values())
    all_keys.update(control_center_keys.values())
    valid_keys = sorted(all_keys - {UNASSIGNED_GEOMETRY, CROSS_FIELD_GEOMETRY})
    if not valid_keys:
        raise ValueError("no endpoints were assigned to an official HSC tract")
    key_to_id = {key: index for index, key in enumerate(valid_keys)}

    def encode(key):
        return int(key_to_id.get(key, -1))

    def encode_rows(pairs, center_map, source_codes):
        center_codes = np.fromiter(
            (encode(center_map[str(value)]) for value in np.asarray(pairs["lens_id"]).astype(str)),
            dtype=np.int32, count=len(pairs["lens_id"]),
        )
        return center_codes, source_codes

    unique_source_codes = np.fromiter(
        (encode(str(value)) for value in unique_source_keys),
        dtype=np.int32, count=len(unique_source_keys),
    )
    lens_source_codes = unique_source_codes[inverse[:lens_n]]
    control_source_codes = unique_source_codes[inverse[lens_n:]]
    lens_center_codes, lens_source_codes = encode_rows(
        lens_pairs, lens_center_keys, lens_source_codes,
    )
    control_center_codes, control_source_codes = encode_rows(
        control_pairs, control_center_keys, control_source_codes,
    )
    return {
        "lens_center_lookup": {key: encode(value) for key, value in lens_center_keys.items()},
        "control_center_lookup": {key: encode(value) for key, value in control_center_keys.items()},
        "lens_center": lens_center_codes,
        "control_center": control_center_codes,
        "lens_source": lens_source_codes,
        "control_source": control_source_codes,
    }, {
        "assignment": "official PDR3 Wide patch polygons; parent tract ID is assigned from the nearest containing patch center",
        "n_unique_source_ids_assigned": int(len(unique_ids)),
        "lens_centers": lens_center_audit,
        "control_centers": control_center_audit,
        "unique_sources": source_audit,
        "region_labels": {str(value): key for key, value in key_to_id.items()},
    }


def run(args) -> dict:
    output = (args.output or (
        DEFAULT_EXACT_OUTPUT if args.geometry_dir else DEFAULT_OUTPUT
    )).resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite output: {output}")
    reference_profiles = read_reference(args.reference.resolve())

    lens_ids, fixed_selection = base.fixed_sample_ids(args.support_dir.resolve())
    lens_pairs, all_controls, center_to_lens, row_provenance = base.load_approved_pairs(
        args.plan_dir.resolve(), args.rows.resolve(),
    )
    controls_by_lens: dict[str, list[str]] = {lens_id: [] for lens_id in lens_ids}
    for center_id, lens_id in center_to_lens.items():
        if lens_id in lens_ids:
            controls_by_lens[lens_id].append(center_id)
    if any(len(centers) != 5 for centers in controls_by_lens.values()):
        raise ValueError("a fixed lens set no longer has exactly five matched control centers")
    control_ids = {center for centers in controls_by_lens.values() for center in centers}
    lens_pairs = base.select_rows(lens_pairs, lens_ids)
    control_pairs = base.select_rows(all_controls, control_ids)
    del all_controls

    source_reuse = base.verify_repeated_source_rows(lens_pairs, control_pairs)
    geometry_summary = None
    geometry_audit = None
    if args.geometry_dir:
        geometry_index, geometry_summary = load_index(args.geometry_dir.resolve(), level="tract")
        exact, geometry_audit = exact_region_maps(geometry_index, lens_pairs, control_pairs)
        lens_region_lookup = exact["lens_center_lookup"]
        control_region_lookup = exact["control_center_lookup"]
        lens_center_region = exact["lens_center"]
        control_center_region = exact["control_center"]
        lens_source_region = exact["lens_source"]
        control_source_region = exact["control_source"]
        region_label_by_id = {
            int(key): str(value) for key, value in geometry_audit["region_labels"].items()
        }
        tract_keys = list(range(geometry_summary["n_patch_parent_tracts"]))
        tract_components = None
    else:
        tract_keys, tract_xyz = tract_centers_from_patch_qa(args.patch_qa.resolve())
        tract_components = tract_region_labels(tract_xyz, args.adjacency_deg)
        lens_region_lookup, lens_center_region = center_region_map(
            lens_pairs, tract_xyz, tract_components,
        )
        control_region_lookup, control_center_region = center_region_map(
            control_pairs, tract_xyz, tract_components,
        )
        lens_source_region = nearest_region_labels(
            lens_pairs["source_ra_deg"], lens_pairs["source_dec_deg"], tract_xyz, tract_components,
        )
        control_source_region = nearest_region_labels(
            control_pairs["source_ra_deg"], control_pairs["source_dec_deg"], tract_xyz, tract_components,
        )
        region_label_by_id = {
            int(value): f"approx_component_{int(value)}"
            for value in np.unique(tract_components)
        }

    nz_z, nz_bins = load_hsc_y3_nz_sacc(args.nz_sacc.resolve())
    cosmology = Cosmology(H0=70.0, Om=0.3)
    lens_cache = base.cache_geometry_and_efficiency(lens_pairs, nz_z, nz_bins, cosmology)
    control_cache = base.cache_geometry_and_efficiency(control_pairs, nz_z, nz_bins, cosmology)
    lens_rows_by_set = base.group_indices(lens_pairs["lens_id"])
    control_rows_by_center = base.group_indices(control_pairs["lens_id"])

    full_by_lens = {}
    full_set_profiles = {"pooled_controls": {}, "equal_control_centers": {}}
    for lens_id in sorted(lens_ids):
        target_indices = lens_rows_by_set.get(lens_id, np.asarray([], dtype=np.int64))
        center_ids = sorted(controls_by_lens[lens_id])
        center_indices = {
            center_id: control_rows_by_center.get(center_id, np.asarray([], dtype=np.int64))
            for center_id in center_ids
        }
        random_indices = np.concatenate([center_indices[center_id] for center_id in center_ids])
        target_result = base.cached_stack(lens_pairs, lens_cache, target_indices, cosmology)
        pooled_result = base.cached_stack(control_pairs, control_cache, random_indices, cosmology)
        equal_result = equal.equal_center_controls(
            control_pairs, control_cache, center_indices, center_ids, cosmology,
        )
        pooled_delta, equal_delta, n_centers = equal.subtracted_profile(
            target_result, pooled_result, equal_result,
        )
        if (np.any(~np.isfinite(pooled_delta)) or np.any(~np.isfinite(equal_delta))
                or np.any(n_centers < base.MIN_CONTROLS_PER_ANNULUS)):
            raise ValueError(f"full fixed-sample support changed for lens {lens_id}")
        full_set_profiles["pooled_controls"][lens_id] = pooled_delta
        full_set_profiles["equal_control_centers"][lens_id] = equal_delta
        full_by_lens[lens_id] = {
            "target_indices": target_indices,
            "random_indices": random_indices,
            "center_indices": center_indices,
            "center_ids": center_ids,
            "target": target_result,
            "pooled": pooled_result,
            "equal": equal_result,
        }

    full_profiles = {
        name: np.mean(np.stack(list(values.values())), axis=0)
        for name, values in full_set_profiles.items()
    }
    full_profiles["equal_minus_pooled"] = (
        full_profiles["equal_control_centers"] - full_profiles["pooled_controls"]
    )
    for name, reference in reference_profiles.items():
        if not np.allclose(full_profiles[name], reference, rtol=2e-9, atol=1e-3):
            raise ValueError(f"recomputed {name} profile does not match its pinned reference")

    active_lens = lens_cache["active"]
    active_control = control_cache["active"]
    if (np.any(lens_center_region[active_lens] < 0)
            or np.any(lens_source_region[active_lens] < 0)):
        raise ValueError("an active lens endpoint is outside the assignable HSC Wide tract geometry")
    if (np.any(control_center_region[active_control] < 0)
            or np.any(control_source_region[active_control] < 0)):
        raise ValueError("an active control endpoint is outside the assignable HSC Wide tract geometry")
    occupied = np.unique(np.concatenate((
        lens_center_region[active_lens], lens_source_region[active_lens],
        control_center_region[active_control], control_source_region[active_control],
    )))
    if len(occupied) < 2:
        raise ValueError("fewer than two spatial regions overlap the active stack")

    region_rows, profile_rows = [], []
    for region in occupied:
        keep_lens = (lens_center_region != region) & (lens_source_region != region)
        keep_control = (control_center_region != region) & (control_source_region != region)
        accepted = {name: {} for name in (
            "pooled_controls", "equal_control_centers", "equal_minus_pooled",
        )}
        min_center_support = []
        affected_sets = 0

        for lens_id in sorted(lens_ids):
            entry = full_by_lens[lens_id]
            target_indices = entry["target_indices"]
            random_indices = entry["random_indices"]
            affected = (
                not np.all(keep_lens[target_indices])
                or not np.all(keep_control[random_indices])
            )
            if affected:
                affected_sets += 1
                target_result = base.cached_stack(
                    lens_pairs, lens_cache,
                    target_indices[keep_lens[target_indices]], cosmology,
                )
                filtered_random = random_indices[keep_control[random_indices]]
                pooled_result = base.cached_stack(
                    control_pairs, control_cache, filtered_random, cosmology,
                )
                equal_result = equal.equal_center_controls(
                    control_pairs, control_cache, entry["center_indices"], entry["center_ids"],
                    cosmology, keep_mask=keep_control,
                )
            else:
                target_result, pooled_result, equal_result = (
                    entry["target"], entry["pooled"], entry["equal"],
                )
            pooled_delta, equal_delta, n_centers = equal.subtracted_profile(
                target_result, pooled_result, equal_result,
            )
            finite = np.all(np.isfinite(pooled_delta), axis=1) & np.all(np.isfinite(equal_delta), axis=1)
            supported = finite & (n_centers >= base.MIN_CONTROLS_PER_ANNULUS)
            if np.all(supported):
                accepted["pooled_controls"][lens_id] = pooled_delta
                accepted["equal_control_centers"][lens_id] = equal_delta
                accepted["equal_minus_pooled"][lens_id] = equal_delta - pooled_delta
                min_center_support.extend(n_centers.tolist())

        common_ids = set(accepted["pooled_controls"])
        if any(set(accepted[name]) != common_ids for name in accepted):
            raise ValueError("regional estimator profiles do not share common lens-set support")
        if not common_ids:
            raise ValueError(f"region deletion {region} leaves no common supported lens sets")
        replicas = {
            name: np.mean(np.stack(list(values.values())), axis=0)
            for name, values in accepted.items()
        }
        status = "influence_only_insufficient_for_covariance"
        for name, replica in replicas.items():
            shift = replica - full_profiles[name]
            for bin_index, (lo, hi) in enumerate(zip(
                    base.RADIAL_EDGES_MPC[:-1], base.RADIAL_EDGES_MPC[1:])):
                full_t, full_x = full_profiles[name][bin_index]
                rep_t, rep_x = replica[bin_index]
                profile_rows.append({
                    "region_id": int(region), "estimator": name,
                    "region_label": region_label_by_id[int(region)],
                    "radial_bin": bin_index, "r_min_mpc": float(lo), "r_max_mpc": float(hi),
                    "common_lens_sets": len(common_ids),
                    "full_delta_sigma_msun_mpc2": float(full_t),
                    "leave_region_out_delta_sigma_msun_mpc2": float(rep_t),
                    "delta_sigma_shift_msun_mpc2": float(shift[bin_index, 0]),
                    "delta_sigma_fractional_shift": float(shift[bin_index, 0] / full_t) if full_t else None,
                    "full_cross_msun_mpc2": float(full_x),
                    "leave_region_out_cross_msun_mpc2": float(rep_x),
                    "cross_shift_msun_mpc2": float(shift[bin_index, 1]),
                    "cross_fractional_shift": float(shift[bin_index, 1] / full_x) if full_x else None,
                })
        lens_set_counts = {
            int(component): int(sum(value == component for value in lens_region_lookup.values()))
            for component in occupied
        }
        control_counts = {
            int(component): int(sum(value == component for value in control_region_lookup.values()))
            for component in occupied
        }
        region_rows.append({
            "region_id": int(region),
            "region_label": region_label_by_id[int(region)],
            "lens_centers_in_region": lens_set_counts[int(region)],
            "control_centers_in_region": control_counts[int(region)],
            "affected_lens_sets": affected_sets,
            "deleted_active_lens_pair_rows": int(np.count_nonzero(~keep_lens & active_lens)),
            "deleted_active_control_pair_rows": int(np.count_nonzero(~keep_control & active_control)),
            "common_supported_lens_sets_after_deletion": len(common_ids),
            "minimum_positive_control_centers_in_common_sets": min(min_center_support),
            "status": status,
        })

    output.mkdir(parents=True, exist_ok=False)
    write_tsv(output / "region_deletions.tsv", region_rows)
    write_tsv(output / "region_profile_shifts.tsv", profile_rows)
    reference_sha = sha(args.reference.resolve())
    patch_sha = sha(args.patch_qa.resolve()) if not args.geometry_dir else None
    exact_mode = args.geometry_dir is not None
    method = ({
        "region_definition": "parent tract IDs from official HSC PDR3 Wide patch polygons",
        "endpoint_assignment": "point-in-polygon over patch polygons; nearest containing same-field patch center resolves overlaps",
        "delete_rule": "remove a lens/control pair when either center or source endpoint is assigned to the deleted tract; recompute calibrations and the four-control radial support rule",
        "estimators": ["pooled_controls", "equal_control_centers", "equal_minus_pooled"],
        "profile_comparison": "each leave-tract-out profile compared with the fixed 52-set full-sample estimator",
    } if exact_mode else {
        "region_definition": "connected components of spherical-mean HSC patch-QA tract centers",
        "tract_adjacency_deg": float(args.adjacency_deg),
        "endpoint_assignment": "nearest tract center, then its connected component",
        "delete_rule": "remove a lens/control pair when either center or source endpoint is assigned to the deleted component; recompute calibrations and the four-control radial support rule",
        "estimators": ["pooled_controls", "equal_control_centers", "equal_minus_pooled"],
        "profile_comparison": "each leave-region-out profile compared with the fixed 52-set full-sample estimator",
    })
    input_provenance = {
        "row_plan_sha256": row_provenance["plan_sha256"],
        "row_status_sha256": row_provenance["row_status_sha256"],
        "row_summary_sha256": row_provenance["row_summary_sha256"],
        "row_batch_sha256": row_provenance["batch_sha256"],
        "fixed_sample_support_summary_sha256": sha(args.support_dir.resolve() / "summary.json"),
        "fixed_sample_support_rows_sha256": sha(args.support_dir.resolve() / "per_lens_radial_support.tsv"),
        "hsc_y3_nz_sacc_sha256": sha(args.nz_sacc.resolve()),
        "reference_spatial_covariance_summary_sha256": reference_sha,
        "fixed_lens_ids_sha256": fixed_selection["fixed_lens_ids_sha256"],
    }
    if exact_mode:
        input_provenance["official_geometry_directory"] = str(args.geometry_dir.resolve())
        input_provenance["official_geometry_files"] = geometry_summary["files"]
        input_provenance["official_geometry_parser_sha256"] = sha(
            Path(__file__).with_name("hsc_pdr3_tract_patch_geometry.py")
        )
    else:
        input_provenance["patch_qa_sha256"] = patch_sha
    summary = {
        "schema": "astrafek-s23b-fixed52-tract-region-influence-v2",
        "status": "complete_exact_tract_influence_only" if exact_mode else "complete_exploratory_regional_influence",
        "science_measurement": False,
        "cas_query_submitted": False,
        "new_data_downloaded": False,
        "images_read": False,
        "covariance_estimated": False,
        "halo_mass_fit_performed": False,
        "sample": {
            **fixed_selection,
            "n_selected_lens_sets": len(lens_ids),
            "n_selected_control_centers": len(control_ids),
            "lens_source_pair_rows": int(len(lens_pairs["lens_id"])),
            "control_source_pair_rows": int(len(control_pairs["lens_id"])),
            "source_reuse_audit": source_reuse,
        },
        "method": method,
        "tracts_total": len(tract_keys),
        "tract_components_total": (None if exact_mode else int(len(np.unique(tract_components)))),
        "occupied_spatial_blocks": len(occupied),
        "lens_region_counts": {str(row["region_label"]): row["lens_centers_in_region"] for row in region_rows},
        "control_region_counts": {str(row["region_label"]): row["control_centers_in_region"] for row in region_rows},
        "official_geometry_assignment_audit": geometry_audit,
        "official_geometry_summary": geometry_summary,
        "full_profiles_msun_mpc2": {name: values.tolist() for name, values in full_profiles.items()},
        "region_deletions": region_rows,
        "inputs": input_provenance,
        "limitations": ([
            "individual tracts are influence blocks, not independent survey regions; covariance is deliberately not estimated",
            "tract IDs are assigned from the nearest containing official patch polygon; patch overlap handling creates a deterministic spatial block label, not a reconstruction of the catalog's owning patch at source-row query time",
        ] if exact_mode else [
            "the approved pair rows have coordinates but no tract or patch IDs; endpoint assignment uses nearest patch-QA tract center",
            "the tract-center graph approximates footprint connectivity; it is not a polygon-overlap map or official HSC field label",
            "only a few connected components overlap the sample, so these are influence checks and cannot estimate a stable survey-region covariance",
        ]) + [
            "the fixed-support sample is selected for radial support; the leave-region-out common sample can change",
            "photo-z uncertainty, reduced shear, magnification, boost/member dilution, shear-selection systematics, and PSF uncertainties remain unpropagated",
        ],
        "outputs": {
            "region_deletions_sha256": sha(output / "region_deletions.tsv"),
            "region_profile_shifts_sha256": sha(output / "region_profile_shifts.tsv"),
        },
        "output_directory": str(output),
    }

    region_word = "tract" if exact_mode else "approximate region"
    region_description = (
        "Official HSC PDR3 Wide patch polygons assign the center and unique source endpoints to their parent tract IDs. Same-field patch overlaps use the nearest containing patch center; active endpoints that are unassigned or cross-field ambiguous stop the run."
        if exact_mode else
        "Connected components of patch-QA tract centers at the configured angular threshold define the approximate regions."
    )
    assessment = [
        f"# S23B fixed-52 {region_word} influence audit",
        "",
        "This CPU analysis reuses only the previously approved 20-file HSC S19A shape-row export and pinned HSC Y3 SACC n(z). It makes no CAS request, additional download, or image read. `science_measurement=false`.",
        "",
        f"{region_description} Each lens/control pair is removed if either endpoint maps to the deleted {region_word}. The calibrated pooled-control and equal-control-center profiles are recomputed under the same four-control floor.",
        "",
        "| Region | Lens centers | Control centers | Affected lens sets | Common supported sets after deletion | Minimum positive controls |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    assessment.extend(
        f"| {row['region_label']} | {row['lens_centers_in_region']} | {row['control_centers_in_region']} | {row['affected_lens_sets']} | {row['common_supported_lens_sets_after_deletion']} | {row['minimum_positive_control_centers_in_common_sets']} |"
        for row in region_rows
    )
    assessment.extend([
        "",
        "The reported regional shifts are influence diagnostics, not errors or significance. Individual adjacent HSC tracts are not independent blocks, so this runner does not estimate covariance. `region_profile_shifts.tsv` gives the per-bin pooled, equal-center, and paired-contrast shifts.",
        "",
        "The fixed sample is selected for radial support, and the retained common lens sets can change after deletion. Other unresolved limits include n(z) uncertainty, reduced shear, magnification, member dilution/boost, shear-selection systematics, and PSF uncertainty. Keep `science_measurement=false`; do not infer a mass or detection.",
        "",
        f"- Runner SHA-256: `{sha(Path(__file__).resolve())}`",
        f"- Approved row-plan SHA-256: `{row_provenance['plan_sha256']}`",
        f"- Reference estimator/covariance summary SHA-256: `{reference_sha}`",
        (f"- Official geometry directory: `{args.geometry_dir.resolve()}`" if exact_mode
         else f"- HSC patch-QA SHA-256: `{patch_sha}`"),
        *( [f"- Official geometry parser SHA-256: `{input_provenance['official_geometry_parser_sha256']}`"]
           if exact_mode else [] ),
        f"- Result directory: `{output}`",
        "",
    ])
    (output / "assessment.md").write_text("\n".join(assessment), encoding="utf-8")
    summary["outputs"]["assessment_sha256"] = sha(output / "assessment.md")
    write_json(output / "summary.json", summary)
    print(json.dumps({
        "status": summary["status"], "science_measurement": False,
        "occupied_spatial_blocks": len(occupied), "region_deletions": region_rows,
        "output_directory": str(output),
        "summary_sha256": sha(output / "summary.json"),
    }, ensure_ascii=False, indent=2), flush=True)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan-dir", type=Path, default=DEFAULT_PLAN_DIR)
    parser.add_argument("--rows", type=Path, default=DEFAULT_RESULTS)
    parser.add_argument("--support-dir", type=Path, default=DEFAULT_SUPPORT)
    parser.add_argument("--nz-sacc", type=Path, default=DEFAULT_NZ)
    parser.add_argument("--patch-qa", type=Path, default=DEFAULT_PATCH_QA)
    parser.add_argument("--reference", type=Path, default=DEFAULT_REFERENCE)
    parser.add_argument("--output", type=Path,
                        help="new output directory; defaults to an approximate or exact-geometry run directory")
    parser.add_argument("--geometry-dir", type=Path,
                        help="directory containing all four local HSC PDR3 Wide tract/patch coordinate files; validates reviewed names and sizes, never downloads")
    parser.add_argument("--adjacency-deg", type=float, default=2.0)
    run(parser.parse_args())


if __name__ == "__main__":
    main()
