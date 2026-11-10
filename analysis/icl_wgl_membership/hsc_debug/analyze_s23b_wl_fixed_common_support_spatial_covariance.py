#!/usr/bin/env python3
"""Estimate paired endpoint spatial covariance for the fixed S23B 52-lens profile.

The estimator keeps the preselected 52 matched lens sets fixed. Within a set,
the target stack is calibrated independently and its matched random centers are
pooled by source-pair lensing weights; the outer set mean is equal-weighted.
Delete-one equal-area sky cells remove a pair whenever either the lens/random
center endpoint or source endpoint falls in the cell. The full HSC S19A
calibration is rerun for every affected set and deletion. Only the already
approved source-row CSVs and pinned Astrafek support products are read.

This is an exploratory spatial covariance diagnostic. Equal-area cells are not
validated contiguous HSC regions, and the output is always
``science_measurement=false``.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
ROOT = HERE.parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from analyze_s23b_wl_random_subtracted_profile import load_approved_pairs
from fetch_ancillary_catalogues import sha
from ogfkit.wlcluster import Cosmology
from run_wl_cluster_spatial_jackknife import equal_area_block_ids
from s19a_shape_stack import (
    load_hsc_y3_nz_sacc,
    pair_geometry,
    stack_s19a_pairs,
    tomographic_nz_lensing_efficiency,
)
from submit_s23b_wl_common_shape_rows import (
    APPROVED_PLAN_SHA256,
    DEFAULT_PLAN_DIR,
    DEFAULT_RESULTS,
)


DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
SUPPORT_DIR = DATA / "wl_s23b_matched_radial_support_2026-11-11_02"
FIXED_PROFILE_DIR = DATA / "wl_s23b_fixed_common_support_profile_2026-11-11_02"
QA_DIR = DEFAULT_PLAN_DIR / "quality_audit_2026-11-11_01"
NZ_PATH = DATA / "hsc_y3_nz_sacc_2026-10-19_01/hsc_y3_fourier_space_data_vector.sacc"
DEFAULT_OUTPUT = DATA / "wl_s23b_fixed_common_support_spatial_covariance_2026-11-11_01"

EXPECTED_SUPPORT_SUMMARY_SHA256 = "c51a4d9268935a0919d0abd3da2f10b6cd16ee7e35f1e6e23a8f3fe8c9e51234"
EXPECTED_PER_LENS_SUPPORT_SHA256 = "6744ece545781faffebef65fc7a0143284da3a234425b6bed540bd1e9d19625c"
EXPECTED_FIXED_PROFILE_SUMMARY_SHA256 = "4d270b81737e5bcb6188941f9f26cd6255614fd7ae8d1e6f2977cf51a4686930"
EXPECTED_QA_SUMMARY_SHA256 = "eb73aa0365813d880d3edf60b12e4039f2575f8a97eab8868fc918121a34d445"
EXPECTED_NZ_SHA256 = "a28f9e2e088e92d96d2d87083a6eeeac4c0c84c6b48e033bd3d1dc1bf0d8958f"
EXPECTED_FIXED_LENS_SHA256 = "18312f27139b3750abdbf71543f22d9b7b18b2d3b24cbd29ed686a377524bb7d"
RADIAL_EDGES_MPC = np.asarray([0.1, 0.2, 0.4, 0.7, 1.0, 1.5, 2.0, 3.0])
SOURCE_BINS = {2, 3, 4}
MIN_CONTROLS_PER_ANNULUS = 4
EXPECTED_FIXED_LENSES = 52
EXPECTED_DIMENSION = 14
REQUIRED_PAIR_COLUMNS = (
    "lens_id", "object_id", "lens_ra_deg", "lens_dec_deg", "lens_z",
    "source_ra_deg", "source_dec_deg", "e1", "e2", "sigma_e", "rms_e",
    "weight", "m", "c1", "c2", "mag_a10", "resolution", "psf_shape11",
    "psf_shape22", "psf_shape12", "b_mode_pass", "source_zbin",
)


def read_json(path: Path) -> dict:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def read_tsv(path: Path) -> list[dict[str, str]]:
    with Path(path).open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def write_json(path: Path, payload: dict) -> None:
    Path(path).write_text(
        json.dumps(payload, ensure_ascii=False, indent=2, allow_nan=False) + "\n",
        encoding="utf-8",
    )


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with Path(path).open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def digest_ids(values) -> str:
    return hashlib.sha256("\n".join(sorted(map(str, values))).encode()).hexdigest()


def group_indices(values: np.ndarray) -> dict[str, np.ndarray]:
    values = np.asarray(values).astype(str)
    if values.size == 0:
        return {}
    order = np.argsort(values, kind="mergesort")
    ordered = values[order]
    starts = np.r_[0, np.flatnonzero(ordered[1:] != ordered[:-1]) + 1, len(ordered)]
    return {
        str(ordered[start]): order[start:end]
        for start, end in zip(starts[:-1], starts[1:])
    }


def select_rows(pairs: dict[str, np.ndarray], centers: set[str]) -> dict[str, np.ndarray]:
    center_ids = np.asarray(pairs["lens_id"]).astype(str)
    bins = np.asarray(pairs["source_zbin"], dtype=int)
    calibration_valid = (
        (np.asarray(pairs["m"], dtype=float) > -1.0)
        & ((np.asarray(pairs["psf_shape11"], dtype=float)
            + np.asarray(pairs["psf_shape22"], dtype=float)) > 0.0)
    )
    keep = np.isin(center_ids, list(centers)) & np.isin(bins, list(SOURCE_BINS)) & calibration_valid
    if not np.any(keep):
        raise ValueError("fixed common-support selection produced no source rows")
    return {key: np.asarray(value)[keep] for key, value in pairs.items()}


def verify_repeated_source_rows(lens_pairs, control_pairs) -> dict:
    fields = ("object_id", "source_ra_deg", "source_dec_deg", "source_zbin")
    combined = {key: np.concatenate((lens_pairs[key], control_pairs[key])) for key in fields}
    order = np.argsort(combined["object_id"].astype(str), kind="mergesort")
    ids = combined["object_id"].astype(str)[order]
    repeated = ids[1:] == ids[:-1]
    if np.any(repeated):
        ra = np.asarray(combined["source_ra_deg"], dtype=float)[order]
        dec = np.asarray(combined["source_dec_deg"], dtype=float)[order]
        bins = np.asarray(combined["source_zbin"], dtype=int)[order]
        bad = repeated & (
            (np.abs(ra[1:] - ra[:-1]) > 1e-10)
            | (np.abs(dec[1:] - dec[:-1]) > 1e-10)
            | (bins[1:] != bins[:-1])
        )
        if np.any(bad):
            raise ValueError("a repeated source ID has inconsistent coordinates or HSC z-bin")
    n_unique = int(len(ids) - np.count_nonzero(repeated))
    lens_unique = np.unique(np.asarray(lens_pairs["object_id"]).astype(str))
    control_unique = np.unique(np.asarray(control_pairs["object_id"]).astype(str))
    n_shared = int(np.intersect1d(lens_unique, control_unique, assume_unique=True).size)
    return {
        "pair_rows": int(len(ids)),
        "unique_source_ids": n_unique,
        "source_ids_in_both_lens_and_control_rows": n_shared,
        "repeated_source_coordinate_and_zbin_mismatches": 0,
    }


def verify_repeated_cells(ids, cells, label: str) -> None:
    ids = np.asarray(ids).astype(str)
    cells = np.asarray(cells, dtype=np.int64)
    order = np.argsort(ids, kind="mergesort")
    ordered_ids, ordered_cells = ids[order], cells[order]
    repeated = ordered_ids[1:] == ordered_ids[:-1]
    if np.any(repeated & (ordered_cells[1:] != ordered_cells[:-1])):
        raise ValueError(f"repeated {label} IDs map to more than one equal-area cell")


def fixed_sample_ids(support_dir: Path) -> tuple[set[str], dict]:
    summary_path = support_dir / "summary.json"
    support_path = support_dir / "per_lens_radial_support.tsv"
    if sha(summary_path) != EXPECTED_SUPPORT_SUMMARY_SHA256:
        raise ValueError("matched radial-support summary hash differs from its pinned artifact")
    if sha(support_path) != EXPECTED_PER_LENS_SUPPORT_SHA256:
        raise ValueError("per-lens radial-support table hash differs from its pinned artifact")
    summary = read_json(summary_path)
    if summary.get("science_measurement") is not False:
        raise ValueError("radial-support input unexpectedly promotes the analysis")
    rule = summary.get("fixed_all_radius_support", {}).get("4", {})
    if rule.get("n_lens_sets_supported_in_all_seven_annuli") != EXPECTED_FIXED_LENSES:
        raise ValueError("pinned fixed ≥4-control sample size changed")
    if rule.get("lens_ids_sha256") != EXPECTED_FIXED_LENS_SHA256:
        raise ValueError("pinned fixed lens-set membership hash changed")

    grouped = {}
    for row in read_tsv(support_path):
        grouped.setdefault(str(row["matched_lens_id"]), []).append(row)
    selected = {
        lens_id for lens_id, rows in grouped.items()
        if len(rows) == 7 and all(
            int(row["n_positive_controls_z234"]) >= MIN_CONTROLS_PER_ANNULUS
            and int(row["lens_pair_rows_z234"]) > 0
            for row in rows
        )
    }
    if len(selected) != EXPECTED_FIXED_LENSES or digest_ids(selected) != EXPECTED_FIXED_LENS_SHA256:
        raise ValueError("fixed lens IDs derived from the per-lens table do not match the pinned sample")
    return selected, {
        "n_fixed_lens_sets": len(selected),
        "fixed_lens_ids_sha256": digest_ids(selected),
        "rule": "positive target support and at least four positive matched controls in every 0.1–3 proper-Mpc annulus, source bins 2–4",
    }


def active_pair_mask(pairs, beta, radius) -> np.ndarray:
    ellipticity_amplitude = np.hypot(
        np.asarray(pairs["e1"], dtype=float), np.asarray(pairs["e2"], dtype=float),
    )
    shape = (
        (ellipticity_amplitude < 2.0)
        & (np.asarray(pairs["sigma_e"], dtype=float) >= 0.0)
        & (np.asarray(pairs["sigma_e"], dtype=float) <= 0.4)
        & (np.asarray(pairs["mag_a10"], dtype=float) <= 25.5)
        & (np.asarray(pairs["resolution"], dtype=float) >= 0.3)
    )
    return (
        shape
        & np.asarray(pairs["b_mode_pass"], dtype=bool)
        & (np.asarray(beta, dtype=float) > 0)
        & (np.asarray(radius, dtype=float) >= RADIAL_EDGES_MPC[0])
        & (np.asarray(radius, dtype=float) <= RADIAL_EDGES_MPC[-1])
    )


def cache_geometry_and_efficiency(pairs, nz_z, nz_bins, cosmology) -> dict:
    beta = tomographic_nz_lensing_efficiency(
        np.asarray(pairs["lens_z"], dtype=float),
        np.asarray(pairs["source_zbin"], dtype=int),
        nz_z, nz_bins, cosmology,
    )
    radius, phi = pair_geometry(pairs, cosmology)
    if not (beta.shape == radius.shape == phi.shape == (len(pairs["lens_id"]),)):
        raise ValueError("cached n(z) geometry arrays do not match selected pair rows")
    return {
        "beta": np.asarray(beta, dtype=float),
        "radius": np.asarray(radius, dtype=float),
        "phi": np.asarray(phi, dtype=float),
        "active": active_pair_mask(pairs, beta, radius),
    }


def cached_stack(pairs, cache, indices, cosmology, *, selection_bias_coefficients=None):
    indices = np.asarray(indices, dtype=np.int64)
    if indices.size == 0:
        return None
    selected = {key: np.asarray(pairs[key])[indices] for key in REQUIRED_PAIR_COLUMNS}
    return stack_s19a_pairs(
        selected, None, None, None, RADIAL_EDGES_MPC, cosmology,
        apply_b_mode_mask=True,
        _pair_lensing_efficiency=cache["beta"][indices],
        _pair_geometry_precomputed=(cache["radius"][indices], cache["phi"][indices]),
        selection_bias_coefficients=selection_bias_coefficients,
    )


def profile_values(lens_result, control_result) -> tuple[np.ndarray | None, np.ndarray]:
    controls = np.zeros(len(RADIAL_EDGES_MPC) - 1, dtype=int)
    if lens_result is None or control_result is None:
        return None, controls
    delta = np.full((len(controls), 2), np.nan, dtype=float)
    for index, (lens_row, control_row) in enumerate(zip(
            lens_result["profile"], control_result["profile"])):
        controls[index] = int(control_row["n_lenses"])
        lens_t = lens_row["delta_sigma_weak_msun_mpc2"]
        lens_x = lens_row["cross_delta_sigma_weak_msun_mpc2"]
        control_t = control_row["delta_sigma_weak_msun_mpc2"]
        control_x = control_row["cross_delta_sigma_weak_msun_mpc2"]
        if None not in (lens_t, lens_x, control_t, control_x):
            delta[index] = (float(lens_t) - float(control_t), float(lens_x) - float(control_x))
    return delta, controls


def covariance_result(replicates: list[np.ndarray], minimum_sets: int,
                      nblocks: int) -> tuple[np.ndarray | None, dict]:
    if len(replicates) != nblocks or nblocks < 3:
        return None, {
            "status": "not_estimable_all_delete_blocks_need_a_finite_common_profile",
            "n_blocks": nblocks,
            "n_finite_replicates": len(replicates),
        }
    matrix = np.asarray(replicates, dtype=float)
    if matrix.shape != (nblocks, EXPECTED_DIMENSION) or not np.all(np.isfinite(matrix)):
        return None, {
            "status": "not_estimable_replicate_matrix_invalid",
            "n_blocks": nblocks,
            "replicate_shape": list(matrix.shape),
        }
    centered = matrix - matrix.mean(axis=0)
    covariance = (nblocks - 1.0) / nblocks * centered.T @ centered
    eigenvalues = np.linalg.eigvalsh(covariance)
    scale = max(float(np.max(np.abs(eigenvalues))), 1.0)
    rank = int(np.linalg.matrix_rank(covariance, tol=scale * 1e-10))
    psd = bool(float(np.min(eigenvalues)) >= -scale * 1e-10)
    min_sets_pass = minimum_sets >= math.ceil(0.8 * EXPECTED_FIXED_LENSES)
    enough_blocks = nblocks >= EXPECTED_DIMENSION + 1
    full_rank = rank == EXPECTED_DIMENSION
    return covariance, {
        "status": "complete_exploratory_delete_one_equal_area_cell_covariance",
        "n_blocks": nblocks,
        "covariance_dimension": EXPECTED_DIMENSION,
        "covariance_rank": rank,
        "rank_deficient": not full_rank,
        "eigenvalues_msun2_mpc4": eigenvalues.tolist(),
        "positive_semidefinite_with_relative_tolerance_1e-10": psd,
        "minimum_common_lens_sets_across_deletions": minimum_sets,
        "minimum_set_support_gate_passed_at_80_percent": min_sets_pass,
        "enough_blocks_for_full_14d_rank": enough_blocks,
        "full_rank_14d": full_rank,
        "exploratory_covariance_quality_gate_passed": bool(psd and min_sets_pass and enough_blocks and full_rank),
        "jackknife_formula": "(B-1)/B sum_b (theta_b - mean(theta)) (theta_b - mean(theta))^T",
    }


def run(args) -> dict:
    output = args.output.resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite existing covariance output: {output}")
    for label, path, expected in (
        ("radial-support summary", args.support_dir / "summary.json", EXPECTED_SUPPORT_SUMMARY_SHA256),
        ("per-lens radial support", args.support_dir / "per_lens_radial_support.tsv", EXPECTED_PER_LENS_SUPPORT_SHA256),
        ("fixed-profile summary", args.fixed_profile_dir / "summary.json", EXPECTED_FIXED_PROFILE_SUMMARY_SHA256),
        ("QA summary", args.qa_dir / "summary.json", EXPECTED_QA_SUMMARY_SHA256),
        ("HSC Y3 n(z)", args.nz_sacc, EXPECTED_NZ_SHA256),
    ):
        observed = sha(path)
        if observed != expected:
            raise ValueError(f"{label} hash mismatch: {observed}")

    qa_summary = read_json(args.qa_dir / "summary.json")
    fixed_summary = read_json(args.fixed_profile_dir / "summary.json")
    if (qa_summary.get("science_measurement") is not False
            or qa_summary.get("network_used") is not False
            or fixed_summary.get("science_measurement") is not False
            or fixed_summary.get("covariance_estimated") is not False):
        raise ValueError("an input artifact unexpectedly promotes the result or contains covariance")

    lens_ids, fixed_selection = fixed_sample_ids(args.support_dir)
    lens_pairs, all_control_pairs, center_to_lens, row_provenance = load_approved_pairs(
        args.plan_dir.resolve(), args.rows.resolve(),
    )
    controls_by_lens: dict[str, list[str]] = {lens_id: [] for lens_id in lens_ids}
    for control_id, matched_lens in center_to_lens.items():
        if matched_lens in lens_ids:
            controls_by_lens[matched_lens].append(control_id)
    for lens_id, controls in controls_by_lens.items():
        if len(controls) != 5:
            raise ValueError(f"fixed lens set {lens_id} does not have exactly five matched controls")
    selected_control_ids = {center for controls in controls_by_lens.values() for center in controls}

    lens_pairs = select_rows(lens_pairs, lens_ids)
    control_pairs = select_rows(all_control_pairs, selected_control_ids)
    del all_control_pairs

    nz_z, nz_bins = load_hsc_y3_nz_sacc(args.nz_sacc)
    cosmology = Cosmology(H0=70.0, Om=0.3)
    lens_cache = cache_geometry_and_efficiency(lens_pairs, nz_z, nz_bins, cosmology)
    control_cache = cache_geometry_and_efficiency(control_pairs, nz_z, nz_bins, cosmology)
    source_reuse = verify_repeated_source_rows(lens_pairs, control_pairs)

    lens_rows_by_set = group_indices(lens_pairs["lens_id"])
    control_rows_by_center = group_indices(control_pairs["lens_id"])
    lens_full = {}
    per_set_rows = []
    per_set_delta = {}
    for lens_id in sorted(lens_ids):
        target_indices = lens_rows_by_set.get(lens_id, np.asarray([], dtype=np.int64))
        random_indices = np.concatenate([
            control_rows_by_center[control_id]
            for control_id in sorted(controls_by_lens[lens_id])
            if control_id in control_rows_by_center
        ]) if any(control_id in control_rows_by_center for control_id in controls_by_lens[lens_id]) else np.asarray([], dtype=np.int64)
        target_result = cached_stack(lens_pairs, lens_cache, target_indices, cosmology)
        random_result = cached_stack(control_pairs, control_cache, random_indices, cosmology)
        delta, n_controls = profile_values(target_result, random_result)
        if delta is None or np.any(~np.isfinite(delta)) or np.any(n_controls < MIN_CONTROLS_PER_ANNULUS):
            raise ValueError(f"full fixed sample does not retain four positive controls in every annulus for {lens_id}")
        lens_full[lens_id] = (target_result, random_result, delta, n_controls,
                              target_indices, random_indices)
        per_set_delta[lens_id] = delta
        for radial_bin, (lo, hi) in enumerate(zip(RADIAL_EDGES_MPC[:-1], RADIAL_EDGES_MPC[1:])):
            target_profile = target_result["profile"][radial_bin]
            random_profile = random_result["profile"][radial_bin]
            per_set_rows.append({
                "matched_lens_id": lens_id,
                "radial_bin": radial_bin,
                "r_min_mpc": float(lo),
                "r_max_mpc": float(hi),
                "n_lens_pairs": int(target_profile["n_pairs"]),
                "n_random_pairs": int(random_profile["n_pairs"]),
                "n_positive_random_centers": int(n_controls[radial_bin]),
                "lens_delta_sigma_msun_mpc2": float(target_profile["delta_sigma_weak_msun_mpc2"]),
                "random_delta_sigma_msun_mpc2": float(random_profile["delta_sigma_weak_msun_mpc2"]),
                "delta_sigma_random_subtracted_msun_mpc2": float(delta[radial_bin, 0]),
                "lens_cross_msun_mpc2": float(target_profile["cross_delta_sigma_weak_msun_mpc2"]),
                "random_cross_msun_mpc2": float(random_profile["cross_delta_sigma_weak_msun_mpc2"]),
                "cross_random_subtracted_msun_mpc2": float(delta[radial_bin, 1]),
            })

    full_delta = np.mean(np.stack(list(per_set_delta.values())), axis=0)
    prior_rows = fixed_summary.get("radial_profile_msun_mpc2", [])
    if len(prior_rows) != len(RADIAL_EDGES_MPC) - 1:
        raise ValueError("pinned fixed-profile summary does not have seven radial bins")
    prior_delta = np.asarray([
        [row["delta_sigma_mean_msun_mpc2"], row["cross_delta_sigma_mean_msun_mpc2"]]
        for row in prior_rows
    ], dtype=float)
    if not np.allclose(full_delta, prior_delta, rtol=2e-9, atol=1e-3):
        raise ValueError("recomputed 52-set full profile differs from the pinned fixed-profile result")

    per_width = []
    replicate_rows, block_support_rows, block_rows, covariance_rows = [], [], [], []
    for width in args.block_widths:
        lens_center_blocks = equal_area_block_ids(
            lens_pairs["lens_ra_deg"], lens_pairs["lens_dec_deg"], cell_width_deg=width,
        )
        lens_source_blocks = equal_area_block_ids(
            lens_pairs["source_ra_deg"], lens_pairs["source_dec_deg"], cell_width_deg=width,
        )
        control_center_blocks = equal_area_block_ids(
            control_pairs["lens_ra_deg"], control_pairs["lens_dec_deg"], cell_width_deg=width,
        )
        control_source_blocks = equal_area_block_ids(
            control_pairs["source_ra_deg"], control_pairs["source_dec_deg"], cell_width_deg=width,
        )
        verify_repeated_cells(
            np.concatenate((lens_pairs["object_id"], control_pairs["object_id"])),
            np.concatenate((lens_source_blocks, control_source_blocks)),
            f"source at {width:g} degrees",
        )
        lens_active = lens_cache["active"]
        control_active = control_cache["active"]
        block_ids = np.unique(np.concatenate((
            lens_center_blocks[lens_active], lens_source_blocks[lens_active],
            control_center_blocks[control_active], control_source_blocks[control_active],
        )))
        if len(block_ids) < 3:
            raise ValueError(f"only {len(block_ids)} occupied {width:g}-degree endpoint cells")

        affected_by_lens = {}
        for lens_id in sorted(lens_ids):
            target_indices = lens_full[lens_id][4]
            random_indices = lens_full[lens_id][5]
            target_active_indices = target_indices[lens_active[target_indices]]
            random_active_indices = random_indices[control_active[random_indices]]
            active_cells = []
            if target_active_indices.size:
                active_cells.extend((
                    lens_center_blocks[target_active_indices],
                    lens_source_blocks[target_active_indices],
                ))
            if random_active_indices.size:
                active_cells.extend((
                    control_center_blocks[random_active_indices],
                    control_source_blocks[random_active_indices],
                ))
            affected_by_lens[lens_id] = set(
                np.unique(np.concatenate(active_cells)).astype(int).tolist()
            ) if active_cells else set()

        block_replicates = []
        n_sets_by_block = []
        min_controls_by_block = []
        print(f"width {width:g} deg: {len(block_ids)} occupied endpoint cells", flush=True)
        for block_index, block_id in enumerate(block_ids, start=1):
            keep_lens = ((lens_center_blocks != block_id) & (lens_source_blocks != block_id))
            keep_control = ((control_center_blocks != block_id) & (control_source_blocks != block_id))
            accepted = {}
            eligible_by_bin = np.zeros(len(RADIAL_EDGES_MPC) - 1, dtype=int)
            valid_pair_sets_by_bin = np.zeros(len(RADIAL_EDGES_MPC) - 1, dtype=int)
            control_ge4_by_bin = np.zeros(len(RADIAL_EDGES_MPC) - 1, dtype=int)
            support_control_counts = []
            deleted_lens_pair_rows = int(np.count_nonzero(~keep_lens))
            deleted_control_pair_rows = int(np.count_nonzero(~keep_control))
            affected_sets = 0

            for lens_id in sorted(lens_ids):
                (full_target, full_random, _, _, target_indices, random_indices) = lens_full[lens_id]
                if int(block_id) in affected_by_lens[lens_id]:
                    affected_sets += 1
                    target_result = cached_stack(
                        lens_pairs, lens_cache, target_indices[keep_lens[target_indices]], cosmology,
                    )
                    random_result = cached_stack(
                        control_pairs, control_cache, random_indices[keep_control[random_indices]], cosmology,
                    )
                else:
                    target_result, random_result = full_target, full_random
                delta, n_controls = profile_values(target_result, random_result)
                if delta is None:
                    continue
                has_values = np.all(np.isfinite(delta), axis=1)
                has_four_controls = n_controls >= MIN_CONTROLS_PER_ANNULUS
                valid_pair_sets_by_bin += has_values.astype(int)
                control_ge4_by_bin += (has_values & has_four_controls).astype(int)
                support_control_counts.extend(n_controls[has_values & has_four_controls].tolist())
                if np.all(has_values & has_four_controls):
                    accepted[lens_id] = delta

            n_common = len(accepted)
            if n_common >= 3:
                values = np.mean(np.stack(list(accepted.values())), axis=0)
                block_replicates.append(values.reshape(-1))
                n_sets_by_block.append(n_common)
                min_controls = int(min(support_control_counts)) if support_control_counts else 0
                min_controls_by_block.append(min_controls)
                status = "complete"
            else:
                values = None
                min_controls = int(min(support_control_counts)) if support_control_counts else 0
                status = "fewer_than_three_common_lens_sets"

            for radial_bin in range(len(RADIAL_EDGES_MPC) - 1):
                block_support_rows.append({
                    "block_width_deg": float(width),
                    "block_id": int(block_id),
                    "radial_bin": radial_bin,
                    "r_min_mpc": float(RADIAL_EDGES_MPC[radial_bin]),
                    "r_max_mpc": float(RADIAL_EDGES_MPC[radial_bin + 1]),
                    "n_lens_sets_with_finite_lens_and_random_profiles": int(valid_pair_sets_by_bin[radial_bin]),
                    "n_lens_sets_with_at_least_four_random_centers": int(control_ge4_by_bin[radial_bin]),
                    "n_lens_sets_common_to_all_seven_annuli": n_common,
                    "minimum_random_centers_among_common_sets": min_controls,
                })
                replicate_rows.append({
                    "block_width_deg": float(width),
                    "block_id": int(block_id),
                    "radial_bin": radial_bin,
                    "r_min_mpc": float(RADIAL_EDGES_MPC[radial_bin]),
                    "r_max_mpc": float(RADIAL_EDGES_MPC[radial_bin + 1]),
                    "replicate_status": status,
                    "n_common_lens_sets": n_common,
                    "delta_sigma_msun_mpc2": float(values[radial_bin, 0]) if values is not None else None,
                    "cross_delta_sigma_msun_mpc2": float(values[radial_bin, 1]) if values is not None else None,
                })
            block_rows.append({
                "block_width_deg": float(width),
                "block_id": int(block_id),
                "affected_lens_sets": affected_sets,
                "deleted_lens_pair_rows": deleted_lens_pair_rows,
                "deleted_control_pair_rows": deleted_control_pair_rows,
                "n_common_lens_sets_all_annuli": n_common,
                "minimum_random_centers_among_common_sets": min_controls,
                "status": status,
            })
            if block_index == 1 or block_index % 5 == 0 or block_index == len(block_ids):
                print(f"  completed {block_index}/{len(block_ids)} cells", flush=True)

        covariance, diagnostics = covariance_result(
            block_replicates,
            min(n_sets_by_block) if n_sets_by_block else 0,
            len(block_ids),
        )
        if covariance is None:
            error_vector = np.full(EXPECTED_DIMENSION, np.nan)
        else:
            error_vector = np.sqrt(np.diag(covariance))
            for i in range(EXPECTED_DIMENSION):
                component_i = "tangential" if i % 2 == 0 else "cross"
                bin_i = i // 2
                for j in range(EXPECTED_DIMENSION):
                    component_j = "tangential" if j % 2 == 0 else "cross"
                    bin_j = j // 2
                    covariance_rows.append({
                        "block_width_deg": float(width),
                        "component_i": component_i,
                        "radial_bin_i": bin_i,
                        "component_j": component_j,
                        "radial_bin_j": bin_j,
                        "covariance_msun2_mpc4": float(covariance[i, j]),
                    })
        profile_rows = []
        for radial_bin in range(len(RADIAL_EDGES_MPC) - 1):
            profile_rows.append({
                "block_width_deg": float(width),
                "radial_bin": radial_bin,
                "r_min_mpc": float(RADIAL_EDGES_MPC[radial_bin]),
                "r_max_mpc": float(RADIAL_EDGES_MPC[radial_bin + 1]),
                "n_full_lens_sets": EXPECTED_FIXED_LENSES,
                "delta_sigma_msun_mpc2": float(full_delta[radial_bin, 0]),
                "delta_sigma_jackknife_error_msun_mpc2": float(error_vector[2 * radial_bin]) if covariance is not None else None,
                "cross_delta_sigma_msun_mpc2": float(full_delta[radial_bin, 1]),
                "cross_jackknife_error_msun_mpc2": float(error_vector[2 * radial_bin + 1]) if covariance is not None else None,
                "n_delete_blocks": len(block_ids),
                "minimum_lens_sets_after_delete": diagnostics.get("minimum_common_lens_sets_across_deletions"),
                "covariance_status": diagnostics["status"],
            })
        per_width.append({
            "block_width_deg": float(width),
            "n_occupied_endpoint_cells": int(len(block_ids)),
            "active_lens_pair_rows": int(np.count_nonzero(lens_active)),
            "active_control_pair_rows": int(np.count_nonzero(control_active)),
            "minimum_common_lens_sets_after_delete": diagnostics.get("minimum_common_lens_sets_across_deletions"),
            "maximum_common_lens_sets_after_delete": max(n_sets_by_block) if n_sets_by_block else 0,
            "covariance": diagnostics,
            "full_profile": profile_rows,
            "covariance_matrix": covariance.tolist() if covariance is not None else None,
        })
    output.mkdir(parents=True, exist_ok=False)
    per_set_fields = list(per_set_rows[0])
    write_tsv(output / "per_lens_set_profile.tsv", per_set_rows, per_set_fields)
    write_tsv(output / "block_deletion_support.tsv", block_rows, list(block_rows[0]))
    write_tsv(output / "block_replicates.tsv", replicate_rows, list(replicate_rows[0]))
    write_tsv(output / "covariance_profile.tsv", [
        row for item in per_width for row in item["full_profile"]
    ], [
        "block_width_deg", "radial_bin", "r_min_mpc", "r_max_mpc", "n_full_lens_sets",
        "delta_sigma_msun_mpc2", "delta_sigma_jackknife_error_msun_mpc2",
        "cross_delta_sigma_msun_mpc2", "cross_jackknife_error_msun_mpc2",
        "n_delete_blocks", "minimum_lens_sets_after_delete", "covariance_status",
    ])
    if covariance_rows:
        write_tsv(output / "covariance_matrix.tsv", covariance_rows, list(covariance_rows[0]))

    covariance_gates = {
        str(item["block_width_deg"]): item["covariance"]["exploratory_covariance_quality_gate_passed"]
        if "exploratory_covariance_quality_gate_passed" in item["covariance"] else False
        for item in per_width
    }
    all_covariances_complete = all(
        item["covariance"]["status"].startswith("complete")
        for item in per_width
    )
    summary = {
        "schema": "astrafek-s23b-wl-fixed-common-support-spatial-covariance-v1",
        "status": (
            "complete_exploratory_fixed_52_paired_endpoint_spatial_covariance"
            if all_covariances_complete
            else "complete_exploratory_fixed_52_endpoint_diagnostics_covariance_gate_failure"
        ),
        "science_measurement": False,
        "network_used": False,
        "cas_query_submitted": False,
        "new_data_downloaded": False,
        "images_read": False,
        "full_photoz_pdfs_included": False,
        "halo_mass_fit_performed": False,
        "scope": "paired endpoint delete-one spatial-cell covariance of the fixed equal-lens 52-set random-subtracted profile",
        "estimator": {
            "source_bins": sorted(SOURCE_BINS),
            "random_control_pooling": "all available matched controls for each lens set pooled by source-pair lensing weights before subtraction",
            "outer_lens_set_weighting": "equal arithmetic mean over lens sets supported in every radial bin after the same delete-one block",
            "minimum_controls_after_deletion": MIN_CONTROLS_PER_ANNULUS,
            "endpoint_rule": "drop a center-source pair if either its lens/random-center endpoint or source endpoint is in the deleted cell",
            "shared_source_rule": "same source IDs and coordinates checked across target/control rows; deterministic equal-area cell applied to both roles",
            "calibration": "existing HSC S19A per-source m/c, radial responsivity, edge-selection response, and pinned HSC Y3 tomographic n(z); recalculated after each affected deletion",
            "radial_edges_proper_mpc": RADIAL_EDGES_MPC.tolist(),
        },
        "sample": {
            **fixed_selection,
            "n_control_centers_selected": len(selected_control_ids),
            "n_lens_source_pair_rows_after_bin_and_calibration_input_selection": int(len(lens_pairs["lens_id"])),
            "n_control_source_pair_rows_after_bin_and_calibration_input_selection": int(len(control_pairs["lens_id"])),
            "source_reuse_audit": source_reuse,
            "full_profile_matches_pinned_fixed_profile": True,
            "full_profile_max_absolute_difference_msun_mpc2": float(np.max(np.abs(full_delta - prior_delta))),
        },
        "cosmology": {"H0_km_s_Mpc": 70.0, "Omega_m": 0.3, "distance_convention": "proper projected radius"},
        "block_method": "equal-area RA/sin(Dec) cells from run_wl_cluster_spatial_jackknife.equal_area_block_ids; widths 2.5 and 5 degrees",
        "block_width_results": per_width,
        "quality_gates": {
            "all_delete_replicates_finite": {str(item["block_width_deg"]): item["covariance"]["status"].startswith("complete") for item in per_width},
            "exploratory_covariance_quality_gate_passed": covariance_gates,
        },
        "limitations": [
            "equal-area RA/sin(Dec) cells are not validated contiguous HSC survey regions and may combine disconnected footprints",
            "block counts and covariance rank limit any 14-dimensional interpretation; rank deficiency is reported",
            "delete blocks can reduce the common set count; every replicate uses the intersection of lens sets with target support and at least four random centers in all seven annuli",
            "fiducial HSC Y3 n(z) uncertainty, photo-z dilution bias, reduced shear, magnification, and calibration uncertainty are not propagated",
            "random-aperture depth, mask, PSF, sky, boost, and cluster-member contamination are not fully matched or corrected",
            "the covariance is exploratory and does not authorize a detection or halo-mass fit",
        ],
        "provenance": {
            "runner_path": str(Path(__file__).resolve()),
            "runner_sha256": sha(Path(__file__).resolve()),
            "approved_row_plan_sha256": row_provenance["plan_sha256"],
            "approved_row_status_sha256": row_provenance["row_status_sha256"],
            "approved_row_summary_sha256": row_provenance["row_summary_sha256"],
            "approved_source_row_batch_sha256": row_provenance["batch_sha256"],
            "radial_support_summary_sha256": sha(args.support_dir / "summary.json"),
            "per_lens_radial_support_sha256": sha(args.support_dir / "per_lens_radial_support.tsv"),
            "fixed_profile_summary_sha256": sha(args.fixed_profile_dir / "summary.json"),
            "qa_summary_sha256": sha(args.qa_dir / "summary.json"),
            "n_z_sacc_sha256": sha(args.nz_sacc),
            "fixed_lens_ids_sha256": fixed_selection["fixed_lens_ids_sha256"],
        },
        "output_directory": str(output),
        "output_sha256": {},
    }
    assessment = [
        "# S23B fixed 52-lens matched-random spatial covariance (2026-11-11)",
        "",
        "The calculation used only the already approved 20 HSC shape-row CSVs, existing HSC Y3 SACC n(z), and pinned Astrafek QA/profile products. It made no CAS request, network query, additional download, image read, or halo-mass fit. `science_measurement=false`.",
        "",
        "The full-data point estimate is the equal arithmetic mean over 52 matched lens sets. Within each set, the target is calibrated separately and the available random-center rows are pooled by source-pair lensing weights before subtraction. The recomputed point profile matches the pinned fixed-support profile within the stated numerical tolerance.",
        "",
        "Each delete-one replicate removes a pair if either endpoint is inside the deleted equal-area RA/sin(Dec) cell. The same common set intersection is used across all seven annuli in that replicate, and each retained set must still have at least four positive random centers per annulus. Both target and random calibrations are recomputed after deletion.",
        "",
        "| Cell width | Occupied cells | Covariance rank / 14 | Minimum common sets | PSD | 14-D quality gate |",
        "|---:|---:|---:|---:|:---:|:---:|",
    ]
    for item in per_width:
        gate = item["covariance"]
        assessment.append(
            f"| {item['block_width_deg']:g}° | {item['n_occupied_endpoint_cells']} | "
            f"{gate.get('covariance_rank', 'n/a')} / 14 | "
            f"{gate.get('minimum_common_lens_sets_across_deletions', 'n/a')} | "
            f"{gate.get('positive_semidefinite_with_relative_tolerance_1e-10', 'n/a')} | "
            f"{gate.get('exploratory_covariance_quality_gate_passed', False)} |"
        )
    assessment.extend([
        "",
        "| R proper (Mpc) | ΔΣt (10¹³ M☉ Mpc⁻²) | Error 2.5° | Error 5° | Cross (10¹³ M☉ Mpc⁻²) | Cross error 2.5° | Cross error 5° |",
        "|---:|---:|---:|---:|---:|---:|---:|",
    ])
    by_width = {item["block_width_deg"]: item["full_profile"] for item in per_width}
    for index, row25 in enumerate(by_width.get(2.5, [])):
        row5 = by_width.get(5.0, [])[index] if len(by_width.get(5.0, [])) > index else {}
        def scaled(value):
            return "n/a" if value is None else f"{float(value) / 1e13:.4f}"
        assessment.append(
            f"| {row25['r_min_mpc']:.1f}–{row25['r_max_mpc']:.1f} | "
            f"{scaled(row25['delta_sigma_msun_mpc2'])} | "
            f"{scaled(row25['delta_sigma_jackknife_error_msun_mpc2'])} | "
            f"{scaled(row5.get('delta_sigma_jackknife_error_msun_mpc2'))} | "
            f"{scaled(row25['cross_delta_sigma_msun_mpc2'])} | "
            f"{scaled(row25['cross_jackknife_error_msun_mpc2'])} | "
            f"{scaled(row5.get('cross_jackknife_error_msun_mpc2'))} |"
        )
    assessment.extend([
        "",
        "A covariance quality gate passes only when all delete blocks are finite, there are at least 15 occupied cells, the 14×14 covariance has full rank and is positive semidefinite to numerical tolerance, and every deletion retains at least 80% of the 52 lens sets. Failure means that this block covariance is a support or footprint diagnostic, not a usable uncertainty matrix.",
        "",
        "Equal-area cells are an approximate fallback, not observed contiguous HSC regions. Even a passed numerical gate does not resolve HSC region modes, random-field matching, n(z) uncertainty, photo-z dilution, reduced shear, boost, or residual member contamination. No significance or mass is inferred.",
        "",
        f"- Fixed lens-set membership SHA-256: `{fixed_selection['fixed_lens_ids_sha256']}`",
        f"- Approved source-row plan SHA-256: `{row_provenance['plan_sha256']}`",
        f"- HSC Y3 n(z) SHA-256: `{EXPECTED_NZ_SHA256}`",
        f"- Runner SHA-256: `{summary['provenance']['runner_sha256']}`",
        f"- Output directory: `{output}`",
        "",
    ])
    (output / "assessment.md").write_text("\n".join(assessment), encoding="utf-8")
    output_names = [
        "assessment.md", "per_lens_set_profile.tsv", "block_deletion_support.tsv",
        "block_replicates.tsv", "covariance_profile.tsv",
    ]
    if covariance_rows:
        output_names.append("covariance_matrix.tsv")
    summary["output_sha256"] = {name: sha(output / name) for name in output_names}
    write_json(output / "summary.json", summary)
    print(json.dumps({
        "status": summary["status"],
        "science_measurement": False,
        "full_profile_matches_pinned_fixed_profile": summary["sample"]["full_profile_matches_pinned_fixed_profile"],
        "block_width_results": [
            {
                "block_width_deg": item["block_width_deg"],
                "occupied_cells": item["n_occupied_endpoint_cells"],
                "covariance_status": item["covariance"]["status"],
                "rank": item["covariance"].get("covariance_rank"),
                "minimum_common_lens_sets": item["covariance"].get("minimum_common_lens_sets_across_deletions"),
                "quality_gate": item["covariance"].get("exploratory_covariance_quality_gate_passed", False),
                "profile": item["full_profile"],
            }
            for item in per_width
        ],
        "output_directory": str(output),
        "summary_sha256": sha(output / "summary.json"),
    }, ensure_ascii=False, indent=2), flush=True)
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan-dir", type=Path, default=DEFAULT_PLAN_DIR)
    parser.add_argument("--rows", type=Path, default=DEFAULT_RESULTS)
    parser.add_argument("--qa-dir", type=Path, default=QA_DIR)
    parser.add_argument("--support-dir", type=Path, default=SUPPORT_DIR)
    parser.add_argument("--fixed-profile-dir", type=Path, default=FIXED_PROFILE_DIR)
    parser.add_argument("--nz-sacc", type=Path, default=NZ_PATH)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--block-widths", type=float, nargs="+", default=[2.5, 5.0])
    args = parser.parse_args()
    if sorted(set(args.block_widths)) != [2.5, 5.0]:
        parser.error("this analysis requires exactly the 2.5- and 5-degree block widths")
    args.block_widths = sorted(args.block_widths)
    run(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
