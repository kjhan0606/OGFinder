#!/usr/bin/env python3
"""Measure descriptive source-shape alignment with the local HSC PSF.

The audit uses only the previously approved S23B lens/control shape-row
export. It selects the profile's 4/5-support sample, verifies every batch,
deduplicates sources within each role, removes sources shared across roles,
and fits a weighted linear relation between calibrated source ellipticity and
PSF ellipticity. Equal-area source cells provide exploratory jackknife errors.
The output remains a diagnostic and always records science_measurement=false.
"""
from __future__ import annotations

import argparse
import hashlib
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

from analyze_s23b_wl_random_subtracted_profile import (
    DEFAULT_NZ,
    DEFAULT_QA,
    DEFAULT_PLAN_DIR,
    DEFAULT_RESULTS,
    EXPECTED_NZ_SHA256,
    load_approved_pairs,
    read_json,
    support_policies,
)
from fetch_ancillary_catalogues import sha
from ogfkit.wlcluster import Cosmology
from run_wl_cluster_spatial_jackknife import equal_area_block_ids
from s19a_shape_stack import (
    load_hsc_y3_nz_sacc,
    pair_geometry,
    tomographic_nz_lensing_efficiency,
)

DEFAULT_OUTPUT = HERE / "wl_s23b_psf_shape_alignment_2026-11-11.json"
SOURCE_BINS = (2, 3, 4)
RADIAL_MIN_MPC = 0.1
RADIAL_MAX_MPC = 3.0
MATCHED_CONTROL_MINIMUM = 4
BLOCK_WIDTHS_DEG = (5.0, 2.5)
MASKED_RAW_VALUE = 0
CONSISTENCY_FIELDS = (
    "source_ra_deg", "source_dec_deg", "e1", "e2", "sigma_e", "rms_e",
    "weight", "m", "c1", "c2", "mag_a10", "resolution",
    "psf_shape11", "psf_shape22", "psf_shape12", "source_zbin",
    "responsivity",
)
PREDICTOR_COUNT = 5
OUTPUT_COUNT = 2


def write_json(path: Path, payload: dict) -> None:
    with Path(path).open("x", encoding="utf-8") as stream:
        json.dump(payload, stream, ensure_ascii=False, indent=2, allow_nan=False)
        stream.write("\n")


def select_centers(pairs: dict[str, np.ndarray], center_ids: set[str],
                   nz_redshift: np.ndarray, nz_by_bin: np.ndarray,
                   cosmology: Cosmology) -> tuple[dict[str, np.ndarray], dict]:
    ids = np.asarray(pairs["lens_id"]).astype(str)
    source_bins = np.asarray(pairs["source_zbin"], dtype=int)
    candidate = np.isin(ids, list(center_ids)) & np.isin(source_bins, SOURCE_BINS)
    if not np.any(candidate):
        raise ValueError("the selected center and source-bin sample is empty")
    selected = {key: np.asarray(value)[candidate] for key, value in pairs.items()}

    radius, _ = pair_geometry(selected, cosmology)
    beta = tomographic_nz_lensing_efficiency(
        selected["lens_z"], selected["source_zbin"],
        nz_redshift, nz_by_bin, cosmology,
    )
    psf_trace = selected["psf_shape11"] + selected["psf_shape22"]
    calibration_valid = (selected["m"] > -1.0) & (psf_trace > 0.0)
    hsc_shape_cuts = (
        (np.hypot(selected["e1"], selected["e2"]) < 2.0)
        & (selected["sigma_e"] >= 0.0)
        & (selected["sigma_e"] <= 0.4)
        & (selected["mag_a10"] <= 25.5)
        & (selected["resolution"] >= 0.3)
        & (selected["weight"] > 0.0)
    )
    unmasked = np.asarray(selected["b_mode_pass"], dtype=bool)
    behind_lens = beta > 0.0
    radial = (radius >= RADIAL_MIN_MPC) & (radius <= RADIAL_MAX_MPC)

    keep = calibration_valid & hsc_shape_cuts & unmasked & behind_lens & radial
    filtered = {key: np.asarray(value)[keep] for key, value in selected.items()}
    if not len(filtered["object_id"]):
        raise ValueError("the common source selection left no rows")
    lens_distance = np.asarray(cosmology.angular_diameter(filtered["lens_z"]), dtype=float)
    geometry_weight = filtered["weight"] * (lens_distance * beta[keep]) ** 2
    responsivity_by_bin = {}
    responsivity = np.full(len(filtered["object_id"]), np.nan, dtype=float)
    for source_bin in SOURCE_BINS:
        in_bin = filtered["source_zbin"] == source_bin
        weight_sum = float(np.sum(geometry_weight[in_bin]))
        if weight_sum <= 0.0:
            raise ValueError(f"source bin {source_bin} has no positive responsivity weight")
        value = float(1.0 - np.sum(
            geometry_weight[in_bin] * filtered["rms_e"][in_bin] ** 2
        ) / weight_sum)
        if not np.isfinite(value) or value <= 0.0:
            raise ValueError(f"source bin {source_bin} has non-positive responsivity")
        responsivity[in_bin] = value
        responsivity_by_bin[str(source_bin)] = value
    filtered["responsivity"] = responsivity
    counts = {
        "candidate_center_sourcebin_pairs": int(len(selected["object_id"])),
        "valid_calibration_input_pairs": int(np.count_nonzero(calibration_valid)),
        "pass_hsc_shape_and_weight_cuts": int(np.count_nonzero(calibration_valid & hsc_shape_cuts)),
        "pass_b_mode_mask": int(np.count_nonzero(calibration_valid & hsc_shape_cuts & unmasked)),
        "positive_nz_lensing_efficiency": int(np.count_nonzero(
            calibration_valid & hsc_shape_cuts & unmasked & behind_lens)),
        "pass_0p1_to_3_proper_mpc": int(np.count_nonzero(keep)),
        "pair_weighted_responsivity_by_source_zbin": responsivity_by_bin,
        "science_measurement": False,
    }
    return filtered, counts


def unique_sources(pairs: dict[str, np.ndarray], role: str) -> dict:
    ids = np.asarray(pairs["object_id"], dtype=object)
    unique_ids, first, inverse, counts = np.unique(
        ids, return_index=True, return_inverse=True, return_counts=True,
    )
    if len(unique_ids) == 0:
        raise ValueError(f"{role} has no unique source IDs")

    for field in CONSISTENCY_FIELDS:
        values = np.asarray(pairs[field])
        reference = values[first][inverse]
        if field == "source_zbin":
            consistent = np.asarray(values == reference)
        else:
            consistent = np.isclose(values, reference, rtol=0.0, atol=1e-9)
        if not np.all(consistent):
            raise ValueError(f"repeated {role} source IDs have inconsistent {field} values")

    return {
        "object_id": unique_ids,
        **{field: np.asarray(pairs[field])[first] for field in CONSISTENCY_FIELDS},
        "pair_rows": int(len(ids)),
        "unique_sources": int(len(unique_ids)),
        "repeated_pair_rows_removed": int(len(ids) - len(unique_ids)),
        "maximum_pairs_per_source": int(np.max(counts)),
    }


def remove_cross_role_overlap(lenses: dict, controls: dict) -> tuple[dict, dict, int]:
    shared = np.intersect1d(lenses["object_id"], controls["object_id"], assume_unique=True)
    if len(shared) == 0:
        return lenses, controls, 0
    lens_keep = ~np.isin(lenses["object_id"], shared, assume_unique=True)
    control_keep = ~np.isin(controls["object_id"], shared, assume_unique=True)

    def subset(data: dict, keep: np.ndarray) -> dict:
        result = {
            key: (np.asarray(value)[keep] if isinstance(value, np.ndarray) else value)
            for key, value in data.items()
        }
        result["unique_sources_after_cross_role_exclusion"] = int(np.count_nonzero(keep))
        return result

    return subset(lenses, lens_keep), subset(controls, control_keep), int(len(shared))


def make_design(data: dict) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    trace = data["psf_shape11"] + data["psf_shape22"]
    if np.any(~np.isfinite(trace)) or np.any(trace <= 0.0):
        raise ValueError("PSF moment traces must be positive and finite")
    psf_e1 = (data["psf_shape11"] - data["psf_shape22"]) / trace
    psf_e2 = 2.0 * data["psf_shape12"] / trace
    zbin = np.asarray(data["source_zbin"], dtype=int)
    x = np.column_stack((zbin == 2, zbin == 3, zbin == 4, psf_e1, psf_e2)).astype(float)
    raw_y = np.column_stack((data["e1"], data["e2"])).astype(float)
    responsivity = data["responsivity"]
    corrected_y = np.column_stack((
        (data["e1"] / (2.0 * responsivity) - data["c1"]) / (1.0 + data["m"]),
        (data["e2"] / (2.0 * responsivity) - data["c2"]) / (1.0 + data["m"]),
    ))
    weights = data["weight"] * (2.0 * responsivity * (1.0 + data["m"])) ** 2
    if (not np.all(np.isfinite(x)) or not np.all(np.isfinite(raw_y))
            or not np.all(np.isfinite(corrected_y)) or not np.all(np.isfinite(weights))
            or np.any(weights <= 0.0)):
        raise ValueError("weighted PSF regression inputs must be finite and positive")
    return x, raw_y, corrected_y, weights


def block_sufficient_statistics(x: np.ndarray, y_raw: np.ndarray, y_corrected: np.ndarray,
                                weights: np.ndarray, block_ids: np.ndarray) -> dict:
    labels, inverse = np.unique(block_ids, return_inverse=True)
    if len(labels) < 3:
        raise ValueError("spatial jackknife requires at least three occupied equal-area cells")
    n_blocks = len(labels)
    xtwx = np.zeros((n_blocks, PREDICTOR_COUNT, PREDICTOR_COUNT), dtype=float)
    xtwy_raw = np.zeros((n_blocks, PREDICTOR_COUNT, OUTPUT_COUNT), dtype=float)
    xtwy_corrected = np.zeros_like(xtwy_raw)

    for i in range(PREDICTOR_COUNT):
        for j in range(i, PREDICTOR_COUNT):
            value = np.bincount(
                inverse, weights=weights * x[:, i] * x[:, j], minlength=n_blocks,
            )
            xtwx[:, i, j] = value
            xtwx[:, j, i] = value
        for output in range(OUTPUT_COUNT):
            xtwy_raw[:, i, output] = np.bincount(
                inverse, weights=weights * x[:, i] * y_raw[:, output], minlength=n_blocks,
            )
            xtwy_corrected[:, i, output] = np.bincount(
                inverse, weights=weights * x[:, i] * y_corrected[:, output], minlength=n_blocks,
            )
    return {
        "labels": labels.astype(int),
        "xtwx": xtwx,
        "xtwy_raw": xtwy_raw,
        "xtwy_corrected": xtwy_corrected,
    }


def fit_coefficients(xtwx: np.ndarray, xtwy: np.ndarray) -> np.ndarray:
    coefficients, _, rank, _ = np.linalg.lstsq(xtwx, xtwy, rcond=1e-12)
    if rank < PREDICTOR_COUNT:
        raise ValueError("PSF regression design lost rank in a spatial jackknife replicate")
    return coefficients[3:5, :].T


def jackknife_error(replicates: np.ndarray) -> np.ndarray:
    mean = np.mean(replicates, axis=0)
    centered = replicates - mean
    return np.sqrt((len(replicates) - 1.0) / len(replicates)
                   * np.sum(centered * centered, axis=0))


def summarize_role(data: dict, *, block_width_deg: float) -> tuple[dict, dict]:
    x, y_raw, y_corrected, weights = make_design(data)
    blocks = equal_area_block_ids(
        data["source_ra_deg"], data["source_dec_deg"],
        cell_width_deg=block_width_deg,
    )
    stats = block_sufficient_statistics(x, y_raw, y_corrected, weights, blocks)
    xtwx_total = np.sum(stats["xtwx"], axis=0)
    xtwy_raw_total = np.sum(stats["xtwy_raw"], axis=0)
    xtwy_corrected_total = np.sum(stats["xtwy_corrected"], axis=0)
    raw_full = fit_coefficients(xtwx_total, xtwy_raw_total)
    corrected_full = fit_coefficients(xtwx_total, xtwy_corrected_total)

    raw_replicates = []
    corrected_replicates = []
    for index in range(len(stats["labels"])):
        matrix = xtwx_total - stats["xtwx"][index]
        raw_replicates.append(fit_coefficients(
            matrix, xtwy_raw_total - stats["xtwy_raw"][index],
        ))
        corrected_replicates.append(fit_coefficients(
            matrix, xtwy_corrected_total - stats["xtwy_corrected"][index],
        ))
    raw_replicates = np.asarray(raw_replicates)
    corrected_replicates = np.asarray(corrected_replicates)

    result = {
        "n_unique_sources": int(len(data["object_id"])),
        "n_occupied_equal_area_cells": int(len(stats["labels"])),
        "source_zbin_counts": {
            str(value): int(np.count_nonzero(data["source_zbin"] == value))
            for value in SOURCE_BINS
        },
        "weighted_least_squares_shape_to_psf_matrix": {
            "raw_catalog_shape": {
                "coefficient_matrix_rows_e1_e2_columns_psf_e1_psf_e2": raw_full.tolist(),
                "equal_area_jackknife_se": jackknife_error(raw_replicates).tolist(),
            },
            "per_object_m_c_corrected_shape": {
                "coefficient_matrix_rows_e1_e2_columns_psf_e1_psf_e2": corrected_full.tolist(),
                "equal_area_jackknife_se": jackknife_error(corrected_replicates).tolist(),
            },
        },
        "block_width_deg": float(block_width_deg),
        "jackknife_method": "delete one equal-area RA/sin(Dec) source cell",
    }
    result["_stats"] = stats
    result["_full_raw"] = raw_full
    result["_full_corrected"] = corrected_full
    return result, stats


def paired_role_difference(lens_result: dict, control_result: dict) -> dict:
    lens_stats = lens_result["_stats"]
    control_stats = control_result["_stats"]
    labels = np.union1d(lens_stats["labels"], control_stats["labels"])

    def totals(stats: dict) -> tuple[np.ndarray, np.ndarray, np.ndarray, dict]:
        return (
            np.sum(stats["xtwx"], axis=0),
            np.sum(stats["xtwy_raw"], axis=0),
            np.sum(stats["xtwy_corrected"], axis=0),
            {int(label): index for index, label in enumerate(stats["labels"])},
        )

    l_m, l_raw, l_corrected, l_index = totals(lens_stats)
    c_m, c_raw, c_corrected, c_index = totals(control_stats)
    full_raw = (fit_coefficients(l_m, l_raw) - fit_coefficients(c_m, c_raw))
    full_corrected = (
        fit_coefficients(l_m, l_corrected) - fit_coefficients(c_m, c_corrected)
    )
    raw_replicates = []
    corrected_replicates = []
    for label in labels:
        li = l_index.get(int(label))
        ci = c_index.get(int(label))
        l_block_m = 0.0 if li is None else lens_stats["xtwx"][li]
        l_block_raw = 0.0 if li is None else lens_stats["xtwy_raw"][li]
        l_block_corrected = 0.0 if li is None else lens_stats["xtwy_corrected"][li]
        c_block_m = 0.0 if ci is None else control_stats["xtwx"][ci]
        c_block_raw = 0.0 if ci is None else control_stats["xtwy_raw"][ci]
        c_block_corrected = 0.0 if ci is None else control_stats["xtwy_corrected"][ci]
        raw_replicates.append(
            fit_coefficients(l_m - l_block_m, l_raw - l_block_raw)
            - fit_coefficients(c_m - c_block_m, c_raw - c_block_raw)
        )
        corrected_replicates.append(
            fit_coefficients(l_m - l_block_m, l_corrected - l_block_corrected)
            - fit_coefficients(c_m - c_block_m, c_corrected - c_block_corrected)
        )
    raw_replicates = np.asarray(raw_replicates)
    corrected_replicates = np.asarray(corrected_replicates)
    return {
        "paired_lens_minus_control_matrix": {
            "raw_catalog_shape": {
                "coefficient_matrix_rows_e1_e2_columns_psf_e1_psf_e2": full_raw.tolist(),
                "equal_area_jackknife_se": jackknife_error(raw_replicates).tolist(),
            },
            "per_object_m_c_corrected_shape": {
                "coefficient_matrix_rows_e1_e2_columns_psf_e1_psf_e2": full_corrected.tolist(),
                "equal_area_jackknife_se": jackknife_error(corrected_replicates).tolist(),
            },
        },
        "paired_occupied_equal_area_cells": int(len(labels)),
    }


def public_role_result(result: dict) -> dict:
    return {key: value for key, value in result.items() if not key.startswith("_")}


def audit(args: argparse.Namespace) -> dict:
    if args.output.exists():
        raise FileExistsError(f"refusing to overwrite output file: {args.output}")
    if sha(args.nz_sacc) != EXPECTED_NZ_SHA256:
        raise ValueError("pinned HSC Y3 n(z) SACC checksum changed")

    lens_pairs, control_pairs, center_to_lens, provenance = load_approved_pairs(
        args.plan_dir, args.rows,
    )
    policies, support_counts = support_policies(args.qa_dir, center_to_lens)
    lens_ids, control_ids = policies[MATCHED_CONTROL_MINIMUM]
    if len(lens_ids) != 74 or len(control_ids) != 367:
        raise ValueError("selected 4/5-support sample differs from the pinned profile sample")

    nz_redshift, nz_by_bin = load_hsc_y3_nz_sacc(args.nz_sacc)
    cosmology = Cosmology(H0=args.h0, Om=args.omega_m)
    selected_lenses, lens_pair_counts = select_centers(
        lens_pairs, lens_ids, nz_redshift, nz_by_bin, cosmology,
    )
    selected_controls, control_pair_counts = select_centers(
        control_pairs, control_ids, nz_redshift, nz_by_bin, cosmology,
    )
    unique_lenses = unique_sources(selected_lenses, "lens")
    unique_controls = unique_sources(selected_controls, "control")
    unique_lenses, unique_controls, shared_source_count = remove_cross_role_overlap(
        unique_lenses, unique_controls,
    )
    if len(unique_lenses["object_id"]) < 100 or len(unique_controls["object_id"]) < 100:
        raise ValueError("cross-role source exclusion left too few unique sources")

    results = {}
    for width in BLOCK_WIDTHS_DEG:
        lens_result, _ = summarize_role(unique_lenses, block_width_deg=width)
        control_result, _ = summarize_role(unique_controls, block_width_deg=width)
        results[f"{width:g}_degree_cells"] = {
            "lenses": public_role_result(lens_result),
            "controls": public_role_result(control_result),
            "paired_difference": paired_role_difference(lens_result, control_result),
        }

    qa_summary = args.qa_dir / "summary.json"
    per_center = args.qa_dir / "per_center.tsv"
    payload = {
        "schema": "astrafek-s23b-wl-psf-shape-alignment-v1",
        "status": "complete_descriptive_psf_alignment_diagnostic",
        "science_measurement": False,
        "network_used": False,
        "additional_data_downloaded": False,
        "image_pixels_read": False,
        "approved_sample": {
            "row_plan_sha256": provenance["plan_sha256"],
            "n_lens_centers": len(lens_ids),
            "n_positive_control_centers": len(control_ids),
            "minimum_positive_controls_per_lens": MATCHED_CONTROL_MINIMUM,
            "support_policy_audit": support_counts,
            "source_redshift_bins": list(SOURCE_BINS),
            "masked_raw_b_mode_value": MASKED_RAW_VALUE,
            "radial_range_proper_mpc": [RADIAL_MIN_MPC, RADIAL_MAX_MPC],
        },
        "pair_selection_counts": {
            "lenses": lens_pair_counts,
            "controls": control_pair_counts,
        },
        "source_deduplication": {
            "lenses": {key: value for key, value in unique_lenses.items()
                       if not isinstance(value, np.ndarray)},
            "controls": {key: value for key, value in unique_controls.items()
                         if not isinstance(value, np.ndarray)},
            "sources_shared_between_lens_and_control_roles_excluded_from_both": shared_source_count,
        },
        "regression": {
            "model": "per-source-z-bin intercept plus linear PSF-ellipticity terms for each galaxy-shape component",
            "psf_ellipticity": "e1=(PSF11-PSF22)/(PSF11+PSF22); e2=2*PSF12/(PSF11+PSF22)",
            "raw_shape_response": "catalog e1/e2",
            "corrected_shape_response": "(e/(2R_z)-c)/(1+m), with pair-weighted sample responsivity R_z per source bin; sample-level edge-selection correction is not applied",
            "weight": "HSC shape_weight*[2R_z*(1+m)]^2, one observation per unique source and role",
            "uncertainty": "delete-one equal-area RA/sin(Dec) source-cell jackknife at 5 and 2.5 degree widths",
            "p_values_reported": False,
        },
        "results": results,
        "inputs": {
            "source_row_plan_sha256": provenance["plan_sha256"],
            "row_status_sha256": provenance["row_status_sha256"],
            "row_summary_sha256": sha(args.rows / "row_summary.json"),
            "source_batch_sha256": provenance["batch_sha256"],
            "support_summary_sha256": sha(qa_summary),
            "per_center_support_sha256": sha(per_center),
            "hsc_y3_nz_sacc_sha256": sha(args.nz_sacc),
        },
        "limitations": [
            "The PSF association is a sample-specific descriptive null and does not isolate PSF leakage from lensing or other spatial correlations.",
            "The regression omits sample-level HSC edge-selection response and responsivity terms.",
            "Equal-area source cells are approximate and are not validated contiguous HSC survey regions.",
            "The audit does not replace the HSC star, random-point, B-mode, or image-level PSF tests.",
            "A small or zero linear coefficient does not establish negligible PSF systematics.",
            "science_measurement remains false.",
        ],
        "analysis_script_sha256": sha(Path(__file__).resolve()),
        "cosmology": {"H0_km_s_mpc": args.h0, "Omega_m": args.omega_m},
    }
    return payload


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan-dir", type=Path, default=DEFAULT_PLAN_DIR)
    parser.add_argument("--rows", type=Path, default=DEFAULT_RESULTS)
    parser.add_argument("--qa-dir", type=Path, default=DEFAULT_QA)
    parser.add_argument("--nz-sacc", type=Path, default=DEFAULT_NZ)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--h0", type=float, default=70.0)
    parser.add_argument("--omega-m", type=float, default=0.3)
    args = parser.parse_args()
    args.output = args.output.resolve()
    result = audit(args)
    write_json(args.output, result)
    print(json.dumps({
        "status": result["status"],
        "science_measurement": result["science_measurement"],
        "selected_lenses": result["approved_sample"]["n_lens_centers"],
        "selected_controls": result["approved_sample"]["n_positive_control_centers"],
        "shared_sources_excluded": result["source_deduplication"][
            "sources_shared_between_lens_and_control_roles_excluded_from_both"],
        "results": result["results"],
        "output": str(args.output),
        "output_sha256": sha(args.output),
    }, ensure_ascii=False, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
