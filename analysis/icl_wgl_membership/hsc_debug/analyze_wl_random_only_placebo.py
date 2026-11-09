#!/usr/bin/env python3
"""Recompute the matched random-center weak-lensing stack as a placebo null.

This offline diagnostic uses only the already approved HSC S19A shape rows,
the existing HSC Y3 n(z) calibration, and the matched random centers used by
the saved lens-minus-random stack. It tests whether the random-only stack is
consistent with zero under the same approximate spatial jackknife.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.stats import chi2

HERE = Path(__file__).resolve().parent
import sys
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
ROOT = HERE.parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from fetch_ancillary_catalogues import sha
from ogfkit.wlcluster import Cosmology
from run_wl_cluster_spatial_jackknife import equal_area_block_ids
from run_wl_random_shape_stack import verify_and_load_pairs as verify_random_rows
from run_wl_shape_stack import verify_and_load_pairs as verify_lens_rows
from s19a_shape_stack import load_hsc_y3_nz_sacc, stack_s19a_tomographic_pairs

DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
STACK_SUMMARY = DATA / "wl_random_shape_stack_lens_reuse_2026-11-04_01/summary.json"
RANDOM_PLAN = DATA / "wl_random_shape_rows_lens_reuse_2026-11-03_01"
RANDOM_ROWS = RANDOM_PLAN / "results"
LENS_ROWS = DATA / "wl_shape_rows_export_2026-10-18_03"
LENS_PREFLIGHT = DATA / "wl_shape_preflight_2026-10-18_01"
LENS_PLAN = DATA / "wl_shape_rows_plan_2026-10-18_04"
DEFAULT_OUTPUT = DATA / "wl_random_only_placebo_2026-11-06_01"


def profile_vector(result: dict) -> np.ndarray:
    values = []
    for row in result["profile"]:
        tangential = row["delta_sigma_weak_msun_mpc2"]
        cross = row["cross_delta_sigma_weak_msun_mpc2"]
        if tangential is None or cross is None or not np.isfinite(tangential) or not np.isfinite(cross):
            raise ValueError("every random-only radial bin must be populated")
        values.extend((float(tangential), float(cross)))
    return np.asarray(values, dtype=float)


def covariance_and_test(covariance: np.ndarray, values: np.ndarray) -> dict:
    eigenvalues, eigenvectors = np.linalg.eigh(covariance)
    max_eigenvalue = float(np.max(eigenvalues))
    tolerance = max_eigenvalue * 1e-10
    positive = eigenvalues > tolerance
    if not np.any(positive):
        raise ValueError("random-only jackknife covariance has no positive eigenmodes")
    projected = eigenvectors[:, positive].T @ values
    statistic = float(np.sum(projected * projected / eigenvalues[positive]))
    dof = int(np.count_nonzero(positive))
    return {
        "chi_square_vs_zero": statistic,
        "effective_dof_covariance_rank": dof,
        "p_value_chi_square_approx": float(chi2.sf(statistic, dof)),
        "covariance_rank_tolerance": tolerance,
        "minimum_eigenvalue": float(eigenvalues[0]),
        "maximum_eigenvalue": max_eigenvalue,
        "condition_number_positive_subspace": float(max_eigenvalue / eigenvalues[positive][0]),
    }


def run(output: Path, cell_width_deg: float) -> dict:
    output = Path(output).resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite {output}")
    for path in (STACK_SUMMARY, RANDOM_ROWS / "row_status.json",
                 RANDOM_ROWS / "row_summary.json", LENS_ROWS / "row_status.json",
                 LENS_PREFLIGHT / "cas_counts/count_summary.json",
                 LENS_PLAN / "row_export_plan.json"):
        if not path.is_file():
            raise FileNotFoundError(path)
    stack_doc = json.loads(STACK_SUMMARY.read_text(encoding="utf-8"))
    if (stack_doc.get("schema") != "hsc-s19a-random-center-stack-run-v1"
            or stack_doc.get("status") != "exploratory"
            or stack_doc.get("science_measurement") is not False
            or stack_doc.get("random_lens_subtraction_applied") is not True
            or stack_doc.get("spatial_block_covariance_applied") is not True):
        raise ValueError("saved paired random stack is incomplete or has an unexpected gate")
    inp = stack_doc["input"]
    nz_path = Path(inp["n_z_sacc_path"])
    if sha(nz_path) != inp["n_z_sacc_sha256"]:
        raise ValueError("saved HSC Y3 n(z) calibration file changed")
    masked_value = int(inp["b_mode_masked_value"])
    random_pairs, random_provenance = verify_random_rows(
        RANDOM_ROWS, RANDOM_PLAN, masked_value)
    lens_pairs, lens_provenance = verify_lens_rows(
        LENS_ROWS, LENS_PREFLIGHT, LENS_PLAN, masked_value)
    if (random_provenance["random_row_plan_sha256"] != inp["random_row_plan_sha256"]
            or lens_provenance["row_plan_sha256"] != inp["lens_row_plan_sha256"]):
        raise ValueError("source-row inputs differ from the saved paired-stack provenance")

    center_to_lens = random_provenance["random_center_matched_lens_ids"]
    center_counts = random_provenance["random_center_expected_source_counts"]
    matched_lens_ids = set(random_provenance["matched_lens_ids"])
    available_lens_ids = set(np.unique(lens_pairs["lens_id"].astype(str)))
    random_support_lens_ids = {
        center_to_lens[center_id] for center_id, count in center_counts.items()
        if int(count) > 0
    }
    common_lens_ids = matched_lens_ids & available_lens_ids & random_support_lens_ids
    expected_common_count = inp["matched_sample"]["n_common_shape_support_lenses"]
    if len(common_lens_ids) != expected_common_count or expected_common_count != 435:
        raise ValueError(f"common lens support mismatch: {len(common_lens_ids)}")
    center_ids = {center_id for center_id, lens_id in center_to_lens.items()
                  if lens_id in common_lens_ids}
    if len(center_ids) != 2175:
        raise ValueError(f"expected five matched random centers for each common lens; got {len(center_ids)}")
    expected_centers_with_rows = {center_id for center_id in center_ids
                                  if int(center_counts[center_id]) > 0}
    mask = np.isin(random_pairs["lens_id"].astype(str), list(expected_centers_with_rows))
    pairs = {key: np.asarray(value)[mask] for key, value in random_pairs.items()}
    if set(np.unique(pairs["lens_id"].astype(str))) != expected_centers_with_rows:
        raise ValueError("loaded random rows do not match common-support centers")

    z_nz, nz_by_bin = load_hsc_y3_nz_sacc(nz_path)
    edges = [float(value) for value in inp["radial_edges_mpc"]]
    cosmology_doc = inp["cosmology"]
    cosmology = Cosmology(H0=float(cosmology_doc["H0_km_s_Mpc"]),
                          Om=float(cosmology_doc["Omega_m"]))
    full = stack_s19a_tomographic_pairs(pairs, z_nz, nz_by_bin, edges, cosmology)
    saved = stack_doc["matched_random_stack"]
    if full["counts"] != saved["counts"] or len(full["profile"]) != len(saved["profile"]):
        raise ValueError("recomputed random-only stack counts differ from the saved paired result")
    for index, (got, expected) in enumerate(zip(full["profile"], saved["profile"], strict=True)):
        for field in ("delta_sigma_weak_msun_mpc2", "cross_delta_sigma_weak_msun_mpc2"):
            if not np.isclose(got[field], expected[field], rtol=1e-12, atol=1e-8):
                raise ValueError(f"random-only stack does not reproduce saved profile bin {index}/{field}")

    lens_cells = equal_area_block_ids(pairs["lens_ra_deg"], pairs["lens_dec_deg"],
                                      cell_width_deg=cell_width_deg)
    source_cells = equal_area_block_ids(pairs["source_ra_deg"], pairs["source_dec_deg"],
                                        cell_width_deg=cell_width_deg)
    blocks = np.unique(np.concatenate((lens_cells, source_cells)))
    if len(blocks) < 3:
        raise ValueError("at least three occupied spatial blocks are required")
    replicas = []
    deletions = []
    for index, block in enumerate(blocks, start=1):
        keep = (lens_cells != block) & (source_cells != block)
        replica = stack_s19a_tomographic_pairs(
            {key: np.asarray(value)[keep] for key, value in pairs.items()},
            z_nz, nz_by_bin, edges, cosmology)
        replicas.append(profile_vector(replica))
        deletions.append({"block_id": int(block),
                          "deleted_pairs": int(np.count_nonzero(~keep)),
                          "remaining_pairs": int(np.count_nonzero(keep))})
        if index % 5 == 0 or index == len(blocks):
            print(f"random-only jackknife {index}/{len(blocks)} blocks", flush=True)
    replicate_array = np.asarray(replicas, dtype=float)
    centered = replicate_array - replicate_array.mean(axis=0)
    covariance = (len(blocks) - 1.0) / len(blocks) * centered.T @ centered
    full_vector = profile_vector(full)
    nbin = len(full["profile"])
    t_idx = np.arange(0, 2 * nbin, 2)
    x_idx = t_idx + 1
    t_cov = covariance[np.ix_(t_idx, t_idx)]
    x_cov = covariance[np.ix_(x_idx, x_idx)]
    tangential_test = covariance_and_test(t_cov, full_vector[t_idx])
    cross_test = covariance_and_test(x_cov, full_vector[x_idx])
    errors = np.sqrt(np.maximum(np.diag(covariance), 0.0))
    profile = []
    for index, row in enumerate(full["profile"]):
        profile.append({
            "bin": index, "r_min_mpc": row["r_min_mpc"], "r_max_mpc": row["r_max_mpc"],
            "random_only_delta_sigma_msun_mpc2": float(full_vector[2 * index]),
            "random_only_cross_delta_sigma_msun_mpc2": float(full_vector[2 * index + 1]),
            "jackknife_error_delta_sigma_msun_mpc2": float(errors[2 * index]),
            "jackknife_error_cross_msun_mpc2": float(errors[2 * index + 1]),
            "n_pairs": int(row["n_pairs"]), "n_lenses": int(row["n_lenses"]),
        })

    output.mkdir(parents=True, exist_ok=False)
    covariance_path = output / "random_only_covariance.npy"
    np.save(covariance_path, covariance, allow_pickle=False)
    summary = {
        "schema": "astrafek-wl-random-only-placebo-v1",
        "status": "complete_exploratory_diagnostic",
        "science_measurement": False,
        "hsc_network_used": False,
        "image_pixels_read": False,
        "data_requested": False,
        "interpretation": "random-center-only tangential and cross profiles are a placebo/null diagnostic, not a cluster lens signal",
        "n_matched_cluster_lenses": len(common_lens_ids),
        "n_random_centers": len(center_ids),
        "n_random_source_pairs": int(len(pairs["lens_id"])),
        "n_random_pairs_used": int(full["counts"]["used_pairs"]),
        "random_center_mapping_sha256": hashlib.sha256("\n".join(
            f"{key}:{center_to_lens[key]}" for key in sorted(center_ids)).encode()).hexdigest(),
        "estimator": "same HSC Y3 tomographic estimator, source cuts, B-mode mask, n(z), bins, cosmology as saved random-subtracted stack",
        "spatial_jackknife": {
            "method": "delete any random lens-source pair whose lens center or source falls in the same equal-area RA/sin(Dec) cell",
            "cell_width_deg": float(cell_width_deg),
            "n_occupied_blocks": int(len(blocks)),
            "deletions": deletions,
            "covariance_file": str(covariance_path),
            "covariance_sha256": sha(covariance_path),
        },
        "global_null_tests": {"tangential": tangential_test, "cross": cross_test},
        "profile": profile,
        "input_sha256": {
            "paired_stack_summary": sha(STACK_SUMMARY),
            "random_row_summary": sha(RANDOM_ROWS / "row_summary.json"),
            "lens_row_summary": sha(LENS_ROWS / "row_summary.json"),
            "nz_sacc": sha(nz_path),
        },
        "random_row_plan_sha256": random_provenance["random_row_plan_sha256"],
        "lens_row_plan_sha256": lens_provenance["row_plan_sha256"],
        "runner_sha256": sha(Path(__file__).resolve()),
        "limitations": [
            "Random-center residuals are conditional on the selected 435-lens matched-control sample.",
            "Equal-area cells are approximate and do not represent validated contiguous HSC survey regions.",
            "Global chi-square p-values are approximate and do not include photo-z, shear-calibration, reduced-shear, magnification, or source-selection uncertainty.",
            "A null random-only stack does not validate the cluster signal; a non-null stack would flag residual systematics/selection.",
        ],
    }
    summary_path = output / "summary.json"
    summary_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2, allow_nan=False) + "\n",
                            encoding="utf-8")
    summary["summary_sha256"] = sha(summary_path)
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--cell-width-deg", type=float, default=5.0)
    args = parser.parse_args()
    if not 0 < args.cell_width_deg <= 180:
        parser.error("--cell-width-deg must be in (0, 180]")
    result = run(args.output, args.cell_width_deg)
    print(json.dumps({"status": result["status"], "science_measurement": False,
                      "n_random_centers": result["n_random_centers"],
                      "n_random_pairs_used": result["n_random_pairs_used"],
                      "n_occupied_blocks": result["spatial_jackknife"]["n_occupied_blocks"],
                      "global_null_tests": result["global_null_tests"],
                      "summary_sha256": result["summary_sha256"]},
                     ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
