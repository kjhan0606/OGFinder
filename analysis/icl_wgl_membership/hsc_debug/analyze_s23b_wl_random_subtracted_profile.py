#!/usr/bin/env python3
"""Measure an exploratory random-subtracted S23B weak-lensing profile.

This offline Astrafek runner consumes only the previously approved 80-lens /
400-control S19A shape-row export and the pinned HSC Y3 tomographic n(z).
It verifies row receipts and the completed support audit, applies the existing
HSC calibrated estimator, and estimates paired endpoint equal-area jackknife
covariances. It makes no network request, reads no image pixels, and always
records ``science_measurement=false``.
"""
from __future__ import annotations

import argparse
import csv
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

from fetch_ancillary_catalogues import catalogue_reader, sha
from audit_s23b_wl_common_shape_rows import read_export_state
from ogfkit.wlcluster import Cosmology
from run_wl_cluster_spatial_jackknife import equal_area_block_ids
from s19a_shape_stack import (
    decode_hsc_b_mode_mask,
    load_hsc_y3_nz_sacc,
    paired_tomographic_endpoint_jackknife,
    stack_s19a_tomographic_pairs,
)
from submit_s23b_wl_common_shape_rows import (
    APPROVED_PLAN_SHA256,
    DEFAULT_PLAN_DIR,
    DEFAULT_RESULTS,
    LENS_COLUMNS,
    RANDOM_COLUMNS,
    read_tsv,
    validate_plan,
)


DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
DEFAULT_QA = DEFAULT_PLAN_DIR / "quality_audit_2026-11-11_01"
DEFAULT_NZ = DATA / "hsc_y3_nz_sacc_2026-10-19_01/hsc_y3_fourier_space_data_vector.sacc"
DEFAULT_OUTPUT = DATA / "wl_s23b_random_subtracted_profile_2026-11-11_01"
EXPECTED_QA_SUMMARY_SHA256 = "eb73aa0365813d880d3edf60b12e4039f2575f8a97eab8868fc918121a34d445"
EXPECTED_PER_CENTER_SHA256 = "e630442cbaeb054763897a46d68d3b4baea7c0884e80ac22f46ece906583f48e"
EXPECTED_NZ_SHA256 = "a28f9e2e088e92d96d2d87083a6eeeac4c0c84c6b48e033bd3d1dc1bf0d8958f"
EXPECTED_LENS_ROWS = 340_256
EXPECTED_CONTROL_ROWS = 1_572_453
MASKED_RAW_VALUE = 0
RADIAL_EDGES_MPC = np.asarray([0.1, 0.2, 0.4, 0.7, 1.0, 1.5, 2.0, 3.0])
REQUIRED_PAIR_COLUMNS = (
    "lens_id", "object_id", "lens_ra_deg", "lens_dec_deg", "lens_z",
    "source_ra_deg", "source_dec_deg", "e1", "e2", "sigma_e", "rms_e",
    "weight", "m", "c1", "c2", "mag_a10", "resolution", "psf_shape11",
    "psf_shape22", "psf_shape12", "b_mode_pass", "source_zbin",
)
NUMERIC_PAIR_COLUMNS = tuple(
    field for field in REQUIRED_PAIR_COLUMNS
    if field not in {
        "lens_id", "object_id", "lens_ra_deg", "lens_dec_deg", "lens_z",
        "b_mode_pass", "source_zbin",
    }
)


def read_json(path: Path) -> dict:
    return json.loads(Path(path).read_text(encoding="utf-8"))


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


def initialize_pair_columns() -> dict[str, list]:
    return {field: [] for field in REQUIRED_PAIR_COLUMNS}


def append_row(columns: dict[str, list], row: dict[str, str], role: str,
               center: dict, center_to_lens: dict[str, str], batch_name: str) -> None:
    if role == "lens":
        center_id = str(row["lens_id"])
        matched_lens_id = center_id
        center_ra = float(center["ra_deg"])
        center_dec = float(center["dec_deg"])
        center_z = float(center["z_lens"])
        comparisons = (
            (float(row["lens_ra_deg"]), center_ra, 1e-9),
            (float(row["lens_dec_deg"]), center_dec, 1e-9),
            (float(row["lens_z"]), center_z, 1e-10),
        )
        if any(abs(left - right) > tol for left, right, tol in comparisons):
            raise ValueError(f"lens center fields changed in {batch_name}: {center_id}")
    else:
        center_id = str(row["random_center_id"])
        matched_lens_id = str(row["matched_lens_id"])
        if matched_lens_id != str(center["matched_lens_id"]):
            raise ValueError(f"matched lens ID changed in {batch_name}: {center_id}")
        center_ra = float(center["ra_deg"])
        center_dec = float(center["dec_deg"])
        center_z = float(center["z_lens"])
        comparisons = (
            (float(row["random_ra_deg"]), center_ra, 1e-9),
            (float(row["random_dec_deg"]), center_dec, 1e-9),
            (float(row["lens_z"]), center_z, 1e-10),
        )
        if any(abs(left - right) > tol for left, right, tol in comparisons):
            raise ValueError(f"control center fields changed in {batch_name}: {center_id}")

    if not str(row["object_id"]):
        raise ValueError(f"empty source object ID in {batch_name}")
    raw_mask_text = str(row["b_mode_mask_raw"]).strip()
    zbin_text = str(row["source_zbin"]).strip()
    raw_mask, source_zbin = int(raw_mask_text), int(zbin_text)
    if str(raw_mask) != raw_mask_text or raw_mask not in (0, 1):
        raise ValueError(f"non-binary raw B-mode mask in {batch_name}")
    if str(source_zbin) != zbin_text or source_zbin not in (0, 1, 2, 3, 4):
        raise ValueError(f"invalid HSC Y3 source-z bin in {batch_name}")

    columns["lens_id"].append(center_id)
    columns["object_id"].append(str(row["object_id"]))
    columns["lens_ra_deg"].append(center_ra)
    columns["lens_dec_deg"].append(center_dec)
    columns["lens_z"].append(center_z)
    for field in NUMERIC_PAIR_COLUMNS:
        value = float(row[field])
        if not np.isfinite(value):
            raise ValueError(f"non-finite {field} value in {batch_name}")
        columns[field].append(value)
    columns["b_mode_pass"].append(raw_mask)
    columns["source_zbin"].append(source_zbin)
    if role == "control":
        previous = center_to_lens.setdefault(center_id, matched_lens_id)
        if previous != matched_lens_id:
            raise ValueError(f"control center maps to multiple lenses: {center_id}")


def freeze_pairs(columns: dict[str, list]) -> dict[str, np.ndarray]:
    arrays = {}
    for field, values in columns.items():
        if field in {"lens_id", "object_id"}:
            arrays[field] = np.asarray(values, dtype=object)
        elif field == "b_mode_pass":
            arrays[field] = decode_hsc_b_mode_mask(np.asarray(values), masked_value=MASKED_RAW_VALUE)
        elif field == "source_zbin":
            arrays[field] = np.asarray(values, dtype=np.int16)
        else:
            arrays[field] = np.asarray(values, dtype=np.float64)
    return arrays


def load_approved_pairs(plan_dir: Path, rows_dir: Path) -> tuple[dict, dict, dict, dict]:
    plan = validate_plan(plan_dir)
    if plan.get("plan_sha256") != APPROVED_PLAN_SHA256:
        raise ValueError("source-row plan hash differs from the approved plan")
    state, row_summary = read_export_state(rows_dir)
    if (row_summary.get("plan_sha256") != APPROVED_PLAN_SHA256
            or row_summary.get("lens_source_pairs") != EXPECTED_LENS_ROWS
            or row_summary.get("random_source_pairs") != EXPECTED_CONTROL_ROWS
            or row_summary.get("combined_source_pairs") != EXPECTED_LENS_ROWS + EXPECTED_CONTROL_ROWS
            or row_summary.get("science_measurement") is not False
            or row_summary.get("credentials_saved") is not False
            or row_summary.get("images_returned") is not False
            or row_summary.get("photoz_pdf_rows_included") is not False):
        raise ValueError("source-row checkpoint does not match the approved export")
    if len(plan["lens_centers"]) != 80 or len(plan["random_centers"]) != 400:
        raise ValueError("approved plan no longer contains 80 lenses and 400 controls")
    if len(state.get("jobs", {})) != 20:
        raise ValueError("approved source-row export does not have exactly 20 completed batches")

    lens_centers = {str(row["lens_id"]): row for row in plan["lens_centers"]}
    control_centers = {str(row["random_center_id"]): row for row in plan["random_centers"]}
    if len(lens_centers) != 80 or len(control_centers) != 400:
        raise ValueError("duplicate center IDs in approved source-row plan")

    lens_columns = initialize_pair_columns()
    control_columns = initialize_pair_columns()
    center_to_lens: dict[str, str] = {
        center_id: str(center["matched_lens_id"])
        for center_id, center in control_centers.items()
    }
    observed = {"lens": {}, "control": {}}
    batch_hashes = {}
    batch_specs = []
    lens_chunk_size = 25
    for index, batch in enumerate(plan["lens_batches"]):
        centers = plan["lens_centers"][index * lens_chunk_size:(index + 1) * lens_chunk_size]
        batch_specs.append(("lens", batch["name"], centers, LENS_COLUMNS))
    random_count_plan = read_json(plan_dir / "random_shape_count_plan.json")
    random_batch_size = int(random_count_plan["batch_size"])
    for index, batch in enumerate(plan["random_batches"]):
        centers = plan["random_centers"][index * random_batch_size:(index + 1) * random_batch_size]
        batch_specs.append(("control", batch["name"], centers, RANDOM_COLUMNS))

    for index, (role, sql_name, centers, expected_columns) in enumerate(batch_specs, start=1):
        base_name = sql_name.removesuffix(".sql")
        receipt = state["jobs"].get(base_name, {})
        path = rows_dir / f"{base_name}.csv"
        if receipt.get("status") != "complete" or not path.is_file():
            raise ValueError(f"source-row batch is missing or incomplete: {base_name}")
        if path.stat().st_size != receipt.get("bytes") or sha(path) != receipt.get("sha256"):
            raise ValueError(f"source-row batch size/hash differs from checkpoint: {base_name}")
        batch_hashes[base_name] = receipt["sha256"]
        center_lookup = lens_centers if role == "lens" else control_centers
        columns = lens_columns if role == "lens" else control_columns
        count_by_center = {str(row["lens_id"] if role == "lens" else row["random_center_id"]): 0
                           for row in centers}
        center_key = "lens_id" if role == "lens" else "random_center_id"
        with path.open(encoding="utf-8-sig", newline="") as stream:
            reader = catalogue_reader(stream)
            if reader.fieldnames != expected_columns:
                raise ValueError(f"unexpected CSV schema in {base_name}")
            for row in reader:
                if None in row or any(value is None for value in row.values()):
                    raise ValueError(f"malformed CSV row in {base_name}")
                center_id = str(row[center_key])
                if center_id not in count_by_center or center_id not in center_lookup:
                    raise ValueError(f"unrequested center ID in {base_name}: {center_id}")
                count_by_center[center_id] += 1
                append_row(columns, row, role, center_lookup[center_id], center_to_lens, base_name)
        if receipt.get("rows") != sum(count_by_center.values()):
            raise ValueError(f"row count differs from checkpoint in {base_name}")
        for center_id, count in count_by_center.items():
            if role == "lens":
                expected_count = int(plan["lens_counts"][center_id])
            else:
                expected_count = int(plan["random_counts"][center_id]["n_shape_sources"])
            if count != expected_count:
                raise ValueError(f"count-only mismatch for center {center_id}: {count} != {expected_count}")
            observed[role][center_id] = count
        print(f"loaded {index}/20: {base_name} ({receipt['rows']:,} verified rows)", flush=True)

    if sum(observed["lens"].values()) != EXPECTED_LENS_ROWS:
        raise ValueError("lens source-row total changed from the approved count-only total")
    if sum(observed["control"].values()) != EXPECTED_CONTROL_ROWS:
        raise ValueError("control source-row total changed from the approved count-only total")
    lens_pairs = freeze_pairs(lens_columns)
    control_pairs = freeze_pairs(control_columns)
    if len(lens_pairs["lens_id"]) != EXPECTED_LENS_ROWS or len(control_pairs["lens_id"]) != EXPECTED_CONTROL_ROWS:
        raise ValueError("frozen pair arrays have unexpected lengths")
    return lens_pairs, control_pairs, center_to_lens, {
        "plan_sha256": APPROVED_PLAN_SHA256,
        "row_status_sha256": sha(rows_dir / "row_status.json"),
        "row_summary_sha256": sha(rows_dir / "row_summary.json"),
        "batch_sha256": batch_hashes,
        "lens_source_pairs": EXPECTED_LENS_ROWS,
        "control_source_pairs": EXPECTED_CONTROL_ROWS,
    }


def support_policies(qa_dir: Path, center_to_lens: dict[str, str]) -> tuple[dict, dict]:
    summary_path = qa_dir / "summary.json"
    center_path = qa_dir / "per_center.tsv"
    if sha(summary_path) != EXPECTED_QA_SUMMARY_SHA256:
        raise ValueError("the pinned radial/tomographic support audit summary changed")
    if sha(center_path) != EXPECTED_PER_CENTER_SHA256:
        raise ValueError("the pinned per-center support table changed")
    summary = read_json(summary_path)
    if (summary.get("science_measurement") is not False
            or summary.get("network_used") is not False
            or summary.get("rows", {}).get("lens_pairs") != EXPECTED_LENS_ROWS
            or summary.get("rows", {}).get("control_pairs") != EXPECTED_CONTROL_ROWS):
        raise ValueError("support audit does not match the verified source-row export")
    bmode = summary.get("b_mode_polarity_validation", {})
    if (bmode.get("status") != "exact_with_published_GAMA09H_rectangle"
            or bmode.get("masked_raw_value") != MASKED_RAW_VALUE
            or bmode.get("mismatch_rows") != 0):
        raise ValueError("raw B-mode mask polarity is not validated for this export")

    center_rows = read_tsv(center_path)
    if len(center_rows) != 480:
        raise ValueError("support table does not contain exactly 80 lens and 400 control centers")
    lens_support, controls_by_lens = {}, {lens_id: [] for lens_id in set(center_to_lens.values())}
    for row in center_rows:
        center_id = str(row["center_id"])
        role = row["role"]
        positive = int(row["n_calibration_valid"]) > 0
        if role == "lens":
            lens_support[center_id] = positive
        elif role == "control":
            matched = str(row["matched_lens_id"])
            if center_to_lens.get(center_id) != matched:
                raise ValueError(f"support table control-to-lens mapping changed: {center_id}")
            if positive:
                controls_by_lens[matched].append(center_id)
        else:
            raise ValueError(f"unexpected support-table role: {role}")
    if len(lens_support) != 80 or sum(len(values) for values in controls_by_lens.values()) != 376:
        raise ValueError("support-table positive-center totals changed")

    result = {}
    for minimum in (3, 4, 5):
        eligible = {
            lens_id for lens_id, has_lens_support in lens_support.items()
            if has_lens_support and len(controls_by_lens[lens_id]) >= minimum
        }
        control_ids = {
            control_id for lens_id in eligible for control_id in controls_by_lens[lens_id]
        }
        result[minimum] = (eligible, control_ids)
    if len(result[3][0]) != 74 or len(result[4][0]) != 74 or len(result[5][0]) != 71:
        raise ValueError("support-policy membership no longer matches the pinned sensitivity audit")
    if result[3][0] != result[4][0] or result[3][1] != result[4][1]:
        raise ValueError("3/5 and 4/5 positive-support policies no longer select the same sample")
    return result, {
        "n_lens_sets_with_positive_lens_support": sum(lens_support.values()),
        "n_lenses_with_any_positive_control": sum(bool(values) for values in controls_by_lens.values()),
        "n_positive_control_centers": sum(len(values) for values in controls_by_lens.values()),
        "policies": {
            str(minimum): {
                "n_eligible_lens_sets": len(values[0]),
                "n_positive_control_centers": len(values[1]),
                "lens_ids_sha256": hashlib.sha256("\n".join(sorted(values[0])).encode()).hexdigest(),
                "control_ids_sha256": hashlib.sha256("\n".join(sorted(values[1])).encode()).hexdigest(),
            }
            for minimum, values in result.items()
        },
    }


def subset_pairs(pairs: dict[str, np.ndarray], centers: set[str], source_bins: set[int]) -> dict[str, np.ndarray]:
    center_values = np.asarray(pairs["lens_id"]).astype(str)
    bins = np.asarray(pairs["source_zbin"], dtype=int)
    valid_calibration_input = (
        (np.asarray(pairs["m"], dtype=float) > -1.0)
        & ((np.asarray(pairs["psf_shape11"], dtype=float)
            + np.asarray(pairs["psf_shape22"], dtype=float)) > 0.0)
    )
    keep = (
        np.isin(center_values, list(centers))
        & np.isin(bins, list(source_bins))
        & valid_calibration_input
    )
    if not np.any(keep):
        raise ValueError("support and source-bin selection produced an empty pair table")
    return {key: np.asarray(value)[keep] for key, value in pairs.items()}


def wald_diagnostic(profile: list[dict], covariance: np.ndarray, component: str) -> dict:
    if component == "tangential":
        indices = np.arange(0, 2 * len(profile), 2)
        key = "delta_sigma_random_subtracted_msun_mpc2"
    else:
        indices = np.arange(1, 2 * len(profile), 2)
        key = "cross_delta_sigma_random_subtracted_msun_mpc2"
    values = np.asarray([float(row[key]) for row in profile], dtype=float)
    matrix = covariance[np.ix_(indices, indices)]
    matrix = 0.5 * (matrix + matrix.T)
    eig = np.linalg.eigvalsh(matrix)
    scale = max(float(np.max(np.abs(eig))), 1.0)
    rank = int(np.count_nonzero(eig > scale * 1e-10))
    if rank == 0:
        return {"chi2": None, "rank": 0, "nominal_p": None}
    statistic = float(values @ np.linalg.pinv(matrix, rcond=1e-10) @ values)
    return {
        "chi2": statistic,
        "rank": rank,
        "nominal_p": float(chi2.sf(statistic, rank)),
        "covariance_eigenvalue_min": float(np.min(eig)),
        "covariance_eigenvalue_max": float(np.max(eig)),
        "interpretation": "conditional diagnostic under the selected equal-area endpoint jackknife; not a calibrated detection probability",
    }


def run_profile_case(lens_pairs: dict[str, np.ndarray], control_pairs: dict[str, np.ndarray],
                     lens_ids: set[str], control_ids: set[str], source_bins: set[int],
                     nz_z: np.ndarray, nz_bins: np.ndarray, cosmology: Cosmology,
                     block_width_deg: float, label: str) -> dict:
    selected_lens = subset_pairs(lens_pairs, lens_ids, source_bins)
    selected_control = subset_pairs(control_pairs, control_ids, source_bins)
    lens_center_blocks = equal_area_block_ids(
        selected_lens["lens_ra_deg"], selected_lens["lens_dec_deg"],
        cell_width_deg=block_width_deg,
    )
    lens_source_blocks = equal_area_block_ids(
        selected_lens["source_ra_deg"], selected_lens["source_dec_deg"],
        cell_width_deg=block_width_deg,
    )
    control_center_blocks = equal_area_block_ids(
        selected_control["lens_ra_deg"], selected_control["lens_dec_deg"],
        cell_width_deg=block_width_deg,
    )
    control_source_blocks = equal_area_block_ids(
        selected_control["source_ra_deg"], selected_control["source_dec_deg"],
        cell_width_deg=block_width_deg,
    )
    result = paired_tomographic_endpoint_jackknife(
        selected_lens, selected_control,
        lens_center_blocks, lens_source_blocks,
        control_center_blocks, control_source_blocks,
        nz_z, nz_bins, RADIAL_EDGES_MPC, cosmology,
        apply_b_mode_mask=True,
        block_method=(
            f"delete-one {block_width_deg:g}-degree equal-area RA/sin(Dec) cells; "
            "drop pairs whose actual lens/random-center or source endpoint lies in the deleted cell"
        ),
    )
    covariance = np.asarray(result["covariance"], dtype=float)
    profile = result["profile"]
    result["profile"] = [
        {**row, "delta_sigma_1e13_msun_mpc2": row["delta_sigma_random_subtracted_msun_mpc2"] / 1e13,
         "delta_sigma_error_1e13_msun_mpc2": row["delta_sigma_jackknife_error_msun_mpc2"] / 1e13,
         "cross_1e13_msun_mpc2": row["cross_delta_sigma_random_subtracted_msun_mpc2"] / 1e13,
         "cross_error_1e13_msun_mpc2": row["cross_jackknife_error_msun_mpc2"] / 1e13}
        for row in profile
    ]
    result["covariance_rank"] = int(np.linalg.matrix_rank(covariance))
    result["tangential_zero_diagnostic"] = wald_diagnostic(result["profile"], covariance, "tangential")
    result["cross_zero_diagnostic"] = wald_diagnostic(result["profile"], covariance, "cross")
    result["n_selected_lenses"] = len(lens_ids)
    result["n_selected_control_centers"] = len(control_ids)
    result["source_bins"] = sorted(source_bins)
    result["block_width_deg"] = float(block_width_deg)
    result["case_label"] = label
    print(
        f"stacked {label}: lenses={len(lens_ids)}, controls={len(control_ids)}, "
        f"bins={sorted(source_bins)}, blocks={result['n_occupied_blocks']}",
        flush=True,
    )
    return result


def flatten_profiles(profile_cases: list[dict], support_policy: str) -> list[dict]:
    rows = []
    for case in profile_cases:
        for profile in case["profile"]:
            rows.append({
                "support_policy": support_policy,
                "case_label": case["case_label"],
                "source_bins": ",".join(map(str, case["source_bins"])),
                "block_width_deg": case["block_width_deg"],
                "n_occupied_blocks": case["n_occupied_blocks"],
                "radial_bin": profile["bin"],
                "r_min_mpc": profile["r_min_mpc"],
                "r_max_mpc": profile["r_max_mpc"],
                "delta_sigma_1e13_msun_mpc2": profile["delta_sigma_1e13_msun_mpc2"],
                "delta_sigma_error_1e13_msun_mpc2": profile["delta_sigma_error_1e13_msun_mpc2"],
                "cross_1e13_msun_mpc2": profile["cross_1e13_msun_mpc2"],
                "cross_error_1e13_msun_mpc2": profile["cross_error_1e13_msun_mpc2"],
                "n_lens_pairs": profile["n_lens_pairs"],
                "n_control_pairs": profile["n_random_pairs"],
            })
    return rows


def make_plot(path: Path, primary: dict, all_bins: dict, strict: dict) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FuncFormatter

    cases = [
        (primary, "Controls ≥4/5; z-bins 2–4", "#1769aa", "o"),
        (all_bins, "Controls ≥4/5; z-bins 1–4", "#d1495b", "s"),
        (strict, "Controls 5/5; z-bins 2–4", "#2a9d8f", "^"),
    ]
    radius = np.sqrt(RADIAL_EDGES_MPC[:-1] * RADIAL_EDGES_MPC[1:])
    figure, axes = plt.subplots(1, 2, figsize=(12.8, 5.4), sharex=True)
    for case, label, color, marker in cases:
        profile = case["profile"]
        tangential = np.asarray([row["delta_sigma_1e13_msun_mpc2"] for row in profile])
        tangential_error = np.asarray([row["delta_sigma_error_1e13_msun_mpc2"] for row in profile])
        cross = np.asarray([row["cross_1e13_msun_mpc2"] for row in profile])
        cross_error = np.asarray([row["cross_error_1e13_msun_mpc2"] for row in profile])
        axes[0].errorbar(radius, tangential, yerr=tangential_error, marker=marker,
                         color=color, lw=1.4, capsize=2, label=label)
        axes[1].errorbar(radius, cross, yerr=cross_error, marker=marker,
                         color=color, lw=1.4, capsize=2, label=label)
    axes[0].axhline(0.0, color="0.3", lw=0.8)
    axes[1].axhline(0.0, color="0.3", lw=0.8)
    for axis, title in zip(axes, (r"Tangential $\Delta\Sigma_t$", r"Cross $\Delta\Sigma_\times$")):
        axis.set_xscale("log")
        axis.set_xticks([0.1, 0.2, 0.4, 0.7, 1.0, 1.5, 2.0, 3.0])
        axis.get_xaxis().set_major_formatter(FuncFormatter(lambda value, _pos: f"{value:g}"))
        axis.set_xlabel("Projected radius R [proper Mpc]")
        axis.set_ylabel(r"$\Delta\Sigma$ [$10^{13}\,M_\odot\,\mathrm{Mpc}^{-2}$]")
        axis.text(0.03, 0.96, title, transform=axis.transAxes, va="top", ha="left", fontsize=11)
        axis.grid(True, which="major", alpha=0.22)
        axis.legend(loc="upper right", frameon=False, fontsize=8)
    figure.text(
        0.5, 0.015,
        "Exploratory weak-shear estimator; error bars are paired endpoint equal-area jackknife 1σ.",
        ha="center", fontsize=9,
    )
    figure.tight_layout(rect=(0, 0.05, 1, 1))
    figure.savefig(path, dpi=180)
    plt.close(figure)


def write_assessment(path: Path, summary: dict) -> None:
    primary = summary["profile_cases"]["controls_ge4_bins_2_3_4_5deg"]
    all_bins = summary["profile_cases"]["controls_ge4_bins_1_2_3_4_5deg"]
    strict = summary["profile_cases"]["controls_5of5_bins_2_3_4_5deg"]
    rows = [
        "# S23B random-subtracted HSC weak-lensing profile (2026-11-11)",
        "",
        "This Astrafek run used only the previously approved 80-lens/400-control S19A shape-row export and the pinned HSC Y3 SACC tomographic n(z). It made no CAS request, network transfer, or image read. Every result remains exploratory and `science_measurement=false`.",
        "",
        "The primary profile requires at least four positive-support controls per lens set and uses HSC Y3 source bins 2–4, following the existing Astrafek S19A ICL/WL protocol that excludes the lowest source bin after earlier tomography and estimator-weight diagnostics. The 3/5 and 4/5 rules selected exactly the same 74 lens sets in the support audit. The strict sensitivity requires 5/5 controls and retains 71 sets. The all-bin case restores source bin 1 as a sensitivity check.",
        "",
        "The calibrated stack applies the existing HSC per-source multiplicative/additive terms, per-radial-bin responsivity, HSC edge-selection correction, fiducial HSC Y3 n(z), and the verified raw B-mode polarity (raw 0 masked, raw 1 retained). It uses the weak-shear approximation and does not fit halo masses.",
        "",
        "The paired endpoint jackknife deletes every lens/random-center–source pair whose actual center or source endpoint lies in the same equal-area RA/sin(Dec) cell. The 2.5° result is a block-width sensitivity. These cells are not validated contiguous HSC survey regions.",
        "",
        "## Primary profile: 4/5 controls, source bins 2–4, 5° cells",
        "",
        "| R proper (Mpc) | ΔΣt (10¹³ M☉ Mpc⁻²) | Cross (10¹³ M☉ Mpc⁻²) | Lens/control source pairs |",
        "|---:|---:|---:|---:|",
    ]
    for row in primary["profile"]:
        rows.append(
            f"| {row['r_min_mpc']:.1f}–{row['r_max_mpc']:.1f} | "
            f"{row['delta_sigma_1e13_msun_mpc2']:.3f} ± {row['delta_sigma_error_1e13_msun_mpc2']:.3f} | "
            f"{row['cross_1e13_msun_mpc2']:.3f} ± {row['cross_error_1e13_msun_mpc2']:.3f} | "
            f"{row['n_lens_pairs']:,} / {row['n_random_pairs']:,} |"
        )
    rows += [
        "",
        f"This case retains {primary['n_selected_lenses']} lens sets and {primary['n_selected_control_centers']} positive-support control centers across {primary['n_occupied_blocks']} occupied shared endpoint cells. Its nominal tangential zero-profile diagnostic is χ²={primary['tangential_zero_diagnostic']['chi2']:.2f} for rank {primary['tangential_zero_diagnostic']['rank']} (p={primary['tangential_zero_diagnostic']['nominal_p']:.4g}); this is conditional on the approximate block covariance and is not a detection probability.",
        "",
        "## Sensitivity summaries",
        "",
        "| Case | Lens sets | Controls | Cells | ΔΣt first / last bin (10¹³ M☉ Mpc⁻²) | Cross null nominal p |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for name, label, case in (
        ("primary", "4/5, bins 2–4, 5°", primary),
        ("all", "4/5, bins 1–4, 5°", all_bins),
        ("strict", "5/5, bins 2–4, 5°", strict),
        ("width", "4/5, bins 2–4, 2.5°", summary["profile_cases"]["controls_ge4_bins_2_3_4_2p5deg"]),
    ):
        profile = case["profile"]
        cross_p = case["cross_zero_diagnostic"]["nominal_p"]
        cross_p_text = f"{cross_p:.4g}" if cross_p is not None else "undefined"
        rows.append(
            f"| {label} | {case['n_selected_lenses']} | {case['n_selected_control_centers']} | "
            f"{case['n_occupied_blocks']} | "
            f"{profile[0]['delta_sigma_1e13_msun_mpc2']:.3f} / {profile[-1]['delta_sigma_1e13_msun_mpc2']:.3f} | "
            f"{cross_p_text} |"
        )
    rows += [
        "",
        "The source-bin-1 result is reported beside the bins-2–4 profile because the previous HSC tomographic audit found strong low-bin estimator-weight excess and inner-bin inconsistency. Excluding bin 1 changes the inner amplitude materially; the bins-2–4 profile is a prespecified robustness comparison from the established Astrafek workflow, not a post-hoc photo-z correction.",
        "",
        "## Interpretation limits",
        "",
        "This produces a reproducible exploratory random-subtracted profile from the newly approved S23B matched sample. It does not yet establish a calibrated weak-lensing detection or halo-mass estimate. The HSC Y3 n(z) posterior and photo-z dilution uncertainty, reduced shear, magnification, cluster-member contamination, boost correction, random-aperture depth/mask/PSF matching, selection-response uncertainty, and a validated survey-region covariance remain unpropagated. The 2.5°/5° equal-area cells can mix disconnected footprints, and the reported Wald p-values are nominal diagnostics only.",
        "",
        "## Reproducibility",
        "",
        f"- Runner SHA-256: `{summary['provenance']['runner_sha256']}`",
        f"- Approved shape-row plan SHA-256: `{summary['provenance']['row_plan_sha256']}`",
        f"- HSC Y3 SACC n(z) SHA-256: `{summary['provenance']['nz_sacc_sha256']}`",
        f"- Output directory: `{summary['output_directory']}`",
        "- Astrafek runner: `analysis/icl_wgl_membership/hsc_debug/analyze_s23b_wl_random_subtracted_profile.py`",
        "",
    ]
    path.write_text("\n".join(rows), encoding="utf-8")


def run(args) -> dict:
    output = args.output.resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite profile output: {output}")
    plan_dir, rows_dir = args.plan_dir.resolve(), args.rows.resolve()
    qa_dir, nz_path = args.qa_dir.resolve(), args.nz_sacc.resolve()
    if sha(nz_path) != EXPECTED_NZ_SHA256:
        raise ValueError("HSC Y3 SACC n(z) file differs from the pinned artifact")
    qa_summary_path = qa_dir / "summary.json"
    qa_summary = read_json(qa_summary_path)
    if sha(qa_summary_path) != EXPECTED_QA_SUMMARY_SHA256:
        raise ValueError("support audit summary differs from the pinned QA artifact")
    if qa_summary.get("network_used") is not False or qa_summary.get("science_measurement") is not False:
        raise ValueError("support audit has an unexpected scope or science flag")

    lens_pairs, control_pairs, center_to_lens, row_provenance = load_approved_pairs(plan_dir, rows_dir)
    support, support_summary = support_policies(qa_dir, center_to_lens)
    nz_z, nz_by_bin = load_hsc_y3_nz_sacc(nz_path)
    cosmology = Cosmology(H0=70.0, Om=0.3)

    scenarios = [
        ("controls_ge4_bins_2_3_4_5deg", 4, {2, 3, 4}, 5.0),
        ("controls_ge4_bins_1_2_3_4_5deg", 4, {1, 2, 3, 4}, 5.0),
        ("controls_5of5_bins_2_3_4_5deg", 5, {2, 3, 4}, 5.0),
        ("controls_ge4_bins_2_3_4_2p5deg", 4, {2, 3, 4}, 2.5),
    ]
    profile_cases = {}
    for label, minimum_controls, source_bins, width in scenarios:
        lens_ids, control_ids = support[minimum_controls]
        case = run_profile_case(
            lens_pairs, control_pairs, lens_ids, control_ids, source_bins,
            nz_z, nz_by_bin, cosmology, width, label,
        )
        profile_cases[label] = case

    input_validity = {}
    for role, pairs in (("lens", lens_pairs), ("control", control_pairs)):
        invalid_m = np.asarray(pairs["m"], dtype=float) <= -1.0
        invalid_psf_trace = (
            np.asarray(pairs["psf_shape11"], dtype=float)
            + np.asarray(pairs["psf_shape22"], dtype=float)
        ) <= 0.0
        input_validity[role] = {
            "input_center_source_pairs": int(len(invalid_m)),
            "m_at_or_below_minus_one": int(np.count_nonzero(invalid_m)),
            "nonpositive_psf_moment_trace": int(np.count_nonzero(invalid_psf_trace)),
            "either_calibration_input_invalid": int(np.count_nonzero(invalid_m | invalid_psf_trace)),
        }

    summary = {
        "schema": "astrafek-s23b-random-subtracted-wl-profile-v1",
        "status": "complete_exploratory",
        "science_measurement": False,
        "network_used": False,
        "cas_query_submitted": False,
        "new_data_downloaded": False,
        "images_read": False,
        "full_photoz_pdfs_included": False,
        "halo_mass_fit_performed": False,
        "scope": "random-subtracted HSC S19A shape stack around the approved S23B CAMIRA lens/control sample",
        "estimator": {
            "name": "Astrafek stack_s19a_tomographic_pairs with paired_tomographic_endpoint_jackknife",
            "calibration": "per-source m/c, radial-bin HSC responsivity, HSC edge-selection response, HSC Y3 tomographic n(z)",
            "b_mode_mask": "raw value 0 excluded; raw value 1 retained after exact GAMA09H coordinate crosscheck",
            "shear_model": "weak-shear DeltaSigma; no reduced-shear model",
            "control_stack": "pooled positive-support matched random centers; lens and random profiles are recomputed independently in each shared endpoint deletion",
            "primary_source_bins": [2, 3, 4],
            "source_bin_1_sensitivity": "repeat with bins 1–4 because the lowest HSC Y3 source bin has earlier estimator-weight and tomographic-consistency flags",
            "primary_support_rule": "nonzero retained lens support and at least four positive-support controls out of five; 3/5 gives the same 74 lens sets in the pinned support audit",
        },
        "sample": {
            "verified_lens_source_pairs_before_cuts": EXPECTED_LENS_ROWS,
            "verified_control_source_pairs_before_cuts": EXPECTED_CONTROL_ROWS,
            "n_lens_sets_by_support_rule": {
                str(minimum): len(support[minimum][0]) for minimum in (3, 4, 5)
            },
            "n_positive_control_centers_by_support_rule": {
                str(minimum): len(support[minimum][1]) for minimum in (3, 4, 5)
            },
            "support_audit": support_summary,
            "calibration_input_validity": input_validity,
        },
        "cosmology": {"H0_km_s_Mpc": 70.0, "Omega_m": 0.3, "distance_convention": "proper projected radius"},
        "radial_edges_proper_mpc": RADIAL_EDGES_MPC.tolist(),
        "profile_cases": profile_cases,
        "limitations": [
            "equal-area endpoint cells are approximate and are not validated contiguous HSC survey regions",
            "fiducial HSC Y3 n(z) posterior/systematic uncertainty and photo-z dilution calibration are not propagated",
            "reduced shear and magnification are not modeled",
            "random-aperture depth, mask, PSF, sky, boost, and cluster-member contamination are not fully matched or corrected",
            "selection-response calibration uncertainty and correlated survey-region modes remain",
            "the profile is exploratory and its nominal covariance diagnostics are not detection probabilities",
        ],
        "provenance": {
            "runner_path": str(Path(__file__).resolve()),
            "runner_sha256": sha(Path(__file__).resolve()),
            "row_plan_sha256": row_provenance["plan_sha256"],
            "row_status_sha256": row_provenance["row_status_sha256"],
            "row_summary_sha256": row_provenance["row_summary_sha256"],
            "row_batch_sha256": row_provenance["batch_sha256"],
            "qa_summary_sha256": sha(qa_summary_path),
            "per_center_support_sha256": sha(qa_dir / "per_center.tsv"),
            "nz_sacc_path": str(nz_path),
            "nz_sacc_sha256": sha(nz_path),
            "science_measurement": False,
        },
        "output_directory": str(output),
    }

    output.mkdir(parents=True, exist_ok=False)
    profile_rows = flatten_profiles(
        list(profile_cases.values()), "support-policy and source-bin sensitivity",
    )
    write_tsv(output / "profile.tsv", profile_rows, list(profile_rows[0].keys()))
    write_json(output / "summary.json", summary)
    write_assessment(output / "assessment.md", summary)
    make_plot(
        output / "profile.png",
        profile_cases["controls_ge4_bins_2_3_4_5deg"],
        profile_cases["controls_ge4_bins_1_2_3_4_5deg"],
        profile_cases["controls_5of5_bins_2_3_4_5deg"],
    )
    summary["output_sha256"] = {
        name: sha(output / name) for name in ("assessment.md", "profile.tsv", "profile.png")
    }
    write_json(output / "summary.json", summary)
    print(json.dumps({
        "status": summary["status"],
        "science_measurement": False,
        "output_directory": str(output),
        "primary_profile": profile_cases["controls_ge4_bins_2_3_4_5deg"]["profile"],
        "summary_sha256": sha(output / "summary.json"),
    }, ensure_ascii=False, indent=2), flush=True)
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan-dir", type=Path, default=DEFAULT_PLAN_DIR)
    parser.add_argument("--rows", type=Path, default=DEFAULT_RESULTS)
    parser.add_argument("--qa-dir", type=Path, default=DEFAULT_QA)
    parser.add_argument("--nz-sacc", type=Path, default=DEFAULT_NZ)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    run(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
