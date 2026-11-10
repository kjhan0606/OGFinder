#!/usr/bin/env python3
"""Audit the approved S23B lens/control shape rows without a lensing stack.

Reads only the already-approved 80-lens/400-control CSV export and the
previously received HSC Y3 SACC n(z) file. It makes no network request, writes
no row-level copy, and always records science_measurement=false.
"""
from __future__ import annotations

import argparse
from array import array
from collections import Counter
import csv
import json
import math
from pathlib import Path
import statistics
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from fetch_ancillary_catalogues import catalogue_reader, sha
from run_wl_shape_stack import Cosmology
from s19a_shape_stack import (
    load_hsc_y3_nz_sacc,
    pair_geometry,
    tomographic_nz_lensing_efficiency,
)
from submit_s23b_wl_common_shape_rows import (
    APPROVED_PLAN_SHA256,
    DEFAULT_PLAN_DIR,
    DEFAULT_RESULTS as DEFAULT_ROWS,
    LENS_COLUMNS,
    RANDOM_COLUMNS,
    validate_plan,
)


DEFAULT_NZ = Path(
    "/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/"
    "hsc_y3_nz_sacc_2026-10-19_01/hsc_y3_fourier_space_data_vector.sacc"
)
DEFAULT_OUTPUT = DEFAULT_ROWS.parent / "quality_audit_2026-11-11_01"
EXPECTED_NZ_SHA256 = "a28f9e2e088e92d96d2d87083a6eeeac4c0c84c6b48e033bd3d1dc1bf0d8958f"
EXPECTED_COUNT_PLAN_SHA256 = "fe19ac426272d4005e3ee74a902dfdb52d1fe83ecb2cc5ec901c10b880d28e8f"
EXPECTED_LENS_COUNT_PLAN_SHA256 = "b713eebded489c280fd9ce1b5840bc8d7d3e42d029035f8dbd414b72fb31e984"
MASKED_RAW_VALUE = 0
GAMA09H_RECTANGLE = (132.5, 140.0, 1.6, 5.0)  # RA min/max, Dec min/max (deg)
RADIAL_EDGES_MPC = np.asarray([0.1, 0.2, 0.4, 0.7, 1.0, 1.5, 2.0, 3.0])
SHAPE_FIELDS = (
    "e_abs", "sigma_e", "rms_e", "weight", "m", "mag_a10", "resolution", "psf_trace",
)
CSV_NUMERIC = (
    "source_ra_deg", "source_dec_deg", "e1", "e2", "sigma_e", "rms_e", "weight",
    "m", "c1", "c2", "mag_a10", "resolution", "psf_shape11", "psf_shape22", "psf_shape12",
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


def effective_n(sum_w: float, sum_w2: float) -> float | None:
    if sum_w <= 0 or sum_w2 <= 0:
        return None
    return float(sum_w * sum_w / sum_w2)


def quantiles(values: list[float], probabilities=(0.16, 0.5, 0.84)) -> dict:
    if not values:
        return {f"q{int(p * 100):02d}": None for p in probabilities}
    result = np.quantile(np.asarray(values, dtype=float), probabilities)
    return {f"q{int(p * 100):02d}": float(value) for p, value in zip(probabilities, result)}


def distribution_summary(values: array) -> dict:
    if not values:
        return {"n": 0, "min": None, "q16": None, "median": None, "q84": None, "max": None}
    data = np.frombuffer(values, dtype=np.float64)
    q16, median, q84 = np.quantile(data, [0.16, 0.5, 0.84])
    return {
        "n": int(len(data)), "min": float(np.min(data)), "q16": float(q16),
        "median": float(median), "q84": float(q84), "max": float(np.max(data)),
    }


def read_export_state(rows_dir: Path) -> tuple[dict, dict]:
    status_path, summary_path = rows_dir / "row_status.json", rows_dir / "row_summary.json"
    if not status_path.is_file() or not summary_path.is_file():
        raise ValueError("completed approved export checkpoint and summary are required")
    state, summary = read_json(status_path), read_json(summary_path)
    if (state.get("status") != "complete" or summary.get("status") != "complete"
            or state.get("plan_sha256") != APPROVED_PLAN_SHA256
            or summary.get("plan_sha256") != APPROVED_PLAN_SHA256
            or state.get("count_only_plan_sha256") != EXPECTED_COUNT_PLAN_SHA256
            or summary.get("count_only_plan_sha256") != EXPECTED_COUNT_PLAN_SHA256):
        raise ValueError("source-row checkpoint is incomplete or belongs to another plan")
    for key in ("credentials_saved", "photoz_pdf_rows_included", "images_returned", "science_measurement"):
        if state.get(key) is not False or summary.get(key) is not False:
            raise ValueError(f"source-row checkpoint has an unexpected scope flag: {key}")
    if len(state.get("jobs", {})) != 20:
        raise ValueError("expected exactly 20 completed source-row batches")
    return state, summary


def new_center_stats(center: dict, role: str) -> dict:
    return {
        "center_id": str(center.get("lens_id", center.get("random_center_id"))),
        "matched_lens_id": str(center.get("lens_id", center.get("matched_lens_id"))),
        "role": role,
        "z_lens": float(center["z_lens"]),
        "n_input": 0,
        "n_shape_quality": 0,
        "n_shape_and_zbin": 0,
        "n_shape_zbin_bmode": 0,
        "n_positive_nz_efficiency": 0,
        "n_in_radial_range": 0,
        "n_calibration_valid": 0,
        "n_bmode_masked_raw0": 0,
        "n_bmode_pass_raw1": 0,
        "zbin_counts": {str(index): 0 for index in range(5)},
        "retained_zbin_counts": {str(index): 0 for index in range(5)},
        "sum_shape_weight": 0.0,
        "sum_shape_weight2": 0.0,
        "sum_sigma_crit_weight_proxy": 0.0,
        "sum_sigma_crit_weight_proxy2": 0.0,
        "n_eff_shape_weight": None,
        "n_eff_sigma_crit_weight_proxy": None,
    }


def load_batch(path: Path, role: str, centers: list[dict], expected_counts: dict[str, int],
               expected_columns: list[str], seen_pairs: dict[str, set[str]],
               global_masks: dict, global_zbins: dict, sample_arrays: dict,
               center_results: dict, radial_results: list[dict],
               radial_tomographic_results: list[dict],
               cosmology: Cosmology, nz_redshift: np.ndarray, nz_by_bin: np.ndarray) -> tuple[int, str]:
    center_key = "lens_id" if role == "lens" else "random_center_id"
    center_index = {str(center[center_key]): index for index, center in enumerate(centers)}
    if len(center_index) != len(centers):
        raise ValueError("duplicate center IDs in one approved SQL batch")
    rows_by_center = np.zeros(len(centers), dtype=np.int64)
    columns = {
        "center_index": [], "object_id": [], "lens_ra_deg": [], "lens_dec_deg": [],
        "lens_z": [], "source_ra_deg": [], "source_dec_deg": [], "e1": [], "e2": [],
        "sigma_e": [], "rms_e": [], "weight": [], "m": [], "c1": [], "c2": [],
        "mag_a10": [], "resolution": [], "psf_shape11": [], "psf_shape22": [],
        "psf_shape12": [], "raw_mask": [], "source_zbin": [],
    }
    with path.open(encoding="utf-8-sig", newline="") as stream:
        reader = catalogue_reader(stream)
        if reader.fieldnames != expected_columns:
            raise ValueError(f"unexpected source-row CSV schema: {path.name}")
        for row in reader:
            if None in row or any(value is None for value in row.values()):
                raise ValueError(f"malformed CSV record in {path.name}")
            center_id = str(row.get(center_key, ""))
            if center_id not in center_index:
                raise ValueError(f"unrequested center ID in {path.name}")
            ci = center_index[center_id]
            center = centers[ci]
            object_id = str(row.get("object_id", ""))
            if not object_id or object_id in seen_pairs[center_id]:
                raise ValueError(f"empty or duplicated center/source pair in {path.name}")
            seen_pairs[center_id].add(object_id)
            rows_by_center[ci] += 1

            if role == "lens":
                lens_ra = float(center["ra_deg"])
                lens_dec = float(center["dec_deg"])
                lens_z = float(center["z_lens"])
                if (str(row["lens_id"]) != center_id
                        or abs(float(row["lens_ra_deg"]) - lens_ra) > 1e-9
                        or abs(float(row["lens_dec_deg"]) - lens_dec) > 1e-9
                        or abs(float(row["lens_z"]) - lens_z) > 1e-10):
                    raise ValueError(f"lens center metadata mismatch in {path.name}")
            else:
                lens_ra = float(center["ra_deg"])
                lens_dec = float(center["dec_deg"])
                lens_z = float(center["z_lens"])
                if (str(row["matched_lens_id"]) != str(center["matched_lens_id"])
                        or abs(float(row["random_ra_deg"]) - lens_ra) > 1e-9
                        or abs(float(row["random_dec_deg"]) - lens_dec) > 1e-9
                        or abs(float(row["lens_z"]) - lens_z) > 1e-10):
                    raise ValueError(f"control center/lens association mismatch in {path.name}")

            numeric = {field: float(row[field]) for field in CSV_NUMERIC}
            if not all(math.isfinite(value) for value in numeric.values()):
                raise ValueError(f"non-finite shape/coordinate value in {path.name}")
            if (not 0 <= numeric["source_ra_deg"] < 360
                    or not -90 <= numeric["source_dec_deg"] <= 90):
                raise ValueError(f"invalid source coordinates in {path.name}")
            mask_text, zbin_text = str(row["b_mode_mask_raw"]).strip(), str(row["source_zbin"]).strip()
            mask, zbin = int(mask_text), int(zbin_text)
            if str(mask) != mask_text or mask not in (0, 1):
                raise ValueError(f"non-binary raw B-mode mask in {path.name}")
            if str(zbin) != zbin_text or not 0 <= zbin <= 4:
                raise ValueError(f"source_zbin outside the documented HSC 0..4 range in {path.name}")
            if numeric["weight"] <= 0 or not 0 <= numeric["rms_e"] < 1:
                raise ValueError(f"invalid HSC weight/rms_e in {path.name}")

            global_masks[role][mask] += 1
            global_zbins[role][zbin] += 1
            values = {
                "e_abs": math.hypot(numeric["e1"], numeric["e2"]),
                "sigma_e": numeric["sigma_e"], "rms_e": numeric["rms_e"],
                "weight": numeric["weight"], "m": numeric["m"],
                "mag_a10": numeric["mag_a10"], "resolution": numeric["resolution"],
                "psf_trace": numeric["psf_shape11"] + numeric["psf_shape22"],
            }
            for field, value in values.items():
                sample_arrays[role]["all"][field].append(value)

            columns["center_index"].append(ci)
            columns["object_id"].append(object_id)
            columns["lens_ra_deg"].append(lens_ra)
            columns["lens_dec_deg"].append(lens_dec)
            columns["lens_z"].append(lens_z)
            for field in CSV_NUMERIC:
                columns[field].append(numeric[field])
            columns["raw_mask"].append(mask)
            columns["source_zbin"].append(zbin)

    if rows_by_center.tolist() != [int(expected_counts[str(center[center_key])]) for center in centers]:
        raise ValueError(f"per-center source rows differ from approved count-only totals: {path.name}")
    a = {key: np.asarray(value) for key, value in columns.items()}
    e_abs = np.hypot(a["e1"], a["e2"])
    shape_quality = (
        (e_abs < 2.0) & (a["sigma_e"] >= 0.0) & (a["sigma_e"] <= 0.4)
        & (a["mag_a10"] <= 25.5) & (a["resolution"] >= 0.3)
    )
    zbin_ok = (a["source_zbin"] >= 1) & (a["source_zbin"] <= 4)
    bmode_ok = a["raw_mask"] != MASKED_RAW_VALUE
    hsc_selected = shape_quality & zbin_ok & bmode_ok
    beta = tomographic_nz_lensing_efficiency(
        a["lens_z"], a["source_zbin"], nz_redshift, nz_by_bin, cosmology,
    )
    radius, _phi = pair_geometry(a, cosmology)
    radial_ok = (radius >= RADIAL_EDGES_MPC[0]) & (radius <= RADIAL_EDGES_MPC[-1])
    beta_ok = beta > 0
    pre_calibration = hsc_selected & beta_ok & radial_ok
    psf_trace = a["psf_shape11"] + a["psf_shape22"]
    calibration_ok = (a["m"] > -1.0) & (psf_trace > 0.0)
    keep = pre_calibration & calibration_ok
    d_l = np.asarray(cosmology.angular_diameter(a["lens_z"]), dtype=float)
    sigma_crit_weight_proxy = a["weight"] * (d_l * beta) ** 2

    n_centers = len(centers)
    ci = a["center_index"].astype(np.int64)
    def count_by_center(mask):
        return np.bincount(ci[mask], minlength=n_centers).astype(np.int64)
    def sum_by_center(values, mask):
        return np.bincount(ci[mask], weights=values[mask], minlength=n_centers)

    n_shape = count_by_center(shape_quality)
    n_shape_z = count_by_center(shape_quality & zbin_ok)
    n_hsc = count_by_center(hsc_selected)
    n_beta = count_by_center(hsc_selected & beta_ok)
    n_radial = count_by_center(hsc_selected & beta_ok & radial_ok)
    n_calibration_invalid = count_by_center(pre_calibration & ~calibration_ok)
    n_keep = count_by_center(keep)
    sum_w = sum_by_center(a["weight"], keep)
    sum_w2 = sum_by_center(a["weight"] ** 2, keep)
    sum_geom = sum_by_center(sigma_crit_weight_proxy, keep)
    sum_geom2 = sum_by_center(sigma_crit_weight_proxy ** 2, keep)
    mask0 = count_by_center(a["raw_mask"] == 0)
    mask1 = count_by_center(a["raw_mask"] == 1)
    zbin_counts = np.vstack([
        count_by_center(a["source_zbin"] == bin_index) for bin_index in range(5)
    ]).T
    retained_zbin_counts = np.vstack([
        count_by_center(keep & (a["source_zbin"] == bin_index)) for bin_index in range(5)
    ]).T

    inside_box = (
        (a["source_ra_deg"] >= GAMA09H_RECTANGLE[0])
        & (a["source_ra_deg"] <= GAMA09H_RECTANGLE[1])
        & (a["source_dec_deg"] >= GAMA09H_RECTANGLE[2])
        & (a["source_dec_deg"] <= GAMA09H_RECTANGLE[3])
    )
    global_masks["rectangle"]["inside_mask0_exact"] += int(np.count_nonzero(inside_box & (a["raw_mask"] == 0)))
    global_masks["rectangle"]["inside_mask1_exact"] += int(np.count_nonzero(inside_box & (a["raw_mask"] == 1)))
    global_masks["rectangle"]["outside_mask0_exact"] += int(np.count_nonzero(~inside_box & (a["raw_mask"] == 0)))
    global_masks["rectangle"]["outside_mask1_exact"] += int(np.count_nonzero(~inside_box & (a["raw_mask"] == 1)))

    for index, center in enumerate(centers):
        center_id = str(center[center_key])
        lens_id = str(center.get("lens_id", center.get("matched_lens_id")))
        stats = center_results[center_id]
        stats.update({
            "n_input": int(rows_by_center[index]),
            "n_shape_quality": int(n_shape[index]),
            "n_shape_and_zbin": int(n_shape_z[index]),
            "n_shape_zbin_bmode": int(n_hsc[index]),
            "n_positive_nz_efficiency": int(n_beta[index]),
            "n_in_radial_range": int(n_radial[index]),
            "n_calibration_invalid_after_cuts": int(n_calibration_invalid[index]),
            "n_calibration_valid": int(n_keep[index]),
            "n_bmode_masked_raw0": int(mask0[index]),
            "n_bmode_pass_raw1": int(mask1[index]),
            "sum_shape_weight": float(sum_w[index]),
            "sum_shape_weight2": float(sum_w2[index]),
            "sum_sigma_crit_weight_proxy": float(sum_geom[index]),
            "sum_sigma_crit_weight_proxy2": float(sum_geom2[index]),
            "n_eff_shape_weight": effective_n(float(sum_w[index]), float(sum_w2[index])),
            "n_eff_sigma_crit_weight_proxy": effective_n(float(sum_geom[index]), float(sum_geom2[index])),
        })
        stats["zbin_counts"] = {str(bin_index): int(zbin_counts[index, bin_index])
                                for bin_index in range(5)}
        stats["retained_zbin_counts"] = {
            str(bin_index): int(retained_zbin_counts[index, bin_index])
            for bin_index in range(5)
        }

    for field, values in (("e_abs", e_abs), ("sigma_e", a["sigma_e"]),
                          ("rms_e", a["rms_e"]), ("weight", a["weight"]),
                          ("m", a["m"]), ("mag_a10", a["mag_a10"]),
                          ("resolution", a["resolution"]), ("psf_trace", psf_trace)):
        sample_arrays[role]["retained"][field].extend(values[keep])

    # Per-center counts and estimator-weight proxies in the established radial bins.
    for radial_index, (lo, hi) in enumerate(zip(RADIAL_EDGES_MPC[:-1], RADIAL_EDGES_MPC[1:])):
        in_annulus = keep & (radius >= lo) & (
            (radius < hi) if radial_index < len(RADIAL_EDGES_MPC) - 2 else (radius <= hi)
        )
        n_annulus = count_by_center(in_annulus)
        sw_annulus = sum_by_center(a["weight"], in_annulus)
        sw2_annulus = sum_by_center(a["weight"] ** 2, in_annulus)
        sg_annulus = sum_by_center(sigma_crit_weight_proxy, in_annulus)
        sg2_annulus = sum_by_center(sigma_crit_weight_proxy ** 2, in_annulus)
        tomo_annulus = {}
        for bin_index in range(1, 5):
            in_annulus_bin = in_annulus & (a["source_zbin"] == bin_index)
            tomo_annulus[bin_index] = (
                count_by_center(in_annulus_bin),
                sum_by_center(a["weight"], in_annulus_bin),
                sum_by_center(a["weight"] ** 2, in_annulus_bin),
                sum_by_center(sigma_crit_weight_proxy, in_annulus_bin),
                sum_by_center(sigma_crit_weight_proxy ** 2, in_annulus_bin),
            )
        for index, center in enumerate(centers):
            center_id = str(center[center_key])
            radial_results.append({
                "role": role, "center_id": center_id, "matched_lens_id": lens_id_for(center, role),
                "radial_bin": radial_index, "r_min_mpc": float(lo), "r_max_mpc": float(hi),
                "n_rows": int(n_annulus[index]),
                "sum_shape_weight": float(sw_annulus[index]),
                "n_eff_shape_weight": effective_n(float(sw_annulus[index]), float(sw2_annulus[index])),
                "sum_sigma_crit_weight_proxy": float(sg_annulus[index]),
                "n_eff_sigma_crit_weight_proxy": effective_n(float(sg_annulus[index]), float(sg2_annulus[index])),
            })
            for bin_index in range(1, 5):
                n_bin, sw_bin, sw2_bin, sg_bin, sg2_bin = tomo_annulus[bin_index]
                radial_tomographic_results.append({
                    "role": role, "center_id": center_id,
                    "matched_lens_id": lens_id_for(center, role),
                    "radial_bin": radial_index, "source_zbin": bin_index,
                    "r_min_mpc": float(lo), "r_max_mpc": float(hi),
                    "n_rows": int(n_bin[index]),
                    "sum_shape_weight": float(sw_bin[index]),
                    "sum_shape_weight2": float(sw2_bin[index]),
                    "n_eff_shape_weight": effective_n(float(sw_bin[index]), float(sw2_bin[index])),
                    "sum_sigma_crit_weight_proxy": float(sg_bin[index]),
                    "sum_sigma_crit_weight_proxy2": float(sg2_bin[index]),
                    "n_eff_sigma_crit_weight_proxy": effective_n(float(sg_bin[index]), float(sg2_bin[index])),
                })
    return int(len(a["center_index"])), sha(path)


def lens_id_for(center: dict, role: str) -> str:
    return str(center["lens_id"] if role == "lens" else center["matched_lens_id"])


def aggregate_role(stats: list[dict]) -> dict:
    fields = (
        "n_input", "n_shape_quality", "n_shape_and_zbin", "n_shape_zbin_bmode",
        "n_positive_nz_efficiency", "n_in_radial_range", "n_calibration_valid",
    )
    out = {field: sum(int(row[field]) for row in stats) for field in fields}
    n_eff_shape = [float(row["n_eff_shape_weight"] or 0.0) for row in stats]
    n_eff_geometry = [float(row["n_eff_sigma_crit_weight_proxy"] or 0.0) for row in stats]
    out.update({
        "n_centers": len(stats),
        "n_zero_after_documented_cuts": sum(row["n_shape_zbin_bmode"] == 0 for row in stats),
        "n_zero_after_geometry_radial_and_calibration_validity": sum(
            row["n_calibration_valid"] == 0 for row in stats
        ),
        "center_n_eff_shape_weight_q16_q50_q84": quantiles(n_eff_shape),
        "center_n_eff_sigma_crit_weight_proxy_q16_q50_q84": quantiles(n_eff_geometry),
        "center_retained_count_q16_q50_q84": quantiles([row["n_calibration_valid"] for row in stats]),
    })
    out["hsc_y3_zbin_retained_after_all_audit_cuts"] = {
        str(index): sum(int(row["retained_zbin_counts"][str(index)]) for row in stats)
        for index in range(5)
    }
    return out


def matched_lens_rows(lens_stats: dict[str, dict], control_stats: dict[str, dict],
                      selected_lens_ids: list[str]) -> tuple[list[dict], dict]:
    controls_by_lens: dict[str, list[dict]] = {lens_id: [] for lens_id in selected_lens_ids}
    for row in control_stats.values():
        controls_by_lens[str(row["matched_lens_id"])].append(row)
    if any(len(rows) != 5 for rows in controls_by_lens.values()):
        raise ValueError("post-cut controls no longer map five-to-one to every selected lens")
    output = []
    count_ratios, weight_ratios, geom_ratios, neff_geometry_ratios = [], [], [], []
    for lens_id in selected_lens_ids:
        lens = lens_stats[lens_id]
        controls = controls_by_lens[lens_id]
        counts = [int(row["n_calibration_valid"]) for row in controls]
        positive_controls = [row for row in controls if int(row["n_calibration_valid"]) > 0]
        sumw = [float(row["sum_shape_weight"]) for row in controls]
        geom = [float(row["sum_sigma_crit_weight_proxy"]) for row in controls]
        neff_geom = [float(row["n_eff_sigma_crit_weight_proxy"] or 0.0) for row in controls]
        lens_count = int(lens["n_calibration_valid"])
        lens_w = float(lens["sum_shape_weight"])
        lens_geom = float(lens["sum_sigma_crit_weight_proxy"])
        lens_neff = float(lens["n_eff_sigma_crit_weight_proxy"] or 0.0)
        mean_control_count = float(np.mean(counts))
        mean_control_w = float(np.mean(sumw))
        mean_control_geom = float(np.mean(geom))
        mean_control_neff = float(np.mean(neff_geom))
        count_ratio = mean_control_count / lens_count if lens_count > 0 else None
        weight_ratio = mean_control_w / lens_w if lens_w > 0 else None
        geom_ratio = mean_control_geom / lens_geom if lens_geom > 0 else None
        neff_ratio = mean_control_neff / lens_neff if lens_neff > 0 else None
        positive_count_mean = (
            float(np.mean([int(row["n_calibration_valid"]) for row in positive_controls]))
            if positive_controls else None
        )
        positive_weight_mean = (
            float(np.mean([float(row["sum_shape_weight"]) for row in positive_controls]))
            if positive_controls else None
        )
        positive_geometry_mean = (
            float(np.mean([float(row["sum_sigma_crit_weight_proxy"]) for row in positive_controls]))
            if positive_controls else None
        )
        if count_ratio is not None:
            count_ratios.append(count_ratio)
        if weight_ratio is not None:
            weight_ratios.append(weight_ratio)
        if geom_ratio is not None:
            geom_ratios.append(geom_ratio)
        if neff_ratio is not None:
            neff_geometry_ratios.append(neff_ratio)
        output.append({
            "lens_id": lens_id, "lens_n_selected": lens_count,
            "control_n_mean": mean_control_count, "control_n_min": min(counts),
            "control_n_median": float(statistics.median(counts)), "control_n_max": max(counts),
            "control_to_lens_count_ratio": count_ratio,
            "lens_sum_shape_weight": lens_w, "control_sum_shape_weight_mean": mean_control_w,
            "control_to_lens_shape_weight_ratio": weight_ratio,
            "lens_n_eff_shape_weight": lens["n_eff_shape_weight"],
            "control_n_eff_shape_weight_mean": float(np.mean([
                float(row["n_eff_shape_weight"] or 0.0) for row in controls
            ])),
            "lens_sigma_crit_weight_proxy": lens_geom,
            "control_sigma_crit_weight_proxy_mean": mean_control_geom,
            "control_to_lens_sigma_crit_weight_ratio": geom_ratio,
            "lens_n_eff_sigma_crit_weight_proxy": lens["n_eff_sigma_crit_weight_proxy"],
            "control_n_eff_sigma_crit_weight_proxy_mean": mean_control_neff,
            "control_to_lens_n_eff_sigma_crit_ratio": neff_ratio,
            "control_n_eff_shape_weight_mean": float(np.mean([
                float(row["n_eff_shape_weight"] or 0.0) for row in controls
            ])),
            "control_to_lens_n_eff_shape_weight_ratio": (
                float(np.mean([float(row["n_eff_shape_weight"] or 0.0) for row in controls]))
                / float(lens["n_eff_shape_weight"])
                if lens["n_eff_shape_weight"] else None
            ),
            "n_zero_control_centers": sum(value == 0 for value in counts),
            "n_positive_control_centers": len(positive_controls),
            "positive_control_count_mean": positive_count_mean,
            "positive_control_to_lens_count_ratio": (
                positive_count_mean / lens_count if lens_count > 0 and positive_count_mean is not None else None
            ),
            "positive_control_to_lens_shape_weight_ratio": (
                positive_weight_mean / lens_w if lens_w > 0 and positive_weight_mean is not None else None
            ),
            "positive_control_to_lens_geometry_proxy_ratio": (
                positive_geometry_mean / lens_geom
                if lens_geom > 0 and positive_geometry_mean is not None else None
            ),
        })
    summary = {
        "lenses_with_zero_retained_sources": sum(row["lens_n_selected"] == 0 for row in output),
        "lenses_with_at_least_one_zero_control": sum(row["n_zero_control_centers"] > 0 for row in output),
        "lenses_with_all_five_zero_controls": sum(row["n_zero_control_centers"] == 5 for row in output),
        "positive_lens_support_pairs": sum(row["lens_n_selected"] > 0 for row in output),
        "mean_control_to_lens_count_ratio_q16_q50_q84": quantiles(count_ratios),
        "mean_control_to_lens_shape_weight_ratio_q16_q50_q84": quantiles(weight_ratios),
        "mean_control_to_lens_sigma_crit_weight_ratio_q16_q50_q84": quantiles(geom_ratios),
        "mean_control_to_lens_sigma_crit_weight_n_eff_ratio_q16_q50_q84": quantiles(neff_geometry_ratios),
        "mean_control_to_lens_shape_weight_n_eff_ratio_q16_q50_q84": quantiles([
            row["control_to_lens_n_eff_shape_weight_ratio"]
            for row in output if row["control_to_lens_n_eff_shape_weight_ratio"] is not None
        ]),
        "ratio_definition": "mean of five matched controls divided by the paired lens, center by center",
    }
    policy_summary = {}
    for minimum_positive_controls in (3, 4, 5):
        eligible = [
            row for row in output
            if row["lens_n_selected"] > 0
            and row["n_positive_control_centers"] >= minimum_positive_controls
        ]
        ratio_sets = {
            "positive_control_count_to_lens_ratio": "positive_control_to_lens_count_ratio",
            "positive_control_shape_weight_to_lens_ratio": "positive_control_to_lens_shape_weight_ratio",
            "positive_control_geometry_proxy_to_lens_ratio": "positive_control_to_lens_geometry_proxy_ratio",
        }
        policy = {
            "minimum_positive_controls": minimum_positive_controls,
            "n_eligible_lens_sets": len(eligible),
            "fraction_of_80_lens_sets": len(eligible) / len(output),
            "lens_pair_rows_retained": sum(int(row["lens_n_selected"]) for row in eligible),
            "control_centers_used": sum(int(row["n_positive_control_centers"]) for row in eligible),
        }
        for summary_key, row_key in ratio_sets.items():
            policy[summary_key + "_q16_q50_q84"] = quantiles([
                float(row[row_key]) for row in eligible if row[row_key] is not None
            ])
        policy_summary[str(minimum_positive_controls)] = policy
    summary["zero_support_policy_sensitivity"] = policy_summary
    return output, summary


def aggregate_radial_tomographic(rows: list[dict]) -> list[dict]:
    totals: dict[tuple[str, int, int], dict] = {}
    for row in rows:
        key = (str(row["role"]), int(row["radial_bin"]), int(row["source_zbin"]))
        aggregate = totals.setdefault(key, {
            "role": key[0], "radial_bin": key[1], "source_zbin": key[2],
            "r_min_mpc": float(row["r_min_mpc"]), "r_max_mpc": float(row["r_max_mpc"]),
            "n_centers_total": 80 if key[0] == "lens" else 400,
            "n_centers_with_support": 0, "n_rows": 0,
            "sum_shape_weight": 0.0, "sum_shape_weight2": 0.0,
            "sum_sigma_crit_weight_proxy": 0.0, "sum_sigma_crit_weight_proxy2": 0.0,
        })
        aggregate["n_centers_with_support"] += int(int(row["n_rows"]) > 0)
        aggregate["n_rows"] += int(row["n_rows"])
        aggregate["sum_shape_weight"] += float(row["sum_shape_weight"])
        aggregate["sum_shape_weight2"] += float(row["sum_shape_weight2"])
        aggregate["sum_sigma_crit_weight_proxy"] += float(row["sum_sigma_crit_weight_proxy"])
        aggregate["sum_sigma_crit_weight_proxy2"] += float(row["sum_sigma_crit_weight_proxy2"])
    output = []
    for key in sorted(totals):
        row = totals[key]
        row["n_centers_without_support"] = row["n_centers_total"] - row["n_centers_with_support"]
        row["n_eff_shape_weight"] = effective_n(row["sum_shape_weight"], row["sum_shape_weight2"])
        row["n_eff_sigma_crit_weight_proxy"] = effective_n(
            row["sum_sigma_crit_weight_proxy"], row["sum_sigma_crit_weight_proxy2"]
        )
        output.append(row)
    return output


def run(args) -> dict:
    plan_dir, rows_dir, output = args.plan_dir.resolve(), args.rows.resolve(), args.output.resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite an existing QA directory: {output}")
    plan = validate_plan(plan_dir)
    state, export_summary = read_export_state(rows_dir)
    if plan["plan_sha256"] != APPROVED_PLAN_SHA256:
        raise ValueError("source-row plan hash differs from the approved export")
    if (plan["manifest"]["count_only_preflight_plan_sha256"] != EXPECTED_COUNT_PLAN_SHA256
            or plan["manifest"]["lens_count_plan_sha256"] != EXPECTED_LENS_COUNT_PLAN_SHA256):
        raise ValueError("one of the count-only inputs has changed from its approved plan")
    if len(plan["lens_centers"]) != 80 or len(plan["random_centers"]) != 400:
        raise ValueError("approved common-support design is not 80 lenses and 400 controls")
    if output == rows_dir or rows_dir in output.parents:
        raise ValueError("QA output must be outside the source-row export directory")

    nz_path = args.nz_sacc.resolve()
    if sha(nz_path) != EXPECTED_NZ_SHA256:
        raise ValueError("local HSC Y3 SACC file differs from the previously approved artifact")
    nz_redshift, nz_by_bin = load_hsc_y3_nz_sacc(nz_path)
    cosmology = Cosmology(H0=70.0, Om=0.3)
    selected_lens_ids = [str(center["lens_id"]) for center in plan["lens_centers"]]
    lens_stats = {center_id: new_center_stats(center, "lens") for center_id, center in (
        (str(center["lens_id"]), center) for center in plan["lens_centers"]
    )}
    control_stats = {center_id: new_center_stats(center, "control") for center_id, center in (
        (str(center["random_center_id"]), center) for center in plan["random_centers"]
    )}
    seen_lens = {center_id: set() for center_id in lens_stats}
    seen_control = {center_id: set() for center_id in control_stats}
    global_masks = {
        "lens": Counter({0: 0, 1: 0}), "control": Counter({0: 0, 1: 0}),
        "rectangle": Counter({
            "inside_mask0_exact": 0, "inside_mask1_exact": 0,
            "outside_mask0_exact": 0, "outside_mask1_exact": 0,
        }),
    }
    global_zbins = {"lens": Counter({i: 0 for i in range(5)}),
                    "control": Counter({i: 0 for i in range(5)})}
    sample_arrays = {
        role: {
            scope: {field: array("d") for field in SHAPE_FIELDS}
            for scope in ("all", "retained")
        }
        for role in ("lens", "control")
    }
    radial_results: list[dict] = []
    radial_tomographic_results: list[dict] = []
    batch_hashes: dict[str, str] = {}
    expected_rows = {"lens": 0, "control": 0}
    row_counts = {"lens": 0, "control": 0}
    batch_specs = []
    for index, batch in enumerate(plan["lens_batches"]):
        chunk = plan["lens_centers"][index * 25:(index + 1) * 25]
        expected = {str(center["lens_id"]): int(plan["lens_counts"][str(center["lens_id"])])
                    for center in chunk}
        batch_specs.append(("lens", batch["name"], chunk, expected, LENS_COLUMNS))
    random_batch_size = int(read_json(plan_dir / "random_shape_count_plan.json")["batch_size"])
    for index, batch in enumerate(plan["random_batches"]):
        chunk = plan["random_centers"][index * random_batch_size:(index + 1) * random_batch_size]
        expected = {center["random_center_id"]: int(plan["random_counts"][center["random_center_id"]]["n_shape_sources"])
                    for center in chunk}
        batch_specs.append(("control", batch["name"], chunk, expected, RANDOM_COLUMNS))

    for index, (role, name, centers, counts, columns) in enumerate(batch_specs, start=1):
        checkpoint_name = name.removesuffix(".sql")
        item = state["jobs"].get(checkpoint_name, {})
        path = rows_dir / f"{checkpoint_name}.csv"
        if item.get("status") != "complete" or not path.is_file():
            raise ValueError(f"source-row batch is missing or incomplete: {checkpoint_name}")
        if path.stat().st_size != item.get("bytes") or sha(path) != item.get("sha256"):
            raise ValueError(f"source-row batch size/hash no longer matches its checkpoint: {checkpoint_name}")
        center_map = lens_stats if role == "lens" else control_stats
        seen = seen_lens if role == "lens" else seen_control
        n_rows, digest = load_batch(
            path, role, centers, counts, columns, seen, global_masks, global_zbins,
            sample_arrays, center_map, radial_results, radial_tomographic_results,
            cosmology, nz_redshift, nz_by_bin,
        )
        if digest != item["sha256"] or n_rows != int(item["rows"]):
            raise ValueError(f"source-row batch row count or digest changed: {checkpoint_name}")
        row_counts[role] += n_rows
        expected_rows[role] += sum(counts.values())
        batch_hashes[checkpoint_name] = digest
        print(f"verified {index}/20: {checkpoint_name} ({n_rows:,} rows)", flush=True)

    if row_counts != {"lens": 340256, "control": 1572453} or row_counts != expected_rows:
        raise ValueError("re-read row totals do not match the approved count-only totals")
    lens_rows = list(lens_stats.values())
    control_rows = list(control_stats.values())
    lens_aggregate = aggregate_role(lens_rows)
    control_aggregate = aggregate_role(control_rows)
    matched_rows, matched_summary = matched_lens_rows(lens_stats, control_stats, selected_lens_ids)
    radial_tomographic_summary = aggregate_radial_tomographic(radial_tomographic_results)

    geometry = global_masks["rectangle"]
    mask_geometry_mismatches = geometry["outside_mask0_exact"] + geometry["inside_mask1_exact"]
    prior_mask_audit = HERE / "wl_bmode_coordinate_validation_2026-10-19.md"
    distribution_summaries = {}
    for role in ("lens", "control"):
        distribution_summaries[role] = {
            scope: {field: distribution_summary(values) for field, values in fields.items()}
            for scope, fields in sample_arrays[role].items()
        }

    if mask_geometry_mismatches:
        mask_geometry_status = "requires_review_nonzero_coordinate_classification_mismatches"
    else:
        mask_geometry_status = "exact_with_published_GAMA09H_rectangle"

    audit_summary = {
        "schema": "astrafek-s23b-wl-common-shape-row-qa-v2",
        "status": "complete_local_shape_quality_radial_tomographic_and_common_support_audit",
        "science_measurement": False,
        "network_used": False,
        "cas_query_submitted": False,
        "additional_data_downloaded": False,
        "credentials_saved": False,
        "images_included": False,
        "full_photoz_pdfs_included": False,
        "row_plan_sha256": APPROVED_PLAN_SHA256,
        "count_only_plan_sha256": EXPECTED_COUNT_PLAN_SHA256,
        "lens_count_plan_sha256": EXPECTED_LENS_COUNT_PLAN_SHA256,
        "n_z_sacc_path": str(nz_path), "n_z_sacc_sha256": EXPECTED_NZ_SHA256,
        "official_hsc_shape_catalog_documentation": "https://hsc-release.mtk.nao.ac.jp/doc/index.php/s19a-shape-catalog-pdr3/",
        "prior_b_mode_polarity_audit_path": str(prior_mask_audit),
        "prior_b_mode_polarity_audit_sha256": sha(prior_mask_audit),
        "source_row_summary_sha256": sha(rows_dir / "row_summary.json"),
        "source_row_status_sha256": sha(rows_dir / "row_status.json"),
        "source_row_batch_sha256": batch_hashes,
        "rows": {
            "lens_pairs": row_counts["lens"], "control_pairs": row_counts["control"],
            "combined_center_source_pairs": sum(row_counts.values()),
            "n_lenses": len(lens_stats), "n_control_centers": len(control_stats),
        },
        "source_selection_contract": {
            "shape": "hypot(e1,e2)<2 and 0<=sigma_e<=0.4",
            "magnitude": "mag_a10<=25.5",
            "resolution": "resolution>=0.3",
            "tomographic": "source_zbin 1..4; bin 0 excluded",
            "b_mode": "exclude stored raw value 0; pass raw value 1",
            "lensing_efficiency": "fiducial HSC Y3 SACC n(z); retain beta=<D_ls/D_s>>0",
            "radial": "0.1<=R_proper<=3.0 Mpc",
            "calibration_input_validity": "m>-1 and psf_shape11+psf_shape22>0; diagnostic gate, not HSC shape cut",
            "identical_for_lenses_and_controls": True,
        },
        "shape_cutflow": {"lens": lens_aggregate, "control": control_aggregate},
        "hsc_raw_mask_counts": {
            role: {str(mask): int(global_masks[role][mask]) for mask in (0, 1)}
            for role in ("lens", "control")
        },
        "hsc_zbin_raw_counts": {
            role: {str(bin_index): int(global_zbins[role][bin_index]) for bin_index in range(5)}
            for role in ("lens", "control")
        },
        "hsc_y3_zbin_retained_after_all_audit_cuts": {
            "lens": lens_aggregate["hsc_y3_zbin_retained_after_all_audit_cuts"],
            "control": control_aggregate["hsc_y3_zbin_retained_after_all_audit_cuts"],
        },
        "b_mode_polarity_validation": {
            "official_HSC_page_defines_mask_purpose_but_not_raw_polarity": True,
            "masked_raw_value": MASKED_RAW_VALUE,
            "coordinate_crosscheck_rectangle_deg": {
                "ra": [GAMA09H_RECTANGLE[0], GAMA09H_RECTANGLE[1]],
                "dec": [GAMA09H_RECTANGLE[2], GAMA09H_RECTANGLE[3]],
            },
            "coordinate_crosscheck_counts": {key: int(geometry[key]) for key in (
                "inside_mask0_exact", "inside_mask1_exact", "outside_mask0_exact", "outside_mask1_exact",
            )},
            "mismatch_rows": int(mask_geometry_mismatches),
            "status": mask_geometry_status,
            "prior_same_release_audit": str(prior_mask_audit),
        },
        "weight_effective_support": {
            "shape_weight_definition": "N_eff=(sum HSC derived weight)^2/sum(weight^2)",
            "sigma_crit_weight_proxy_definition": "HSC weight*(D_l*<D_ls/D_s>)^2; common 4*pi*G/c^2 factor cancels N_eff",
            "calibration_note": "does not include per-object (1+m), HSC responsivity, selection-bias response, covariance, or n(z) uncertainty; not final estimator weight",
            "lens": distribution_summaries["lens"],
            "control": distribution_summaries["control"],
        },
        "matched_control_support": matched_summary,
        "zero_support_policy_sensitivity": matched_summary["zero_support_policy_sensitivity"],
        "radial_tomographic_support": radial_tomographic_summary,
        "limitations": [
            "Rows are center/source pairs; sources repeated across apertures are intentionally retained and not independent.",
            "The HSC Y3 SACC curves are fiducial n(z); their systematic/posterior uncertainty is not propagated.",
            "No shear profile, random-subtracted estimator, mass fit, PSF spatial model, or covariance estimate was calculated.",
            "The five controls are matched positions, not a validated uncontaminated environmental null; prior CAMIRA and footprint limits remain.",
        ],
    }

    output.mkdir(parents=True, exist_ok=False)
    write_json(output / "summary.json", audit_summary)
    center_fields = [
        "center_id", "matched_lens_id", "role", "z_lens", "n_input", "n_shape_quality",
        "n_shape_and_zbin", "n_shape_zbin_bmode", "n_positive_nz_efficiency", "n_in_radial_range",
        "n_calibration_invalid_after_cuts", "n_calibration_valid", "n_bmode_masked_raw0",
        "n_bmode_pass_raw1", "sum_shape_weight", "sum_shape_weight2", "n_eff_shape_weight",
        "sum_sigma_crit_weight_proxy", "sum_sigma_crit_weight_proxy2",
        "n_eff_sigma_crit_weight_proxy",
        "zbin_counts", "retained_zbin_counts",
    ]
    write_tsv(output / "per_center.tsv", [
        {
            **row,
            "zbin_counts": json.dumps(row["zbin_counts"], sort_keys=True),
            "retained_zbin_counts": json.dumps(row["retained_zbin_counts"], sort_keys=True),
        }
        for row in [*lens_rows, *control_rows]
    ], center_fields)
    match_fields = list(matched_rows[0].keys())
    write_tsv(output / "per_lens_five_control_support.tsv", matched_rows, match_fields)
    radial_fields = list(radial_results[0].keys())
    write_tsv(output / "radial_support.tsv", radial_results, radial_fields)
    radial_tomographic_fields = list(radial_tomographic_results[0].keys())
    write_tsv(
        output / "per_center_radial_tomographic_support.tsv",
        radial_tomographic_results, radial_tomographic_fields,
    )
    aggregate_fields = list(radial_tomographic_summary[0].keys())
    write_tsv(
        output / "radial_tomographic_support.tsv",
        radial_tomographic_summary, aggregate_fields,
    )
    policy_rows = []
    for policy in matched_summary["zero_support_policy_sensitivity"].values():
        flat = {key: value for key, value in policy.items() if not key.endswith("_q16_q50_q84")}
        for key, value in policy.items():
            if key.endswith("_q16_q50_q84"):
                for quantile in ("q16", "q50", "q84"):
                    flat[f"{key.removesuffix('_q16_q50_q84')}_{quantile}"] = value[quantile]
        policy_rows.append(flat)
    policy_fields = list(policy_rows[0].keys())
    write_tsv(output / "zero_support_policy_sensitivity.tsv", policy_rows, policy_fields)
    write_assessment(output / "assessment.md", audit_summary, output)
    return audit_summary


def write_assessment(path: Path, summary: dict, output: Path) -> None:
    lens, control = summary["shape_cutflow"]["lens"], summary["shape_cutflow"]["control"]
    bmode = summary["b_mode_polarity_validation"]
    matched = summary["matched_control_support"]
    lens_retained = summary["weight_effective_support"]["lens"]["retained"]
    control_retained = summary["weight_effective_support"]["control"]["retained"]
    def fmt_quantiles(values: dict) -> str:
        return "/".join(f"{values[key]:.3f}" for key in ("q16", "q50", "q84"))
    def fmt_distribution(values: dict) -> str:
        return "/".join(f"{values[key]:.3f}" for key in ("q16", "median", "q84"))
    shape_rows = [
        ("|e|", "e_abs"), ("sigma_e", "sigma_e"), ("rms_e", "rms_e"),
        ("HSC weight", "weight"), ("multiplicative m", "m"),
        ("mag_a10", "mag_a10"), ("resolution", "resolution"),
        ("PSF moment trace", "psf_trace"),
    ]
    distribution_table = [
        "| Retained quantity | Lenses q16 / median / q84 | Controls q16 / median / q84 |",
        "|---|---:|---:|",
    ]
    for label, key in shape_rows:
        distribution_table.append(
            f"| {label} | {fmt_distribution(lens_retained[key])} | {fmt_distribution(control_retained[key])} |"
        )
    selected_bins = summary["hsc_y3_zbin_retained_after_all_audit_cuts"]
    mask_counts = bmode["coordinate_crosscheck_counts"]
    support_policies = summary["zero_support_policy_sensitivity"]
    policy_lines = [
        "| Minimum positive controls | Eligible lens sets | Lens sets retained | Control/lens count ratio q16 / median / q84 |",
        "|---:|---:|---:|---:|",
    ]
    for minimum, policy in support_policies.items():
        policy_lines.append(
            f"| {minimum}/5 | {policy['n_eligible_lens_sets']} | "
            f"{100 * policy['fraction_of_80_lens_sets']:.1f}% | "
            f"{fmt_quantiles(policy['positive_control_count_to_lens_ratio_q16_q50_q84'])} |"
        )
    radial_support = summary["radial_tomographic_support"]
    radial_table_lines = [
        "| Role | R proper (Mpc) | Pair rows | z-bin pair rows (1/2/3/4) | Centers with support (1/2/3/4) |",
        "|---|---:|---:|---:|---:|",
    ]
    for role in ("lens", "control"):
        for radial_bin in range(len(RADIAL_EDGES_MPC) - 1):
            cells = [
                row for row in radial_support
                if row["role"] == role and row["radial_bin"] == radial_bin
            ]
            cells.sort(key=lambda row: row["source_zbin"])
            total = sum(int(row["n_rows"]) for row in cells)
            bins = "/".join(f"{int(row['n_rows']):,}" for row in cells)
            coverage = "/".join(str(row["n_centers_with_support"]) for row in cells)
            radial_table_lines.append(
                f"| {role} | {cells[0]['r_min_mpc']:.1f}–{cells[0]['r_max_mpc']:.1f} | "
                f"{total:,} | {bins} | {coverage} |"
            )
    lines = [
        "# S23B weak-lensing source-row quality and common-support audit",
        "",
        f"Approved row-plan SHA-256: `{summary['row_plan_sha256']}`  ",
        f"Input HSC Y3 n(z) SHA-256: `{summary['n_z_sacc_sha256']}`  ",
        f"Output: `{output}`  ",
        "No network request, additional download, or lensing stack was run.",
        "",
        "## Source cuts and semantics",
        "",
        f"The same source cuts were applied to target lenses and matched controls: `|e|<2`, `0<=sigma_e<=0.4`, `mag_a10<=25.5`, `resolution>=0.3`, HSC Y3 `source_zbin=1..4`, exclude the verified raw B-mode value 0, positive fiducial HSC Y3 `n(z)` lensing efficiency, and proper radius 0.1–3 Mpc. The official HSC S19A page defines `b_mode_mask` as a large-B-mode region cut but does not document the stored 0/1 direction. Prior coordinate validation of the same release identified raw 0 as masked; the new sample's coordinate comparison to the published GAMA09H rectangle found **{bmode['mismatch_rows']:,} mismatches** ({bmode['status']}).",
        f"The exact coordinate cross-check found {mask_counts['inside_mask0_exact']:,} raw-0 pairs inside the rectangle, {mask_counts['inside_mask1_exact']:,} raw-1 pairs inside, {mask_counts['outside_mask0_exact']:,} raw-0 pairs outside, and {mask_counts['outside_mask1_exact']:,} raw-1 pairs outside.",
        "",
        "Official definitions: [HSC S19A Shape Catalog (PDR3)](https://hsc-release.mtk.nao.ac.jp/doc/index.php/s19a-shape-catalog-pdr3/). Prior independent polarity audit: `wl_bmode_coordinate_validation_2026-10-19.md`.",
        "",
        "## Retained pair support",
        "",
        f"Lens apertures: {lens['n_input']:,} input pairs; {lens['n_shape_zbin_bmode']:,} pass the HSC shape, magnitude, resolution, tomographic-bin, and B-mode cuts; {lens['n_positive_nz_efficiency']:,} also have positive fiducial n(z) efficiency; {lens['n_calibration_valid']:,} remain after 0.1–3 Mpc and calibration-input-validity checks.",
        f"Five-control apertures: {control['n_input']:,} input pairs; {control['n_shape_zbin_bmode']:,} pass the same HSC cuts; {control['n_positive_nz_efficiency']:,} also have positive fiducial n(z) efficiency; {control['n_calibration_valid']:,} remain after the same radial and calibration-input-validity checks.",
        f"Retained HSC Y3 z-bin counts (1/2/3/4): lenses {selected_bins['lens']['1']:,}/{selected_bins['lens']['2']:,}/{selected_bins['lens']['3']:,}/{selected_bins['lens']['4']:,}; controls {selected_bins['control']['1']:,}/{selected_bins['control']['2']:,}/{selected_bins['control']['3']:,}/{selected_bins['control']['4']:,}.",
        f"At the paired-lens level, {matched['lenses_with_zero_retained_sources']} lenses have zero retained lens rows; {matched['lenses_with_at_least_one_zero_control']} have at least one zero-support control; {matched['lenses_with_all_five_zero_controls']} have all five controls at zero.",
        f"The q16/median/q84 ratios of mean-five-control to lens are {fmt_quantiles(matched['mean_control_to_lens_count_ratio_q16_q50_q84'])} by count, {fmt_quantiles(matched['mean_control_to_lens_shape_weight_ratio_q16_q50_q84'])} by HSC weight sum, and {fmt_quantiles(matched['mean_control_to_lens_sigma_crit_weight_ratio_q16_q50_q84'])} by geometry-weight proxy sum.",
        "The table below shows how many nonzero-lens sets remain under a minimum number of positive-support controls. Ratios use only positive-support controls for each set; the table is a support sensitivity, not a selected estimator.",
        "",
        *policy_lines,
        "",
        f"Center-level `N_eff` q16/median/q84 is {fmt_quantiles(lens['center_n_eff_shape_weight_q16_q50_q84'])} (lens) and {fmt_quantiles(control['center_n_eff_shape_weight_q16_q50_q84'])} (control) for catalog weight; geometry-weight proxy `N_eff` is {fmt_quantiles(lens['center_n_eff_sigma_crit_weight_proxy_q16_q50_q84'])} and {fmt_quantiles(control['center_n_eff_sigma_crit_weight_proxy_q16_q50_q84'])}.",
        "",
        "Retained shape and weight distributions (q16 / median / q84):",
        "",
        *distribution_table,
        "",
        "Radial/tomographic support is pair-counted, with the number of apertures having at least one retained source in each cell:",
        "",
        *radial_table_lines,
        "",
        "Per-aperture retained counts and effective-support measures are in `per_center.tsv`; paired five-control comparisons are in `per_lens_five_control_support.tsv`; 0.1–3 Mpc annular counts and weights are in `radial_support.tsv`. The 13,440-row `per_center_radial_tomographic_support.tsv` and its 56-row `radial_tomographic_support.tsv` aggregate the same support by radial annulus and HSC Y3 source bin.",
        "",
        "## Interpretation limits",
        "",
        "`N_eff_shape_weight` uses only the HSC catalog weight. The reported `shape_weight*(D_l beta)^2` is a geometry-weighted support proxy from fiducial HSC Y3 n(z); it omits per-object shear calibration, responsivity, HSC selection-bias correction, and n(z) uncertainty. These summaries do not constitute a shear profile or measurement. Repeated source IDs across apertures are retained pair rows and are not independent samples. PSF spatial matching, control contamination, source-reuse/field covariance, and injection/recovery remain unvalidated. `science_measurement=false`.",
        "",
    ]
    Path(path).write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan-dir", type=Path, default=DEFAULT_PLAN_DIR)
    parser.add_argument("--rows", type=Path, default=DEFAULT_ROWS)
    parser.add_argument("--nz-sacc", type=Path, default=DEFAULT_NZ)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    summary = run(args)
    print(json.dumps({
        "status": summary["status"], "science_measurement": summary["science_measurement"],
        "rows": summary["rows"], "cutflow": summary["shape_cutflow"],
        "matched_control_support": summary["matched_control_support"],
        "zero_support_policy_sensitivity": summary["zero_support_policy_sensitivity"],
        "radial_tomographic_cells": len(summary["radial_tomographic_support"]),
        "b_mode_polarity_validation": summary["b_mode_polarity_validation"],
        "output": str(args.output.resolve()),
    }, ensure_ascii=False, indent=2), flush=True)
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as error:
        message = str(error) if isinstance(error, ValueError) else type(error).__name__
        print(f"audit_s23b_wl_common_shape_rows: {message}", file=sys.stderr)
        sys.exit(1)
