#!/usr/bin/env python3
"""Propagate measured coadd lag covariance through the S23B ICL estimator.

The runner replays the saved S23B target/control masks and the production
bilinear metric-grid warp. It constructs the exact linear input-pixel weights
for the common-support 50--300, 100--300, and 200--300 minus 300--800 kpc
contrasts, then contracts those weights with each frame's measured short-lag
residual covariance kernel and delivered VARIANCE plane.

The covariance kernels were measured in outer-field patches. Applying them to
the cluster annuli assumes local stationarity, so the propagated values are a
conditional pixel-noise diagnostic, not a calibrated science covariance or an
error-bar correction. The analysis reads only already received files, makes
no network request, and keeps science_measurement=false.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import estimate_s23b_icl_coadd_noise_correlation as covariance
import stack_s23b_icl_target_control_2d as s23b_stack
from download_pilot import validate_fits
from validate_control_images import nJy_factor


DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
TARGET_MANIFEST = covariance.TARGET_PLAN / "target_manifest.json"
CONTROL_MANIFEST = covariance.CONTROL_PLAN / "control_manifest.json"
PAIR_PROFILES = covariance.STACK / "paired_radial_profiles.tsv"
STACK_SUMMARY = covariance.STACK / "summary.json"
COVARIANCE_DIR = DATA / "s23b_icl_residual_covariance_2026-11-07_01"
COVARIANCE_SUMMARY = COVARIANCE_DIR / "summary.json"
COVARIANCE_FRAMES = COVARIANCE_DIR / "frame_noise_correlation.tsv"
COVARIANCE_LAGS = COVARIANCE_DIR / "lag_correlation.tsv"
RELATION_DIR = DATA / "s23b_icl_relative_radial_contrasts_2026-11-07_02"
RELATION_SUMMARY = RELATION_DIR / "summary.json"
RELATION_ROWS = RELATION_DIR / "paired_relative_contrasts.tsv"
DEFAULT_OUTPUT = DATA / "s23b_icl_relative_covariance_propagation_2026-11-07_02"

EXPECTED_PAIR_PROFILE_SHA256 = "285906b03360e9d137727e948c15419767c356d3343e61f409562567ac4abc43"
EXPECTED_RELATION_ROWS_SHA256 = "c5de5f32b51c6899c8d17fe2b3e4b0a3966d810afd8db5ee27ba649fc0642800"
EXPECTED_CONTRASTS = (
    (50, 300, "50-300_minus_300-800"),
    (100, 300, "100-300_minus_300-800"),
    (200, 300, "200-300_minus_300-800"),
)
OUTER = (300, 800)
LAG_RADIUS = 4
N_SIDE = 401
KPC_MAX = 1000.0
MASK_FACTOR = 2.5
FILTERS = ("i", "i2")
SCHEMA = "astrafek-s23b-icl-relative-contrast-covariance-v1"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(8 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_json(path: Path) -> dict:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def read_tsv(path: Path) -> list[dict[str, str]]:
    with Path(path).open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f"refusing to write an empty table: {path}")
    with path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t",
                                lineterminator="\n", extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def write_json_atomic(path: Path, payload: dict) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True,
                                    allow_nan=False) + "\n", encoding="utf-8")
    temporary.replace(path)


def verify_inputs() -> tuple[list[dict], dict, dict]:
    records, record_meta = covariance.build_records()
    if len(records) != 2044 or record_meta["n_pairs"] != 1022:
        raise ValueError("the audited S23B covariance records do not contain 1,022 pairs")

    stack_summary = read_json(STACK_SUMMARY)
    if (stack_summary.get("status") != "complete_exploratory"
            or stack_summary.get("science_measurement") is not False
            or stack_summary.get("n_pairs_stacked") != 1022
            or float(stack_summary.get("mask_factor", -1)) != MASK_FACTOR
            or stack_summary.get("grid", {}).get("n_side") != N_SIDE
            or float(stack_summary.get("grid", {}).get("kpc_max", -1)) != KPC_MAX
            or stack_summary.get("outputs", {}).get("paired_profiles_sha256")
            != EXPECTED_PAIR_PROFILE_SHA256
            or sha256(PAIR_PROFILES) != EXPECTED_PAIR_PROFILE_SHA256):
        raise ValueError("the pinned S23B stack, grid, or mask definition changed")

    covariance_summary = read_json(COVARIANCE_SUMMARY)
    if (covariance_summary.get("status") != "complete_exploratory_residual_covariance_diagnostic"
            or covariance_summary.get("science_measurement") is not False
            or covariance_summary.get("sample", {}).get("n_pairs_in_input_stack") != 1022
            or covariance_summary.get("method", {}).get("maximum_lag_pixels") != LAG_RADIUS
            or sha256(COVARIANCE_FRAMES)
            != covariance_summary.get("outputs", {}).get("frame_noise_correlation.tsv")
            or sha256(COVARIANCE_LAGS)
            != covariance_summary.get("outputs", {}).get("lag_correlation.tsv")):
        raise ValueError("the pinned S23B lag-covariance diagnostic is incomplete or changed")

    relation_summary = read_json(RELATION_SUMMARY)
    if (relation_summary.get("status") != "complete_exploratory_relative_contrast"
            or relation_summary.get("science_measurement") is not False
            or relation_summary.get("n_pairs_by_filter") != {"i": 241, "i2": 781}
            or sha256(RELATION_ROWS) != EXPECTED_RELATION_ROWS_SHA256
            or relation_summary.get("outputs", {}).get("paired_relative_contrasts.tsv", {}).get("sha256")
            != EXPECTED_RELATION_ROWS_SHA256):
        raise ValueError("the pinned S23B relative-contrast table is incomplete or changed")

    hashes = {
        str(path.resolve()): sha256(path)
        for path in (TARGET_MANIFEST, CONTROL_MANIFEST, STACK_SUMMARY, PAIR_PROFILES,
                     COVARIANCE_SUMMARY, COVARIANCE_FRAMES, COVARIANCE_LAGS,
                     RELATION_SUMMARY, RELATION_ROWS)
    }
    hashes.update(record_meta["provenance"])
    return records, stack_summary, {
        "input_sha256": hashes,
        "record_meta": record_meta,
        "covariance_summary": covariance_summary,
        "relation_summary": relation_summary,
    }


def read_frame_kernels() -> tuple[dict, dict]:
    frame_rows = read_tsv(COVARIANCE_FRAMES)
    lag_rows = read_tsv(COVARIANCE_LAGS)
    frame_key = lambda row: (row["pair_id"], row["cluster_id"], row["role"], row["filter"])
    frame_meta = {frame_key(row): row for row in frame_rows}
    if len(frame_meta) != len(frame_rows) or len(frame_rows) != 2044:
        raise ValueError("the per-frame covariance table has duplicate or missing records")

    lag_values: dict[tuple, dict[tuple[int, int], float]] = {}
    expected_lags = {(dy, dx) for dy in range(-LAG_RADIUS, LAG_RADIUS + 1)
                     for dx in range(-LAG_RADIUS, LAG_RADIUS + 1)}
    for row in lag_rows:
        key = frame_key(row)
        lag = (int(row["dy"]), int(row["dx"]))
        if key not in frame_meta or lag in lag_values.setdefault(key, {}):
            raise ValueError(f"duplicate or unknown lag row: {key}/{lag}")
        value = float(row["median_variance_normalized_covariance"])
        if not math.isfinite(value):
            raise ValueError(f"non-finite measured lag covariance: {key}/{lag}")
        lag_values[key][lag] = value

    complete: dict[tuple, np.ndarray] = {}
    for key, values in lag_values.items():
        if set(values) != expected_lags:
            raise ValueError(f"incomplete lag grid for {key}")
        kernel = np.empty((2 * LAG_RADIUS + 1, 2 * LAG_RADIUS + 1), dtype=np.float64)
        for dy, dx in expected_lags:
            kernel[dy + LAG_RADIUS, dx + LAG_RADIUS] = 0.5 * (
                values[(dy, dx)] + values[(-dy, -dx)])
        complete[key] = kernel

    missing = set(frame_meta) - set(complete)
    if len(missing) != 1:
        raise ValueError(f"expected exactly one known insufficient-patch frame, found {len(missing)}")
    fallback_kernels = {}
    for band in FILTERS:
        for role in ("target", "matched_control"):
            peers = [kernel for key, kernel in complete.items()
                     if key[2] == role and key[3] == band]
            if not peers:
                raise ValueError(f"no measured covariance kernels for {band}/{role}")
            fallback_kernels[(band, role)] = np.median(np.stack(peers), axis=0)
    for key in missing:
        complete[key] = fallback_kernels[(key[3], key[2])]
    return complete, {"missing_kernel_keys": [list(key) for key in sorted(missing)],
                      "frame_count": len(frame_meta), "measured_kernel_count": len(lag_values),
                      "imputation": "per-filter and per-role median of measured frame kernels"}


def load_mask(path: Path, expected_sha: str) -> np.ndarray:
    if sha256(path) != expected_sha:
        raise ValueError(f"prepared mask checksum mismatch: {path}")
    with fits.open(path, memmap=True, checksum=True) as hdus:
        if "MASK_K2P5" not in hdus:
            raise ValueError(f"prepared mask lacks MASK_K2P5: {path}")
        return np.asarray(hdus["MASK_K2P5"].data, dtype=bool)


def broad_values(target_map: np.ndarray, control_map: np.ndarray,
                 common: np.ndarray, radius: np.ndarray,
                 pair_profile: list[dict], lo: int, hi: int) -> tuple[int, float | None]:
    selected = (radius >= lo) & (radius < hi) & common
    count = int(selected.sum())
    observed = (float(np.mean(target_map[selected] - control_map[selected], dtype=np.float64))
                if count else None)
    expected_count = sum(int(row["n_valid_pixels"]) for row in pair_profile
                         if float(row["r_min_kpc"]) >= lo
                         and float(row["r_max_kpc"]) <= hi)
    expected_sum = sum(float(row["delta_mean_njy_arcsec2"]) * int(row["n_valid_pixels"])
                       for row in pair_profile
                       if float(row["r_min_kpc"]) >= lo
                       and float(row["r_max_kpc"]) <= hi
                       and int(row["n_valid_pixels"]) > 0)
    if count != expected_count:
        raise ValueError(f"replayed {lo}-{hi} kpc common-pixel count {count} != {expected_count}")
    if count == 0:
        if expected_sum != 0.0:
            raise ValueError(f"empty {lo}-{hi} kpc support has a non-zero saved profile sum")
        return count, None
    expected_mean = expected_sum / expected_count
    if observed is None or not math.isclose(observed, expected_mean, rel_tol=0, abs_tol=2e-5):
        raise ValueError(f"replayed {lo}-{hi} kpc profile mean differs from saved stack")
    return count, observed


def output_estimator(common: np.ndarray, radius: np.ndarray,
                     inner: tuple[int, int], outer: tuple[int, int]) -> np.ndarray:
    inner_mask = common & (radius >= inner[0]) & (radius < inner[1])
    outer_mask = common & (radius >= outer[0]) & (radius < outer[1])
    n_inner, n_outer = int(inner_mask.sum()), int(outer_mask.sum())
    if n_inner < 1 or n_outer < 1:
        raise ValueError(f"missing support in {inner} / {outer} kpc estimator")
    coefficient = np.zeros(common.shape, dtype=np.float64)
    coefficient[inner_mask] = 1.0 / n_inner
    coefficient[outer_mask] = -1.0 / n_outer
    return coefficient


def scatter_bilinear_coefficients(output_coefficient: np.ndarray, *, axis_kpc: np.ndarray,
                                  center_xy: tuple[float, float], kpc_per_pixel: float,
                                  native_to_sb: float, native_shape: tuple[int, int]) -> tuple[np.ndarray, tuple[int, int]]:
    """Apply the transpose of warp_metric's order-1 interpolation to an aperture."""
    output_y, output_x = np.nonzero(output_coefficient)
    if output_y.size == 0:
        raise ValueError("empty output-grid estimator")
    cx, cy = center_xy
    source_x = cx + axis_kpc[output_x] / kpc_per_pixel
    source_y = cy + axis_kpc[output_y] / kpc_per_pixel
    x0 = np.floor(source_x).astype(np.int64)
    y0 = np.floor(source_y).astype(np.int64)
    fx, fy = source_x - x0, source_y - y0
    out_weight = output_coefficient[output_y, output_x]
    height, width = native_shape
    x_min = max(0, int(x0.min()))
    y_min = max(0, int(y0.min()))
    x_max = min(width, int(x0.max()) + 2)
    y_max = min(height, int(y0.max()) + 2)
    if x_max <= x_min or y_max <= y_min:
        raise ValueError("bilinear source footprint is outside the input image")
    coefficient = np.zeros((y_max - y_min, x_max - x_min), dtype=np.float64)
    for dx, dy, weight in (
        (0, 0, (1.0 - fx) * (1.0 - fy)),
        (1, 0, fx * (1.0 - fy)),
        (0, 1, (1.0 - fx) * fy),
        (1, 1, fx * fy),
    ):
        ix, iy = x0 + dx, y0 + dy
        keep = weight > 0.0
        if np.any(keep & ((ix < 0) | (ix >= width) | (iy < 0) | (iy >= height))):
            raise ValueError("valid output sample used an out-of-bounds input pixel")
        if np.any(keep):
            np.add.at(coefficient, (iy[keep] - y_min, ix[keep] - x_min),
                      out_weight[keep] * weight[keep] * native_to_sb)
    if not np.any(coefficient):
        raise ValueError("bilinear transpose produced no input weights")
    return coefficient, (y_min, x_min)


def contract_kernel(coefficient: np.ndarray, variance: np.ndarray,
                    kernel: np.ndarray) -> tuple[float, float]:
    """Return diagonal-plane and full finite-kernel estimator variances."""
    y, x = np.nonzero(coefficient)
    weights = coefficient[y, x]
    variances = variance[y, x].astype(np.float64, copy=False)
    if (not np.isfinite(weights).all() or not np.isfinite(variances).all()
            or np.any(variances <= 0)):
        raise ValueError("estimator uses a non-finite or non-positive VARIANCE pixel")
    diagonal = float(np.dot(weights * weights, variances))
    full = 0.0
    height, width = coefficient.shape
    for dy in range(-LAG_RADIUS, LAG_RADIUS + 1):
        iy = y + dy
        valid_y = (iy >= 0) & (iy < height)
        for dx in range(-LAG_RADIUS, LAG_RADIUS + 1):
            rho = float(kernel[dy + LAG_RADIUS, dx + LAG_RADIUS])
            if rho == 0.0:
                continue
            ix = x + dx
            keep = valid_y & (ix >= 0) & (ix < width)
            if not np.any(keep):
                continue
            next_weights = coefficient[iy[keep], ix[keep]]
            nonzero = next_weights != 0.0
            if not np.any(nonzero):
                continue
            y_keep, x_keep = iy[keep][nonzero], ix[keep][nonzero]
            left_var = variances[keep][nonzero]
            right_var = variance[y_keep, x_keep].astype(np.float64, copy=False)
            if not np.isfinite(right_var).all() or np.any(right_var <= 0):
                raise ValueError("lag pair touches a non-finite or non-positive VARIANCE pixel")
            full += rho * float(np.sum(
                weights[keep][nonzero] * next_weights[nonzero]
                * np.sqrt(left_var * right_var), dtype=np.float64))
    if not math.isfinite(diagonal) or diagonal <= 0:
        raise ValueError("non-positive diagonal propagated variance")
    return diagonal, full


def read_frame(record: dict, frame_meta: dict, role_meta: dict,
               z: float, cosmology) -> tuple[dict, np.ndarray, np.ndarray, np.ndarray]:
    image_path, mask_path = Path(record["image_path"]), Path(record["mask_path"])
    if image_path.stat().st_size != int(record["image_bytes"]):
        raise ValueError(f"image byte count changed: {image_path}")
    image_sha = sha256(image_path)
    if image_sha != record["image_sha256"]:
        raise ValueError(f"image SHA-256 changed: {image_path}")
    info = validate_fits(image_path, frame_meta)
    if info["coadd_header_filter"] != record["filter"]:
        raise ValueError(f"filter mismatch: {image_path}")
    if not np.allclose(info["center_pixel_xy"], [record["center_x"], record["center_y"]],
                       rtol=0, atol=1e-6):
        raise ValueError(f"center pixel differs from audited frame record: {image_path}")
    with fits.open(image_path, memmap=True, checksum=True) as hdus:
        image = np.asarray(hdus[info["science_hdu"]].data, dtype=np.float32)
        variance = np.asarray(hdus[info["variance_hdu"]].data, dtype=np.float32)
    source_mask = load_mask(mask_path, record["mask_sha256"])
    if image.shape != variance.shape or image.shape != source_mask.shape:
        raise ValueError(f"image, variance, and saved mask grids differ: {image_path}")
    if (str(frame_meta.get("id")) != str(record["cluster_id"])
            or float(frame_meta["z"]) != z
            or float(role_meta["z"]) != z):
        raise ValueError(f"unexpected frame metadata: {image_path}")
    if record["role"] == "target":
        if frame_meta.get("sample_role") != "target":
            raise ValueError(f"target frame metadata has the wrong sample role: {image_path}")
    elif (role_meta.get("frame_role") != "matched_control"
          or role_meta.get("science_measurement") is not False
          or str(role_meta.get("target_cluster_id")) != str(record["pair_id"])):
        raise ValueError(f"control frame metadata is outside its matched pair: {image_path}")
    warped = s23b_stack.stack.target_and_control_map(
        image, source_mask, info, z, cosmology, KPC_MAX, N_SIDE)
    return info, image, variance, warped


def process_pair(pair: dict[str, str], records_by_role: dict[str, dict],
                 target_meta: dict, control_meta: dict, profile_rows: list[dict],
                 relation_lookup: dict, kernels: dict, imputed_kernel_keys: set,
                 cosmology, radius: np.ndarray) -> list[dict]:
    target_record = records_by_role["target"]
    control_record = records_by_role["matched_control"]
    target_info, target_image, target_variance, target_warp = read_frame(
        target_record, target_meta, target_meta, float(pair["z"]), cosmology)
    control_info, control_image, control_variance, control_warp = read_frame(
        control_record, control_meta, control_meta, float(pair["z"]), cosmology)
    common = ((target_warp["weight"] > 0) & (control_warp["weight"] > 0)
              & np.isfinite(target_warp["sb"]) & np.isfinite(control_warp["sb"]))
    if not np.any(common):
        raise ValueError(f"empty common warped support: {pair['target_id']}")
    target_cov_key = (target_record["pair_id"], target_record["cluster_id"],
                      target_record["role"], target_record["filter"])
    control_cov_key = (control_record["pair_id"], control_record["cluster_id"],
                       control_record["role"], control_record["filter"])
    target_kernel, control_kernel = kernels[target_cov_key], kernels[control_cov_key]
    profile = json.loads(pair["profile_json"])
    output_rows = []

    for inner_low, inner_high, contrast_name in EXPECTED_CONTRASTS:
        n_inner, inner_mean = broad_values(
            target_warp["sb"], control_warp["sb"], common, radius, profile,
            inner_low, inner_high)
        n_outer, outer_mean = broad_values(
            target_warp["sb"], control_warp["sb"], common, radius, profile,
            *OUTER)
        saved = relation_lookup[(pair["target_id"], pair["filter"], contrast_name)]
        saved_inner = (float(saved["inner_mean_njy_arcsec2"])
                       if saved["inner_mean_njy_arcsec2"] not in (None, "") else None)
        saved_outer = (float(saved["outer_mean_njy_arcsec2"])
                       if saved["outer_mean_njy_arcsec2"] not in (None, "") else None)
        saved_contrast = (float(saved["contrast_njy_arcsec2"])
                          if saved["contrast_njy_arcsec2"] not in (None, "") else None)
        if (n_inner != int(saved["inner_common_pixels"])
                or n_outer != int(saved["outer_common_pixels"])
                or ((inner_mean is None) != (saved_inner is None))
                or ((outer_mean is None) != (saved_outer is None))
                or ((inner_mean is None or saved_inner is None)
                    and n_inner != 0)
                or ((outer_mean is None or saved_outer is None)
                    and n_outer != 0)
                or (inner_mean is not None and saved_inner is not None
                    and not math.isclose(inner_mean, saved_inner, rel_tol=0, abs_tol=2e-5))
                or (outer_mean is not None and saved_outer is not None
                    and not math.isclose(outer_mean, saved_outer, rel_tol=0, abs_tol=2e-5))
                or ((n_inner == 0 or n_outer == 0) != (saved_contrast is None))
                or (saved_contrast is not None and inner_mean is not None and outer_mean is not None
                    and not math.isclose(inner_mean - outer_mean, saved_contrast,
                                         rel_tol=0, abs_tol=3e-5))):
            raise ValueError(f"exact estimator replay disagrees with saved contrast: {pair['target_id']}/{contrast_name}")

        if n_inner == 0 or n_outer == 0:
            output_rows.append({
                "target_id": pair["target_id"], "control_id": pair["control_id"],
                "filter": pair["filter"], "z": float(pair["z"]), "contrast": contrast_name,
                "inner_r_min_kpc": inner_low, "inner_r_max_kpc": inner_high,
                "outer_r_min_kpc": OUTER[0], "outer_r_max_kpc": OUTER[1],
                "n_inner_common_pixels": n_inner, "n_outer_common_pixels": n_outer,
                "observed_contrast_njy_arcsec2": saved_contrast,
                "target_diagonal_variance": None, "control_diagonal_variance": None,
                "target_zero_lag_variance": None, "control_zero_lag_variance": None,
                "target_finite_lag_variance": None, "control_finite_lag_variance": None,
                "pair_diagonal_variance": None, "pair_zero_lag_variance": None,
                "pair_finite_lag_variance": None, "pair_finite_lag_sigma_njy_arcsec2": None,
                "target_kernel_zero_lag": float(target_kernel[LAG_RADIUS, LAG_RADIUS]),
                "control_kernel_zero_lag": float(control_kernel[LAG_RADIUS, LAG_RADIUS]),
                "target_kernel_imputed": target_cov_key in imputed_kernel_keys,
                "control_kernel_imputed": control_cov_key in imputed_kernel_keys,
                "status": "insufficient_common_support",
                "science_measurement": False,
            })
            continue

        output_weights = output_estimator(common, radius,
                                          (inner_low, inner_high), OUTER)
        target_coeff, target_origin = scatter_bilinear_coefficients(
            output_weights, axis_kpc=target_warp["kpc"],
            center_xy=tuple(target_info["center_pixel_xy"]),
            kpc_per_pixel=float(target_warp["kpc_per_pixel"]),
            native_to_sb=nJy_factor(target_info), native_shape=target_image.shape)
        control_coeff, control_origin = scatter_bilinear_coefficients(
            output_weights, axis_kpc=control_warp["kpc"],
            center_xy=tuple(control_info["center_pixel_xy"]),
            kpc_per_pixel=float(control_warp["kpc_per_pixel"]),
            native_to_sb=nJy_factor(control_info), native_shape=control_image.shape)
        target_y, target_x = target_origin
        control_y, control_x = control_origin
        target_var_crop = target_variance[target_y:target_y + target_coeff.shape[0],
                                          target_x:target_x + target_coeff.shape[1]]
        control_var_crop = control_variance[control_y:control_y + control_coeff.shape[0],
                                            control_x:control_x + control_coeff.shape[1]]
        if target_var_crop.shape != target_coeff.shape or control_var_crop.shape != control_coeff.shape:
            raise ValueError("coefficient and cropped VARIANCE shapes differ")
        target_diag, target_full = contract_kernel(target_coeff, target_var_crop, target_kernel)
        control_diag, control_full = contract_kernel(control_coeff, control_var_crop, control_kernel)
        target_zero = target_diag * float(target_kernel[LAG_RADIUS, LAG_RADIUS])
        control_zero = control_diag * float(control_kernel[LAG_RADIUS, LAG_RADIUS])
        full_total = target_full + control_full
        if full_total <= 0 or not math.isfinite(full_total):
            raise ValueError(f"non-positive full covariance variance: {pair['target_id']}/{contrast_name}")
        output_rows.append({
            "target_id": pair["target_id"], "control_id": pair["control_id"],
            "filter": pair["filter"], "z": float(pair["z"]), "contrast": contrast_name,
            "inner_r_min_kpc": inner_low, "inner_r_max_kpc": inner_high,
            "outer_r_min_kpc": OUTER[0], "outer_r_max_kpc": OUTER[1],
            "n_inner_common_pixels": n_inner, "n_outer_common_pixels": n_outer,
            "observed_contrast_njy_arcsec2": float(saved["contrast_njy_arcsec2"]),
            "target_diagonal_variance": target_diag,
            "control_diagonal_variance": control_diag,
            "target_zero_lag_variance": target_zero,
            "control_zero_lag_variance": control_zero,
            "target_finite_lag_variance": target_full,
            "control_finite_lag_variance": control_full,
            "pair_diagonal_variance": target_diag + control_diag,
            "pair_zero_lag_variance": target_zero + control_zero,
            "pair_finite_lag_variance": full_total,
            "pair_finite_lag_sigma_njy_arcsec2": math.sqrt(full_total),
            "target_kernel_zero_lag": float(target_kernel[LAG_RADIUS, LAG_RADIUS]),
            "control_kernel_zero_lag": float(control_kernel[LAG_RADIUS, LAG_RADIUS]),
            "target_kernel_imputed": target_cov_key in imputed_kernel_keys,
            "control_kernel_imputed": control_cov_key in imputed_kernel_keys,
            "status": "complete_conditional_covariance",
            "science_measurement": False,
        })
    return output_rows


def summarize(rows: list[dict], relation_summary: dict) -> list[dict]:
    summaries = []
    block_rows = read_tsv(RELATION_DIR / "contrast_block_uncertainty.tsv")
    block_lookup = {(row["filter"], row["contrast"], float(row["block_width_deg"]),
                     float(row["grid_phase_fraction"])): row for row in block_rows}
    for band in FILTERS:
        for _, _, contrast in EXPECTED_CONTRASTS:
            sample = [row for row in rows
                      if row["filter"] == band and row["contrast"] == contrast
                      and row["status"] == "complete_conditional_covariance"]
            if not sample:
                raise ValueError(f"no propagated pair results for {band}/{contrast}")
            n = len(sample)
            observed = np.asarray([float(row["observed_contrast_njy_arcsec2"])
                                   for row in sample])
            variances = {
                "delivered_variance_diagonal": np.asarray([float(row["pair_diagonal_variance"])
                                                            for row in sample]),
                "measured_zero_lag_only": np.asarray([float(row["pair_zero_lag_variance"])
                                                       for row in sample]),
                "measured_finite_lag_kernel": np.asarray([float(row["pair_finite_lag_variance"])
                                                           for row in sample]),
            }
            reference = block_lookup[(band, contrast, 5.0, 0.0)]
            expected_supported = int(reference["n_supported_pairs"])
            if n != expected_supported:
                raise ValueError(f"propagated support for {band}/{contrast} is {n}, expected {expected_supported}")
            sigma_pixel = {name: math.sqrt(float(values.sum()) / (n * n))
                           for name, values in variances.items()}
            pair_sem = float(observed.std(ddof=1) / math.sqrt(n)) if n > 1 else None
            ci_low = (float(reference["block_bootstrap_ci95_low"])
                      if reference.get("block_bootstrap_ci95_low") not in (None, "") else None)
            ci_high = (float(reference["block_bootstrap_ci95_high"])
                       if reference.get("block_bootstrap_ci95_high") not in (None, "") else None)
            half_width = ((ci_high - ci_low) / 2.0
                          if ci_low is not None and ci_high is not None else None)
            summaries.append({
                "filter": band, "contrast": contrast, "n_pairs_propagated": n,
                "n_pairs_in_saved_stack": int(relation_summary["n_pairs_by_filter"][band]),
                "n_pairs_in_contrast": expected_supported,
                "mean_observed_contrast_njy_arcsec2": float(observed.mean()),
                "pair_sem_njy_arcsec2": pair_sem,
                "block_5deg_phase0_ci95_low": ci_low,
                "block_5deg_phase0_ci95_high": ci_high,
                "independent_pixel_noise_sigma_of_equal_pair_mean": sigma_pixel["delivered_variance_diagonal"],
                "measured_zero_lag_noise_sigma_of_equal_pair_mean": sigma_pixel["measured_zero_lag_only"],
                "measured_finite_lag_noise_sigma_of_equal_pair_mean": sigma_pixel["measured_finite_lag_kernel"],
                "finite_lag_to_diagonal_sigma_ratio": (
                    sigma_pixel["measured_finite_lag_kernel"]
                    / sigma_pixel["delivered_variance_diagonal"]),
                "finite_lag_noise_over_pair_sem": (
                    sigma_pixel["measured_finite_lag_kernel"] / pair_sem
                    if pair_sem and pair_sem > 0 else None),
                "finite_lag_noise_over_block_ci_half_width": (
                    sigma_pixel["measured_finite_lag_kernel"] / half_width
                    if half_width and half_width > 0 else None),
                "science_measurement": False,
            })
    return summaries


def render_plot(path: Path, summaries: list[dict]) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 5.1), sharey=True)
    colors = {"delivered_variance_diagonal": "#4976a8",
              "measured_zero_lag_noise": "#d08b2c",
              "measured_finite_lag_noise": "#9f4f51"}
    keys = (("independent_pixel_noise_sigma_of_equal_pair_mean", "delivered_variance_diagonal",
             "VARIANCE diagonal"),
            ("measured_zero_lag_noise_sigma_of_equal_pair_mean", "measured_zero_lag_noise",
             "Measured zero lag"),
            ("measured_finite_lag_noise_sigma_of_equal_pair_mean", "measured_finite_lag_noise",
             "Measured ±4 lag"))
    labels = ["50–300", "100–300", "200–300"]
    x = np.arange(3, dtype=float)
    width = 0.22
    for ax, band in zip(axes, FILTERS):
        selected = [row for row in summaries if row["filter"] == band]
        selected.sort(key=lambda row: next(i for i, item in enumerate(EXPECTED_CONTRASTS)
                                           if item[2] == row["contrast"]))
        for j, (key, color_key, label) in enumerate(keys):
            values = [float(row[key]) for row in selected]
            ax.bar(x + (j - 1) * width, values, width=width,
                   color=colors[color_key], label=label)
        ax.plot(x, [float(row["block_5deg_phase0_ci95_high"] - row["block_5deg_phase0_ci95_low"])
                    / 2.0 for row in selected], color="#222222", marker="o",
                linestyle="none", label="Block-CI half-width")
        ax.set_xticks(x, labels)
        ax.set_title(f"HSC-{band}", loc="left")
        ax.set_xlabel("Inner annulus [proper kpc]")
        ax.grid(axis="y", alpha=0.2)
        ax.legend(loc="upper right", fontsize=8, frameon=True)
    axes[0].set_ylabel("Noise σ or block-CI half-width [nJy arcsec⁻²]")
    fig.tight_layout()
    fig.savefig(path, dpi=170)
    plt.close(fig)


def run(output: Path, *, resume: bool) -> dict:
    output = output.resolve()
    records, stack_summary, provenance = verify_inputs()
    kernels, kernel_audit = read_frame_kernels()
    imputed_kernel_keys = {tuple(key) for key in kernel_audit["missing_kernel_keys"]}
    target_doc, control_doc = read_json(TARGET_MANIFEST), read_json(CONTROL_MANIFEST)
    targets = {str(row["id"]): row for row in target_doc["clusters"]}
    controls = {str(row["target_cluster_id"]): row for row in control_doc["clusters"]}
    if len(targets) != 1043 or len(controls) != 1041:
        raise ValueError("S23B target/control manifests changed")
    profiles = read_tsv(PAIR_PROFILES)
    relations = read_tsv(RELATION_ROWS)
    relation_lookup = {(row["target_id"], row["filter"], row["contrast"]): row
                       for row in relations}
    expected_relation_rows = 1022 * len(EXPECTED_CONTRASTS)
    if len(profiles) != 1022 or len(relation_lookup) != expected_relation_rows:
        raise ValueError("saved pair/profile contrast tables do not have the expected rows")
    by_pair_role: dict[str, dict[str, dict]] = {}
    for record in records:
        role = str(record["role"])
        if role not in ("target", "matched_control"):
            raise ValueError(f"unknown frame role in covariance records: {role}")
        by_pair_role.setdefault(str(record["pair_id"]), {})[role] = record
    if len(by_pair_role) != 1022 or any(set(pair) != {"target", "matched_control"}
                                         for pair in by_pair_role.values()):
        raise ValueError("audited records do not contain one target and control per pair")

    manifest = {
        "schema": SCHEMA, "output_grid": {"kpc_max": KPC_MAX, "n_side": N_SIDE},
        "mask_factor": MASK_FACTOR, "lag_radius_pixels": LAG_RADIUS,
        "contrasts": [list(item) for item in EXPECTED_CONTRASTS], "outer_kpc": list(OUTER),
        "input_sha256": provenance["input_sha256"],
        "runner_sha256": sha256(Path(__file__).resolve()),
        "kernel_assumption": "per-frame outer-field normalized kernel is stationary when applied to the full relative annulus estimator",
    }
    manifest_path, progress_path = output / "run_manifest.json", output / "progress.json"
    jsonl_path = output / "pair_measurements.jsonl"
    if output.exists():
        if not resume:
            raise FileExistsError(f"refusing to overwrite {output}; pass --resume for an incomplete run")
        if not manifest_path.is_file() or read_json(manifest_path) != manifest:
            raise ValueError("resume inputs or runner differ from the saved run manifest")
        if (output / "summary.json").is_file():
            raise ValueError("completed output cannot be resumed")
    else:
        if resume:
            raise FileNotFoundError("--resume requested but output does not exist")
        output.mkdir(parents=True, exist_ok=False)
        write_json_atomic(manifest_path, manifest)

    completed: dict[str, list[dict]] = {}
    if jsonl_path.is_file():
        with jsonl_path.open(encoding="utf-8") as stream:
            for line in stream:
                if line.strip():
                    row = json.loads(line)
                    if row["target_id"] in completed:
                        raise ValueError(f"duplicate completed pair checkpoint: {row['target_id']}")
                    completed[row["target_id"]] = row["measurements"]
    if not resume and completed:
        raise ValueError("output unexpectedly contains prior pair checkpoints")

    axis = np.linspace(-KPC_MAX, KPC_MAX, N_SIDE, dtype=np.float64)
    yy, xx = np.meshgrid(axis, axis, indexing="ij")
    radius = np.hypot(xx, yy)
    cosmology = s23b_stack.stack.Cosmology(
        H0=float(target_doc["h0"]), Om=float(target_doc["omega_m"]))
    target_by_id = targets
    control_by_target = controls
    append_mode = "a" if jsonl_path.exists() else "x"
    done = len(completed)
    with jsonl_path.open(append_mode, encoding="utf-8") as checkpoint:
        for index, pair in enumerate(profiles, 1):
            target_id = str(pair["target_id"])
            if target_id in completed:
                continue
            if target_id not in target_by_id or target_id not in control_by_target:
                raise ValueError(f"pair lacks manifest records: {target_id}")
            role_records = by_pair_role[target_id]
            if (role_records["target"]["filter"] != pair["filter"]
                    or role_records["matched_control"]["filter"] != pair["filter"]
                    or role_records["matched_control"]["cluster_id"] != pair["control_id"]):
                raise ValueError(f"pair image records differ from saved stack: {target_id}")
            measurements = process_pair(
                pair, role_records, target_by_id[target_id],
                control_by_target[target_id], json.loads(pair["profile_json"]),
                relation_lookup, kernels, imputed_kernel_keys, cosmology, radius)
            checkpoint.write(json.dumps({"target_id": target_id,
                                         "measurements": measurements}, allow_nan=False) + "\n")
            checkpoint.flush()
            completed[target_id] = measurements
            done += 1
            if done % 10 == 0 or done == len(profiles):
                write_json_atomic(progress_path, {
                    "status": "running", "completed_pairs": done,
                    "total_pairs": len(profiles), "current_pair_index": index,
                    "science_measurement": False,
                })
                print(f"propagated covariance for {done}/{len(profiles)} pairs", flush=True)

    rows = [measurement for pair_rows in completed.values() for measurement in pair_rows]
    if len(rows) != len(profiles) * len(EXPECTED_CONTRASTS):
        raise ValueError("completed pair checkpoint count does not close")
    summary_rows = summarize(rows, provenance["relation_summary"])
    pair_table = output / "pair_measurements.tsv"
    summary_table = output / "filter_summary.tsv"
    ordered_rows = sorted(rows, key=lambda row: (row["filter"], row["target_id"], row["contrast"]))
    write_tsv(pair_table, ordered_rows)
    write_tsv(summary_table, summary_rows)
    figure_path = output / "relative_contrast_noise_budget.png"
    render_plot(figure_path, summary_rows)

    result = {
        "schema": SCHEMA,
        "status": "complete_conditional_covariance_propagation",
        "science_measurement": False,
        "images_ready_for_science": False,
        "network_used": False,
        "new_data_downloaded": False,
        "science_pixels_read": True,
        "sample": {"n_pairs": len(profiles), "by_filter": {band: sum(
            row["filter"] == band for row in profiles) for band in FILTERS}},
        "estimator": {
            "operation": "saved stack's exact bilinear warp, per-pair common-valid mask, pixel-count weighted broad-annulus means, and inner-minus-outer contrast",
            "grid": {"kpc_max": KPC_MAX, "n_side": N_SIDE,
                     "step_kpc": float(axis[1] - axis[0])},
            "mask_factor": MASK_FACTOR,
            "contrasts": [list(item) for item in EXPECTED_CONTRASTS],
            "outer_kpc": list(OUTER),
        },
        "covariance_model": {
            "kernel_source": str(COVARIANCE_LAGS),
            "measured_lag_radius_pixels": LAG_RADIUS,
            "kernel_application": "per-frame median variance-normalized residual kernel from outer-field patches; assumed stationary over each annulus",
            "variance_plane": "native HSC coadd VARIANCE; converted through the exact bilinear transpose weights",
            "target_control_cross_covariance": "assumed zero because target and random-control images are separate fields",
            "missing_kernel_audit": kernel_audit,
            "population_noise_floor": "independent pair noise; sum per-frame target+control propagated variances and divide by pair count squared",
        },
        "filter_summary": summary_rows,
        "input_sha256": provenance["input_sha256"],
        "analysis_script_sha256": sha256(Path(__file__).resolve()),
        "outputs": {
            "pair_measurements.tsv": sha256(pair_table),
            "filter_summary.tsv": sha256(summary_table),
            "relative_contrast_noise_budget.png": sha256(figure_path),
        },
        "limitations": [
            "The covariance kernel is estimated outside a 60-arcsec center exclusion, then transferred to cluster annuli; spatial stationarity is unverified.",
            "Only short lags through four native coadd pixels are propagated; long-range sky covariance is not included.",
            "The delivered VARIANCE plane does not include source Poisson, sky-subtraction uncertainty, PSF wings, extinction, centering, or mask-selection uncertainty.",
            "Pair pixel-noise terms are treated as independent across clusters; overlapping coadd footprints and large-scale covariance are not included.",
            "The estimates describe the variance of the saved exploratory annular contrast under this model and are not applied as profile corrections or final error bars.",
            "The existing connected sky-block intervals remain the population uncertainty reference; the pixel-noise term is a conditional diagnostic component.",
        ],
        "science_measurement_gate": False,
    }
    write_json_atomic(output / "summary.json", result)
    progress_path.unlink(missing_ok=True)
    print(json.dumps({"status": result["status"], "output": str(output),
                      "filter_summary": summary_rows,
                      "summary_sha256": sha256(output / "summary.json")}, indent=2), flush=True)
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    run(args.output, resume=args.resume)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
