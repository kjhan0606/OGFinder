#!/usr/bin/env python3
"""Estimate ICL/satellite covariance from shared estimator pixels.

The satellite estimator includes its accepted azimuthal control annuli with
negative weights. For each matched cluster/filter, this script reconstructs
the actual saved target and control supports, combines equal-satellite weights,
and computes the covariance implied by the HSC variance map. An optional
measured short-range lag kernel estimates local coadd/residual covariance; this
remains an exploratory noise diagnostic, not a total science error.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from astropy.cosmology import FlatLambdaCDM
from astropy.io import fits
from astropy import units as u
from scipy.spatial import cKDTree
from scipy.signal import fftconvolve

import measure_satellite_outer_profiles as satellite
from audit_satellite_icl_pixel_overlap import (
    DEFAULT_ICL, DEFAULT_MASKS, DEFAULT_PREP, DEFAULT_PROFILES, ROOT,
    ICL_INNER, ICL_OUTER, SAT_INNER, SAT_OUTER, SAT_COVERAGE_FLOOR,
    ICL_COVERAGE_FLOOR, read_tsv, selected_pairs, sha256,
)

DEFAULT_OUTPUT = ROOT / "satellite_icl_shared_pixel_covariance_2026-10-18_01"
DEFAULT_LAG_COVARIANCE = ROOT / "icl_coadd_noise_correlation_2026-10-19_04"
BIN_INDEX = len(satellite.R_EDGES_KPC) - 2


def supports_outer_annulus(profile: dict | None) -> bool:
    """Match the per-annulus estimator: only 50--75 kpc coverage gates support."""
    return bool(profile is not None
                and profile["coverage"][BIN_INDEX] >= SAT_COVERAGE_FLOOR)


def annulus_pixel_indices(center_xy: np.ndarray, target_number: int | None,
                          image_shape: tuple[int, int], base: np.ndarray,
                          segmentation: np.ndarray, detections: dict,
                          fallback_tree, pixel_scale: float, kpc_per_pixel: float,
                          finite: np.ndarray) -> np.ndarray:
    """Return valid 50--75 kpc pixel indices using the profile estimator mask."""
    cx, cy = map(float, center_xy)
    half = int(math.ceil(SAT_OUTER / kpc_per_pixel)) + 2
    x0, x1 = max(0, int(math.floor(cx-half))), min(image_shape[1], int(math.ceil(cx+half))+1)
    y0, y1 = max(0, int(math.floor(cy-half))), min(image_shape[0], int(math.ceil(cy+half))+1)
    if x0 >= x1 or y0 >= y1:
        return np.empty(0, dtype=np.int64)
    mask = satellite.make_stamp_mask(base, segmentation, detections, fallback_tree,
                                     (cx, cy), x0, y0, target_number,
                                     pixel_scale, half, satellite.KRON_FACTOR)
    yy, xx = np.indices((y1-y0, x1-x0), dtype=np.float32)
    radius = np.hypot(xx+x0-cx, yy+y0-cy)*kpc_per_pixel
    keep = ((radius >= SAT_INNER) & (radius < SAT_OUTER) & ~mask
            & finite[y0:y1, x0:x1])
    ly, lx = np.nonzero(keep)
    return np.ravel_multi_index((ly+y0, lx+x0), image_shape).astype(np.int64)


def diagonal_covariance(pixel_weights: dict[int, float], icl_indices: np.ndarray,
                        variance_flat: np.ndarray, conversion: float) -> dict:
    """Compute ICL-target/satellite covariance under diagonal pixel variance."""
    keys = np.fromiter(pixel_weights.keys(), dtype=np.int64)
    weights = np.fromiter(pixel_weights.values(), dtype=np.float64)
    active = np.abs(weights) > 1e-15
    keys, weights = keys[active], weights[active]
    icl_valid = np.zeros(len(variance_flat), dtype=bool)
    icl_valid[np.asarray(icl_indices, dtype=np.int64)] = True
    shared = icl_valid[keys]
    n_icl = len(icl_indices)
    var_sat = float(np.sum(variance_flat[keys] * weights**2) * conversion**2)
    var_icl = float(np.sum(variance_flat[icl_indices]) / n_icl**2 * conversion**2)
    cov = float(np.sum(variance_flat[keys[shared]] * weights[shared])
                / n_icl * conversion**2)
    corr = cov / math.sqrt(var_icl*var_sat) if var_icl > 0 and var_sat > 0 else None
    return {"n_satellite_weight_pixels": int(len(keys)),
            "n_shared_nonzero_weight_pixels": int(shared.sum()),
            "diagonal_var_icl_target_njy2_arcsec4": var_icl,
            "diagonal_var_satellite_unweighted_njy2_arcsec4": var_sat,
            "diagonal_covariance_njy2_arcsec4": cov,
            "diagonal_pixel_noise_correlation": corr}


def lag_weighted_covariance(weight_a: np.ndarray, weight_b: np.ndarray,
                            variance: np.ndarray, kernel: np.ndarray,
                            conversion: float = 1.0) -> float:
    """Contract signed aperture weights with an empirical stationary lag kernel.

    The kernel is the covariance divided by ``sqrt(var_i * var_j)``. Its
    central pixel is at ``kernel.shape // 2``. Zero padding treats pixels
    outside the image as absent.
    """
    wa, wb = np.asarray(weight_a, dtype=np.float64), np.asarray(weight_b, dtype=np.float64)
    var, cov_kernel = np.asarray(variance, dtype=np.float64), np.asarray(kernel, dtype=np.float64)
    if wa.ndim != 2 or wa.shape != wb.shape or wa.shape != var.shape:
        raise ValueError("weight maps and variance must be matching 2D arrays")
    if (cov_kernel.ndim != 2 or cov_kernel.shape[0] != cov_kernel.shape[1]
            or cov_kernel.shape[0] % 2 != 1 or not np.all(np.isfinite(cov_kernel))):
        raise ValueError("lag kernel must be a finite odd-sized square array")
    weighted = (wa != 0) | (wb != 0)
    if (not np.all(np.isfinite(wa)) or not np.all(np.isfinite(wb))
            or np.any((~np.isfinite(var) | (var <= 0)) & weighted)
            or not math.isfinite(float(conversion)) or conversion <= 0):
        raise ValueError("nonfinite weights or nonpositive variance on weighted pixels")
    sqrt_variance = np.sqrt(np.where(np.isfinite(var) & (var > 0), var, 0.0))
    a, b = wa * sqrt_variance, wb * sqrt_variance
    return float(np.sum(a * fftconvolve(b, cov_kernel, mode="same"), dtype=np.float64)
                 * float(conversion) ** 2)


def load_lag_covariance(directory: Path) -> tuple[dict, dict]:
    """Read and checksum-verify the matched target/control lag-kernel audit."""
    summary_path = directory / "summary.json"
    table_path = directory / "lag_correlation.tsv"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    if (summary.get("schema") != "astrafek-hsc-coadd-noise-correlation-v1"
            or summary.get("status") != "complete_exploratory"
            or summary.get("science_measurement") is not False
            or summary.get("n_target_frames") != 57
            or summary.get("n_matched_control_frames") != 57):
        raise ValueError("expected the audited matched 57-frame target/control lag kernel")
    if sha256(table_path) != summary.get("lag_results_sha256"):
        raise ValueError("lag-correlation table checksum mismatch")
    for path, expected in summary.get("input_sha256", {}).items():
        if sha256(Path(path)) != expected:
            raise ValueError(f"lag-kernel input checksum mismatch: {path}")
    analysis_script = Path(__file__).with_name("estimate_hsc_coadd_noise_correlation.py")
    if sha256(analysis_script) != summary.get("analysis_script_sha256"):
        raise ValueError("lag-kernel analysis script changed after measurement")
    grouped: dict[tuple[str, str], dict[tuple[int, int], float]] = {}
    for row in read_tsv(table_path):
        if row["role"] != "target":
            continue
        key = (str(row["cluster_id"]), str(row["filter"]))
        lag = (int(row["dy"]), int(row["dx"]))
        value = float(row["median_variance_normalized_covariance"])
        if not math.isfinite(value) or lag in grouped.setdefault(key, {}):
            raise ValueError(f"invalid or duplicate lag covariance for {key} {lag}")
        grouped[key][lag] = value
    if len(grouped) != 57:
        raise ValueError(f"expected 57 target kernels, got {len(grouped)}")
    width = int(summary["max_lag_pixels"])
    expected_lags = {(dy, dx) for dy in range(-width, width + 1)
                     for dx in range(-width, width + 1)}
    kernels = {}
    for key, values in grouped.items():
        if set(values) != expected_lags:
            raise ValueError(f"incomplete lag grid for {key}")
        kernel = np.empty((2 * width + 1, 2 * width + 1), dtype=np.float64)
        for dy, dx in expected_lags:
            kernel[dy + width, dx + width] = 0.5 * (values[(dy, dx)] + values[(-dy, -dx)])
        kernels[key] = kernel
    provenance = {"summary_path": str(summary_path.resolve()),
                  "summary_sha256": sha256(summary_path),
                  "lag_table_path": str(table_path.resolve()),
                  "lag_table_sha256": sha256(table_path)}
    return kernels, provenance


def process_cluster(ident: str, rows: list[dict[str, str]], data: Path,
                    masks_dir: Path, prep_dir: Path, manifest: dict,
                    mask_summary: dict, mask_records: list[dict[str, str]],
                    icl_pairs: dict[tuple[str, str], dict[str, str]],
                    lag_kernel: np.ndarray | None = None) -> dict:
    mask_matches = [row for row in mask_records if row["id"] == ident]
    if len(mask_matches) != 1:
        raise ValueError(f"{ident}: expected one saved ICL mask record")
    mask_record = mask_matches[0]
    filt = mask_record["primary_filter"]
    rows = [row for row in rows if row["filter"] == filt]
    if not rows:
        raise ValueError(f"{ident}: no matched satellite rows for the ICL filter")
    pair = icl_pairs.get((ident, filt))
    if pair is None or float(pair["coverage_200_300"]) < ICL_COVERAGE_FLOOR:
        raise ValueError(f"{ident}: missing supported ICL annulus pair")
    frames = [row for row in manifest["clusters"] if str(row["id"]) == ident]
    if len(frames) != 1:
        raise ValueError(f"{ident}: manifest identity is not unique")
    frame = dict(frames[0]); image_path = Path(mask_record["image"]); frame["image"] = str(image_path)
    icl_mask_path = Path(mask_record["mask_file"])
    if sha256(image_path) != mask_record["image_sha256"]:
        raise ValueError(f"{ident}: image checksum mismatch")
    if (sha256(icl_mask_path) != mask_record["mask_sha256"]
            or mask_summary.get("mask_output_sha256", {}).get(str(icl_mask_path))
            != mask_record["mask_sha256"]):
        raise ValueError(f"{ident}: ICL mask checksum mismatch")
    prep = prep_dir / f"cluster_{ident}"
    detections_list = satellite.parse_detections(satellite.read_tsv(prep/"detections.tsv"))
    member_rows = satellite.read_tsv(prep/"members.tsv")
    preparation = satellite.read_json(prep/"preparation.json")
    if preparation.get("science_measurement") is not False:
        raise ValueError(f"{ident}: unexpected satellite preparation science gate")
    info = satellite.validate_fits(image_path, frame)
    if info["coadd_header_filter"] != filt:
        raise ValueError(f"{ident}: filter mismatch")
    matched, _ = satellite.match_members(member_rows, detections_list, info["pixel_scale_arcsec"])
    matched_by_index = {item["member_index"]: item for item in matched}
    with fits.open(image_path, memmap=False, checksum=True) as hdul:
        image = np.asarray(hdul[info["science_hdu"]].data, dtype=np.float32)
        variance = np.asarray(hdul[info["variance_hdu"]].data, dtype=np.float32)
    with fits.open(icl_mask_path, memmap=False, checksum=True) as hdul:
        icl_excluded = np.asarray(hdul["MASK_K2P5"].data, dtype=bool)
    with fits.open(prep/"masks.fits.gz", memmap=False, checksum=True) as hdul:
        segmentation = np.asarray(hdul["SEGMENTATION"].data, dtype=np.int32)
        base = np.asarray(hdul["BASE_MASK"].data, dtype=bool)
    if not (image.shape == variance.shape == icl_excluded.shape == segmentation.shape == base.shape):
        raise ValueError(f"{ident}: image/mask shapes differ")
    finite = np.isfinite(image) & np.isfinite(variance) & (variance > 0)
    if np.any(~finite & ~base):
        raise ValueError(f"{ident}: BASE_MASK does not include all invalid pixels")
    pixel_scale = float(info["pixel_scale_arcsec"])
    center_xy = np.asarray(info["center_pixel_xy"], dtype=float)
    kpc_per_pixel = float(FlatLambdaCDM(
        H0=float(manifest.get("h0", 70.)), Om0=float(manifest.get("omega_m", .3)),
        Tcmb0=0*u.K).kpc_proper_per_arcmin(float(frame["z"])).value*pixel_scale/60.)
    dets = satellite.detection_arrays(detections_list, pixel_scale)
    member_in_frame = [r for r in member_rows if r["IN_FRAME"] == "True"]
    member_xy = np.asarray([[float(r["X_IMAGE"])-1., float(r["Y_IMAGE"])-1.]
                            for r in member_in_frame], dtype=float).reshape((-1,2))
    matched_indices = set(matched_by_index)
    fallback_xy = np.asarray([[float(r["X_IMAGE"])-1., float(r["Y_IMAGE"])-1.]
                              for r in member_in_frame if int(r["INDEX"]) not in matched_indices],
                             dtype=float).reshape((-1,2))
    fallback_tree = cKDTree(fallback_xy) if len(fallback_xy) else None
    member_tree = cKDTree(member_xy) if len(member_xy) else None
    yy, xx = np.ogrid[:image.shape[0], :image.shape[1]]
    radius_cluster = np.hypot((xx-center_xy[0])*kpc_per_pixel,
                              (yy-center_xy[1])*kpc_per_pixel)
    icl_valid = ~icl_excluded & finite
    icl_annulus = (radius_cluster >= ICL_INNER) & (radius_cluster < ICL_OUTER) & icl_valid
    icl_indices = np.flatnonzero(icl_annulus)
    if not len(icl_indices):
        raise ValueError(f"{ident}: ICL annulus has no valid pixels")

    pixel_weights: dict[int, float] = {}
    n_satellites = len(rows)
    control_overlap_incidence = 0
    target_overlap_incidence = 0
    for row in rows:
        member_index, number = int(row["member_index"]), int(row["detection_number"])
        match = matched_by_index.get(member_index)
        if match is None or int(detections_list[match["detection_index"]]["NUMBER"]) != number:
            raise ValueError(f"{ident}: profile/member/detection identity mismatch")
        target_det = detections_list[match["detection_index"]]
        target_xy = np.array([target_det["x"], target_det["y"]])
        target_profile = satellite.profile_at(
            image, variance, base, segmentation, dets, fallback_tree,
            target_xy, pixel_scale, kpc_per_pixel, 1.0, number, SAT_COVERAGE_FLOOR)
        if not supports_outer_annulus(target_profile):
            raise ValueError(f"{ident}/{member_index}: target no longer passes saved support")
        target_ids = annulus_pixel_indices(
            target_xy, number, image.shape, base, segmentation, dets, fallback_tree,
            pixel_scale, kpc_per_pixel, finite)
        if len(target_ids) != int(target_profile["n_pixels"][BIN_INDEX]) \
                or len(target_ids) != int(row["target_n_pixels"]):
            raise ValueError(f"{ident}/{member_index}: reconstructed target support count mismatch")
        weight = 1.0/(n_satellites*len(target_ids))
        for pixel in target_ids:
            pixel_weights[int(pixel)] = pixel_weights.get(int(pixel), 0.0)+weight
        target_overlap_incidence += int(np.count_nonzero(icl_annulus.ravel()[target_ids]))

        cluster_radius_kpc = float(row["cluster_radius_kpc"])
        points = satellite.control_centers(
            target_xy, center_xy, cluster_radius_kpc/kpc_per_pixel, dets,
            member_tree, pixel_scale, kpc_per_pixel, image.shape)
        accepted = []
        for point in points:
            profile = satellite.profile_at(
                image, variance, base, segmentation, dets, fallback_tree,
                point, pixel_scale, kpc_per_pixel, 1.0, None, SAT_COVERAGE_FLOOR)
            if supports_outer_annulus(profile):
                accepted.append((point, profile))
        if len(accepted) != int(row["n_controls"]):
            raise ValueError(f"{ident}/{member_index}: reconstructed {len(accepted)} controls, "
                             f"saved table records {row['n_controls']}")
        for point, profile in accepted:
            control_ids = annulus_pixel_indices(
                point, None, image.shape, base, segmentation, dets, fallback_tree,
                pixel_scale, kpc_per_pixel, finite)
            if len(control_ids) != int(profile["n_pixels"][BIN_INDEX]):
                raise ValueError(f"{ident}/{member_index}: reconstructed control support mismatch")
            cweight = -1.0/(n_satellites*len(accepted)*len(control_ids))
            for pixel in control_ids:
                pixel_weights[int(pixel)] = pixel_weights.get(int(pixel), 0.0)+cweight
            control_overlap_incidence += int(np.count_nonzero(icl_annulus.ravel()[control_ids]))

    variance_flat = variance.ravel().astype(np.float64)
    conversion = float(satellite.nJy_factor(info))
    result = diagonal_covariance(pixel_weights, icl_indices, variance_flat, conversion)
    if lag_kernel is not None:
        icl_weights = np.zeros(image.shape, dtype=np.float64)
        icl_weights.ravel()[icl_indices] = 1.0 / len(icl_indices)
        satellite_weights = np.zeros(image.shape, dtype=np.float64)
        sat_keys = np.fromiter(pixel_weights.keys(), dtype=np.int64)
        sat_values = np.fromiter(pixel_weights.values(), dtype=np.float64)
        satellite_weights.ravel()[sat_keys] = sat_values
        variance_image = np.asarray(variance, dtype=np.float64)
        lag_var_icl = lag_weighted_covariance(
            icl_weights, icl_weights, variance_image, lag_kernel, conversion,
        )
        lag_var_sat = lag_weighted_covariance(
            satellite_weights, satellite_weights, variance_image, lag_kernel, conversion,
        )
        lag_cov = lag_weighted_covariance(
            icl_weights, satellite_weights, variance_image, lag_kernel, conversion,
        )
        if lag_var_icl <= 0 or lag_var_sat <= 0:
            raise ValueError(f"{ident}: measured lag kernel gives nonpositive aperture variance")
        result.update({
            "lag_kernel_var_icl_target_njy2_arcsec4": lag_var_icl,
            "lag_kernel_var_satellite_unweighted_njy2_arcsec4": lag_var_sat,
            "lag_kernel_covariance_njy2_arcsec4": lag_cov,
            "lag_kernel_pixel_noise_correlation": lag_cov / math.sqrt(lag_var_icl * lag_var_sat),
            "lag_kernel_icl_sigma_inflation_vs_diagonal": math.sqrt(
                lag_var_icl / result["diagonal_var_icl_target_njy2_arcsec4"]),
            "lag_kernel_satellite_sigma_inflation_vs_diagonal": math.sqrt(
                lag_var_sat / result["diagonal_var_satellite_unweighted_njy2_arcsec4"]),
        })
    result.update({"cluster_id": ident, "filter": filt,
                   "n_satellite_profiles": n_satellites,
                   "n_icl_valid_pixels": int(len(icl_indices)),
                   "target_shared_pixel_incidence": target_overlap_incidence,
                   "control_shared_pixel_incidence": control_overlap_incidence,
                   "icl_delta_sb_observed_njy_arcsec2": float(pair["delta_sb_200_300"]),
                   "satellite_unweighted_delta_observed_njy_arcsec2": float(np.mean([
                       float(row["delta_sb_njy_arcsec2"]) for row in rows])),
                   "input_image_sha256": sha256(image_path),
                   "icl_mask_sha256": sha256(icl_mask_path),
                   "satellite_mask_sha256": sha256(prep/"masks.fits.gz")})
    return result


def run(data: Path, masks_dir: Path, prep_dir: Path, profile_dir: Path,
        icl_pairs_path: Path, output_dir: Path,
        cluster_ids: list[str] | None = None,
        lag_covariance_dir: Path | None = None) -> dict:
    if output_dir.exists():
        raise FileExistsError(f"refusing to overwrite {output_dir}")
    manifest_path = data/"sample_manifest.json"
    mask_summary_path, mask_records_path = masks_dir/"summary.json", masks_dir/"records.tsv"
    profile_summary_path = profile_dir/"summary.json"
    profile_table_path = profile_dir/"satellite_profiles.tsv"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    mask_summary = json.loads(mask_summary_path.read_text(encoding="utf-8"))
    profile_summary = json.loads(profile_summary_path.read_text(encoding="utf-8"))
    if (mask_summary.get("science_measurement") is not False
            or profile_summary.get("science_measurement") is not False
            or profile_summary.get("n_frames_completed") != 800):
        raise ValueError("expected complete exploratory 800-frame inputs")
    grouped = selected_pairs(icl_pairs_path, profile_table_path)
    if cluster_ids:
        missing = sorted(set(cluster_ids)-set(grouped))
        if missing:
            raise ValueError(f"requested clusters are outside the matched support: {missing}")
        grouped = {ident: grouped[ident] for ident in cluster_ids}
    mask_records = read_tsv(mask_records_path)
    primary_filter = {str(row["id"]): str(row["primary_filter"]) for row in mask_records}
    lag_kernels = None
    lag_provenance = None
    if lag_covariance_dir is not None:
        lag_kernels, lag_provenance = load_lag_covariance(lag_covariance_dir)
        supported_keys = {(ident, primary_filter[ident]) for ident in grouped}
        selected_kernel_keys = supported_keys & set(lag_kernels)
        if len(selected_kernel_keys) != 57:
            raise ValueError(f"expected exactly 57 supported target/kernel pairs; got {len(selected_kernel_keys)}")
        grouped = {ident: rows for ident, rows in grouped.items()
                   if (ident, primary_filter[ident]) in selected_kernel_keys}
    icl_pairs = {(r["id"], r["filter"]): r for r in read_tsv(icl_pairs_path)}
    results = []
    for index, (ident, rows) in enumerate(sorted(grouped.items()), 1):
        print(f"[{index}/{len(grouped)}] pixel covariance cluster {ident}", flush=True)
        filt = primary_filter[ident]
        if not any(row["filter"] == filt for row in rows):
            raise ValueError(f"{ident}: selected satellite rows lack primary filter {filt}")
        kernel = lag_kernels[(ident, filt)] if lag_kernels is not None else None
        results.append(process_cluster(ident, rows, data, masks_dir, prep_dir, manifest,
                                       mask_summary, mask_records, icl_pairs, kernel))
    output_dir.mkdir(parents=True, exist_ok=False)
    table_path = output_dir/"pair_pixel_covariance.tsv"
    with table_path.open("x", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(results[0]), delimiter="\t")
        writer.writeheader(); writer.writerows(results)
    by_filter = {}
    for filt in sorted({r["filter"] for r in results}):
        subset = [r for r in results if r["filter"] == filt]
        corrs = [r["diagonal_pixel_noise_correlation"] for r in subset]
        covs = [r["diagonal_covariance_njy2_arcsec4"] for r in subset]
        by_filter[filt] = {"n_cluster_filter_pairs": len(subset),
                           "median_diagonal_pixel_noise_correlation": float(np.median(corrs)),
                           "q16_q84_diagonal_pixel_noise_correlation": np.quantile(corrs,[.16,.84]).tolist(),
                           "median_diagonal_covariance_njy2_arcsec4": float(np.median(covs)),
                           "negative_covariance_pairs": int(np.count_nonzero(np.asarray(covs)<0))}
    input_paths = [manifest_path, mask_summary_path, mask_records_path, profile_summary_path,
                   profile_table_path, icl_pairs_path, Path(__file__).resolve()]
    if lag_provenance is not None:
        input_paths.extend((Path(lag_provenance["summary_path"]),
                            Path(lag_provenance["lag_table_path"])))
    summary = {
        "schema": ("hsc-icl-satellite-shared-pixel-lag-covariance-v1"
                   if lag_kernels is not None else "hsc-icl-satellite-shared-pixel-diagonal-covariance-v1"),
        "status": "complete_exploratory", "science_measurement": False,
        "coadd_covariance_available": False,
        "short_range_empirical_lag_kernel_applied": lag_kernels is not None,
        "method": ("Reconstruct actual target and accepted azimuthal-control pixel weights. Contract the ICL 200-300 kpc and satellite 50-75 kpc signed estimators with per-frame symmetrized 9x9 empirical covariance kernels measured from matched target/control image patches."
                   if lag_kernels is not None else "Reconstruct target and accepted azimuthal-control 50-75 kpc masks using the saved satellite estimator; build equal-satellite signed pixel weights and contract them with the ICL target 200-300 kpc weights and HSC variance map under diagonal independent-pixel noise."),
        "selection": {"n_cluster_filter_pairs": len(results),
                      "n_satellite_profiles": sum(r["n_satellite_profiles"] for r in results),
                      "icl_annulus_kpc": [ICL_INNER, ICL_OUTER],
                      "satellite_annulus_kpc": [SAT_INNER, SAT_OUTER],
                      "icl_coverage_floor": ICL_COVERAGE_FLOOR,
                      "satellite_target_control_coverage_floor": SAT_COVERAGE_FLOOR},
        "by_filter": by_filter,
        "lag_kernel_provenance": lag_provenance,
        "lag_covariance_by_filter": ({filt: {
            "n_cluster_filter_pairs": sum(r["filter"] == filt for r in results),
            "median_icl_sigma_inflation_vs_diagonal": float(np.median([
                r["lag_kernel_icl_sigma_inflation_vs_diagonal"] for r in results if r["filter"] == filt])),
            "median_satellite_sigma_inflation_vs_diagonal": float(np.median([
                r["lag_kernel_satellite_sigma_inflation_vs_diagonal"] for r in results if r["filter"] == filt])),
            "median_lag_kernel_pixel_noise_correlation": float(np.median([
                r["lag_kernel_pixel_noise_correlation"] for r in results if r["filter"] == filt])),
            "q16_q84_lag_kernel_pixel_noise_correlation": np.quantile([
                r["lag_kernel_pixel_noise_correlation"] for r in results if r["filter"] == filt], [.16,.84]).tolist(),
        } for filt in sorted({r["filter"] for r in results})} if lag_kernels is not None else {}),
        "input_sha256": {str(p): sha256(p) for p in input_paths},
        "output_sha256": {str(table_path): sha256(table_path)},
        "limitations": [
            "The measured lag kernel is a stationary, finite-lag patch estimate; it does not capture the exact position-dependent full coadd covariance.",
            "The ICL random-control term is not included; this is covariance with the target-side ICL annulus only.",
            "The diagnostic is conditional on fixed detections, masks, support cuts and equal cluster/satellite weighting; it does not correct the observed ICL-satellite relation.",
            "Pixel covariance is a noise-only sensitivity and excludes source-Poisson, sky-transfer, PSF, calibration and astrophysical covariance.",
        ],
    }
    (output_dir/"summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True)+"\n",
                                             encoding="utf-8")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=ROOT)
    parser.add_argument("--icl-masks", type=Path, default=DEFAULT_MASKS)
    parser.add_argument("--satellite-prep", type=Path, default=DEFAULT_PREP)
    parser.add_argument("--satellite-profiles", type=Path, default=DEFAULT_PROFILES)
    parser.add_argument("--icl-pairs", type=Path, default=DEFAULT_ICL)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--lag-covariance-dir", type=Path,
                        help="apply a checksum-verified measured short-range coadd lag kernel")
    parser.add_argument("--cluster-id", action="append", default=[])
    args = parser.parse_args()
    summary = run(args.data.resolve(), args.icl_masks.resolve(), args.satellite_prep.resolve(),
                  args.satellite_profiles.resolve(), args.icl_pairs.resolve(), args.output_dir.resolve(),
                  args.cluster_id or None,
                  args.lag_covariance_dir.resolve() if args.lag_covariance_dir else None)
    print(json.dumps({"status": summary["status"], "science_measurement": False,
                      "selection": summary["selection"], "by_filter": summary["by_filter"],
                      "output": str(args.output_dir.resolve())}, indent=2))


if __name__ == "__main__":
    main()
