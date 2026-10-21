#!/usr/bin/env python3
"""Inject synthetic ICL/satellite light convolved with a public wide HSC PSF.

This is a separate sensitivity path. It applies the published i-band reference
PSF only to synthetic components; it does not claim a target-specific PSF or
alter the observed coadd.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
import platform
import statistics
import time
from pathlib import Path

import numpy as np
from astropy import units as u
from astropy.cosmology import FlatLambdaCDM
from astropy.io import fits

import measure_satellite_outer_profiles as satellite
import prestack_core
from joint_icl_satellite_injection import (
    DEFAULT_DATA, DEFAULT_ICL_MASKS, DEFAULT_ICL_PAIRS, DEFAULT_SAT_PREP,
    add_surface_brightness_injection, redetect_masks_and_targets,
    satellite_halo_template, satellite_outer_measurement, sha256,
)
from position_psf import (PositionPSFBank, build_core_wing_hybrid_psf,
                          convolve_surface_brightness,
                          satellite_halos_with_psfs)
from run_joint_icl_satellite_batch import (
    DEFAULT_SATELLITE_TABLE, select_jointly_supported_pairs,
)


HERE = Path(__file__).resolve().parent
DEFAULT_REFERENCE_DIR = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/extended_psf_reference_2026-10-14_01")
DEFAULT_PLAN_DIR = DEFAULT_DATA / "icl_psf_picker_injection_plan_2026-10-19_02"
DEFAULT_PSFS = DEFAULT_DATA / "icl_psf_picker_injection_2026-10-19_01"
DEFAULT_TEST_CENTER = DEFAULT_DATA / "icl_psf_picker_test_2026-10-19_01/cluster10_i2_center_psf.fits"


def load_public_reference(path: Path, reference_dir: Path) -> tuple[np.ndarray, dict, Path]:
    """Load a checksum-verified HSC PDR3 i-band extended PSF without renormalizing."""
    path, reference_dir = Path(path).resolve(), Path(reference_dir).resolve()
    summary_path = reference_dir / "summary.json"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    if summary.get("status") != "complete" or summary.get("license") != "CC-BY-4.0":
        raise ValueError("published PSF reference metadata or license is unexpected")
    if summary.get("psf_approved") is not False:
        raise ValueError("this sensitivity runner requires the reference to remain unapproved")
    rows = [row for row in summary.get("records", [])
            if Path(row.get("path", "")).resolve() == path]
    if len(rows) != 1:
        raise ValueError("PSF file is not uniquely listed in the audited reference summary")
    record = rows[0]
    if record.get("shape") != [4001, 4001] or not math.isclose(
            float(record.get("pixel_scale_arcsec", 0)), 0.168, rel_tol=0, abs_tol=1e-9):
        raise ValueError("unexpected extended PSF shape or pixel scale")
    digest = sha256(path)
    if digest != record.get("sha256"):
        raise ValueError("extended PSF checksum differs from its audited reference record")
    with fits.open(path, memmap=False, checksum=True) as hdus:
        hdu = int(record["hdu"])
        if len(hdus) <= hdu or hdus[hdu].data is None:
            raise ValueError("extended PSF HDU is missing")
        kernel = np.asarray(hdus[hdu].data, dtype=np.float32)
    total = float(kernel.sum(dtype=np.float64))
    if (kernel.shape != (4001, 4001) or not np.isfinite(kernel).all()
            or abs(total - float(record["finite_grid_sum"])) > 2e-6):
        raise ValueError("extended PSF array differs from audited finite-grid properties")
    if not np.any(kernel < 0):
        raise ValueError("signed interpolation values were unexpectedly clipped")
    return kernel, {**record, "path": str(path), "sha256": digest, "license": summary["license"],
                    "citation": "Garate-Nuñez et al. 2024; Zenodo DOI 10.5281/zenodo.14915578",
                    "finite_grid_sum_loaded_float32": total,
                    "summary_sha256": sha256(summary_path),
                    "reference_dir": str(reference_dir)}, summary_path


def satellite_halos_truncated(shape: tuple[int, int], centers_xy: np.ndarray,
                              kpc_per_pixel: float, amplitude_at_50: float,
                              scale_kpc: float, truncation_scale_lengths: float = 15.0) -> np.ndarray:
    """Reproduce the local-core run's 15-scale-length satellite-envelope support."""
    centers = np.asarray(centers_xy, dtype=float).reshape((-1, 2))
    if kpc_per_pixel <= 0 or scale_kpc <= 0 or amplitude_at_50 < 0:
        raise ValueError("scales must be positive and amplitude non-negative")
    radius = int(np.ceil(truncation_scale_lengths * scale_kpc / kpc_per_pixel))
    axis = np.arange(-radius, radius + 1, dtype=np.float32)
    yy, xx = np.meshgrid(axis, axis, indexing="ij")
    output = np.zeros(shape, dtype=np.float32)
    for x, y in centers:
        cx, cy = int(round(float(x))), int(round(float(y)))
        frac_x, frac_y = float(x) - cx, float(y) - cy
        radial = np.hypot(xx - frac_x, yy - frac_y) * kpc_per_pixel
        halo = amplitude_at_50 * np.exp((50.0 - radial) / scale_kpc)
        x0, x1 = max(0, cx - radius), min(shape[1], cx + radius + 1)
        y0, y1 = max(0, cy - radius), min(shape[0], cy + radius + 1)
        if x0 >= x1 or y0 >= y1:
            continue
        sx0, sx1 = x0 - (cx - radius), x1 - (cx - radius)
        sy0, sy1 = y0 - (cy - radius), y1 - (cy - radius)
        output[y0:y1, x0:x1] += halo[sy0:sy1, sx0:sx1]
    return output


def unit_point_source_surface_brightness(shape: tuple[int, int], center_xy,
                                         pixel_scale_arcsec: float,
                                         total_flux_njy: float = 1.0) -> np.ndarray:
    """Build a bilinearly centered point source in nJy/arcsec², conserving flux."""
    x, y = map(float, center_xy)
    if (len(shape) != 2 or min(shape) <= 0
            or not np.isfinite([x, y, pixel_scale_arcsec, total_flux_njy]).all()
            or pixel_scale_arcsec <= 0 or total_flux_njy <= 0):
        raise ValueError("point-source geometry and flux must be finite and positive")
    x0, y0 = int(math.floor(x)), int(math.floor(y))
    if x0 < 0 or y0 < 0 or x0 + 1 >= shape[1] or y0 + 1 >= shape[0]:
        raise ValueError("point-source center must leave room for bilinear centering")
    fx, fy = x - x0, y - y0
    weights = ((y0, x0, (1-fx)*(1-fy)),
               (y0, x0+1, fx*(1-fy)),
               (y0+1, x0, (1-fx)*fy),
               (y0+1, x0+1, fx*fy))
    source = np.zeros(shape, dtype=np.float32)
    pixel_area = pixel_scale_arcsec**2
    for yy, xx, weight in weights:
        source[yy, xx] += total_flux_njy * weight / pixel_area
    return source


def measure_point_source_wing_profile(surface_brightness: np.ndarray,
                                      excluded: np.ndarray, variance: np.ndarray,
                                      radius_kpc: np.ndarray, kpc_per_pixel: float,
                                      pixel_scale_arcsec: float,
                                      total_flux_njy: float) -> dict:
    """Measure masked point-source leakage in fixed physical annuli."""
    image = np.asarray(surface_brightness, dtype=float)
    mask = np.asarray(excluded, dtype=bool)
    var = np.asarray(variance, dtype=float)
    radius = np.asarray(radius_kpc, dtype=float)
    if not (image.shape == mask.shape == var.shape == radius.shape):
        raise ValueError("point-source image, mask, variance, and radius shapes differ")
    if (not np.isfinite(image).all() or not np.isfinite(radius).all()
            or not np.isfinite([kpc_per_pixel, pixel_scale_arcsec, total_flux_njy]).all()
            or kpc_per_pixel <= 0 or pixel_scale_arcsec <= 0 or total_flux_njy <= 0):
        raise ValueError("point-source profile inputs must be finite and positive")
    valid = ~mask & np.isfinite(var) & (var > 0)
    result = {}
    pixel_area_arcsec2 = pixel_scale_arcsec**2
    for label, (lo, hi) in {"50_100_kpc": (50., 100.),
                            "200_300_kpc": (200., 300.)}.items():
        annulus = (radius >= lo) & (radius < hi)
        use = annulus & valid
        count = int(use.sum())
        expected = math.pi * (hi**2 - lo**2) / kpc_per_pixel**2
        sb_sum = float(image[use].sum(dtype=np.float64))
        sb_mean = sb_sum / count if count else None
        result[label] = {
            "n_valid_pixels": count,
            "geometric_pixel_count": expected,
            "mask_coverage": min(1.0, count / expected) if expected > 0 else 0.0,
            "mean_surface_brightness_njy_arcsec2_per_njy": (
                sb_mean / total_flux_njy if sb_mean is not None else None),
            "masked_annulus_flux_fraction": sb_sum * pixel_area_arcsec2 / total_flux_njy,
        }
    return {"source_flux_njy": float(total_flux_njy),
            "pixel_scale_arcsec": float(pixel_scale_arcsec),
            "kpc_per_pixel": float(kpc_per_pixel),
            "annuli": result}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def select_reference_frames(pairs, supported_frames, filter_name: str):
    """Select the paired-injection frames with an approved local core PSF mapping."""
    if filter_name not in ("i", "i2"):
        raise ValueError("filter_name must be i or i2")
    expected_count = {"i": 13, "i2": 44}[filter_name]
    selected = [(str(ident), str(filt)) for ident, filt in pairs
                if str(filt) == filter_name and (str(ident), str(filt)) in supported_frames]
    if len(selected) != expected_count:
        raise ValueError(
            f"expected {expected_count} {filter_name}-band frames with approved injection positions, got {len(selected)}"
        )
    return selected


def calculate_transfer_response(measurements: dict, diffuse_amplitude: float,
                                satellite_amplitude: float) -> dict:
    """Calculate unit responses and additivity residuals for four scenarios."""
    if (not math.isfinite(diffuse_amplitude) or diffuse_amplitude <= 0
            or not math.isfinite(satellite_amplitude) or satellite_amplitude <= 0):
        raise ValueError("injection amplitudes must be finite and positive")
    base_icl = float(measurements["baseline"]["icl_raw_200_300_njy_arcsec2"])
    base_sat = float(measurements["baseline"]["satellite_delta_50_75_njy_arcsec2"])
    icl_only_icl = float(measurements["icl_only"]["icl_raw_200_300_njy_arcsec2"])
    icl_only_sat = float(measurements["icl_only"]["satellite_delta_50_75_njy_arcsec2"])
    sat_only_icl = float(measurements["satellite_only"]["icl_raw_200_300_njy_arcsec2"])
    sat_only_sat = float(measurements["satellite_only"]["satellite_delta_50_75_njy_arcsec2"])
    both_icl = float(measurements["both"]["icl_raw_200_300_njy_arcsec2"])
    both_sat = float(measurements["both"]["satellite_delta_50_75_njy_arcsec2"])
    return {
        "icl_response_to_diffuse_unit_injection": (icl_only_icl - base_icl) / diffuse_amplitude,
        "icl_response_to_all_satellite_envelopes": (sat_only_icl - base_icl) / satellite_amplitude,
        "satellite_response_to_smooth_diffuse": (icl_only_sat - base_sat) / diffuse_amplitude,
        "satellite_response_to_satellite_unit_injection": (sat_only_sat - base_sat) / satellite_amplitude,
        "both_minus_sum_residual_icl_njy_arcsec2":
            both_icl - base_icl - (icl_only_icl - base_icl) - (sat_only_icl - base_icl),
        "both_minus_sum_residual_satellite_njy_arcsec2":
            both_sat - base_sat - (icl_only_sat - base_sat) - (sat_only_sat - base_sat),
    }


def run_frame(data: Path, ident: str, output: Path, reference: dict,
              kernel: np.ndarray, icl_masks_dir: Path, satellite_prep_dir: Path,
              icl_pairs_path: Path, psf_bank: PositionPSFBank,
              diffuse_amplitude: float = 1.0, satellite_amplitude: float = 1.0,
              diffuse_scale_kpc: float = 100.0, satellite_scale_kpc: float = 15.0,
              unresolved_bcg_policy: str = "mask_all", filter_name: str = "i",
              hybrid_parameters: dict | None = None) -> dict:
    if filter_name not in ("i", "i2"):
        raise ValueError("filter_name must be i or i2")
    if output.exists():
        raise FileExistsError(f"refusing to overwrite {output}")
    if (not platform.node().lower().startswith("lageunha")
            and not os.environ.get("SLURM_JOB_ID")):
        raise RuntimeError("run full-image extended-PSF convolution on lageunha or through Slurm")
    manifest_path = data / "sample_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    frames = [row for row in manifest["clusters"] if str(row["id"]) == str(ident)]
    if len(frames) != 1:
        raise ValueError(f"expected one sample frame for cluster {ident}")
    frame = dict(frames[0])
    mask_records_path = icl_masks_dir / "records.tsv"
    mask_summary_path = icl_masks_dir / "summary.json"
    mask_summary = json.loads(mask_summary_path.read_text(encoding="utf-8"))
    mask_rows = [row for row in read_tsv(mask_records_path) if row["id"] == str(ident)]
    if len(mask_rows) != 1:
        raise ValueError("missing or duplicate saved ICL mask record")
    mask_record = mask_rows[0]
    frame["image"] = mask_record["image"]
    if str(mask_record["primary_filter"]).lower() != filter_name:
        raise ValueError(f"frame mask does not identify the requested {filter_name}-band coadd")
    pair_rows = [row for row in read_tsv(icl_pairs_path)
                 if row["id"] == str(ident) and row["filter"] == filter_name]
    if len(pair_rows) != 1 or float(pair_rows[0]["coverage_200_300"]) < 0.40:
        raise ValueError("frame is not in the paired K2.5/40% ICL support sample")

    prep = satellite_prep_dir / f"cluster_{ident}"
    detections_path, members_path = prep / "detections.tsv", prep / "members.tsv"
    sat_mask_path, sat_prep_path = prep / "masks.fits.gz", prep / "preparation.json"
    for path in (manifest_path, mask_records_path, mask_summary_path,
                 detections_path, members_path, sat_mask_path, sat_prep_path,
                 icl_pairs_path):
        if not path.is_file():
            raise FileNotFoundError(path)
    image_path = Path(frame["image"])
    if sha256(image_path) != mask_record["image_sha256"]:
        raise ValueError("image checksum differs from the saved ICL-mask input")
    icl_mask_path = Path(mask_record["mask_file"])
    if (sha256(icl_mask_path) != mask_record["mask_sha256"]
            or mask_summary.get("mask_output_sha256", {}).get(str(icl_mask_path))
            != mask_record["mask_sha256"]):
        raise ValueError("saved ICL mask checksum mismatch")

    info = satellite.validate_fits(image_path, frame)
    if str(info["coadd_header_filter"]).lower() != filter_name:
        raise ValueError(f"HSC image header does not identify a {filter_name}-band coadd")
    preparation = satellite.read_json(sat_prep_path)
    if preparation.get("science_measurement") is not False:
        raise ValueError("expected explicitly exploratory satellite preparation")
    detections_list = satellite.parse_detections(satellite.read_tsv(detections_path))
    member_rows = satellite.read_tsv(members_path)
    matched, match_stats = satellite.match_members(
        member_rows, detections_list, info["pixel_scale_arcsec"])
    center_xy = np.asarray(info["center_pixel_xy"], dtype=float)
    pixel_scale = float(info["pixel_scale_arcsec"])
    cosmology = FlatLambdaCDM(H0=float(manifest.get("h0", 70.)),
                              Om0=float(manifest.get("omega_m", .3)), Tcmb0=0*u.K)
    kpc_per_pixel = float(cosmology.kpc_proper_per_arcmin(float(frame["z"])).value
                           * pixel_scale / 60.0)
    conversion = float(satellite.nJy_factor(info))

    with fits.open(image_path, memmap=False, checksum=True) as hdus:
        original = np.asarray(hdus[info["science_hdu"]].data, dtype=np.float32)
        variance = np.asarray(hdus[info["variance_hdu"]].data, dtype=np.float32)
        bits = np.asarray(hdus[info["mask_hdu"]].data, dtype=np.int32)
    with fits.open(icl_mask_path, memmap=False, checksum=True) as hdus:
        icl_excluded = np.asarray(hdus["MASK_K2P5"].data, dtype=bool)
    with fits.open(sat_mask_path, memmap=False, checksum=True) as hdus:
        segmentation = np.asarray(hdus["SEGMENTATION"].data, dtype=np.int32)
        base = np.asarray(hdus["BASE_MASK"].data, dtype=bool)
    if not (original.shape == variance.shape == bits.shape == icl_excluded.shape
            == segmentation.shape == base.shape):
        raise ValueError("image and saved masks have different shapes")
    invalid = ~np.isfinite(original) | ~np.isfinite(variance) | (variance <= 0)
    bad_bits = sum(1 << bit for bit in satellite.BAD_BIT_POSITIONS) | satellite.BRIGHT_OBJECT_BIT
    if not np.array_equal(base, ((bits & bad_bits) != 0) | invalid):
        raise ValueError("satellite base mask no longer matches the source image")

    detections = satellite.detection_arrays(detections_list, pixel_scale)
    in_frame = [row for row in member_rows if row["IN_FRAME"] == "True"]
    member_xy = np.asarray([[float(row["X_IMAGE"])-1., float(row["Y_IMAGE"])-1.]
                            for row in in_frame], dtype=float).reshape((-1, 2))
    matched_indices = {item["member_index"] for item in matched}
    fallback_xy = np.asarray([[float(row["X_IMAGE"])-1., float(row["Y_IMAGE"])-1.]
                              for row in in_frame if int(row["INDEX"]) not in matched_indices],
                             dtype=float).reshape((-1, 2))
    from scipy.spatial import cKDTree
    fallback_tree = cKDTree(fallback_xy) if len(fallback_xy) else None
    member_tree = cKDTree(member_xy) if len(member_xy) else None
    targets, sat_centers = [], []
    for item in matched:
        det = detections_list[item["detection_index"]]
        radius = math.hypot(det["x"]-center_xy[0], det["y"]-center_xy[1]) * kpc_per_pixel
        if satellite.MIN_CLUSTER_RADIUS_KPC <= radius < satellite.MAX_CLUSTER_RADIUS_KPC:
            targets.append((item, det, radius))
            sat_centers.append([det["x"], det["y"]])
    if not targets:
        raise ValueError("no matched satellite targets in 50–500 kpc")
    expected_satellites = sum(
        any(role["cluster_id"] == str(ident) and role["frame_filter"] == filter_name
            and role["role"] == "satellite_injection_center" for role in row["roles"])
        for row in json.loads(psf_bank.mapping_path.read_text(encoding="utf-8")))
    if len(targets) != expected_satellites:
        raise ValueError("satellite target count differs from the approved PSF position plan")

    hybrid_metadata = None
    hybrid_core_hashes = {}
    hybrid_satellite_kernels = None
    hybrid_satellite_metadata = []
    if hybrid_parameters is not None:
        required = {"splice_radius_arcsec", "blend_width_arcsec", "wing_amplitude"}
        if set(hybrid_parameters) != required:
            raise ValueError(f"hybrid parameters must contain exactly {sorted(required)}")
        center_core, satellite_cores, psf_names = psf_bank.frame(
            ident, filter_name, len(targets))
        kernel, hybrid_metadata = build_core_wing_hybrid_psf(
            center_core, kernel, pixel_scale_arcsec=float(reference["pixel_scale_arcsec"]),
            **hybrid_parameters)
        if len(satellite_cores) != len(targets) or len(psf_names) != len(targets) + 1:
            raise ValueError("position PSF bank returned an incomplete hybrid mapping")
        hybrid_core_hashes = psf_bank.file_hashes(psf_names)
        hashes_by_name = {Path(path).name: digest
                          for path, digest in hybrid_core_hashes.items()}
        hybrid_satellite_kernels = []
        for index, (local_core, psf_name) in enumerate(
                zip(satellite_cores, psf_names[1:])):
            satellite_kernel, satellite_metadata = build_core_wing_hybrid_psf(
                local_core, kernel,
                pixel_scale_arcsec=float(reference["pixel_scale_arcsec"]),
                **hybrid_parameters)
            hybrid_satellite_kernels.append(satellite_kernel)
            hybrid_satellite_metadata.append({
                **satellite_metadata,
                "target_index": index,
                "local_core_filename": psf_name,
                "local_core_sha256": hashes_by_name[psf_name],
            })

    icl_geom = prestack_core.geometry(original.shape, center_xy, kpc_per_pixel)
    dx, dy, _, _, _ = icl_geom
    diffuse = prestack_core.diffuse_model(dx, dy, diffuse_scale_kpc, 1.0,
                                           diffuse_amplitude).astype(np.float32)
    diffuse_convolved = convolve_surface_brightness(diffuse, kernel)
    point_source_flux_njy = 1.0
    point_source = unit_point_source_surface_brightness(
        original.shape, center_xy, pixel_scale, point_source_flux_njy)
    point_source_convolved = convolve_surface_brightness(point_source, kernel)
    point_source_wing_profile = measure_point_source_wing_profile(
        point_source_convolved, icl_excluded, variance, icl_geom[2],
        kpc_per_pixel, pixel_scale, point_source_flux_njy)
    if hybrid_satellite_kernels is not None:
        satellite_convolved = satellite_halos_with_psfs(
            original.shape, np.asarray(sat_centers), kpc_per_pixel,
            satellite_amplitude, satellite_scale_kpc,
            hybrid_satellite_kernels)
    else:
        satellites = satellite_halos_truncated(
            original.shape, np.asarray(sat_centers), kpc_per_pixel,
            satellite_amplitude, satellite_scale_kpc)
        satellite_convolved = convolve_surface_brightness(satellites, kernel)
    empty = np.zeros_like(original, dtype=np.float32)
    scenarios = {"baseline": empty, "icl_only": diffuse_convolved,
                 "satellite_only": satellite_convolved,
                 "both": diffuse_convolved + satellite_convolved}
    measurements, fixed_mask_measurements, redetection = {}, {}, {}
    fixed_detected, fixed_targets = None, None
    for name, injected in scenarios.items():
        synthetic = add_surface_brightness_injection(original, injected, conversion)
        detected = redetect_masks_and_targets(
            synthetic, variance, bits, member_rows, center_xy, pixel_scale,
            unresolved_bcg_policy=unresolved_bcg_policy)
        targets_now = []
        for match, det, radius_pix in detected["targets_pixel_radius"]:
            radius_kpc = radius_pix * kpc_per_pixel
            if satellite.MIN_CLUSTER_RADIUS_KPC <= radius_kpc < satellite.MAX_CLUSTER_RADIUS_KPC:
                targets_now.append((match, det, radius_kpc))
        if not targets_now:
            raise ValueError("redetection left no matched satellite targets in 50–500 kpc")
        if name == "baseline":
            fixed_detected = detected
            fixed_targets = targets_now
        icl_result = prestack_core.measure(synthetic, variance, detected["icl_mask"],
                                           icl_geom, center_xy, kpc_per_pixel)
        sat_metric = satellite_outer_measurement(
            synthetic, variance, detected["base"], detected["segmentation"],
            detected["detections"], detected["fallback_tree"], detected["member_tree"],
            targets_now, center_xy, pixel_scale, kpc_per_pixel, conversion)
        measurements[name] = {
            "icl_raw_200_300_njy_arcsec2": float(icl_result["raw"][3] * conversion),
            "icl_plane_corrected_200_300_njy_arcsec2": float(icl_result["corrected"][3] * conversion),
            "satellite_delta_50_75_njy_arcsec2": sat_metric["delta_sb_njy_arcsec2"],
            "satellite_support": {key: sat_metric[key]
                                  for key in ("n_satellites", "n_control_profiles")},
            "icl_coverage_200_300": float(icl_result["coverage"][3]),
        }
        fixed_icl_result = prestack_core.measure(
            synthetic, variance, fixed_detected["icl_mask"], icl_geom,
            center_xy, kpc_per_pixel)
        fixed_sat_metric = satellite_outer_measurement(
            synthetic, variance, fixed_detected["base"], fixed_detected["segmentation"],
            fixed_detected["detections"], fixed_detected["fallback_tree"],
            fixed_detected["member_tree"], fixed_targets, center_xy,
            pixel_scale, kpc_per_pixel, conversion)
        fixed_mask_measurements[name] = {
            "icl_raw_200_300_njy_arcsec2": float(fixed_icl_result["raw"][3] * conversion),
            "icl_plane_corrected_200_300_njy_arcsec2": float(fixed_icl_result["corrected"][3] * conversion),
            "satellite_delta_50_75_njy_arcsec2": fixed_sat_metric["delta_sb_njy_arcsec2"],
            "satellite_support": {key: fixed_sat_metric[key]
                                  for key in ("n_satellites", "n_control_profiles")},
            "icl_coverage_200_300": float(fixed_icl_result["coverage"][3]),
        }
        redetection[name] = detected["diagnostics"]

    transfer_response = calculate_transfer_response(
        measurements, diffuse_amplitude, satellite_amplitude)
    fixed_mask_transfer_response = calculate_transfer_response(
        fixed_mask_measurements, diffuse_amplitude, satellite_amplitude)
    redetection_minus_fixed_mask_response = {
        key: transfer_response[key] - fixed_mask_transfer_response[key]
        for key in transfer_response
    }
    image_hashes = {str(path): sha256(path) for path in (
        manifest_path, mask_records_path, mask_summary_path, image_path, icl_mask_path,
        detections_path, members_path, sat_mask_path, sat_prep_path, icl_pairs_path)}
    result = {
        "schema": "astrafek-extended-psf-joint-injection-v3",
        "status": "complete",
        "science_measurement": False,
        "observational_validation_passed": False,
        "flux_covariance_available": False,
        "cluster_id": str(ident), "filter": filter_name, "redshift": float(frame["z"]),
        "injection": {
            "diffuse_profile": "circular exponential, normalized at R=100 kpc",
            "satellite_profile": "sum of circular exponentials, normalized at r=50 kpc; truncated at 15 scale lengths",
            "diffuse_amplitude_at_100_njy_arcsec2": diffuse_amplitude,
            "satellite_amplitude_at_50_per_galaxy_njy_arcsec2": satellite_amplitude,
            "diffuse_scale_kpc": diffuse_scale_kpc,
            "satellite_scale_kpc": satellite_scale_kpc,
            "separate_unit_point_source_flux_njy": point_source_flux_njy,
            "point_source_center_xy_pixel": center_xy.tolist(),
            "point_source_profile_uses_saved_icl_mask": True,
            "point_source_added_to_redetection_scenarios": False,
            "psf_convolution_applied": True,
            "psf_model": "published HSC PDR3 intermediate coadd/bg HSC-I extended reference; position/epoch/seeing unmatched",
            "frame_filter": filter_name,
            "HSC_I_to_i2_is_proxy_sensitivity_only": filter_name == "i2",
            "psf_reference_name": reference["name"],
            "psf_reference_sha256": reference["sha256"],
            "psf_finite_grid_sum": reference["finite_grid_sum"],
            "psf_negative_pixels_preserved": True,
            "source_detection_and_masks_rerun": True,
            "fixed_baseline_masks_also_measured": True,
            "reselection_effect_definition": "redetected response minus fixed-baseline-mask response",
            "unresolved_bcg_policy": unresolved_bcg_policy,
            "injected_source_photon_noise_added": False,
            "original_fits_modified": False,
            "n_injected_satellite_centers": len(targets),
        },
        "support": {"n_matched_members": len(matched),
                    "n_ambiguous_matches_rejected": match_stats["ambiguous"],
                    "n_injected_targets_in_50_500_kpc": len(targets)},
        "measurements_by_scenario": measurements,
        "measurements_by_scenario_fixed_masks": fixed_mask_measurements,
        "redetection_by_scenario": redetection,
        "point_source_wing_profile": point_source_wing_profile,
        "transfer_response": transfer_response,
        "fixed_mask_transfer_response": fixed_mask_transfer_response,
        "redetection_minus_fixed_mask_response": redetection_minus_fixed_mask_response,
        "input_sha256": image_hashes,
        "code_sha256": {
            str(Path(__file__).resolve()): sha256(Path(__file__)),
            str(Path(__file__).with_name("joint_icl_satellite_injection.py").resolve()):
                sha256(Path(__file__).with_name("joint_icl_satellite_injection.py")),
            str(Path(__file__).with_name("position_psf.py").resolve()):
                sha256(Path(__file__).with_name("position_psf.py")),
        },
        "limitations": [
            "The public HSC-I reference is not matched to this frame's position, epoch, or seeing.",
            "For i2 frames, HSC-I is used only as an exploratory proxy sensitivity; the resulting transfer is not an i2 PSF calibration.",
            "The reference symmetric/asymmetric variants share input stars and are not independent measurements.",
            "Only synthetic components are convolved; the observed coadd's original PSF and upstream sky are unchanged.",
            "No injected-source photon noise or correlated coadd noise is added.",
            "This is a selected-frame transfer sensitivity, not an ICL detection, completeness, or covariance measurement.",
        ],
    }
    result["input_sha256"].update({
        str(Path(reference["path"]).resolve()): reference["sha256"],
        str(Path(reference["summary_path"]).resolve()): reference["summary_sha256"],
        str(psf_bank.mapping_path): sha256(psf_bank.mapping_path),
    })
    if hybrid_metadata is not None:
        result["input_sha256"].update(hybrid_core_hashes)
        result["hybrid_psf"] = hybrid_metadata
        result["hybrid_satellite_psfs"] = hybrid_satellite_metadata
        result["injection"]["psf_model"] = (
            "position-specific center and satellite PSF cores smoothly joined to public HSC-I extended wing; exploratory sensitivity only"
        )
        result["injection"]["hybrid_parameters"] = dict(hybrid_parameters)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n")
    return result


def run_batch(data: Path, output_dir: Path, reference_path: Path,
              reference_dir: Path, plan_dir: Path, psf_dir: Path,
              test_center_psf: Path, *, only_cluster: str | None = None,
              filter_name: str = "i", hybrid_parameters: dict | None = None) -> dict:
    if filter_name not in ("i", "i2"):
        raise ValueError("filter_name must be i or i2")
    if output_dir.exists():
        raise FileExistsError(f"refusing to overwrite {output_dir}")
    bank = PositionPSFBank(plan_dir, psf_dir, test_center_psf)
    kernel, reference, reference_summary_path = load_public_reference(reference_path, reference_dir)
    reference["summary_path"] = str(reference_summary_path)
    pairs = select_jointly_supported_pairs(
        DEFAULT_ICL_PAIRS, DEFAULT_SATELLITE_TABLE)
    selected = select_reference_frames(pairs, bank.supported_frames, filter_name)
    if only_cluster is not None:
        selected = [pair for pair in selected if pair[0] == str(only_cluster)]
        if len(selected) != 1:
            raise ValueError(f"cluster {only_cluster} is not a unique selected {filter_name}-band frame")

    output_dir.mkdir(parents=True, exist_ok=False)
    started = time.time()
    records, errors = [], []
    for index, (ident, filt) in enumerate(selected, 1):
        output = output_dir / f"cluster_{ident}.json"
        print(f"[{index}/{len(selected)}] extended PSF injection cluster {ident} ({filt})", flush=True)
        try:
            result = run_frame(data, ident, output, reference, kernel,
                               DEFAULT_ICL_MASKS, DEFAULT_SAT_PREP, DEFAULT_ICL_PAIRS,
                               bank, filter_name=filter_name,
                               hybrid_parameters=hybrid_parameters)
            records.append({"cluster_id": ident, "filter": filt, "path": str(output),
                            "sha256": sha256(output),
                            "transfer_response": result["transfer_response"],
                            "fixed_mask_transfer_response": result["fixed_mask_transfer_response"],
                            "redetection_minus_fixed_mask_response": result["redetection_minus_fixed_mask_response"],
                            "hybrid_psf": result.get("hybrid_psf"),
                            "hybrid_satellite_psfs": result.get("hybrid_satellite_psfs"),
                            "support": result["support"]})
        except Exception as exc:
            errors.append({"cluster_id": ident, "filter": filt,
                           "error_type": type(exc).__name__, "error": str(exc)})
            print(f"  ERROR {type(exc).__name__}: {exc}", flush=True)
    summary = {
        "schema": "astrafek-extended-psf-joint-injection-batch-v3",
        "status": "complete" if not errors else "partial",
        "science_measurement": False,
        "observational_validation_passed": False,
        "flux_covariance_available": False,
        "n_selected": len(selected), "n_completed": len(records), "n_errors": len(errors),
        "frame_filter": filter_name,
        "public_reference_transfer_mode": ("exploratory_HSC_I_proxy_for_i2" if filter_name == "i2"
                                            else "HSC_I_reference_for_i_coadd"),
        "selected_pairs": [{"cluster_id": ident, "filter": filt} for ident, filt in selected],
        "only_cluster_pilot": only_cluster,
        "psf_model_mode": "position_core_public_wing_hybrid_sensitivity" if hybrid_parameters else "public_extended_reference_sensitivity",
        "hybrid_parameters": hybrid_parameters,
        "reference": {k: reference[k] for k in (
            "name", "path", "sha256", "summary_sha256", "license",
            "pixel_scale_arcsec", "finite_grid_sum", "negative_pixels",
            "psf_approved", "citation")},
        "records": records, "errors": errors,
        "reference_summary_path": str(reference_summary_path),
        "reference_summary_sha256": sha256(reference_summary_path),
        "position_mapping_sha256": sha256(plan_dir / "position_mapping.json"),
        "code_sha256": {
            str(Path(__file__).resolve()): sha256(Path(__file__)),
            str(Path(__file__).with_name("joint_icl_satellite_injection.py").resolve()):
                sha256(Path(__file__).with_name("joint_icl_satellite_injection.py")),
            str(Path(__file__).with_name("position_psf.py").resolve()):
                sha256(Path(__file__).with_name("position_psf.py")),
        },
        "elapsed_seconds": time.time() - started,
        "limitations": [
            f"Public extended reference is not position-, epoch-, or seeing-matched to the selected {filter_name}-band target frames.",
            "For i2, this is a proxy sensitivity only, not a filter-specific PSF calibration.",
            "The symmetric and asymmetric reference products share the same star sample and are sensitivity variants, not independent PSFs.",
            "Only synthetic light is convolved; observed data, source Poisson noise, coadd covariance, and upstream sky transfer are not changed.",
            "This is a model-dependence test on a support-selected subset, not science validation or a correction.",
        ],
    }
    summary_path = output_dir / "summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=DEFAULT_DATA)
    parser.add_argument("--reference-dir", type=Path, default=DEFAULT_REFERENCE_DIR)
    parser.add_argument("--reference-psf", type=Path, required=True)
    parser.add_argument("--plan-dir", type=Path, default=DEFAULT_PLAN_DIR)
    parser.add_argument("--psf-dir", type=Path, default=DEFAULT_PSFS)
    parser.add_argument("--test-center-psf", type=Path, default=DEFAULT_TEST_CENTER)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--cluster-id")
    parser.add_argument("--filter", choices=("i", "i2"), default="i",
                        help="analyze one coadd filter; i2 uses HSC-I only as a proxy sensitivity")
    parser.add_argument("--hybrid-splice-radius-arcsec", type=float)
    parser.add_argument("--hybrid-blend-width-arcsec", type=float)
    parser.add_argument("--hybrid-wing-amplitude", type=float)
    args = parser.parse_args()
    if (not platform.node().lower().startswith("lageunha")
            and not os.environ.get("SLURM_JOB_ID")):
        parser.error("run 4k-image extended-PSF convolution on lageunha or through Slurm")
    hybrid_values = (args.hybrid_splice_radius_arcsec,
                     args.hybrid_blend_width_arcsec,
                     args.hybrid_wing_amplitude)
    if any(value is not None for value in hybrid_values) and not all(
            value is not None for value in hybrid_values):
        parser.error("all three --hybrid-* options must be supplied together")
    hybrid_parameters = None
    if all(value is not None for value in hybrid_values):
        hybrid_parameters = {
            "splice_radius_arcsec": args.hybrid_splice_radius_arcsec,
            "blend_width_arcsec": args.hybrid_blend_width_arcsec,
            "wing_amplitude": args.hybrid_wing_amplitude,
        }
    result = run_batch(args.data.resolve(), args.output_dir.resolve(),
                       args.reference_psf.resolve(), args.reference_dir.resolve(),
                       args.plan_dir.resolve(), args.psf_dir.resolve(),
                       args.test_center_psf.resolve(), only_cluster=args.cluster_id,
                       filter_name=args.filter, hybrid_parameters=hybrid_parameters)
    print(json.dumps({k: result[k] for k in (
        "status", "n_selected", "n_completed", "n_errors", "elapsed_seconds",
        "only_cluster_pilot", "reference")}, indent=2))


if __name__ == "__main__":
    main()
