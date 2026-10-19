#!/usr/bin/env python3
"""Joint ICL/satellite injection test on one existing HSC frame.

The runner adds synthetic components only to an in-memory science array. It
reuses the stored K2.5 ICL mask and satellite segmentation by default. An
optional redetection mode reruns SEP detection and K2.5 masks for each injected
scenario. Neither mode is a science measurement or models HSC upstream sky,
PSF convolution, or injected-source photon noise.
"""
from __future__ import annotations

import argparse
import ast
import csv
import hashlib
import inspect
import json
import math
import os
import platform
import time
from pathlib import Path

import numpy as np
from astropy import units as u
from astropy.cosmology import FlatLambdaCDM
from astropy.io import fits
from scipy.spatial import cKDTree

import measure_full_sample_profiles as icl
import measure_satellite_outer_profiles as satellite
import prestack_core
import prepare_detection_psf as detection_pipeline
from position_psf import PositionPSFBank, convolve_surface_brightness, satellite_halos_with_psfs


DEFAULT_DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
DEFAULT_ICL_MASKS = DEFAULT_DATA / "prestack_cluster_masks_2026-10-16_01"
DEFAULT_SAT_PREP = DEFAULT_DATA / "satellite_detection_psf_2026-10-18_01"
DEFAULT_ICL_PAIRS = DEFAULT_DATA / "prestack_flux_profiles_2026-10-18_02/pair_metrics_k2p5_c40_r50_300.tsv"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


SCRIPT_SHA256_AT_IMPORT = sha256(Path(__file__))
POSITION_PSF_SCRIPT = Path(__file__).with_name("position_psf.py")
POSITION_PSF_SHA256_AT_IMPORT = sha256(POSITION_PSF_SCRIPT)


def runtime_code_hashes() -> dict[str, str]:
    """Hash project sources that implement this injection measurement."""
    paths = {
        Path(__file__).resolve(), POSITION_PSF_SCRIPT.resolve(),
        Path(satellite.__file__).resolve(), Path(prestack_core.__file__).resolve(),
        Path(detection_pipeline.__file__).resolve(),
        Path(satellite.validate_fits.__code__.co_filename).resolve(),
        Path(satellite.nJy_factor.__code__.co_filename).resolve(),
        Path(inspect.getsourcefile(detection_pipeline.Cosmology)).resolve(),
    }
    return {str(path): sha256(path) for path in sorted(paths)}


RUNTIME_CODE_SHA256_AT_IMPORT = runtime_code_hashes()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def satellite_halo_template(shape: tuple[int, int], centers_xy: np.ndarray,
                            kpc_per_pixel: float, amplitude_at_50: float,
                            scale_kpc: float = 15.0) -> np.ndarray:
    """Sum exponential satellite envelopes, normalized at 50 projected kpc."""
    if kpc_per_pixel <= 0 or scale_kpc <= 0 or amplitude_at_50 < 0:
        raise ValueError("scales must be positive and amplitude non-negative")
    yy, xx = np.ogrid[:shape[0], :shape[1]]
    result = np.zeros(shape, dtype=np.float64)
    for x, y in np.asarray(centers_xy, dtype=float).reshape((-1, 2)):
        radius = np.hypot((xx - x) * kpc_per_pixel,
                          (yy - y) * kpc_per_pixel)
        result += amplitude_at_50 * np.exp((50.0 - radius) / scale_kpc)
    return result.astype(np.float32)


def add_surface_brightness_injection(native_image: np.ndarray,
                                     injected_njy_arcsec2: np.ndarray,
                                     native_to_njy_arcsec2: float) -> np.ndarray:
    """Convert an injected SB model back to the image's native flux units."""
    if not math.isfinite(native_to_njy_arcsec2) or native_to_njy_arcsec2 <= 0:
        raise ValueError("native-to-nJy conversion must be finite and positive")
    if native_image.shape != injected_njy_arcsec2.shape:
        raise ValueError("native image and injection template shapes differ")
    return native_image + injected_njy_arcsec2 / native_to_njy_arcsec2


def satellite_outer_measurement(image: np.ndarray, variance: np.ndarray,
                                base: np.ndarray, segmentation: np.ndarray,
                                detections: dict, fallback_tree,
                                member_tree, targets: list[tuple[dict, dict, float]],
                                center_xy: np.ndarray, pixel_scale: float,
                                kpc_per_pixel: float, conversion: float,
                                coverage_floor: float = 0.60,
                                min_controls: int = 3) -> dict:
    """Match the saved per-annulus estimator for the 50--75 kpc bin."""
    target_values, control_values = [], []
    n_controlled = 0
    bin_index = len(satellite.R_EDGES_KPC) - 2
    for match, target_det, cluster_radius in targets:
        target_xy = np.asarray([target_det["x"], target_det["y"]], dtype=float)
        profile = satellite.profile_at(
            image, variance, base, segmentation, detections, fallback_tree,
            target_xy, pixel_scale, kpc_per_pixel, conversion,
            target_det["NUMBER"], coverage_floor)
        if profile is None or profile["coverage"][bin_index] < coverage_floor:
            continue
        points = satellite.control_centers(
            target_xy, center_xy, cluster_radius / kpc_per_pixel, detections,
            member_tree, pixel_scale, kpc_per_pixel, image.shape)
        controls = []
        for point in points:
            item = satellite.profile_at(
                image, variance, base, segmentation, detections, fallback_tree,
                point, pixel_scale, kpc_per_pixel, conversion, None,
                coverage_floor)
            if item is not None and item["coverage"][bin_index] >= coverage_floor:
                controls.append(float(item["sb"][bin_index]))
        if len(controls) < min_controls:
            continue
        target_values.append(float(profile["sb"][bin_index]))
        control_values.append(float(np.mean(controls)))
        n_controlled += len(controls)
    if not target_values:
        return {
            "delta_sb_njy_arcsec2": None,
            "n_satellites": 0,
            "n_control_profiles": 0,
            "status": "no_satellite_passed_minimum_support",
        }
    return {
        "delta_sb_njy_arcsec2": float(np.mean(np.asarray(target_values)
                                                - np.asarray(control_values))),
        "n_satellites": len(target_values),
        "n_control_profiles": n_controlled,
        "status": "supported",
    }


def redetect_masks_and_targets(image: np.ndarray, variance: np.ndarray,
                               bits: np.ndarray, member_rows: list[dict],
                               center_xy: np.ndarray, pixel_scale: float,
                               unresolved_bcg_policy: str = "error"):
    """Re-run SEP detection and regenerate K2.5 masks on an injected image."""
    if unresolved_bcg_policy not in {"error", "mask_all"}:
        raise ValueError("unresolved_bcg_policy must be 'error' or 'mask_all'")
    invalid = ~np.isfinite(image) | ~np.isfinite(variance) | (variance <= 0)
    bad_bits = sum(1 << bit for bit in satellite.BAD_BIT_POSITIONS)
    bad = ((bits & bad_bits) != 0) | invalid
    bright = (bits & satellite.BRIGHT_OBJECT_BIT) != 0
    base = bad | bright
    in_frame = [row for row in member_rows if row["IN_FRAME"] == "True"]
    member_xy = np.asarray([[float(row["X_IMAGE"]) - 1.,
                             float(row["Y_IMAGE"]) - 1.] for row in in_frame], dtype=float).reshape((-1, 2))
    objects, segmentation, kron, kron_flags, bcg, _, diagnostics = detection_pipeline.detect(
        image, variance, bad, center_xy, member_xy, pixel_scale,
    )
    bcg_candidates = int(diagnostics["bcg_matches"])
    bcg_for_mask = bcg
    if bcg < 0 and unresolved_bcg_policy == "error":
        raise ValueError(
            f"redetection found {bcg_candidates} BCG candidates at the supplied center "
            f"(nearest={diagnostics['nearest_bcg_detection_arcsec']})")
    if bcg < 0:
        # Conservatively mask every segmentation/source rather than guessing which
        # detection is the BCG. This is a separately recorded sensitivity mode.
        bcg_for_mask = -1
    masks, fallback_count = detection_pipeline.source_masks(
        objects, segmentation, kron, kron_flags, base, bcg_for_mask, member_xy, pixel_scale,
    )
    rows = []
    for index, obj in enumerate(objects):
        rows.append({
            "NUMBER": str(index + 1),
            "X_IMAGE": str(float(obj["x"]) + 1.), "Y_IMAGE": str(float(obj["y"]) + 1.),
            "A_IMAGE": str(float(obj["a"])), "B_IMAGE": str(float(obj["b"])),
            "THETA_RAD": str(float(obj["theta"])), "KRON_RADIUS": str(float(kron[index])),
            "KRON_FLAGS": str(int(kron_flags[index])),
        })
    detections = satellite.parse_detections(rows)
    matched, match_stats = satellite.match_members(member_rows, detections, pixel_scale)
    matched_indices = {item["member_index"] for item in matched}
    fallback_xy = np.asarray([[float(row["X_IMAGE"]) - 1., float(row["Y_IMAGE"]) - 1.]
                              for row in in_frame if int(row["INDEX"]) not in matched_indices],
                             dtype=float).reshape((-1, 2))
    fallback_tree = cKDTree(fallback_xy) if len(fallback_xy) else None
    member_tree = cKDTree(member_xy) if len(member_xy) else None
    det_arrays = satellite.detection_arrays(detections, pixel_scale)
    targets = []
    for item in matched:
        det = detections[item["detection_index"]]
        radius = math.hypot(det["x"] - center_xy[0], det["y"] - center_xy[1])
        # The caller converts this pixel radius to proper kpc using its frame cosmology.
        targets.append((item, det, radius))
    return {
        "base": base, "segmentation": segmentation, "detections": det_arrays,
        "fallback_tree": fallback_tree, "member_tree": member_tree,
        "targets_pixel_radius": targets, "icl_mask": masks[0],
        "diagnostics": {**diagnostics, "n_detections": len(rows),
                         "n_bcg_candidates": bcg_candidates,
                         "bcg_mask_policy": ("preserve_unique_bcg" if bcg >= 0 else "mask_all_sources"),
                         "n_fallback_member_masks": fallback_count,
                         "n_matched_members": len(matched),
                         "n_ambiguous_member_matches": match_stats["ambiguous"]},
    }


def run_frame(data: Path, ident: str, output: Path,
              icl_masks_dir: Path = DEFAULT_ICL_MASKS,
              satellite_prep_dir: Path = DEFAULT_SAT_PREP,
              icl_pairs_path: Path = DEFAULT_ICL_PAIRS,
              diffuse_amplitude: float = 1.0,
              satellite_amplitude: float = 1.0,
              diffuse_scale_kpc: float = 100.0,
              satellite_scale_kpc: float = 15.0,
              rerun_detection_masks: bool = False,
              unresolved_bcg_policy: str = "error",
              psf_bank: PositionPSFBank | None = None) -> dict:
    if output.exists():
        raise FileExistsError(f"refusing to overwrite {output}")
    if (not platform.node().lower().startswith("lageunha")
            and not os.environ.get("SLURM_JOB_ID")):
        raise RuntimeError("Run image-level injection on lageunha, not the syntax login node")
    if min(diffuse_amplitude, satellite_amplitude) < 0:
        raise ValueError("injected amplitudes must be non-negative")

    manifest_path = data / "sample_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    frames = [frame for frame in manifest["clusters"] if str(frame["id"]) == str(ident)]
    if len(frames) != 1:
        raise ValueError(f"expected exactly one frame for cluster {ident}")
    frame = dict(frames[0])
    mask_records_path = icl_masks_dir / "records.tsv"
    mask_summary_path = icl_masks_dir / "summary.json"
    mask_summary = json.loads(mask_summary_path.read_text(encoding="utf-8"))
    mask_records = [r for r in read_tsv(mask_records_path) if r["id"] == str(ident)]
    if len(mask_records) != 1:
        raise ValueError("missing or duplicate saved ICL mask record")
    mask_record = mask_records[0]
    frame["image"] = mask_record["image"]
    pair_rows = [r for r in read_tsv(icl_pairs_path) if r["id"] == str(ident)]
    pair_rows = [r for r in pair_rows if r["filter"] == mask_record["primary_filter"]]
    if len(pair_rows) != 1 or float(pair_rows[0]["coverage_200_300"]) < 0.40:
        raise ValueError("frame is not in the paired K2.5/40% ICL 200--300 kpc support sample")
    prep = satellite_prep_dir / f"cluster_{ident}"
    detections_path, members_path = prep / "detections.tsv", prep / "members.tsv"
    sat_mask_path, sat_prep_path = prep / "masks.fits.gz", prep / "preparation.json"
    for path in (manifest_path, mask_records_path, mask_summary_path, detections_path,
                 members_path, sat_mask_path, sat_prep_path, icl_pairs_path):
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
    if info["coadd_header_filter"] != mask_record["primary_filter"]:
        raise ValueError("filter mismatch between image and ICL mask record")
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
                           * pixel_scale / 60.)
    conversion = float(satellite.nJy_factor(info))

    with fits.open(image_path, memmap=False, checksum=True) as hdul:
        original = np.asarray(hdul[info["science_hdu"]].data, dtype=np.float32)
        variance = np.asarray(hdul[info["variance_hdu"]].data, dtype=np.float32)
        bits = np.asarray(hdul[info["mask_hdu"]].data, dtype=np.int32)
    with fits.open(icl_mask_path, memmap=False, checksum=True) as hdul:
        icl_excluded = np.asarray(hdul["MASK_K2P5"].data, dtype=bool)
    with fits.open(sat_mask_path, memmap=False, checksum=True) as hdul:
        segmentation = np.asarray(hdul["SEGMENTATION"].data, dtype=np.int32)
        base = np.asarray(hdul["BASE_MASK"].data, dtype=bool)
    if not (original.shape == variance.shape == bits.shape == icl_excluded.shape
            == segmentation.shape == base.shape):
        raise ValueError("image and saved masks have different shapes")
    invalid = ~np.isfinite(original) | ~np.isfinite(variance) | (variance <= 0)
    bad_bits = sum(1 << bit for bit in satellite.BAD_BIT_POSITIONS) | satellite.BRIGHT_OBJECT_BIT
    if not np.array_equal(base, ((bits & bad_bits) != 0) | invalid):
        raise ValueError("satellite base mask no longer matches the source image")
    if np.any(icl_excluded & ~base & ~np.isfinite(original)):
        raise ValueError("unexpected non-finite ICL-mask relation")

    detections = satellite.detection_arrays(detections_list, pixel_scale)
    in_frame_members = [r for r in member_rows if r["IN_FRAME"] == "True"]
    member_xy = np.asarray([[float(r["X_IMAGE"])-1., float(r["Y_IMAGE"])-1.]
                            for r in in_frame_members], dtype=float).reshape((-1, 2))
    matched_indices = {item["member_index"] for item in matched}
    fallback_xy = np.asarray([[float(r["X_IMAGE"])-1., float(r["Y_IMAGE"])-1.]
                              for r in in_frame_members
                              if int(r["INDEX"]) not in matched_indices], dtype=float).reshape((-1, 2))
    fallback_tree = cKDTree(fallback_xy) if len(fallback_xy) else None
    member_tree = cKDTree(member_xy) if len(member_xy) else None
    targets = []
    sat_centers = []
    for item in matched:
        det = detections_list[item["detection_index"]]
        radius = math.hypot(det["x"]-center_xy[0], det["y"]-center_xy[1]) * kpc_per_pixel
        if satellite.MIN_CLUSTER_RADIUS_KPC <= radius < satellite.MAX_CLUSTER_RADIUS_KPC:
            targets.append((item, det, radius))
            sat_centers.append([det["x"], det["y"]])
    if not targets:
        raise ValueError("no matched satellite targets in the configured cluster-radius range")

    icl_geom = prestack_core.geometry(original.shape, center_xy, kpc_per_pixel)
    dx, dy, _, _, _ = icl_geom
    diffuse_template = prestack_core.diffuse_model(
        dx, dy, diffuse_scale_kpc, 1.0, diffuse_amplitude).astype(np.float32)
    psf_names: list[str] = []
    if psf_bank is None:
        satellite_template = satellite_halo_template(
            original.shape, np.asarray(sat_centers), kpc_per_pixel,
            satellite_amplitude, satellite_scale_kpc)
    else:
        center_psf, satellite_psfs, psf_names = psf_bank.frame(
            ident, str(mask_record["primary_filter"]), len(sat_centers))
        diffuse_template = convolve_surface_brightness(diffuse_template, center_psf)
        satellite_template = satellite_halos_with_psfs(
            original.shape, np.asarray(sat_centers), kpc_per_pixel,
            satellite_amplitude, satellite_scale_kpc, satellite_psfs)
    empty = np.zeros_like(original, dtype=np.float32)
    scenarios = {
        "baseline": empty,
        "icl_only": diffuse_template,
        "satellite_only": satellite_template,
        "both": diffuse_template + satellite_template,
    }
    measurements = {}
    redetection_diagnostics = {}
    for name, injected in scenarios.items():
        synthetic = add_surface_brightness_injection(original, injected, conversion)
        if rerun_detection_masks:
            redetected = redetect_masks_and_targets(
                synthetic, variance, bits, member_rows, center_xy, pixel_scale,
                unresolved_bcg_policy,
            )
            targets_now = []
            for match, det, radius_pix in redetected["targets_pixel_radius"]:
                radius_kpc = radius_pix * kpc_per_pixel
                if satellite.MIN_CLUSTER_RADIUS_KPC <= radius_kpc < satellite.MAX_CLUSTER_RADIUS_KPC:
                    targets_now.append((match, det, radius_kpc))
            if not targets_now:
                raise ValueError("redetection left no matched satellite targets in the configured radius range")
            icl_result = prestack_core.measure(
                synthetic, variance, redetected["icl_mask"], icl_geom,
                center_xy, kpc_per_pixel)
            sat_metric = satellite_outer_measurement(
                synthetic, variance, redetected["base"], redetected["segmentation"],
                redetected["detections"], redetected["fallback_tree"],
                redetected["member_tree"], targets_now, center_xy, pixel_scale,
                kpc_per_pixel, conversion)
            redetection_diagnostics[name] = {
                **redetected["diagnostics"],
                "n_targets_in_50_500_kpc": len(targets_now),
                "icl_mask_pixels": int(np.count_nonzero(redetected["icl_mask"])),
            }
        else:
            icl_result = prestack_core.measure(
                synthetic, variance, icl_excluded, icl_geom, center_xy, kpc_per_pixel)
            sat_metric = satellite_outer_measurement(
                synthetic, variance, base, segmentation, detections, fallback_tree,
                member_tree, targets, center_xy, pixel_scale, kpc_per_pixel, conversion)
        icl_metric = float(icl_result["raw"][3] * conversion)
        measurements[name] = {"icl_raw_200_300_njy_arcsec2": icl_metric,
                              "icl_plane_corrected_200_300_njy_arcsec2": float(icl_result["corrected"][3] * conversion),
                              "satellite_delta_50_75_njy_arcsec2": sat_metric["delta_sb_njy_arcsec2"],
                              "satellite_support": {k: sat_metric[k] for k in
                                                    ("n_satellites", "n_control_profiles", "status")},
                              "icl_coverage_200_300": float(icl_result["coverage"][3])}

    base_icl = measurements["baseline"]["icl_raw_200_300_njy_arcsec2"]
    inc_icl = measurements["icl_only"]["icl_raw_200_300_njy_arcsec2"] - base_icl
    icl_from_sat = measurements["satellite_only"]["icl_raw_200_300_njy_arcsec2"] - base_icl
    both_icl = measurements["both"]["icl_raw_200_300_njy_arcsec2"] - base_icl

    def satellite_delta(scenario: str) -> float | None:
        value = measurements[scenario]["satellite_delta_50_75_njy_arcsec2"]
        return None if value is None else float(value)

    def satellite_response(scenario: str, amplitude: float) -> float | None:
        baseline = satellite_delta("baseline")
        value = satellite_delta(scenario)
        if baseline is None or value is None or amplitude == 0:
            return None
        return (value - baseline) / amplitude

    sat_from_icl = satellite_response("icl_only", diffuse_amplitude)
    inc_sat = satellite_response("satellite_only", satellite_amplitude)
    both_sat = satellite_response("both", satellite_amplitude)
    satellite_residual = (
        None if both_sat is None or sat_from_icl is None or inc_sat is None
        else (both_sat * satellite_amplitude
              - sat_from_icl * diffuse_amplitude
              - inc_sat * satellite_amplitude)
    )
    satellite_support_gaps = [
        name for name, row in measurements.items()
        if row["satellite_support"]["status"] != "supported"
    ]
    result = {
        "schema": ("hsc-joint-icl-satellite-redetection-injection-v1" if rerun_detection_masks
                   else "hsc-joint-icl-satellite-fixed-mask-injection-v1"),
        "status": "complete",
        "support_status": "partial_satellite_support" if satellite_support_gaps else "complete",
        "satellite_support_gaps": satellite_support_gaps,
        "science_measurement": False,
        "observational_validation_passed": False,
        "flux_covariance_available": False,
        "cluster_id": str(ident), "filter": info["coadd_header_filter"],
        "z": float(frame["z"]), "cluster_radius_selection_kpc": [50.0, 500.0],
        "satellite_profile_bin_kpc": [50.0, 75.0],
        "icl_profile_bin_kpc": [200.0, 300.0],
        "injection": {
            "diffuse_profile": "circular exponential, normalized to amplitude at R=100 kpc",
            "satellite_profile": "sum of circular exponentials, each normalized at r=50 kpc",
            "diffuse_amplitude_at_100_njy_arcsec2": diffuse_amplitude,
            "diffuse_scale_kpc": diffuse_scale_kpc,
            "satellite_amplitude_at_50_per_galaxy_njy_arcsec2": satellite_amplitude,
            "satellite_scale_kpc": satellite_scale_kpc,
            "template_unit": "nJy/arcsec^2, divided by FLUXMAG0-derived conversion before adding to native image",
            "n_injected_satellite_centers": len(targets),
            "source_detection_and_masks_rerun": rerun_detection_masks,
            "unresolved_bcg_policy": unresolved_bcg_policy if rerun_detection_masks else None,
            "psf_convolution_applied": psf_bank is not None,
            "psf_model": ("position-dependent HSC-I coadd model PSF core; center PSF for diffuse component, local PSF per satellite"
                          if psf_bank is not None else None),
            "psf_file_count": len(psf_names),
            "psf_stamp_pixels": ("position-dependent odd arrays, 41–43 pixels"
                                 if psf_bank is not None else None),
            "psf_negative_pixels_preserved": psf_bank is not None,
            "variance_and_poisson_noise_changed": False,
            "original_fits_modified": False,
        },
        "transfer_response": {
            "icl_response_to_diffuse_unit_injection": inc_icl / diffuse_amplitude if diffuse_amplitude else None,
            "icl_plane_corrected_response_to_diffuse_unit_injection": (
                (measurements["icl_only"]["icl_plane_corrected_200_300_njy_arcsec2"]
                 - measurements["baseline"]["icl_plane_corrected_200_300_njy_arcsec2"])
                / diffuse_amplitude if diffuse_amplitude else None),
            "icl_response_to_all_satellite_envelopes": icl_from_sat / satellite_amplitude if satellite_amplitude else None,
            "icl_plane_corrected_response_to_all_satellite_envelopes": (
                (measurements["satellite_only"]["icl_plane_corrected_200_300_njy_arcsec2"]
                 - measurements["baseline"]["icl_plane_corrected_200_300_njy_arcsec2"])
                / satellite_amplitude if satellite_amplitude else None),
            "satellite_response_to_smooth_diffuse": sat_from_icl,
            "satellite_response_to_satellite_unit_injection": inc_sat,
            "both_minus_sum_residual_icl_njy_arcsec2": both_icl - (inc_icl + icl_from_sat),
            "both_minus_sum_residual_satellite_njy_arcsec2": satellite_residual,
        },
        "measurements_by_scenario": measurements,
        "redetection_by_scenario": redetection_diagnostics,
        "support": {"n_matched_members": len(matched),
                    "n_ambiguous_matches_rejected": match_stats["ambiguous"],
                    "n_injected_targets_in_50_500_kpc": len(targets)},
        "input_sha256": {
            str(manifest_path): sha256(manifest_path), str(image_path): sha256(image_path),
            str(icl_mask_path): sha256(icl_mask_path), str(mask_records_path): sha256(mask_records_path),
            str(mask_summary_path): sha256(mask_summary_path), str(detections_path): sha256(detections_path),
            str(members_path): sha256(members_path), str(sat_mask_path): sha256(sat_mask_path),
            str(sat_prep_path): sha256(sat_prep_path), str(icl_pairs_path): sha256(icl_pairs_path),
        },
        "code_sha256": RUNTIME_CODE_SHA256_AT_IMPORT,
        "limitations": [
            ("This is a support-selected, redetection-and-mask-response sensitivity, not a population calibration."
             if rerun_detection_masks else "Single-frame fixed-mask measurement-transfer pilot, not a population calibration."),
            "Synthetic satellite envelopes are added to every positional member match; they are not fitted galaxy light profiles.",
            "The ICL metric is the target-side raw annular mean; target-minus-control stack injection is modeled by leaving controls unmodified.",
            (("SEP source detection and K2.5 masks are regenerated, but injected-source photon noise and HSC upstream sky processing are not modeled."
              if rerun_detection_masks else "Masks and source detections remain fixed, so detection completeness, segmentation/Kron response and mask growth are not tested.")
             if psf_bank is not None else
             ("SEP source detection and K2.5 masks are regenerated, but injected-source photon noise, PSF convolution and HSC upstream sky processing are not modeled."
              if rerun_detection_masks else "Masks and source detections remain fixed, so detection completeness, segmentation/Kron response, mask growth and PSF convolution are not tested.")),
            "The result is conditional on a single image, diagonal variance, assumed injection shapes, and existing automated centers/masks.",
            ("It uses only the compact model-PSF core; extended PSF wings, upstream sky processing, source Poisson noise, extinction and passband calibration remain unmodeled."
             if psf_bank is not None else "It does not include PSF convolution, upstream sky processing, source Poisson noise, extinction or passband calibration."),
        ],
        "elapsed_seconds": None,
    }
    if psf_bank is not None:
        result["input_sha256"].update({
            str(psf_bank.plan_path): sha256(psf_bank.plan_path),
            str(psf_bank.mapping_path): sha256(psf_bank.mapping_path),
            str(psf_bank.audit_path): sha256(psf_bank.audit_path),
            **psf_bank.file_hashes(psf_names),
        })
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n",
                      encoding="utf-8")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=DEFAULT_DATA)
    parser.add_argument("--icl-masks", type=Path, default=DEFAULT_ICL_MASKS)
    parser.add_argument("--satellite-prep", type=Path, default=DEFAULT_SAT_PREP)
    parser.add_argument("--icl-pairs", type=Path, default=DEFAULT_ICL_PAIRS)
    parser.add_argument("--cluster-id", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--diffuse-amplitude", type=float, default=1.0)
    parser.add_argument("--satellite-amplitude", type=float, default=1.0)
    parser.add_argument("--diffuse-scale-kpc", type=float, default=100.0)
    parser.add_argument("--satellite-scale-kpc", type=float, default=15.0)
    parser.add_argument("--rerun-detection-masks", action="store_true",
                        help="rerun SEP detection and regenerate source masks for each injection scenario")
    parser.add_argument("--unresolved-bcg-policy", choices=("error", "mask_all"), default="error",
                        help="on ambiguous/missing central detection, fail or mask every detected source")
    parser.add_argument("--psf-plan-dir", type=Path,
                        help="apply audited position-dependent HSC model PSF cores from this approved plan")
    parser.add_argument("--psf-download-dir", type=Path,
                        help="directory containing the audited PSF batch")
    parser.add_argument("--psf-test-center", type=Path,
                        default=Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_psf_picker_test_2026-10-19_01/cluster10_i2_center_psf.fits"))
    args = parser.parse_args()
    started = time.time()
    if bool(args.psf_plan_dir) != bool(args.psf_download_dir):
        parser.error("--psf-plan-dir and --psf-download-dir must be supplied together")
    psf_bank = (PositionPSFBank(args.psf_plan_dir.resolve(), args.psf_download_dir.resolve(),
                                args.psf_test_center.resolve())
                if args.psf_plan_dir else None)
    result = run_frame(args.data.resolve(), args.cluster_id, args.output.resolve(),
                       args.icl_masks.resolve(), args.satellite_prep.resolve(),
                       args.icl_pairs.resolve(),
                       args.diffuse_amplitude, args.satellite_amplitude,
                       args.diffuse_scale_kpc, args.satellite_scale_kpc,
                       args.rerun_detection_masks, args.unresolved_bcg_policy, psf_bank)
    result["elapsed_seconds"] = time.time() - started
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n",
                           encoding="utf-8")
    print(json.dumps({"cluster_id": result["cluster_id"], "filter": result["filter"],
                      "transfer_response": result["transfer_response"],
                      "support": result["support"], "output": str(args.output)}, indent=2))


if __name__ == "__main__":
    main()
