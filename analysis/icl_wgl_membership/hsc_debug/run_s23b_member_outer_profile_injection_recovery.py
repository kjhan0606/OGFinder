#!/usr/bin/env python3
"""Measure source-aware CAMIRA-member halo injection recovery on the full support.

One already-profiled 50--75 proper-kpc member is selected per eligible
target/filter frame. A circular exponential halo normalized at 50 kpc is
added in memory at amplitudes 0.5, 1, and 2 nJy arcsec^-2. SEP detection,
segmentation, member matching, source masks, and rotated controls are rebuilt
for every injected image. The run is a deterministic pipeline-transfer and
mask-selection diagnostic, not a science measurement; no data are downloaded.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import sys
import time
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from astropy.cosmology import FlatLambdaCDM
from astropy.io import fits
from scipy.spatial import cKDTree

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import measure_satellite_outer_profiles as profile_core
import measure_s23b_camira_member_outer_profiles as full_profile
import prepare_detection_psf as detection_pipeline
from position_psf import (S23BCamiraMemberProfilePSFBank,
                           satellite_halos_with_psfs)
from download_pilot import BAD_BIT_POSITIONS, BRIGHT_OBJECT_BIT, validate_fits

DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
TARGET_MANIFEST = DATA / "s23b_icl_target_image_plan_2026-11-07_01/target_manifest.json"
MEMBER_TABLE = DATA / "s23b_camira_member_catalog_2026-11-07_01/import_2026-11-07_02/target_member_positions.tsv"
MASK_RECORDS = DATA / "s23b_icl_mask_support_2026-11-07_01/target/records.tsv"
BASELINE_DIR = DATA / "s23b_camira_member_outer_profiles_full_2026-11-07_01"
BASELINE_SUMMARY = BASELINE_DIR / "summary.json"
BASELINE_FRAMES = BASELINE_DIR / "selected_target_frames.tsv"
BASELINE_PROFILES = BASELINE_DIR / "member_profile_bins.tsv"
DEFAULT_OUTPUT = DATA / "s23b_camira_member_outer_injection_recovery_2026-11-07_01"
DEFAULT_PSF_OUTPUT = DATA / "s23b_camira_member_outer_injection_recovery_psf_2026-11-07_01"

EXPECTED = {
    str(TARGET_MANIFEST): "e487582ba40cd3e58da67704242db4d7219e381fb08f4b33fed47ce4ff268b44",
    str(MEMBER_TABLE): "d1eb604922146c4d953ba4cad2f92d4c515fb4d4f10cd553f62b683535bfa6ae",
    str(MASK_RECORDS): "609354bfa7777088e995a487d0c2eaaed717b80e94589f1c5e0c3758f3e32ff8",
    str(BASELINE_SUMMARY): "e6b91fd7fe07caf89a6012053a46bba38be17fdb5298eb265b931ccc5e0c4cdd",
    str(BASELINE_FRAMES): "3fe2f2240a750be2ea597d7ceeeefef3593da2cd9cf9cb97f0cb193b1723457e",
    str(BASELINE_PROFILES): "56ff7feaaf2c2a7dd138c62f919df3d631d7798c9ff44cdc9482564045b963e1",
}
EXPECTED_BASELINE_RUNNER = "cc91779a92b3e6cfa86a564d9b9072f17f519292f95fca17c7186f1cbe3f4c87"
AMPLITUDES = (0.5, 1.0, 2.0)
SCALE_KPC = 15.0
PROFILE_BIN = "50-75"
PROFILE_BIN_INDEX = len(profile_core.R_EDGES_KPC) - 2
COVERAGE_FLOOR = 0.60
MIN_CONTROLS = 3
SEED = 20261107
N_BOOTSTRAP = 10_000
SCHEMA = "astrafek-s23b-member-outer-profile-injection-recovery-v1"


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


def write_json_atomic(path: Path, payload: dict) -> None:
    temporary = Path(str(path) + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True,
                                    allow_nan=False) + "\n", encoding="utf-8")
    temporary.replace(path)


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    temporary = Path(str(path) + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t",
                                extrasaction="raise", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def finite(value):
    if value is None:
        return None
    result = float(value)
    return result if math.isfinite(result) else None


def runtime_hashes() -> dict[str, str]:
    paths = {
        Path(__file__).resolve(), Path(full_profile.__file__).resolve(),
        Path(profile_core.__file__).resolve(), Path(detection_pipeline.__file__).resolve(),
        Path(__file__).with_name("position_psf.py").resolve(),
        Path(validate_fits.__code__.co_filename).resolve(),
        Path(profile_core.nJy_factor.__code__.co_filename).resolve(),
        Path(detection_pipeline.Cosmology.__module__.replace(".", "/") + ".py"),
    }
    paths = {p.resolve() for p in paths if p.is_file()}
    return {str(path): sha256(path) for path in sorted(paths)}


def verify_inputs() -> tuple[dict, list[dict], dict[str, list[dict]], dict[tuple[str, str], dict], dict]:
    for path, expected in EXPECTED.items():
        source = Path(path)
        if not source.is_file() or sha256(source) != expected:
            raise ValueError(f"pinned input is missing or changed: {source}")
    baseline = read_json(BASELINE_SUMMARY)
    if (baseline.get("status") != "complete_exploratory_source_aware_profile_stack"
            or baseline.get("science_measurement") is not False
            or baseline.get("analysis_script_sha256") != EXPECTED_BASELINE_RUNNER):
        raise ValueError("saved full-profile baseline is incomplete or has a different runner")

    target_manifest = read_json(TARGET_MANIFEST)
    frame_by_id = {str(row["id"]): row for row in target_manifest["clusters"]}
    selected_frames = read_tsv(BASELINE_FRAMES)
    record_rows = read_tsv(MASK_RECORDS)
    member_rows = read_tsv(MEMBER_TABLE)
    profile_rows = read_tsv(BASELINE_PROFILES)
    record_by_id = {str(row["id"]): row for row in record_rows}
    members_by_target: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in member_rows:
        members_by_target[str(row["target_cluster_id"])].append(row)

    frame_keys = {(str(row["target_id"]), str(row["filter"])) for row in selected_frames}
    if len(frame_keys) != len(selected_frames):
        raise ValueError("baseline selection contains duplicate target/filter frames")
    for ident, band in frame_keys:
        frame, record = frame_by_id.get(ident), record_by_id.get(ident)
        if frame is None or record is None or str(record["primary_filter"]).lower() != band:
            raise ValueError(f"baseline frame does not map to its manifest/mask record: {ident}/{band}")
        if str(record.get("image_sha256", "")) == "" or not Path(record["image"]).is_file():
            raise ValueError(f"missing audited image for {ident}")

    outer_rows = [row for row in profile_rows
                  if row["radial_bin_kpc"] == PROFILE_BIN]
    by_frame: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in outer_rows:
        by_frame[(str(row["target_id"]), str(row["filter"]))].append(row)
    if set(by_frame) - frame_keys:
        raise ValueError("outer-profile rows refer to a frame outside the pinned baseline")

    rng = np.random.default_rng(SEED)
    selected: list[dict] = []
    profile_lookup: dict[tuple[str, str], dict] = {}
    for key in sorted(by_frame):
        candidates = sorted(by_frame[key], key=lambda row: (int(row["member_index"]), row["source_row"]))
        if len({row["source_row"] for row in candidates}) != len(candidates):
            raise ValueError(f"duplicate 50--75 kpc member row for {key}")
        chosen = candidates[int(rng.integers(len(candidates)))]
        ident, band = key
        group = members_by_target.get(ident, [])
        index = int(chosen["member_index"])
        if index < 0 or index >= len(group) or group[index]["source_row"] != chosen["source_row"]:
            raise ValueError(f"baseline member index/source-row mapping differs for {key}")
        member = group[index]
        if member.get("inside_image_pixel_centers") != "True":
            raise ValueError(f"profiled member is not in image: {key}/{chosen['source_row']}")
        if int(chosen["n_valid_controls"]) < MIN_CONTROLS:
            raise ValueError(f"baseline profile lacks the required rotated controls: {key}")
        selected.append({
            "target_id": ident, "filter": band, "member_index": index,
            "source_row": chosen["source_row"],
            "cluster_radius_kpc": float(chosen["cluster_radius_kpc"]),
            "baseline_target_sb": float(chosen["target_sb_njy_arcsec2"]),
            "baseline_control_sb": float(chosen["control_sb_njy_arcsec2"]),
            "baseline_delta_sb": float(chosen["delta_sb_njy_arcsec2"]),
            "baseline_target_coverage": float(chosen["target_coverage"]),
            "baseline_control_coverage": float(chosen["control_coverage_mean"]),
            "baseline_n_controls": int(chosen["n_valid_controls"]),
            "baseline_detection_number": int(chosen["detection_number"]),
            "member_x": float(member["image_x_pixel"]),
            "member_y": float(member["image_y_pixel"]),
            "member_ra_deg": float(member["member_ra_deg"]),
            "member_dec_deg": float(member["member_dec_deg"]),
            "frame": frame_by_id[ident], "record": record_by_id[ident],
            "members": group,
        })
        profile_lookup[key] = chosen

    by_filter = {band: sum(item["filter"] == band for item in selected)
                 for band in ("i", "i2")}
    if by_filter != {"i": 62, "i2": 365}:
        raise ValueError(f"eligible member-profile frame counts changed: {by_filter}")
    return target_manifest, selected, members_by_target, profile_lookup, baseline


def make_context(image: np.ndarray, variance: np.ndarray, bits: np.ndarray,
                 frame: dict, member_group: list[dict[str, str]], info: dict,
                 cosmology: FlatLambdaCDM, selected: dict,
                 image_path: Path, image_sha: str) -> dict:
    bad_bitmask = sum(1 << bit for bit in BAD_BIT_POSITIONS)
    bad = (bits & bad_bitmask) != 0
    base = ((bits & (bad_bitmask | BRIGHT_OBJECT_BIT)) != 0)
    invalid = ~np.isfinite(image) | ~np.isfinite(variance) | (variance <= 0)
    base |= invalid

    member_rows = full_profile.build_member_rows(member_group)
    member_index = int(selected["member_index"])
    member_row = member_rows[member_index]
    member_xy = np.asarray([float(member_row["X_IMAGE"]) - 1.0,
                            float(member_row["Y_IMAGE"]) - 1.0], dtype=float)
    if not np.allclose(member_xy, [selected["member_x"], selected["member_y"]], rtol=0, atol=1e-9):
        raise ValueError(f"selected member coordinates changed for {selected['target_id']}")

    pixel_scale = float(info["pixel_scale_arcsec"])
    kpc_per_pixel = float(cosmology.kpc_proper_per_arcmin(float(frame["z"])).value
                           * pixel_scale / 60.0)
    conversion = float(profile_core.nJy_factor(info))
    if not math.isfinite(conversion) or conversion <= 0:
        raise ValueError(f"invalid native-to-SB conversion for {selected['target_id']}")
    center_xy = np.asarray(info["center_pixel_xy"], dtype=float)
    yy, xx = np.indices(image.shape, dtype=np.float32)
    radius_kpc = np.hypot((xx - member_xy[0]) * kpc_per_pixel,
                          (yy - member_xy[1]) * kpc_per_pixel)
    unit_profile = np.exp((50.0 - radius_kpc) / SCALE_KPC).astype(np.float32)
    return {
        "image": image, "variance": variance, "base": base, "bad": bad,
        "member_rows": member_rows, "member_xy": member_xy,
        "member_tree_points": np.asarray([
            [float(row["X_IMAGE"]) - 1.0, float(row["Y_IMAGE"]) - 1.0]
            for row in member_rows if row["IN_FRAME"] == "True"], dtype=float),
        "frame": frame, "info": info, "cosmology": cosmology,
        "selected": selected, "center_xy": center_xy,
        "pixel_scale": pixel_scale, "kpc_per_pixel": kpc_per_pixel,
        "conversion": conversion, "unit_profile_sb": unit_profile,
        "image_path": str(image_path), "image_sha256": image_sha,
    }


def detect_and_measure(context: dict, image: np.ndarray,
                       model_native: np.ndarray | None = None) -> dict:
    objects, segmentation, kron, kron_flags, bcg, snr, detection_audit = detection_pipeline.detect(
        image, context["variance"], context["bad"], context["center_xy"],
        context["member_tree_points"], context["pixel_scale"],
    )
    detections = full_profile.detection_rows(objects, kron, kron_flags)
    if not detections:
        return {"status": "no_sep_detections", "n_detections": 0,
                "segmentation_pixels": 0, "bcg_matches": detection_audit["bcg_matches"]}
    arrays = profile_core.detection_arrays(detections, context["pixel_scale"])
    member_rows = context["member_rows"]
    matched, match_stats = profile_core.match_members(member_rows, detections,
                                                       context["pixel_scale"])
    matched_by_index = {int(row["member_index"]): row for row in matched}
    selected = context["selected"]
    match = matched_by_index.get(int(selected["member_index"]))
    frame_record = {
        "n_detections": len(detections),
        "segmentation_pixels": int(np.count_nonzero(segmentation)),
        "n_member_matches": len(matched),
        "match_stats": match_stats,
        "bcg_matches": detection_audit["bcg_matches"],
    }
    if match is None:
        return {"status": "member_not_uniquely_recovered", **frame_record}

    in_frame_indices = [int(row["INDEX"]) for row in member_rows
                        if row["IN_FRAME"] == "True"]
    matched_indices = set(matched_by_index)
    fallback_xy = np.asarray([
        [float(member_rows[index]["X_IMAGE"]) - 1.0,
         float(member_rows[index]["Y_IMAGE"]) - 1.0]
        for index in in_frame_indices if index not in matched_indices
    ], dtype=float)
    fallback_tree = cKDTree(fallback_xy) if len(fallback_xy) else None
    member_xy = context["member_tree_points"]
    member_tree = cKDTree(member_xy) if len(member_xy) else None
    target_xy = np.asarray([match["member_x"], match["member_y"]], dtype=float)
    detection_number = int(match["detection_number"])
    target = profile_core.profile_at(
        image, context["variance"], context["base"], segmentation, arrays,
        fallback_tree, target_xy, context["pixel_scale"], context["kpc_per_pixel"],
        context["conversion"], detection_number, COVERAGE_FLOOR,
    )
    if target is None or target["coverage"][PROFILE_BIN_INDEX] < COVERAGE_FLOOR:
        return {"status": "target_annulus_below_coverage", **frame_record,
                "detection_number": detection_number,
                "target_coverage": None if target is None else finite(target["coverage"][PROFILE_BIN_INDEX]),
                "match_arcsec": finite(match["match_arcsec"])}

    cluster_radius_pix = float(selected["cluster_radius_kpc"]) / context["kpc_per_pixel"]
    control_xy = profile_core.control_centers(
        target_xy, context["center_xy"], cluster_radius_pix, arrays, member_tree,
        context["pixel_scale"], context["kpc_per_pixel"], image.shape,
    )
    controls = []
    model_controls = []
    for point in control_xy:
        profile = profile_core.profile_at(
            image, context["variance"], context["base"], segmentation, arrays,
            fallback_tree, point, context["pixel_scale"], context["kpc_per_pixel"],
            context["conversion"], None, COVERAGE_FLOOR,
        )
        if profile is not None and profile["coverage"][PROFILE_BIN_INDEX] >= COVERAGE_FLOOR:
            controls.append(profile)
        if model_native is not None:
            model_profile = profile_core.profile_at(
                model_native, context["variance"], context["base"], segmentation, arrays,
                fallback_tree, point, context["pixel_scale"], context["kpc_per_pixel"],
                context["conversion"], None, COVERAGE_FLOOR,
            )
            if (model_profile is not None
                    and model_profile["coverage"][PROFILE_BIN_INDEX] >= COVERAGE_FLOOR):
                model_controls.append(model_profile)
    if len(controls) < MIN_CONTROLS:
        return {"status": "too_few_rotated_controls", **frame_record,
                "detection_number": detection_number,
                "target_coverage": finite(target["coverage"][PROFILE_BIN_INDEX]),
                "n_valid_controls": len(controls),
                "match_arcsec": finite(match["match_arcsec"])}

    index = PROFILE_BIN_INDEX
    result = {
        "status": "supported",
        **frame_record,
        "detection_number": detection_number,
        "match_arcsec": finite(match["match_arcsec"]),
        "target_sb_njy_arcsec2": finite(target["sb"][index]),
        "control_sb_njy_arcsec2": finite(np.mean([row["sb"][index] for row in controls])),
        "delta_sb_njy_arcsec2": finite(
            target["sb"][index] - np.mean([row["sb"][index] for row in controls])),
        "target_coverage": finite(target["coverage"][index]),
        "control_coverage_mean": finite(np.mean([row["coverage"][index] for row in controls])),
        "n_valid_controls": len(controls),
    }
    if model_native is not None:
        model_target = profile_core.profile_at(
            model_native, context["variance"], context["base"], segmentation, arrays,
            fallback_tree, target_xy, context["pixel_scale"], context["kpc_per_pixel"],
            context["conversion"], detection_number, COVERAGE_FLOOR,
        )
        if model_target is None or model_target["coverage"][index] < COVERAGE_FLOOR:
            result["model_status"] = "target_model_below_coverage"
        elif len(model_controls) < MIN_CONTROLS:
            result["model_status"] = "too_few_model_controls"
        else:
            target_model_sb = float(model_target["sb"][index])
            control_model_sb = float(np.mean([row["sb"][index] for row in model_controls]))
            result.update({
                "model_status": "supported",
                "model_target_sb_njy_arcsec2": finite(target_model_sb),
                "model_control_sb_njy_arcsec2": finite(control_model_sb),
                "model_delta_sb_njy_arcsec2": finite(target_model_sb - control_model_sb),
                "model_target_coverage": finite(model_target["coverage"][index]),
                "model_n_valid_controls": len(model_controls),
            })
    return result


def nominal_annulus_mean(unit_profile: np.ndarray, member_xy: np.ndarray,
                         kpc_per_pixel: float, amplitude: float) -> float:
    half_size = int(math.ceil(profile_core.R_EDGES_KPC[-1] / kpc_per_pixel)) + 2
    cx, cy = map(float, member_xy)
    x0, x1 = max(0, int(math.floor(cx - half_size))), min(unit_profile.shape[1], int(math.ceil(cx + half_size)) + 1)
    y0, y1 = max(0, int(math.floor(cy - half_size))), min(unit_profile.shape[0], int(math.ceil(cy + half_size)) + 1)
    yy, xx = np.indices((y1 - y0, x1 - x0), dtype=np.float32)
    radius = np.hypot(xx + x0 - cx, yy + y0 - cy) * kpc_per_pixel
    annulus = (radius >= profile_core.R_EDGES_KPC[-2]) & (radius < profile_core.R_EDGES_KPC[-1])
    if not np.any(annulus):
        return 0.0
    return float(np.mean(unit_profile[y0:y1, x0:x1][annulus], dtype=np.float64) * amplitude)


def process_frame(task: dict) -> dict:
    selected = task["selected"]
    frame = dict(selected["frame"])
    record = selected["record"]
    image_path = Path(record["image"])
    if image_path.stat().st_size != int(record["image_bytes"]):
        raise ValueError(f"image byte count changed: {image_path}")
    image_sha = sha256(image_path)
    if image_sha != record["image_sha256"]:
        raise ValueError(f"image checksum differs from audited record: {image_path}")
    info = validate_fits(image_path, frame)
    if info["coadd_header_filter"] != selected["filter"]:
        raise ValueError(f"primary filter mismatch for {selected['target_id']}")
    with fits.open(image_path, mode="readonly", memmap=False, checksum=True) as hdus:
        image = np.asarray(hdus[info["science_hdu"]].data, dtype=np.float32)
        variance = np.asarray(hdus[info["variance_hdu"]].data, dtype=np.float32)
        bits = np.asarray(hdus[info["mask_hdu"]].data, dtype=np.int32)
    if image.shape != variance.shape or image.shape != bits.shape:
        raise ValueError(f"FITS image/VARIANCE/MASK dimensions differ: {selected['target_id']}")

    manifest = task["manifest"]
    cosmology = FlatLambdaCDM(H0=float(manifest["h0"]), Om0=float(manifest["omega_m"]))
    context = make_context(image, variance, bits, frame, selected["members"],
                           info, cosmology, selected, image_path, image_sha)
    baseline = detect_and_measure(context, image)
    if baseline["status"] != "supported":
        raise ValueError(f"baseline member profile unsupported for {selected['target_id']}: {baseline['status']}")

    baseline_fields = {
        "target_sb_njy_arcsec2": selected["baseline_target_sb"],
        "control_sb_njy_arcsec2": selected["baseline_control_sb"],
        "delta_sb_njy_arcsec2": selected["baseline_delta_sb"],
        "target_coverage": selected["baseline_target_coverage"],
        "control_coverage_mean": selected["baseline_control_coverage"],
    }
    baseline_replay_differences = {
        key: float(baseline[key]) - float(expected)
        for key, expected in baseline_fields.items()
    }
    if any(abs(value) > 2e-5 for value in baseline_replay_differences.values()):
        raise ValueError(
            f"baseline replay differs from saved full profile for {selected['target_id']}: "
            f"{baseline_replay_differences}")
    if int(baseline["detection_number"]) != int(selected["baseline_detection_number"]):
        raise ValueError(f"baseline matched detection number changed for {selected['target_id']}")

    member_psf = task.get("member_psf")
    if member_psf is None:
        unit_model_sb = context["unit_profile_sb"]
    else:
        unit_model_sb = satellite_halos_with_psfs(
            image.shape, context["member_xy"][None, :], context["kpc_per_pixel"],
            1.0, SCALE_KPC, [member_psf], truncation_scale_lengths=15.0,
        )
    nominal_unit = nominal_annulus_mean(unit_model_sb, context["member_xy"],
                                        context["kpc_per_pixel"], 1.0)
    rows = []
    for amplitude in AMPLITUDES:
        model_sb = unit_model_sb * np.float32(amplitude)
        model_native = np.asarray(model_sb / np.float32(context["conversion"]), dtype=np.float32)
        injected_image = np.asarray(image + model_native, dtype=np.float32)
        measured = detect_and_measure(context, injected_image, model_native)
        nominal_input = float(nominal_unit * amplitude)
        row = {
            "target_id": selected["target_id"], "filter": selected["filter"],
            "member_index": int(selected["member_index"]), "source_row": selected["source_row"],
            "cluster_radius_kpc": float(selected["cluster_radius_kpc"]),
            "member_ra_deg": selected["member_ra_deg"],
            "member_dec_deg": selected["member_dec_deg"],
            "amplitude_at_50_kpc_njy_arcsec2": float(amplitude),
            "exponential_scale_kpc": SCALE_KPC,
            "nominal_unmasked_50_75_mean_njy_arcsec2": nominal_input,
            "psf_convolution_applied": member_psf is not None,
            "member_psf_filename": task.get("member_psf_filename"),
            "member_psf_sha256": task.get("member_psf_sha256"),
            "member_psf_shape": (list(member_psf.shape) if member_psf is not None else None),
            "member_psf_negative_pixels_preserved": (
                bool(np.any(member_psf < 0)) if member_psf is not None else None),
            "baseline_delta_sb_njy_arcsec2": float(baseline["delta_sb_njy_arcsec2"]),
            "baseline_target_sb_njy_arcsec2": float(baseline["target_sb_njy_arcsec2"]),
            "baseline_control_sb_njy_arcsec2": float(baseline["control_sb_njy_arcsec2"]),
            "baseline_target_coverage": float(baseline["target_coverage"]),
            "baseline_control_coverage_mean": float(baseline["control_coverage_mean"]),
            "baseline_n_controls": int(baseline["n_valid_controls"]),
            "baseline_n_detections": int(baseline["n_detections"]),
            "baseline_segmentation_pixels": int(baseline["segmentation_pixels"]),
            "baseline_replay_delta_abs_njy_arcsec2": abs(baseline_replay_differences["delta_sb_njy_arcsec2"]),
            "injected_status": measured["status"],
            "injected_n_detections": measured.get("n_detections"),
            "injected_segmentation_pixels": measured.get("segmentation_pixels"),
            "injected_detection_count_change": (
                int(measured["n_detections"]) - int(baseline["n_detections"])
                if measured.get("n_detections") is not None else None),
            "injected_segmentation_pixel_change": (
                int(measured["segmentation_pixels"]) - int(baseline["segmentation_pixels"])
                if measured.get("segmentation_pixels") is not None else None),
            "science_measurement": False,
        }
        if measured["status"] == "supported":
            injected_delta = float(measured["delta_sb_njy_arcsec2"])
            recovered = injected_delta - float(baseline["delta_sb_njy_arcsec2"])
            row.update({
                "injected_target_sb_njy_arcsec2": measured["target_sb_njy_arcsec2"],
                "injected_control_sb_njy_arcsec2": measured["control_sb_njy_arcsec2"],
                "injected_delta_sb_njy_arcsec2": injected_delta,
                "recovered_delta_increment_njy_arcsec2": recovered,
                "recovered_increment_over_nominal_input": (
                    recovered / nominal_input if nominal_input > 0 else None),
                "injected_target_coverage": measured["target_coverage"],
                "injected_control_coverage_mean": measured["control_coverage_mean"],
                "injected_n_controls": measured["n_valid_controls"],
                "injected_member_match_arcsec": measured["match_arcsec"],
                "injected_detection_number": measured["detection_number"],
                "model_status": measured.get("model_status", "not_available"),
                "masked_model_delta_njy_arcsec2": measured.get("model_delta_sb_njy_arcsec2"),
                "masked_model_over_nominal_input": (
                    measured["model_delta_sb_njy_arcsec2"] / nominal_input
                    if measured.get("model_delta_sb_njy_arcsec2") is not None and nominal_input > 0
                    else None),
                "recovered_minus_masked_model_njy_arcsec2": (
                    recovered - measured["model_delta_sb_njy_arcsec2"]
                    if measured.get("model_delta_sb_njy_arcsec2") is not None else None),
            })
        else:
            row.update({key: None for key in (
                "injected_target_sb_njy_arcsec2", "injected_control_sb_njy_arcsec2",
                "injected_delta_sb_njy_arcsec2", "recovered_delta_increment_njy_arcsec2",
                "recovered_increment_over_nominal_input", "injected_target_coverage",
                "injected_control_coverage_mean", "injected_n_controls",
                "injected_member_match_arcsec", "injected_detection_number",
                "model_status", "masked_model_delta_njy_arcsec2",
                "masked_model_over_nominal_input", "recovered_minus_masked_model_njy_arcsec2")})
        rows.append(row)

    return {
        "status": "complete",
        "target_id": selected["target_id"], "filter": selected["filter"],
        "source_row": selected["source_row"],
        "image_sha256": image_sha,
        "psf_convolution_applied": member_psf is not None,
        "member_psf_filename": task.get("member_psf_filename"),
        "member_psf_sha256": task.get("member_psf_sha256"),
        "baseline_replay_differences": baseline_replay_differences,
        "baseline_n_detections": int(baseline["n_detections"]),
        "baseline_segmentation_pixels": int(baseline["segmentation_pixels"]),
        "injections": rows,
    }


def bootstrap_summary(rows: list[dict]) -> dict:
    usable = [row for row in rows if row.get("injected_status") == "supported"
              and row.get("recovered_increment_over_nominal_input") is not None]
    by_filter: dict[str, list[dict]] = {
        band: [row for row in usable if row["filter"] == band]
        for band in ("i", "i2")
    }
    rng = np.random.default_rng(SEED + 17)
    result = {}
    for band, band_rows in by_filter.items():
        frame_keys = sorted({(row["target_id"], row["filter"], row["source_row"])
                             for row in band_rows})
        indexed = {amp: [] for amp in AMPLITUDES}
        for key in frame_keys:
            frame_rows = {float(row["amplitude_at_50_kpc_njy_arcsec2"]): row
                          for row in band_rows
                          if (row["target_id"], row["filter"], row["source_row"]) == key}
            for amp in AMPLITUDES:
                if amp in frame_rows:
                    indexed[amp].append(frame_rows[amp])
        band_summary = {}
        for amp in AMPLITUDES:
            amp_rows = indexed[amp]
            response = np.asarray([float(row["recovered_increment_over_nominal_input"])
                                   for row in amp_rows], dtype=float)
            transfer = np.asarray([
                float(row["masked_model_over_nominal_input"])
                for row in amp_rows if row.get("masked_model_over_nominal_input") is not None
            ], dtype=float)
            if response.size:
                boot_median = np.empty(N_BOOTSTRAP, dtype=float)
                for draw in range(N_BOOTSTRAP):
                    sample = rng.integers(0, response.size, response.size)
                    boot_median[draw] = float(np.median(response[sample]))
                band_summary[str(amp)] = {
                    "n_supported_frames": int(response.size),
                    "median_recovered_response": float(np.median(response)),
                    "response_q16_q84": [float(x) for x in np.quantile(response, [0.16, 0.84])],
                    "cluster_bootstrap_median_ci95": [float(x) for x in np.quantile(boot_median, [0.025, 0.975])],
                    "median_masked_model_transfer": float(np.median(transfer)) if transfer.size else None,
                    "n_masked_model_rows": int(transfer.size),
                    "fraction_detection_count_changed": float(np.mean([
                        int(row["injected_detection_count_change"] or 0) != 0
                        for row in amp_rows])),
                    "median_detection_count_change": float(np.median([
                        int(row["injected_detection_count_change"] or 0)
                        for row in amp_rows])),
                    "median_segmentation_pixel_change": float(np.median([
                        int(row["injected_segmentation_pixel_change"] or 0)
                        for row in amp_rows])),
                }
            else:
                band_summary[str(amp)] = {"n_supported_frames": 0}

        ref_rows = { (row["target_id"], row["source_row"]): row
                     for row in indexed[1.0] }
        for amp in (0.5, 2.0):
            diffs = []
            for row in indexed[amp]:
                key = (row["target_id"], row["source_row"])
                ref = ref_rows.get(key)
                if ref is not None:
                    diffs.append(float(row["recovered_increment_over_nominal_input"])
                                 - float(ref["recovered_increment_over_nominal_input"]))
            if diffs:
                values = np.asarray(diffs, dtype=float)
                draws = np.empty(N_BOOTSTRAP, dtype=float)
                for index in range(N_BOOTSTRAP):
                    draws[index] = float(np.median(values[rng.integers(0, len(values), len(values))]))
                band_summary[f"paired_{amp:g}_minus_1"] = {
                    "n_paired_frames": len(values),
                    "median_response_difference": float(np.median(values)),
                    "cluster_bootstrap_median_ci95": [float(x) for x in np.quantile(draws, [0.025, 0.975])],
                }
        result[band] = band_summary
    return result


def render_plot(path: Path, summary: dict, *, psf_convolved: bool = False) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.8), sharey=True)
    colors = {"i": "#3568a8", "i2": "#b65b30"}
    for ax, band, label in zip(axes, ("i", "i2"), ("HSC-i", "HSC-i2")):
        stats = summary.get(band, {})
        xs, medians, lows, highs = [], [], [], []
        for amp in AMPLITUDES:
            row = stats.get(str(amp), {})
            if row.get("n_supported_frames", 0) == 0:
                continue
            xs.append(amp)
            median = float(row["median_recovered_response"])
            lo, hi = row["cluster_bootstrap_median_ci95"]
            medians.append(median)
            lows.append(median - lo)
            highs.append(hi - median)
        if xs:
            ax.errorbar(xs, medians, yerr=np.asarray([lows, highs]), fmt="o-",
                        capsize=3, color=colors[band], label=f"{label}: median, cluster-bootstrap 95% CI")
        ax.axhline(1.0, color="0.35", lw=0.9, ls="--", label="unit recovery")
        ax.set_xscale("log", base=2)
        ax.set_xticks(AMPLITUDES, [str(x) for x in AMPLITUDES])
        ax.set_xlabel("Injected SB at 50 kpc [nJy arcsec⁻²]")
        ax.set_title(label)
        ax.legend(frameon=False, fontsize=9, loc="best")
    axes[0].set_ylabel("Recovered target-minus-rotated-control / unmasked input")
    title = ("Member-halo injection transfer; model-PSF core, no photon noise"
             if psf_convolved else
             "Member-halo profile injection transfer; no PSF or photon-noise injection")
    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(path, dpi=160)
    plt.close(fig)


def manifest_payload(input_hashes: dict[str, str]) -> dict:
    return {
        "schema": SCHEMA,
        "inputs_sha256": input_hashes,
        "runtime_code_sha256": runtime_hashes(),
        "selection": "one seeded-random CAMIRA member per target/filter frame with an existing supported 50--75 kpc profile",
        "n_frames_by_filter": {"i": 62, "i2": 365},
        "amplitude_at_50_kpc_njy_arcsec2": list(AMPLITUDES),
        "satellite_profile": "circular exponential exp((50 kpc-r)/15 kpc), no PSF convolution",
        "measurement": "baseline reproduces saved full profile; inject in target image; rerun SEP segmentation, one-to-one CAMIRA match, own-source-preserving masks, and same-radius rotated controls",
        "coverage_floor": COVERAGE_FLOOR,
        "minimum_rotated_controls": MIN_CONTROLS,
        "redetection_repeated_for_each_amplitude": True,
        "source_poisson_added": False,
        "variance_updated": False,
        "science_measurement": False,
        "images_ready_for_science": False,
        "network_used": False,
        "data_downloaded": False,
        "seed": SEED,
        "bootstrap_replicas": N_BOOTSTRAP,
    }


def run(output: Path, *, workers: int, resume: bool,
        psf_bank: S23BCamiraMemberProfilePSFBank | None = None) -> dict:
    output = output.resolve()
    target_manifest, selected, members_by_target, _, baseline_summary = verify_inputs()
    input_hashes = {path: sha256(Path(path)) for path in EXPECTED}
    run_manifest = manifest_payload(input_hashes)
    run_manifest["psf_input"] = psf_bank.manifest() if psf_bank is not None else None
    run_manifest["satellite_profile"] = (
        "circular exponential exp((50 kpc-r)/15 kpc), convolved with the member-position HSC-I coadd model PSF core"
        if psf_bank is not None else
        "circular exponential exp((50 kpc-r)/15 kpc), no PSF convolution")
    run_manifest["source_poisson_added"] = False
    run_manifest["script_sha256"] = sha256(Path(__file__).resolve())
    manifest_path = output / "run_manifest.json"
    progress_path = output / "progress.json"
    checkpoint_path = output / "frame_results.jsonl"

    if output.exists():
        if not resume:
            raise FileExistsError(f"refusing to overwrite existing output: {output}")
        if not manifest_path.is_file() or not checkpoint_path.is_file():
            raise ValueError("--resume requires the existing run manifest and frame checkpoint")
        saved_manifest = read_json(manifest_path)
        if saved_manifest != run_manifest:
            raise ValueError("existing output manifest differs from current inputs or script")
        if (output / "summary.json").is_file():
            raise ValueError("completed output cannot be resumed")
    else:
        if resume:
            raise FileNotFoundError("--resume requested but output directory does not exist")
        output.mkdir(parents=True, exist_ok=False)
        write_json_atomic(manifest_path, run_manifest)
        checkpoint_path.touch(exist_ok=False)

    completed = {}
    with checkpoint_path.open(encoding="utf-8") as stream:
        for line in stream:
            if not line.strip():
                continue
            record = json.loads(line)
            if record.get("status") == "complete":
                completed[(record["target_id"], record["filter"])] = record
    tasks = []
    frame_by_id = {str(row["id"]): row for row in target_manifest["clusters"]}
    record_by_id = {str(row["id"]): row for row in read_tsv(MASK_RECORDS)}
    for selected_item in selected:
        key = (selected_item["target_id"], selected_item["filter"])
        if key in completed:
            continue
        selected_item["members"] = members_by_target[selected_item["target_id"]]
        selected_item["frame"] = frame_by_id[selected_item["target_id"]]
        selected_item["record"] = record_by_id[selected_item["target_id"]]
        task = {"selected": selected_item, "manifest": target_manifest}
        if psf_bank is not None:
            psf, psf_filename, psf_sha256 = psf_bank.member(
                selected_item["target_id"], selected_item["filter"],
                selected_item["member_index"], selected_item["source_row"],
            )
            task.update({"member_psf": psf, "member_psf_filename": psf_filename,
                         "member_psf_sha256": psf_sha256})
        tasks.append(task)

    progress = {"status": "running", "completed_frames": len(completed),
                "total_frames": len(selected), "workers": workers,
                "science_measurement": False}
    write_json_atomic(progress_path, progress)
    started = time.time()
    new_completed = len(completed)
    failures = 0
    with checkpoint_path.open("a", encoding="utf-8") as stream:
        if tasks:
            with ProcessPoolExecutor(max_workers=workers) as pool:
                futures = {pool.submit(process_frame, task): task["selected"] for task in tasks}
                for future in as_completed(futures):
                    selected_item = futures[future]
                    try:
                        record = future.result()
                    except Exception as error:
                        failures += 1
                        record = {
                            "status": "failed", "target_id": selected_item["target_id"],
                            "filter": selected_item["filter"],
                            "source_row": selected_item["source_row"],
                            "error_type": type(error).__name__, "error": str(error),
                        }
                    stream.write(json.dumps(record, sort_keys=True, allow_nan=False) + "\n")
                    stream.flush()
                    new_completed += int(record.get("status") == "complete")
                    done_records = len(completed) + len([1])
                    progress.update({
                        "completed_frames": new_completed,
                        "finished_tasks": progress.get("finished_tasks", 0) + 1,
                        "total_frames": len(selected),
                        "last_target_id": selected_item["target_id"],
                        "last_status": record.get("status"),
                        "worker_failures": failures,
                        "science_measurement": False,
                    })
                    write_json_atomic(progress_path, progress)
                    if progress["finished_tasks"] % 10 == 0 or progress["finished_tasks"] == len(tasks):
                        print(f"member injection recovery {progress['finished_tasks']}/{len(tasks)} tasks; "
                              f"{new_completed}/{len(selected)} frames complete", flush=True)

    all_records = []
    with checkpoint_path.open(encoding="utf-8") as stream:
        for line in stream:
            if line.strip():
                all_records.append(json.loads(line))
    by_key = {}
    for record in all_records:
        by_key[(record.get("target_id"), record.get("filter"))] = record
    ordered_records = [by_key[(item["target_id"], item["filter"])]
                       for item in selected if (item["target_id"], item["filter"]) in by_key]
    if len(ordered_records) != len(selected):
        raise ValueError("checkpoint does not contain exactly one latest record per selected frame")
    response_rows = [row for record in ordered_records
                     for row in record.get("injections", [])]
    fields = [
        "target_id", "filter", "member_index", "source_row", "cluster_radius_kpc",
        "member_ra_deg", "member_dec_deg", "amplitude_at_50_kpc_njy_arcsec2",
        "exponential_scale_kpc", "nominal_unmasked_50_75_mean_njy_arcsec2",
        "psf_convolution_applied", "member_psf_filename", "member_psf_sha256",
        "member_psf_shape", "member_psf_negative_pixels_preserved",
        "baseline_delta_sb_njy_arcsec2", "baseline_target_sb_njy_arcsec2",
        "baseline_control_sb_njy_arcsec2", "baseline_target_coverage",
        "baseline_control_coverage_mean", "baseline_n_controls", "baseline_n_detections",
        "baseline_segmentation_pixels", "baseline_replay_delta_abs_njy_arcsec2",
        "injected_status", "injected_n_detections", "injected_segmentation_pixels",
        "injected_detection_count_change", "injected_segmentation_pixel_change",
        "injected_target_sb_njy_arcsec2", "injected_control_sb_njy_arcsec2",
        "injected_delta_sb_njy_arcsec2", "recovered_delta_increment_njy_arcsec2",
        "recovered_increment_over_nominal_input", "injected_target_coverage",
        "injected_control_coverage_mean", "injected_n_controls",
        "injected_member_match_arcsec", "injected_detection_number", "model_status",
        "masked_model_delta_njy_arcsec2", "masked_model_over_nominal_input",
        "recovered_minus_masked_model_njy_arcsec2", "science_measurement",
    ]
    write_tsv(output / "injection_responses.tsv", response_rows, fields)
    write_tsv(output / "selected_members.tsv", [
        {key: value for key, value in item.items()
         if key not in ("frame", "record", "members")}
        for item in selected
    ], ["target_id", "filter", "member_index", "source_row", "cluster_radius_kpc",
        "baseline_target_sb", "baseline_control_sb", "baseline_delta_sb",
        "baseline_target_coverage", "baseline_control_coverage", "baseline_n_controls",
        "baseline_detection_number", "member_x", "member_y", "member_ra_deg", "member_dec_deg"])

    response_summary = bootstrap_summary(response_rows)
    plot_path = output / "member_halo_injection_recovery.png"
    render_plot(plot_path, response_summary, psf_convolved=psf_bank is not None)
    n_frame_failures = sum(record.get("status") != "complete" for record in ordered_records)
    status = "complete_exploratory_member_halo_injection_recovery" if n_frame_failures == 0 else "complete_with_frame_failures"
    limitations = [
        "One already-profiled member per target/filter frame is sampled; the experiment does not estimate completeness for unmatched or faint CAMIRA candidates.",
        "Injected source Poisson noise is not added and the delivered VARIANCE plane is held fixed.",
        "The injection is made only in the target image; the rotated controls remain fixed sky positions in the same target image.",
        "The local rotated controls can contain cluster light and subtract some injected halo wings.",
        "The test calibrates measurement transfer for this model; it does not establish a physical satellite profile or an ICL association.",
        "The 50--75 kpc profile is conditional on the original detection, coverage, and rotated-control support selection.",
    ]
    if psf_bank is None:
        limitations.append("The circular exponential halo is not convolved with a position-dependent PSF.")
    else:
        limitations.append("Only the local signed coadd model-PSF core is applied; extended PSF wings, coadd-noise effects, and source-Poisson noise remain unmodeled.")
    output_hashes = {
        name: sha256(output / name)
        for name in ("run_manifest.json", "frame_results.jsonl", "injection_responses.tsv",
                     "selected_members.tsv", "member_halo_injection_recovery.png")
    }
    summary = {
        "schema": SCHEMA,
        "status": status,
        "science_measurement": False,
        "images_ready_for_science": False,
        "network_used": False,
        "data_downloaded": False,
        "n_frames_selected": len(selected),
        "n_frames_complete": sum(record.get("status") == "complete" for record in ordered_records),
        "n_frames_failed": n_frame_failures,
        "n_frames_by_filter": {band: sum(record.get("filter") == band and record.get("status") == "complete"
                                          for record in ordered_records) for band in ("i", "i2")},
        "n_member_amplitude_rows": len(response_rows),
        "n_supported_injection_rows": sum(row.get("injected_status") == "supported" for row in response_rows),
        "n_baseline_profile_replays": sum(record.get("status") == "complete" for record in ordered_records),
        "amplitude_response_by_filter": response_summary,
        "psf_convolution_applied": psf_bank is not None,
        "method": run_manifest,
        "limitations": limitations,
        "input_sha256": input_hashes,
        "runtime_code_sha256": run_manifest["runtime_code_sha256"],
        "analysis_script_sha256": run_manifest["script_sha256"],
        "output_sha256": output_hashes,
        "elapsed_seconds_this_invocation": time.time() - started,
        "science_measurement_gate": False,
    }
    write_json_atomic(output / "summary.json", summary)
    progress.update({"status": status, "completed_frames": summary["n_frames_complete"],
                     "failed_frames": n_frame_failures, "science_measurement": False})
    write_json_atomic(progress_path, progress)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--workers", type=int, default=12)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--psf-plan-dir", type=Path,
                        help="approved S23B member-profile PSF request plan directory")
    parser.add_argument("--psf-download-dir", type=Path,
                        help="audited S23B member-profile PSF response directory")
    parser.add_argument("--approved-psf-plan-sha256",
                        help="exact plan SHA-256 explicitly approved for PSF retrieval")
    args = parser.parse_args()
    if not 1 <= args.workers <= 24:
        parser.error("--workers must be between 1 and 24")
    if bool(args.psf_plan_dir) != bool(args.psf_download_dir):
        parser.error("--psf-plan-dir and --psf-download-dir must be supplied together")
    if args.psf_plan_dir and not args.approved_psf_plan_sha256:
        parser.error("PSF convolution requires the explicit approved plan SHA-256")
    if not args.psf_plan_dir and args.approved_psf_plan_sha256:
        parser.error("--approved-psf-plan-sha256 requires both PSF directories")
    psf_bank = (S23BCamiraMemberProfilePSFBank(
        args.psf_plan_dir.resolve(), args.psf_download_dir.resolve(),
        args.approved_psf_plan_sha256,
    ) if args.psf_plan_dir else None)
    output = args.output or (DEFAULT_PSF_OUTPUT if psf_bank is not None else DEFAULT_OUTPUT)
    result = run(output, workers=args.workers, resume=args.resume, psf_bank=psf_bank)
    print(json.dumps({
        "status": result["status"], "output": str(output.resolve()),
        "n_frames_complete": result["n_frames_complete"],
        "n_frames_failed": result["n_frames_failed"],
        "n_supported_injection_rows": result["n_supported_injection_rows"],
        "psf_convolution_applied": result["psf_convolution_applied"],
        "science_measurement": result["science_measurement"],
        "summary_sha256": sha256(output.resolve() / "summary.json"),
    }, indent=2), flush=True)


if __name__ == "__main__":
    main()
