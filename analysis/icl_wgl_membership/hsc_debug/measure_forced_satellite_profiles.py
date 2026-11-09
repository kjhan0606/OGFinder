#!/usr/bin/env python3
"""Measure satellite-position profiles where CAMIRA positions lack SEP detections.

The runner uses the existing HSC images and preparation masks. It measures
only CAMIRA positions with no unique SEP match within 1 arcsec and no SEP
detection within 3 arcsec. The central position is forced, while all detected
sources remain masked and azimuthal controls estimate the local baseline.
This expands a positional stack; it does not identify individual satellite
stellar envelopes or certify a science measurement.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import platform
import sys
import time

import numpy as np
from astropy import units as u
from astropy.cosmology import FlatLambdaCDM
from astropy.io import fits
from scipy.spatial import cKDTree

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import measure_satellite_outer_profiles as core
from download_pilot import BAD_BIT_POSITIONS, BRIGHT_OBJECT_BIT
from evaluate_prestack_pair_support import pair_block_labels

DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
TARGET_MANIFEST = DATA / "download_01" / "download_manifest.json"
SUPPORT_AUDIT_DIR = DATA / "satellite_forced_position_support_2026-11-06_03"
SUPPORT_TABLE = SUPPORT_AUDIT_DIR / "member_position_support.tsv"
SUPPORT_SUMMARY = SUPPORT_AUDIT_DIR / "summary.json"
DETECTED_PROFILE_DIR = DATA / "satellite_outer_profile_perbin_full_2026-10-18_01"
DETECTED_PROFILE_TABLE = DETECTED_PROFILE_DIR / "satellite_profiles.tsv"
ICL_PROFILE_TABLE = DATA / "icl_annulus_support_2026-10-18_01" / "annulus_pair_profiles.tsv"
CONTROL_MANIFEST = DATA / "control_download_01" / "download_manifest.json"
DEFAULT_OUTPUT = DATA / "satellite_forced_outer_profile_2026-11-06_01"
R_EDGES_KPC = core.R_EDGES_KPC
COVERAGE_FLOOR = 0.60
MIN_CONTROLS = 3
N_BOOT = 5000
SEED = 20261106


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty table: {path}")
    with path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t",
                                lineterminator="\n", extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def read_fits_products(frame: dict, prep_dir: Path):
    image_path = Path(frame["image"])
    mask_path = prep_dir / f"cluster_{frame['id']}" / "masks.fits.gz"
    info = core.validate_fits(image_path, frame)
    with fits.open(image_path, memmap=False, checksum=True) as hdul:
        image = np.asarray(hdul[info["science_hdu"]].data, dtype=np.float32)
        variance = np.asarray(hdul[info["variance_hdu"]].data, dtype=np.float32)
        bits = np.asarray(hdul[info["mask_hdu"]].data, dtype=np.int32)
    with fits.open(mask_path, memmap=False, checksum=True) as hdul:
        segmentation = np.asarray(hdul["SEGMENTATION"].data, dtype=np.int32)
        base = np.asarray(hdul["BASE_MASK"].data, dtype=bool)
    if not (image.shape == variance.shape == bits.shape == segmentation.shape == base.shape):
        raise ValueError(f"image/variance/mask shapes differ for cluster {frame['id']}")
    bad_bits = sum(1 << bit for bit in BAD_BIT_POSITIONS) | BRIGHT_OBJECT_BIT
    invalid = ~np.isfinite(image) | ~np.isfinite(variance) | (variance <= 0)
    if not np.array_equal(base, ((bits & bad_bits) != 0) | invalid):
        raise ValueError(f"base mask differs from original HSC quality bits for cluster {frame['id']}")
    return image, variance, base, segmentation, info, mask_path


def forced_profiles_for_frame(frame: dict, prep_dir: Path,
                              candidates: list[dict], *,
                              mask_other_unmatched_camira: bool = True,
                              coverage_floor: float = COVERAGE_FLOOR,
                              control_coverage_floor: float = COVERAGE_FLOOR
                              ) -> tuple[list[dict], dict, dict]:
    ident = str(frame["id"])
    prep = prep_dir / f"cluster_{ident}"
    detections_path, members_path = prep / "detections.tsv", prep / "members.tsv"
    prep_json_path = prep / "preparation.json"
    for path in (detections_path, members_path, prep_json_path):
        if not path.is_file():
            raise FileNotFoundError(path)
    preparation = core.read_json(prep_json_path)
    if preparation.get("science_measurement") is not False:
        raise ValueError(f"unexpected preparation science gate for {ident}")
    image, variance, base, segmentation, info, mask_path = read_fits_products(frame, prep_dir)
    detections = core.parse_detections(core.read_tsv(detections_path))
    members = core.read_tsv(members_path)
    matched, _ = core.match_members(members, detections, info["pixel_scale_arcsec"])
    matched_ids = {int(row["member_index"]) for row in matched}
    in_frame = [row for row in members if row["IN_FRAME"] == "True"]
    all_member_xy = np.asarray([[float(row["X_IMAGE"]) - 1., float(row["Y_IMAGE"]) - 1.]
                                for row in in_frame], dtype=float)
    member_indices = np.asarray([int(row["INDEX"]) for row in in_frame], dtype=int)
    member_tree = cKDTree(all_member_xy) if len(all_member_xy) else None
    unmatched = np.asarray([i for i, index in enumerate(member_indices)
                            if int(index) not in matched_ids], dtype=int)
    fallback_xy = all_member_xy[unmatched] if len(unmatched) else np.empty((0, 2))
    fallback_indices = member_indices[unmatched] if len(unmatched) else np.asarray([], dtype=int)
    fallback_tree_all = cKDTree(fallback_xy) if len(fallback_xy) else None
    dets = core.detection_arrays(detections, float(info["pixel_scale_arcsec"]))

    center_xy = np.asarray(info["center_pixel_xy"], dtype=float)
    pixel_scale = float(info["pixel_scale_arcsec"])
    kpc_per_pixel = float(COSMOLOGY.kpc_proper_per_arcmin(float(frame["z"])).value
                          * pixel_scale / 60.)
    conversion = float(core.nJy_factor(info))
    member_by_index = {int(row["INDEX"]): row for row in in_frame}
    result_rows, diagnostics = [], []
    for candidate in candidates:
        member_index = int(candidate["member_index"])
        member = member_by_index.get(member_index)
        if member is None:
            raise ValueError(f"forced member {ident}/{member_index} is no longer in frame")
        separation = float(member["DETECTION_SEPARATION_ARCSEC"])
        if member_index in matched_ids or not np.isfinite(separation) or separation <= 3.0:
            raise ValueError(f"forced sample definition changed for {ident}/{member_index}")
        xy = np.array([float(member["X_IMAGE"]) - 1., float(member["Y_IMAGE"]) - 1.])
        radius = math.hypot(xy[0] - center_xy[0], xy[1] - center_xy[1]) * kpc_per_pixel
        if not 50.0 <= radius < 500.0:
            raise ValueError(f"forced position left 50-500 kpc range for {ident}/{member_index}")

        keep = fallback_indices != member_index
        fallback_tree = cKDTree(fallback_xy[keep]) if np.any(keep) else None
        # Other unmatched CAMIRA positions may themselves be the faint signal
        # this forced-position stack is intended to recover. Masking them in
        # every target stamp can erase most of the annulus in crowded clusters.
        target_fallback_tree = fallback_tree if mask_other_unmatched_camira else None
        control_fallback_tree = fallback_tree_all if mask_other_unmatched_camira else None
        target = core.profile_at(image, variance, base, segmentation, dets,
                                 target_fallback_tree, xy, pixel_scale, kpc_per_pixel,
                                 conversion, None, coverage_floor)
        controls_xy = core.control_centers(xy, center_xy, radius / kpc_per_pixel,
                                           dets, member_tree, pixel_scale,
                                           kpc_per_pixel, image.shape)
        controls = []
        for point in controls_xy:
            control = core.profile_at(image, variance, base, segmentation, dets,
                                      control_fallback_tree, point, pixel_scale,
                                      kpc_per_pixel, conversion, None, control_coverage_floor)
            if control is not None:
                controls.append(control)

        candidate_rows = []
        for j, (lo, hi) in enumerate(zip(R_EDGES_KPC[:-1], R_EDGES_KPC[1:])):
            valid_controls = [row for row in controls
                              if row["coverage"][j] >= control_coverage_floor]
            if target["coverage"][j] < coverage_floor or len(valid_controls) < MIN_CONTROLS:
                continue
            control_sb = float(np.mean([row["sb"][j] for row in valid_controls]))
            candidate_rows.append({
                "cluster_id": ident, "filter": info["coadd_header_filter"],
                "redshift": float(frame["z"]), "richness": float(frame["richness"]),
                "member_index": member_index, "cluster_radius_kpc": radius,
                "nearest_sep_detection_arcsec": separation,
                "radial_bin_kpc": f"{lo:g}-{hi:g}", "n_controls": len(valid_controls),
                "target_sb_njy_arcsec2": float(target["sb"][j]),
                "control_sb_njy_arcsec2": control_sb,
                "delta_sb_njy_arcsec2": float(target["sb"][j] - control_sb),
                "target_err_nominal_njy_arcsec2": float(target["err"][j]),
                "target_coverage": float(target["coverage"][j]),
                "control_coverage": float(np.mean([row["coverage"][j] for row in valid_controls])),
                "target_n_pixels": int(target["n_pixels"][j]),
                "science_measurement": False,
            })
        result_rows.extend(candidate_rows)
        diagnostics.append({
            "cluster_id": ident, "filter": info["coadd_header_filter"],
            "member_index": member_index, "cluster_radius_kpc": radius,
            "nearest_sep_detection_arcsec": separation,
            "candidate_controls": len(controls_xy), "valid_control_profiles": len(controls),
            "retained_annuli": len(candidate_rows),
            "retained_bins": ",".join(row["radial_bin_kpc"] for row in candidate_rows),
            "target_coverage_by_bin": ",".join(f"{x:.4f}" for x in target["coverage"]),
            "mask_other_unmatched_camira": mask_other_unmatched_camira,
            "status": "retained" if candidate_rows else "no_annulus_passed_support",
        })
    provenance = {
        "image_sha256": core.sha256(Path(frame["image"])),
        "mask_sha256": core.sha256(mask_path),
        "detections_sha256": core.sha256(detections_path),
        "members_sha256": core.sha256(members_path),
        "preparation_sha256": core.sha256(prep_json_path),
    }
    frame_result = {
        "cluster_id": ident, "filter": info["coadd_header_filter"],
        "n_candidates_requested": len(candidates),
        "n_profiles_retained": sum(row["cluster_id"] == ident for row in result_rows),
        **provenance,
    }
    return result_rows, {"frame": frame_result, "candidates": diagnostics}, provenance


def partial_corr(x, y, z, log_richness):
    design = np.column_stack([np.ones(len(x)), z, log_richness])
    rx = x - design @ np.linalg.lstsq(design, x, rcond=None)[0]
    ry = y - design @ np.linalg.lstsq(design, y, rcond=None)[0]
    if np.std(rx) == 0 or np.std(ry) == 0:
        return float("nan")
    return float(np.corrcoef(rx, ry)[0, 1])


def block_bootstrap_association(rows: list[dict], band: str, satellite_field: str,
                                target_manifest: dict, control_manifest: dict,
                                replicas=N_BOOT, seed=SEED) -> dict:
    """Partial correlation with shared target/control connected sky blocks."""
    target_by_id = {str(row["id"]): row for row in target_manifest["clusters"]}
    control_by_id = {str(row["target_cluster_id"]): row for row in control_manifest["clusters"]}
    selected = [row for row in rows
                if row["filter"] == band and row[satellite_field] is not None]
    if len(selected) < 12:
        return {"status": "insufficient_support", "filter": band,
                "metric": satellite_field, "n_pairs": len(selected)}
    target_coords, control_coords = [], []
    for row in selected:
        target = target_by_id[row["cluster_id"]]
        control = control_by_id[row["cluster_id"]]
        target_coords.append({"ra_deg": float(target["ra_deg"]), "dec_deg": float(target["dec_deg"])})
        control_coords.append({"ra_deg": float(control["ra_deg"]), "dec_deg": float(control["dec_deg"])})
    _, labels, n_target_cells, n_control_cells, n_blocks, largest = pair_block_labels(
        target_coords, control_coords, np.ones(len(selected), dtype=bool))
    groups = np.unique(labels)
    by_block = {group: np.flatnonzero(labels == group) for group in groups}
    x = np.asarray([float(row["icl_delta_200_300"]) for row in selected])
    y = np.asarray([float(row[satellite_field]) for row in selected])
    z = np.asarray([float(row["redshift"]) for row in selected])
    richness = np.log(np.asarray([float(row["richness"]) for row in selected]))
    estimate = partial_corr(x, y, z, richness)
    rng = np.random.default_rng(seed)
    draws = []
    for _ in range(replicas):
        chosen = rng.choice(groups, size=len(groups), replace=True)
        indices = np.concatenate([by_block[group] for group in chosen])
        value = partial_corr(x[indices], y[indices], z[indices], richness[indices])
        if np.isfinite(value):
            draws.append(value)
    if len(draws) < replicas * 0.9:
        raise ValueError(f"too few valid {band} block-bootstrap draws: {len(draws)}/{replicas}")
    return {
        "status": "complete_exploratory", "filter": band,
        "metric": satellite_field, "n_pairs": len(selected),
        "n_connected_5deg_blocks": n_blocks,
        "n_target_cells": n_target_cells, "n_control_cells": n_control_cells,
        "largest_block_pairs": largest,
        "partial_r_z_logrichness": estimate,
        "block_bootstrap_ci95": np.percentile(draws, [2.5, 97.5]).tolist(),
        "bootstrap_fraction_positive": float(np.mean(np.asarray(draws) > 0)),
        "replicas_requested": replicas, "replicas_valid": len(draws),
        "seed": seed,
    }


def summarize_stacks(forced_rows: list[dict], detected_rows: list[dict]) -> tuple[list[dict], dict]:
    detected_by_key = {}
    for row in detected_rows:
        if row["radial_bin_kpc"] == "50-75":
            key = (str(row["cluster_id"]), row["filter"])
            detected_by_key.setdefault(key, []).append(float(row["delta_sb_njy_arcsec2"]))
    forced_by_key = {}
    for row in forced_rows:
        if row["radial_bin_kpc"] == "50-75":
            key = (str(row["cluster_id"]), row["filter"])
            forced_by_key.setdefault(key, []).append(float(row["delta_sb_njy_arcsec2"]))
    all_keys = set(detected_by_key) | set(forced_by_key)
    by_cluster = []
    for key in sorted(all_keys, key=lambda item: (item[1], int(item[0]))):
        detected = detected_by_key.get(key, [])
        forced = forced_by_key.get(key, [])
        combined = detected + forced
        by_cluster.append({
            "cluster_id": key[0], "filter": key[1],
            "n_detected_positions": len(detected), "n_forced_positions": len(forced),
            "n_combined_positions": len(combined),
            "detected_mean_delta_50_75": float(np.mean(detected)) if detected else None,
            "forced_mean_delta_50_75": float(np.mean(forced)) if forced else None,
            "combined_mean_delta_50_75": float(np.mean(combined)) if combined else None,
            "science_measurement": False,
        })
    return by_cluster, {"detected_by_key": detected_by_key, "forced_by_key": forced_by_key}


def build_associations(cluster_rows, icl_rows):
    icl = {(str(row["id"]), row["filter"]): row for row in icl_rows
           if row["radial_bin_kpc"] == "200-300"}
    output = []
    for row in cluster_rows:
        key = (row["cluster_id"], row["filter"])
        if key not in icl:
            continue
        reference = icl[key]
        record = {
            "cluster_id": key[0], "filter": key[1],
            "redshift": float(reference["z"]), "richness": float(reference["richness"]),
            "icl_delta_200_300": float(reference["delta_sb_njy_arcsec2"]),
            "satellite_delta_detected_50_75": row["detected_mean_delta_50_75"],
            "satellite_delta_forced_50_75": row["forced_mean_delta_50_75"],
            "satellite_delta_combined_50_75": row["combined_mean_delta_50_75"],
            "n_detected_positions": row["n_detected_positions"],
            "n_forced_positions": row["n_forced_positions"],
            "icl_common_coverage": float(reference["common_coverage"]),
            "science_measurement": False,
        }
        output.append(record)
    return output


def run(output: Path, replicas: int, seed: int,
        mask_other_unmatched_camira: bool = True,
        coverage_floor: float = COVERAGE_FLOOR,
        control_coverage_floor: float = COVERAGE_FLOOR) -> dict:
    if output.exists():
        raise FileExistsError(f"refusing to overwrite {output}")
    target_manifest = json.loads(TARGET_MANIFEST.read_text(encoding="utf-8"))
    support_summary = json.loads(SUPPORT_SUMMARY.read_text(encoding="utf-8"))
    if support_summary.get("status") != "complete_exploratory" or support_summary.get("science_measurement") is not False:
        raise ValueError("forced-position support audit is incomplete or invalid")
    if support_summary.get("output_sha256", {}).get("member_position_support.tsv") != sha256(SUPPORT_TABLE):
        raise ValueError("forced-position support table hash differs from its audit summary")
    support_rows = [row for row in read_tsv(SUPPORT_TABLE)
                    if row["category"] == "no_sep_detection_within_3arcsec"]
    if len(support_rows) != 3197:
        raise ValueError(f"unexpected forced-position support: {len(support_rows)}")
    by_cluster = {}
    for row in support_rows:
        by_cluster.setdefault(str(row["cluster_id"]), []).append(row)
    frames = core.resolve_image_paths(target_manifest["clusters"], TARGET_MANIFEST)
    prep_dir = DATA / "satellite_detection_psf_2026-10-18_01"
    all_rows, all_diagnostics, frame_records = [], [], []
    hashes = {
        str(TARGET_MANIFEST): sha256(TARGET_MANIFEST),
        str(SUPPORT_SUMMARY): sha256(SUPPORT_SUMMARY),
        str(SUPPORT_TABLE): sha256(SUPPORT_TABLE),
        str(DETECTED_PROFILE_TABLE): sha256(DETECTED_PROFILE_TABLE),
        str(ICL_PROFILE_TABLE): sha256(ICL_PROFILE_TABLE),
        str(CONTROL_MANIFEST): sha256(CONTROL_MANIFEST),
    }
    start = time.time()
    n_frames_processed = 0
    for frame in frames:
        ident = str(frame["id"])
        candidates = by_cluster.get(ident, [])
        if not candidates:
            continue
        print(f"[{n_frames_processed + 1}] forced profiles for cluster {ident} ({len(candidates)} positions)", flush=True)
        rows, details, provenance = forced_profiles_for_frame(
            frame, prep_dir, candidates,
            mask_other_unmatched_camira=mask_other_unmatched_camira,
            coverage_floor=coverage_floor,
            control_coverage_floor=control_coverage_floor)
        all_rows.extend(rows)
        all_diagnostics.extend(details["candidates"])
        frame_records.append(details["frame"])
        for key, value in provenance.items():
            hashes[f"cluster_{ident}/{key}"] = value
        n_frames_processed += 1

    if len(all_diagnostics) != len(support_rows):
        raise ValueError(f"processed {len(all_diagnostics)} of {len(support_rows)} forced candidates")
    detected_rows = read_tsv(DETECTED_PROFILE_TABLE)
    cluster_profiles, _ = summarize_stacks(all_rows, detected_rows)
    association_rows = build_associations(cluster_profiles, read_tsv(ICL_PROFILE_TABLE))
    association_metrics = {"forced_only": {}, "combined": {}}
    control_manifest = json.loads(CONTROL_MANIFEST.read_text(encoding="utf-8"))
    for band in ("i", "i2"):
        association_metrics["forced_only"][band] = block_bootstrap_association(
            association_rows, band, "satellite_delta_forced_50_75",
            target_manifest, control_manifest, replicas,
            seed + (0 if band == "i" else 10))
        association_metrics["combined"][band] = block_bootstrap_association(
            association_rows, band, "satellite_delta_combined_50_75",
            target_manifest, control_manifest, replicas,
            seed + (1 if band == "i" else 11))

    results = {
        "schema": "hsc-forced-satellite-outer-profile-v1",
        "status": "complete_exploratory",
        "science_measurement": False,
        "n_candidate_positions": len(support_rows),
        "coverage_floor": coverage_floor,
        "n_target_frames_processed": n_frames_processed,
        "n_forced_candidate_bins_retained": len(all_rows),
        "n_candidates_with_any_retained_bin": sum(bool(row["retained_annuli"]) for row in all_diagnostics),
        "n_candidates_retained_50_75": sum(row["radial_bin_kpc"] == "50-75" for row in all_rows),
        "forced_profiles_by_filter_bin": {},
        "association_support_by_filter": {
            band: {
                "forced_only": sum(row["filter"] == band and row["satellite_delta_forced_50_75"] is not None for row in association_rows),
                "combined": sum(row["filter"] == band and row["satellite_delta_combined_50_75"] is not None for row in association_rows),
            } for band in ("i", "i2")
        },
        "association_metrics": association_metrics,
        "method": {
            "forced_position_selection": "in-frame CAMIRA member, 50-500 proper kpc, no unique SEP match within 1 arcsec and nearest SEP detection farther than 3 arcsec",
            "target_measurement": (
                "forced at CAMIRA coordinates; mask all SEP detections and original HSC quality/bright-star bits; "
                + ("also mask other unmatched CAMIRA positions" if mask_other_unmatched_camira
                   else "leave other unmatched CAMIRA positions unmasked so candidate satellite light is not self-subtracted")
            ),
            "mask_other_unmatched_camira": mask_other_unmatched_camira,
            "local_baseline": "up to eight azimuthally rotated positions at the same cluster-centric radius, excluding positions near detections or CAMIRA members; average valid control profiles by annulus",
            "target_coverage_floor": coverage_floor,
            "control_coverage_floor": control_coverage_floor,
            "annulus_support": f"target coverage >= {coverage_floor:.2f}; at least {MIN_CONTROLS} controls with coverage >= {control_coverage_floor:.2f} per annulus",
            "stacking": "mean member profile within cluster/filter/radius, then equal-weight cluster; combine forced positions with saved detected-position profiles only after cluster means are formed",
            "association": "200-300 proper-kpc ICL contrast versus 50-75 proper-kpc satellite-position contrast; adjust redshift and log richness; connected target/control 5-degree block bootstrap",
        },
        "limitations": [
            "Forced positions without nearby SEP detections are not guaranteed true/isolated satellites; CAMIRA selection and blends remain.",
            "Their centers are not independently verified galaxy centroids, so forced stacks measure mean light around CAMIRA positions.",
            "The profile mask, BCG center, PSF wings, sky transfer, source-Poisson, and total coadd covariance are not science-certified.",
            "Association samples remain coverage-selected; the missing-at-random assumption is unverified and no familywise adjustment across analyses is claimed.",
            "A 5-degree connected-block partition is approximate rather than a validated HSC field covariance.",
        ],
        "elapsed_seconds": time.time() - start,
        "input_hashes": hashes,
    }
    for band in ("i", "i2"):
        for radial_bin in sorted({row["radial_bin_kpc"] for row in all_rows},
                                 key=lambda value: (float(value.split("-")[0]), value)):
            selected = [row for row in all_rows if row["filter"] == band and row["radial_bin_kpc"] == radial_bin]
            results["forced_profiles_by_filter_bin"].setdefault(band, {})[radial_bin] = {
                "n_member_profiles": len(selected),
                "n_clusters": len({row["cluster_id"] for row in selected}),
            }
    output.mkdir(parents=True, exist_ok=False)
    write_tsv(output / "forced_satellite_profiles.tsv", all_rows)
    write_tsv(output / "candidate_diagnostics.tsv", all_diagnostics)
    write_tsv(output / "cluster_profiles.tsv", cluster_profiles)
    write_tsv(output / "icl_satellite_association.tsv", association_rows)
    results["output_sha256"] = {
        name: sha256(output / name) for name in (
            "forced_satellite_profiles.tsv", "candidate_diagnostics.tsv",
            "cluster_profiles.tsv", "icl_satellite_association.tsv")
    }
    results["script_sha256"] = sha256(Path(__file__).resolve())
    (output / "summary.json").write_text(json.dumps(results, indent=2, ensure_ascii=False) + "\n",
                                         encoding="utf-8")
    return results


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--bootstrap-replicas", type=int, default=N_BOOT)
    parser.add_argument("--seed", type=int, default=SEED)
    parser.add_argument("--coverage-floor", type=float, default=COVERAGE_FLOOR,
                        help="minimum usable area fraction per target/control annulus (default 0.60)")
    parser.add_argument("--control-coverage-floor", type=float, default=COVERAGE_FLOOR,
                        help="minimum usable area fraction for each rotated control annulus (default 0.60)")
    parser.add_argument("--leave-other-unmatched-camira-unmasked", action="store_false",
                        dest="mask_other_unmatched_camira",
                        help="sensitivity mode: do not mask neighboring unmatched CAMIRA positions in target or control stamps")
    parser.set_defaults(mask_other_unmatched_camira=True)
    args = parser.parse_args()
    if not platform.node().lower().startswith("lageunha"):
        parser.error("Run CPU image processing on lageunha")
    if args.bootstrap_replicas < 100:
        parser.error("at least 100 bootstrap replicates are required")
    if not (0.0 < args.coverage_floor <= 1.0
            and 0.0 < args.control_coverage_floor <= 1.0):
        parser.error("coverage floors must be in (0, 1]")
    result = run(args.output, args.bootstrap_replicas, args.seed,
                 args.mask_other_unmatched_camira, args.coverage_floor,
                 args.control_coverage_floor)
    print(json.dumps({key: result[key] for key in (
        "status", "output", "n_candidate_positions", "coverage_floor", "n_candidates_with_any_retained_bin",
        "n_candidates_retained_50_75", "forced_profiles_by_filter_bin",
        "association_support_by_filter", "association_metrics", "science_measurement")
        if key in result}, indent=2))


COSMOLOGY = FlatLambdaCDM(H0=70.0, Om0=0.3, Tcmb0=0 * u.K)

if __name__ == "__main__":
    main()
