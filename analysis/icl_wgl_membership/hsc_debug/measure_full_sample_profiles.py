#!/usr/bin/env python3
"""Measure exploratory matched-control profiles for the existing 800 HSC-I cutouts.

Uses the already downloaded global-sky coadd/bg images, the reviewed automatic
mask candidates, and existing matched-control profiles. It never downloads
data, changes source images, or activates a science-measurement flag.
"""
import argparse
import ast
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import sys
import time

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
from astropy.io import fits

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from analyze_satellite_association_pilot import association_stats
from download_pilot import file_sha256 as sha256, validate_fits
from evaluate_prestack_pair_support import pair_block_labels
from prestack_core import EDGES, FACTORS, covariance_diagnostics, geometry, measure
from validate_control_images import nJy_factor
from validate_prestack import surface_result


SCOPES = {"50-300": np.arange(4), "50-500": np.arange(6)}
CONFIGS = [(2.5, 0.40, "50-300"), (2.5, 0.40, "50-500"),
           (2.5, 0.50, "50-300"), (3.5, 0.50, "50-300")]
REPLICAS = 2000


def read_json(path):
    with Path(path).open(encoding="utf-8") as stream:
        return json.load(stream)


def read_tsv(path):
    with Path(path).open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def write_json(path, value):
    target = Path(path)
    temporary = target.with_suffix(target.suffix + ".tmp")
    temporary.write_text(json.dumps(value, ensure_ascii=False, indent=2,
                                    sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    temporary.replace(target)


def write_tsv(path, rows, fields):
    target = Path(path)
    temporary = Path(str(target) + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(target)


def safe(value):
    if isinstance(value, dict):
        return {str(k): safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [safe(v) for v in value]
    if isinstance(value, np.ndarray):
        return safe(value.tolist())
    if isinstance(value, (np.bool_, bool)):
        return bool(value)
    if isinstance(value, (np.integer, int)):
        return int(value)
    if isinstance(value, (np.floating, float)):
        return float(value) if np.isfinite(value) else None
    return value


def plain_profile(result, factor):
    out = dict(result)
    for key in ("raw", "corrected", "plane"):
        out[key] = result[key] * factor
    for key in ("nominal_raw", "nominal_corrected"):
        out[key] = result[key] * factor ** 2
    return out


def process_cluster(frame, mask_record, mask_summary, cosmo):
    ident = str(frame["id"])
    image_path = Path(mask_record["image"])
    mask_path = Path(mask_record["mask_file"])
    if not image_path.is_file() or image_path.stat().st_size != int(mask_record["image_bytes"]):
        raise ValueError(f"image missing or byte count changed for {ident}")
    image_hash = sha256(image_path)
    if image_hash != mask_record["image_sha256"]:
        raise ValueError(f"image SHA-256 changed for {ident}")
    expected_mask = mask_summary.get("mask_output_sha256", {}).get(str(mask_path))
    if not mask_path.is_file() or expected_mask != mask_record["mask_sha256"]:
        raise ValueError(f"mask provenance mismatch for {ident}")
    mask_hash = sha256(mask_path)
    if mask_hash != expected_mask:
        raise ValueError(f"mask SHA-256 changed for {ident}")

    info = validate_fits(image_path, frame)
    if info["coadd_header_filter"] != mask_record["primary_filter"]:
        raise ValueError(f"FITS filter differs from mask record for {ident}")
    factor = nJy_factor(info)
    center = info["center_pixel_xy"]
    scale = info["pixel_scale_arcsec"]
    kpp = cosmo.kpc_proper_per_arcmin(float(frame["z"])).value * scale / 60.0
    with fits.open(image_path, memmap=False, checksum=True) as hdul:
        image = np.asarray(hdul[info["science_hdu"]].data, dtype=np.float32)
        variance = np.asarray(hdul[info["variance_hdu"]].data, dtype=np.float32)
    with fits.open(mask_path, memmap=False, checksum=True) as hdul:
        masks = [np.asarray(hdul["MASK_K" + str(k).replace(".", "P")].data, dtype=bool)
                 for k in FACTORS]
    if any(mask.shape != image.shape for mask in masks) or variance.shape != image.shape:
        raise ValueError(f"image/variance/mask shapes differ for {ident}")
    geom = geometry(image.shape, center, kpp)
    profiles = [plain_profile(measure(image, variance, mask, geom, center, kpp), factor)
                for mask in masks]
    arrays = {key: np.asarray([profile[key] for profile in profiles]) for key in
              ("raw", "corrected", "coverage", "count", "nominal_raw",
               "nominal_corrected", "plane", "n_tiles")}
    row = {
        "id": ident,
        "ra_deg": float(frame["ra_deg"]),
        "dec_deg": float(frame["dec_deg"]),
        "z": float(frame["z"]),
        "richness": float(frame["richness"]),
        "filter": info["coadd_header_filter"],
        "n_member_positions": len(frame.get("members_radec_deg", [])),
        "auto_unique_center_candidate": mask_record["auto_unique_center_candidate"] == "True",
        "bcg_match_count": int(mask_record["bcg_match_count"]),
        "center_masked_by_kron": ast.literal_eval(mask_record["center_masked_by_kron"]),
        "image_sha256": image_hash,
        "mask_sha256": mask_hash,
        "native_to_njy_arcsec2": float(factor),
        "pixel_scale_arcsec": float(scale),
        "kpc_per_pixel": float(kpp),
    }
    return arrays, row


def member_counts(frame):
    positions = np.asarray(frame.get("members_radec_deg", []), dtype=float)
    if positions.ndim != 2 or positions.shape[1] != 2:
        raise ValueError(f"member coordinates malformed for {frame['id']}")
    if not len(positions):
        return np.zeros(6, dtype=int)
    centre = SkyCoord(float(frame["ra_deg"]) * u.deg, float(frame["dec_deg"]) * u.deg)
    members = SkyCoord(positions[:, 0] * u.deg, positions[:, 1] * u.deg)
    cosmo = FlatLambdaCDM(H0=70.0, Om0=0.3, Tcmb0=0)
    radii = centre.separation(members).rad * cosmo.angular_diameter_distance(float(frame["z"])).to_value(u.kpc)
    return np.asarray([np.count_nonzero((radii >= lo) & (radii < hi))
                       for lo, hi in zip(EDGES[:-1], EDGES[1:])], dtype=int)


def pair_selection(cluster, control, factor_index, floor, bins):
    coverage = np.minimum(cluster["coverage"][:, factor_index, bins],
                          control["coverage"][:, factor_index, bins])
    finite = np.isfinite(cluster["raw"][:, factor_index, bins]).all(axis=1)
    finite &= np.isfinite(control["raw"][:, factor_index, bins]).all(axis=1)
    return np.all(coverage >= floor, axis=1) & finite


def make_stack(population, rows, control_rows, cluster, control, factor_index,
               floor, scope, mode, band, seed):
    bins = SCOPES[scope]
    band_selection = np.asarray([row["filter"] == band and control_rows[i]["filter"] == band
                                 for i, row in enumerate(rows)], dtype=bool)
    selected = pair_selection(cluster, control, factor_index, floor, bins) & band_selection
    target_values = cluster[mode][:, factor_index, bins]
    control_values = control[mode][:, factor_index, bins]
    differences = target_values - control_values
    nominal = cluster["nominal_" + mode][:, factor_index][:, bins][:, :, bins]
    nominal = nominal + control["nominal_" + mode][:, factor_index][:, bins][:, :, bins]
    selected &= np.isfinite(differences).all(axis=1) & np.isfinite(nominal).all(axis=(1, 2))
    if not np.any(selected):
        return None, selected
    indices, labels, target_cells, control_cells, blocks, largest = pair_block_labels(
        rows, control_rows, selected)
    info, covariance = covariance_diagnostics(
        differences[selected], labels, nominal[selected], seed=seed, replicas=REPLICAS)
    mean_target = np.mean(target_values[selected], axis=0)
    mean_control = np.mean(control_values[selected], axis=0)
    mean_difference = np.mean(differences[selected], axis=0)
    cov = covariance.get("spatial_bootstrap_covariance")
    standard_error = np.sqrt(np.diag(cov)).tolist() if cov is not None and np.isfinite(cov).all() else [None] * len(bins)
    info.update(selected_indices=indices.tolist(), n_target_sky_cells=target_cells,
                n_control_sky_cells=control_cells, n_joint_blocks=blocks,
                largest_joint_block=largest, band=band, kron_factor=FACTORS[factor_index],
                coverage_floor=floor, radial_scope_kpc=scope, extraction=mode,
                mean_target_sb=mean_target.tolist(), mean_control_sb=mean_control.tolist(),
                mean_pair_difference=mean_difference.tolist(),
                spatial_bootstrap_standard_error=standard_error)
    for key, value in covariance.items():
        covariance[key] = value
    return (info, covariance, differences, selected), selected


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", type=Path, required=True)
    parser.add_argument("--work", type=Path, required=True)
    args = parser.parse_args()
    if not os.environ.get("SLURM_JOB_ID"):
        parser.error("Run this full-sample measurement inside a Slurm allocation")
    sample = args.sample.resolve()
    work = args.work.resolve()
    if work.exists():
        parser.error(f"Refusing to overwrite existing output: {work}")

    sample_manifest_path = sample / "sample_manifest.json"
    mask_dir = sample / "prestack_cluster_masks_2026-10-16_01"
    mask_summary_path = mask_dir / "summary.json"
    mask_records_path = mask_dir / "records.tsv"
    control_dir = sample / "prestack_extraction_2026-10-15_01"
    control_summary_path = control_dir / "extraction_summary.json"
    area_dir = sample / "prestack_full_pair_area_support_2026-10-16_01"
    area_summary_path = area_dir / "summary.json"
    area_support_path = area_dir / "support_counts.tsv"
    paths = [sample_manifest_path, mask_summary_path, mask_records_path,
             control_summary_path, area_summary_path, area_support_path]
    if any(not path.is_file() for path in paths):
        parser.error("A required manifest or support artifact is missing")
    sample_doc = read_json(sample_manifest_path)
    mask_summary = read_json(mask_summary_path)
    control_summary = read_json(control_summary_path)
    area_summary = read_json(area_summary_path)
    if len(sample_doc.get("clusters", [])) != 800 or mask_summary.get("n_success") != 800:
        parser.error("Expected the complete 800-target sample and mask set")
    if mask_summary.get("status") != "complete" or mask_summary.get("mask_ready_for_science") is not False:
        parser.error("Mask provenance/status is unexpected")
    if control_summary.get("status") != "complete" or control_summary.get("n_measured") != 798:
        parser.error("Expected the complete 798-field control-profile set")
    if area_summary.get("status") != "complete" or area_summary.get("science_measurement") is not False:
        parser.error("Full-sample area-support provenance is unexpected")
    if sha256(mask_records_path) != mask_summary["output_sha256"][str(mask_records_path)]:
        parser.error("Mask record table SHA-256 does not match its summary")

    frames = sample_doc["clusters"]
    by_id = {str(row["id"]): row for row in frames}
    if len(by_id) != 800:
        parser.error("Sample manifest contains duplicate cluster identifiers")
    masks_by_id = {str(row["id"]): row for row in read_tsv(mask_records_path)}
    if len(masks_by_id) != 800:
        parser.error("Mask record table does not contain exactly 800 rows")
    control_records = control_summary.get("records", [])
    controls_by_target = {}
    for row in control_records:
        target_id = str(row["target_cluster_id"])
        if target_id in controls_by_target:
            parser.error(f"Multiple controls linked to target {target_id}")
        controls_by_target[target_id] = row
    if len(controls_by_target) != 798:
        parser.error("Control extraction does not contain 798 unique target links")

    work.mkdir(parents=True, exist_ok=False)
    start = time.monotonic()
    source_hashes = {str(path): sha256(path) for path in paths}
    profile_rows = []
    profile_parts = {key: [] for key in ("raw", "corrected", "coverage", "count",
                                         "nominal_raw", "nominal_corrected", "plane", "n_tiles")}
    error_file = work / "errors.json"
    progress_file = work / "progress.json"
    cosmo = FlatLambdaCDM(H0=float(sample_doc["h0"]), Om0=float(sample_doc["omega_m"]), Tcmb0=0)
    try:
        for index, frame in enumerate(frames, 1):
            ident = str(frame["id"])
            arrays, row = process_cluster(frame, masks_by_id[ident], mask_summary, cosmo)
            for key in profile_parts:
                profile_parts[key].append(arrays[key])
            profile_rows.append(row)
            source_hashes[row["image"] if "image" in row else masks_by_id[ident]["image"]] = row["image_sha256"]
            source_hashes[masks_by_id[ident]["mask_file"]] = row["mask_sha256"]
            if index % 10 == 0 or index == len(frames):
                write_json(progress_file, {"status": "running", "completed": index,
                                          "expected": len(frames), "elapsed_seconds": time.monotonic() - start})
                print(f"Measured {index}/{len(frames)} targets", flush=True)
        cluster = {key: np.asarray(value) for key, value in profile_parts.items()}
        if cluster["raw"].shape != (800, 3, 6):
            raise ValueError(f"Unexpected target profile cube shape {cluster['raw'].shape}")

        control_rows = []
        control_parts = {key: [] for key in ("raw", "corrected", "coverage", "nominal_raw", "nominal_corrected")}
        for frame in frames:
            ident = str(frame["id"])
            record = controls_by_target.get(ident)
            if record is None:
                control_rows.append({"filter": "", "ra_deg": frame["ra_deg"], "dec_deg": frame["dec_deg"]})
                for key in control_parts:
                    shape = (3, 6) if key not in ("nominal_raw", "nominal_corrected") else (3, 6, 6)
                    control_parts[key].append(np.full(shape, np.nan))
                continue
            control_id = record["control_id"]
            baseline = control_dir / control_id / "baseline.npz"
            expected = record["baseline_sha256"]
            if not baseline.is_file() or sha256(baseline) != expected:
                raise ValueError(f"Control baseline checksum failed for {control_id}")
            source_hashes[str(baseline)] = expected
            with np.load(baseline) as archive:
                for key in control_parts:
                    control_parts[key].append(archive[key].copy())
            control_rows.append({"filter": record["header_filter"], "ra_deg": record["ra_deg"],
                                 "dec_deg": record["dec_deg"], "control_id": control_id})
        control = {key: np.asarray(value) for key, value in control_parts.items()}
        if control["raw"].shape != (800, 3, 6):
            raise ValueError(f"Unexpected control profile cube shape {control['raw'].shape}")

        # Satellite-position counts use proper projected radii; no member flux or
        # membership-probability weights are available in this manifest.
        for frame, row in zip(frames, profile_rows):
            row["satellite_counts"] = member_counts(frame).tolist()

        stack_rows = []
        covariance_archive = {}
        candidate_metrics = []
        for config_index, (factor, floor, scope) in enumerate(CONFIGS):
            factor_index = FACTORS.index(factor)
            bins = SCOPES[scope]
            for band_index, band in enumerate(("i", "i2")):
                for mode_index, mode in enumerate(("raw", "corrected")):
                    result, selected = make_stack(
                        profile_rows, profile_rows, control_rows, cluster, control,
                        factor_index, floor, scope, mode, band,
                        seed=20261018 + config_index * 100 + band_index * 10 + mode_index)
                    if result is None:
                        stack_rows.append({"config": f"k{factor}_c{int(floor*100)}_r{scope}",
                                           "filter": band, "extraction": mode, "n_pairs": 0,
                                           "status": "no_supported_pairs"})
                        continue
                    info, cov, differences, selected = result
                    name = f"k{factor}_c{int(floor*100)}_r{scope}_{band}_{mode}"
                    for key in ("sample_covariance", "independent_stack_covariance",
                                "spatial_bootstrap_covariance", "delete_cell_covariance",
                                "mean"):
                        if key in cov:
                            covariance_archive[name + "__" + key] = cov[key]
                    covariance_archive[name + "__selected_indices"] = np.asarray(info["selected_indices"], int)
                    annuli = list(bins)
                    for local, bin_index in enumerate(annuli):
                        stack_rows.append({
                            "config": f"k{factor}_c{int(floor*100)}_r{scope}",
                            "kron_factor": factor, "coverage_floor": floor, "radial_scope_kpc": scope,
                            "filter": band, "extraction": mode, "n_pairs": info["n"],
                            "n_joint_blocks": info["n_joint_blocks"],
                            "r_min_kpc": EDGES[bin_index], "r_max_kpc": EDGES[bin_index + 1],
                            "mean_target_sb_njy_arcsec2": info["mean_target_sb"][local],
                            "mean_control_sb_njy_arcsec2": info["mean_control_sb"][local],
                            "mean_pair_difference_njy_arcsec2": info["mean_pair_difference"][local],
                            "spatial_bootstrap_se_njy_arcsec2": info["spatial_bootstrap_standard_error"][local],
                            "covariance_rank": info.get("rank"), "status": info["status"],
                        })
                    if factor == 2.5 and floor == 0.40 and scope == "50-300" and mode == "raw":
                        for index in np.flatnonzero(selected):
                            frame_row = profile_rows[index]
                            ctrl_row = control_rows[index]
                            pair_metrics_row = {
                                "id": frame_row["id"], "filter": band,
                                "control_id": ctrl_row.get("control_id", ""),
                                "richness": frame_row["richness"], "z": frame_row["z"],
                                "auto_unique_center_candidate": frame_row["auto_unique_center_candidate"],
                                "center_masked_k2p5": frame_row["center_masked_by_kron"]["2.5"],
                                "satellite_counts": frame_row["satellite_counts"][:4],
                                "satellite_density_mpc2": (np.asarray(frame_row["satellite_counts"][:4]) /
                                    (np.pi * np.diff(EDGES[:5] ** 2) / 1.0e6)).tolist(),
                                "delta_sb_njy_arcsec2": (cluster["raw"][index, factor_index, :4] -
                                                         control["raw"][index, factor_index, :4]).tolist(),
                                "common_coverage": np.minimum(cluster["coverage"][index, factor_index, :4],
                                                              control["coverage"][index, factor_index, :4]).tolist(),
                            }
                            candidate_metrics.append(pair_metrics_row)

        profile_fields = ["id", "ra_deg", "dec_deg", "z", "richness", "filter", "n_member_positions",
                          "auto_unique_center_candidate", "bcg_match_count", "center_masked_by_kron",
                          "native_to_njy_arcsec2", "pixel_scale_arcsec", "kpc_per_pixel",
                          "coverage_k2p5_50_100", "coverage_k2p5_100_150", "coverage_k2p5_150_200",
                          "coverage_k2p5_200_300", "coverage_k2p5_300_400", "coverage_k2p5_400_500",
                          "members_50_100", "members_100_150", "members_150_200", "members_200_300",
                          "members_300_400", "members_400_500"]
        profile_table = []
        for index, row in enumerate(profile_rows):
            flat = {key: row[key] for key in profile_fields[:13]}
            flat["center_masked_by_kron"] = json.dumps(row["center_masked_by_kron"], sort_keys=True)
            factor_index = FACTORS.index(2.5)
            for b, name in enumerate(("50_100", "100_150", "150_200", "200_300", "300_400", "400_500")):
                flat["coverage_k2p5_" + name] = float(cluster["coverage"][index, factor_index, b])
                flat["members_" + name] = int(row["satellite_counts"][b])
            profile_table.append(flat)
        write_tsv(work / "cluster_profiles.tsv", profile_table, profile_fields)
        stack_fields = ["config", "kron_factor", "coverage_floor", "radial_scope_kpc", "filter", "extraction",
                        "n_pairs", "n_joint_blocks", "r_min_kpc", "r_max_kpc", "mean_target_sb_njy_arcsec2",
                        "mean_control_sb_njy_arcsec2", "mean_pair_difference_njy_arcsec2",
                        "spatial_bootstrap_se_njy_arcsec2", "covariance_rank", "status"]
        write_tsv(work / "stack_profiles.tsv", stack_rows, stack_fields)
        pair_fields = ["id", "filter", "control_id", "richness", "z", "auto_unique_center_candidate",
                       "center_masked_k2p5", "members_50_100", "members_100_150", "members_150_200",
                       "members_200_300", "density_mpc2_50_100", "density_mpc2_100_150",
                       "density_mpc2_150_200", "density_mpc2_200_300", "delta_sb_50_100",
                       "delta_sb_100_150", "delta_sb_150_200", "delta_sb_200_300",
                       "coverage_50_100", "coverage_100_150", "coverage_150_200", "coverage_200_300"]
        flat_pairs = []
        for row in candidate_metrics:
            flat = {key: row[key] for key in ("id", "filter", "control_id", "richness", "z",
                                               "auto_unique_center_candidate", "center_masked_k2p5")}
            for b, name in enumerate(("50_100", "100_150", "150_200", "200_300")):
                flat["members_" + name] = row["satellite_counts"][b]
                flat["density_mpc2_" + name] = row["satellite_density_mpc2"][b]
                flat["delta_sb_" + name] = row["delta_sb_njy_arcsec2"][b]
                flat["coverage_" + name] = row["common_coverage"][b]
            flat_pairs.append(flat)
        write_tsv(work / "pair_metrics_k2p5_c40_r50_300.tsv", flat_pairs, pair_fields)
        np.savez_compressed(work / "cluster_profiles.npz", **cluster)
        np.savez_compressed(work / "stack_covariances.npz", **covariance_archive)

        # Add empirical profile covariance and sample diagnostics for the principal exploratory cut.
        assoc_rows = []
        for row in candidate_metrics:
            assoc_rows.append({"filter": row["filter"], "richness": row["richness"], "z": row["z"],
                               "satellite_density_mpc2": row["satellite_density_mpc2"],
                               "delta_sb_njy_arcsec2": row["delta_sb_njy_arcsec2"]})
        association = {band: association_stats(assoc_rows, band) for band in ("i", "i2", "all")}
        summary = {
            "schema": "hsc-icl-full-sample-profiles-v1", "status": "complete",
            "n_targets_measured": len(profile_rows), "n_matched_controls": len(control_records),
            "n_matched_filter_pairs": sum(row["filter"] == control_rows[i]["filter"]
                                           for i, row in enumerate(profile_rows)),
            "filter_counts_targets": {band: sum(row["filter"] == band for row in profile_rows) for band in ("i", "i2")},
            "filter_counts_controls": {band: sum(row["header_filter"] == band for row in control_records) for band in ("i", "i2")},
            "pixel_calibration": "FLUXMAG0-derived nJy/arcsec2; no Galactic extinction correction or i/i2 response cross-calibration",
            "measurement": "Kron source-mask native-pixel annular means; target minus matched control; equal pair weighting",
            "radial_bins_kpc": EDGES.tolist(), "stack_configurations": CONFIGS,
            "stack_profiles": stack_rows, "member_position_association": association,
            "science_gates": {"science_measurement": False, "sky_correction_approved": False,
                              "psf_approved": False, "mask_ready_for_science": False,
                              "extinction_crosscheck_passed": False, "weak_lensing_measurement": False},
            "limitations": [
                "This is the first full-sample exploratory radial-profile stack, not a science-approved ICL measurement.",
                "Only 798 targets have a matched control; exact filter mismatch further reduces usable pairs.",
                "The K2.5/40% cut is an exploratory small-mask/area sensitivity branch, not an adopted science threshold.",
                "The source masks are automated, BCG confirmation is incomplete, and residual source wings may bias the profile.",
                "Global-sky coadd/bg does not validate upstream sky transfer or extended PSF-wing transfer.",
                "FLUXMAG0 zero-point calibration is at the image level; Galactic extinction and i/i2 passband response are not corrected.",
                "The CAMIRA association uses projected member positions and counts only, not satellite photometry, stellar masses, or membership weights.",
                "The member-light correlation is descriptive; its nominal tests do not propagate spatial covariance or multiple-bin selection.",
                "No satellite-centered outer-light stack or S19A background-shape lensing stack is made by this run.",
            ],
            "source_sha256": source_hashes,
            "output_sha256": {str(path): sha256(path) for path in
                              (work / "cluster_profiles.tsv", work / "stack_profiles.tsv",
                               work / "pair_metrics_k2p5_c40_r50_300.tsv", work / "cluster_profiles.npz",
                               work / "stack_covariances.npz")},
            "script_sha256": {str(Path(__file__).resolve()): sha256(__file__)},
            "elapsed_seconds": time.monotonic() - start,
        }
        write_json(work / "summary.json", safe(summary))
        write_json(progress_file, {"status": "complete", "completed": len(profile_rows),
                                  "expected": len(frames), "elapsed_seconds": time.monotonic() - start})
        print(json.dumps({"status": "complete", "n_targets": len(profile_rows),
                          "n_pairs": summary["n_matched_controls"],
                          "candidate_pairs": len(candidate_metrics), "output": str(work)}, ensure_ascii=False), flush=True)
    except Exception as exc:
        write_json(error_file, {"status": "failed", "exception": type(exc).__name__, "error": str(exc),
                                "completed_profiles": len(profile_rows),
                                "elapsed_seconds": time.monotonic() - start})
        raise


if __name__ == "__main__":
    main()
