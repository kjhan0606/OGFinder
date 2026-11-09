#!/usr/bin/env python3
"""Propagate measured local coadd lag kernels through the full 200--300 kpc ICL apertures.

For every annulus-supported same-filter target/control pair, this contracts the
existing per-frame 9x9 residual-covariance kernel with the actual saved K=2.5
target mask or the established matched-control source mask. It validates the
reconstructed aperture mean and coverage against the saved ICL profile table.
The result is an approximate finite-lag pixel-noise floor, not total stack
uncertainty or a correction to the measured ICL profile.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import time

import numpy as np
from astropy.io import fits
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u

import estimate_hsc_coadd_noise_correlation as noise
from audit_icl_full_sample_coadd_covariance import build_records, read_json, read_tsv, sha256
from download_pilot import validate_fits
from validate_control_images import nJy_factor


SAMPLE = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
PROFILE = SAMPLE / "icl_annulus_support_2026-10-18_01"
COVARIANCE = SAMPLE / "icl_coadd_noise_correlation_full_2026-11-03_01"
DEFAULT_OUTPUT = SAMPLE / "icl_aperture_covariance_200_300_2026-11-03_01"
BIN_NAME = "200-300"
INNER_KPC, OUTER_KPC = 200.0, 300.0
MAX_LAG = 4


def atomic_json(path: Path, value: dict) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, ensure_ascii=False, indent=2,
                                    sort_keys=True, allow_nan=False) + "\n",
                         encoding="utf-8")
    temporary.replace(path)


def load_profiles() -> tuple[dict[tuple[str, str], dict], dict]:
    summary_path = PROFILE / "summary.json"
    pair_path = PROFILE / "annulus_pair_profiles.tsv"
    stack_path = PROFILE / "annulus_stack_summary.tsv"
    summary = read_json(summary_path)
    if (summary.get("status") != "complete" or summary.get("science_gates", {}).get("science_measurement") is not False
            or summary.get("n_200_300_pairs") != 224):
        raise ValueError("200--300 kpc ICL support profile is incomplete or outside its exploratory gate")
    for path in (pair_path, stack_path):
        expected = summary.get("output_sha256", {}).get(str(path))
        if not expected or expected != sha256(path):
            raise ValueError(f"saved ICL profile checksum mismatch: {path}")
    rows = [r for r in read_tsv(pair_path) if r["radial_bin_kpc"] == BIN_NAME]
    if len(rows) != 224:
        raise ValueError(f"expected 224 supported annulus pairs, found {len(rows)}")
    profiles = {}
    for row in rows:
        key = (str(row["id"]), str(row["filter"]))
        if key in profiles:
            raise ValueError(f"duplicate saved ICL pair: {key}")
        profiles[key] = row
    stack = {(r["filter"], r["radial_bin_kpc"]): r for r in read_tsv(stack_path)}
    return profiles, {"summary": summary, "stack": stack,
                      "paths": (summary_path, pair_path, stack_path)}


def load_kernels() -> tuple[dict[tuple[str, str, str], np.ndarray], dict]:
    summary_path = COVARIANCE / "summary.json"
    audit_path = COVARIANCE / "independent_audit.json"
    frame_path = COVARIANCE / "frame_noise_correlation.tsv"
    lag_path = COVARIANCE / "lag_correlation.tsv"
    summary, audit = read_json(summary_path), read_json(audit_path)
    if (summary.get("status") != "complete_exploratory"
            or summary.get("science_measurement") is not False
            or summary.get("sample", {}).get("same_filter_pairs") != 783
            or audit.get("status") != "verified"
            or audit.get("science_measurement") is not False
            or audit.get("frame_and_lag_arithmetic_reproduced") is not True
            or audit.get("input_acquisition_hashes_matched") is not True
            or audit.get("n_frame_rows") != 1566
            or audit.get("n_lag_rows") != 1566 * 81
            or audit.get("frame_table_sha256") != sha256(frame_path)
            or audit.get("lag_table_sha256") != sha256(lag_path)
            or audit.get("summary_sha256") != sha256(summary_path)):
        raise ValueError("coadd-lag input lacks complete independent provenance")
    if (summary.get("output_sha256", {}).get(str(frame_path)) != sha256(frame_path)
            or summary.get("output_sha256", {}).get(str(lag_path)) != sha256(lag_path)):
        raise ValueError("coadd-lag output hashes differ from the saved summary")
    frame_rows = read_tsv(frame_path)
    frames = {(r["cluster_id"], r["role"], r["filter"]): r for r in frame_rows}
    if len(frames) != 1566:
        raise ValueError("duplicate frame key in the 1,566-row coadd-covariance table")
    raw: dict[tuple[str, str, str], dict[tuple[int, int], float]] = {}
    for row in read_tsv(lag_path):
        key = (row["cluster_id"], row["role"], row["filter"])
        lag = (int(row["dy"]), int(row["dx"]))
        value = float(row["median_variance_normalized_covariance"])
        if not math.isfinite(value) or lag in raw.setdefault(key, {}):
            raise ValueError(f"invalid or duplicate covariance lag {key}/{lag}")
        raw[key][lag] = value
    if set(raw) != set(frames):
        raise ValueError("frame and lag-covariance key sets differ")
    expected_lags = {(dy, dx) for dy in range(-MAX_LAG, MAX_LAG + 1)
                     for dx in range(-MAX_LAG, MAX_LAG + 1)}
    kernels = {}
    for key, values in raw.items():
        if set(values) != expected_lags:
            raise ValueError(f"incomplete 9x9 lag kernel for {key}")
        kernel = np.empty((9, 9), dtype=np.float64)
        for dy, dx in expected_lags:
            kernel[dy + MAX_LAG, dx + MAX_LAG] = 0.5 * (
                values[(dy, dx)] + values[(-dy, -dx)])
        if not np.all(np.isfinite(kernel)) or kernel[4, 4] <= 0:
            raise ValueError(f"invalid symmetrized covariance kernel for {key}")
        kernels[key] = kernel
    return kernels, {"summary_path": summary_path, "audit_path": audit_path,
                     "frame_path": frame_path, "lag_path": lag_path,
                     "frames": frames}


def aperture_measurement(record: dict, frame: dict, profile: dict,
                         kernel: np.ndarray) -> dict:
    image_path = Path(record["path"])
    if image_path.stat().st_size <= 0:
        raise ValueError(f"empty input image: {image_path}")
    info = validate_fits(image_path, frame)
    if info["coadd_header_filter"] != record["filter"]:
        raise ValueError(f"FILTER changed for {record['cluster_id']}/{record['role']}")
    with fits.open(image_path, memmap=True, checksum=True) as hdus:
        image = np.asarray(hdus[info["science_hdu"]].data, dtype=np.float32)
        variance = np.asarray(hdus[info["variance_hdu"]].data, dtype=np.float32)
        image_mask = np.asarray(hdus[info["mask_hdu"]].data, dtype=np.int64)
        if record["role"] == "target":
            mask_path = Path(record["target_mask_path"])
            if sha256(mask_path) != record["target_mask_sha256"]:
                raise ValueError(f"saved K=2.5 mask hash changed for {record['cluster_id']}")
            with fits.open(mask_path, memmap=True, checksum=True) as mask_hdus:
                excluded = np.asarray(mask_hdus["MASK_K2P5"].data, dtype=bool)
        else:
            excluded = noise._matched_control_source_mask(
                image, variance, image_mask, float(record["pixel_scale_arcsec"]),
            )
    if image.shape != variance.shape or image.shape != excluded.shape:
        raise ValueError(f"image/variance/mask grid mismatch for {record['cluster_id']}/{record['role']}")
    valid = np.isfinite(image) & np.isfinite(variance) & (variance > 0) & ~excluded
    z = float(frame["z"])
    kpc_per_pixel = float(FlatLambdaCDM(H0=70., Om0=0.3, Tcmb0=0*u.K)
                          .kpc_proper_per_arcmin(z).value * info["pixel_scale_arcsec"] / 60.)
    center_x, center_y = map(float, info["center_pixel_xy"])
    half = int(math.ceil(OUTER_KPC / kpc_per_pixel)) + MAX_LAG + 2
    x0, x1 = max(0, int(math.floor(center_x-half))), min(image.shape[1], int(math.ceil(center_x+half))+1)
    y0, y1 = max(0, int(math.floor(center_y-half))), min(image.shape[0], int(math.ceil(center_y+half))+1)
    yy, xx = np.ogrid[y0:y1, x0:x1]
    radius = np.hypot((xx-center_x)*kpc_per_pixel, (yy-center_y)*kpc_per_pixel)
    aperture = ((radius >= INNER_KPC) & (radius < OUTER_KPC)
                & valid[y0:y1, x0:x1])
    n_pixel = int(aperture.sum())
    if n_pixel < 1:
        raise ValueError(f"empty ICL aperture for {record['cluster_id']}/{record['role']}")
    area = math.pi * (OUTER_KPC**2 - INNER_KPC**2) / kpc_per_pixel**2
    coverage = min(1., n_pixel / area)
    coverage_key = "target_coverage" if record["role"] == "target" else "control_coverage"
    expected_coverage = float(profile[coverage_key])
    if not math.isclose(coverage, expected_coverage, rel_tol=0, abs_tol=2e-6):
        raise ValueError(f"aperture coverage mismatch for {record['cluster_id']}/{record['role']}: {coverage} vs {expected_coverage}")
    conversion = float(nJy_factor(info))
    mean = float(np.mean(image[y0:y1, x0:x1][aperture], dtype=np.float64) * conversion)
    mean_key = "target_sb_njy_arcsec2" if record["role"] == "target" else "control_sb_njy_arcsec2"
    expected_mean = float(profile[mean_key])
    if not math.isclose(mean, expected_mean, rel_tol=0, abs_tol=1e-5):
        raise ValueError(f"reconstructed annulus mean mismatch for {record['cluster_id']}/{record['role']}: {mean} vs {expected_mean}")

    local_var = variance[y0:y1, x0:x1].astype(np.float64, copy=False)
    field = np.zeros(local_var.shape, dtype=np.float64)
    field[aperture] = np.sqrt(local_var[aperture])
    diagonal_sum = float(np.sum(local_var[aperture], dtype=np.float64))
    covariance_sum = 0.0
    height, width = field.shape
    for dy in range(-MAX_LAG, MAX_LAG+1):
        ya0, ya1 = max(0, -dy), min(height, height-dy)
        yb0, yb1 = ya0+dy, ya1+dy
        for dx in range(-MAX_LAG, MAX_LAG+1):
            xa0, xa1 = max(0, -dx), min(width, width-dx)
            xb0, xb1 = xa0+dx, xa1+dx
            a = field[ya0:ya1, xa0:xa1]
            b = field[yb0:yb1, xb0:xb1]
            covariance_sum += float(kernel[dy+MAX_LAG, dx+MAX_LAG]) * float(
                np.sum(a*b, dtype=np.float64))
    diagonal_variance = diagonal_sum * conversion**2 / n_pixel**2
    lag_variance = covariance_sum * conversion**2 / n_pixel**2
    if not math.isfinite(diagonal_variance) or diagonal_variance <= 0 or not math.isfinite(lag_variance) or lag_variance <= 0:
        raise ValueError(f"non-positive aperture variance for {record['cluster_id']}/{record['role']}")
    return {
        "cluster_id": record["cluster_id"], "filter": record["filter"], "role": record["role"],
        "n_valid_aperture_pixels": n_pixel, "coverage_reconstructed": coverage,
        "surface_brightness_reconstructed_njy_arcsec2": mean,
        "diagonal_aperture_variance_njy2_arcsec4": diagonal_variance,
        "finite_9x9_lag_aperture_variance_njy2_arcsec4": lag_variance,
        "finite_9x9_lag_sigma_inflation": math.sqrt(lag_variance/diagonal_variance),
        "kernel_center": float(kernel[4, 4]),
        "source_image_sha256_from_verified_acquisition": record["source_image_sha256"],
        "target_mask_sha256": record["target_mask_sha256"] if record["role"] == "target" else None,
    }


def run(output: Path, allow_local_full_sample: bool = False, resume: bool = False) -> dict:
    if not os.environ.get("SLURM_JOB_ID") and not allow_local_full_sample:
        raise RuntimeError("use Slurm for the full 224-pair run; --allow-local-full-sample is only for an explicitly managed Lageunha CPU run")
    profiles, profile_info = load_profiles()
    kernels, kernel_info = load_kernels()
    records, setup = build_records(SAMPLE)
    frame_by_key = {(r["cluster_id"], r["role"], r["filter"]): r for r in records}
    sample_doc = read_json(SAMPLE / "sample_manifest.json")
    control_doc = read_json(SAMPLE / "control_download_01/download_manifest.json")
    targets = {str(r["id"]): r for r in sample_doc["clusters"]}
    controls = {str(r["target_cluster_id"]): r for r in control_doc["clusters"]}
    if len(targets) != 800 or len(controls) != 798 or len(frame_by_key) != 1566:
        raise ValueError("full-sample target/control metadata or covariance records are incomplete")
    input_paths = [*profile_info["paths"], kernel_info["summary_path"], kernel_info["audit_path"],
                   kernel_info["frame_path"], kernel_info["lag_path"],
                   SAMPLE/"sample_manifest.json", SAMPLE/"control_download_01/download_manifest.json",
                   SAMPLE/"prestack_cluster_masks_2026-10-16_01/records.tsv"]
    input_hashes = {str(path.resolve()): sha256(path) for path in input_paths}
    script_hash = sha256(Path(__file__).resolve())
    manifest = {"schema": "astrafek-icl-200-300-aperture-covariance-run-v1",
                "n_supported_pairs": len(profiles), "input_sha256": input_hashes,
                "analysis_script_sha256": script_hash, "science_measurement": False}
    records_dir = output / "pair_records"
    if output.exists():
        if not resume:
            raise FileExistsError(f"refusing to overwrite {output}; pass --resume only for a matching partial run")
        manifest_path = output / "run_manifest.json"
        if not manifest_path.is_file() or read_json(manifest_path) != manifest:
            raise ValueError("partial output run manifest does not match current inputs or analysis script")
        if not records_dir.is_dir() or (output / "summary.json").is_file():
            raise ValueError("output is not a resumable partial run")
    else:
        if resume:
            raise FileNotFoundError("--resume requested but output directory does not exist")
        output.mkdir(parents=True, exist_ok=False)
        records_dir.mkdir()
        atomic_json(output / "run_manifest.json", manifest)

    pair_records = {}
    for path in records_dir.glob("pair_*.json"):
        record = read_json(path)
        key = (str(record.get("cluster_id")), str(record.get("filter")))
        if (key not in profiles or record.get("schema") != "astrafek-icl-aperture-covariance-pair-v1"
                or record.get("science_measurement") is not False
                or len(record.get("frame_rows", [])) != 2
                or key in pair_records):
            raise ValueError(f"invalid or duplicate partial pair record: {path}")
        pair_records[key] = record
    if len(pair_records) > len(profiles):
        raise ValueError("partial output contains more pairs than the approved profile sample")

    frame_rows, pair_rows = [], []
    for key, record in pair_records.items():
        frame_rows.extend(record["frame_rows"])
        pair_rows.append(record["pair_row"])
    expected = len(profiles) * 2
    progress = output / "progress.json"
    atomic_json(progress, {"status": "running", "frames_done": len(frame_rows),
                           "frames_expected": expected, "pairs_done": len(pair_rows),
                           "pairs_expected": len(profiles), "science_measurement": False})
    started = time.time()
    try:
        for index, (key, profile) in enumerate(sorted(profiles.items()), 1):
            ident, band = key
            if key in pair_records:
                continue
            target_frame, control_frame = targets[ident], controls.get(ident)
            if control_frame is None:
                raise ValueError(f"missing matched control for {ident}")
            for role, source_frame in (("target", target_frame), ("matched_control", control_frame)):
                record = frame_by_key.get((ident, role, band))
                kernel = kernels.get((ident, role, band))
                if record is None or kernel is None:
                    raise ValueError(f"missing frame/kernel input for {ident}/{role}/{band}")
                measurement = aperture_measurement(record, source_frame, profile, kernel)
                frame_rows.append(measurement)
            trow, crow = frame_rows[-2], frame_rows[-1]
            pair_rows.append({
                "cluster_id": ident, "filter": band,
                "target_n_valid_pixels": trow["n_valid_aperture_pixels"],
                "control_n_valid_pixels": crow["n_valid_aperture_pixels"],
                "target_diagonal_variance": trow["diagonal_aperture_variance_njy2_arcsec4"],
                "control_diagonal_variance": crow["diagonal_aperture_variance_njy2_arcsec4"],
                "target_finite_9x9_lag_variance": trow["finite_9x9_lag_aperture_variance_njy2_arcsec4"],
                "control_finite_9x9_lag_variance": crow["finite_9x9_lag_aperture_variance_njy2_arcsec4"],
                "target_minus_control_delta_njy_arcsec2": float(profile["delta_sb_njy_arcsec2"]),
            })
            pair_record = {"schema": "astrafek-icl-aperture-covariance-pair-v1",
                           "cluster_id": ident, "filter": band,
                           "frame_rows": frame_rows[-2:], "pair_row": pair_rows[-1],
                           "science_measurement": False}
            atomic_json(records_dir / f"pair_{ident}_{band}.json", pair_record)
            atomic_json(progress, {"status": "running", "frames_done": len(frame_rows),
                                   "frames_expected": expected, "pairs_done": len(pair_rows),
                                   "last_pair": [ident, band], "science_measurement": False})
            print(f"[{index}/{len(profiles)}] aperture covariance {ident} {band}", flush=True)

        frame_path, pair_path = output / "frame_aperture_covariance.tsv", output / "pair_aperture_covariance.tsv"
        for path, rows in ((frame_path, frame_rows), (pair_path, pair_rows)):
            with path.open("x", encoding="utf-8", newline="") as stream:
                writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t", extrasaction="raise")
                writer.writeheader(); writer.writerows(rows)
        stack_rows = profile_info["stack"]
        by_filter = {}
        for band in ("i", "i2"):
            selected = [r for r in pair_rows if r["filter"] == band]
            n = len(selected)
            if n == 0:
                raise ValueError(f"no supported annulus pairs for {band}")
            diag_stack_var = sum(r["target_diagonal_variance"]+r["control_diagonal_variance"] for r in selected)/n**2
            corr_stack_var = sum(r["target_finite_9x9_lag_variance"]+r["control_finite_9x9_lag_variance"] for r in selected)/n**2
            stack = stack_rows.get((band, BIN_NAME))
            if stack is None or int(stack["n_pairs"]) != n:
                raise ValueError(f"saved stack support disagrees with aperture audit for {band}")
            ci = [float(stack["spatial_block_bootstrap_95_low"]),
                  float(stack["spatial_block_bootstrap_95_high"])]
            halfwidth = (ci[1]-ci[0])/2
            by_filter[band] = {
                "n_pairs": n,
                "median_target_sigma_inflation": float(np.median([r["finite_9x9_lag_sigma_inflation"] for r in frame_rows if r["filter"] == band and r["role"] == "target"])),
                "median_control_sigma_inflation": float(np.median([r["finite_9x9_lag_sigma_inflation"] for r in frame_rows if r["filter"] == band and r["role"] == "matched_control"])),
                "diagonal_independent_frame_noise_stack_sigma_njy_arcsec2": math.sqrt(diag_stack_var),
                "finite_9x9_lag_independent_frame_noise_stack_sigma_njy_arcsec2": math.sqrt(corr_stack_var),
                "finite_9x9_lag_to_diagonal_stack_sigma_ratio": math.sqrt(corr_stack_var/diag_stack_var),
                "existing_spatial_block_bootstrap_95_interval_njy_arcsec2": ci,
                "existing_spatial_interval_halfwidth_njy_arcsec2": halfwidth,
                "finite_9x9_lag_noise_sigma_to_spatial_interval_halfwidth": math.sqrt(corr_stack_var)/halfwidth if halfwidth > 0 else None,
                "stack_mean_pair_delta_njy_arcsec2": float(stack["mean_pair_delta_njy_arcsec2"]),
            }
        summary = {
            "schema": "astrafek-icl-200-300-aperture-covariance-v1",
            "status": "complete_exploratory_diagnostic",
            "science_measurement": False,
            "n_supported_pairs": len(pair_rows), "n_frames": len(frame_rows),
            "n_pairs_by_filter": {band: sum(r["filter"] == band for r in pair_rows) for band in ("i", "i2")},
            "method": "contract each saved per-frame symmetric 9x9 variance-normalized lag-covariance kernel with the actual valid 200-300 kpc annulus mean weights; propagate target+control per-pair pixel variances to the equal-weight supported stack assuming independent frames",
            "filters": by_filter,
            "limitations": [
                "The finite 9x9 kernels are medians over sampled image patches and assume local stationarity across each aperture.",
                "The calculation omits covariance beyond 4 pixels, upstream sky transfer, source Poisson noise, PSF wings, mask uncertainty, and cross-frame spatial covariance.",
                "Target and control frame pixel-noise terms are combined as independent; larger-scale field covariance remains in the spatial block bootstrap and is not replaced here.",
                "The value is an approximate pixel-noise floor for the uncorrected aperture mean, not a correction to the observed ICL profile or total uncertainty.",
            ],
            "inputs": input_hashes,
            "analysis_script_sha256": script_hash,
            "coadd_kernel_independent_audit_sha256": sha256(kernel_info["audit_path"]),
            "hsc_network_used": False, "data_requested": False,
            "images_modified": False, "covariance_correction_applied": False,
            "elapsed_seconds": time.time()-started,
        }
        summary_path = output / "summary.json"
        atomic_json(summary_path, summary)
        output_hashes = {"summary_sha256": sha256(summary_path),
                         "frame_aperture_covariance_sha256": sha256(frame_path),
                         "pair_aperture_covariance_sha256": sha256(pair_path)}
        atomic_json(output / "output_hashes.json", output_hashes)
        atomic_json(progress, {"status": summary["status"], "frames_done": len(frame_rows),
                               "frames_expected": expected, "pairs_done": len(pair_rows),
                               "science_measurement": False, "output_hashes": output_hashes})
        print(json.dumps({"status": summary["status"], "pairs": summary["n_pairs_by_filter"],
                          "filters": by_filter, "output_hashes": output_hashes,
                          "science_measurement": False}, ensure_ascii=False, indent=2), flush=True)
        return summary
    except BaseException as exc:
        atomic_json(progress, {"status": "interrupted" if isinstance(exc, KeyboardInterrupt) else "stopped_on_error",
                               "frames_done": len(frame_rows), "frames_expected": expected,
                               "pairs_done": len(pair_rows), "error_type": type(exc).__name__,
                               "science_measurement": False})
        raise


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--allow-local-full-sample", action="store_true",
                        help="run a single CPU on Lageunha when Slurm is unavailable")
    parser.add_argument("--resume", action="store_true",
                        help="resume a matching incomplete run from hash-checked per-pair records")
    args = parser.parse_args()
    run(args.output.resolve(), args.allow_local_full_sample, args.resume)


if __name__ == "__main__":
    main()
