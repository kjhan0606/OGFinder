#!/usr/bin/env python3
"""Propagate measured ±8-pixel covariance kernels through supported ICL annuli.

This is a diagnostic extension of analyze_icl_annulus_aperture_covariance.py.
It contracts the new stationary 33x33 per-frame residual kernel with each
actual valid 200--300 proper-kpc annulus mask and compares the resulting
independent-frame pixel-noise floor with the saved spatial-block interval.
No HSC data are requested and no profile flux is altered.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import time

import numpy as np
from astropy.io import fits
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u

import estimate_hsc_coadd_noise_correlation as noise
from audit_icl_full_sample_coadd_covariance import (
    DEFAULT_SAMPLE, build_records, read_json, read_tsv, sha256,
)
from analyze_icl_annulus_aperture_covariance import (
    load_profiles, SAMPLE, PROFILE, INNER_KPC, OUTER_KPC, BIN_NAME,
)
from download_pilot import validate_fits
from validate_control_images import nJy_factor

INPUT = SAMPLE / "icl_supported_covariance_to_16px_2026-11-06_01"
DEFAULT_OUTPUT = SAMPLE / "icl_aperture_covariance_200_300_16px_2026-11-06_01"
MAX_LAG = 16


def atomic_json(path: Path, value: dict) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, ensure_ascii=False, indent=2,
                                    sort_keys=True, allow_nan=False) + "\n",
                         encoding="utf-8")
    temporary.replace(path)


def load_long_lag_kernels() -> tuple[dict, dict]:
    summary_path = INPUT / "summary.json"
    progress_path = INPUT / "progress.json"
    lag_path = INPUT / "lag_correlation.tsv"
    frame_path = INPUT / "frame_noise_correlation.tsv"
    summary = read_json(summary_path)
    progress = read_json(progress_path)
    if (summary.get("status") != "complete_exploratory_diagnostic"
            or summary.get("science_measurement") is not False
            or summary.get("sample", {}).get("supported_200_300_kpc_pairs") != 224
            or summary.get("method", {}).get("max_lag_pixels") != MAX_LAG
            or progress.get("status") != "complete_exploratory_diagnostic"
            or progress.get("summary_sha256") != sha256(summary_path)
            or progress.get("output_sha256", {}).get("lag_correlation.tsv") != sha256(lag_path)
            or progress.get("output_sha256", {}).get("frame_noise_correlation.tsv") != sha256(frame_path)):
        raise ValueError("long-lag covariance output is incomplete or has changed")
    expected = {(dy, dx) for dy in range(-MAX_LAG, MAX_LAG + 1)
                for dx in range(-MAX_LAG, MAX_LAG + 1)}
    raw: dict[tuple[str, str, str], dict[tuple[int, int], float]] = {}
    for row in read_tsv(lag_path):
        key = (str(row["cluster_id"]), str(row["role"]), str(row["filter"]))
        lag = (int(row["dy"]), int(row["dx"]))
        values = raw.setdefault(key, {})
        if lag in values:
            raise ValueError(f"duplicate lag row for {key}/{lag}")
        value = float(row["median_variance_normalized_covariance"])
        if not math.isfinite(value):
            raise ValueError(f"non-finite lag value for {key}/{lag}")
        values[lag] = value
    kernels = {}
    for key, values in raw.items():
        if set(values) != expected:
            raise ValueError(f"incomplete 33x33 kernel for {key}")
        kernel = np.empty((2 * MAX_LAG + 1, 2 * MAX_LAG + 1), dtype=np.float64)
        for dy, dx in expected:
            kernel[dy + MAX_LAG, dx + MAX_LAG] = 0.5 * (
                values[(dy, dx)] + values[(-dy, -dx)])
        if not np.isfinite(kernel).all() or kernel[MAX_LAG, MAX_LAG] <= 0:
            raise ValueError(f"invalid symmetrized kernel for {key}")
        kernels[key] = kernel
    return kernels, {"summary_path": summary_path, "lag_path": lag_path,
                     "frame_path": frame_path, "summary": summary}


def aperture_measurement(record: dict, frame: dict, profile: dict,
                         kernel: np.ndarray) -> dict:
    image_path = Path(record["path"])
    if sha256(image_path) != record["source_image_sha256"]:
        raise ValueError(f"source image changed for {record['cluster_id']}/{record['role']}")
    info = validate_fits(image_path, frame)
    if info["coadd_header_filter"] != record["filter"]:
        raise ValueError(f"FILTER mismatch for {record['cluster_id']}/{record['role']}")
    with fits.open(image_path, memmap=True, checksum=True) as hdus:
        image = np.asarray(hdus[info["science_hdu"]].data, dtype=np.float32)
        variance = np.asarray(hdus[info["variance_hdu"]].data, dtype=np.float32)
        image_mask = np.asarray(hdus[info["mask_hdu"]].data, dtype=np.int64)
        if record["role"] == "target":
            mask_path = Path(record["target_mask_path"])
            if sha256(mask_path) != record["target_mask_sha256"]:
                raise ValueError(f"target mask changed for {record['cluster_id']}")
            with fits.open(mask_path, memmap=True, checksum=True) as mask_hdus:
                excluded = np.asarray(mask_hdus["MASK_K2P5"].data, dtype=bool)
        else:
            excluded = noise._matched_control_source_mask(
                image, variance, image_mask, float(record["pixel_scale_arcsec"]),
            )
    if image.shape != variance.shape or image.shape != excluded.shape:
        raise ValueError(f"image/variance/mask shape mismatch for {record['cluster_id']}")
    valid = np.isfinite(image) & np.isfinite(variance) & (variance > 0) & ~excluded
    redshift = float(frame["z"])
    kpc_per_pixel = float(FlatLambdaCDM(H0=70., Om0=0.3, Tcmb0=0*u.K)
                          .kpc_proper_per_arcmin(redshift).value
                          * info["pixel_scale_arcsec"] / 60.)
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
        raise ValueError(f"empty annulus for {record['cluster_id']}/{record['role']}")
    area = math.pi * (OUTER_KPC**2-INNER_KPC**2) / kpc_per_pixel**2
    coverage = min(1., n_pixel/area)
    coverage_key = "target_coverage" if record["role"] == "target" else "control_coverage"
    if not math.isclose(coverage, float(profile[coverage_key]), rel_tol=0, abs_tol=2e-6):
        raise ValueError(f"saved profile coverage mismatch for {record['cluster_id']}/{record['role']}")
    conversion = float(nJy_factor(info))
    mean = float(np.mean(image[y0:y1, x0:x1][aperture], dtype=np.float64) * conversion)
    mean_key = "target_sb_njy_arcsec2" if record["role"] == "target" else "control_sb_njy_arcsec2"
    if not math.isclose(mean, float(profile[mean_key]), rel_tol=0, abs_tol=1e-5):
        raise ValueError(f"saved profile mean mismatch for {record['cluster_id']}/{record['role']}")

    local_var = variance[y0:y1, x0:x1].astype(np.float64, copy=False)
    weighted_aperture = np.zeros(local_var.shape, dtype=np.float64)
    weighted_aperture[aperture] = np.sqrt(local_var[aperture])
    diagonal_sum = float(np.sum(local_var[aperture], dtype=np.float64))
    covariance_sum = 0.0
    height, width = weighted_aperture.shape
    for dy in range(-MAX_LAG, MAX_LAG + 1):
        ay0, ay1 = max(0, -dy), min(height, height-dy)
        by0, by1 = ay0+dy, ay1+dy
        for dx in range(-MAX_LAG, MAX_LAG + 1):
            ax0, ax1 = max(0, -dx), min(width, width-dx)
            bx0, bx1 = ax0+dx, ax1+dx
            a = weighted_aperture[ay0:ay1, ax0:ax1]
            b = weighted_aperture[by0:by1, bx0:bx1]
            covariance_sum += kernel[dy+MAX_LAG, dx+MAX_LAG] * float(np.sum(a*b, dtype=np.float64))
    diagonal_variance = diagonal_sum * conversion**2 / n_pixel**2
    lag_variance = covariance_sum * conversion**2 / n_pixel**2
    if not math.isfinite(diagonal_variance) or diagonal_variance <= 0:
        raise ValueError("non-positive diagonal aperture variance")
    if not math.isfinite(lag_variance) or lag_variance <= 0:
        raise ValueError("non-positive 33x33 aperture variance")
    return {
        "cluster_id": str(record["cluster_id"]), "filter": str(record["filter"]),
        "role": str(record["role"]), "n_valid_aperture_pixels": n_pixel,
        "coverage_reconstructed": coverage,
        "surface_brightness_reconstructed_njy_arcsec2": mean,
        "diagonal_aperture_variance_njy2_arcsec4": diagonal_variance,
        "finite_33x33_lag_aperture_variance_njy2_arcsec4": lag_variance,
        "finite_33x33_lag_sigma_inflation": math.sqrt(lag_variance/diagonal_variance),
        "kernel_center": float(kernel[MAX_LAG, MAX_LAG]),
        "source_image_sha256": record["source_image_sha256"],
        "target_mask_sha256": record["target_mask_sha256"] if record["role"] == "target" else None,
    }


def run(output: Path, *, resume: bool = False) -> dict:
    output = Path(output).resolve()
    profiles, profile_info = load_profiles()
    kernels, kernel_info = load_long_lag_kernels()
    records, setup = build_records(DEFAULT_SAMPLE.resolve())
    frame_by_key = {(str(r["cluster_id"]), str(r["role"]), str(r["filter"])): r for r in records}
    sample_doc = read_json(SAMPLE / "sample_manifest.json")
    control_doc = read_json(SAMPLE / "control_download_01/download_manifest.json")
    targets = {str(r["id"]): r for r in sample_doc["clusters"]}
    controls = {str(r["target_cluster_id"]): r for r in control_doc["clusters"]}
    if len(profiles) != 224 or len(frame_by_key) != 1566 or len(targets) != 800 or len(controls) != 798:
        raise ValueError("saved ICL profile or acquisition sample is incomplete")
    expected_keys = {(ident, role, band) for ident, band in profiles
                     for role in ("target", "matched_control")}
    if not expected_keys <= set(kernels) or not expected_keys <= set(frame_by_key):
        raise ValueError("33x33 lag kernels or image records do not cover all supported frames")

    input_paths = [*profile_info["paths"], kernel_info["summary_path"], kernel_info["lag_path"],
                   kernel_info["frame_path"], SAMPLE/"sample_manifest.json",
                   SAMPLE/"control_download_01/download_manifest.json",
                   SAMPLE/"prestack_cluster_masks_2026-10-16_01/records.tsv"]
    input_hashes = {str(path.resolve()): sha256(path) for path in input_paths}
    manifest = {"schema": "astrafek-icl-33x33-aperture-covariance-v1",
                "n_supported_pairs": len(profiles), "max_lag_pixels": MAX_LAG,
                "input_sha256": input_hashes,
                "analysis_script_sha256": sha256(Path(__file__).resolve()),
                "science_measurement": False}
    records_dir = output / "pair_records"
    if output.exists():
        if not resume:
            raise FileExistsError(f"refusing to overwrite {output}; pass --resume for matching checkpoints")
        if not (output / "run_manifest.json").is_file() or read_json(output / "run_manifest.json") != manifest:
            raise ValueError("resume manifest differs from current inputs or runner")
        if not records_dir.is_dir() or (output / "summary.json").is_file():
            raise ValueError("output is not an incomplete, resumable run")
    else:
        if resume:
            raise FileNotFoundError("--resume requested but output does not exist")
        output.mkdir(parents=True, exist_ok=False)
        records_dir.mkdir()
        atomic_json(output / "run_manifest.json", manifest)

    pair_records = {}
    for path in records_dir.glob("pair_*.json"):
        item = read_json(path)
        key = (str(item.get("cluster_id")), str(item.get("filter")))
        if (key not in profiles or item.get("schema") != "astrafek-icl-33x33-aperture-pair-v1"
                or item.get("science_measurement") is not False or key in pair_records
                or len(item.get("frame_rows", [])) != 2):
            raise ValueError(f"invalid saved pair checkpoint: {path}")
        pair_records[key] = item

    target_frames = {str(r["id"]): r for r in sample_doc["clusters"]}
    control_by_target = {str(r["target_cluster_id"]): r for r in control_doc["clusters"]}
    started = time.time()
    for index, (ident, band) in enumerate(sorted(profiles), 1):
        key = (ident, band)
        if key in pair_records:
            continue
        profile = profiles[key]
        target_source = target_frames.get(ident)
        control_source = control_by_target.get(ident)
        if target_source is None or control_source is None:
            raise ValueError(f"missing matched source frame metadata for {ident}")
        frame_rows = []
        for role, source_frame in (("target", target_source), ("matched_control", control_source)):
            frame = frame_by_key[(ident, role, band)]
            measurement = aperture_measurement(frame, source_frame, profile,
                                               kernels[(ident, role, band)])
            frame_rows.append(measurement)
        row = {
            "cluster_id": ident, "filter": band,
            "target_diagonal_variance": frame_rows[0]["diagonal_aperture_variance_njy2_arcsec4"],
            "control_diagonal_variance": frame_rows[1]["diagonal_aperture_variance_njy2_arcsec4"],
            "target_finite_33x33_lag_variance": frame_rows[0]["finite_33x33_lag_aperture_variance_njy2_arcsec4"],
            "control_finite_33x33_lag_variance": frame_rows[1]["finite_33x33_lag_aperture_variance_njy2_arcsec4"],
            "target_minus_control_delta_njy_arcsec2": float(profile["delta_sb_njy_arcsec2"]),
        }
        item = {"schema": "astrafek-icl-33x33-aperture-pair-v1", "cluster_id": ident,
                "filter": band, "frame_rows": frame_rows, "pair_row": row,
                "science_measurement": False}
        atomic_json(records_dir / f"pair_{ident}_{band}.json", item)
        pair_records[key] = item
        atomic_json(output / "progress.json", {
            "status": "running", "pairs_done": len(pair_records),
            "pairs_expected": len(profiles), "last_pair": [ident, band],
            "science_measurement": False,
        })
        print(f"[{index}/{len(profiles)}] 33x33 ICL aperture covariance {ident} {band}", flush=True)

    frame_rows = [row for item in pair_records.values() for row in item["frame_rows"]]
    pair_rows = [item["pair_row"] for item in pair_records.values()]
    if len(pair_rows) != 224 or len(frame_rows) != 448:
        raise ValueError("not all supported pairs have two aperture variance records")
    frame_path, pair_path = output/"frame_aperture_covariance.tsv", output/"pair_aperture_covariance.tsv"
    for path, rows in ((frame_path, frame_rows), (pair_path, pair_rows)):
        with path.open("x", encoding="utf-8", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t",
                                    extrasaction="raise", lineterminator="\n")
            writer.writeheader(); writer.writerows(rows)
    old_summary = read_json(SAMPLE/"icl_aperture_covariance_200_300_2026-11-06_01/summary.json")
    old_pairs = read_tsv(SAMPLE/"icl_aperture_covariance_200_300_2026-11-06_01/pair_aperture_covariance.tsv")
    old_by_key = {(str(r["cluster_id"]), str(r["filter"])): r for r in old_pairs}
    if len(old_by_key) != 224:
        raise ValueError("existing 17x17 exact-aperture comparison table is incomplete")
    stacks = profile_info["stack"]
    by_filter = {}
    for band in ("i", "i2"):
        selected = [r for r in pair_rows if r["filter"] == band]
        n = len(selected)
        diag_var = sum(r["target_diagonal_variance"]+r["control_diagonal_variance"] for r in selected)/n**2
        lag_var = sum(r["target_finite_33x33_lag_variance"]+r["control_finite_33x33_lag_variance"] for r in selected)/n**2
        old_diag = sum(float(old_by_key[(r["cluster_id"], band)]["target_diagonal_variance"])
                       + float(old_by_key[(r["cluster_id"], band)]["control_diagonal_variance"])
                       for r in selected)/n**2
        old_lag = sum(float(old_by_key[(r["cluster_id"], band)]["target_finite_17x17_lag_variance"])
                      + float(old_by_key[(r["cluster_id"], band)]["control_finite_17x17_lag_variance"])
                      for r in selected)/n**2
        stack = stacks[(band, BIN_NAME)]
        ci = [float(stack["spatial_block_bootstrap_95_low"]),
              float(stack["spatial_block_bootstrap_95_high"])]
        halfwidth = (ci[1]-ci[0])/2
        by_filter[band] = {
            "n_pairs": n,
            "target_median_sigma_inflation": float(np.median([r["finite_33x33_lag_sigma_inflation"] for r in frame_rows if r["filter"] == band and r["role"] == "target"])),
            "control_median_sigma_inflation": float(np.median([r["finite_33x33_lag_sigma_inflation"] for r in frame_rows if r["filter"] == band and r["role"] == "matched_control"])),
            "diagonal_independent_frame_noise_stack_sigma_njy_arcsec2": math.sqrt(diag_var),
            "finite_33x33_lag_independent_frame_noise_stack_sigma_njy_arcsec2": math.sqrt(lag_var),
            "finite_33x33_to_diagonal_stack_sigma_ratio": math.sqrt(lag_var/diag_var),
            "previous_finite_17x17_lag_stack_sigma_njy_arcsec2": math.sqrt(old_lag),
            "33x33_to_previous_17x17_stack_sigma_ratio": math.sqrt(lag_var/old_lag),
            "existing_spatial_block_bootstrap_95_interval_njy_arcsec2": ci,
            "33x33_pixel_noise_sigma_to_spatial_interval_halfwidth": math.sqrt(lag_var)/halfwidth if halfwidth > 0 else None,
            "stack_mean_pair_delta_njy_arcsec2": float(stack["mean_pair_delta_njy_arcsec2"]),
        }
    summary = {
        "schema": "astrafek-icl-200-300-33x33-aperture-covariance-v1",
        "status": "complete_exploratory_diagnostic", "science_measurement": False,
        "n_supported_pairs": len(pair_rows), "n_frames": len(frame_rows),
        "n_pairs_by_filter": {band: sum(r["filter"] == band for r in pair_rows) for band in ("i", "i2")},
        "method": "contract each saved per-frame symmetric 33x33 variance-normalized lag kernel with the actual valid 200-300 kpc annulus mean weights; propagate target+control variances to equal-weight supported stack assuming independent frames",
        "filters": by_filter, "input_sha256": input_hashes,
        "long_lag_kernel_summary_sha256": sha256(kernel_info["summary_path"]),
        "previous_17x17_summary_sha256": sha256(SAMPLE/"icl_aperture_covariance_200_300_2026-11-06_01/summary.json"),
        "hsc_network_used": False, "data_requested": False,
        "images_modified": False, "covariance_correction_applied": False,
        "limitations": [
            "The 33x33 kernels are patch medians and assume local stationarity throughout the aperture.",
            "This finite-lag pixel-noise floor omits longer-range covariance, source Poisson, sky transfer, PSF wings, mask uncertainty, and cross-frame spatial covariance.",
            "Target and control terms are propagated as independent frames; the saved spatial-block interval remains the broader empirical uncertainty reference.",
            "No covariance correction is applied to the measured ICL surface-brightness contrast.",
        ],
        "analysis_script_sha256": sha256(Path(__file__).resolve()),
        "elapsed_seconds": time.time()-started,
    }
    summary_path = output/"summary.json"
    atomic_json(summary_path, summary)
    hashes = {"summary_sha256": sha256(summary_path),
              "frame_aperture_covariance_sha256": sha256(frame_path),
              "pair_aperture_covariance_sha256": sha256(pair_path)}
    atomic_json(output/"output_hashes.json", hashes)
    atomic_json(output/"progress.json", {"status": summary["status"],
                "pairs_done": len(pair_rows), "pairs_expected": len(profiles),
                "science_measurement": False, "output_hashes": hashes})
    print(json.dumps({"status": summary["status"], "filters": by_filter,
                      "output": str(output), "output_hashes": hashes,
                      "science_measurement": False}, ensure_ascii=False, indent=2), flush=True)
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    run(args.output, resume=args.resume)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
