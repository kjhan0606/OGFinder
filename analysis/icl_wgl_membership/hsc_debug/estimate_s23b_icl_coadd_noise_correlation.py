#!/usr/bin/env python3
"""Measure local HSC coadd residual covariance for the S23B ICL stack.

This CPU-only diagnostic reads the already approved S23B target and matched
random-control images, their independent FITS/WCS audits, the exact 1,022-pair
stack table, and the saved K=2.5 masks. It estimates variance-normalized
pixel-lag covariance outside a fixed 60-arcsec center exclusion. It does not
modify images, contact HSC, or provide the full covariance of an ICL estimator.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import re
import zlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits

import estimate_hsc_coadd_noise_correlation as noise


DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
TARGET_PLAN = DATA / "s23b_icl_target_image_plan_2026-11-07_01"
CONTROL_PLAN = DATA / "s23b_icl_control_image_plan_2026-11-07_02"
CONTROL_DOWNLOAD = DATA / "s23b_icl_control_image_plan_2026-11-07_01/control_download_01"
TARGET_AUDIT = DATA / "s23b_icl_target_image_audit_2026-11-07_01"
CONTROL_AUDIT = DATA / "s23b_icl_control_image_audit_2026-11-07_01"
MASK_ROOT = DATA / "s23b_icl_mask_support_2026-11-07_01"
STACK = DATA / "s23b_icl_target_control_stack_2026-11-07_01"
DEFAULT_OUTPUT = DATA / "s23b_icl_residual_covariance_2026-11-07_01"
TARGET_MANIFEST_SHA256 = "e487582ba40cd3e58da67704242db4d7219e381fb08f4b33fed47ce4ff268b44"
CONTROL_MANIFEST_SHA256 = "ae9eb406cfc179d9980031bb67f85f8ee71c4580ba59e97713278823db5c06df"
EXPECTED_PAIRS = 1022
EXPECTED_FILTERS = {"i": 241, "i2": 781}
SCHEMA = "astrafek-s23b-icl-coadd-residual-covariance-v1"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def write_json_atomic(path: Path, value: dict) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, ensure_ascii=False, indent=2,
                                    sort_keys=True, allow_nan=False) + "\n",
                         encoding="utf-8")
    temporary.replace(path)


def verify_audit(summary_path: Path, table_path: Path, *, expected_count: int) -> tuple[dict, list[dict[str, str]]]:
    summary = read_json(summary_path)
    rows = read_tsv(table_path)
    if (summary.get("status") != "complete_header_wcs_file_integrity_audit"
            or summary.get("n_expected_frames") != expected_count
            or summary.get("n_verified_frames") != expected_count
            or summary.get("science_measurement") is not False
            or summary.get("network_used_by_audit") is not False):
        raise ValueError(f"incomplete or unexpected image audit: {summary_path}")
    recorded = summary.get("output_sha256", {}).get(table_path.name)
    if not recorded or sha256(table_path) != recorded:
        raise ValueError(f"frame audit table hash mismatch: {table_path}")
    for raw_path, digest in summary.get("inputs_sha256", {}).items():
        path = Path(raw_path)
        if not path.is_file() or sha256(path) != digest:
            raise ValueError(f"image-audit input changed: {path}")
    if len(rows) != expected_count or len({row["id"] for row in rows}) != expected_count:
        raise ValueError(f"image audit table does not contain {expected_count} unique frames")
    return summary, rows


def read_mask_records(role: str) -> tuple[dict, dict[str, dict[str, str]]]:
    directory = MASK_ROOT / role
    summary_path = directory / "summary.json"
    records_path = directory / "records.tsv"
    summary = read_json(summary_path)
    rows = read_tsv(records_path)
    if summary.get("science_measurement") is not False:
        raise ValueError(f"{role} mask summary has an unexpected science gate")
    expected = summary.get("output_sha256", {}).get(str(records_path))
    if expected is None:
        expected = summary.get("output_sha256", {}).get(records_path.name)
    if expected is None or sha256(records_path) != expected:
        raise ValueError(f"{role} mask records checksum mismatch")
    indexed = {str(row["id"]): row for row in rows}
    if len(indexed) != len(rows):
        raise ValueError(f"duplicate IDs in {role} mask records")
    return summary, indexed


def build_records() -> tuple[list[dict], dict]:
    target_manifest_path = TARGET_PLAN / "target_manifest.json"
    control_manifest_path = CONTROL_PLAN / "control_manifest.json"
    target_download = TARGET_PLAN / "target_download_01"
    control_download = CONTROL_DOWNLOAD
    target_status_path = target_download / "download_status.json"
    control_status_path = control_download / "download_status.json"
    target_download_manifest_path = target_download / "download_manifest.json"
    control_download_manifest_path = control_download / "download_manifest.json"
    target_audit_summary_path = TARGET_AUDIT / "summary.json"
    target_audit_table_path = TARGET_AUDIT / "verified_target_frames.tsv"
    control_audit_summary_path = CONTROL_AUDIT / "summary.json"
    control_audit_table_path = CONTROL_AUDIT / "verified_control_frames.tsv"
    stack_summary_path = STACK / "summary.json"
    profile_path = STACK / "paired_radial_profiles.tsv"
    required = (
        target_manifest_path, control_manifest_path, target_status_path,
        control_status_path, target_download_manifest_path,
        control_download_manifest_path, target_audit_summary_path,
        target_audit_table_path, control_audit_summary_path,
        control_audit_table_path, stack_summary_path, profile_path,
    )
    for path in required:
        if not path.is_file():
            raise FileNotFoundError(path)

    target_manifest_sha = sha256(target_manifest_path)
    control_manifest_sha = sha256(control_manifest_path)
    if (target_manifest_sha != TARGET_MANIFEST_SHA256
            or control_manifest_sha != CONTROL_MANIFEST_SHA256):
        raise ValueError("approved S23B image manifest hash changed")
    target_status = read_json(target_status_path)
    control_status = read_json(control_status_path)
    target_downloaded = read_json(target_download_manifest_path)
    control_downloaded = read_json(control_download_manifest_path)
    if (target_status.get("status") != "complete"
            or target_status.get("manifest_sha256") != target_manifest_sha
            or target_status.get("n_expected") != 1043
            or target_status.get("credentials_saved") is not False
            or control_status.get("status") != "complete"
            or control_status.get("manifest_sha256") != control_manifest_sha
            or control_status.get("n_expected") != 1041
            or control_status.get("credentials_saved") is not False
            or target_downloaded.get("science_measurement") is not False
            or control_downloaded.get("science_measurement") is not False):
        raise ValueError("S23B image receipts are incomplete or outside their approved manifests")

    target_audit, target_rows = verify_audit(
        target_audit_summary_path, target_audit_table_path, expected_count=1043)
    control_audit, control_rows = verify_audit(
        control_audit_summary_path, control_audit_table_path, expected_count=1041)
    target_frames = {row["id"]: row for row in target_rows}
    control_frames = {row["id"]: row for row in control_rows}
    target_mask_summary, target_masks = read_mask_records("target")
    control_mask_summary, control_masks = read_mask_records("control")

    stack_summary = read_json(stack_summary_path)
    profiles = read_tsv(profile_path)
    stack_profile_sha = stack_summary.get("outputs", {}).get("paired_profiles_sha256")
    if (stack_summary.get("status") != "complete_exploratory"
            or stack_summary.get("science_measurement") is not False
            or stack_summary.get("n_pairs_stacked") != EXPECTED_PAIRS
            or stack_profile_sha != sha256(profile_path)
            or len(profiles) != EXPECTED_PAIRS):
        raise ValueError("S23B stack table or summary is incomplete")
    stack_filters = {band: int(count) for band, count in stack_summary.get("n_by_filter", {}).items()}
    if stack_filters != EXPECTED_FILTERS:
        raise ValueError(f"unexpected paired filter counts: {stack_filters}")
    if len({(row["target_id"], row["control_id"]) for row in profiles}) != EXPECTED_PAIRS:
        raise ValueError("duplicate target/control pairs in the saved S23B stack table")

    target_status_rows = {str(row["id"]): row for row in target_status.get("completed", [])}
    control_status_rows = {str(row["id"]): row for row in control_status.get("completed", [])}
    records: list[dict] = []
    for pair in sorted(profiles, key=lambda row: (row["filter"], row["target_id"])):
        target_id = str(pair["target_id"])
        control_id = str(pair["control_id"])
        band = str(pair["filter"])
        target = target_frames.get(target_id)
        control = control_frames.get(control_id)
        target_mask = target_masks.get(target_id)
        control_mask = control_masks.get(control_id)
        if any(value is None for value in (target, control, target_mask, control_mask)):
            raise ValueError(f"missing audited image or mask for {target_id}/{control_id}")
        if (target["primary_filter"] != band or control["filter"] != band
                or target_mask["primary_filter"] != band
                or control_mask["primary_filter"] != band):
            raise ValueError(f"filter mismatch in saved stack pair {target_id}/{control_id}")
        target_saved = target_status_rows.get(target_id)
        control_saved = control_status_rows.get(control_id)
        if (target_saved is None or control_saved is None
                or target_saved["sha256"] != target["sha256"]
                or control_saved["sha256"] != control["sha256"]):
            raise ValueError(f"download receipt differs from image audit for {target_id}/{control_id}")
        for role, ident, frame, mask, base_dir, status_record, mask_summary in (
            ("target", target_id, target, target_mask, target_download / "images", target_saved, target_mask_summary),
            ("matched_control", control_id, control, control_mask, control_download / "images", control_saved, control_mask_summary),
        ):
            image_path = base_dir / frame["file"]
            mask_path = Path(mask["mask_file"]).resolve()
            if (not image_path.is_file() or image_path.stat().st_size != int(frame["bytes"])
                    or not mask_path.is_file()
                    or mask.get("mask_sha256") != mask_summary.get("mask_output_sha256", {}).get(str(mask_path))):
                raise ValueError(f"image/mask output missing or inconsistent for {ident}")
            records.append({
                "pair_id": target_id,
                "cluster_id": ident,
                "role": role,
                "filter": band,
                "image_path": str(image_path.resolve()),
                "image_sha256": frame["sha256"],
                "image_bytes": int(frame["bytes"]),
                "mask_path": str(mask_path),
                "mask_sha256": mask["mask_sha256"],
                "pixel_scale_arcsec": float(frame["pixel_scale_arcsec"]),
                "center_x": float(frame["center_x"]),
                "center_y": float(frame["center_y"]),
            })

    if len(records) != 2 * EXPECTED_PAIRS:
        raise ValueError("exact S23B stack pair list did not produce two image rows per pair")
    for band, expected in EXPECTED_FILTERS.items():
        target_n = sum(row["filter"] == band and row["role"] == "target" for row in records)
        control_n = sum(row["filter"] == band and row["role"] == "matched_control" for row in records)
        if target_n != expected or control_n != expected:
            raise ValueError(f"unexpected {band} pair frame counts")

    provenance_paths = [*required,
        TARGET_AUDIT / "summary.json", TARGET_AUDIT / "verified_target_frames.tsv",
        CONTROL_AUDIT / "summary.json", CONTROL_AUDIT / "verified_control_frames.tsv",
        MASK_ROOT / "target" / "summary.json", MASK_ROOT / "target" / "records.tsv",
        MASK_ROOT / "control" / "summary.json", MASK_ROOT / "control" / "records.tsv",
    ]
    provenance = {str(path.resolve()): sha256(path) for path in sorted(set(provenance_paths))}
    provenance["analysis_helper"] = sha256(Path(noise.__file__).resolve())
    return records, {
        "provenance": provenance,
        "n_pairs": EXPECTED_PAIRS,
        "filter_counts": EXPECTED_FILTERS,
        "target_audit_sha256": sha256(target_audit_summary_path),
        "control_audit_sha256": sha256(control_audit_summary_path),
        "stack_summary_sha256": sha256(stack_summary_path),
        "paired_profile_sha256": sha256(profile_path),
        "target_mask_summary_sha256": sha256(MASK_ROOT / "target" / "summary.json"),
        "control_mask_summary_sha256": sha256(MASK_ROOT / "control" / "summary.json"),
    }


def patch_origins(valid: np.ndarray, *, center_x: float, center_y: float,
                  pixel_scale_arcsec: float, patch_size: int, max_patches: int,
                  min_valid_fraction: float, outer_radius_arcsec: float,
                  seed: int) -> list[tuple[int, int]]:
    height, width = valid.shape
    candidates: list[tuple[int, int]] = []
    radius_pixels = outer_radius_arcsec / pixel_scale_arcsec
    for y0 in range(0, height - patch_size + 1, patch_size):
        cy = y0 + 0.5 * (patch_size - 1)
        for x0 in range(0, width - patch_size + 1, patch_size):
            if float(valid[y0:y0 + patch_size, x0:x0 + patch_size].mean()) < min_valid_fraction:
                continue
            cx = x0 + 0.5 * (patch_size - 1)
            if math.hypot(cx - center_x, cy - center_y) < radius_pixels:
                continue
            candidates.append((y0, x0))
    if len(candidates) <= max_patches:
        return candidates
    rng = np.random.default_rng(seed)
    chosen = rng.choice(len(candidates), size=max_patches, replace=False)
    return [candidates[int(index)] for index in sorted(chosen)]


def analyze_frame(record: dict, index: int, *, patch_size: int, max_patches: int,
                  max_lag: int, min_valid_fraction: float,
                  outer_radius_arcsec: float) -> dict:
    image_path = Path(record["image_path"])
    mask_path = Path(record["mask_path"])
    if sha256(image_path) != record["image_sha256"]:
        raise ValueError(f"image hash changed: {image_path}")
    if sha256(mask_path) != record["mask_sha256"]:
        raise ValueError(f"mask hash changed: {mask_path}")
    with fits.open(image_path, memmap=True, checksum=True) as hdus:
        if len(hdus) < 4 or any(hdus[index].data is None for index in (1, 2, 3)):
            raise ValueError(f"unexpected HSC image HDUs: {image_path}")
        image = np.asarray(hdus[1].data, dtype=np.float32)
        variance = np.asarray(hdus[3].data, dtype=np.float32)
        image_mask = np.asarray(hdus[2].data)
    with fits.open(mask_path, memmap=True, checksum=True) as hdus:
        if "BASE_MASK" not in hdus or "MASK_K2P5" not in hdus:
            raise ValueError(f"saved source mask lacks required planes: {mask_path}")
        base_mask = np.asarray(hdus["BASE_MASK"].data, dtype=bool)
        source_mask = np.asarray(hdus["MASK_K2P5"].data, dtype=bool)
    if (image.ndim != 2 or image.shape != variance.shape or image.shape != image_mask.shape
            or image.shape != base_mask.shape or image.shape != source_mask.shape):
        raise ValueError(f"science, variance, quality, and source-mask grids differ: {record['cluster_id']}")
    if np.any(base_mask & ~source_mask):
        raise ValueError(f"Kron mask does not include all base exclusions: {record['cluster_id']}")
    valid = np.isfinite(image) & np.isfinite(variance) & (variance > 0) & ~source_mask
    seed = zlib.crc32(f"20261107:{record['pair_id']}:{record['role']}:{record['filter']}".encode())
    origins = patch_origins(
        valid, center_x=record["center_x"], center_y=record["center_y"],
        pixel_scale_arcsec=record["pixel_scale_arcsec"], patch_size=patch_size,
        max_patches=max_patches, min_valid_fraction=min_valid_fraction,
        outer_radius_arcsec=outer_radius_arcsec, seed=seed,
    )
    lag_rows = []
    frame_corr: dict[tuple[int, int], float | None] = {}
    for dy in range(-max_lag, max_lag + 1):
        for dx in range(-max_lag, max_lag + 1):
            frame_corr[(dy, dx)] = None
    if len(origins) >= 8:
        patches = []
        for y0, x0 in origins:
            section = np.s_[y0:y0 + patch_size, x0:x0 + patch_size]
            patches.append(noise.lag_correlations(
                image[section], variance[section], valid[section], max_lag=max_lag,
            ))
        for dy in range(-max_lag, max_lag + 1):
            for dx in range(-max_lag, max_lag + 1):
                values = np.asarray([patch[(dy, dx)][0] for patch in patches], dtype=float)
                values = values[np.isfinite(values)]
                median = float(np.median(values)) if len(values) else None
                q16, q84 = (np.quantile(values, [0.16, 0.84]).tolist()
                            if len(values) else (None, None))
                lag_rows.append({
                    "pair_id": record["pair_id"], "cluster_id": record["cluster_id"],
                    "role": record["role"], "filter": record["filter"],
                    "dy": dy, "dx": dx, "n_patches": len(values),
                    "median_variance_normalized_covariance": median,
                    "q16_variance_normalized_covariance": q16,
                    "q84_variance_normalized_covariance": q84,
                })
                frame_corr[(dy, dx)] = median
    factor = sum(value for value in frame_corr.values() if value is not None)
    center_ratio = frame_corr[(0, 0)]
    frame = {
        "record_index": index,
        "pair_id": record["pair_id"], "cluster_id": record["cluster_id"],
        "role": record["role"], "filter": record["filter"],
        "image_path": record["image_path"], "image_sha256": record["image_sha256"],
        "mask_path": record["mask_path"], "mask_sha256": record["mask_sha256"],
        "image_shape": list(image.shape), "pixel_scale_arcsec": record["pixel_scale_arcsec"],
        "center_exclusion_arcsec": outer_radius_arcsec,
        "valid_fraction_whole_image": float(valid.mean()),
        "n_patches": len(origins),
        "analysis_status": "complete" if len(origins) >= 8 else "insufficient_outer_patches",
        "patch_size": patch_size, "max_lag_pixels": max_lag,
        "variance_plane_ratio_at_zero_lag": center_ratio,
        "top_hat_variance_inflation_approx_square_lag_window": float(factor) if len(origins) >= 8 else None,
        "top_hat_sigma_inflation_approx_square_lag_window": math.sqrt(factor) if len(origins) >= 8 and factor > 0 else None,
    }
    return {"frame": frame, "lag_rows": lag_rows}


def write_table(path: Path, rows: list[dict], fields: list[str]) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t",
                                extrasaction="raise", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def summarize_frames(frames: list[dict], n_pairs: int, filters: dict) -> dict:
    by_filter: dict[str, dict] = {}
    for band in ("i", "i2"):
        role_summary = {}
        frame_index = {
            role: {row["pair_id"]: row for row in frames
                   if row["filter"] == band and row["role"] == role
                   and row["analysis_status"] == "complete"}
            for role in ("target", "matched_control")
        }
        for role in ("target", "matched_control"):
            subset = list(frame_index[role].values())
            factors = np.asarray([row["top_hat_variance_inflation_approx_square_lag_window"]
                                  for row in subset], dtype=float)
            centers = np.asarray([row["variance_plane_ratio_at_zero_lag"] for row in subset], dtype=float)
            role_summary[role] = {
                "n_frames_complete": len(subset),
                "n_frames_insufficient_patches": sum(
                    row["filter"] == band and row["role"] == role
                    and row["analysis_status"] != "complete" for row in frames),
                "median_variance_inflation_factor": float(np.median(factors)) if len(factors) else None,
                "q16_q84_variance_inflation_factor": np.quantile(factors, [0.16, 0.84]).tolist() if len(factors) else None,
                "median_sigma_inflation_factor": float(np.median(np.sqrt(factors))) if len(factors) else None,
                "median_variance_plane_ratio_at_zero_lag": float(np.median(centers)) if len(centers) else None,
            }
        paired = [(frame_index["target"][pair_id], frame_index["matched_control"][pair_id])
                  for pair_id in sorted(frame_index["target"].keys() & frame_index["matched_control"].keys())]
        differences = np.asarray([
            target["top_hat_variance_inflation_approx_square_lag_window"]
            - control["top_hat_variance_inflation_approx_square_lag_window"]
            for target, control in paired], dtype=float)
        role_summary["paired_target_minus_control"] = {
            "n_pairs": int(len(differences)),
            "median": float(np.median(differences)) if len(differences) else None,
            "q16_q84_frame_pair_spread": np.quantile(differences, [0.16, 0.84]).tolist() if len(differences) else None,
            "interval_type": "descriptive paired frame-pair spread, not a confidence interval",
        }
        by_filter[band] = role_summary
    return {
        "n_pairs_in_input_stack": n_pairs,
        "filter_pair_counts": filters,
        "n_image_frames_expected": 2 * n_pairs,
        "n_image_frames_complete": sum(row["analysis_status"] == "complete" for row in frames),
        "n_image_frames_insufficient_patches": sum(row["analysis_status"] != "complete" for row in frames),
        "by_filter": by_filter,
    }


def make_figure(frames: list[dict], path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.4), sharey=True)
    for ax, band in zip(axes, ("i", "i2")):
        values = [
            [row["top_hat_variance_inflation_approx_square_lag_window"] for row in frames
             if row["filter"] == band and row["role"] == role and row["analysis_status"] == "complete"]
            for role in ("target", "matched_control")
        ]
        ax.boxplot(values, labels=["Target", "Random control"], showfliers=False,
                   whis=(16, 84), widths=0.58)
        ax.axhline(1.0, color="0.45", linewidth=0.8, linestyle="--")
        ax.set_title(f"HSC-{band}", fontsize=10, loc="left")
        ax.set_ylabel("Approximate 9×9 variance factor")
        ax.tick_params(axis="x", labelrotation=10)
        ax.grid(axis="y", alpha=0.25)
    fig.suptitle("S23B outer-field coadd residual covariance", fontsize=10)
    fig.tight_layout()
    fig.savefig(path, dpi=220)
    plt.close(fig)


def run(output: Path, *, resume: bool, patch_size: int, max_patches: int,
        max_lag: int, min_valid_fraction: float, outer_radius_arcsec: float) -> dict:
    output = output.resolve()
    if patch_size < 16 or max_lag < 0 or 2 * max_lag + 8 > patch_size:
        raise ValueError("invalid patch_size/max_lag combination")
    if max_patches < 8 or not 0 < min_valid_fraction <= 1 or outer_radius_arcsec < 0:
        raise ValueError("invalid patch-selection settings")
    records, setup = build_records()
    config = {
        "schema": SCHEMA, "target_manifest_sha256": TARGET_MANIFEST_SHA256,
        "control_manifest_sha256": CONTROL_MANIFEST_SHA256,
        "stacked_pair_count": EXPECTED_PAIRS, "filter_pair_counts": EXPECTED_FILTERS,
        "patch_size": patch_size, "max_patches_per_frame": max_patches,
        "max_lag_pixels": max_lag, "minimum_valid_pixel_fraction": min_valid_fraction,
        "excluded_center_radius_arcsec": outer_radius_arcsec,
        "variance_normalization": "delivered VARIANCE plane",
        "mask_plane": "saved MASK_K2P5, verified to include BASE_MASK",
        "analysis_script_sha256": sha256(Path(__file__).resolve()),
        "input_provenance": setup["provenance"],
    }
    config_hash = hashlib.sha256(json.dumps(config, sort_keys=True,
                                          separators=(",", ":")).encode()).hexdigest()
    config_path = output / "config.json"
    state_path = output / "state.json"
    frames_dir = output / "frame_results"
    summary_path = output / "summary.json"
    if output.exists():
        if not resume:
            raise FileExistsError(f"refusing to reuse {output}; pass --resume to continue its exact run")
        if not config_path.is_file() or read_json(config_path) != config:
            raise ValueError("resume configuration differs from the pinned run")
        if summary_path.exists():
            raise FileExistsError(f"run is already complete: {summary_path}")
    else:
        if resume:
            raise FileNotFoundError(f"cannot resume a run that does not exist: {output}")
        output.mkdir(parents=True, exist_ok=False)
        frames_dir.mkdir(parents=True, exist_ok=False)
        write_json_atomic(config_path, config)
    frames_dir.mkdir(parents=True, exist_ok=True)
    completed = 0
    frames: list[dict] = []
    lags: list[dict] = []
    write_json_atomic(state_path, {
        "status": "running", "config_sha256": config_hash,
        "n_frames_expected": len(records), "n_frames_complete": 0,
        "science_measurement": False,
    })
    for index, record in enumerate(records):
        safe_id = re.sub(r"[^A-Za-z0-9_.-]+", "_", record["pair_id"])
        result_path = frames_dir / f"{index:04d}_{safe_id}_{record['role']}.json"
        if result_path.is_file():
            result = read_json(result_path)
            frame = result.get("frame", {})
            if (result.get("config_sha256") != config_hash
                    or frame.get("record_index") != index
                    or frame.get("image_sha256") != record["image_sha256"]
                    or frame.get("mask_sha256") != record["mask_sha256"]
                    or frame.get("pair_id") != record["pair_id"]
                    or frame.get("role") != record["role"]):
                raise ValueError(f"resume frame result does not match pinned input: {result_path}")
            if sha256(Path(record["image_path"])) != record["image_sha256"] or sha256(Path(record["mask_path"])) != record["mask_sha256"]:
                raise ValueError(f"resume input changed since frame result was saved: {record['pair_id']}")
        else:
            print(f"[{index + 1}/{len(records)}] {record['pair_id']} {record['role']} {record['filter']}", flush=True)
            result = analyze_frame(
                record, index, patch_size=patch_size, max_patches=max_patches,
                max_lag=max_lag, min_valid_fraction=min_valid_fraction,
                outer_radius_arcsec=outer_radius_arcsec,
            )
            result["config_sha256"] = config_hash
            write_json_atomic(result_path, result)
        frames.append(result["frame"])
        lags.extend(result["lag_rows"])
        completed += 1
        if completed % 25 == 0 or completed == len(records):
            write_json_atomic(state_path, {
                "status": "running", "config_sha256": config_hash,
                "n_frames_expected": len(records), "n_frames_complete": completed,
                "last_pair_id": record["pair_id"], "last_role": record["role"],
                "science_measurement": False,
            })

    frame_path = output / "frame_noise_correlation.tsv"
    lag_path = output / "lag_correlation.tsv"
    frame_fields = list(frames[0])
    lag_fields = list(lags[0]) if lags else [
        "pair_id", "cluster_id", "role", "filter", "dy", "dx", "n_patches",
        "median_variance_normalized_covariance", "q16_variance_normalized_covariance",
        "q84_variance_normalized_covariance",
    ]
    write_table(frame_path, frames, frame_fields)
    write_table(lag_path, lags, lag_fields)
    figure_path = output / "variance_inflation_by_filter.png"
    make_figure(frames, figure_path)
    summary = {
        "schema": SCHEMA,
        "status": "complete_exploratory_residual_covariance_diagnostic",
        "science_measurement": False,
        "images_ready_for_science": False,
        "sample": summarize_frames(frames, EXPECTED_PAIRS, EXPECTED_FILTERS),
        "method": {
            "pair_selection": "exact 1,022 target/control pairs in the saved exploratory S23B flux stack",
            "mask": "saved K=2.5 mask for targets and matched controls; BASE_MASK inclusion checked per frame",
            "patch_selection": f"non-overlapping {patch_size}x{patch_size} native-pixel patches with at least {min_valid_fraction:.2f} valid fraction; patch centers beyond {outer_radius_arcsec:g} arcsec from the target/random center",
            "patches_per_frame": max_patches,
            "maximum_lag_pixels": max_lag,
            "residual_model": "local affine plane with iterative 4-sigma clipping",
            "normalization": "residual divided by square root of the delivered HSC VARIANCE plane",
            "variance_factor": "sum of the finite variance-normalized lag covariance grid; approximate square top-hat factor only",
            "no_network_or_new_data": True,
        },
        "inputs_sha256": setup["provenance"],
        "config_sha256": config_hash,
        "outputs": {
            frame_path.name: sha256(frame_path), lag_path.name: sha256(lag_path),
            figure_path.name: sha256(figure_path),
        },
        "limitations": [
            "The finite local lag kernel is not propagated through the Lanczos metric warp or exact signed annular estimator weights.",
            "The outer-field patch diagnostic can still contain unresolved sources and diffuse cluster emission.",
            "Frame-pair spreads are descriptive, not confidence intervals or a substitute for connected sky-block covariance.",
            "This analysis does not include source Poisson noise, residual sky-transfer uncertainty, PSF-wing uncertainty, or member-selection uncertainty.",
            "No correction is applied to the ICL stack; science_measurement and images_ready_for_science remain false.",
        ],
    }
    write_json_atomic(summary_path, summary)
    write_json_atomic(state_path, {
        "status": "complete", "config_sha256": config_hash,
        "n_frames_expected": len(records), "n_frames_complete": len(frames),
        "summary_sha256": sha256(summary_path), "science_measurement": False,
    })
    print(json.dumps({
        "status": summary["status"], "n_pairs": EXPECTED_PAIRS,
        "n_frames_expected": len(records),
        "n_frames_complete": summary["sample"]["n_image_frames_complete"],
        "n_frames_insufficient_patches": summary["sample"]["n_image_frames_insufficient_patches"],
        "by_filter": summary["sample"]["by_filter"],
        "science_measurement": False, "output": str(output),
    }, indent=2), flush=True)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--resume", action="store_true")
    parser.add_argument("--patch-size", type=int, default=64)
    parser.add_argument("--max-patches", type=int, default=16)
    parser.add_argument("--max-lag", type=int, default=4)
    parser.add_argument("--min-valid-fraction", type=float, default=0.80)
    parser.add_argument("--outer-radius-arcsec", type=float, default=60.0)
    args = parser.parse_args()
    run(args.output, resume=args.resume, patch_size=args.patch_size,
        max_patches=args.max_patches, max_lag=args.max_lag,
        min_valid_fraction=args.min_valid_fraction,
        outer_radius_arcsec=args.outer_radius_arcsec)


if __name__ == "__main__":
    main()
