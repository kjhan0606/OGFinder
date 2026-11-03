#!/usr/bin/env python3
"""Measure matched target/control coadd residual covariance for the full ICL stack.

This CPU audit reuses the acquired 800 target images, 798 matched controls,
and saved K=2.5 target masks. It retains the 783 same-filter pairs used by the
full target-control stack and measures each image's variance-normalized local
residual lag kernel. It does not modify images, submit HSC queries, or provide
a final ICL uncertainty model.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import zlib

import numpy as np
from astropy.io import fits

import estimate_hsc_coadd_noise_correlation as noise
from evaluate_prestack_pair_support import pair_block_labels


DEFAULT_SAMPLE = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
DEFAULT_OUTPUT = DEFAULT_SAMPLE / "icl_coadd_noise_correlation_full_2026-11-03_01"


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


def write_json(path: Path, value: dict) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, ensure_ascii=False, indent=2,
                                    sort_keys=True, allow_nan=False) + "\n",
                         encoding="utf-8")
    temporary.replace(path)


def frame_correlations(record: dict, *, patch_size: int, max_patches: int,
                       max_lag: int, min_valid_fraction: float) -> tuple[dict, list[dict]]:
    """Use the archived K=2.5 mask on targets and the established SEP mask on controls."""
    image_path = Path(record["path"])
    if sha256(image_path) != record["source_image_sha256"]:
        raise ValueError(f"{record['cluster_id']}/{record['role']}: image hash differs from acquisition record")
    with fits.open(image_path, memmap=True, checksum=True) as hdus:
        image = np.asarray(hdus[1].data, dtype=np.float32)
        image_mask = np.asarray(hdus[2].data, dtype=np.int64)
        variance = np.asarray(hdus[3].data, dtype=np.float32)
    if record["role"] == "target":
        mask_path = Path(record["target_mask_path"])
        if sha256(mask_path) != record["target_mask_sha256"]:
            raise ValueError(f"{record['cluster_id']}: target mask checksum mismatch")
        with fits.open(mask_path, memmap=True, checksum=True) as hdus:
            excluded = np.asarray(hdus["MASK_K2P5"].data, dtype=bool)
    else:
        excluded = noise._matched_control_source_mask(
            image, variance, image_mask, float(record["pixel_scale_arcsec"]),
        )
    if image.shape != variance.shape or image.shape != excluded.shape:
        raise ValueError(f"{record['cluster_id']}/{record['role']}: image, variance, and mask grids differ")
    valid = np.isfinite(image) & np.isfinite(variance) & (variance > 0) & ~excluded
    seed = zlib.crc32(f"20261103:{record['cluster_id']}:{record['filter']}:{record['role']}".encode())
    origins = noise._sample_patch_origins(
        valid, patch_size=patch_size, max_patches=max_patches, seed=seed,
        min_valid_fraction=min_valid_fraction,
    )
    if len(origins) < 8:
        raise ValueError(f"{record['cluster_id']}/{record['role']}: only {len(origins)} usable patches")
    patches = []
    for y0, x0 in origins:
        section = np.s_[y0:y0 + patch_size, x0:x0 + patch_size]
        patches.append(noise.lag_correlations(image[section], variance[section],
                                              valid[section], max_lag=max_lag))
    lag_rows, frame_corr = [], {}
    for dy in range(-max_lag, max_lag + 1):
        for dx in range(-max_lag, max_lag + 1):
            values = np.asarray([patch[(dy, dx)][0] for patch in patches], dtype=float)
            values = values[np.isfinite(values)]
            median = float(np.median(values)) if len(values) else None
            q16, q84 = np.quantile(values, [0.16, 0.84]).tolist() if len(values) else (None, None)
            lag_rows.append({"cluster_id": record["cluster_id"], "role": record["role"],
                             "filter": record["filter"], "dy": dy, "dx": dx,
                             "n_patches": len(values),
                             "median_variance_normalized_covariance": median,
                             "q16_variance_normalized_covariance": q16,
                             "q84_variance_normalized_covariance": q84})
            frame_corr[(dy, dx)] = median
    factor = sum(value for value in frame_corr.values() if value is not None)
    center = frame_corr[(0, 0)]
    frame = {
        "cluster_id": record["cluster_id"], "role": record["role"], "filter": record["filter"],
        "image_id": record["image_id"], "image_sha256": record["source_image_sha256"],
        "mask_path": record["target_mask_path"] if record["role"] == "target" else "SEP/Kron 2.5 mask derived with prepare_detection_psf.py",
        "mask_sha256": record["target_mask_sha256"] if record["role"] == "target" else None,
        "image_shape": list(image.shape), "usable_pixels_fraction": float(valid.mean()),
        "patch_size": patch_size, "n_patches": len(origins), "max_lag_pixels": max_lag,
        "variance_plane_ratio_at_zero_lag": center,
        "top_hat_variance_inflation_approx_square_lag_window": float(factor),
        "top_hat_sigma_inflation_approx_square_lag_window": float(np.sqrt(factor)) if factor > 0 else None,
        "connected_5deg_block": record["connected_5deg_block"],
    }
    return frame, lag_rows


def build_records(sample: Path) -> tuple[list[dict], dict]:
    target_manifest_path = sample / "download_01/download_manifest.json"
    target_status_path = sample / "download_01/download_status.json"
    control_manifest_path = sample / "control_download_01/download_manifest.json"
    control_status_path = sample / "control_download_01/download_status.json"
    mask_summary_path = sample / "prestack_cluster_masks_2026-10-16_01/summary.json"
    mask_records_path = sample / "prestack_cluster_masks_2026-10-16_01/records.tsv"
    required = (target_manifest_path, target_status_path, control_manifest_path,
                control_status_path, mask_summary_path, mask_records_path)
    for path in required:
        if not path.is_file():
            raise FileNotFoundError(path)

    target_manifest, target_status = read_json(target_manifest_path), read_json(target_status_path)
    control_manifest, control_status = read_json(control_manifest_path), read_json(control_status_path)
    mask_summary = read_json(mask_summary_path)
    target_masks = read_tsv(mask_records_path)
    if (target_status.get("status") != "complete" or target_status.get("n_expected") != 800
            or target_status.get("downloaded") != 800
            or control_status.get("status") != "complete" or control_status.get("n_expected") != 798
            or control_status.get("downloaded") != 798
            or control_manifest.get("frame_role") != "matched_control"
            or mask_summary.get("status") != "complete" or mask_summary.get("n_success") != 800
            or len(target_manifest.get("clusters", [])) != 800
            or len(control_manifest.get("clusters", [])) != 798
            or len(target_masks) != 800):
        raise ValueError("full target/control image or K=2.5 mask provenance is incomplete")
    mask_records_expected = mask_summary.get("output_sha256", {}).get(str(mask_records_path))
    if mask_records_expected and mask_records_expected != sha256(mask_records_path):
        raise ValueError("K=2.5 mask record table differs from its saved checksum")

    target_downloaded = {str(row["id"]): row for row in target_status["completed"]}
    control_downloaded = {str(row["id"]): row for row in control_status["completed"]}
    target_by_id = {str(row["id"]): row for row in target_manifest["clusters"]}
    control_by_target = {str(row["target_cluster_id"]): row for row in control_manifest["clusters"]}
    mask_by_id = {str(row["id"]): row for row in target_masks}
    if (len(target_downloaded) != 800 or len(control_downloaded) != 798
            or len(target_by_id) != 800 or len(control_by_target) != 798
            or len(mask_by_id) != 800):
        raise ValueError("duplicate or missing target/control IDs")

    records = []
    pair_rows = []
    for ident, target in target_by_id.items():
        control = control_by_target.get(ident)
        target_saved = target_downloaded.get(ident)
        mask = mask_by_id.get(ident)
        if target_saved is None or mask is None or control is None:
            continue
        target_review = target_saved.get("review", {})
        control_saved = control_downloaded.get(str(control["id"]))
        if control_saved is None:
            continue
        control_review = control_saved.get("review", {})
        target_filter = str(target_review.get("coadd_header_filter") or target_review.get("filter"))
        control_filter = str(control_review.get("coadd_header_filter") or control_review.get("filter"))
        if target_filter != control_filter:
            continue
        if (str(mask.get("primary_filter")) != target_filter
                or str(control.get("coadd_header_filter")) != control_filter
                or str(target_saved.get("file")) != Path(str(mask["image"])).name):
            raise ValueError(f"filter/image/mask provenance mismatch for {ident}")
        target_path = sample / "download_01/images" / str(target_saved["file"])
        control_path = sample / "control_download_01/images" / str(control_saved["file"])
        mask_path = Path(mask["mask_file"])
        if not target_path.is_file() or not control_path.is_file() or not mask_path.is_file():
            raise FileNotFoundError(f"missing input for {ident}")
        if target_saved.get("bytes") != target_path.stat().st_size or control_saved.get("bytes") != control_path.stat().st_size:
            raise ValueError(f"input size differs from acquisition record for {ident}")
        if mask_summary.get("mask_output_sha256", {}).get(str(mask_path)) != mask.get("mask_sha256"):
            raise ValueError(f"target mask checksum index mismatch for {ident}")

        center = target
        pair_rows.append({
            "cluster_id": ident, "filter": target_filter,
            "target_ra_deg": float(center["ra_deg"]), "target_dec_deg": float(center["dec_deg"]),
            "control_ra_deg": float(control["ra_deg"]), "control_dec_deg": float(control["dec_deg"]),
        })
        for role, image_path, saved, mask_file in (
                ("target", target_path, target_saved, mask_path),
                ("matched_control", control_path, control_saved, None)):
            records.append({
                "cluster_id": ident, "filter": target_filter, "role": role,
                "image_id": ident if role == "target" else str(control["id"]),
                "path": str(image_path.resolve()),
                "source_image_sha256": str(saved["sha256"]),
                "pixel_scale_arcsec": float(saved["review"]["pixel_scale_arcsec"]),
                "target_mask_path": str(mask_file) if mask_file else None,
                "target_mask_sha256": str(mask["mask_sha256"]),
            })

    if len(pair_rows) != 783:
        raise ValueError(f"expected 783 same-filter pairs, found {len(pair_rows)}")
    filter_counts = {band: sum(row["filter"] == band for row in pair_rows) for band in ("i", "i2")}
    if filter_counts != {"i": 302, "i2": 481}:
        raise ValueError(f"unexpected same-filter pair counts: {filter_counts}")

    target_coords = [{"ra_deg": row["target_ra_deg"], "dec_deg": row["target_dec_deg"]}
                     for row in pair_rows]
    control_coords = [{"ra_deg": row["control_ra_deg"], "dec_deg": row["control_dec_deg"]}
                      for row in pair_rows]
    _, labels, n_target_cells, n_control_cells, n_blocks, largest_block = pair_block_labels(
        target_coords, control_coords, np.ones(len(pair_rows), dtype=bool))
    for row, label in zip(pair_rows, labels):
        row["connected_5deg_block"] = str(label)
    for record in records:
        record["connected_5deg_block"] = next(
            row["connected_5deg_block"] for row in pair_rows if row["cluster_id"] == record["cluster_id"])
    provenance_paths = [*required]
    provenance = {str(path): sha256(path) for path in provenance_paths}
    provenance.update({
        "noise_helper_script": sha256(Path(noise.__file__).resolve()),
        "target_control_pairing": "783 same-filter pairs selected by saved target/control FILTER headers",
        "target_mask": "saved proposed_masks.fits.gz MASK_K2P5 plane; controls use SEP/Kron-2.5 masks",
        "same_filter_pair_counts": filter_counts,
        "connected_5deg_blocks": n_blocks,
        "target_sky_cells": n_target_cells,
        "control_sky_cells": n_control_cells,
        "largest_block_pairs": largest_block,
    })
    return records, {"provenance": provenance, "pair_rows": pair_rows}


def run(sample: Path, output: Path, patch_size: int, max_patches: int,
        max_lag: int, min_valid_fraction: float) -> dict:
    sample, output = sample.resolve(), output.resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite {output}")
    if not os.environ.get("SLURM_JOB_ID"):
        raise RuntimeError("full image audit must run inside a Slurm allocation on syntax")
    records, setup = build_records(sample)
    output.mkdir(parents=True, exist_ok=False)
    progress_path = output / "progress.json"
    details, lag_rows = [], []
    write_json(progress_path, {"status": "running", "n_frames_done": 0,
                               "n_frames_expected": len(records),
                               "slurm_job_id": os.environ["SLURM_JOB_ID"]})
    try:
        for index, record in enumerate(records, 1):
            print(f"[{index}/{len(records)}] covariance {record['cluster_id']} "
                  f"{record['role']} {record['filter']}", flush=True)
            frame, rows = frame_correlations(
                record, patch_size=patch_size, max_patches=max_patches,
                max_lag=max_lag, min_valid_fraction=min_valid_fraction,
            )
            frame["connected_5deg_block"] = record["connected_5deg_block"]
            details.append(frame)
            lag_rows.extend(rows)
            write_json(progress_path, {"status": "running", "n_frames_done": index,
                                       "n_frames_expected": len(records),
                                       "last_cluster_id": record["cluster_id"],
                                       "last_role": record["role"],
                                       "slurm_job_id": os.environ["SLURM_JOB_ID"]})
        detail_path = output / "frame_noise_correlation.tsv"
        lag_path = output / "lag_correlation.tsv"
        for path, rows in ((detail_path, details), (lag_path, lag_rows)):
            with path.open("x", encoding="utf-8", newline="") as stream:
                writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t",
                                        extrasaction="raise", lineterminator="\n")
                writer.writeheader()
                writer.writerows(rows)
        summary = {
            "schema": "astrafek-icl-full-sample-coadd-residual-covariance-v1",
            "status": "complete_exploratory", "science_measurement": False,
            "sample": {"same_filter_pairs": 783, "image_frames": len(details),
                       "filter_counts": {band: sum(row["filter"] == band and row["role"] == "target"
                                                    for row in details) for band in ("i", "i2")}},
            "method": {
                "patch_size": patch_size, "max_patches_per_image": max_patches,
                "max_lag_pixels": max_lag, "min_valid_fraction": min_valid_fraction,
                "target_mask": "saved K=2.5 mask; includes source masks and base invalid pixels",
                "control_mask": "SEP/Kron-2.5 and bad/bright-star bit masks derived with the existing Astrafek helper",
                "residual": "local affine plane fit with iterative 4-sigma clipping; normalize by delivered variance plane",
                "factor": "sum of the finite 9x9 variance-normalized lag covariance grid; approximate large-top-hat variance factor only",
                "image_content_rehashed": True,
                "no_new_hsc_query_or_download": True,
            },
            "by_filter_and_role": {}, "paired_target_control": {},
            "connected_block_bootstrap": False,
            "inputs": setup["provenance"],
            "outputs": {"frame_table": str(detail_path), "lag_table": str(lag_path)},
            "limitations": [
                "A finite local lag kernel is not the exact covariance of the radial ICL stack.",
                "This sample-wide diagnostic does not include source Poisson noise, sky-transfer uncertainty, PSF wings, or mask-growth response.",
                "Per-frame distributions are descriptive and do not replace cluster-level spatial-block uncertainty.",
                "Spatial blocks are recorded for later paired covariance analysis but are not resampled in this job.",
            ],
        }
        for band in ("i", "i2"):
            band_summary = {}
            target = {row["cluster_id"]: row for row in details
                      if row["filter"] == band and row["role"] == "target"}
            control = {row["cluster_id"]: row for row in details
                       if row["filter"] == band and row["role"] == "matched_control"}
            if set(target) != set(control):
                raise ValueError(f"target/control covariance coverage differs in {band}")
            for role, values in (("target", list(target.values())),
                                 ("matched_control", list(control.values()))):
                factors = np.asarray([row["top_hat_variance_inflation_approx_square_lag_window"]
                                      for row in values], dtype=float)
                centers = np.asarray([row["variance_plane_ratio_at_zero_lag"] for row in values], dtype=float)
                band_summary[role] = {
                    "n_images": len(values), "median_variance_factor": float(np.median(factors)),
                    "q16_q84_variance_factor": np.quantile(factors, [0.16, 0.84]).tolist(),
                    "median_sigma_inflation": float(np.median(np.sqrt(factors))),
                    "median_zero_lag_variance_ratio": float(np.median(centers)),
                }
            delta = np.asarray([
                target[ident]["top_hat_variance_inflation_approx_square_lag_window"]
                - control[ident]["top_hat_variance_inflation_approx_square_lag_window"]
                for ident in sorted(target)], dtype=float)
            band_summary["paired_target_minus_control"] = {
                "median": float(np.median(delta)),
                "q16_q84": np.quantile(delta, [0.16, 0.84]).tolist(),
                "n_pairs": len(delta),
            }
            summary["by_filter_and_role"][band] = band_summary
        summary["output_sha256"] = {str(path): sha256(path) for path in (detail_path, lag_path)}
        summary_path = output / "summary.json"
        write_json(summary_path, summary)
        summary["summary_sha256"] = sha256(summary_path)
        write_json(progress_path, {"status": "complete", "n_frames_done": len(details),
                                   "n_frames_expected": len(records),
                                   "summary": str(summary_path),
                                   "summary_sha256": summary["summary_sha256"],
                                   "slurm_job_id": os.environ["SLURM_JOB_ID"]})
        print(json.dumps({"status": summary["status"], "science_measurement": False,
                          "sample": summary["sample"], "by_filter_and_role": summary["by_filter_and_role"],
                          "summary_sha256": summary["summary_sha256"]}, ensure_ascii=False, indent=2))
        return summary
    except Exception as exc:
        write_json(progress_path, {"status": "stopped_on_error", "n_frames_done": len(details),
                                   "n_frames_expected": len(records), "error": repr(exc),
                                   "slurm_job_id": os.environ["SLURM_JOB_ID"]})
        raise


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", type=Path, default=DEFAULT_SAMPLE)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--patch-size", type=int, default=64)
    parser.add_argument("--max-patches", type=int, default=96)
    parser.add_argument("--max-lag", type=int, default=4)
    parser.add_argument("--min-valid-fraction", type=float, default=0.80)
    args = parser.parse_args()
    run(args.sample, args.output, args.patch_size, args.max_patches,
        args.max_lag, args.min_valid_fraction)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
