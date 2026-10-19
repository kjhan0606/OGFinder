#!/usr/bin/env python3
"""Estimate short-range coadd pixel correlation from existing HSC images.

This is a technical diagnostic, not a source-Poisson model or science error
bar. It uses masked pixels in the already downloaded ICL target coadds, removes
a local affine surface in each sampled patch, normalizes by the delivered
variance plane, and measures robust lag correlations. No network or HSC query
is used.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import zlib

import numpy as np
from astropy.io import fits


ROOT = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
DEFAULT_METADATA = ROOT / "icl_coadd_noise_metadata_2026-10-19_01/summary.json"
DEFAULT_PREP = ROOT / "satellite_detection_psf_2026-10-18_01"
DEFAULT_OUTPUT = ROOT / "icl_coadd_noise_correlation_2026-10-19_01"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _fit_plane_residual(image: np.ndarray, variance: np.ndarray,
                        valid: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Subtract a clipped affine plane and return normalized residuals/mask."""
    height, width = image.shape
    yy, xx = np.indices((height, width), dtype=np.float64)
    x = (xx - (width - 1) / 2) / max(width, 1)
    y = (yy - (height - 1) / 2) / max(height, 1)
    design = np.column_stack((np.ones(valid.sum()), x[valid], y[valid]))
    values = image[valid].astype(np.float64, copy=False)
    keep = np.ones(len(values), dtype=bool)
    coefficients = np.zeros(3, dtype=np.float64)
    for _ in range(4):
        if keep.sum() < 64:
            raise ValueError("too few pixels remain after local clipping")
        coefficients, *_ = np.linalg.lstsq(design[keep], values[keep], rcond=None)
        residual = values - design @ coefficients
        median = np.median(residual[keep])
        scale = 1.4826 * np.median(np.abs(residual[keep] - median))
        if not np.isfinite(scale) or scale <= 0:
            break
        updated = np.abs(residual - median) <= 4.0 * scale
        if np.array_equal(updated, keep):
            break
        keep = updated

    plane = coefficients[0] + coefficients[1] * x + coefficients[2] * y
    normalized = (image.astype(np.float64) - plane) / np.sqrt(variance)
    accepted = valid & np.isfinite(normalized) & (np.abs(normalized) <= 5.0)
    return normalized, accepted


def lag_correlations(image: np.ndarray, variance: np.ndarray, valid: np.ndarray,
                     *, max_lag: int = 4) -> dict[tuple[int, int], tuple[float, int]]:
    """Measure normalized residual covariance at all lags in a patch."""
    if image.ndim != 2 or image.shape != variance.shape or image.shape != valid.shape:
        raise ValueError("image, variance, and valid mask must be matching 2D arrays")
    if min(image.shape) < 2 * max_lag + 8:
        raise ValueError("patch is too small for the requested maximum lag")
    if valid.sum() < 64:
        raise ValueError("patch has too few valid pixels")
    normalized, accepted = _fit_plane_residual(image, variance, valid)
    result = {}
    height, width = image.shape
    for dy in range(-max_lag, max_lag + 1):
        for dx in range(-max_lag, max_lag + 1):
            y0, y1 = max(0, -dy), min(height, height - dy)
            x0, x1 = max(0, -dx), min(width, width - dx)
            a = normalized[y0:y1, x0:x1]
            b = normalized[y0 + dy:y1 + dy, x0 + dx:x1 + dx]
            pair_valid = (accepted[y0:y1, x0:x1]
                          & accepted[y0 + dy:y1 + dy, x0 + dx:x1 + dx])
            count = int(pair_valid.sum())
            if count < 32:
                result[(dy, dx)] = (float("nan"), count)
            else:
                result[(dy, dx)] = (float(np.mean(a[pair_valid] * b[pair_valid])), count)
    return result


def _sample_patch_origins(valid: np.ndarray, *, patch_size: int, max_patches: int,
                          seed: int, min_valid_fraction: float) -> list[tuple[int, int]]:
    height, width = valid.shape
    candidates = []
    for y0 in range(0, height - patch_size + 1, patch_size):
        for x0 in range(0, width - patch_size + 1, patch_size):
            fraction = float(valid[y0:y0 + patch_size, x0:x0 + patch_size].mean())
            if fraction >= min_valid_fraction:
                candidates.append((y0, x0))
    if len(candidates) <= max_patches:
        return candidates
    rng = np.random.default_rng(seed)
    chosen = rng.choice(len(candidates), size=max_patches, replace=False)
    return [candidates[int(i)] for i in sorted(chosen)]


def frame_correlations(record: dict, prep_dir: Path, *, patch_size: int,
                       max_patches: int, max_lag: int,
                       min_valid_fraction: float) -> tuple[dict, list[dict]]:
    ident, band = str(record["cluster_id"]), str(record["filter"])
    image_path = Path(record["path"])
    if sha256(image_path) != record["source_image_sha256_from_audited_run"]:
        raise ValueError(f"{ident}/{band}: image hash differs from audited input")
    prep_path = prep_dir / f"cluster_{ident}" / "masks.fits.gz"
    if not prep_path.is_file():
        raise FileNotFoundError(prep_path)
    with fits.open(image_path, memmap=True, checksum=True) as hdus:
        image = np.asarray(hdus[1].data, dtype=np.float32)
        variance = np.asarray(hdus[3].data, dtype=np.float32)
    with fits.open(prep_path, memmap=True, checksum=True) as hdus:
        base = np.asarray(hdus["BASE_MASK"].data, dtype=bool)
    if image.shape != variance.shape or image.shape != base.shape:
        raise ValueError(f"{ident}/{band}: science, variance, and mask shapes differ")
    valid = np.isfinite(image) & np.isfinite(variance) & (variance > 0) & ~base
    seed = zlib.crc32(f"20261019:{ident}:{band}".encode("utf-8"))
    origins = _sample_patch_origins(
        valid, patch_size=patch_size, max_patches=max_patches, seed=seed,
        min_valid_fraction=min_valid_fraction,
    )
    if len(origins) < 8:
        raise ValueError(f"{ident}/{band}: only {len(origins)} usable patches")
    patches = []
    for y0, x0 in origins:
        sl = np.s_[y0:y0 + patch_size, x0:x0 + patch_size]
        patches.append(lag_correlations(image[sl], variance[sl], valid[sl], max_lag=max_lag))
    rows = []
    frame_corr = {}
    for dy in range(-max_lag, max_lag + 1):
        for dx in range(-max_lag, max_lag + 1):
            values = np.asarray([p[(dy, dx)][0] for p in patches], dtype=float)
            values = values[np.isfinite(values)]
            median = float(np.median(values)) if len(values) else None
            q16, q84 = (np.quantile(values, [0.16, 0.84]).tolist()
                        if len(values) else (None, None))
            rows.append({"cluster_id": ident, "filter": band, "dy": dy, "dx": dx,
                         "n_patches": len(values), "median_rho": median,
                         "q16_rho": q16, "q84_rho": q84})
            frame_corr[(dy, dx)] = median
    factor = sum(value for value in frame_corr.values() if value is not None)
    center = frame_corr[(0, 0)]
    result = {
        "cluster_id": ident, "filter": band, "image_sha256": record["source_image_sha256_from_audited_run"],
        "mask_path": str(prep_path), "mask_sha256": sha256(prep_path),
        "image_shape": list(image.shape), "usable_pixels_fraction": float(valid.mean()),
        "patch_size": patch_size, "n_patches": len(origins),
        "max_lag_pixels": max_lag, "central_rho": center,
        "top_hat_variance_inflation_approx_square_lag_window": float(factor),
        "top_hat_sigma_inflation_approx_square_lag_window": math.sqrt(factor) if factor > 0 else None,
    }
    return result, rows


def run(metadata_path: Path, prep_dir: Path, output_dir: Path, *, patch_size: int = 64,
        max_patches: int = 96, max_lag: int = 4,
        min_valid_fraction: float = 0.80) -> dict:
    if output_dir.exists():
        raise FileExistsError(f"refusing to overwrite {output_dir}")
    if patch_size < 16 or max_lag < 0 or 2 * max_lag + 8 > patch_size:
        raise ValueError("require patch_size >= 16 and 0 <= 2*max_lag+8 <= patch_size")
    if max_patches < 8 or not 0 < min_valid_fraction <= 1:
        raise ValueError("max_patches must be >= 8 and min_valid_fraction in (0, 1]")
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    if (metadata.get("schema") != "astrafek-icl-coadd-noise-metadata-audit-v1"
            or metadata.get("status") != "verified_metadata_only"
            or metadata.get("science_measurement") is not False
            or metadata.get("n_frames") != 57
            or metadata.get("variance_plane_present_in_all_frames") is not True):
        raise ValueError("expected the audited 57-frame coadd-metadata input")
    details, lag_rows = [], []
    for index, record in enumerate(metadata["frame_headers"], 1):
        print(f"[{index}/57] coadd noise correlation {record['cluster_id']} {record['filter']}", flush=True)
        frame, rows = frame_correlations(
            record, prep_dir, patch_size=patch_size, max_patches=max_patches,
            max_lag=max_lag, min_valid_fraction=min_valid_fraction,
        )
        details.append(frame)
        lag_rows.extend(rows)
    by_filter = {}
    for band in sorted({row["filter"] for row in details}):
        subset = [row for row in details if row["filter"] == band]
        factors = np.asarray([row["top_hat_variance_inflation_approx_square_lag_window"]
                              for row in subset], dtype=float)
        central = np.asarray([row["central_rho"] for row in subset], dtype=float)
        by_filter[band] = {
            "n_frames": len(subset),
            "median_square_window_variance_inflation_approx": float(np.median(factors)),
            "q16_q84_square_window_variance_inflation_approx": np.quantile(factors, [0.16, 0.84]).tolist(),
            "median_normalized_neighbor_covariance": float(np.median(central)),
        }
    output_dir.mkdir(parents=True, exist_ok=False)
    frame_path, lag_path = output_dir / "frame_noise_correlation.tsv", output_dir / "lag_correlation.tsv"
    for path, rows in ((frame_path, details), (lag_path, lag_rows)):
        with path.open("x", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)
    summary = {
        "schema": "astrafek-hsc-coadd-noise-correlation-v1",
        "status": "complete_exploratory", "science_measurement": False,
        "n_frames": len(details), "patch_size_pixels": patch_size,
        "max_patches_per_frame": max_patches, "max_lag_pixels": max_lag,
        "minimum_valid_fraction_per_patch": min_valid_fraction,
        "method": "random non-overlapping 64-pixel patches; reject source/bad-mask pixels; fit a local affine plane with 4-sigma iterative clipping; normalize residuals by sqrt(delivered variance); estimate robust median lag covariance across patches",
        "by_filter": by_filter,
        "input_sha256": {str(metadata_path): sha256(metadata_path)},
        "frame_results_sha256": sha256(frame_path), "lag_results_sha256": sha256(lag_path),
        "limitations": [
            "patch-scale residuals can include unresolved sources and real diffuse structure despite masking and local plane removal",
            "the square-lag-window top-hat variance factor is an approximate noise-correlation diagnostic, not the exact covariance of the ICL or satellite estimators",
            "the normalized kernel uses delivered variance and does not identify source-photon noise or an effective gain",
            "these 57 target frames do not establish the covariance for all 800 cluster or matched-control frames",
            "no upstream sky-transfer or PSF-wing uncertainty is measured",
        ],
        "science_measurement": False,
    }
    summary_path = output_dir / "summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True, allow_nan=False) + "\n",
                            encoding="utf-8")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metadata", type=Path, default=DEFAULT_METADATA)
    parser.add_argument("--satellite-prep", type=Path, default=DEFAULT_PREP)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--patch-size", type=int, default=64)
    parser.add_argument("--max-patches", type=int, default=96)
    parser.add_argument("--max-lag", type=int, default=4)
    parser.add_argument("--min-valid-fraction", type=float, default=0.80)
    args = parser.parse_args()
    result = run(args.metadata.resolve(), args.satellite_prep.resolve(), args.output.resolve(),
                 patch_size=args.patch_size, max_patches=args.max_patches,
                 max_lag=args.max_lag, min_valid_fraction=args.min_valid_fraction)
    print(json.dumps({"status": result["status"], "science_measurement": False,
                      "n_frames": result["n_frames"], "by_filter": result["by_filter"],
                      "output": str(args.output.resolve())}, indent=2))


if __name__ == "__main__":
    main()
