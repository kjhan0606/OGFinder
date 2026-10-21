#!/usr/bin/env python3
"""Measure whether local HSC PSF cores match the public PSF in an overlap annulus.

This evaluates a potential core-to-wing stitching region. It does not create a
hybrid PSF or approve the public model as a correction.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from astropy.io import fits


EXPECTED_PUBLIC_SHA256 = "8897af89767e5fbb8d5f9c70965126bc89409f6c65ecd8d117c7995f34be2d85"
EXPECTED_LOCAL_SUMMARY_SHA256 = "2a1a0f63f7cbabadaf6f1a99f27d53ea9b87d8df2232305c73a89e8ce75807ad"
EXPECTED_LOCAL_BATCH_AUDIT_SHA256 = "6375bc58df6d5a15d837cf5b13d3db79ebdc6944dbf7510164f465b6cbe72f7e"
EXPECTED_LOCAL_MODELS = 744
PIXEL_SCALE_ARCSEC = 0.168
DEFAULT_RADIAL_EDGES_ARCSEC = np.array([1.50, 1.80, 2.10, 2.40, 2.70, 3.00, 3.30])


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def radial_mean_profile(image: np.ndarray, radial_edges_arcsec=DEFAULT_RADIAL_EDGES_ARCSEC,
                        pixel_scale_arcsec=PIXEL_SCALE_ARCSEC):
    data = np.asarray(image, dtype=np.float64)
    edges = np.asarray(radial_edges_arcsec, dtype=float)
    if data.ndim != 2 or any(size % 2 != 1 for size in data.shape):
        raise ValueError("PSF core must be an odd-sized 2-D image")
    if not np.isfinite(data).all():
        raise ValueError("PSF core contains non-finite pixels")
    if edges.ndim != 1 or len(edges) < 3 or not np.isfinite(edges).all() or np.any(np.diff(edges) <= 0):
        raise ValueError("at least two increasing finite radial bins are required")
    y, x = np.indices(data.shape, dtype=float)
    radius = np.hypot(x - (data.shape[1] - 1) / 2, y - (data.shape[0] - 1) / 2)
    radius *= pixel_scale_arcsec
    profile, counts = [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        values = data[(radius >= lo) & (radius < hi)]
        counts.append(int(len(values)))
        profile.append(float(np.mean(values)) if len(values) else np.nan)
    return np.asarray(profile), np.asarray(counts)


def fit_profile_scale(local_profile: np.ndarray, public_profile: np.ndarray,
                      pixel_counts: np.ndarray) -> dict:
    local = np.asarray(local_profile, dtype=float)
    public = np.asarray(public_profile, dtype=float)
    counts = np.asarray(pixel_counts, dtype=float)
    valid = np.isfinite(local) & np.isfinite(public) & (counts > 0)
    if valid.sum() < 3:
        raise ValueError("fewer than three populated radial bins in the overlap region")
    local, public, counts = local[valid], public[valid], counts[valid]
    denominator = float(np.sum(counts * public * public))
    if denominator <= 0 or not math.isfinite(denominator):
        raise ValueError("public overlap profile has no positive weighted power")
    scale = float(np.sum(counts * local * public) / denominator)
    residual = local - scale * public
    local_power = float(np.sum(counts * local * local))
    if local_power <= 0 or not math.isfinite(local_power):
        raise ValueError("local overlap profile has no finite power")
    return {
        "amplitude_scale_local_over_public": scale,
        "weighted_fractional_shape_residual_rms": float(np.sqrt(np.sum(counts * residual * residual) / local_power)),
        "max_abs_residual_over_weighted_local_rms": float(
            np.max(np.abs(residual)) / np.sqrt(local_power / np.sum(counts))),
        "n_radial_bins_fit": int(valid.sum()),
        "n_pixels_used": int(counts.sum()),
        "local_profile_values": local.tolist(),
        "public_profile_values": public.tolist(),
    }


def _read_psf(path: Path, hdu: int) -> np.ndarray:
    with fits.open(path, memmap=True, checksum=True) as hdus:
        if hdu >= len(hdus) or hdus[hdu].data is None:
            raise ValueError(f"missing PSF image HDU {hdu}: {path}")
        image = np.asarray(hdus[hdu].data, dtype=np.float64)
    if image.ndim != 2 or any(size % 2 != 1 for size in image.shape) or not np.isfinite(image).all():
        raise ValueError(f"PSF image must be finite and odd-sized: {path}")
    return image


def analyze(public_path: Path, local_summary_path: Path, radial_edges_arcsec=DEFAULT_RADIAL_EDGES_ARCSEC):
    public_path, local_summary_path = Path(public_path).resolve(), Path(local_summary_path).resolve()
    if public_path.stat().st_size != 9809280 or sha256(public_path) != EXPECTED_PUBLIC_SHA256:
        raise ValueError("approved symmetric public PSF file mismatch")
    if sha256(local_summary_path) != EXPECTED_LOCAL_SUMMARY_SHA256:
        raise ValueError("approved local PSF summary checksum mismatch")
    local_summary = json.loads(local_summary_path.read_text(encoding="utf-8"))
    if (local_summary.get("schema") != "astrafek-hsc-position-psf-core-variation-v1"
            or local_summary.get("status") != "complete_exploratory_psf_core_audit"
            or local_summary.get("n_total_models_including_test") != 745
            or local_summary.get("science_measurement") is not False):
        raise ValueError("local PSF summary is incomplete or science-promoted")
    batch_audit = Path(local_summary["inputs"]["batch_audit"]).resolve()
    if sha256(batch_audit) != EXPECTED_LOCAL_BATCH_AUDIT_SHA256:
        raise ValueError("approved local PSF batch audit checksum mismatch")
    local_download_dir = batch_audit.parent
    local_records = [r for r in local_summary["records"]
                     if r.get("source") == "approved_744_position_batch"]
    if len(local_records) != EXPECTED_LOCAL_MODELS:
        raise ValueError("approved local PSF batch must contain exactly 744 records")

    public_full = _read_psf(public_path, hdu=1)
    if public_full.shape != (4001, 4001):
        raise ValueError("public PSF must be 4001x4001 pixels")
    public_core = public_full[1979:2022, 1979:2022]
    public_profile, counts = radial_mean_profile(public_core, radial_edges_arcsec)
    frames = []
    for row in local_records:
        local_path = local_download_dir / row["file"]
        if not local_path.is_file() or sha256(local_path) != row["sha256"]:
            raise ValueError(f"local PSF file missing or changed: {local_path.name}")
        local_image = _read_psf(local_path, hdu=0)
        local_profile, local_counts = radial_mean_profile(local_image, radial_edges_arcsec)
        if not np.array_equal(local_counts, counts):
            raise ValueError(f"overlap radial support differs from public PSF: {local_path.name}")
        fitted = fit_profile_scale(local_profile, public_profile, counts)
        role_names = sorted({str(role["role"]) for role in row.get("roles", [])})
        frame_filters = sorted({str(role.get("frame_filter", "unknown")) for role in row.get("roles", [])})
        frames.append({"file": local_path.name, "sha256": row["sha256"],
                       "roles": role_names, "frame_filters": frame_filters, **fitted})
    scales = np.asarray([r["amplitude_scale_local_over_public"] for r in frames])
    residuals = np.asarray([r["weighted_fractional_shape_residual_rms"] for r in frames])
    by_filter = {}
    for frame_filter in sorted({f for row in frames for f in row["frame_filters"]}):
        subset = [row for row in frames if frame_filter in row["frame_filters"]]
        by_filter[frame_filter] = {
            "n_positions": len(subset),
            "scale_median_q16_q84": np.quantile(
                [r["amplitude_scale_local_over_public"] for r in subset], [0.16, 0.5, 0.84]).tolist(),
            "shape_residual_median_q16_q84": np.quantile(
                [r["weighted_fractional_shape_residual_rms"] for r in subset], [0.16, 0.5, 0.84]).tolist(),
            "max_residual_over_local_rms_median_q16_q84": np.quantile(
                [r["max_abs_residual_over_weighted_local_rms"] for r in subset], [0.16, 0.5, 0.84]).tolist(),
        }
    return {
        "schema": "astrafek-hsc-public-psf-picker-core-overlap-audit-v1",
        "status": "complete_overlap_transfer_diagnostic",
        "public_psf": {"path": str(public_path), "sha256": sha256(public_path),
                       "overlap_radius_arcsec": [float(radial_edges_arcsec[0]), float(radial_edges_arcsec[-1])],
                       "radial_edges_arcsec": [float(v) for v in radial_edges_arcsec],
                       "annular_mean_flux_per_pixel": public_profile.tolist(),
                       "pixel_counts_per_annulus": counts.tolist()},
        "local_psf_provenance": {"summary_path": str(local_summary_path),
                                 "summary_sha256": EXPECTED_LOCAL_SUMMARY_SHA256,
                                 "batch_audit_sha256": EXPECTED_LOCAL_BATCH_AUDIT_SHA256,
                                 "n_positions": len(frames)},
        "fit_summary": {
            "scale_median_q16_q84": np.quantile(scales, [0.16, 0.5, 0.84]).tolist(),
            "scale_p05_p95": np.quantile(scales, [0.05, 0.95]).tolist(),
            "shape_residual_median_q16_q84": np.quantile(residuals, [0.16, 0.5, 0.84]).tolist(),
            "shape_residual_p05_p95": np.quantile(residuals, [0.05, 0.95]).tolist(),
            "max_residual_over_local_rms_median_q16_q84": np.quantile(
                [r["max_abs_residual_over_weighted_local_rms"] for r in frames], [0.16, 0.5, 0.84]).tolist(),
            "by_frame_filter": by_filter,
            "interpretation": "Amplitude and shape consistency in the overlap annulus are diagnostics; no hybrid PSF is constructed.",
        },
        "records": frames,
        "psf_approved": False,
        "science_measurement": False,
        "limitations": [
            "This compares model PSF cores, not observed star wings or cluster/control image transfer.",
            "A stable overlap fit does not validate the external wing amplitude beyond 3.3 arcsec.",
            "No local-core/public-wing hybrid is produced or applied.",
        ],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--public-psf", type=Path, required=True)
    parser.add_argument("--local-core-summary", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    output = args.output.resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite existing output: {output}")
    result = analyze(args.public_psf, args.local_core_summary)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({"status": result["status"], "n_positions": len(result["records"]),
                      "output": str(output), "sha256": sha256(output)}, indent=2))


if __name__ == "__main__":
    main()
