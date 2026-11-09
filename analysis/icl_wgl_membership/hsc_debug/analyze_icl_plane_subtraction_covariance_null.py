#!/usr/bin/env python3
"""Compare measured long-lag kernels with ideal white noise after plane removal."""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
from estimate_hsc_coadd_noise_correlation import lag_correlations

ROOT = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
KERNEL_DIR = ROOT / "icl_supported_covariance_to_16px_2026-11-06_01"
APERTURE_DIR = ROOT / "icl_aperture_covariance_200_300_16px_2026-11-06_01"
OUTPUT = ROOT / "icl_plane_subtraction_covariance_null_2026-11-06_02"
PATCH = 64
MAX_LAG = 16
MC_REPLICATES = 100
MC_SEED = 413


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def read_lag_rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def affine_projection_kernel(size: int, max_lag: int) -> dict[tuple[int, int], float]:
    """Exact E[r_i r_j] lag means for unit white noise after unmasked OLS plane removal."""
    coordinate = (np.arange(size, dtype=np.float64) - (size - 1) / 2) / size
    n_pixel = size * size
    sum_x2 = size * float(np.sum(coordinate**2))
    sum_y2 = sum_x2
    result = {}
    for dy in range(-max_lag, max_lag + 1):
        y0 = coordinate[max(0, -dy):min(size, size-dy)]
        y1 = coordinate[max(0, -dy)+dy:min(size, size-dy)+dy]
        mean_y_product = float(np.mean(y0*y1))
        for dx in range(-max_lag, max_lag + 1):
            x0 = coordinate[max(0, -dx):min(size, size-dx)]
            x1 = coordinate[max(0, -dx)+dx:min(size, size-dx)+dx]
            mean_h = (1.0/n_pixel + float(np.mean(x0*x1))/sum_x2
                      + mean_y_product/sum_y2)
            value = -mean_h
            if dx == 0 and dy == 0:
                value += 1.0
            result[(dy, dx)] = value
    return result


def summarize_empirical(rows: list[dict[str, str]], maximum: int) -> dict:
    kernels: dict[tuple[str, str, str], dict[tuple[int, int], float]] = {}
    for row in rows:
        key = (str(row["cluster_id"]), str(row["role"]), str(row["filter"]))
        lag = (int(row["dy"]), int(row["dx"]))
        value = float(row["median_variance_normalized_covariance"])
        if not np.isfinite(value):
            raise ValueError(f"non-finite measured lag value: {key}/{lag}")
        kernels.setdefault(key, {})[lag] = value
    expected = {(dy, dx) for dy in range(-maximum, maximum+1)
                for dx in range(-maximum, maximum+1)}
    if len(kernels) != 448 or any(set(kernel) != expected for kernel in kernels.values()):
        raise ValueError("measured kernels do not cover all 448 frames and 33x33 lag cells")
    result = {}
    for band in ("i", "i2"):
        result[band] = {}
        for role in ("target", "matched_control"):
            selected = [kernel for (ident, this_role, this_band), kernel in kernels.items()
                        if this_band == band and this_role == role]
            if len(selected) != (71 if band == "i" else 153):
                raise ValueError(f"unexpected frame count for {band}/{role}")
            by_frame = []
            for kernel in selected:
                c8 = sum(value for (dy, dx), value in kernel.items()
                         if max(abs(dx), abs(dy)) <= 8)
                c16 = sum(kernel.values())
                by_frame.append((c8, c16, c16-c8))
            arr = np.asarray(by_frame)
            result[band][role] = {
                "n_frames": len(selected),
                "radius_8_median_q16_q84": [float(x) for x in
                    (np.quantile(arr[:, 0], .16), np.median(arr[:, 0]), np.quantile(arr[:, 0], .84))],
                "radius_16_median_q16_q84": [float(x) for x in
                    (np.quantile(arr[:, 1], .16), np.median(arr[:, 1]), np.quantile(arr[:, 1], .84))],
                "shell_9_to_16_median_q16_q84": [float(x) for x in
                    (np.quantile(arr[:, 2], .16), np.median(arr[:, 2]), np.quantile(arr[:, 2], .84))],
            }
    return result


def simulate_estimator_null(replicates: int, seed: int) -> dict:
    rng = np.random.default_rng(seed)
    variance = np.ones((PATCH, PATCH), dtype=np.float64)
    valid = np.ones((PATCH, PATCH), dtype=bool)
    draws = {radius: np.empty(replicates) for radius in (8, 16)}
    for index in range(replicates):
        image = rng.normal(size=(PATCH, PATCH))
        kernel = lag_correlations(image, variance, valid, max_lag=MAX_LAG)
        for radius in draws:
            draws[radius][index] = sum(
                value for (dy, dx), (value, count) in kernel.items()
                if max(abs(dx), abs(dy)) <= radius and count >= 32
            )
    return {
        str(radius): {"mean": float(np.mean(values)),
                      "q16_median_q84": [float(x) for x in
                          (np.quantile(values, .16), np.median(values), np.quantile(values, .84))]}
        for radius, values in draws.items()
    }


def main() -> None:
    if OUTPUT.exists():
        raise FileExistsError(f"refusing to overwrite {OUTPUT}")
    kernel_summary_path = KERNEL_DIR / "summary.json"
    kernel_progress_path = KERNEL_DIR / "progress.json"
    lag_path = KERNEL_DIR / "lag_correlation.tsv"
    aperture_summary_path = APERTURE_DIR / "summary.json"
    kernel_summary = read_json(kernel_summary_path)
    progress = read_json(kernel_progress_path)
    aperture_summary = read_json(aperture_summary_path)
    if (kernel_summary.get("status") != "complete_exploratory_diagnostic"
            or kernel_summary.get("science_measurement") is not False
            or kernel_summary.get("method", {}).get("max_lag_pixels") != MAX_LAG
            or progress.get("status") != "complete_exploratory_diagnostic"
            or progress.get("summary_sha256") != sha256(kernel_summary_path)
            or progress.get("output_sha256", {}).get("lag_correlation.tsv") != sha256(lag_path)
            or aperture_summary.get("science_measurement") is not False
            or "33x33" not in aperture_summary.get("method", "")
            or aperture_summary.get("long_lag_kernel_summary_sha256") != sha256(kernel_summary_path)):
        raise ValueError("kernel or aperture inputs are incomplete or changed")
    measured = summarize_empirical(read_lag_rows(lag_path), MAX_LAG)
    theoretical = affine_projection_kernel(PATCH, MAX_LAG)
    cumulative = {}
    previous = 0.0
    for radius in (0, 1, 2, 4, 8, 12, 16):
        value = sum(v for (dy, dx), v in theoretical.items()
                    if max(abs(dx), abs(dy)) <= radius)
        cumulative[str(radius)] = {"factor": value, "shell_from_previous_requested_radius": value-previous}
        previous = value
    summary = {
        "schema": "astrafek-icl-affine-plane-covariance-null-v1",
        "status": "complete_exploratory_diagnostic", "science_measurement": False,
        "model": {
            "noise": "independent unit-variance Gaussian pixels",
            "patch_shape": [PATCH, PATCH], "mask": "none",
            "fit": "unweighted ordinary least-squares affine plane, no clipping",
            "covariance": "exact expectation of (I-H) epsilon for the finite patch",
            "purpose": "estimate the long-lag anticorrelation induced by local affine-plane subtraction alone",
        },
        "idealized_affine_projection_cumulative_factors": cumulative,
        "same_estimator_white_noise_monte_carlo": {
            "replicates": MC_REPLICATES, "seed": MC_SEED,
            "cumulative_factors": simulate_estimator_null(MC_REPLICATES, MC_SEED),
        },
        "measured_frame_kernel_factors": measured,
        "limits": [
            "The idealized null has no masks, heteroscedastic variance, iterative clipping, astrophysical residuals, or real resampling covariance.",
            "Observed kernels use a median across patches at each lag; the resulting kernel need not be positive-semidefinite.",
            "This comparison diagnoses estimator behavior and does not estimate or correct the physical coadd covariance.",
        ],
        "input_sha256": {str(p.resolve()): sha256(p) for p in
                         (kernel_summary_path, kernel_progress_path, lag_path, aperture_summary_path)},
        "hsc_network_used": False, "data_requested": False,
    }
    OUTPUT.mkdir(parents=True)
    out = OUTPUT / "summary.json"
    out.write_text(json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True)+"\n",
                   encoding="utf-8")
    print(json.dumps({"output": str(out), "summary_sha256": sha256(out),
                      "idealized": cumulative, "measured": measured,
                      "science_measurement": False}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
