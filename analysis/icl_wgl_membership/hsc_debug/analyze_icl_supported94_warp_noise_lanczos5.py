#!/usr/bin/env python3
"""Propagate diagonal warp variance through Lanczos-5 cutout-to-coadd mapping.

This is an offline resampling sensitivity for the already approved warp set.
It uses the HSC PDR3 fifth-order Lanczos choice for the fractional WCS mapping,
but does not model correlations already present inside each resampled warp.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import tarfile
from collections import defaultdict
from pathlib import Path

import numpy as np
from astropy import units as u
from astropy.cosmology import FlatLambdaCDM
from astropy.io import fits
from astropy.wcs import WCS
from numba import njit

import analyze_icl_supported94_warp_noise as baseline


ROOT = baseline.ROOT
DOWNLOAD_DIR = baseline.DOWNLOAD_DIR
DEFAULT_OUTPUT = ROOT / "icl_supported94_warp_noise_lanczos5_2026-11-06_02"
KPC_BOUNDS = (200.0, 300.0)
LANCZOS_ORDER = 5


@njit(cache=False)
def lanczos_weight(distance: float) -> float:
    ad = abs(distance)
    if ad >= LANCZOS_ORDER:
        return 0.0
    if ad < 1e-14:
        return 1.0
    nearest = round(distance)
    if abs(distance - nearest) < 1e-13:
        return 0.0
    pix = math.pi * distance
    taper = pix / LANCZOS_ORDER
    return (math.sin(pix) / pix) * (math.sin(taper) / taper)


@njit(cache=False)
def lanczos5_variance_at_positions(image: np.ndarray, mask: np.ndarray,
                                   variance: np.ndarray, sx: np.ndarray,
                                   sy: np.ndarray, no_data_bit: np.int64,
                                   excluded_bits: np.int64):
    """Resample diagonal variance with squared separable Lanczos-5 weights.

    Pixels with any nonzero kernel coefficient outside the warp or on a
    no-data/nonfinite input are invalid. The clean flag additionally requires
    every contributing input mask to avoid the common exclusion bits.
    """
    n = sx.size
    out = np.full(n, np.nan, dtype=np.float64)
    valid = np.zeros(n, dtype=np.bool_)
    clean = np.zeros(n, dtype=np.bool_)
    for p in range(n):
        ix0 = math.floor(sx[p])
        iy0 = math.floor(sy[p])
        wx = np.zeros(11, dtype=np.float64)
        wy = np.zeros(11, dtype=np.float64)
        sum_x = 0.0
        sum_y = 0.0
        for j in range(11):
            ix = ix0 - 4 + j
            iy = iy0 - 4 + j
            wx[j] = lanczos_weight(sx[p] - ix)
            wy[j] = lanczos_weight(sy[p] - iy)
            sum_x += wx[j]
            sum_y += wy[j]
        norm = sum_x * sum_y
        if abs(norm) < 1e-12:
            continue
        variance_sum = 0.0
        is_valid = True
        is_clean = True
        for jy in range(11):
            if wy[jy] == 0.0:
                continue
            iy = iy0 - 4 + jy
            for jx in range(11):
                if wx[jx] == 0.0:
                    continue
                ix = ix0 - 4 + jx
                weight = wx[jx] * wy[jy]
                if ix < 0 or ix >= variance.shape[1] or iy < 0 or iy >= variance.shape[0]:
                    is_valid = False
                    break
                pixel_mask = mask[iy, ix]
                pixel_var = variance[iy, ix]
                if ((pixel_mask & no_data_bit) != 0
                        or not math.isfinite(pixel_var) or pixel_var <= 0.0
                        or not math.isfinite(image[iy, ix])):
                    is_valid = False
                    break
                if (pixel_mask & excluded_bits) != 0:
                    is_clean = False
                variance_sum += weight * weight * pixel_var
            if not is_valid:
                break
        if is_valid:
            value = variance_sum / (norm * norm)
            if math.isfinite(value) and value > 0.0:
                out[p] = value
                valid[p] = True
                clean[p] = is_clean
    return out, valid, clean


def fit_mapping(warp_wcs: WCS, coadd_wcs: WCS,
                shape: tuple[int, int]) -> tuple[np.ndarray, np.ndarray, float]:
    sy_grid, sx_grid = np.meshgrid(
        np.linspace(0, shape[0] - 1, 7),
        np.linspace(0, shape[1] - 1, 7), indexing="ij")
    ra, dec = warp_wcs.pixel_to_world_values(sx_grid, sy_grid)
    cx, cy = coadd_wcs.world_to_pixel_values(ra, dec)
    design = np.column_stack((np.ones(cx.size), sx_grid.ravel(), sy_grid.ravel()))
    fit_x = np.linalg.lstsq(design, cx.ravel(), rcond=None)[0]
    fit_y = np.linalg.lstsq(design, cy.ravel(), rcond=None)[0]
    residual = float(np.max(np.hypot(
        design @ fit_x - cx.ravel(), design @ fit_y - cy.ravel())))
    if residual > 0.05:
        raise ValueError(f"visit affine WCS residual exceeds 0.05 pixel: {residual}")
    return fit_x, fit_y, residual


def quantiles(values: np.ndarray) -> dict:
    x = np.asarray(values, dtype=np.float64)
    x = x[np.isfinite(x) & (x > 0)]
    if not x.size:
        return {"n_pixels": 0, "q16_median_q84": None,
                "fraction_within_1pct": None, "fraction_within_2pct": None}
    q16, median, q84 = np.quantile(x, [0.16, 0.5, 0.84])
    return {"n_pixels": int(x.size),
            "q16_median_q84": [float(q16), float(median), float(q84)],
            "fraction_within_1pct": float(np.mean(np.abs(x - 1) <= 0.01)),
            "fraction_within_2pct": float(np.mean(np.abs(x - 1) <= 0.02))}


def analyze_center(request: dict, info: dict, archive_path: Path,
                   record: dict) -> dict:
    if not archive_path.is_file() or baseline.sha256(archive_path) != record.get("archive_sha256"):
        raise ValueError(f"warp archive missing or hash-mismatched: {archive_path}")
    expected_visits = set(request["expected_visits"])
    if not expected_visits.issubset(set(record["audit"]["actual_visits"])):
        raise ValueError(f"metadata-expected visits are incomplete: {request['request_id']}")

    with fits.open(info["image_path"], memmap=True, checksum=True) as coadd:
        if len(coadd) < 4:
            raise ValueError("coadd lacks image/mask/variance planes")
        image = np.asarray(coadd[1].data, dtype=np.float32)
        coadd_mask = np.asarray(coadd[2].data, dtype=np.int64)
        coadd_variance = np.asarray(coadd[3].data, dtype=np.float32)
        if image.shape != coadd_mask.shape or image.shape != coadd_variance.shape:
            raise ValueError("coadd image/mask/variance shape mismatch")
        coadd_wcs = WCS(coadd[1].header).celestial
        center_x, center_y = baseline.scalar_pair(*coadd_wcs.world_to_pixel_values(
            info["ra_deg"], info["dec_deg"]))
        coadd_bits = baseline.mask_bits(coadd[2].header)
        if baseline.NO_DATA_NAME not in coadd_bits:
            raise ValueError("coadd mask lacks MP_NO_DATA")
        no_data_bit = np.int64(1 << coadd_bits[baseline.NO_DATA_NAME])
        coadd_available = (np.isfinite(image) & np.isfinite(coadd_variance)
                           & (coadd_variance > 0)
                           & ((coadd_mask & no_data_bit) == 0))
        common_names = set(baseline.MASK_EXCLUDE_NAMES) & set(coadd_bits)
        coadd_excluded_bits = np.int64(baseline.bitmask(sorted(common_names), coadd_bits))

        kpc_per_pixel = float(
            FlatLambdaCDM(H0=70.0, Om0=0.3, Tcmb0=0*u.K)
            .angular_diameter_distance(info["redshift"]).to_value(u.kpc)
            * 0.168 / (180.0 * 3600.0 / math.pi)
        )
        outer = KPC_BOUNDS[1]
        half = int(math.ceil(outer / kpc_per_pixel)) + 2
        x0, x1 = max(0, int(math.floor(center_x - half))), min(image.shape[1], int(math.ceil(center_x + half)) + 1)
        y0, y1 = max(0, int(math.floor(center_y - half))), min(image.shape[0], int(math.ceil(center_y + half)) + 1)
        yy, xx = np.mgrid[y0:y1, x0:x1]
        radius = np.hypot((xx + 0.5 - center_x) * kpc_per_pixel,
                          (yy + 0.5 - center_y) * kpc_per_pixel)
        annulus = (radius >= KPC_BOUNDS[0]) & (radius < KPC_BOUNDS[1])
        local_y, local_x = np.nonzero(annulus)
        gy, gx = local_y + y0, local_x + x0
        flat_n = len(gx)
        if flat_n == 0:
            raise ValueError("no geometric annulus pixels on the coadd")
        coadd_var = coadd_variance[gy, gx].astype(np.float64)
        coadd_ok = coadd_available[gy, gx]
        coadd_clean = coadd_ok & ((coadd_mask[gy, gx] & coadd_excluded_bits) == 0)
        inverse_variance_sum = np.zeros(flat_n, dtype=np.float64)
        n_valid = np.zeros(flat_n, dtype=np.uint16)
        n_clean = np.zeros(flat_n, dtype=np.uint16)
        actual_visits: set[int] = set()
        expected_visits_found: set[int] = set()
        mapping_residuals = []
        jacobians = []
        variance_support_loss = 0

        with tarfile.open(archive_path, mode="r:*") as archive:
            members = [m for m in archive.getmembers()
                       if m.isfile() and m.name.endswith(".fits")]
            for member in members:
                visit = int(Path(member.name).stem.removeprefix("warp-"))
                if visit in actual_visits:
                    raise ValueError(f"duplicate visit {visit}")
                actual_visits.add(visit)
                expected_visits_found.add(visit) if visit in expected_visits else None
                stream = archive.extractfile(member)
                if stream is None:
                    raise ValueError(f"cannot read {member.name}")
                with fits.open(stream, memmap=False) as warp:
                    if len(warp) < 4:
                        raise ValueError(f"visit {visit} lacks image/mask/variance")
                    w_image = np.asarray(warp[1].data, dtype=np.float32)
                    w_mask = np.asarray(warp[2].data, dtype=np.int64)
                    w_variance = np.asarray(warp[3].data, dtype=np.float32)
                    if w_image.shape != w_mask.shape or w_image.shape != w_variance.shape:
                        raise ValueError(f"visit {visit} plane shape mismatch")
                    warp_wcs = WCS(warp[1].header).celestial
                    fit_x, fit_y, residual = fit_mapping(warp_wcs, coadd_wcs, w_image.shape)
                    affine = np.array([[fit_x[1], fit_x[2]],
                                       [fit_y[1], fit_y[2]]], dtype=np.float64)
                    inverse = np.linalg.inv(affine)
                    source = inverse @ np.vstack((
                        gx.astype(np.float64) - fit_x[0],
                        gy.astype(np.float64) - fit_y[0]))
                    sx, sy = source[0], source[1]
                    w_bits = baseline.mask_bits(warp[2].header)
                    if w_bits.get(baseline.NO_DATA_NAME) != coadd_bits.get(baseline.NO_DATA_NAME):
                        raise ValueError(f"visit {visit} MP_NO_DATA bit differs from coadd")
                    warp_no_data = np.int64(1 << w_bits[baseline.NO_DATA_NAME])
                    names = tuple(name for name in baseline.MASK_EXCLUDE_NAMES if name in w_bits)
                    warp_excluded = np.int64(baseline.bitmask(names, w_bits))
                    local_var, valid, clean = lanczos5_variance_at_positions(
                        w_image, w_mask, w_variance, sx, sy, warp_no_data, warp_excluded)
                    inverse_variance_sum[valid] += 1.0 / local_var[valid]
                    n_valid[valid] += 1
                    n_clean[clean] += 1
                    nearest_x = np.rint(sx).astype(int)
                    nearest_y = np.rint(sy).astype(int)
                    in_nearest = ((nearest_x >= 0) & (nearest_x < w_image.shape[1])
                                  & (nearest_y >= 0) & (nearest_y < w_image.shape[0]))
                    variance_support_loss += int(np.count_nonzero(~valid & in_nearest))
                    mapping_residuals.append(residual)
                    jacobians.append(float(abs(np.linalg.det(affine))))

        if expected_visits_found != expected_visits:
            raise ValueError(f"expected visit set mismatch for {request['request_id']}")
        predicted = np.divide(1.0, inverse_variance_sum,
                              out=np.full(flat_n, np.nan), where=inverse_variance_sum > 0)
        ratio = np.divide(coadd_var, predicted,
                          out=np.full(flat_n, np.nan),
                          where=coadd_ok & np.isfinite(predicted) & (predicted > 0))
        primary = coadd_ok & (n_valid >= 3) & np.isfinite(ratio) & (ratio > 0)
        strict = coadd_clean & (n_clean >= 3) & np.isfinite(ratio) & (ratio > 0)
        geometric_area = math.pi * (KPC_BOUNDS[1]**2 - KPC_BOUNDS[0]**2) / kpc_per_pixel**2
        primary_q, strict_q = quantiles(ratio[primary]), quantiles(ratio[strict])
        return {
            "request_id": request["request_id"],
            "cluster_id": info["pair"]["cluster_id"],
            "filter": request["filter_code"], "role": info["role"],
            "connected_5deg_block": info["block"],
            "inner_kpc": KPC_BOUNDS[0], "outer_kpc": KPC_BOUNDS[1],
            "geometric_pixels_in_coadd": flat_n,
            "primary_pixels": int(primary.sum()),
            "primary_coverage_fraction_full_annulus": min(1.0, float(primary.sum() / geometric_area)),
            "primary_ratio": primary_q,
            "strict_pixels": int(strict.sum()),
            "strict_coverage_fraction_full_annulus": min(1.0, float(strict.sum() / geometric_area)),
            "strict_ratio": strict_q,
            "median_n_valid_visits": float(np.median(n_valid[primary])) if primary.any() else None,
            "median_n_clean_visits": float(np.median(n_clean[strict])) if strict.any() else None,
            "n_expected_visits": len(expected_visits),
            "n_returned_visits": len(actual_visits),
            "n_additional_visits": len(actual_visits - expected_visits),
            "max_affine_residual_pixel": float(max(mapping_residuals)),
            "median_warp_to_coadd_pixel_area_ratio": float(np.median(jacobians)),
            "kernel_complete_support_losses_inside_nearest_coverage": variance_support_loss,
            "variance_resampling": "separable Lanczos-5; diagonal input variance propagated with squared normalized kernel coefficients",
        }


def run(output_dir: Path, only_cluster: str | None = None) -> dict:
    if output_dir.exists():
        raise FileExistsError(f"refusing to overwrite analysis output: {output_dir}")
    plan, download, checkpoint, manifests, input_hashes = baseline.load_inputs()
    pairs = baseline.read_tsv(baseline.PAIR_TABLE)
    inventory = baseline.pair_inventory(plan, pairs, manifests)
    all_requests = list(plan["requests"])
    if only_cluster is not None:
        all_requests = [r for r in all_requests if str(only_cluster) in r["pair_ids"]]
        if not all_requests:
            raise ValueError(f"cluster {only_cluster} is not in the approved warp plan")
    rows, failures = [], []
    for request in all_requests:
        index = next(i for i, original in enumerate(plan["requests"])
                     if original["request_id"] == request["request_id"])
        record = checkpoint["requests"].get(f"request_{index:03d}", {})
        if record.get("status") == "http_404_permanent":
            failures.append({"request_id": request["request_id"], "status": "http_404_permanent"})
            continue
        archive_path = DOWNLOAD_DIR / f"request_{index:03d}.tar"
        try:
            rows.append(analyze_center(request, inventory[request["request_id"]],
                                       archive_path, record))
        except Exception as exc:
            failures.append({"request_id": request["request_id"],
                             "status": "analysis_error",
                             "error": f"{type(exc).__name__}: {exc}"})

    pair_bootstrap = baseline.paired_block_bootstrap(rows, pairs) if only_cluster is None else {}
    summary = {
        "schema": "astrafek-icl-supported94-warp-noise-lanczos5-v1",
        "status": "complete_exploratory" if not failures else "complete_with_data_gaps",
        "science_measurement": False, "hsc_network_used": False, "gpu_used": False,
        "parent_plan_sha256": baseline.PLAN_SHA256,
        "parent_download_summary_sha256": baseline.DOWNLOAD_SUMMARY_SHA256,
        "planned_pairs": 94, "planned_centers": 188,
        "requested_centers_this_run": len(all_requests),
        "analyzed_centers": len(rows), "failed_count": len(failures),
        "failures": failures,
        "radial_region_proper_kpc": list(KPC_BOUNDS),
        "kernel": {
            "name": "Lanczos-5",
            "official_hsc_pdr3_reference": "https://hsc-release.mtk.nao.ac.jp/doc/index.php/processing__pdr3/",
            "application": "secondary fractional-WCS mapping of downloaded warp cutouts onto paired coadd pixels",
            "variance_rule": "sum of squared normalized separable interpolation coefficients times input diagonal variances",
            "complete_support": "every nonzero kernel tap must be in-bounds, finite, positive-variance, and not MP_NO_DATA",
        },
        "assumptions_and_limits": [
            "warp-pixel errors are treated as diagonal for the secondary cutout-to-coadd mapping",
            "does not recover off-diagonal covariance already induced when HSC PDR3 created the warp from CCD pixels",
            "does not reproduce exact HSC coadd weights, clipping, flux-scale factors, or mask interpolation",
            "result is an approximation sensitivity, not an exact coadd-noise reconstruction or ICL correction",
        ],
        "paired_target_minus_control": pair_bootstrap,
        "input_sha256": input_hashes,
        "analysis_script_sha256": baseline.sha256(Path(__file__)),
    }
    output_dir.mkdir(parents=True, exist_ok=False)
    table = output_dir / "per_center_lanczos5.tsv"
    columns = list(rows[0]) if rows else ["request_id", "cluster_id", "filter", "role"]
    with table.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    summary["output_sha256"] = {table.name: baseline.sha256(table)}
    (output_dir / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8")
    print(json.dumps({
        "status": summary["status"], "science_measurement": False,
        "requested_centers": summary["requested_centers_this_run"],
        "analyzed_centers": summary["analyzed_centers"],
        "failed": summary["failed_count"], "output_dir": str(output_dir),
        "paired_target_minus_control": pair_bootstrap,
        "center_medians": [{"cluster": r["cluster_id"], "filter": r["filter"],
                            "role": r["role"],
                            "primary_ratio": (r["primary_ratio"]["q16_median_q84"][1]
                                              if r["primary_ratio"]["q16_median_q84"] else None),
                            "coverage": r["primary_coverage_fraction_full_annulus"]}
                           for r in rows],
    }, indent=2))
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--only-cluster", help="offline pilot on one already downloaded cluster")
    args = parser.parse_args()
    run(args.output, args.only_cluster)


if __name__ == "__main__":
    main()
