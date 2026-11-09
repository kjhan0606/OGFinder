#!/usr/bin/env python3
"""Compare HSC coadd VARIANCE with visit-warp variance over the supported set.

All network acquisition is complete before this offline analysis. Inputs are
hash-pinned; this runner does not access HSC, does not treat pixels as
independent draws, and never promotes science_measurement.
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


ROOT = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
PLAN_DIR = ROOT / "icl_supported94_warp_noise_plan_2026-11-06_01"
DOWNLOAD_DIR = ROOT / "icl_supported94_warp_noise_download_2026-11-06_01"
TARGET_MANIFEST = ROOT / "download_01/download_manifest.json"
CONTROL_MANIFEST = ROOT / "control_download_01/download_manifest.json"
PAIR_TABLE = ROOT / "icl_satellite_supported_94_association_2026-11-05_01/supported_pairs.tsv"
DEFAULT_OUTPUT = ROOT / "icl_supported94_warp_noise_analysis_2026-11-06_01"
PLAN_SHA256 = "bd936a94bcdbff83a607acbfca66616c396db526ada50cc89c9e3e312d37007a"
DOWNLOAD_SUMMARY_SHA256 = "0039321462f8158d1c006a1df64bb3e96ae988744ed780c845f74a60c79b95b4"
NO_DATA_NAME = "MP_NO_DATA"
MASK_EXCLUDE_NAMES = (
    "MP_BAD", "MP_SAT", "MP_INTRP", "MP_EDGE", "MP_CR", "MP_SUSPECT",
    "MP_NO_DATA", "MP_BRIGHT_OBJECT", "MP_CROSSTALK", "MP_UNMASKEDNAN",
    "MP_CLIPPED", "MP_REJECTED", "MP_SENSOR_EDGE", "MP_INEXACT_PSF",
)
RADIAL_BINS_KPC = ((0, 50), (50, 100), (100, 200), (200, 300), (300, 400))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def mask_bits(header) -> dict[str, int]:
    return {key: int(value) for key, value in header.items() if key.startswith("MP_")}


def bitmask(names, mapping: dict[str, int]) -> int:
    return sum(1 << mapping[name] for name in names)


def scalar_pair(x, y) -> tuple[float, float]:
    return float(np.asarray(x)), float(np.asarray(y))


def quantiles(values: np.ndarray) -> dict:
    values = np.asarray(values, dtype=np.float64)
    values = values[np.isfinite(values) & (values > 0)]
    if not values.size:
        return {"n_pixels": 0, "q16_median_q84": None,
                "fraction_within_1pct": None, "fraction_within_2pct": None}
    q16, median, q84 = np.quantile(values, [0.16, 0.5, 0.84])
    return {
        "n_pixels": int(values.size),
        "q16_median_q84": [float(q16), float(median), float(q84)],
        "fraction_within_1pct": float(np.mean(np.abs(values - 1) <= 0.01)),
        "fraction_within_2pct": float(np.mean(np.abs(values - 1) <= 0.02)),
    }


def load_inputs() -> tuple[dict, dict, dict, dict, dict]:
    paths = [PLAN_DIR / "plan.json", DOWNLOAD_DIR / "summary.json",
             DOWNLOAD_DIR / "checkpoint.json", TARGET_MANIFEST, CONTROL_MANIFEST,
             PAIR_TABLE]
    for path in paths:
        if not path.is_file():
            raise FileNotFoundError(path)
    hashes = {str(path): sha256(path) for path in paths}
    if hashes[str(PLAN_DIR / "plan.json")] != PLAN_SHA256:
        raise ValueError("approved plan hash differs from the fixed analysis scope")
    if hashes[str(DOWNLOAD_DIR / "summary.json")] != DOWNLOAD_SUMMARY_SHA256:
        raise ValueError("download summary changed after the acquisition audit")
    plan = json.loads((PLAN_DIR / "plan.json").read_text(encoding="utf-8"))
    summary = json.loads((DOWNLOAD_DIR / "summary.json").read_text(encoding="utf-8"))
    checkpoint = json.loads((DOWNLOAD_DIR / "checkpoint.json").read_text(encoding="utf-8"))
    target_manifest = json.loads(TARGET_MANIFEST.read_text(encoding="utf-8"))
    control_manifest = json.loads(CONTROL_MANIFEST.read_text(encoding="utf-8"))
    if (plan.get("status") != "prepared_not_submitted"
            or plan.get("query_submitted") is not False
            or plan.get("science_measurement") is not False
            or plan.get("scope", {}).get("n_pairs") != 94
            or summary.get("plan_sha256") != PLAN_SHA256
            or summary.get("status") != "complete_with_data_gaps"
            or summary.get("science_measurement") is not False
            or summary.get("completed_archive_count") != 184
            or summary.get("http_404_request_count") != 4
            or checkpoint.get("plan_sha256") != PLAN_SHA256
            or checkpoint.get("science_measurement") is not False):
        raise ValueError("plan/download status or research guard state is inconsistent")
    return plan, summary, checkpoint, {"target": target_manifest, "control": control_manifest}, hashes


def pair_inventory(plan: dict, pair_rows: list[dict], manifests: dict) -> dict[str, dict]:
    targets = {str(row["id"]): row for row in manifests["target"]["clusters"]}
    controls = {str(row["id"]): row for row in manifests["control"]["clusters"]}
    by_request = {}
    for pair in pair_rows:
        band = pair["filter"]
        cid = str(pair["cluster_id"])
        control_id = str(pair["control_id"])
        target_request = f"{band}:target:{cid}"
        control_request = f"{band}:control:{control_id}"
        if cid not in targets or control_id not in controls:
            raise ValueError(f"paired coadd manifest is missing cluster/control {cid}")
        for request_id, role, row, ra_key, dec_key in (
            (target_request, "target", targets[cid], "target_ra_deg", "target_dec_deg"),
            (control_request, "control", controls[control_id], "control_ra_deg", "control_dec_deg"),
        ):
            pair_value = float(pair[ra_key]), float(pair[dec_key])
            manifest_value = float(row["ra_deg"]), float(row["dec_deg"])
            if (abs(pair_value[0] - manifest_value[0]) > 1e-6
                    or abs(pair_value[1] - manifest_value[1]) > 1e-6
                    or row.get("science_hdu") != 1
                    or row.get("band") != "HSC-I"):
                raise ValueError(f"paired/coadd manifest geometry or band mismatch: {request_id}")
            image_dir = ROOT / ("download_01/images" if role == "target" else "control_download_01/images")
            image_path = image_dir / row["image"]
            if not image_path.is_file():
                raise FileNotFoundError(image_path)
            by_request[request_id] = {
                "pair": pair,
                "role": role,
                "image_path": image_path,
                "ra_deg": pair_value[0],
                "dec_deg": pair_value[1],
                "redshift": float(pair["redshift"]),
                "block": pair["connected_5deg_block"],
            }
    expected_ids = {row["request_id"] for row in plan["requests"]}
    if set(by_request) != expected_ids or len(by_request) != 188:
        raise ValueError("pair/coadd inventory does not exactly cover the approved request list")
    return by_request


def analyze_center(request: dict, info: dict, archive_path: Path,
                    record: dict, plan: dict) -> list[dict]:
    if not archive_path.is_file() or sha256(archive_path) != record.get("archive_sha256"):
        raise ValueError(f"warp archive is missing or hash-mismatched: {archive_path}")
    expected_visits = set(request["expected_visits"])
    returned_visits = set(record["audit"]["actual_visits"])
    if not expected_visits.issubset(returned_visits):
        raise ValueError(f"one or more expected visits are absent in {request['request_id']}")

    rows = []
    coadd_path = info["image_path"]
    with fits.open(coadd_path, memmap=True) as coadd:
        if len(coadd) < 4 or str(coadd[0].header.get("FILTER", "")).lower() not in ("i", "i2"):
            raise ValueError(f"unexpected target/control coadd product: {coadd_path}")
        image = np.asarray(coadd[1].data)
        coadd_mask = np.asarray(coadd[2].data)
        coadd_variance = np.asarray(coadd[3].data)
        if (image.ndim != 2 or image.shape != coadd_mask.shape
                or image.shape != coadd_variance.shape):
            raise ValueError("coadd IMAGE/MASK/VARIANCE arrays are not aligned")
        coadd_wcs = WCS(coadd[1].header).celestial
        center_x, center_y = scalar_pair(*coadd_wcs.world_to_pixel_values(
            info["ra_deg"], info["dec_deg"]))
        coadd_bits = mask_bits(coadd[2].header)
        if NO_DATA_NAME not in coadd_bits:
            raise ValueError("coadd MASK lacks MP_NO_DATA")
        coadd_no_data = 1 << coadd_bits[NO_DATA_NAME]
        coadd_available = (np.isfinite(image) & np.isfinite(coadd_variance)
                           & (coadd_variance > 0) & ((coadd_mask & coadd_no_data) == 0))

        variance_sum = np.zeros(image.shape, dtype=np.float64)
        n_valid = np.zeros(image.shape, dtype=np.uint16)
        n_clean = np.zeros(image.shape, dtype=np.uint16)
        footprint = np.zeros(image.shape, dtype=bool)
        offsets = []
        visit_mask_maps = []
        per_visit = []
        with tarfile.open(archive_path, mode="r:*") as archive:
            members = [item for item in archive.getmembers()
                       if item.isfile() and item.name.endswith(".fits")]
            found = set()
            for member in members:
                visit = int(Path(member.name).stem.removeprefix("warp-"))
                if visit in found:
                    raise ValueError(f"duplicate visit FITS in {request['request_id']}: {visit}")
                found.add(visit)
                stream = archive.extractfile(member)
                if stream is None:
                    raise ValueError(f"cannot read {member.name}")
                with fits.open(stream, memmap=False) as warp:
                    if len(warp) < 4:
                        raise ValueError(f"visit {visit} lacks IMAGE/MASK/VARIANCE")
                    w_image = np.asarray(warp[1].data)
                    w_mask = np.asarray(warp[2].data)
                    w_variance = np.asarray(warp[3].data)
                    if (w_image.ndim != 2 or w_image.shape != w_mask.shape
                            or w_image.shape != w_variance.shape):
                        raise ValueError(f"visit {visit} planes are misaligned")
                    warp_wcs = WCS(warp[1].header).celestial
                    # Warp cutouts have their own fractional-pixel WCS offsets.
                    # Fit the smooth warp-to-coadd mapping on a 7x7 grid and
                    # nearest-sample it onto coadd pixel centers.  The measured
                    # affine residual is retained as a QA field; this is a
                    # variance consistency diagnostic, not HSC's exact resampler.
                    sx_grid, sy_grid = np.meshgrid(
                        np.linspace(0, w_image.shape[1]-1, 7),
                        np.linspace(0, w_image.shape[0]-1, 7))
                    wx, wy = warp_wcs.pixel_to_world_values(sx_grid, sy_grid)
                    cx_grid, cy_grid = coadd_wcs.world_to_pixel_values(wx, wy)
                    design = np.column_stack((np.ones(cx_grid.size),
                                              sx_grid.ravel(), sy_grid.ravel()))
                    fit_x = np.linalg.lstsq(design, cx_grid.ravel(), rcond=None)[0]
                    fit_y = np.linalg.lstsq(design, cy_grid.ravel(), rcond=None)[0]
                    residual = float(np.max(np.hypot(
                        design @ fit_x - cx_grid.ravel(),
                        design @ fit_y - cy_grid.ravel())))
                    if residual > 0.05:
                        raise ValueError(f"visit {visit} affine WCS residual exceeds 0.05 pixel: {residual}")
                    corners = np.array([[1, 0, 0], [1, w_image.shape[1]-1, 0],
                                        [1, 0, w_image.shape[0]-1],
                                        [1, w_image.shape[1]-1, w_image.shape[0]-1]])
                    corner_x, corner_y = corners @ fit_x, corners @ fit_y
                    x0 = max(0, int(np.floor(corner_x.min())))
                    y0 = max(0, int(np.floor(corner_y.min())))
                    x1 = min(image.shape[1], int(np.ceil(corner_x.max())) + 1)
                    y1 = min(image.shape[0], int(np.ceil(corner_y.max())) + 1)
                    if x1 <= x0 or y1 <= y0:
                        raise ValueError(f"visit {visit} footprint misses paired coadd")
                    yy_local, xx_local = np.mgrid[y0:y1, x0:x1]
                    inverse = np.linalg.inv(np.array([[fit_x[1], fit_x[2]],
                                                       [fit_y[1], fit_y[2]]]))
                    source = inverse @ np.vstack((
                        xx_local.ravel() - fit_x[0], yy_local.ravel() - fit_y[0]))
                    sx = np.rint(source[0]).astype(int).reshape(xx_local.shape)
                    sy = np.rint(source[1]).astype(int).reshape(yy_local.shape)
                    inside = ((source[0].reshape(xx_local.shape) >= -0.5)
                              & (source[0].reshape(xx_local.shape) < w_image.shape[1]-0.5)
                              & (source[1].reshape(yy_local.shape) >= -0.5)
                              & (source[1].reshape(yy_local.shape) < w_image.shape[0]-0.5))
                    sx = np.clip(sx, 0, w_image.shape[1]-1)
                    sy = np.clip(sy, 0, w_image.shape[0]-1)
                    w_bits = mask_bits(warp[2].header)
                    if w_bits.get(NO_DATA_NAME) != coadd_bits.get(NO_DATA_NAME):
                        raise ValueError(f"visit {visit} MP_NO_DATA plane bit differs from coadd")
                    no_data = 1 << w_bits[NO_DATA_NAME]
                    excludes = tuple(name for name in MASK_EXCLUDE_NAMES if name in w_bits)
                    excluded = bitmask(excludes, w_bits)
                    sample_image = w_image[sy, sx]
                    sample_mask = w_mask[sy, sx]
                    sample_variance = w_variance[sy, sx]
                    valid = (inside & np.isfinite(sample_image) & np.isfinite(sample_variance)
                             & (sample_variance > 0) & ((sample_mask & no_data) == 0))
                    clean = valid & ((sample_mask & excluded) == 0)
                    target_slice = np.s_[y0:y1, x0:x1]
                    local_inverse = np.zeros(valid.shape, dtype=np.float64)
                    local_inverse[valid] = 1.0 / sample_variance[valid]
                    variance_sum[target_slice] += local_inverse
                    n_valid[target_slice] += valid.astype(np.uint16)
                    n_clean[target_slice] += clean.astype(np.uint16)
                    footprint[target_slice] |= inside
                    offsets.append([float(fit_x[0]), float(fit_y[0])])
                    visit_mask_maps.append(set(w_bits))
                    per_visit.append({
                        "visit": visit,
                        "is_metadata_expected": visit in expected_visits,
                        "member_bytes": int(member.size),
                        "valid_image_variance_fraction": float(valid.sum() / max(inside.sum(), 1)),
                        "clean_pixel_fraction": float(clean.sum() / max(inside.sum(), 1)),
                        "median_positive_variance": float(np.median(sample_variance[valid])) if valid.any() else None,
                        "coadd_affine_offset_xy": [float(fit_x[0]), float(fit_y[0])],
                        "wcs_affine_grid_max_residual_pixel": residual,
                        "variance_sampling": "nearest_neighbor_on_coadd_grid",
                    })
        if not expected_visits.issubset(found):
            raise ValueError(f"expected visits are missing in {request['request_id']}")
        predicted = np.divide(1.0, variance_sum, out=np.full(image.shape, np.nan),
                              where=variance_sum > 0)
        ratio = np.divide(coadd_variance, predicted, out=np.full(image.shape, np.nan),
                          where=coadd_available & np.isfinite(predicted) & (predicted > 0))

        common_names = set(MASK_EXCLUDE_NAMES) & set(coadd_bits)
        for present in visit_mask_maps:
            common_names.intersection_update(present)
        coadd_excluded = bitmask(sorted(common_names), coadd_bits)
        clean_coadd = coadd_available & ((coadd_mask & coadd_excluded) == 0)

        kpc_per_pixel = float(
            FlatLambdaCDM(H0=70.0, Om0=0.3, Tcmb0=0*u.K)
            .angular_diameter_distance(info["redshift"]).to_value(u.kpc)
            * 0.168 / (180.0 * 3600.0 / math.pi)
        )
        yy, xx = np.indices(image.shape)
        radius_kpc = np.hypot((xx + 0.5 - center_x) * kpc_per_pixel,
                              (yy + 0.5 - center_y) * kpc_per_pixel)
        output = []
        for inner, outer in RADIAL_BINS_KPC:
            annulus = (radius_kpc >= inner) & (radius_kpc < outer)
            annulus_footprint = annulus & footprint
            primary = annulus & footprint & coadd_available & (n_valid >= 3) & np.isfinite(ratio) & (ratio > 0)
            strict = annulus & footprint & clean_coadd & (n_clean >= 3) & np.isfinite(ratio) & (ratio > 0)
            n_geometry = int(annulus.sum())
            n_in_footprint = int(annulus_footprint.sum())
            output.append({
                "request_id": request["request_id"],
                "cluster_id": info["pair"]["cluster_id"],
                "filter": request["filter_code"],
                "role": info["role"],
                "connected_5deg_block": info["block"],
                "inner_kpc": inner,
                "outer_kpc": outer,
                "geometric_pixels": n_geometry,
                "footprint_pixels": n_in_footprint,
                "footprint_coverage_fraction": float(n_in_footprint / n_geometry) if n_geometry else 0.0,
                "primary_pixels": int(primary.sum()),
                "primary_coverage_fraction": float(primary.sum() / n_geometry) if n_geometry else 0.0,
                "primary_coverage_within_footprint": float(primary.sum() / n_in_footprint) if n_in_footprint else 0.0,
                "median_n_valid_visits": float(np.median(n_valid[primary])) if primary.any() else None,
                "primary_ratio": quantiles(ratio[primary]),
                "strict_pixels": int(strict.sum()),
                "strict_coverage_fraction": float(strict.sum() / n_geometry) if n_geometry else 0.0,
                "strict_coverage_within_footprint": float(strict.sum() / n_in_footprint) if n_in_footprint else 0.0,
                "strict_median_n_clean_visits": float(np.median(n_clean[strict])) if strict.any() else None,
                "strict_ratio": quantiles(ratio[strict]),
                "common_mask_exclusions": sorted(common_names),
                "n_metadata_expected_visits": len(expected_visits),
                "n_returned_visits": len(found),
                "n_additional_returned_visits": len(found - expected_visits),
                "n_valid_visits_center_median": float(np.median(n_valid[footprint])) if footprint.any() else None,
                "n_clean_visits_center_median": float(np.median(n_clean[footprint])) if footprint.any() else None,
                "valid_visit_count_footprint_fraction": float(np.mean(n_valid[footprint] >= 3)) if footprint.any() else 0.0,
                "offsets_xy_unique": len({tuple(item) for item in offsets}),
                "warp_footprint_pixels": int(footprint.sum()),
                "warp_archive_bytes": int(archive_path.stat().st_size),
            })
        return output


def paired_block_bootstrap(center_rows: list[dict], pair_rows: list[dict],
                           *, n_replicates: int = 5000, seed: int = 20261109) -> dict:
    by_pair = defaultdict(dict)
    for row in center_rows:
        if row["inner_kpc"] == 200 and row["outer_kpc"] == 300:
            value = row["primary_ratio"]["q16_median_q84"]
            if value:
                by_pair[(row["filter"], row["cluster_id"])][row["role"]] = value[1]
    pair_block = {(row["filter"], str(row["cluster_id"])): row["connected_5deg_block"]
                  for row in pair_rows}
    results = {}
    for band in ("i", "i2"):
        band_seed = seed + (0 if band == "i" else 1)
        rng = np.random.default_rng(band_seed)
        paired = [(cid, roles["target"] - roles["control"],
                   pair_block[(filter_code, cid)])
                  for (filter_code, cid), roles in by_pair.items()
                  if filter_code == band and {"target", "control"}.issubset(roles)]
        blocks = sorted({item[2] for item in paired})
        deltas = np.asarray([item[1] for item in paired], dtype=float)
        if not deltas.size:
            results[band] = {"n_complete_pairs": 0, "median_target_minus_control": None,
                             "block_bootstrap_ci95": None}
            continue
        index_by_block = {block: np.asarray([i for i, item in enumerate(paired) if item[2] == block], dtype=int)
                          for block in blocks}
        estimates = np.empty(n_replicates, dtype=float)
        for replica in range(n_replicates):
            chosen = rng.choice(blocks, size=len(blocks), replace=True)
            selected = np.concatenate([index_by_block[block] for block in chosen])
            estimates[replica] = float(np.median(deltas[selected]))
        results[band] = {
            "n_complete_pairs": int(len(paired)),
            "n_connected_blocks": int(len(blocks)),
            "median_target_minus_control": float(np.median(deltas)),
            "block_bootstrap_ci95": [float(x) for x in np.quantile(estimates, [0.025, 0.975])],
            "replicates": n_replicates,
            "seed": band_seed,
        }
    return results


def run(output_dir: Path) -> dict:
    output_dir = Path(output_dir).resolve()
    if output_dir.exists():
        raise FileExistsError(f"refusing to overwrite analysis output: {output_dir}")
    plan, download, checkpoint, manifests, input_hashes = load_inputs()
    pairs = read_tsv(PAIR_TABLE)
    inventory = pair_inventory(plan, pairs, manifests)
    records = checkpoint["requests"]
    centers = []
    failures = []
    archive_index = {request["request_id"]: index for index, request in enumerate(plan["requests"])}
    for request in plan["requests"]:
        request_id = request["request_id"]
        index = archive_index[request_id]
        record = records.get(f"request_{index:03d}", {})
        if record.get("status") == "http_404_permanent":
            failures.append({"request_id": request_id, "cluster_id": request.get("pair_ids", [None])[0],
                             "filter": request["filter_code"], "role": request["role"],
                             "status": "http_404_permanent"})
            continue
        archive_path = DOWNLOAD_DIR / f"request_{index:03d}.tar"
        try:
            centers.extend(analyze_center(request, inventory[request_id], archive_path, record, plan))
        except Exception as exc:
            failures.append({"request_id": request_id,
                             "cluster_id": request.get("pair_ids", [None])[0],
                             "filter": request["filter_code"], "role": request["role"],
                             "status": "analysis_error", "error": f"{type(exc).__name__}: {exc}"})

    radial_200_300 = [row for row in centers if row["inner_kpc"] == 200 and row["outer_kpc"] == 300]
    pair_comparison = paired_block_bootstrap(centers, pairs)
    summary = {
        "schema": "astrafek-icl-supported94-warp-noise-analysis-v1",
        "status": "complete_exploratory" if not failures else "complete_with_data_gaps",
        "science_measurement": False,
        "hsc_network_used": False,
        "gpu_used": False,
        "parent_plan_sha256": PLAN_SHA256,
        "parent_download_summary_sha256": DOWNLOAD_SUMMARY_SHA256,
        "planned_pair_count": 94,
        "planned_center_count": 188,
        "analyzed_center_count": len({row["request_id"] for row in centers}),
        "failed_center_count": len(failures),
        "failed_centers": failures,
        "radial_200_300_center_count": len(radial_200_300),
        "radial_200_300_full_annulus_coverage_median": float(np.median([r["primary_coverage_fraction"] for r in radial_200_300])) if radial_200_300 else None,
        "radial_200_300_primary_ratio_median_across_centers": float(np.median([
            r["primary_ratio"]["q16_median_q84"][1] for r in radial_200_300
            if r["primary_ratio"]["q16_median_q84"]])) if any(r["primary_ratio"]["q16_median_q84"] for r in radial_200_300) else None,
        "radial_200_300_strict_coverage_median": float(np.median([r["strict_coverage_fraction"] for r in radial_200_300])) if radial_200_300 else None,
        "paired_target_minus_control_primary_median_ratio": pair_comparison,
        "interpretation": [
            "The pixelwise coadd/predicted variance ratio is a product-consistency diagnostic, not an independent detector-noise calibration.",
            "Pixel quantiles are descriptive; pixels are spatially correlated and are not used as independent inferential replicates.",
            "A 2-arcmin square does not fully cover every low-redshift 200-300 kpc annulus; report per-center geometric support.",
            "This analysis does not establish source-Poisson/effective gain, resampling covariance, residual-sky covariance, PSF-wing uncertainty, or total ICL-estimator covariance.",
            "Four approved center requests returned HTTP 404; center-level data gaps remain visible in the paired sample.",
        ],
        "input_sha256": input_hashes,
    }
    output_dir.mkdir(parents=True, exist_ok=False)
    radial_path = output_dir / "per_center_radial.tsv"
    columns = list(centers[0]) if centers else ["request_id", "cluster_id", "filter", "role", "inner_kpc", "outer_kpc"]
    with radial_path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(centers)
    summary_path = output_dir / "summary.json"
    summary["per_center_radial_sha256"] = sha256(radial_path)
    summary_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2,
                                       sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({"status": summary["status"],
                      "science_measurement": False,
                      "analyzed_centers": summary["analyzed_center_count"],
                      "failed_centers": summary["failed_center_count"],
                      "summary": str(summary_path),
                      "per_center_radial": str(radial_path)}, ensure_ascii=False, indent=2))
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    run(args.output_dir)


if __name__ == "__main__":
    main()
