#!/usr/bin/env python3
"""Pilot a same-visit split-stack noise check for one supported ICL pair.

This offline runner uses only the already approved HSC-I warp cutouts for
cluster 101 and its matched control. It does not reproduce the HSC coadd
resampler or promote any ICL measurement.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
import tarfile
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
DEFAULT_OUTPUT = ROOT / "icl_visit_split_covariance_pilot_2026-11-06_03"
PLAN_SHA256 = "bd936a94bcdbff83a607acbfca66616c396db526ada50cc89c9e3e312d37007a"
DOWNLOAD_SHA256 = "0039321462f8158d1c006a1df64bb3e96ae988744ed780c845f74a60c79b95b4"
EXCLUDE_NAMES = (
    "MP_BAD", "MP_SAT", "MP_INTRP", "MP_EDGE", "MP_CR", "MP_SUSPECT",
    "MP_NO_DATA", "MP_BRIGHT_OBJECT", "MP_CROSSTALK", "MP_UNMASKEDNAN",
    "MP_CLIPPED", "MP_REJECTED", "MP_SENSOR_EDGE", "MP_INEXACT_PSF",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def finite_list(array) -> list:
    values = np.asarray(array, dtype=np.float64)
    return [[float(x) if np.isfinite(x) else None for x in row] for row in values]


def mask_bits(header) -> dict[str, int]:
    return {key: int(value) for key, value in header.items() if key.startswith("MP_")}


def build_visit_map(warp_wcs, coadd_wcs, warp_shape, coadd_shape):
    sx_grid, sy_grid = np.meshgrid(np.linspace(0, warp_shape[1] - 1, 7),
                                   np.linspace(0, warp_shape[0] - 1, 7))
    ra, dec = warp_wcs.pixel_to_world_values(sx_grid, sy_grid)
    cx, cy = coadd_wcs.world_to_pixel_values(ra, dec)
    design = np.column_stack((np.ones(cx.size), sx_grid.ravel(), sy_grid.ravel()))
    fit_x = np.linalg.lstsq(design, cx.ravel(), rcond=None)[0]
    fit_y = np.linalg.lstsq(design, cy.ravel(), rcond=None)[0]
    residual = float(np.max(np.hypot(design @ fit_x - cx.ravel(),
                                     design @ fit_y - cy.ravel())))
    if residual > 0.05:
        raise ValueError(f"affine warp-to-coadd residual {residual:.4f} pixel exceeds 0.05")
    corners = np.array([[1, 0, 0], [1, warp_shape[1]-1, 0],
                        [1, 0, warp_shape[0]-1], [1, warp_shape[1]-1, warp_shape[0]-1]])
    corner_x, corner_y = corners @ fit_x, corners @ fit_y
    x0 = max(0, int(np.floor(corner_x.min())))
    y0 = max(0, int(np.floor(corner_y.min())))
    x1 = min(coadd_shape[1], int(np.ceil(corner_x.max())) + 1)
    y1 = min(coadd_shape[0], int(np.ceil(corner_y.max())) + 1)
    yy, xx = np.mgrid[y0:y1, x0:x1]
    inverse = np.linalg.inv(np.array([[fit_x[1], fit_x[2]], [fit_y[1], fit_y[2]]]))
    source = inverse @ np.vstack((xx.ravel() - fit_x[0], yy.ravel() - fit_y[0]))
    sx_float = source[0].reshape(xx.shape)
    sy_float = source[1].reshape(yy.shape)
    inside = ((sx_float >= -0.5) & (sx_float < warp_shape[1] - 0.5)
              & (sy_float >= -0.5) & (sy_float < warp_shape[0] - 0.5))
    sx = np.clip(np.rint(sx_float).astype(int), 0, warp_shape[1] - 1)
    sy = np.clip(np.rint(sy_float).astype(int), 0, warp_shape[0] - 1)
    return (x0, y0, x1, y1, sx, sy, inside, residual)


def center_covariance(request, request_index, pair, role, coadd_row, download_record,
                      output_dir):
    archive_path = DOWNLOAD_DIR / f"request_{request_index:03d}.tar"
    if sha256(archive_path) != download_record["archive_sha256"]:
        raise ValueError(f"archive hash mismatch for {request['request_id']}")
    coadd_path = ROOT / ("download_01/images" if role == "target" else "control_download_01/images") / coadd_row["image"]
    with fits.open(coadd_path, memmap=True) as coadd:
        image = np.asarray(coadd[1].data)
        coadd_mask = np.asarray(coadd[2].data)
        coadd_var = np.asarray(coadd[3].data)
        wcs = WCS(coadd[1].header).celestial
        ra = float(pair["target_ra_deg"] if role == "target" else pair["control_ra_deg"])
        dec = float(pair["target_dec_deg"] if role == "target" else pair["control_dec_deg"])
        cx, cy = (float(np.asarray(v)) for v in wcs.world_to_pixel_values(ra, dec))
        c_bits = mask_bits(coadd[2].header)
        if "MP_NO_DATA" not in c_bits:
            raise ValueError("coadd mask lacks MP_NO_DATA")
        c_no_data = 1 << c_bits["MP_NO_DATA"]
        excluded = sum(1 << c_bits[name] for name in EXCLUDE_NAMES if name in c_bits)
        c_valid = (np.isfinite(image) & np.isfinite(coadd_var) & (coadd_var > 0)
                   & ((coadd_mask & c_no_data) == 0) & ((coadd_mask & excluded) == 0))
        s_a = np.zeros(image.shape, dtype=np.float64)
        w_a = np.zeros(image.shape, dtype=np.float64)
        s_b = np.zeros(image.shape, dtype=np.float64)
        w_b = np.zeros(image.shape, dtype=np.float64)
        visit_ids = sorted(request["expected_visits"])
        split_a = set(visit_ids[::2])
        split_b = set(visit_ids[1::2])
        # Sorted alternating visits give each half repeated coverage across the
        # same observing seasons, while target/control use the identical split.
        wcs_residuals = []
        split_visit_valid = {"A": [], "B": []}
        with tarfile.open(archive_path, mode="r:*") as archive:
            for member in archive.getmembers():
                if not member.isfile() or not member.name.endswith(".fits"):
                    continue
                visit = int(Path(member.name).stem.removeprefix("warp-"))
                if visit not in split_a | split_b:
                    continue
                stream = archive.extractfile(member)
                if stream is None:
                    raise ValueError(f"cannot extract {member.name}")
                with fits.open(stream, memmap=False) as hdus:
                    wi, wm, wv = (np.asarray(hdus[i].data) for i in (1, 2, 3))
                    warp_wcs = WCS(hdus[1].header).celestial
                    x0, y0, x1, y1, sx, sy, inside, residual = build_visit_map(
                        warp_wcs, wcs, wi.shape, image.shape)
                    wcs_residuals.append(residual)
                    bits = mask_bits(hdus[2].header)
                    if bits.get("MP_NO_DATA") != c_bits.get("MP_NO_DATA"):
                        raise ValueError("warp and coadd MP_NO_DATA bit differ")
                    invalid = sum(1 << bits[name] for name in EXCLUDE_NAMES if name in bits)
                    sample_i, sample_m, sample_v = wi[sy, sx], wm[sy, sx], wv[sy, sx]
                    valid = (inside & np.isfinite(sample_i) & np.isfinite(sample_v)
                             & (sample_v > 0) & ((sample_m & (1 << bits["MP_NO_DATA"])) == 0)
                             & ((sample_m & invalid) == 0))
                    local_w = np.zeros(valid.shape, dtype=np.float64)
                    local_w[valid] = 1.0 / sample_v[valid]
                    target = np.s_[y0:y1, x0:x1]
                    if visit in split_a:
                        s_a[target] += np.where(valid, sample_i * local_w, 0.0)
                        w_a[target] += local_w
                        split_visit_valid["A"].append(int(valid.sum()))
                    else:
                        s_b[target] += np.where(valid, sample_i * local_w, 0.0)
                        w_b[target] += local_w
                        split_visit_valid["B"].append(int(valid.sum()))
        mean_a = np.divide(s_a, w_a, out=np.full_like(s_a, np.nan), where=w_a > 0)
        mean_b = np.divide(s_b, w_b, out=np.full_like(s_b, np.nan), where=w_b > 0)
        var_sum = np.divide(1.0, w_a, out=np.full_like(w_a, np.nan), where=w_a > 0)
        var_sum += np.divide(1.0, w_b, out=np.full_like(w_b, np.nan), where=w_b > 0)
        common = c_valid & (w_a > 0) & (w_b > 0) & np.isfinite(var_sum) & (var_sum > 0)
        z = np.divide(mean_a - mean_b, np.sqrt(var_sum), out=np.full_like(mean_a, np.nan), where=common)
        kpc_per_pixel = float(FlatLambdaCDM(H0=70, Om0=0.3, Tcmb0=0*u.K)
                               .angular_diameter_distance(float(pair["redshift"])).to_value(u.kpc)
                               * 0.168 / (180*3600/math.pi))
        yy, xx = np.indices(image.shape)
        radius = np.hypot((xx + 0.5 - cx) * kpc_per_pixel, (yy + 0.5 - cy) * kpc_per_pixel)
        annulus = (radius >= 200) & (radius < 300) & common
        kernel_sum = np.zeros((17, 17), dtype=np.float64)
        kernel_n = np.zeros((17, 17), dtype=np.int64)
        half = 32
        ny, nx = image.shape
        patch_count = 0
        for y0p in range(0, ny - half, half):
            for x0p in range(0, nx - half, half):
                sl = np.s_[y0p:y0p+64, x0p:x0p+64]
                valid = annulus[sl] & np.isfinite(z[sl])
                if int(valid.sum()) < 500:
                    continue
                py, px = np.mgrid[:64, :64]
                design = np.column_stack((np.ones(valid.sum()), px[valid], py[valid]))
                vals = z[sl][valid]
                keep = np.ones(vals.shape, dtype=bool)
                for _ in range(4):
                    if keep.sum() < 300:
                        break
                    beta = np.linalg.lstsq(design[keep], vals[keep], rcond=None)[0]
                    resid = vals - design @ beta
                    med = np.median(resid[keep])
                    mad = 1.4826 * np.median(np.abs(resid[keep] - med))
                    if not np.isfinite(mad) or mad <= 0:
                        break
                    new_keep = np.abs(resid - med) <= 4 * mad
                    if np.array_equal(new_keep, keep):
                        break
                    keep = new_keep
                resid_grid = np.full((64, 64), np.nan)
                valid_grid = np.zeros((64, 64), dtype=bool)
                all_y, all_x = np.where(valid)
                resid_grid[all_y[keep], all_x[keep]] = resid[keep]
                valid_grid[all_y[keep], all_x[keep]] = True
                if valid_grid.sum() < 300:
                    continue
                patch_count += 1
                for dy in range(-8, 9):
                    for dx in range(-8, 9):
                        y_src = slice(max(0, -dy), min(64, 64-dy))
                        x_src = slice(max(0, -dx), min(64, 64-dx))
                        y_dst = slice(max(0, dy), min(64, 64+dy))
                        x_dst = slice(max(0, dx), min(64, 64+dx))
                        ok = valid_grid[y_src, x_src] & valid_grid[y_dst, x_dst]
                        if ok.any():
                            iy, ix = dy + 8, dx + 8
                            kernel_sum[iy, ix] += np.sum(resid_grid[y_src, x_src][ok] * resid_grid[y_dst, x_dst][ok])
                            kernel_n[iy, ix] += int(ok.sum())
        kernel = np.divide(kernel_sum, kernel_n, out=np.full_like(kernel_sum, np.nan), where=kernel_n > 0)
        center_cov = float(kernel[8, 8]) if np.isfinite(kernel[8, 8]) else None
        finite_kernel = np.nan_to_num(kernel, nan=0.0)
        sym = 0.5 * (finite_kernel + finite_kernel[::-1, ::-1])
        lag_sums = {str(rad): float(sym[8-rad:9+rad, 8-rad:9+rad].sum())
                    for rad in (0, 1, 2, 4, 8)}
        return {
            "request_id": request["request_id"], "role": role,
            "planned_visits": len(visit_ids), "split_a_visits": sorted(split_a),
            "split_b_visits": sorted(split_b), "valid_pixels_200_300_kpc": int(annulus.sum()),
            "finite_normalized_difference_pixels_200_300_kpc": int(np.isfinite(z[annulus]).sum()),
            "coverage_in_geometric_annulus": float(annulus.sum() / max(((radius >= 200) & (radius < 300)).sum(), 1)),
            "split_difference_rms_normalized_before_plane": (float(np.nanstd(z[annulus]))
                                                               if annulus.any() and np.isfinite(z[annulus]).any()
                                                               else None),
            "kernel_patch_count": patch_count, "kernel_lag_pair_counts_center_and_8px": {
                "zero": int(kernel_n[8, 8]), "8px_x": int(kernel_n[8, 16])},
            "normalized_covariance_zero_lag": center_cov,
            "symmetric_cumulative_lag_sum": lag_sums,
            "kernel": finite_list(kernel), "kernel_pair_counts": kernel_n.tolist(),
            "wcs_affine_residual_max_pixel": max(wcs_residuals) if wcs_residuals else None,
            "valid_pixels_per_visit_by_split": split_visit_valid,
            "coadd_fits_sha256": sha256(coadd_path),
            "warp_archive_sha256": sha256(archive_path),
        }


def run(output_dir: Path) -> dict:
    output_dir = Path(output_dir).resolve()
    if output_dir.exists():
        raise FileExistsError(output_dir)
    plan_path = PLAN_DIR / "plan.json"
    summary_path = DOWNLOAD_DIR / "summary.json"
    checkpoint_path = DOWNLOAD_DIR / "checkpoint.json"
    for path in (plan_path, summary_path, checkpoint_path, TARGET_MANIFEST, CONTROL_MANIFEST, PAIR_TABLE):
        if not path.is_file():
            raise FileNotFoundError(path)
    hashes = {str(p): sha256(p) for p in (plan_path, summary_path, checkpoint_path,
                                          TARGET_MANIFEST, CONTROL_MANIFEST, PAIR_TABLE)}
    if hashes[str(plan_path)] != PLAN_SHA256 or hashes[str(summary_path)] != DOWNLOAD_SHA256:
        raise ValueError("approved plan or download summary hash mismatch")
    plan = json.loads(plan_path.read_text())
    download = json.loads(summary_path.read_text())
    checkpoint = json.loads(checkpoint_path.read_text())
    if (plan["status"] != "prepared_not_submitted" or plan["query_submitted"] is not False
            or download["completed_archive_count"] != 184 or download["http_404_request_count"] != 4
            or download["science_measurement"] is not False):
        raise ValueError("download gate, count, or science flag mismatch")
    pairs = list(csv.DictReader(PAIR_TABLE.open(encoding="utf-8"), delimiter="\t"))
    pair = next(row for row in pairs if row["cluster_id"] == "101" and row["filter"] == "i2")
    manifests = {
        "target": json.loads(TARGET_MANIFEST.read_text()),
        "control": json.loads(CONTROL_MANIFEST.read_text()),
    }
    outputs = []
    for role, center_id in (("target", "101"), ("control", pair["control_id"])):
        manifest_rows = manifests[role]["clusters"]
        coadd_row = next(r for r in manifest_rows if str(r["id"]) == center_id)
        req_id = f"i2:{role}:{center_id}"
        index, request = next((i, r) for i, r in enumerate(plan["requests"]) if r["request_id"] == req_id)
        download_record = checkpoint["requests"][f"request_{index:03d}"]
        if download_record.get("status") == "http_404_permanent":
            raise ValueError(f"approved pilot center is unavailable: {req_id}")
        outputs.append(center_covariance(request, index, pair, role, coadd_row,
                                         download_record, output_dir))
    result = {
        "schema": "astrafek-icl-visit-split-covariance-pilot-v1",
        "status": "complete_exploratory",
        "science_measurement": False,
        "hsc_network_used": False,
        "gpu_used": False,
        "cluster_id": "101", "filter": "i2", "matched_pairs_used": 1,
        "method": "same expected visit IDs alternated after sorting; inverse-variance half-stacks; normalized A-B difference; 64x64 local affine-plane residual covariance in 200-300 proper-kpc annulus",
        "limitations": [
            "The nearest-neighbor affine warp-to-coadd mapping is not the HSC coadd resampler.",
            "Visit-split stacks are a single-pair diagnostic and do not support population inference.",
            "Time-variable sources, seeing-dependent source masking, sky residuals, and visit-dependent PSF remain in the split difference.",
            "The normalized kernel is conditional on the delivered warp variances and the selected common mask; no correction is applied to the science stack.",
        ],
        "input_sha256": hashes,
        "centers": outputs,
    }
    output_dir.mkdir(parents=True, exist_ok=False)
    out = output_dir / "summary.json"
    out.write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    result["summary_sha256"] = sha256(out)
    print(json.dumps({"status": result["status"], "science_measurement": False,
                      "output": str(out), "centers": [
                          {k: v for k, v in c.items() if k not in ("kernel", "kernel_pair_counts")}
                          for c in outputs]}, indent=2))
    return result


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    run(parser.parse_args().output_dir)
