#!/usr/bin/env python3
"""Propagate warp-level heteroscedastic white noise through the split estimator."""
from __future__ import annotations

import argparse
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

import analyze_icl_visit_split_covariance_pilot as pilot

ROOT = pilot.ROOT
PLAN_DIR = pilot.PLAN_DIR
DOWNLOAD_DIR = pilot.DOWNLOAD_DIR
TARGET_MANIFEST = pilot.TARGET_MANIFEST
CONTROL_MANIFEST = pilot.CONTROL_MANIFEST
PAIR_TABLE = pilot.PAIR_TABLE
SPLIT_PILOT = ROOT / "icl_visit_split_covariance_pilot_2026-11-06_03/summary.json"
SPLIT_PILOT_SHA256 = "09e557de604d4161cd0281065a59f32669744ab480b180029d192f79ba192359"
DEFAULT_OUTPUT = ROOT / "icl_visit_split_covariance_null_2026-11-06_01"
N_REALIZATIONS = 100
MASK_NAMES = pilot.EXCLUDE_NAMES


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def prepare_center(pair: dict, role: str, plan: dict, checkpoint: dict,
                   manifest: dict, plan_index: dict) -> dict:
    center_id = str(pair["cluster_id"] if role == "target" else pair["control_id"])
    request_id = f"{pair['filter']}:{role}:{center_id}"
    index = plan_index[request_id]
    request = plan["requests"][index]
    record = checkpoint["requests"][f"request_{index:03d}"]
    archive_path = DOWNLOAD_DIR / f"request_{index:03d}.tar"
    if not archive_path.is_file() or sha256(archive_path) != record.get("archive_sha256"):
        raise ValueError(f"archive missing or hash mismatch: {request_id}")
    coadd_row = next(row for row in manifest["clusters"] if str(row["id"]) == center_id)
    image_dir = ROOT / ("download_01/images" if role == "target" else "control_download_01/images")
    coadd_path = image_dir / coadd_row["image"]
    with fits.open(coadd_path, memmap=True) as hdu:
        coadd_image = np.asarray(hdu[1].data)
        coadd_mask = np.asarray(hdu[2].data)
        coadd_variance = np.asarray(hdu[3].data)
        coadd_wcs = WCS(hdu[1].header).celestial
        ra = float(pair["target_ra_deg"] if role == "target" else pair["control_ra_deg"])
        dec = float(pair["target_dec_deg"] if role == "target" else pair["control_dec_deg"])
        center_x, center_y = (float(np.asarray(v)) for v in coadd_wcs.world_to_pixel_values(ra, dec))
        coadd_bits = pilot.mask_bits(hdu[2].header)
        no_data_bit = 1 << coadd_bits["MP_NO_DATA"]
        excluded = sum(1 << coadd_bits[name] for name in MASK_NAMES if name in coadd_bits)
        coadd_valid = (np.isfinite(coadd_image) & np.isfinite(coadd_variance)
                       & (coadd_variance > 0) & ((coadd_mask & no_data_bit) == 0)
                       & ((coadd_mask & excluded) == 0))
        shape = coadd_image.shape
        zred = float(pair["redshift"])
        kpc_per_pixel = float(FlatLambdaCDM(H0=70, Om0=0.3, Tcmb0=0*u.K)
                              .angular_diameter_distance(zred).to_value(u.kpc)
                              * 0.168 / (180*3600/math.pi))
        yy, xx = np.indices(shape)
        radius = np.hypot((xx + 0.5 - center_x) * kpc_per_pixel,
                          (yy + 0.5 - center_y) * kpc_per_pixel)
        annulus = (radius >= 200) & (radius < 300)

        wsum = {"A": np.zeros(shape, dtype=np.float64), "B": np.zeros(shape, dtype=np.float64)}
        visits = []
        visit_ids = sorted(request["expected_visits"])
        split_a = set(visit_ids[::2])
        with tarfile.open(archive_path, mode="r:*") as archive:
            by_visit = {}
            for member in archive.getmembers():
                if member.isfile() and member.name.endswith(".fits"):
                    visit = int(Path(member.name).stem.removeprefix("warp-"))
                    if visit in visit_ids:
                        by_visit[visit] = member
            if set(by_visit) != set(visit_ids):
                raise ValueError(f"warp visits differ from the plan: {request_id}")
            for visit in visit_ids:
                member = by_visit[visit]
                stream = archive.extractfile(member)
                if stream is None:
                    raise ValueError(f"cannot read {member.name}")
                with fits.open(stream, memmap=False) as whdu:
                    w_image = np.asarray(whdu[1].data)
                    w_mask = np.asarray(whdu[2].data)
                    w_variance = np.asarray(whdu[3].data)
                    wcs = WCS(whdu[1].header).celestial
                    x0, y0, x1, y1, sx, sy, inside, residual = pilot.build_visit_map(
                        wcs, coadd_wcs, w_image.shape, shape)
                    if residual > 0.05:
                        raise ValueError("affine mapping residual exceeds pilot threshold")
                    bits = pilot.mask_bits(whdu[2].header)
                    if bits.get("MP_NO_DATA") != coadd_bits.get("MP_NO_DATA"):
                        raise ValueError("MP_NO_DATA bit differs between warp and coadd")
                    warp_excluded = sum(1 << bits[name] for name in MASK_NAMES if name in bits)
                    image = w_image[sy, sx]
                    mask = w_mask[sy, sx]
                    variance = w_variance[sy, sx]
                    valid_local = (inside & np.isfinite(image) & np.isfinite(variance)
                                   & (variance > 0)
                                   & ((mask & (1 << bits["MP_NO_DATA"])) == 0)
                                   & ((mask & warp_excluded) == 0))
                    source_index = np.full(shape, -1, dtype=np.int32)
                    weight = np.zeros(shape, dtype=np.float32)
                    target_slice = np.s_[y0:y1, x0:x1]
                    local_index = sy * w_image.shape[1] + sx
                    source_index[target_slice] = np.where(valid_local, local_index, -1).astype(np.int32)
                    weight[target_slice] = np.where(valid_local, 1.0 / variance, 0.0).astype(np.float32)
                    split = "A" if visit in split_a else "B"
                    wsum[split] += weight
                    visits.append({"visit": visit, "split": split,
                                   "source_shape": tuple(w_image.shape),
                                   "source_index": source_index,
                                   "source_variance": np.asarray(w_variance, dtype=np.float32),
                                   "weight": weight,
                                   "valid_source_fraction": float(valid_local.mean())})
    var_sum = np.divide(1.0, wsum["A"], out=np.full(shape, np.nan), where=wsum["A"] > 0)
    var_sum += np.divide(1.0, wsum["B"], out=np.full(shape, np.nan), where=wsum["B"] > 0)
    common = coadd_valid & annulus & (wsum["A"] > 0) & (wsum["B"] > 0)
    common &= np.isfinite(var_sum) & (var_sum > 0)
    if int(common.sum()) == 0:
        raise ValueError(f"no common annulus pixels: {request_id}")
    return {"request_id": request_id, "role": role, "shape": shape,
            "visits": visits, "wsum": wsum, "normalizer": np.sqrt(var_sum),
            "common": common, "annulus": annulus, "coadd_path": str(coadd_path),
            "archive_path": str(archive_path), "visit_ids": visit_ids,
            "split_a_visits": sorted(split_a), "n_common": int(common.sum()),
            "coverage": float(common.sum() / max(annulus.sum(), 1))}


def estimate_kernel(z: np.ndarray, valid_annulus: np.ndarray) -> tuple[float, float | None, int]:
    h, w = z.shape
    sums = np.zeros((17, 17), dtype=np.float64)
    counts = np.zeros((17, 17), dtype=np.int64)
    n_patches = 0
    for y0 in range(0, h - 63, 32):
        for x0 in range(0, w - 63, 32):
            sl = np.s_[y0:y0+64, x0:x0+64]
            valid = valid_annulus[sl] & np.isfinite(z[sl])
            if int(valid.sum()) < 500:
                continue
            py, px = np.mgrid[:64, :64]
            design = np.column_stack((np.ones(valid.sum()), px[valid], py[valid]))
            values = z[sl][valid]
            keep = np.ones(values.shape, dtype=bool)
            beta = np.zeros(3)
            residual = values.copy()
            for _ in range(4):
                if int(keep.sum()) < 300:
                    break
                beta = np.linalg.lstsq(design[keep], values[keep], rcond=None)[0]
                residual = values - design @ beta
                median = np.median(residual[keep])
                mad = 1.4826 * np.median(np.abs(residual[keep] - median))
                if not np.isfinite(mad) or mad <= 0:
                    break
                new_keep = np.abs(residual - median) <= 4 * mad
                if np.array_equal(new_keep, keep):
                    break
                keep = new_keep
            if int(keep.sum()) < 300:
                continue
            iy, ix = np.where(valid)
            residual_grid = np.full((64, 64), np.nan)
            keep_grid = np.zeros((64, 64), dtype=bool)
            residual_grid[iy[keep], ix[keep]] = residual[keep]
            keep_grid[iy[keep], ix[keep]] = True
            n_patches += 1
            for dy in range(-8, 9):
                for dx in range(-8, 9):
                    ys = slice(max(0, -dy), min(64, 64-dy))
                    xs = slice(max(0, -dx), min(64, 64-dx))
                    yd = slice(max(0, dy), min(64, 64+dy))
                    xd = slice(max(0, dx), min(64, 64+dx))
                    ok = keep_grid[ys, xs] & keep_grid[yd, xd]
                    if ok.any():
                        i, j = dy + 8, dx + 8
                        sums[i, j] += np.sum(residual_grid[ys, xs][ok] * residual_grid[yd, xd][ok])
                        counts[i, j] += int(ok.sum())
    kernel = np.divide(sums, counts, out=np.full_like(sums, np.nan), where=counts > 0)
    zero = float(kernel[8, 8]) if np.isfinite(kernel[8, 8]) else float("nan")
    if np.all(counts > 0) and np.all(np.isfinite(kernel)):
        symmetric = 0.5 * (kernel + kernel[::-1, ::-1])
        lag_sum = float(symmetric.sum())
    else:
        lag_sum = None
    return zero, lag_sum, n_patches


def simulate_center(prepared: dict, n_real: int, seed: int) -> list[dict]:
    rng = np.random.default_rng(seed)
    shape = prepared["shape"]
    out = []
    for realization in range(n_real):
        accum = {"A": np.zeros(shape, dtype=np.float64), "B": np.zeros(shape, dtype=np.float64)}
        for visit in prepared["visits"]:
            source_index = visit["source_index"]
            valid = source_index >= 0
            source_size = int(np.prod(visit["source_shape"]))
            noise = rng.standard_normal(source_size, dtype=np.float32)
            noise *= np.sqrt(np.clip(np.nan_to_num(visit["source_variance"], nan=0.0), 0, None).ravel())
            sampled = np.zeros(shape, dtype=np.float32)
            sampled[valid] = noise[source_index[valid]]
            accum[visit["split"]] += sampled * visit["weight"]
        mean_a = np.divide(accum["A"], prepared["wsum"]["A"], out=np.zeros(shape),
                           where=prepared["wsum"]["A"] > 0)
        mean_b = np.divide(accum["B"], prepared["wsum"]["B"], out=np.zeros(shape),
                           where=prepared["wsum"]["B"] > 0)
        z = np.divide(mean_a - mean_b, prepared["normalizer"],
                      out=np.full(shape, np.nan), where=prepared["common"])
        zero, lag8, n_patches = estimate_kernel(z, prepared["common"])
        raw_rms = float(np.sqrt(np.mean(np.square(z[prepared["common"]]))))
        out.append({"realization": realization, "raw_normalized_rms": raw_rms,
                    "zero_lag_covariance": zero if np.isfinite(zero) else None,
                    "cumulative_covariance_sum_8px": lag8,
                    "kernel_patch_count": n_patches})
        if (realization + 1) % 10 == 0:
            print(f"{prepared['request_id']}: {realization+1}/{n_real} realizations", flush=True)
    return out


def summarize(values: list[dict], key: str) -> dict:
    array = np.asarray([row[key] for row in values if row.get(key) is not None], dtype=np.float64)
    if not array.size:
        return {"n": 0, "mean": None, "sd": None, "q16_median_q84": None}
    return {"n": int(array.size), "mean": float(np.mean(array)), "sd": float(np.std(array, ddof=1)),
            "q16_median_q84": [float(x) for x in np.quantile(array, [0.16, 0.5, 0.84])]}


def run(output: Path, n_real: int) -> dict:
    output = Path(output).resolve()
    if output.exists():
        raise FileExistsError(output)
    files = [PLAN_DIR / "plan.json", DOWNLOAD_DIR / "summary.json", DOWNLOAD_DIR / "checkpoint.json",
             TARGET_MANIFEST, CONTROL_MANIFEST, PAIR_TABLE, SPLIT_PILOT]
    for path in files:
        if not path.is_file():
            raise FileNotFoundError(path)
    hashes = {str(path): sha256(path) for path in files}
    if (hashes[str(PLAN_DIR / "plan.json")] != pilot.PLAN_SHA256
            or hashes[str(DOWNLOAD_DIR / "summary.json")] != pilot.DOWNLOAD_SHA256
            or hashes[str(SPLIT_PILOT)] != SPLIT_PILOT_SHA256):
        raise ValueError("parent plan/download/pilot hashes mismatch")
    plan = json.loads((PLAN_DIR / "plan.json").read_text(encoding="utf-8"))
    download = json.loads((DOWNLOAD_DIR / "summary.json").read_text(encoding="utf-8"))
    checkpoint = json.loads((DOWNLOAD_DIR / "checkpoint.json").read_text(encoding="utf-8"))
    if download.get("completed_archive_count") != 184 or download.get("science_measurement") is not False:
        raise ValueError("approved acquisition checkpoint changed")
    pair = next(r for r in read_tsv(PAIR_TABLE) if r["cluster_id"] == "101" and r["filter"] == "i2")
    manifests = {"target": json.loads(TARGET_MANIFEST.read_text(encoding="utf-8")),
                 "control": json.loads(CONTROL_MANIFEST.read_text(encoding="utf-8"))}
    plan_index = {row["request_id"]: i for i, row in enumerate(plan["requests"])}
    simulations = {}
    observed = {}
    for role in ("target", "control"):
        request_id = f"i2:{role}:{pair['cluster_id'] if role == 'target' else pair['control_id']}"
        index = plan_index[request_id]
        record = checkpoint["requests"][f"request_{index:03d}"]
        if record.get("status") == "http_404_permanent":
            raise ValueError(f"approved input archive missing: {request_id}")
        prepared = prepare_center(pair, role, plan, checkpoint, manifests[role], plan_index)
        sims = simulate_center(prepared, n_real, 9100611 + (0 if role == "target" else 1))
        simulations[role] = sims
        observed[role] = {"request_id": prepared["request_id"],
                          "clean_annulus_pixels": prepared["n_common"],
                          "clean_annulus_coverage": prepared["coverage"],
                          "planned_visit_ids": prepared["visit_ids"],
                          "split_a_visit_ids": prepared["split_a_visits"],
                          "split_b_visit_ids": sorted(set(prepared["visit_ids"]) - set(prepared["split_a_visits"]))}
    observed_pilot = json.loads(SPLIT_PILOT.read_text(encoding="utf-8"))
    observed_by_role = {x["role"]: x for x in observed_pilot["centers"]}
    summaries = {}
    for role in ("target", "control"):
        obs = observed_by_role[role]
        summaries[role] = {
            "observed": {
                "raw_normalized_rms": obs["split_difference_rms_normalized_before_plane"],
                "zero_lag_covariance": obs["normalized_covariance_zero_lag"],
                "cumulative_covariance_sum_8px": obs["symmetric_cumulative_lag_sum"]["8"],
            },
            "noise_only_simulation": {
                key: summarize(simulations[role], key)
                for key in ("raw_normalized_rms", "zero_lag_covariance", "cumulative_covariance_sum_8px")
            },
        }
    result = {
        "schema": "astrafek-icl-visit-split-covariance-noise-null-v1",
        "status": "complete_exploratory", "science_measurement": False,
        "hsc_network_used": False, "gpu_used": False,
        "n_realizations_per_center": n_real,
        "purpose": "Propagate independent heteroscedastic white noise from each warp's supplied variance plane through the exact saved visit masks, visit split, nearest-neighbor mapping, inverse-variance half-stacks, affine-plane clipping, and lag estimator.",
        "limitations": [
            "The mapping is the project's nearest-neighbor approximation and not the HSC production resampling operator.",
            "The null excludes residual sky, astrophysical scene variability, correlated resampling noise present in the input warp, and variance-plane calibration errors.",
            "Only one matched i2 target/control pair is simulated; this calibrates the pilot estimator, not the full sample or ICL signal.",
            "Monte Carlo quantiles are conditional on the observed masks, visit split, and supplied variance maps.",
        ],
        "input_sha256": hashes,
        "centers": {role: {**observed[role], **summaries[role]} for role in ("target", "control")},
        "realizations": simulations,
    }
    output.mkdir(parents=True, exist_ok=False)
    path = output / "summary.json"
    path.write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({"status": result["status"], "science_measurement": False,
                      "summary": str(path), "centers": summaries}, indent=2), flush=True)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--realizations", type=int, default=N_REALIZATIONS)
    args = parser.parse_args()
    run(args.output_dir, args.realizations)
