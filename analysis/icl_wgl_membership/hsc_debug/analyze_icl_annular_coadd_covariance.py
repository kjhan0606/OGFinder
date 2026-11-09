#!/usr/bin/env python3
"""Compare masked coadd covariance in the ICL annulus and outer sky annulus."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import zlib
from collections import defaultdict
from pathlib import Path

import numpy as np
from astropy import units as u
from astropy.cosmology import FlatLambdaCDM
from astropy.io import fits
from astropy.wcs import WCS

import audit_icl_full_sample_coadd_covariance as full
import estimate_hsc_coadd_noise_correlation as noise

ROOT = full.DEFAULT_SAMPLE
PAIR_TABLE = ROOT / "icl_satellite_supported_94_association_2026-11-05_01/supported_pairs.tsv"
DEFAULT_OUTPUT = ROOT / "icl_annular_coadd_covariance_2026-11-06_01"
PATCH_SIZE = 64
MAX_LAG = 8
MAX_PATCHES = 32
MIN_VALID_FRACTION = 0.80
MIN_PATCHES = 4
REPLICATES = 5000
SEED = 20261112
ANNULI = {"icl_200_300": (200.0, 300.0), "outer_350_500": (350.0, 500.0)}


def read_tsv(path: Path) -> list[dict]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def annular_origins(valid: np.ndarray, center_xy: tuple[float, float], kpc_per_pixel: float,
                    inner: float, outer: float) -> tuple[list[tuple[int, int]], int, int]:
    height, width = valid.shape
    cx, cy = center_xy
    candidates, accepted = [], []
    for y0 in range(0, height - PATCH_SIZE + 1, PATCH_SIZE):
        for x0 in range(0, width - PATCH_SIZE + 1, PATCH_SIZE):
            nearest_x = min(max(cx, x0), x0 + PATCH_SIZE - 1)
            nearest_y = min(max(cy, y0), y0 + PATCH_SIZE - 1)
            r_min = math.hypot((nearest_x - cx) * kpc_per_pixel,
                               (nearest_y - cy) * kpc_per_pixel)
            radii = [math.hypot((x - cx) * kpc_per_pixel, (y - cy) * kpc_per_pixel)
                     for x, y in ((x0, y0), (x0 + PATCH_SIZE - 1, y0),
                                  (x0, y0 + PATCH_SIZE - 1),
                                  (x0 + PATCH_SIZE - 1, y0 + PATCH_SIZE - 1))]
            if r_min >= inner and max(radii) < outer:
                candidates.append((y0, x0))
                if float(valid[y0:y0+PATCH_SIZE, x0:x0+PATCH_SIZE].mean()) >= MIN_VALID_FRACTION:
                    accepted.append((y0, x0))
    if len(accepted) > MAX_PATCHES:
        seed = zlib.crc32(f"annular-v1:{inner}:{outer}:{valid.shape}".encode())
        rng = np.random.default_rng(seed)
        keep = rng.choice(len(accepted), size=MAX_PATCHES, replace=False)
        accepted = [accepted[int(i)] for i in sorted(keep)]
    return accepted, len(candidates), sum(
        float(valid[y:y+PATCH_SIZE, x:x+PATCH_SIZE].mean()) >= MIN_VALID_FRACTION
        for y, x in candidates)


def prepare_center(record: dict, pair: dict) -> dict:
    image_path = Path(record["path"])
    if sha256(image_path) != record["source_image_sha256"]:
        raise ValueError(f"coadd hash mismatch: {record['cluster_id']}/{record['role']}")
    with fits.open(image_path, memmap=True, checksum=True) as hdus:
        image = np.asarray(hdus[1].data, dtype=np.float32)
        mask = np.asarray(hdus[2].data, dtype=np.int64)
        variance = np.asarray(hdus[3].data, dtype=np.float32)
        wcs = WCS(hdus[1].header).celestial
    if record["role"] == "target":
        mask_path = Path(record["target_mask_path"])
        if sha256(mask_path) != record["target_mask_sha256"]:
            raise ValueError(f"target K=2.5 mask hash mismatch: {record['cluster_id']}")
        with fits.open(mask_path, memmap=True, checksum=True) as hdus:
            excluded = np.asarray(hdus["MASK_K2P5"].data, dtype=bool)
    else:
        excluded = noise._matched_control_source_mask(
            image, variance, mask, float(record["pixel_scale_arcsec"]),
        )
    if image.shape != variance.shape or image.shape != excluded.shape:
        raise ValueError(f"image/variance/mask shape mismatch: {record['cluster_id']}/{record['role']}")
    valid = np.isfinite(image) & np.isfinite(variance) & (variance > 0) & ~excluded
    ra = float(pair["target_ra_deg"] if record["role"] == "target" else pair["control_ra_deg"])
    dec = float(pair["target_dec_deg"] if record["role"] == "target" else pair["control_dec_deg"])
    cx, cy = map(float, wcs.world_to_pixel_values(ra, dec))
    scale = float(WCS(fits.getheader(image_path, 1)).celestial.proj_plane_pixel_scales()[0].to_value(u.arcsec))
    kpc_per_pixel = float(FlatLambdaCDM(H0=70, Om0=0.3, Tcmb0=0*u.K)
                          .kpc_proper_per_arcmin(float(pair["redshift"])).value * scale / 60.0)
    return {"record": record, "pair": pair, "image": image, "variance": variance,
            "valid": valid, "center_xy": (cx, cy), "kpc_per_pixel": kpc_per_pixel}


def center_annulus(prepared: dict, annulus_name: str,
                   inner: float, outer: float) -> tuple[dict, list[dict]]:
    record, pair = prepared["record"], prepared["pair"]
    image, variance, valid = prepared["image"], prepared["variance"], prepared["valid"]
    cx, cy = prepared["center_xy"]
    kpc_per_pixel = prepared["kpc_per_pixel"]
    origins, n_geometric_patches, n_mask_eligible_patches = annular_origins(
        valid, (cx, cy), kpc_per_pixel, inner, outer)
    annulus_area_pixels = math.pi * (outer**2 - inner**2) / kpc_per_pixel**2
    y, x = np.ogrid[:image.shape[0], :image.shape[1]]
    radial = np.hypot((x - cx) * kpc_per_pixel, (y - cy) * kpc_per_pixel)
    annulus = (radial >= inner) & (radial < outer)
    geometry_pixels = int(annulus.sum())
    valid_pixels = int((annulus & valid).sum())
    summary = {
        "cluster_id": record["cluster_id"], "filter": record["filter"],
        "role": record["role"], "annulus": annulus_name,
        "inner_kpc": inner, "outer_kpc": outer,
        "geometric_pixels_in_coadd": geometry_pixels,
        "analytic_annulus_area_pixels": float(annulus_area_pixels),
        "coadd_geometry_coverage": float(geometry_pixels / max(annulus_area_pixels, 1)),
        "valid_annulus_coverage": float(valid_pixels / max(geometry_pixels, 1)),
        "n_geometric_whole_patches": n_geometric_patches,
        "n_mask_eligible_whole_patches": n_mask_eligible_patches,
        "n_patches_used": len(origins),
        "image_sha256": record["source_image_sha256"],
        "mask_sha256": record["target_mask_sha256"] if record["role"] == "target" else None,
        "connected_5deg_block": pair["connected_5deg_block"],
        "status": "ok" if len(origins) >= MIN_PATCHES else "insufficient_patch_support",
    }
    if len(origins) < MIN_PATCHES:
        return summary, []
    patch_kernels = []
    for y0, x0 in origins:
        sl = np.s_[y0:y0+PATCH_SIZE, x0:x0+PATCH_SIZE]
        patch_kernels.append(noise.lag_correlations(
            image[sl], variance[sl], valid[sl], max_lag=MAX_LAG))
    lag_rows = []
    medians = {}
    for dy in range(-MAX_LAG, MAX_LAG + 1):
        for dx in range(-MAX_LAG, MAX_LAG + 1):
            values = np.asarray([entry[(dy, dx)][0] for entry in patch_kernels], dtype=np.float64)
            values = values[np.isfinite(values)]
            medians[(dy, dx)] = float(np.median(values)) if values.size else None
            q16, q84 = (np.quantile(values, [0.16, 0.84]).tolist()
                        if values.size else (None, None))
            lag_rows.append({"cluster_id": record["cluster_id"], "filter": record["filter"],
                             "role": record["role"], "annulus": annulus_name,
                             "dy": dy, "dx": dx, "n_patches": int(values.size),
                             "median_variance_normalized_covariance": medians[(dy, dx)],
                             "q16_variance_normalized_covariance": q16,
                             "q84_variance_normalized_covariance": q84})
    kernel = np.empty((2*MAX_LAG+1, 2*MAX_LAG+1), dtype=np.float64)
    for dy in range(-MAX_LAG, MAX_LAG+1):
        for dx in range(-MAX_LAG, MAX_LAG+1):
            a, b = medians[(dy, dx)], medians[(-dy, -dx)]
            kernel[dy+MAX_LAG, dx+MAX_LAG] = 0.5 * (a+b) if a is not None and b is not None else np.nan
    if np.all(np.isfinite(kernel)):
        summary["zero_lag_variance_ratio"] = float(kernel[MAX_LAG, MAX_LAG])
        summary["square_covariance_sum_8px"] = float(kernel.sum())
    else:
        summary["zero_lag_variance_ratio"] = medians[(0, 0)]
        summary["square_covariance_sum_8px"] = None
    return summary, lag_rows


def block_bootstrap(values: list[tuple[str, float]], seed: int) -> dict:
    if not values:
        return {"n_pairs": 0, "n_blocks": 0, "median": None, "ci95": None}
    by_block = defaultdict(list)
    for block, value in values:
        by_block[block].append(float(value))
    blocks = sorted(by_block)
    observed = float(np.median([v for _, v in values]))
    if len(blocks) < 3:
        return {"n_pairs": len(values), "n_blocks": len(blocks), "median": observed, "ci95": None}
    rng = np.random.default_rng(seed)
    estimates = np.empty(REPLICATES)
    for i in range(REPLICATES):
        selected = rng.choice(blocks, size=len(blocks), replace=True)
        sample = [value for block in selected for value in by_block[block]]
        estimates[i] = np.median(sample)
    return {"n_pairs": len(values), "n_blocks": len(blocks), "median": observed,
            "ci95": [float(x) for x in np.quantile(estimates, [0.025, 0.975])],
            "replicates": REPLICATES, "seed": seed}


def write_tsv(path: Path, rows: list[dict]) -> str:
    if not rows:
        raise ValueError(f"refusing empty table: {path}")
    with path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t",
                                lineterminator="\n", extrasaction="raise")
        writer.writeheader(); writer.writerows(rows)
    return sha256(path)


def run(output: Path) -> dict:
    output = Path(output).resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite {output}")
    if not PAIR_TABLE.is_file():
        raise FileNotFoundError(PAIR_TABLE)
    pairs = read_tsv(PAIR_TABLE)
    if len(pairs) != 94:
        raise ValueError(f"expected 94 supported target/control pairs, got {len(pairs)}")
    records, setup = full.build_records(ROOT)
    supported = {(str(row["cluster_id"]), str(row["filter"])) for row in pairs}
    selected = [row for row in records if (str(row["cluster_id"]), str(row["filter"])) in supported]
    expected = {(ident, band, role) for ident, band in supported
                for role in ("target", "matched_control")}
    if len(selected) != 188 or {(str(r["cluster_id"]), str(r["filter"]), r["role"]) for r in selected} != expected:
        raise ValueError("coadd/mask records do not exactly cover the supported 94 pairs")
    pair_lookup = {(str(row["cluster_id"]), str(row["filter"])): row for row in pairs}
    rows, lag_rows, failures = [], [], []
    for index, record in enumerate(selected, 1):
        pair = pair_lookup[(str(record["cluster_id"]), str(record["filter"]))]
        print(f"[{index}/{len(selected)}] prepare {record['filter']} "
              f"{record['cluster_id']} {record['role']}", flush=True)
        try:
            prepared = prepare_center(record, pair)
        except Exception as exc:
            for annulus_name in ANNULI:
                failures.append({"cluster_id": record["cluster_id"], "filter": record["filter"],
                                 "role": record["role"], "annulus": annulus_name,
                                 "error": f"{type(exc).__name__}: {exc}"})
            continue
        for annulus_name, (inner, outer) in ANNULI.items():
            print(f"  {annulus_name} {record['filter']} {record['cluster_id']} {record['role']}", flush=True)
            try:
                result, lags = center_annulus(prepared, annulus_name, inner, outer)
                rows.append(result); lag_rows.extend(lags)
            except Exception as exc:
                failures.append({"cluster_id": record["cluster_id"], "filter": record["filter"],
                                 "role": record["role"], "annulus": annulus_name,
                                 "error": f"{type(exc).__name__}: {exc}"})
    result_lookup = {(r["cluster_id"], r["filter"], r["role"], r["annulus"]): r for r in rows}
    comparisons = []
    for pair in pairs:
        ident, band, block = str(pair["cluster_id"]), str(pair["filter"]), pair["connected_5deg_block"]
        row = {"cluster_id": ident, "filter": band, "block": block}
        quartet = True
        for annulus_name in ANNULI:
            for role, key in (("target", "target"), ("matched_control", "control")):
                result = result_lookup.get((ident, band, role, annulus_name))
                if result is None or result["status"] != "ok":
                    quartet = False
                    row[f"{key}_{annulus_name}_zero_lag"] = None
                    row[f"{key}_{annulus_name}_sum8"] = None
                    row[f"{key}_{annulus_name}_coverage"] = None
                else:
                    row[f"{key}_{annulus_name}_zero_lag"] = result["zero_lag_variance_ratio"]
                    row[f"{key}_{annulus_name}_sum8"] = result["square_covariance_sum_8px"]
                    row[f"{key}_{annulus_name}_coverage"] = result["valid_annulus_coverage"]
        row["complete_four_center_annulus_support"] = quartet
        comparisons.append(row)
    metrics = {}
    for band_index, band in enumerate(("i", "i2")):
        band_summary = {}
        band_pairs = [row for row in comparisons if row["filter"] == band]
        for metric in ("zero_lag", "sum8", "coverage"):
            band_summary[metric] = {}
            for annulus_name in ANNULI:
                values = []
                for row in band_pairs:
                    target, control = row[f"target_{annulus_name}_{metric}"], row[f"control_{annulus_name}_{metric}"]
                    if target is not None and control is not None:
                        values.append((row["block"], target-control))
                band_summary[metric][annulus_name] = block_bootstrap(
                    values, SEED + 1000*band_index + 100*list(ANNULI).index(annulus_name) + list(("zero_lag", "sum8", "coverage")).index(metric))
            # Difference-in-differences for transfer of outer annulus to the ICL annulus.
            did = []
            for row in band_pairs:
                vals = [row.get(f"{role}_{ring}_{metric}") for role in ("target", "control")
                        for ring in ANNULI]
                if all(v is not None for v in vals):
                    # (target ICL - target outer) - (control ICL - control outer)
                    d = ((row[f"target_icl_200_300_{metric}"] - row[f"target_outer_350_500_{metric}"])
                         - (row[f"control_icl_200_300_{metric}"] - row[f"control_outer_350_500_{metric}"]))
                    did.append((row["block"], d))
            band_summary[metric]["target_control_radial_difference_in_differences"] = block_bootstrap(
                did, SEED + 50000 + 1000*band_index + list(("zero_lag", "sum8", "coverage")).index(metric))
        metrics[band] = band_summary
    output.mkdir(parents=True, exist_ok=False)
    center_path = output / "per_center_annular_covariance.tsv"
    lag_path = output / "annular_lag_covariance.tsv"
    pair_path = output / "paired_annular_covariance_comparison.tsv"
    center_hash = write_tsv(center_path, rows)
    lag_hash = write_tsv(lag_path, lag_rows)
    pair_hash = write_tsv(pair_path, comparisons)
    summary = {
        "schema": "astrafek-icl-annular-coadd-covariance-v1",
        "status": "complete_exploratory" if not failures else "complete_with_data_gaps",
        "science_measurement": False, "hsc_network_used": False,
        "input_pair_count": len(pairs), "analyzed_center_annulus_count": len(rows),
        "successful_center_annulus_count": sum(r["status"] == "ok" for r in rows),
        "failed_inputs": failures,
        "annuli_proper_kpc": ANNULI,
        "method": {
            "patch_size_pixels": PATCH_SIZE, "max_lag_pixels": MAX_LAG,
            "max_patches_per_center_annulus": MAX_PATCHES,
            "minimum_valid_fraction_per_patch": MIN_VALID_FRACTION,
            "minimum_patches_per_center_annulus": MIN_PATCHES,
            "patch_placement": "non-overlapping 64x64 patches fully contained inside the radial annulus; K=2.5 target mask and matched-control SEP/Kron source mask",
            "lag_estimator": "existing Astrafek HSC coadd estimator: local affine-plane removal, 4-sigma fit clipping, variance normalization, 5-sigma residual clipping",
            "bootstrap": "paired target-control differences resampled by saved connected 5-degree sky blocks; 5,000 draws",
        },
        "paired_block_bootstrap": metrics,
        "limitations": [
            "The 350-500 kpc outer annulus may contain real cluster-envelope light; it is not assumed to be a blank sky mean.",
            "This measures covariance of delivered coadds after source masking; it does not identify correlated warp noise versus sky/source residuals or variance calibration.",
            "The kernel is a median of local patches and the ±8 square sum is an estimator diagnostic, not the exact radial-aperture covariance.",
            "No empirical covariance correction is transferred to the ICL or weak-lensing science profiles.",
        ],
        "input_sha256": {str(PAIR_TABLE): sha256(PAIR_TABLE), **setup["provenance"]},
        "outputs_sha256": {str(center_path): center_hash, str(lag_path): lag_hash, str(pair_path): pair_hash},
    }
    summary_path = output / "summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True, allow_nan=False) + "\n",
                            encoding="utf-8")
    print(json.dumps({"status": summary["status"], "science_measurement": False,
                      "input_pairs": len(pairs), "successful_center_annuli": summary["successful_center_annulus_count"],
                      "output": str(summary_path), "metrics": metrics}, indent=2), flush=True)
    return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    run(args.output_dir)
