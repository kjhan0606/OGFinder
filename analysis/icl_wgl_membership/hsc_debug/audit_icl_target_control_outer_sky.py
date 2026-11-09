#!/usr/bin/env python3
"""Compare observed outer-annulus residuals in existing HSC target/control pairs.

This is an offline diagnostic, not a sky subtraction. It uses the 600--900 kpc
annulus and the same K=2.5 source-mask family as the target-control ICL stack.
The target annulus can contain cluster light, ICL, PSF wings, cirrus and missed
sources; the measured intercept is therefore a residual indicator, not a pure
sky estimate. Its association sensitivity is reported without activating a
science measurement.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import platform
from pathlib import Path

import numpy as np
from astropy.cosmology import FlatLambdaCDM
from astropy.io import fits

from download_pilot import validate_fits
from evaluate_prestack_pair_support import pair_block_labels
from validate_control_images import control_masks, image_path as control_image_path, load_controls, nJy_factor
from diagnose_outer_sky import tile_diagnostics


DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
DEFAULT_OUTPUT = DATA / "icl_target_control_outer_sky_2026-11-03_01"
ASSOCIATION_TABLE = DATA / "icl_pixel_extinction_2026-10-19_03/pixel_extinction_satellite_association.tsv"


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


def write_tsv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty table: {path}")
    with path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t", extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def measure_outer_frame(path: Path, frame: dict, mask_path: Path | None,
                        cosmology: FlatLambdaCDM, expected_bytes: int) -> dict:
    if not path.is_file() or path.stat().st_size != int(expected_bytes):
        raise ValueError(f"missing image or delivered file-size mismatch: {path}")
    info = validate_fits(path, frame)
    with fits.open(path, memmap=False, checksum=True) as hdus:
        image = np.asarray(hdus[info["science_hdu"]].data, dtype=np.float32)
        variance = np.asarray(hdus[info["variance_hdu"]].data, dtype=np.float32)
        bits = np.asarray(hdus[info["mask_hdu"]].data, dtype=np.int32)
        if image.shape != variance.shape or image.shape != bits.shape:
            raise ValueError(f"image/variance/mask shape mismatch: {path}")
        if mask_path is None:
            center = np.asarray(info["center_pixel_xy"], dtype=float)
            base, variants, _, _, _ = control_masks(
                image, variance, bits, center, info["pixel_scale_arcsec"],
                info["mask_planes"])
            mask = np.asarray(variants[0], dtype=bool)
        else:
            with fits.open(mask_path, memmap=False, checksum=True) as masks:
                if "MASK_K2P5" not in masks:
                    raise ValueError(f"target mask lacks MASK_K2P5: {mask_path}")
                mask = np.asarray(masks["MASK_K2P5"].data, dtype=bool)
            if mask.shape != image.shape:
                raise ValueError(f"target image/mask shape mismatch: {path}")
            center = np.asarray(info["center_pixel_xy"], dtype=float)

    factor = nJy_factor(info)
    image = image.astype(np.float64) * factor
    variance = variance.astype(np.float64) * (factor * factor)
    kpc_per_pixel = float(cosmology.kpc_proper_per_arcmin(float(frame["z"])).value) \
        * float(info["pixel_scale_arcsec"]) / 60.0
    yy, xx = np.ogrid[:image.shape[0], :image.shape[1]]
    radius = np.hypot(xx - center[0], yy - center[1]) * kpc_per_pixel
    region = (radius >= 600.0) & (radius < 900.0)
    tiles, lag_rows, diagnostic = tile_diagnostics(
        image, variance, mask, region, center, kpc_per_pixel)
    if diagnostic["status"] != "diagnostic_only":
        raise ValueError(f"fewer than ten usable sky tiles: {path}")
    return {
        "n_tiles": diagnostic["n_tiles"],
        "valid_fraction": diagnostic["outer_valid_fraction"],
        "median_tile_mean_njy_arcsec2": diagnostic["median_tile_mean"],
        "plane_constant_njy_arcsec2": diagnostic["plane_constant"],
        "plane_x_njy_arcsec2_per_mpc": diagnostic["plane_x_per_mpc"],
        "plane_y_njy_arcsec2_per_mpc": diagnostic["plane_y_per_mpc"],
        "plane_residual_mad_njy_arcsec2": diagnostic["plane_residual_mad"],
        "plane_residual_mad_over_nominal": diagnostic["plane_residual_mad_over_nominal"],
        "median_nominal_tile_mean_sigma_njy_arcsec2": diagnostic["median_nominal_mean_sigma"],
        "lag1_x_correlation": next((r["correlation"] for r in lag_rows
                                     if r["lag_pixels"] == 1 and r["axis"] == "x"), None),
        "lag1_y_correlation": next((r["correlation"] for r in lag_rows
                                     if r["lag_pixels"] == 1 and r["axis"] == "y"), None),
        "pixel_scale_arcsec": info["pixel_scale_arcsec"],
        "filter": info["coadd_header_filter"],
    }


def partial_pearson(x, y, covariates) -> float:
    arrays = [np.asarray(v, dtype=float) for v in (x, y, *covariates)]
    if len(arrays[0]) < 6 or any(a.ndim != 1 or len(a) != len(arrays[0]) for a in arrays):
        return float("nan")
    if any(not np.isfinite(a).all() for a in arrays):
        return float("nan")
    design = np.column_stack((np.ones(len(arrays[0])), *arrays[2:]))
    if np.linalg.matrix_rank(design) != design.shape[1]:
        return float("nan")
    rx = arrays[0] - design @ np.linalg.lstsq(design, arrays[0], rcond=None)[0]
    ry = arrays[1] - design @ np.linalg.lstsq(design, arrays[1], rcond=None)[0]
    denominator = np.linalg.norm(rx) * np.linalg.norm(ry)
    return float(rx @ ry / denominator) if denominator > 0 else float("nan")


def block_bootstrap_association(rows: list[dict], labels: np.ndarray,
                                replicas: int, seed: int) -> dict:
    blocks = np.unique(labels)
    if len(blocks) < 6:
        raise ValueError("at least six connected target/control blocks are required")
    by_block = {label: np.flatnonzero(labels == label) for label in blocks}
    arrays = {key: np.asarray([float(row[key]) for row in rows]) for key in
              ("icl_delta_pixel_sfd_njy_arcsec2", "satellite_delta_center_sfd_njy_arcsec2",
               "z", "richness", "outer_plane_delta_njy_arcsec2")}
    base_cov = [arrays["z"], arrays["richness"]]
    sky_cov = [*base_cov, arrays["outer_plane_delta_njy_arcsec2"]]
    baseline = partial_pearson(arrays["icl_delta_pixel_sfd_njy_arcsec2"],
                               arrays["satellite_delta_center_sfd_njy_arcsec2"], base_cov)
    adjusted = partial_pearson(arrays["icl_delta_pixel_sfd_njy_arcsec2"],
                               arrays["satellite_delta_center_sfd_njy_arcsec2"], sky_cov)
    rng = np.random.default_rng(seed)
    base_draws, adjusted_draws = [], []
    for _ in range(replicas):
        chosen = rng.choice(blocks, size=len(blocks), replace=True)
        idx = np.concatenate([by_block[b] for b in chosen])
        base_draws.append(partial_pearson(arrays["icl_delta_pixel_sfd_njy_arcsec2"][idx],
                                          arrays["satellite_delta_center_sfd_njy_arcsec2"][idx],
                                          [a[idx] for a in base_cov]))
        adjusted_draws.append(partial_pearson(arrays["icl_delta_pixel_sfd_njy_arcsec2"][idx],
                                              arrays["satellite_delta_center_sfd_njy_arcsec2"][idx],
                                              [a[idx] for a in sky_cov]))
    base_draws, adjusted_draws = np.asarray(base_draws), np.asarray(adjusted_draws)
    valid_base, valid_adjusted = np.isfinite(base_draws), np.isfinite(adjusted_draws)
    paired = valid_base & valid_adjusted
    delta = adjusted_draws[paired] - base_draws[paired]
    q = lambda x: [float(np.percentile(x, 2.5)), float(np.percentile(x, 97.5))] if len(x) else None
    return {
        "n_pairs": len(rows), "n_connected_blocks": int(len(blocks)),
        "n_valid_base_replicas": int(valid_base.sum()),
        "n_valid_sky_adjusted_replicas": int(valid_adjusted.sum()),
        "baseline_partial_r_z_richness": baseline,
        "baseline_block_bootstrap_ci95": q(base_draws[valid_base]),
        "outer_residual_adjusted_partial_r_z_richness_sky": adjusted,
        "outer_residual_adjusted_block_bootstrap_ci95": q(adjusted_draws[valid_adjusted]),
        "paired_adjusted_minus_baseline_r": adjusted - baseline,
        "paired_delta_ci95": q(delta), "bootstrap_replicas": replicas, "seed": seed,
        "interpretation": "diagnostic covariate sensitivity only; outer annulus is not clean sky and no measurement was subtracted",
    }


def run(sample: Path, output: Path, replicas: int, seed: int) -> dict:
    if not platform.node().lower().startswith("syn") and not __import__("os").environ.get("SLURM_JOB_ID"):
        raise RuntimeError("Run the full-sample FITS analysis on syntax Slurm, not a login node")
    if output.exists():
        raise FileExistsError(f"refusing to overwrite {output}")
    sample_doc = read_json(sample / "sample_manifest.json")
    target_masks = sample / "prestack_cluster_masks_2026-10-16_01"
    mask_summary_path = target_masks / "summary.json"
    mask_records_path = target_masks / "records.tsv"
    mask_summary = read_json(mask_summary_path)
    mask_records = read_tsv(mask_records_path)
    if (mask_summary.get("status") != "complete" or mask_summary.get("n_success") != 800
            or mask_summary.get("mask_ready_for_science") is not False):
        raise ValueError("target mask bank is incomplete or has an unexpected science gate")
    target_by_id = {str(row["id"]): row for row in sample_doc["clusters"]}
    mask_by_id = {str(row["id"]): row for row in mask_records}
    if len(target_by_id) != 800 or len(mask_by_id) != 800:
        raise ValueError("expected exactly 800 targets and target masks")
    controls_doc, control_downloads, control_base, control_hashes = load_controls(sample)
    if len(controls_doc["clusters"]) != 798:
        raise ValueError("expected the existing 798 matched controls")
    target_status = read_json(sample / "download_01/download_status.json")
    target_downloads = {str(row["id"]): row for row in target_status["completed"]}
    controls_by_target = {str(row["target_cluster_id"]): row for row in controls_doc["clusters"]}
    control_status_by_id = control_downloads
    cosmo = FlatLambdaCDM(H0=float(sample_doc["h0"]), Om0=float(sample_doc["omega_m"]), Tcmb0=0)

    selected = []
    for ident, target in target_by_id.items():
        control = controls_by_target.get(ident)
        mask_record = mask_by_id[ident]
        if not control or control.get("coadd_header_filter") != mask_record["primary_filter"]:
            continue
        if str(control.get("target_cluster_id")) != ident:
            raise ValueError(f"control association mismatch for target {ident}")
        selected.append((target, mask_record, control))
    if len(selected) != 783:
        raise ValueError(f"expected 783 same-filter pairs, got {len(selected)}")

    output.mkdir(parents=True, exist_ok=False)
    frame_rows, pair_rows = [], []
    for index, (target, mask_record, control) in enumerate(selected, 1):
        ident, control_id = str(target["id"]), str(control["id"])
        try:
            target_path = Path(mask_record["image"])
            target_mask_path = Path(mask_record["mask_file"])
            control_path = control_image_path(control_base, control)
            target_status = target_downloads[ident]
            control_status = control_status_by_id[control_id]
            if (target_path.stat().st_size != int(target_status["bytes"])
                    or not target_mask_path.is_file() or target_mask_path.stat().st_size <= 0
                    or control_path.stat().st_size != int(control_status["bytes"])):
                raise ValueError("current file size differs from the completed download/mask record")
            target_result = measure_outer_frame(target_path, target, target_mask_path,
                                                cosmo, target_status["bytes"])
            control_result = measure_outer_frame(control_path, control, None,
                                                 cosmo, control_status["bytes"])
            if target_result["filter"] != control_result["filter"]:
                raise ValueError("target/control FITS filters differ")
            common = {"target_id": ident, "control_id": control_id,
                      "filter": target_result["filter"], "z": float(target["z"]),
                      "target_ra_deg": float(target["ra_deg"]),
                      "target_dec_deg": float(target["dec_deg"]),
                      "control_ra_deg": float(control["ra_deg"]),
                      "control_dec_deg": float(control["dec_deg"])}
            target_row = {**common, "role": "target", **target_result,
                          "input_sha256_from_download_status": target_status["sha256"],
                          "mask_sha256_from_mask_record": mask_record["mask_sha256"]}
            control_row = {**common, "role": "matched_control", **control_result,
                           "input_sha256_from_download_status": control_status["sha256"]}
            frame_rows.extend((target_row, control_row))
            pair_rows.append({**common,
                "target_minus_control_plane_constant_njy_arcsec2": (
                    target_result["plane_constant_njy_arcsec2"] - control_result["plane_constant_njy_arcsec2"]),
                "target_minus_control_median_tile_mean_njy_arcsec2": (
                    target_result["median_tile_mean_njy_arcsec2"] - control_result["median_tile_mean_njy_arcsec2"]),
                "target_plane_residual_mad_over_nominal": target_result["plane_residual_mad_over_nominal"],
                "control_plane_residual_mad_over_nominal": control_result["plane_residual_mad_over_nominal"],
                "target_valid_fraction": target_result["valid_fraction"],
                "control_valid_fraction": control_result["valid_fraction"],
                "target_n_tiles": target_result["n_tiles"], "control_n_tiles": control_result["n_tiles"]})
        except Exception as exc:
            raise RuntimeError(f"pair {index}/{len(selected)} target={ident} control={control_id}: {exc}") from exc
        if index % 25 == 0 or index == len(selected):
            print(f"outer-sky pair diagnostics {index}/{len(selected)}", flush=True)

    target_coord_rows = [{"ra_deg": r["target_ra_deg"], "dec_deg": r["target_dec_deg"]} for r in pair_rows]
    control_coord_rows = [{"ra_deg": r["control_ra_deg"], "dec_deg": r["control_dec_deg"]} for r in pair_rows]
    _, labels, n_target_cells, n_control_cells, n_blocks, max_block = pair_block_labels(
        target_coord_rows, control_coord_rows, np.ones(len(pair_rows), dtype=bool))
    for row, label in zip(pair_rows, labels):
        row["connected_5deg_block"] = str(label)

    association_rows = read_tsv(ASSOCIATION_TABLE)
    by_key = {(str(r["cluster_id"]), str(r["filter"])): r for r in association_rows}
    association_by_filter = {}
    association_block_meta = {}
    for band in ("i2", "i"):
        merged = []
        for pair in pair_rows:
            if pair["filter"] != band:
                continue
            assoc = by_key.get((pair["target_id"], band))
            if assoc is None:
                continue
            merged.append({
                "target_id": pair["target_id"], "filter": band,
                "icl_delta_pixel_sfd_njy_arcsec2": float(assoc["icl_delta_pixel_sfd_njy_arcsec2"]),
                "satellite_delta_center_sfd_njy_arcsec2": float(assoc["satellite_delta_center_sfd_njy_arcsec2"]),
                "z": float(assoc["z"]), "richness": float(assoc["richness"]),
                "outer_plane_delta_njy_arcsec2": float(pair["target_minus_control_plane_constant_njy_arcsec2"]),
                "target_ra_deg": float(pair["target_ra_deg"]),
                "target_dec_deg": float(pair["target_dec_deg"]),
                "control_ra_deg": float(pair["control_ra_deg"]),
                "control_dec_deg": float(pair["control_dec_deg"]),
            })
        if merged:
            # Define the covariance partition on the actual association sample
            # within this filter. The full 783-pair graph can bridge otherwise
            # separate association regions through unrelated clusters.
            target_coords = [{"ra_deg": r["target_ra_deg"], "dec_deg": r["target_dec_deg"]}
                             for r in merged]
            control_coords = [{"ra_deg": r["control_ra_deg"], "dec_deg": r["control_dec_deg"]}
                              for r in merged]
            _, assoc_labels, n_target_cells_assoc, n_control_cells_assoc, n_blocks_assoc, max_assoc = \
                pair_block_labels(target_coords, control_coords, np.ones(len(merged), dtype=bool))
            for row, label in zip(merged, assoc_labels):
                row["connected_5deg_block"] = str(label)
            association_block_meta[band] = {
                "target_cells": n_target_cells_assoc, "control_cells": n_control_cells_assoc,
                "connected_blocks": n_blocks_assoc, "largest_block_pairs": max_assoc,
                "scope": "recomputed within the fixed ICL-satellite association subset for this filter",
            }
            association_by_filter[band] = block_bootstrap_association(
                merged, np.asarray([r["connected_5deg_block"] for r in merged]), replicas,
                seed + (0 if band == "i2" else 1))
            association_by_filter[band]["spatial_block_assignment"] = association_block_meta[band]

    write_tsv(output / "frame_outer_sky_diagnostics.tsv", frame_rows)
    write_tsv(output / "paired_outer_sky_diagnostics.tsv", pair_rows)
    hashes = {
        str((sample / "sample_manifest.json").resolve()): sha256(sample / "sample_manifest.json"),
        str((sample / "download_01/download_status.json").resolve()): sha256(sample / "download_01/download_status.json"),
        str((sample / "control_download_01/download_manifest.json").resolve()): sha256(sample / "control_download_01/download_manifest.json"),
        str((sample / "control_download_01/download_status.json").resolve()): sha256(sample / "control_download_01/download_status.json"),
        str(mask_summary_path.resolve()): sha256(mask_summary_path),
        str(mask_records_path.resolve()): sha256(mask_records_path),
        str(ASSOCIATION_TABLE.resolve()): sha256(ASSOCIATION_TABLE),
        str(Path(__file__).resolve()): sha256(Path(__file__).resolve()),
    }
    result = {
        "schema": "astrafek-icl-target-control-outer-sky-pair-audit-v1",
        "status": "complete_diagnostic",
        "science_measurement": False,
        "n_target_control_pairs": len(pair_rows),
        "n_target_frames": len(frame_rows) // 2,
        "filter_counts": {band: sum(r["filter"] == band for r in pair_rows) for band in ("i", "i2")},
        "outer_annulus_kpc": [600.0, 900.0], "tile_pixels": 64,
        "source_mask": "target saved MASK_K2P5; matched control regenerated SEP/Kron-2.5; native bad/bright pixels excluded",
        "target_control_spatial_blocks": {"target_cells": n_target_cells, "control_cells": n_control_cells,
                                           "connected_blocks": n_blocks, "largest_block_pairs": max_block,
                                           "rule": "connected pairs sharing either endpoint's equal-area 5-degree sky cell"},
        "association_outer_residual_sensitivity_by_filter": association_by_filter,
        "association_spatial_blocks_by_filter": association_block_meta,
        "input_sha256": hashes,
        "image_content_rehashed": False,
        "image_content_integrity": "current file sizes matched completed download records; recorded image SHA-256 values were inherited from prior verified acquisitions, not recomputed here",
        "hsc_network_used": False,
        "data_requested": False,
        "sky_subtracted": False,
        "limitations": [
            "The 600--900 kpc target annulus can contain real cluster/ICL light, PSF wings, cirrus, and undetected sources; the fitted plane intercept is a residual indicator, not pure sky.",
            "The tile residual MAD divided by nominal mean sigma is descriptive and is not a covariance correction.",
            "The association adjustment is a paired spatial-block covariate sensitivity on the existing exploratory sample; no flux is subtracted and it is not a causal or confirmatory test.",
            "Image contents were not rehashed in this run; acquisition SHA-256 values and file sizes were cross-referenced from prior completed records.",
            "PSF-wing transfer, source completeness, extinction/passband uncertainty, upstream sky transfer, and full pixel covariance remain unresolved.",
        ],
    }
    with (output / "summary.json").open("x", encoding="utf-8") as stream:
        json.dump(result, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    print(json.dumps({"status": result["status"], "pairs": len(pair_rows),
                      "filters": result["filter_counts"], "blocks": n_blocks,
                      "association": association_by_filter,
                      "science_measurement": False, "output": str(output)}, ensure_ascii=False), flush=True)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", type=Path, default=DATA)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--replicas", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=20261103)
    args = parser.parse_args()
    if args.replicas < 1000:
        parser.error("at least 1,000 spatial block bootstrap replicas are required")
    try:
        run(args.sample.resolve(), args.output.resolve(), args.replicas, args.seed)
    except (OSError, ValueError, RuntimeError) as exc:
        parser.error(str(exc))


if __name__ == "__main__":
    main()
