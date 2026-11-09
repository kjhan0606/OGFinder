#!/usr/bin/env python3
"""Recompute the supported ICL--satellite association under warp-mask support.

Uses only the already approved/downloaded 94-pair HSC-I warp cutouts and local
coadd/mask products. The strict pixel support is the coadd's finite variance
intersected with at least three valid visit warps and mask planes shared by all
returned warps. The baseline and strict ICL apertures use the same saved K=2.5
source masks. This is a mask-selection sensitivity, not a corrected science
measurement.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import sys
import tarfile
from pathlib import Path

import numpy as np
from astropy import units as u
from astropy.cosmology import FlatLambdaCDM
from astropy.io import fits
from astropy.wcs import WCS

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import analyze_icl_supported94_warp_noise as warp_noise
from evaluate_prestack_pair_support import pair_block_labels
from audit_icl_target_control_outer_sky import partial_pearson
from prestack_core import geometry, measure

ROOT = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
PLAN_DIR = ROOT / "icl_supported94_warp_noise_plan_2026-11-06_01"
DOWNLOAD_DIR = ROOT / "icl_supported94_warp_noise_download_2026-11-06_01"
PAIR_TABLE = ROOT / "icl_satellite_supported_94_association_2026-11-05_01/supported_pairs.tsv"
PROFILE_TABLE = ROOT / "icl_annulus_sensitivity_2026-10-18_02/selected_pair_profiles.tsv"
TARGET_MASK_RECORDS = ROOT / "prestack_cluster_masks_2026-10-16_01/records.tsv"
CONTROL_EXTRACTION = ROOT / "prestack_extraction_2026-10-15_01"
DEFAULT_OUTPUT = ROOT / "icl_supported94_warp_mask_association_2026-11-06_01"
NO_DATA = "MP_NO_DATA"
EXCLUDE_NAMES = warp_noise.MASK_EXCLUDE_NAMES


def read_tsv(path: Path) -> list[dict]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty table: {path}")
    with path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t",
                                lineterminator="\n", extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def strict_warp_support(request: dict, info: dict, archive_path: Path,
                        checkpoint_record: dict) -> tuple[np.ndarray, list[str], dict]:
    if not archive_path.is_file() or warp_noise.sha256(archive_path) != checkpoint_record.get("archive_sha256"):
        raise ValueError(f"archive checksum mismatch: {archive_path}")
    with fits.open(info["image_path"], memmap=True, checksum=True) as coadd:
        image = np.asarray(coadd[1].data)
        mask = np.asarray(coadd[2].data, dtype=np.int64)
        variance = np.asarray(coadd[3].data)
        wcs = WCS(coadd[1].header).celestial
        cbits = warp_noise.mask_bits(coadd[2].header)
        if NO_DATA not in cbits:
            raise ValueError("coadd mask lacks MP_NO_DATA")
        no_data = 1 << cbits[NO_DATA]
        available = (np.isfinite(image) & np.isfinite(variance) & (variance > 0)
                     & ((mask & no_data) == 0))
        shape = image.shape

    counts = np.zeros(shape, dtype=np.uint8)
    found, per_visit_bits, residual_max = set(), [], 0.0
    expected = set(request["expected_visits"])
    with tarfile.open(archive_path, mode="r:*") as archive:
        members = [m for m in archive.getmembers()
                   if m.isfile() and m.name.endswith(".fits")]
        for member in members:
            visit = int(Path(member.name).stem.removeprefix("warp-"))
            if visit in found:
                raise ValueError(f"duplicate visit {visit} in {request['request_id']}")
            found.add(visit)
            stream = archive.extractfile(member)
            if stream is None:
                raise ValueError(f"cannot extract {member.name}")
            with fits.open(stream, memmap=False) as warp:
                wi = np.asarray(warp[1].data)
                wm = np.asarray(warp[2].data, dtype=np.int64)
                wv = np.asarray(warp[3].data)
                if wi.shape != wm.shape or wi.shape != wv.shape:
                    raise ValueError(f"unaligned warp planes: {member.name}")
                bits = warp_noise.mask_bits(warp[2].header)
                if bits.get(NO_DATA) != cbits.get(NO_DATA):
                    raise ValueError(f"MP_NO_DATA bit mismatch: {member.name}")
                w_no_data = 1 << bits[NO_DATA]
                warp_wcs = WCS(warp[1].header).celestial
                sx_grid, sy_grid = np.meshgrid(
                    np.linspace(0, wi.shape[1]-1, 7), np.linspace(0, wi.shape[0]-1, 7))
                wx, wy = warp_wcs.pixel_to_world_values(sx_grid, sy_grid)
                cx, cy = wcs.world_to_pixel_values(wx, wy)
                design = np.column_stack((np.ones(cx.size), sx_grid.ravel(), sy_grid.ravel()))
                fx = np.linalg.lstsq(design, cx.ravel(), rcond=None)[0]
                fy = np.linalg.lstsq(design, cy.ravel(), rcond=None)[0]
                residual = float(np.max(np.hypot(design @ fx-cx.ravel(), design @ fy-cy.ravel())))
                residual_max = max(residual_max, residual)
                if residual > 0.05:
                    raise ValueError(f"affine WCS residual {residual:.6f} px exceeds 0.05")
                corners = np.array([[1, 0, 0], [1, wi.shape[1]-1, 0],
                                    [1, 0, wi.shape[0]-1], [1, wi.shape[1]-1, wi.shape[0]-1]])
                xcorn, ycorn = corners @ fx, corners @ fy
                x0, y0 = max(0, int(np.floor(xcorn.min()))), max(0, int(np.floor(ycorn.min())))
                x1 = min(shape[1], int(np.ceil(xcorn.max()))+1)
                y1 = min(shape[0], int(np.ceil(ycorn.max()))+1)
                if x1 <= x0 or y1 <= y0:
                    raise ValueError(f"warp misses coadd footprint: {member.name}")
                yy, xx = np.mgrid[y0:y1, x0:x1]
                inverse = np.linalg.inv(np.array([[fx[1], fx[2]], [fy[1], fy[2]]]))
                source = inverse @ np.vstack((xx.ravel()-fx[0], yy.ravel()-fy[0]))
                wxp = source[0].reshape(xx.shape); wyp = source[1].reshape(yy.shape)
                inside = ((wxp >= -0.5) & (wxp < wi.shape[1]-0.5)
                          & (wyp >= -0.5) & (wyp < wi.shape[0]-0.5))
                ix = np.clip(np.rint(wxp).astype(int), 0, wi.shape[1]-1)
                iy = np.clip(np.rint(wyp).astype(int), 0, wi.shape[0]-1)
                sample_mask, sample_image, sample_var = wm[iy, ix], wi[iy, ix], wv[iy, ix]
                valid = (inside & np.isfinite(sample_image) & np.isfinite(sample_var)
                         & (sample_var > 0) & ((sample_mask & w_no_data) == 0))
                shared = set(EXCLUDE_NAMES) & set(cbits) & set(bits)
                excluded = warp_noise.bitmask(shared, bits)
                clean = valid & ((sample_mask & excluded) == 0)
                local = counts[y0:y1, x0:x1]
                local[clean] = np.minimum(local[clean].astype(np.uint16)+1, 255).astype(np.uint8)
                per_visit_bits.append(set(bits))
    if not expected.issubset(found):
        raise ValueError(f"missing metadata-expected visits for {request['request_id']}")
    shared_names = set(EXCLUDE_NAMES) & set(cbits)
    for bits in per_visit_bits:
        shared_names.intersection_update(bits)
    coadd_excluded = warp_noise.bitmask(sorted(shared_names), cbits)
    strict = available & ((mask & coadd_excluded) == 0) & (counts >= 3)
    return strict, sorted(shared_names), {
        "n_expected_visits": len(expected), "n_returned_visits": len(found),
        "n_additional_visits": len(found-expected), "max_affine_wcs_residual_pixel": residual_max,
        "strict_valid_pixel_fraction_of_coadd": float(strict.mean()),
    }


def load_source_mask(role: str, cluster_id: str, control_id: str,
                     image_shape: tuple[int, int], mask_records: dict) -> tuple[np.ndarray, str, str]:
    if role == "target":
        record = mask_records[cluster_id]
        path = Path(record["mask_file"])
        expected = str(record["mask_sha256"])
        if warp_noise.sha256(path) != expected:
            raise ValueError(f"K=2.5 source-mask checksum mismatch: {path}")
        with fits.open(path, memmap=True, checksum=True) as hdus:
            mask = np.asarray(hdus["MASK_K2P5"].data, dtype=bool)
        return mask, str(path), expected
    path = CONTROL_EXTRACTION / control_id / "baseline.npz"
    with np.load(path) as data:
        shape = tuple(int(v) for v in data["image_shape"])
        bits = np.asarray(data["mask_bits"][0], dtype=np.uint8)
        mask = np.unpackbits(bits, axis=-1, count=shape[1], bitorder="little").astype(bool).reshape(shape)
    return mask, str(path), warp_noise.sha256(path)


def reconstruct_profile(image_path: Path, source_mask: np.ndarray, center_ra: float,
                        center_dec: float, redshift: float, sample: dict) -> dict:
    with fits.open(image_path, memmap=True, checksum=True) as hdus:
        image = np.asarray(hdus[1].data, dtype=np.float32)
        variance = np.asarray(hdus[3].data, dtype=np.float32)
        primary = hdus[0].header
        wcs = WCS(hdus[1].header).celestial
        cx, cy = wcs.world_to_pixel_values(center_ra, center_dec)
        scale = 3600.0 * math.sqrt(abs(float(np.linalg.det(wcs.pixel_scale_matrix))))
        factor = 3.631e12 / float(primary["FLUXMAG0"]) / scale**2
    kpp = float(FlatLambdaCDM(H0=float(sample["h0"]), Om0=float(sample["omega_m"]), Tcmb0=0*u.K)
                             .kpc_proper_per_arcmin(redshift).value * scale / 60.0)
    geom = geometry(image.shape, (cx, cy), kpp)
    result = measure(image*factor, variance*factor**2, source_mask, geom, (cx, cy), kpp)
    return {"flux_njy_arcsec2": float(result["raw"][3]),
            "coverage": float(result["coverage"][3]), "n_pixels": int(result["count"][3]),
            "center_pixel_xy": [float(cx), float(cy)], "pixel_scale_arcsec": scale,
            "kpc_per_pixel": kpp}


def block_bootstrap(rows: list[dict], labels: np.ndarray, band: str,
                    replicas: int, seed: int) -> dict:
    selected = [row for row in rows if row["filter"] == band]
    indices = np.asarray([i for i,row in enumerate(rows) if row["filter"] == band], dtype=int)
    labels = np.asarray(labels, dtype=str)[indices]
    blocks = np.unique(labels)
    if len(blocks) < 3:
        raise ValueError(f"too few connected spatial blocks for {band}: {len(blocks)}")
    by_block = {b: np.flatnonzero(labels == b) for b in blocks}
    baseline = np.asarray([row["baseline_icl_delta"] for row in selected], dtype=float)
    strict = np.asarray([row["strict_icl_delta"] for row in selected], dtype=float)
    satellite = np.asarray([row["satellite_delta_50_75"] for row in selected], dtype=float)
    z = np.asarray([row["redshift"] for row in selected], dtype=float)
    richness = np.asarray([row["richness"] for row in selected], dtype=float)
    point_base = partial_pearson(baseline, satellite, [z, richness])
    point_strict = partial_pearson(strict, satellite, [z, richness])
    rng = np.random.default_rng(seed)
    draws = np.full((replicas, 3), np.nan)
    for j in range(replicas):
        chosen = rng.choice(blocks, size=len(blocks), replace=True)
        idx = np.concatenate([by_block[b] for b in chosen])
        draws[j, 0] = partial_pearson(baseline[idx], satellite[idx], [z[idx], richness[idx]])
        draws[j, 1] = partial_pearson(strict[idx], satellite[idx], [z[idx], richness[idx]])
        draws[j, 2] = draws[j, 1] - draws[j, 0]
    good = np.isfinite(draws).all(axis=1)
    if good.sum() < replicas*0.9:
        raise ValueError(f"too few finite bootstrap replicates for {band}")
    return {"n_pairs": len(selected), "n_connected_blocks": len(blocks),
            "baseline_partial_r": float(point_base), "strict_partial_r": float(point_strict),
            "strict_minus_baseline_partial_r": float(point_strict-point_base),
            "baseline_ci95": np.quantile(draws[good,0],[.025,.975]).tolist(),
            "strict_ci95": np.quantile(draws[good,1],[.025,.975]).tolist(),
            "paired_delta_ci95": np.quantile(draws[good,2],[.025,.975]).tolist(),
            "replicas_requested": replicas, "replicas_valid": int(good.sum()), "seed": seed}


def run(output: Path, replicas: int, seed: int) -> dict:
    output = Path(output).resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite {output}")
    for path in (PAIR_TABLE, PROFILE_TABLE, TARGET_MASK_RECORDS,
                 ROOT/"sample_manifest.json", ROOT/"control_download_01/download_manifest.json"):
        if not path.is_file():
            raise FileNotFoundError(path)
    plan, download, checkpoint, manifests, input_hashes = warp_noise.load_inputs()
    pair_rows = read_tsv(PAIR_TABLE)
    profile_rows = read_tsv(PROFILE_TABLE)
    mask_records = {r["id"]: r for r in read_tsv(TARGET_MASK_RECORDS)}
    sample = json.loads((ROOT/"sample_manifest.json").read_text(encoding="utf-8"))
    base_profiles = {(r["id"], r["filter"]): r for r in profile_rows
                     if r["kron_factor"] == "2.5" and r["coverage_floor"] == "0.4"
                     and r["radial_bin_kpc"] == "200-300"}
    manifests_map = warp_noise.pair_inventory(plan, pair_rows, manifests)
    request_by_id = {r["request_id"]: r for r in plan["requests"]}
    request_index = {r["request_id"]: i for i,r in enumerate(plan["requests"])}
    strict_maps, support_meta, failures = {}, {}, []
    for request_id, request in request_by_id.items():
        ix = request_index[request_id]
        record = checkpoint["requests"].get(f"request_{ix:03d}", {})
        if record.get("status") == "http_404_permanent":
            failures.append(request_id)
            continue
        try:
            strict, common_names, meta = strict_warp_support(
                request, manifests_map[request_id], DOWNLOAD_DIR/f"request_{ix:03d}.tar", record)
            strict_maps[request_id] = strict
            support_meta[request_id] = {**meta, "common_mask_planes": common_names,
                                        "archive_sha256": record.get("archive_sha256")}
        except Exception as exc:
            failures.append({"request_id": request_id, "error": f"{type(exc).__name__}: {exc}"})

    # Only compare pairs with both complete target and control warp support.
    rows, mask_hashes = [], {}
    target_manifest = {str(r["id"]): r for r in manifests["target"]["clusters"]}
    control_manifest = {str(r["id"]): r for r in manifests["control"]["clusters"]}
    for pair in pair_rows:
        band, ident, control_id = pair["filter"], str(pair["cluster_id"]), str(pair["control_id"])
        tid, cid = f"{band}:target:{ident}", f"{band}:control:{control_id}"
        if tid not in strict_maps or cid not in strict_maps:
            continue
        original = base_profiles.get((ident, band))
        if original is None or str(original["control_id"]) != control_id:
            raise ValueError(f"baseline selected-profile row missing or mismatched: {ident}/{band}")
        role_results = {}
        for role, req_id, frame, center_ra, center_dec, source_id in (
            ("target", tid, target_manifest[ident], float(pair["target_ra_deg"]),
             float(pair["target_dec_deg"]), ident),
            ("control", cid, control_manifest[control_id], float(pair["control_ra_deg"]),
             float(pair["control_dec_deg"]), control_id),
        ):
            image_dir = ROOT/("download_01/images" if role == "target" else "control_download_01/images")
            image_path = image_dir/frame["image"]
            source_mask, mask_path, mask_sha = load_source_mask(role, ident, control_id,
                                                                   tuple(strict_maps[req_id].shape), mask_records)
            mask_hashes[mask_path] = mask_sha
            if source_mask.shape != strict_maps[req_id].shape:
                raise ValueError(f"source mask shape mismatch for {req_id}")
            baseline = reconstruct_profile(image_path, source_mask, center_ra, center_dec,
                                           float(pair["redshift"]), sample)
            expected_flux = float(original["target_flux"] if role == "target" else original["control_flux"])
            expected_cov = float(original["target_coverage"] if role == "target" else original["control_coverage"])
            if (not math.isclose(baseline["flux_njy_arcsec2"], expected_flux, rel_tol=0, abs_tol=2e-6)
                    or not math.isclose(baseline["coverage"], expected_cov, rel_tol=0, abs_tol=2e-6)):
                raise ValueError(f"baseline flux/coverage does not reproduce for {req_id}")
            strict_mask = source_mask | ~strict_maps[req_id]
            strict_profile = reconstruct_profile(image_path, strict_mask, center_ra, center_dec,
                                                 float(pair["redshift"]), sample)
            role_results[role] = {"baseline": baseline, "strict": strict_profile,
                                  "strict_support": support_meta[req_id]}
        t, c = role_results["target"], role_results["control"]
        baseline_delta = t["baseline"]["flux_njy_arcsec2"] - c["baseline"]["flux_njy_arcsec2"]
        strict_delta = t["strict"]["flux_njy_arcsec2"] - c["strict"]["flux_njy_arcsec2"]
        if not math.isclose(baseline_delta, float(pair["icl_delta_200_300"]), rel_tol=0, abs_tol=2e-6):
            raise ValueError(f"baseline pair delta does not reproduce: {ident}/{band}")
        rows.append({"cluster_id": ident, "filter": band, "control_id": control_id,
                     "baseline_icl_delta": baseline_delta, "strict_icl_delta": strict_delta,
                     "strict_minus_baseline_icl_delta": strict_delta-baseline_delta,
                     "satellite_delta_50_75": float(pair["satellite_delta_50_75"]),
                     "redshift": float(pair["redshift"]), "richness": float(pair["richness"]),
                     "target_baseline_coverage": t["baseline"]["coverage"],
                     "control_baseline_coverage": c["baseline"]["coverage"],
                     "target_strict_coverage": t["strict"]["coverage"],
                     "control_strict_coverage": c["strict"]["coverage"],
                     "target_strict_support_fraction": t["strict_support"]["strict_valid_pixel_fraction_of_coadd"],
                     "control_strict_support_fraction": c["strict_support"]["strict_valid_pixel_fraction_of_coadd"],
                     "target_ra_deg": float(pair["target_ra_deg"]),
                     "target_dec_deg": float(pair["target_dec_deg"]),
                     "control_ra_deg": float(pair["control_ra_deg"]),
                     "control_dec_deg": float(pair["control_dec_deg"]),
                     "connected_5deg_block": pair["connected_5deg_block"],
                     "science_measurement": False})

    # Keep exactly the original >=40% annulus-coverage rule after strict masking.
    accepted = [r for r in rows if r["target_strict_coverage"] >= .40
                and r["control_strict_coverage"] >= .40]
    labels = np.empty(len(accepted), dtype=object)
    strict_blocks = {}
    for band in ("i", "i2"):
        indices = np.asarray([i for i, row in enumerate(accepted)
                              if row["filter"] == band], dtype=int)
        selected = [accepted[i] for i in indices]
        target_coords = [{"ra_deg": r["target_ra_deg"], "dec_deg": r["target_dec_deg"]}
                         for r in selected]
        control_coords = [{"ra_deg": r["control_ra_deg"], "dec_deg": r["control_dec_deg"]}
                          for r in selected]
        _, band_labels, target_cells, control_cells, nblocks, largest = pair_block_labels(
            target_coords, control_coords, np.ones(len(selected), dtype=bool))
        labels[indices] = np.asarray(band_labels, dtype=str)
        strict_blocks[band] = {"n_connected_5deg_blocks": int(nblocks),
                               "n_target_cells": int(target_cells),
                               "n_control_cells": int(control_cells),
                               "largest_block_pairs": int(largest)}
        for index, label in zip(indices, band_labels):
            accepted[int(index)]["strict_mask_connected_5deg_block"] = str(label)
    by_filter = {band: block_bootstrap(accepted, labels, band, replicas,
                                      seed+(0 if band=="i" else 1))
                 for band in ("i", "i2")}
    per_band_coverage = {}
    for band in ("i", "i2"):
        group = [r for r in accepted if r["filter"] == band]
        per_band_coverage[band] = {
            "n_warp_complete_pairs": sum(r["filter"] == band for r in rows),
            "n_strict_pairs_after_original_40pct_floor": len(group),
            "median_target_strict_coverage": float(np.median([r["target_strict_coverage"] for r in group])) if group else None,
            "median_control_strict_coverage": float(np.median([r["control_strict_coverage"] for r in group])) if group else None,
            "median_target_minus_control_icl_change": float(np.median([r["strict_minus_baseline_icl_delta"] for r in group])) if group else None,
        }
    output.mkdir(parents=True, exist_ok=False)
    pair_path = output/"strict_mask_association_pairs.tsv"
    write_tsv(pair_path, rows)
    summary = {
        "schema": "astrafek-icl-supported94-warp-mask-association-v1",
        "status": "complete_exploratory_sensitivity" if not failures else "complete_with_data_gaps",
        "science_measurement": False, "hsc_network_used": False, "images_read": True,
        "approved_plan_sha256": warp_noise.PLAN_SHA256,
        "download_summary_sha256": warp_noise.DOWNLOAD_SUMMARY_SHA256,
        "planned_pairs": 94, "planned_centers": 188,
        "warp_complete_pair_count_before_mask_coverage_floor": len(rows),
        "strict_pair_count_after_original_40pct_floor": len(accepted),
        "failed_requests_or_centers": failures,
        "per_filter_coverage": per_band_coverage,
        "strict_mask_blocks_by_filter": strict_blocks,
        "association_by_filter": by_filter,
        "method": "Recompute raw 200-300 proper-kpc ICL target-minus-control from coadd IMAGE under original saved K=2.5 masks, intersected with coadd pixels having at least three valid visits and only mask planes shared by coadd and all returned warps; retain the original 40% target/control annular coverage floor; compare baseline and strict masks on identical retained pairs and use connected 5-degree block bootstrap.",
        "limitations": [
            "Warp-to-coadd mapping uses a local affine fit and nearest-neighbor mask sampling, not the exact HSC resampling kernel.",
            "The strict mask is a quality/support sensitivity, not a full propagated covariance, source-Poisson, sky-transfer, PSF-wing, or completeness correction.",
            "HTTP 404 center gaps and strict-mask coverage selection can change the retained sample; report both counts and coverages.",
            "The satellite profile remains the saved 50-75 kpc contrast; only the ICL annulus pixel support is recomputed.",
            "Keep science_measurement=false; this analysis cannot establish causality or validate the 800-parent population.",
        ],
        "input_sha256": {**input_hashes,
                         str(PAIR_TABLE): warp_noise.sha256(PAIR_TABLE),
                         str(PROFILE_TABLE): warp_noise.sha256(PROFILE_TABLE),
                         str(TARGET_MASK_RECORDS): warp_noise.sha256(TARGET_MASK_RECORDS),
                         **mask_hashes},
        "analysis_script_sha256": warp_noise.sha256(Path(__file__).resolve()),
    }
    summary_path = output/"summary.json"
    summary_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2,
                                       sort_keys=True, allow_nan=False)+"\n", encoding="utf-8")
    summary["summary_sha256"] = warp_noise.sha256(summary_path)
    print(json.dumps({"status": summary["status"], "science_measurement": False,
                      "warp_complete_pairs": len(rows), "strict_pairs": len(accepted),
                      "by_filter": by_filter, "summary": str(summary_path),
                      "pair_table": str(pair_path)}, ensure_ascii=False, indent=2))
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--replicas", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=20261106)
    args = parser.parse_args()
    run(args.output, args.replicas, args.seed)


if __name__ == "__main__":
    main()
