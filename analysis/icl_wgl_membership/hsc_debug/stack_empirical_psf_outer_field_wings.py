#!/usr/bin/env python3
"""Measure extended stellar wings in existing HSC target/control coadds.

This is a matched, frame-level diagnostic only. It does not certify an effective
coadd PSF or produce science-corrected ICL photometry.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import sys

import numpy as np
from astropy.cosmology import FlatLambdaCDM
from astropy.io import fits

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
if str(HERE.parents[2]) not in sys.path:
    sys.path.insert(0, str(HERE.parents[2]))

import prepare_detection_psf as prep
from diagnose_empirical_psf import objects_from_rows, read_table
from fetch_ancillary_catalogues import sha

BASE = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
ARCSEC_PER_RADIAN = 206264.80624709636
RADIAL_EDGES_KPC = np.arange(0.0, 325.0, 25.0)
SKY_ANNULI_PIX = ((120.0, 150.0), (160.0, 190.0))
N_SECTORS = 8


def read_tsv(path: Path) -> list[dict[str, str]]:
    with Path(path).open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def kpc_per_pixel(z: float, pixel_scale_arcsec: float, cosmology: FlatLambdaCDM) -> float:
    if not np.isfinite([z, pixel_scale_arcsec]).all() or z <= 0 or pixel_scale_arcsec <= 0:
        raise ValueError("redshift and pixel scale must be finite and positive")
    return float(cosmology.angular_diameter_distance(z).to_value("kpc")
                 * pixel_scale_arcsec / ARCSEC_PER_RADIAN)


def edge_radius(shape: tuple[int, int], x: float, y: float) -> float:
    """Inscribed circular radius, in pixels, for a fractional image position."""
    ny, nx = shape
    if not np.isfinite([x, y]).all() or not (0 <= x < nx and 0 <= y < ny):
        return -1.0
    return float(min(x, y, nx - 1 - x, ny - 1 - y))


def select_candidates(sample: Path, *, outer_radius_kpc: float = 800.0) -> tuple[list[dict], dict]:
    """Select isolated Gaia stars with full 300-kpc target/control footprints."""
    sample = Path(sample).resolve()
    target_manifest_path = sample / "download_01/download_manifest.json"
    target_status_path = sample / "download_01/download_status.json"
    control_manifest_path = sample / "control_download_01/download_manifest.json"
    control_status_path = sample / "control_download_01/download_status.json"
    ledger_path = sample / "empirical_psf_gaia_full_sample_2026-10-20_02/star_ledger.tsv"
    quality_path = sample / "empirical_psf_star_stamp_quality_audit_2026-10-20_02/star_stamp_quality.tsv"
    support_path = sample / "empirical_psf_matched_control_support_2026-10-20_01/matched_control_support.tsv"
    source_paths = (target_manifest_path, target_status_path, control_manifest_path,
                    control_status_path, ledger_path, quality_path, support_path)
    if not all(path.is_file() for path in source_paths):
        raise FileNotFoundError("one or more existing HSC/PSF input ledgers are missing")
    target_manifest = json.loads(target_manifest_path.read_text(encoding="utf-8"))
    target_status = json.loads(target_status_path.read_text(encoding="utf-8"))
    control_manifest = json.loads(control_manifest_path.read_text(encoding="utf-8"))
    control_status = json.loads(control_status_path.read_text(encoding="utf-8"))
    if target_manifest.get("science_measurement") is not False or control_manifest.get("science_measurement") is not False:
        raise ValueError("input manifests have unexpected science-measurement status")
    h0, omega_m = float(target_manifest["h0"]), float(target_manifest["omega_m"])
    cosmology = FlatLambdaCDM(H0=h0, Om0=omega_m, Tcmb0=0)
    targets = {str(row["id"]): row for row in target_manifest["clusters"]}
    target_done = {str(row["id"]): row for row in target_status["completed"]}
    controls = {str(row["id"]): row for row in control_manifest["clusters"]}
    control_done = {str(row["id"]): row for row in control_status["completed"]}
    ledger = read_tsv(ledger_path)
    quality = {str(row["candidate_id"]): row for row in read_tsv(quality_path)}
    support = {str(row["candidate_id"]): row for row in read_tsv(support_path)}
    if len(targets) != 800 or len(target_done) != 800 or len(controls) != 798 or len(control_done) != 798:
        raise ValueError("the pinned 800-target/798-control data set is incomplete")
    if len(ledger) != 480 or sum(bool(row.get("stamp")) for row in ledger) != 333:
        raise ValueError("unexpected Gaia full-sample ledger cardinality")
    if len(quality) != 333 or len(support) != 333:
        raise ValueError("quality/support tables do not cover all 333 complete Gaia stamps")
    counts = {key: 0 for key in ("complete_stamps", "known_nonvariable_or_unclassified", "outer_field",
                                  "target_300kpc_footprint", "same_filter_control", "control_300kpc_footprint",
                                  "k2p5_core_99pct", "selected")}
    selected = []
    for star in ledger:
        if not star.get("stamp"):
            continue
        counts["complete_stamps"] += 1
        # No known Gaia VARIABLE source; UNKNOWN/NOT_AVAILABLE remain flagged in provenance.
        if str(star.get("gaia_variability", "")).upper() == "VARIABLE":
            continue
        counts["known_nonvariable_or_unclassified"] += 1
        cid = str(star["cluster_id"])
        target = targets[cid]
        target_record = target_done[cid]
        target_review = target_record["review"]
        z = float(target["z"])
        pixel_scale = float(target_review["pixel_scale_arcsec"])
        kpp = kpc_per_pixel(z, pixel_scale, cosmology)
        x, y = float(star["native_x"]), float(star["native_y"])
        cx, cy = (float(v) for v in target_review["center_pixel_xy"])
        cluster_radius_kpc = math.hypot(x - cx, y - cy) * kpp
        if cluster_radius_kpc <= outer_radius_kpc:
            continue
        counts["outer_field"] += 1
        target_shape = tuple(int(v) for v in target_review["shape"])
        if edge_radius(target_shape, x, y) * kpp < 300.0:
            continue
        counts["target_300kpc_footprint"] += 1
        pair = support[str(star["candidate_id"])]
        if (pair["pair_status"] != "same_filter_offset_available"
                or pair["target_filter"] != pair["control_filter"]):
            continue
        counts["same_filter_control"] += 1
        control_id = str(pair["control_id"])
        control = controls[control_id]
        control_record = control_done[control_id]
        control_review = control_record["review"]
        if control.get("target_cluster_id") != cid:
            raise ValueError(f"matched control belongs to a different target for star {star['candidate_id']}")
        if control_review["coadd_header_filter"] != target_review["coadd_header_filter"]:
            raise ValueError(f"matched pair filter mismatch for star {star['candidate_id']}")
        if abs(float(control_review["pixel_scale_arcsec"]) - pixel_scale) > 1e-5:
            raise ValueError(f"matched pair pixel scales differ for star {star['candidate_id']}")
        control_x, control_y = float(pair["control_x"]), float(pair["control_y"])
        control_shape = tuple(int(v) for v in control_review["shape"])
        if edge_radius(control_shape, control_x, control_y) * kpp < 300.0:
            continue
        counts["control_300kpc_footprint"] += 1
        q = quality[str(star["candidate_id"])]
        if int(q["detection_flags"]) != 0 or int(q["hsc_bright_object_pixels_core_0_12"]) != 0:
            continue
        if float(q["mask_k2p5_core_0_12_open_fraction"]) < 0.99:
            continue
        counts["k2p5_core_99pct"] += 1
        selected.append({
            "candidate_id": str(star["candidate_id"]), "cluster_id": cid,
            "detection_number": int(star["number"]), "source_id": str(star["source_id"]),
            "gaia_variability": str(star["gaia_variability"]), "epoch_review_required": True,
            "target_image": str((sample / "download_01/images" / target_record["file"]).resolve()),
            "control_id": control_id,
            "control_image": str((sample / "control_download_01/images" / control_record["file"]).resolve()),
            "target_x": x, "target_y": y, "control_x": control_x, "control_y": control_y,
            "z_lens": z, "cluster_centric_radius_kpc": cluster_radius_kpc,
            "pixel_scale_arcsec": pixel_scale, "kpc_per_pixel": kpp,
            "target_filter": str(target_review["coadd_header_filter"]),
            "k2p5_core_open_fraction": float(q["mask_k2p5_core_0_12_open_fraction"]),
            "target_shape": target_shape, "control_shape": control_shape,
            "target_status_row": target_record, "control_status_row": control_record,
        })
    if len({row["cluster_id"] for row in selected}) != len(selected):
        raise ValueError("selection contains multiple stars in one HSC target frame")
    counts["selected"] = len(selected)
    if not selected:
        raise ValueError("no matched target/control Gaia star candidates passed the geometry gates")
    provenance = {
        "sample_root": str(sample), "cosmology": {"H0_km_s_Mpc": h0, "Omega_m": omega_m, "Tcmb0_K": 0.0},
        "input_sha256": {str(path.resolve()): sha(path) for path in source_paths},
        "counts": counts,
    }
    return selected, provenance


def crop_bounds(shape: tuple[int, int], x: float, y: float, radius_pix: float) -> tuple[int, int, int, int]:
    ny, nx = shape
    if edge_radius(shape, x, y) < radius_pix:
        raise ValueError("requested radial profile extends beyond the image edge")
    x0, x1 = max(0, int(math.floor(x - radius_pix))), min(nx, int(math.ceil(x + radius_pix)) + 1)
    y0, y1 = max(0, int(math.floor(y - radius_pix))), min(ny, int(math.ceil(y + radius_pix)) + 1)
    return x0, x1, y0, y1


def radial_profile(image: np.ndarray, variance: np.ndarray, mask: np.ndarray, center: tuple[float, float],
                   kpc_per_pix: float, core_flux: float | dict[tuple[float, float], float] | None = None) -> list[dict]:
    """Return 25-kpc annular surface brightness relative to the target core flux."""
    image = np.asarray(image, dtype=float)
    variance = np.asarray(variance, dtype=float)
    mask = np.asarray(mask, dtype=bool)
    if image.ndim != 2 or variance.shape != image.shape or mask.shape != image.shape:
        raise ValueError("science, variance, and mask arrays must have matching 2-D shapes")
    if not np.isfinite(kpc_per_pix) or kpc_per_pix <= 0:
        raise ValueError("physical pixel scale must be finite and positive")
    if core_flux is not None and not isinstance(core_flux, dict) and (not np.isfinite(core_flux) or core_flux <= 0):
        raise ValueError("target core flux must be finite and positive")
    yy, xx = np.indices(image.shape, dtype=float)
    dx, dy = xx - float(center[0]), yy - float(center[1])
    radius_pix = np.hypot(dx, dy)
    radius_kpc = radius_pix * kpc_per_pix
    angle = np.arctan2(dy, dx)
    sector = np.clip(np.floor((angle + np.pi) * (N_SECTORS / (2 * np.pi))).astype(int), 0, N_SECTORS - 1)
    good = (~mask & np.isfinite(image) & np.isfinite(variance) & (variance > 0))
    output = []
    for sky_lo, sky_hi in SKY_ANNULI_PIX:
        sky_region = (radius_pix >= sky_lo) & (radius_pix < sky_hi)
        sky_good = sky_region & good
        sector_means, sector_fractions, sector_counts = [], [], []
        for sid in range(N_SECTORS):
            in_sector = sky_region & (sector == sid)
            valid = in_sector & good
            count = int(valid.sum())
            sector_counts.append(count)
            sector_fractions.append(count / max(1, int(in_sector.sum())))
            if count:
                sector_means.append(float(np.mean(image[valid])))
        sky_fraction = float(sky_good.sum() / max(1, sky_region.sum()))
        sky_ok = (int(sky_good.sum()) >= 300 and sky_fraction >= 0.30
                  and min(sector_counts) >= 20 and min(sector_fractions) >= 0.25)
        background = float(np.mean(sector_means)) if sky_ok else None
        core_region = radius_pix < 12.0
        core_good = core_region & good
        core_fraction = float(core_good.sum() / max(1, core_region.sum()))
        target_flux = (core_flux.get((sky_lo, sky_hi)) if isinstance(core_flux, dict) else core_flux)
        if target_flux is None and core_flux is None and sky_ok and core_fraction >= 0.99:
            target_flux = float(np.sum(image[core_good] - background))
        for lo, hi in zip(RADIAL_EDGES_KPC[:-1], RADIAL_EDGES_KPC[1:]):
            ring = (radius_kpc >= lo) & (radius_kpc < hi)
            valid = ring & good
            area = int(ring.sum())
            n_valid = int(valid.sum())
            area_fraction = n_valid / max(1, area)
            n_sectors = sum(
                np.count_nonzero(valid & (sector == sid))
                >= max(1, 0.25 * np.count_nonzero(ring & (sector == sid)))
                for sid in range(N_SECTORS)
            )
            status = "ok"
            if not sky_ok:
                status = "blocked_background_coverage"
            elif lo < 12.0 * kpc_per_pix and core_fraction < 0.99:
                status = "blocked_core_coverage"
            elif target_flux is None or not np.isfinite(target_flux) or target_flux <= 0:
                status = "blocked_target_core_flux"
            elif area_fraction < 0.50 or n_sectors != N_SECTORS:
                status = "blocked_radial_coverage"
            value = None
            raw_mean = float(np.mean(image[valid])) if n_valid else None
            if status == "ok":
                value = (raw_mean - background) / target_flux
            output.append({
                "sky_r_min_pix": sky_lo, "sky_r_max_pix": sky_hi,
                "r_min_kpc": float(lo), "r_max_kpc": float(hi),
                "r_mid_kpc": float((lo + hi) / 2), "n_valid": n_valid,
                "area_fraction": area_fraction, "n_covered_sectors": n_sectors,
                "core_open_fraction": core_fraction, "background": background,
                "target_core_flux": target_flux, "raw_mean": raw_mean,
                "profile_per_target_core_per_pixel": value, "status": status,
                "science_measurement": False, "psf_approved": False,
            })
    return output


def _cutout_plane(hdu, bounds: tuple[int, int, int, int]) -> np.ndarray:
    x0, x1, y0, y1 = bounds
    return np.asarray(hdu.section[y0:y1, x0:x1])


def _control_bad_bits(manifest_row: dict) -> int:
    bits = set(manifest_row["bad_mask"]["exclude_bits"])
    bits.update(manifest_row["bright_star_mask"]["exclude_bits"])
    value = 0
    for bit in bits:
        value |= int(bit)
    return value


def measure_candidate(candidate: dict, sample: Path, detection_dir: Path) -> list[dict]:
    """Read and measure one existing target/control pair, verifying source hashes."""
    target_row = candidate["target_status_row"]
    control_row = candidate["control_status_row"]
    if sha(Path(candidate["target_image"])) != target_row["sha256"]:
        raise ValueError(f"target image SHA-256 mismatch for {candidate['cluster_id']}")
    if sha(Path(candidate["control_image"])) != control_row["sha256"]:
        raise ValueError(f"control image SHA-256 mismatch for {candidate['control_id']}")
    rmax_pix = max(300.0 / candidate["kpc_per_pixel"], SKY_ANNULI_PIX[-1][1]) + 2
    target_bounds = crop_bounds(tuple(candidate["target_shape"]), candidate["target_x"], candidate["target_y"], rmax_pix)
    control_bounds = crop_bounds(tuple(candidate["control_shape"]), candidate["control_x"], candidate["control_y"], rmax_pix)
    record_rows = []
    target_path = Path(candidate["target_image"])
    control_path = Path(candidate["control_image"])
    det_path = detection_dir / f"cluster_{candidate['cluster_id']}"
    detections_path = det_path / "detections.tsv"
    masks_path = det_path / "masks.fits.gz"
    if not detections_path.is_file() or not masks_path.is_file():
        raise FileNotFoundError(f"missing prepared target masks for cluster {candidate['cluster_id']}")
    detections = read_table(detections_path)
    number = candidate["detection_number"]
    if number < 1 or number > len(detections) or int(detections[number - 1]["NUMBER"]) != number:
        raise ValueError(f"HSC detection identity changed for candidate {candidate['candidate_id']}")
    det_row = detections[number - 1]
    detected_xy = (float(det_row["X_IMAGE"]) - 1.0, float(det_row["Y_IMAGE"]) - 1.0)
    if not np.allclose(detected_xy, [candidate["target_x"], candidate["target_y"]], rtol=0, atol=1e-6):
        raise ValueError(f"target star position differs from its HSC detection: {candidate['candidate_id']}")
    kpp = candidate["kpc_per_pixel"]
    with fits.open(target_path, memmap=True) as target_hdus, fits.open(masks_path, memmap=False) as mask_hdus:
        review = target_row["review"]
        sci = target_hdus[int(review["science_hdu"])]
        var = target_hdus[int(review["variance_hdu"])]
        bits_hdu = target_hdus[int(review["mask_hdu"])]
        labels = np.asarray(mask_hdus["SEGMENTATION"].data)
        if tuple(sci.shape) != tuple(candidate["target_shape"]) or labels.shape != sci.shape:
            raise ValueError(f"target image/segmentation dimensions differ for {candidate['cluster_id']}")
        bad_bits = sum(1 << int(bit) for bit in prep.BAD_BIT_POSITIONS) | int(prep.BRIGHT_OBJECT_BIT)
        x0, x1, y0, y1 = target_bounds
        target_image = _cutout_plane(sci, target_bounds).astype(float)
        target_variance = _cutout_plane(var, target_bounds).astype(float)
        target_bits = _cutout_plane(bits_hdu, target_bounds)
        target_base = (target_bits & bad_bits) != 0
        target_base |= ~np.isfinite(target_image) | ~np.isfinite(target_variance) | (target_variance <= 0)
        target_labels = labels[y0:y1, x0:x1]
        objects, kron, kron_flags = objects_from_rows(detections, x0, y0)
        masks, _ = prep.source_masks(objects, target_labels, kron, kron_flags, target_base,
                                     number - 1, np.empty((0, 2)), candidate["pixel_scale_arcsec"])
        target_mask = masks[0]
        target_center = (candidate["target_x"] - x0, candidate["target_y"] - y0)
        target_profile = radial_profile(target_image, target_variance, target_mask, target_center, kpp)
    core_flux_by_sky = {}
    for row in target_profile:
        key = (row["sky_r_min_pix"], row["sky_r_max_pix"])
        if row["r_min_kpc"] == 0.0 and row["status"] == "ok":
            core_flux_by_sky[key] = row["target_core_flux"]
    with fits.open(control_path, memmap=True) as control_hdus:
        review = control_row["review"]
        sci = control_hdus[int(review["science_hdu"])]
        var = control_hdus[int(review["variance_hdu"])]
        bits_hdu = control_hdus[int(review["mask_hdu"])]
        control_image = _cutout_plane(sci, control_bounds).astype(float)
        control_variance = _cutout_plane(var, control_bounds).astype(float)
        control_bits = _cutout_plane(bits_hdu, control_bounds)
    control_manifest_path = sample / "control_download_01/download_manifest.json"
    control_manifest = json.loads(control_manifest_path.read_text(encoding="utf-8"))
    control_manifest_row = next(row for row in control_manifest["clusters"] if row["id"] == candidate["control_id"])
    control_base = (control_bits & _control_bad_bits(control_manifest_row)) != 0
    control_base |= ~np.isfinite(control_image) | ~np.isfinite(control_variance) | (control_variance <= 0)
    control_center = (candidate["control_x"] - control_bounds[0], candidate["control_y"] - control_bounds[2])
    # Keep control detections local to the cutout; all detected objects are masked.
    control_mask_status = "ok"
    try:
        objects, segmentation, kron, kron_flags, _, _, _ = prep.detect(
            control_image, control_variance, control_base, np.asarray(control_center),
            np.empty((0, 2)), candidate["pixel_scale_arcsec"],
        )
    except ValueError as error:
        if "Too few valid pixels" not in str(error):
            raise
        control_masks = [np.ones(control_image.shape, dtype=bool)] * len(prep.FACTORS)
        control_mask_status = "blocked_control_too_few_quality_pixels"
    else:
        control_masks, _ = prep.source_masks(
            objects, segmentation, kron, kron_flags, control_base, -1,
            np.empty((0, 2)), candidate["pixel_scale_arcsec"],
        )
    control_mask = control_masks[0]
    control_profile = radial_profile(control_image, control_variance, control_mask,
                                     control_center, kpp, core_flux=core_flux_by_sky)
    for role, rows in (("target", target_profile), ("control", control_profile)):
        for row in rows:
            output = {
                "candidate_id": candidate["candidate_id"], "cluster_id": candidate["cluster_id"],
                "control_id": candidate["control_id"], "filter": candidate["target_filter"],
                "role": role, "z_lens": candidate["z_lens"],
                "cluster_centric_radius_kpc": candidate["cluster_centric_radius_kpc"],
                "gaia_variability": candidate["gaia_variability"], "epoch_review_required": True,
                "kpc_per_pixel": kpp, **row,
            }
            if role == "control" and control_mask_status != "ok" and output["status"] == "ok":
                output["status"] = control_mask_status
                output["profile_per_target_core_per_pixel"] = None
            record_rows.append(output)
    return record_rows


def stack_pairs(rows: list[dict], *, n_bootstrap: int = 10000, seed: int = 20261020) -> list[dict]:
    """Frame-equal paired target-minus-control stack with frame bootstrap ranges."""
    grouped: dict[tuple, dict[str, dict[str, float]]] = {}
    for row in rows:
        if row["status"] != "ok" or row["profile_per_target_core_per_pixel"] is None:
            continue
        key = (int(row["sky_r_min_pix"]), int(row["sky_r_max_pix"]),
               float(row["r_min_kpc"]), float(row["r_max_kpc"]))
        grouped.setdefault(key, {}).setdefault(row["candidate_id"], {})[row["role"]] = float(row["profile_per_target_core_per_pixel"])
    rng = np.random.default_rng(seed)
    summary = []
    keys = sorted({(int(r["sky_r_min_pix"]), int(r["sky_r_max_pix"]),
                    float(r["r_min_kpc"]), float(r["r_max_kpc"])) for r in rows})
    for key in keys:
        pairs = grouped.get(key, {})
        complete = [(cid, values) for cid, values in sorted(pairs.items()) if {"target", "control"} <= values.keys()]
        target = np.asarray([v["target"] for _, v in complete], dtype=float)
        control = np.asarray([v["control"] for _, v in complete], dtype=float)
        difference = target - control
        if len(difference) and n_bootstrap > 0:
            draws = rng.integers(0, len(difference), size=(n_bootstrap, len(difference)))
            boot = np.mean(difference[draws], axis=1)
            low, high = (float(v) for v in np.quantile(boot, (0.16, 0.84)))
        else:
            low = high = None
        summary.append({
            "sky_r_min_pix": key[0], "sky_r_max_pix": key[1],
            "r_min_kpc": key[2], "r_max_kpc": key[3], "r_mid_kpc": (key[2] + key[3]) / 2,
            "n_paired_frames": int(len(difference)),
            "target_frame_equal_mean": float(np.mean(target)) if len(target) else None,
            "control_frame_equal_mean": float(np.mean(control)) if len(control) else None,
            "paired_difference_frame_equal_mean": float(np.mean(difference)) if len(difference) else None,
            "paired_difference_frame_median": float(np.median(difference)) if len(difference) else None,
            "paired_difference_bootstrap_q16": low, "paired_difference_bootstrap_q84": high,
            "bootstrap_unit": "matched HSC target/control frame pair",
            "science_measurement": False, "psf_approved": False,
        })
    return summary


def write_tsv(path: Path, rows: list[dict]) -> None:
    if not rows:
        Path(path).write_text("", encoding="utf-8")
        return
    with Path(path).open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t", extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def run(sample: Path, output: Path, *, bootstrap_draws: int = 10000, seed: int = 20261020,
        outer_radius_kpc: float = 800.0) -> dict:
    sample, output = Path(sample).resolve(), Path(output).resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite {output}")
    candidates, provenance = select_candidates(sample, outer_radius_kpc=outer_radius_kpc)
    ledger_summary_path = sample / "empirical_psf_gaia_full_sample_2026-10-20_02/summary.json"
    ledger_summary = json.loads(ledger_summary_path.read_text(encoding="utf-8"))
    detection_dir = Path(ledger_summary["detection_dir"])
    profiles = []
    input_hashes = dict(provenance["input_sha256"])
    input_hashes[str(ledger_summary_path.resolve())] = sha(ledger_summary_path)
    input_hashes[str((sample / "empirical_psf_star_stamp_quality_audit_2026-10-20_02/summary.json").resolve())] = sha(
        sample / "empirical_psf_star_stamp_quality_audit_2026-10-20_02/summary.json")
    for candidate in candidates:
        measured = measure_candidate(candidate, sample, detection_dir)
        profiles.extend(measured)
        input_hashes[candidate["target_image"]] = candidate["target_status_row"]["sha256"]
        input_hashes[candidate["control_image"]] = candidate["control_status_row"]["sha256"]
        frame_dir = detection_dir / f"cluster_{candidate['cluster_id']}"
        for path in (frame_dir / "detections.tsv", frame_dir / "masks.fits.gz"):
            input_hashes[str(path.resolve())] = sha(path)
        print(f"measured candidate={candidate['candidate_id']} cluster={candidate['cluster_id']} filter={candidate['target_filter']}", flush=True)
    stack = stack_pairs(profiles, n_bootstrap=bootstrap_draws, seed=seed)
    output.mkdir(parents=True, exist_ok=False)
    candidate_rows = [{key: value for key, value in row.items()
                       if key not in ("target_status_row", "control_status_row")} for row in candidates]
    write_tsv(output / "selected_candidates.tsv", candidate_rows)
    write_tsv(output / "matched_radial_profiles.tsv", profiles)
    write_tsv(output / "frame_equal_paired_stack.tsv", stack)
    files = ("selected_candidates.tsv", "matched_radial_profiles.tsv", "frame_equal_paired_stack.tsv")
    summary = {
        "schema": "astrafek-empirical-psf-outer-field-wings-v1", "status": "complete_exploratory_diagnostic",
        "counts": provenance["counts"], "n_profile_rows": len(profiles),
        "n_complete_pair_bins": sum(row["n_paired_frames"] > 0 for row in stack),
        "n_bootstrap_draws": int(bootstrap_draws), "seed": int(seed),
        "radial_edges_kpc": RADIAL_EDGES_KPC.tolist(), "sky_annuli_pix": [list(v) for v in SKY_ANNULI_PIX],
        "selected_candidate_ids": [row["candidate_id"] for row in candidates],
        "input_sha256": input_hashes,
        "output_sha256": {name: sha(output / name) for name in files},
        "science_measurement": False, "psf_approved": False,
        "limitations": [
            "Only 18 isolated outer-field Gaia-supported stars on 18 HSC frames are selected; this is an exploratory pilot, not a representative PSF model.",
            "Gaia-to-coadd epoch propagation is unresolved; all selected Gaia variability flags require review.",
            "The local sky annuli and subtraction may remove real extended light or retain gradients; matched controls do not establish a PSF-only signal.",
            "The frame bootstrap describes resampling of this small selected sample and omits coadd covariance, selection uncertainty, and systematic calibration.",
            "No ICL surface-brightness correction or scientific measurement is authorized by this diagnostic.",
        ],
    }
    (output / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", type=Path, default=BASE)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--outer-radius-kpc", type=float, default=800.0)
    parser.add_argument("--bootstrap-draws", type=int, default=10000)
    parser.add_argument("--seed", type=int, default=20261020)
    parser.add_argument("--preflight-only", action="store_true", help="inspect metadata and selection only; no image pixels are opened")
    args = parser.parse_args()
    candidates, provenance = select_candidates(args.sample, outer_radius_kpc=args.outer_radius_kpc)
    if args.preflight_only:
        print(json.dumps({**provenance, "selected_candidate_ids": [row["candidate_id"] for row in candidates]}, indent=2))
        return
    if args.output is None:
        parser.error("--output is required unless --preflight-only is used")
    if args.bootstrap_draws < 100:
        parser.error("--bootstrap-draws must be at least 100")
    print(json.dumps(run(args.sample, args.output, bootstrap_draws=args.bootstrap_draws,
                         seed=args.seed, outer_radius_kpc=args.outer_radius_kpc), indent=2))


if __name__ == "__main__":
    main()
