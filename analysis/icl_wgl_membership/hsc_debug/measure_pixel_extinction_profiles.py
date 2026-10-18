#!/usr/bin/env python3
"""Re-measure selected 200--300 kpc ICL annuli with pixel-varying SFD correction.

Uses only already-downloaded HSC target/control images, saved masks and cached SFD
maps. No network access or source-image modification occurs. This remains an
exploratory calibration sensitivity, not a science-approved measurement.
"""
import argparse
import ast
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import sys
import time

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.cosmology import FlatLambdaCDM
from astropy.io import fits
from astropy.wcs import WCS
import astropy.units as u
from scipy.ndimage import map_coordinates

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from download_pilot import file_sha256 as sha256, validate_fits
from evaluate_prestack_pair_support import pair_block_labels
from prepare_ancillary import COEFFICIENTS
from prestack_core import FACTORS, geometry, measure
from validate_control_images import image_path, load_controls, nJy_factor


DEFAULT_DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
GRID_SIZE = 128
RING_INDEX = 3  # EDGES[3:5] = 200--300 kpc


def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def read_tsv(path):
    with Path(path).open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def write_json(path, value):
    Path(path).write_text(json.dumps(value, ensure_ascii=False, indent=2,
                                    sort_keys=True, allow_nan=False) + "\n",
                          encoding="utf-8")


def write_tsv(path, rows, fields):
    with Path(path).open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t",
                                extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def verify_map(path, expected_sha, expected_bytes):
    path = Path(path)
    if not path.is_file() or path.stat().st_size != int(expected_bytes):
        raise ValueError(f"SFD map missing or wrong size: {path}")
    digest = sha256(path)
    if digest != expected_sha:
        raise ValueError(f"SFD map SHA-256 mismatch: {path}")
    with fits.open(path, memmap=True, checksum=True) as hdul:
        if hdul[0].shape != (4096, 4096):
            raise ValueError(f"Unexpected SFD map shape: {path}")
        header = hdul[0].header.copy()
    return path, WCS(header), digest


def ebv_cutout(science_wcs, shape, maps, rng):
    """Approximate dense SFD cutout by oversampled WCS grid, then audit samples."""
    height, width = shape
    gx = np.linspace(0, width - 1, GRID_SIZE)
    gy = np.linspace(0, height - 1, GRID_SIZE)
    grid_x, grid_y = np.meshgrid(gx, gy)
    sky = science_wcs.pixel_to_world(grid_x, grid_y).galactic
    gal_l = sky.l.deg.ravel()
    gal_b = sky.b.deg.ravel()
    grid_values = np.empty(gal_l.shape, dtype=np.float64)
    for pole, selected in (("ngp", gal_b >= 0), ("sgp", gal_b < 0)):
        if not np.any(selected):
            continue
        map_path, map_wcs, map_data = maps[pole]
        px, py = map_wcs.all_world2pix(gal_l[selected], gal_b[selected], 0)
        if (np.any(~np.isfinite(px)) or np.any(~np.isfinite(py))
                or np.any(px < 0) or np.any(py < 0)
                or np.any(px > map_data.shape[1] - 1)
                or np.any(py > map_data.shape[0] - 1)):
            raise ValueError(f"SFD grid leaves map footprint: {map_path}")
        grid_values[selected] = map_coordinates(map_data, [py, px], order=1, prefilter=False)
    coarse = grid_values.reshape(GRID_SIZE, GRID_SIZE)
    if not np.all(np.isfinite(coarse)) or np.any(coarse < 0):
        raise ValueError("Nonfinite or negative SFD values in cutout grid")
    ycoord = np.linspace(0, GRID_SIZE - 1, height)
    xcoord = np.linspace(0, GRID_SIZE - 1, width)
    dense_y = np.broadcast_to(ycoord[:, None], (height, width))
    dense_x = np.broadcast_to(xcoord[None, :], (height, width))
    ebv = map_coordinates(coarse, [dense_y, dense_x],
                           order=1, prefilter=False)

    # Compare the coarse interpolation to direct SFD WCS interpolation at a
    # randomized sample of science-image pixels on every frame.
    ncheck = min(512, height * width)
    ys = rng.integers(0, height, size=ncheck)
    xs = rng.integers(0, width, size=ncheck)
    exact_sky = science_wcs.pixel_to_world(xs, ys).galactic
    exact = np.empty(ncheck, dtype=np.float64)
    for pole, selected in (("ngp", exact_sky.b.deg >= 0), ("sgp", exact_sky.b.deg < 0)):
        if not np.any(selected):
            continue
        _, map_wcs, map_data = maps[pole]
        px, py = map_wcs.all_world2pix(exact_sky.l.deg[selected], exact_sky.b.deg[selected], 0)
        if (np.any(~np.isfinite(px)) or np.any(~np.isfinite(py))
                or np.any(px < 0) or np.any(py < 0)
                or np.any(px > map_data.shape[1] - 1)
                or np.any(py > map_data.shape[0] - 1)):
            raise ValueError(f"SFD validation sample leaves map footprint: {maps[pole][0]}")
        exact[selected] = map_coordinates(map_data, [py, px], order=1, prefilter=False)
    approx = ebv[ys, xs]
    max_abs_error = float(np.max(np.abs(exact - approx)))
    if not np.all(np.isfinite(ebv)) or np.any(ebv < 0):
        raise ValueError("Invalid dense SFD cutout")
    return ebv, max_abs_error


def profile_pixel_corrected(image_path, frame, expected_image_sha,
                            mask, coefficient, maps, rng, h0, omega_m, expected_raw=None,
                            expected_coverage=None, center_ebv=None, control_center=None,
                            control_kpp=None):
    image_path = Path(image_path)
    if not image_path.is_file() or sha256(image_path) != expected_image_sha:
        raise ValueError(f"Image hash mismatch: {image_path}")
    info = validate_fits(image_path, frame)
    if info["coadd_header_filter"] not in COEFFICIENTS:
        raise ValueError(f"Unrecognized HSC-I filter: {info['coadd_header_filter']}")
    factor = nJy_factor(info)
    scale = info["pixel_scale_arcsec"]
    cosmo = FlatLambdaCDM(H0=float(h0), Om0=float(omega_m), Tcmb0=0)
    kpp = cosmo.kpc_proper_per_arcmin(float(frame["z"])).value * scale / 60.0
    center = info["center_pixel_xy"]
    if control_center is not None and not np.allclose(center, control_center, rtol=0, atol=1e-4):
        raise ValueError("Control profile center differs from the reviewed FITS center")
    if control_kpp is not None and not math.isclose(kpp, control_kpp, rel_tol=0, abs_tol=1e-6):
        raise ValueError("Control profile kpc/pixel differs from the assigned target redshift")
    with fits.open(image_path, memmap=False, checksum=True) as hdul:
        image = np.asarray(hdul[info["science_hdu"]].data, dtype=np.float32)
        variance = np.asarray(hdul[info["variance_hdu"]].data, dtype=np.float32)
        header = hdul[info["science_hdu"]].header.copy()
        if header.get("INHERIT") is True:
            header.extend(hdul[0].header, unique=True)
    if mask.shape != image.shape or variance.shape != image.shape:
        raise ValueError("Image, variance and saved mask shapes differ")
    science_wcs = WCS(header).celestial
    ebv, map_interpolation_error = ebv_cutout(science_wcs, image.shape, maps, rng)
    ext_factor = np.power(10.0, 0.4 * coefficient * ebv, dtype=np.float64)
    geom = geometry(image.shape, center, kpp)
    native_variance = variance * factor**2
    image_njy = image * factor
    raw_measurement = measure(image_njy, native_variance, mask, geom, center, kpp)
    if expected_raw is not None and not math.isclose(
            float(raw_measurement["raw"][RING_INDEX]), float(expected_raw),
            rel_tol=0, abs_tol=2e-6):
        raise ValueError(f"Raw profile does not reproduce saved result: {image_path}")
    if expected_coverage is not None and not math.isclose(
            float(raw_measurement["coverage"][RING_INDEX]), float(expected_coverage),
            rel_tol=0, abs_tol=2e-6):
        raise ValueError(f"Raw coverage does not reproduce saved result: {image_path}")
    in_ring = geom[3] == RING_INDEX
    valid = (in_ring & ~mask & np.isfinite(image) & np.isfinite(variance) &
             (variance > 0))
    if not np.any(valid):
        raise ValueError(f"No valid pixels in 200--300 kpc annulus: {image_path}")
    plane_corrected = float(raw_measurement["corrected"][RING_INDEX])
    plane = np.asarray(raw_measurement["plane"], dtype=float)
    if np.all(np.isfinite(plane)):
        plane_image = plane[0] + plane[1] * (geom[0] / 1000.0) + plane[2] * (geom[1] / 1000.0)
        pixel_ext_plane_corrected = float(np.mean((image_njy[valid] - plane_image[valid]) * ext_factor[valid]))
    else:
        pixel_ext_plane_corrected = None
    center_factor = None if center_ebv is None else 10 ** (0.4 * coefficient * float(center_ebv))
    return {
        "raw_njy_arcsec2": float(raw_measurement["raw"][RING_INDEX]),
        "pixel_corrected_njy_arcsec2": float(np.mean(image_njy[valid] * ext_factor[valid])),
        "plane_corrected_njy_arcsec2": plane_corrected if math.isfinite(plane_corrected) else None,
        "pixel_extinction_after_plane_njy_arcsec2": pixel_ext_plane_corrected,
        "center_sfd_corrected_raw_njy_arcsec2": float(raw_measurement["raw"][RING_INDEX] * center_factor) if center_factor is not None else None,
        "center_sfd_corrected_plane_njy_arcsec2": plane_corrected * center_factor if center_factor is not None and math.isfinite(plane_corrected) else None,
        "coverage": float(raw_measurement["coverage"][RING_INDEX]),
        "n_valid_pixels": int(valid.sum()),
        "mean_ebv_valid_pixels": float(np.mean(ebv[valid])) if np.any(valid) else None,
        "mean_pixel_extinction_mag": float(np.mean(coefficient * ebv[valid])) if np.any(valid) else None,
        "p05_pixel_extinction_mag": float(np.percentile(coefficient * ebv[valid], 5)) if np.any(valid) else None,
        "p95_pixel_extinction_mag": float(np.percentile(coefficient * ebv[valid], 95)) if np.any(valid) else None,
        "mean_flux_factor_valid_pixels": float(np.mean(ext_factor[valid])) if np.any(valid) else None,
        "map_interpolation_max_abs_ebv_error": map_interpolation_error,
        "native_to_njy_arcsec2": float(factor), "center_pixel_xy": [float(v) for v in center],
        "kpc_per_pixel": float(kpp), "shape": list(image.shape),
    }


def bootstrap_mean(values, labels, replicas, seed):
    blocks = np.unique(labels)
    by_block = {block: np.flatnonzero(labels == block) for block in blocks}
    rng = np.random.default_rng(seed)
    draws = np.empty(replicas, dtype=float)
    values = np.asarray(values, dtype=float)
    for i in range(replicas):
        selected = rng.choice(blocks, size=len(blocks), replace=True)
        index = np.concatenate([by_block[block] for block in selected])
        draws[i] = np.mean(values[index])
    return [float(v) for v in np.percentile(draws, [2.5, 97.5])]


def partial_pearson(x, y, z, richness):
    arrays = [np.asarray(value, dtype=float) for value in (x, y, z, richness)]
    if len(arrays[0]) < 5 or any(not np.all(np.isfinite(value)) for value in arrays):
        return np.nan
    design = np.column_stack((np.ones(len(arrays[0])), arrays[2], arrays[3]))
    if np.linalg.matrix_rank(design) < 3:
        return np.nan
    xr = arrays[0] - design @ np.linalg.lstsq(design, arrays[0], rcond=None)[0]
    yr = arrays[1] - design @ np.linalg.lstsq(design, arrays[1], rcond=None)[0]
    denominator = float(np.linalg.norm(xr) * np.linalg.norm(yr))
    return float(xr @ yr / denominator) if denominator > 0 else np.nan


def bootstrap_association(rows, labels, replicas, seed):
    blocks = np.unique(labels)
    per_block = {block: np.flatnonzero(labels == block) for block in blocks}
    raw_x = np.asarray([row["icl_delta_raw_njy_arcsec2"] for row in rows], dtype=float)
    center_x = np.asarray([row["icl_delta_center_sfd_njy_arcsec2"] for row in rows], dtype=float)
    pixel_x = np.asarray([row["icl_delta_pixel_sfd_njy_arcsec2"] for row in rows], dtype=float)
    plane_x = np.asarray([row["icl_delta_pixel_sfd_after_plane_njy_arcsec2"] for row in rows], dtype=float)
    satellite_y = np.asarray([row["satellite_delta_center_sfd_njy_arcsec2"] for row in rows], dtype=float)
    z = np.asarray([row["z"] for row in rows], dtype=float)
    richness = np.asarray([row["richness"] for row in rows], dtype=float)
    inputs = {"raw": raw_x, "center_sfd": center_x, "pixel_sfd_raw": pixel_x,
              "pixel_sfd_after_plane": plane_x}
    point = {key: partial_pearson(value, satellite_y, z, richness) for key, value in inputs.items()}
    draws = {key: np.full(replicas, np.nan) for key in inputs}
    rng = np.random.default_rng(seed)
    for i in range(replicas):
        selected = rng.choice(blocks, size=len(blocks), replace=True)
        index = np.concatenate([per_block[block] for block in selected])
        for key, value in inputs.items():
            draws[key][i] = partial_pearson(value[index], satellite_y[index], z[index], richness[index])
    result = {}
    for key in inputs:
        finite = draws[key][np.isfinite(draws[key])]
        result[key] = {
            "partial_pearson_r": float(point[key]) if math.isfinite(point[key]) else None,
            "block_bootstrap_ci95": [float(x) for x in np.percentile(finite, [2.5, 97.5])] if len(finite) else None,
            "valid_replicas": int(len(finite)),
        }
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=DEFAULT_DATA)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--replicas", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=20261020)
    parser.add_argument("--limit", type=int)
    parser.add_argument("--filter", choices=("i", "i2"))
    args = parser.parse_args()
    if not os.environ.get("SLURM_JOB_ID"):
        parser.error("Run this selected-image I/O inside a Slurm allocation")
    data, work = args.data.resolve(), args.work.resolve()
    if work.exists():
        parser.error(f"Refusing to overwrite existing output: {work}")

    sample_path = data / "sample_manifest.json"
    mask_dir = data / "prestack_cluster_masks_2026-10-16_01"
    mask_summary_path = mask_dir / "summary.json"
    mask_records_path = mask_dir / "records.tsv"
    control_dir = data / "prestack_extraction_2026-10-15_01"
    control_summary_path = control_dir / "extraction_summary.json"
    pair_path = data / "icl_annulus_sensitivity_2026-10-18_02/selected_pair_profiles.tsv"
    pair_summary_path = data / "icl_annulus_sensitivity_2026-10-18_02/summary.json"
    center_sens_dir = data / "icl_extinction_sensitivity_2026-10-19_03"
    center_sens_summary_path = center_sens_dir / "summary.json"
    center_sens_path = center_sens_dir / "annulus_extinction_corrected.tsv"
    association_dir = data / "icl_satellite_cut_sensitivity_2026-10-18_01"
    association_summary_path = association_dir / "summary.json"
    association_path = association_dir / "matched_cluster_values.tsv"
    extinction_dir = data / "ancillary_2026-10-13_01"
    candidates_path = extinction_dir / "extinction_candidates.tsv"
    prep_summary_path = extinction_dir / "preparation_status.json"
    control_manifest_path = data / "control_download_01/download_manifest.json"
    control_status_path = data / "control_download_01/download_status.json"
    required = (sample_path, mask_summary_path, mask_records_path, control_summary_path,
                pair_path, pair_summary_path, candidates_path, prep_summary_path,
                control_manifest_path, control_status_path, center_sens_summary_path,
                center_sens_path, association_summary_path, association_path)
    absent = [str(path) for path in required if not path.is_file()]
    if absent:
        parser.error("Missing inputs: " + ", ".join(absent))

    sample = read_json(sample_path)
    masks_summary = read_json(mask_summary_path)
    controls_summary = read_json(control_summary_path)
    pair_summary = read_json(pair_summary_path)
    center_sens_summary = read_json(center_sens_summary_path)
    association_summary = read_json(association_summary_path)
    ext_preparation = read_json(prep_summary_path)
    if len(sample.get("clusters", [])) != 800 or masks_summary.get("n_success") != 800:
        parser.error("Expected all 800 target frames and mask records")
    if masks_summary.get("status") != "complete" or masks_summary.get("mask_ready_for_science") is not False:
        parser.error("Unexpected mask-candidate status or science gate")
    if controls_summary.get("status") != "complete" or controls_summary.get("n_measured") != 798:
        parser.error("Expected the complete exploratory control profiles")
    if pair_summary.get("status") != "complete" or pair_summary.get("science_measurement") is not False:
        parser.error("Unexpected annulus profile status or science gate")
    if pair_summary.get("output_sha256", {}).get(str(pair_path)) != sha256(pair_path):
        parser.error("Selected-pair input checksum mismatch")
    if center_sens_summary.get("science_measurement") is not False or center_sens_summary.get("status") != "complete_exploratory":
        parser.error("Unexpected center-extinction sensitivity status")
    if center_sens_summary.get("output_sha256", {}).get(str(center_sens_path)) != sha256(center_sens_path):
        parser.error("Center-extinction sensitivity table hash mismatch")
    if association_summary.get("science_measurement") is not False or association_summary.get("status") != "complete":
        parser.error("Unexpected satellite association status")
    if association_summary.get("output_sha256", {}).get(str(association_path)) != sha256(association_path):
        parser.error("Satellite association input hash mismatch")
    if ext_preparation.get("status") != "extinction_candidates_ready_controls_await_catalogues":
        parser.error("Unexpected SFD map provenance status")

    map_provenance = ext_preparation.get("maps", {})
    maps = {}
    for pole in ("ngp", "sgp"):
        item = map_provenance.get(pole)
        if not item:
            parser.error(f"Missing SFD map provenance for {pole}")
        path, map_wcs, digest = verify_map(item["file"], item["sha256"], item["bytes"])
        with fits.open(path, memmap=True, checksum=True) as hdul:
            map_data = hdul[0].data
        maps[pole] = (path, map_wcs, map_data)

    sample_rows = {str(row["id"]): row for row in sample["clusters"]}
    mask_rows = {str(row["id"]): row for row in read_tsv(mask_records_path)}
    control_rows_by_id = {str(row["control_id"]): row
                          for row in controls_summary.get("records", [])}
    if len(sample_rows) != 800 or len(mask_rows) != 800 or len(control_rows_by_id) != 798:
        parser.error("Duplicate or missing input manifest IDs")
    candidates = {str(row["cluster_id"]): row for row in read_tsv(candidates_path)}
    if len(candidates) != 800:
        parser.error("Expected 800 SFD center candidates")
    annulus = [row for row in read_tsv(pair_path)
               if float(row["kron_factor"]) == 2.5
               and float(row["coverage_floor"]) == 0.40
               and row["radial_bin_kpc"] == "200-300"]
    if len(annulus) != 224 or len({(row["id"], row["filter"]) for row in annulus}) != 224:
        parser.error("Unexpected fixed-cut annulus sample")
    if args.limit is not None and args.limit < 1:
        parser.error("--limit must be positive")
    requested_source = annulus if args.filter is None else [row for row in annulus if row["filter"] == args.filter]
    requested = requested_source if args.limit is None else requested_source[:args.limit]
    full_run = args.limit is None and args.filter is None
    center_rows = read_tsv(center_sens_path)
    center_by_key = {(str(row["cluster_id"]), row["filter"]): row for row in center_rows}
    if len(center_by_key) != 224:
        parser.error("Unexpected center-extinction sensitivity table size")

    control_doc, control_download, control_base, control_source_hashes = load_controls(data)
    control_frames = {str(row["id"]): row for row in control_doc["clusters"]}
    rng = np.random.default_rng(args.seed)
    work.mkdir(parents=True, exist_ok=False)
    start = time.monotonic()
    result_rows = []
    max_map_error = 0.0
    source_hashes = {str(path): sha256(path) for path in required}
    source_hashes.update(control_source_hashes)
    source_hashes.update({str(path): sha256(path) for path in
                          (Path(map_provenance[pole]["file"]) for pole in ("ngp", "sgp"))})

    for index, pair in enumerate(requested, 1):
        ident, band, control_id = str(pair["id"]), pair["filter"], str(pair["control_id"])
        target = sample_rows[ident]
        target_mask_record = mask_rows[ident]
        target_path = Path(target_mask_record["image"])
        target_mask_path = Path(target_mask_record["mask_file"])
        if mask_dir.resolve() not in target_mask_path.resolve().parents:
            raise ValueError(f"Target mask escapes prepared mask directory: {target_mask_path}")
        expected_mask_sha = masks_summary.get("mask_output_sha256", {}).get(str(target_mask_path))
        if expected_mask_sha != target_mask_record.get("mask_sha256") or sha256(target_mask_path) != expected_mask_sha:
            raise ValueError(f"Target mask checksum mismatch for {ident}")
        with fits.open(target_mask_path, memmap=False, checksum=True) as hdul:
            target_mask = np.asarray(hdul["MASK_K2P5"].data, dtype=bool)
        target_product = profile_pixel_corrected(
            target_path, target, target_mask_record["image_sha256"],
            target_mask, COEFFICIENTS[band], maps, rng,
            sample["h0"], sample["omega_m"],
            expected_raw=float(pair["target_flux"]),
            expected_coverage=float(pair["target_coverage"]),
            center_ebv=float(candidates[ident]["ebv_sfd"]))
        target_product["image_sha256"] = target_mask_record["image_sha256"]
        target_product["mask_sha256"] = expected_mask_sha
        source_hashes[str(target_path)] = target_mask_record["image_sha256"]
        source_hashes[str(target_mask_path)] = expected_mask_sha

        control_profile = control_rows_by_id.get(control_id)
        control_frame = control_frames.get(control_id)
        if (control_profile is None or control_frame is None
                or str(control_profile["target_cluster_id"]) != ident
                or str(control_frame["target_cluster_id"]) != ident
                or control_profile["header_filter"] != band):
            raise ValueError(f"Target/control/filter provenance mismatch for {ident}")
        control_baseline = control_dir / control_id / "baseline.npz"
        if sha256(control_baseline) != control_profile["baseline_sha256"]:
            raise ValueError(f"Control profile checksum mismatch for {control_id}")
        with np.load(control_baseline) as archive:
            shape = tuple(int(v) for v in archive["image_shape"])
            mask_bits = archive["mask_bits"][0]
            control_mask = np.unpackbits(mask_bits, axis=-1, count=shape[1],
                                         bitorder="little").astype(bool)
            control_mask = control_mask.reshape(shape)
            center = archive["center"].astype(float).tolist()
            kpp = float(archive["kpp"])
            control_raw_profile = float(archive["raw"][0, RING_INDEX])
            control_coverage_profile = float(archive["coverage"][0, RING_INDEX])
        control_img_path = Path(control_profile["source_image"])
        control_product = profile_pixel_corrected(
            control_img_path, control_frame, control_profile["source_sha256"],
            control_mask, COEFFICIENTS[band], maps, rng,
            sample["h0"], sample["omega_m"],
            expected_raw=float(pair["control_flux"]),
            expected_coverage=float(pair["control_coverage"]),
            center_ebv=float(control_frame["ebv_sfd"]),
            control_center=center, control_kpp=kpp)
        if not math.isclose(control_raw_profile, float(pair["control_flux"]), rel_tol=0, abs_tol=2e-6):
            raise ValueError(f"Control saved profile does not match pair table for {control_id}")
        if not math.isclose(control_coverage_profile, float(pair["control_coverage"]), rel_tol=0, abs_tol=2e-6):
            raise ValueError(f"Control saved coverage does not match pair table for {control_id}")
        source_hashes[str(control_img_path)] = control_profile["source_sha256"]
        source_hashes[str(control_baseline)] = control_profile["baseline_sha256"]
        max_map_error = max(max_map_error, target_product["map_interpolation_max_abs_ebv_error"],
                            control_product["map_interpolation_max_abs_ebv_error"])
        center_row = center_by_key.get((ident, band))
        if center_row is None:
            raise ValueError(f"Missing center-corrected comparison row for {(ident, band)}")
        if not math.isclose(float(pair["delta_flux"]), float(center_row["icl_delta_raw_njy_arcsec2"]),
                            rel_tol=0, abs_tol=1e-7):
            raise ValueError(f"Center-corrected table raw profile mismatch for {(ident, band)}")
        result_rows.append({
            "cluster_id": ident, "filter": band, "control_id": control_id,
            "target_raw_njy_arcsec2": target_product["raw_njy_arcsec2"],
            "control_raw_njy_arcsec2": control_product["raw_njy_arcsec2"],
            "raw_delta_njy_arcsec2": target_product["raw_njy_arcsec2"] - control_product["raw_njy_arcsec2"],
            "target_pixel_ext_corrected_njy_arcsec2": target_product["pixel_corrected_njy_arcsec2"],
            "control_pixel_ext_corrected_njy_arcsec2": control_product["pixel_corrected_njy_arcsec2"],
            "pixel_ext_corrected_delta_njy_arcsec2": target_product["pixel_corrected_njy_arcsec2"] - control_product["pixel_corrected_njy_arcsec2"],
            "target_pixel_ext_after_plane_njy_arcsec2": target_product["pixel_extinction_after_plane_njy_arcsec2"],
            "control_pixel_ext_after_plane_njy_arcsec2": control_product["pixel_extinction_after_plane_njy_arcsec2"],
            "pixel_ext_after_plane_delta_njy_arcsec2": (
                target_product["pixel_extinction_after_plane_njy_arcsec2"] -
                control_product["pixel_extinction_after_plane_njy_arcsec2"]
                if target_product["pixel_extinction_after_plane_njy_arcsec2"] is not None and
                control_product["pixel_extinction_after_plane_njy_arcsec2"] is not None else None),
            "center_sfd_raw_delta_njy_arcsec2": target_product["center_sfd_corrected_raw_njy_arcsec2"] - control_product["center_sfd_corrected_raw_njy_arcsec2"],
            "target_center_sfd_plane_njy_arcsec2": target_product["center_sfd_corrected_plane_njy_arcsec2"],
            "control_center_sfd_plane_njy_arcsec2": control_product["center_sfd_corrected_plane_njy_arcsec2"],
            "center_sfd_plane_delta_njy_arcsec2": (
                target_product["center_sfd_corrected_plane_njy_arcsec2"] -
                control_product["center_sfd_corrected_plane_njy_arcsec2"]
                if target_product["center_sfd_corrected_plane_njy_arcsec2"] is not None and
                control_product["center_sfd_corrected_plane_njy_arcsec2"] is not None else None),
            "center_sfd_corrected_delta_njy_arcsec2": float(center_row["icl_delta_sfd_corrected_njy_arcsec2"]),
            "target_coverage": target_product["coverage"], "control_coverage": control_product["coverage"],
            "common_coverage": min(target_product["coverage"], control_product["coverage"]),
            "target_n_valid_pixels": target_product["n_valid_pixels"],
            "control_n_valid_pixels": control_product["n_valid_pixels"],
            "target_mean_ebv_in_annulus": target_product["mean_ebv_valid_pixels"],
            "control_mean_ebv_in_annulus": control_product["mean_ebv_valid_pixels"],
            "target_mean_pixel_extinction_mag": target_product["mean_pixel_extinction_mag"],
            "control_mean_pixel_extinction_mag": control_product["mean_pixel_extinction_mag"],
            "target_pixel_extinction_p05_mag": target_product["p05_pixel_extinction_mag"],
            "target_pixel_extinction_p95_mag": target_product["p95_pixel_extinction_mag"],
            "control_pixel_extinction_p05_mag": control_product["p05_pixel_extinction_mag"],
            "control_pixel_extinction_p95_mag": control_product["p95_pixel_extinction_mag"],
            "target_mean_pixel_flux_factor": target_product["mean_flux_factor_valid_pixels"],
            "control_mean_pixel_flux_factor": control_product["mean_flux_factor_valid_pixels"],
            "target_image_sha256": target_product["image_sha256"],
            "target_mask_sha256": target_product["mask_sha256"],
            "control_image_sha256": control_profile["source_sha256"],
            "control_baseline_sha256": control_profile["baseline_sha256"],
            "target_ra_deg": float(pair["target_ra_deg"]),
            "target_dec_deg": float(pair["target_dec_deg"]),
            "control_ra_deg": float(pair["control_ra_deg"]),
            "control_dec_deg": float(pair["control_dec_deg"]),
        })
        if index % 25 == 0 or index == len(requested):
            print(f"Pixel-corrected {index}/{len(requested)} pairs", flush=True)

    result_path = work / "pixel_extinction_annulus_profiles.tsv"
    write_tsv(result_path, result_rows, list(result_rows[0]))
    summary = {
        "schema": "hsc-pixel-extinction-annulus-v1",
        "status": "complete_exploratory" if full_run and len(requested) == 224 else "pilot_incomplete",
        "science_measurement": False, "n_pairs": len(result_rows),
        "filter_counts": {band: sum(row["filter"] == band for row in result_rows) for band in ("i", "i2")},
        "radial_annulus_kpc": [200, 300], "kron_factor": 2.5, "coverage_floor": 0.40,
        "extinction": {"map": "SFD dust maps, bilinear interpolation in native map WCS",
                       "application": "pixel-level factor 10^(0.4 * A_lambda), applied separately to each unmasked target/control pixel",
                       "coefficients": COEFFICIENTS, "coarse_map_grid": [GRID_SIZE, GRID_SIZE],
                       "dense_interpolation_validation_pixels_per_frame": 512,
                       "max_abs_ebv_interpolation_error": max_map_error,
                       "map_provenance": map_provenance},
        "selection": "fixed existing K2.5 / common target-control coverage >=40% / 200-300 kpc sample; selection not recomputed after correction",
        "bootstrap": {"replicas": args.replicas, "seed": args.seed,
                      "blocks": "connected joint 5-degree target/control blocks"},
        "limitations": [
            "SFD map is sampled at native resolution; coarse-grid interpolation has a recorded direct-WCS error audit but does not increase the map's intrinsic angular resolution.",
            "Existing source masks are held fixed; the analysis does not rerun source detection after pixel correction.",
            "The existing selected-pair table fixes the sample before correction; this is a sensitivity audit, not a new science selection.",
            "Upstream sky subtraction, large-scale PSF wings, mask leakage and the full ICL/satellite measurement covariance remain open.",
        ],
        "input_sha256": source_hashes, "program_sha256": sha256(__file__),
        "output_sha256": {str(result_path): sha256(result_path)},
        "elapsed_seconds": time.monotonic() - start,
    }
    rows_by_band = {}
    for band in ("i", "i2"):
        selected = [row for row in result_rows if row["filter"] == band]
        if not selected:
            continue
        target_coords = [{"ra_deg": row["target_ra_deg"], "dec_deg": row["target_dec_deg"]} for row in selected]
        control_coords = [{"ra_deg": row["control_ra_deg"], "dec_deg": row["control_dec_deg"]} for row in selected]
        _, labels, nt, nc, nblocks, largest = pair_block_labels(
            target_coords, control_coords, np.ones(len(selected), dtype=bool))
        raw = np.asarray([row["raw_delta_njy_arcsec2"] for row in selected])
        pixel = np.asarray([row["pixel_ext_corrected_delta_njy_arcsec2"] for row in selected])
        plane_values = np.asarray([row["pixel_ext_after_plane_delta_njy_arcsec2"]
                                   if row["pixel_ext_after_plane_delta_njy_arcsec2"] is not None
                                   else np.nan for row in selected], dtype=float)
        center_plane_values = np.asarray([row["center_sfd_plane_delta_njy_arcsec2"]
                                          if row["center_sfd_plane_delta_njy_arcsec2"] is not None
                                          else np.nan for row in selected], dtype=float)
        rows_by_band[band] = {
            "n_pairs": len(selected), "n_joint_5deg_blocks": nblocks,
            "n_target_cells": nt, "n_control_cells": nc, "largest_joint_block": largest,
            "mean_raw_delta_njy_arcsec2": float(np.mean(raw)),
            "mean_pixel_corrected_delta_njy_arcsec2": float(np.mean(pixel)),
            "mean_pixel_corrected_shift_njy_arcsec2": float(np.mean(pixel - raw)),
            "mean_center_sfd_delta_njy_arcsec2": float(np.mean([
                row["center_sfd_corrected_delta_njy_arcsec2"] for row in selected])),
            "mean_pixel_minus_center_sfd_raw_delta_njy_arcsec2": float(np.mean(
                pixel - np.asarray([row["center_sfd_corrected_delta_njy_arcsec2"] for row in selected]))),
            "n_plane_first_pixel_corrected_pairs": int(np.isfinite(plane_values).sum()),
            "mean_plane_first_pixel_corrected_delta_njy_arcsec2": (
                float(np.nanmean(plane_values)) if np.isfinite(plane_values).any() else None),
            "mean_center_sfd_plane_corrected_delta_njy_arcsec2": (
                float(np.nanmean(center_plane_values)) if np.isfinite(center_plane_values).any() else None),
            "raw_block_bootstrap_ci95": bootstrap_mean(raw, labels, args.replicas, args.seed),
            "pixel_corrected_block_bootstrap_ci95": bootstrap_mean(pixel, labels, args.replicas, args.seed + 1),
            "plane_first_pixel_corrected_block_bootstrap_ci95": (
                bootstrap_mean(plane_values, labels, args.replicas, args.seed + 2)
                if np.all(np.isfinite(plane_values)) else None),
        }
    summary["profile_means_by_filter"] = rows_by_band
    if full_run and len(requested) == 224:
        pixel_by_key = {(row["cluster_id"], row["filter"]): row for row in result_rows}
        association_rows = [row for row in read_tsv(association_path)
                            if float(row["kron_factor"]) == 2.5
                            and float(row["coverage_floor"]) == 0.40]
        association_output = []
        for row in association_rows:
            key = (str(row["cluster_id"]), row["filter"])
            pixel_row = pixel_by_key.get(key)
            if pixel_row is None:
                raise ValueError(f"Association row lacks pixel-corrected ICL profile: {key}")
            if not math.isclose(float(row["icl_delta_200_300"]), pixel_row["raw_delta_njy_arcsec2"],
                                rel_tol=0, abs_tol=1e-6):
                raise ValueError(f"Association raw ICL value does not reproduce for {key}")
            a_center = float(candidates[key[0]]["ebv_sfd"]) * COEFFICIENTS[key[1]]
            satellite_raw = float(row["satellite_delta_50_75"])
            satellite_center = satellite_raw * 10 ** (0.4 * a_center)
            association_output.append({
                "cluster_id": key[0], "filter": key[1],
                "icl_delta_raw_njy_arcsec2": pixel_row["raw_delta_njy_arcsec2"],
                "icl_delta_center_sfd_njy_arcsec2": pixel_row["center_sfd_corrected_delta_njy_arcsec2"],
                "icl_delta_pixel_sfd_njy_arcsec2": pixel_row["pixel_ext_corrected_delta_njy_arcsec2"],
                "icl_delta_pixel_sfd_after_plane_njy_arcsec2": pixel_row["pixel_ext_after_plane_delta_njy_arcsec2"],
                "satellite_delta_raw_njy_arcsec2": satellite_raw,
                "satellite_delta_center_sfd_njy_arcsec2": satellite_center,
                "z": float(row["z"]), "richness": float(row["richness"]),
                "target_ra_deg": float(row["target_ra_deg"]), "target_dec_deg": float(row["target_dec_deg"]),
                "control_ra_deg": float(row["control_ra_deg"]), "control_dec_deg": float(row["control_dec_deg"]),
            })
        association_path_out = work / "pixel_extinction_satellite_association.tsv"
        write_tsv(association_path_out, association_output, list(association_output[0]))
        assoc_metrics = {}
        for index, band in enumerate(("i", "i2")):
            selected = [row for row in association_output if row["filter"] == band]
            target_coords = [{"ra_deg": row["target_ra_deg"], "dec_deg": row["target_dec_deg"]} for row in selected]
            control_coords = [{"ra_deg": row["control_ra_deg"], "dec_deg": row["control_dec_deg"]} for row in selected]
            _, labels, _, _, nblocks, largest = pair_block_labels(
                target_coords, control_coords, np.ones(len(selected), dtype=bool))
            assoc_metrics[band] = {"n_pairs": len(selected), "n_joint_5deg_blocks": nblocks,
                                   "largest_joint_block": largest,
                                   "partial_pearson_vs_center_corrected_satellite": bootstrap_association(
                                       selected, labels, args.replicas, args.seed + 200 + index)}
        summary["pixel_corrected_icl_satellite_association"] = assoc_metrics
        summary["output_sha256"][str(association_path_out)] = sha256(association_path_out)
    summary_path = work / "summary.json"
    write_json(summary_path, summary)
    print(json.dumps({"status": summary["status"], "science_measurement": False,
                      "n_pairs": summary["n_pairs"], "filters": summary["filter_counts"],
                      "max_map_interpolation_error": max_map_error,
                      "output": str(work)}, ensure_ascii=False), flush=True)


if __name__ == "__main__":
    main()
