#!/usr/bin/env python3
"""Audit saturation, nearby detections, and source-mask support for Gaia PSF stamps.

Uses existing HSC detection tables, mask HDUs, and saved star-stamp masks. It
does not read science-image arrays, request data, select a final PSF sample, or
approve an empirical PSF.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import time

import numpy as np
from astropy.io import fits


DEFAULT_SAMPLE = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
DEFAULT_LEDGER = DEFAULT_SAMPLE / "empirical_psf_gaia_full_sample_2026-10-20_02/star_ledger.tsv"
RINGS = {"core_0_12": (0.0, 12.0), "wing_48_64": (48.0, 64.0), "wing_64_96": (64.0, 96.0)}
NEIGHBOR_RADII = (12.0, 24.0, 48.0, 64.0, 96.0)
SAT_BIT = 1 << 1
BRIGHT_OBJECT_BIT = 1 << 9


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with Path(path).open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def precompute_rings(shape: tuple[int, int], x: float, y: float) -> dict[str, dict]:
    """Build annulus and sector geometry once for reuse across masks."""
    if len(shape) != 2 or not np.isfinite([x, y]).all():
        raise ValueError("invalid stamp geometry")
    yy, xx = np.indices(shape)
    dx, dy = xx - x, yy - y
    radius = np.hypot(dx, dy)
    angle = (np.arctan2(dy, dx) + np.pi) % (2.0 * np.pi)
    sector = np.minimum(7, np.floor(angle * (8.0 / (2.0 * np.pi))).astype(int))
    result = {"core_ring": radius < 12.0}
    for name, (r_min, r_max) in RINGS.items():
        ring = (radius >= r_min) & (radius < r_max)
        if not np.any(ring):
            raise ValueError(f"ring has no pixels in stamp: {name}")
        result[name] = {"ring": ring, "sector": sector}
    return result


def ring_support(mask: np.ndarray, x: float, y: float, r_min: float, r_max: float) -> dict:
    """Summarize unmasked pixel fraction and azimuthal support in a stamp ring."""
    if mask.ndim != 2 or not np.isfinite([x, y, r_min, r_max]).all() or r_min < 0 or r_max <= r_min:
        raise ValueError("invalid mask or ring geometry")
    yy, xx = np.indices(mask.shape)
    dx, dy = xx - x, yy - y
    radius = np.hypot(dx, dy)
    ring = (radius >= r_min) & (radius < r_max)
    if not np.any(ring):
        raise ValueError("ring has no pixels in stamp")
    angle = (np.arctan2(dy, dx) + np.pi) % (2.0 * np.pi)
    sector = np.minimum(7, np.floor(angle * (8.0 / (2.0 * np.pi))).astype(int))
    return _ring_support_from_geometry(mask, ring, sector)


def _ring_support_from_geometry(mask: np.ndarray, ring: np.ndarray, sector: np.ndarray) -> dict:
    available = ring & ~mask.astype(bool)
    if not np.any(ring):
        raise ValueError("ring has no pixels in stamp")
    sector_fractions = []
    for index in range(8):
        total = int(np.count_nonzero(ring & (sector == index)))
        valid = int(np.count_nonzero(available & (sector == index)))
        sector_fractions.append(valid / total if total else 0.0)
    return {
        "ring_pixels": int(np.count_nonzero(ring)),
        "available_pixels": int(np.count_nonzero(available)),
        "available_fraction": float(np.count_nonzero(available) / np.count_nonzero(ring)),
        "sectors_at_least_25pct_open": int(sum(value >= 0.25 for value in sector_fractions)),
    }


def _flagged_pixels(flag: np.ndarray, ring: np.ndarray) -> int:
    return int(np.count_nonzero(flag & ring))


def neighbor_summary(detections: list[dict[str, str]], number: int, x: float, y: float,
                     target_flux: float) -> dict:
    candidates = []
    for row in detections:
        if int(row["NUMBER"]) == number:
            continue
        nx, ny = float(row["X_IMAGE"]) - 1.0, float(row["Y_IMAGE"]) - 1.0
        distance = math.hypot(nx - x, ny - y)
        if distance < max(NEIGHBOR_RADII):
            try:
                flux = float(row["FLUX_ISO"])
            except (KeyError, TypeError, ValueError):
                flux = float("nan")
            candidates.append((distance, flux))
    distances = [value[0] for value in candidates]
    positive_flux = [flux / target_flux for _distance, flux in candidates
                     if math.isfinite(flux) and flux > 0 and math.isfinite(target_flux) and target_flux > 0]
    result = {
        "nearest_neighbor_distance_pix": min(distances) if distances else None,
        "max_neighbor_to_target_flux_iso_ratio_within_96pix": max(positive_flux) if positive_flux else None,
    }
    for radius in NEIGHBOR_RADII:
        result[f"n_neighbors_within_{int(radius)}pix"] = sum(distance < radius for distance in distances)
    return result


def _load_detection_table(path: Path) -> list[dict[str, str]]:
    rows = read_tsv(path)
    if not rows or not {"NUMBER", "X_IMAGE", "Y_IMAGE", "FLUX_ISO", "FLAGS"}.issubset(rows[0]):
        raise ValueError(f"unexpected detection-table schema: {path}")
    if [int(row["NUMBER"]) for row in rows] != list(range(1, len(rows) + 1)):
        raise ValueError(f"detection numbering is not contiguous: {path}")
    return rows


def run(sample_root: Path, ledger_path: Path, output: Path) -> dict:
    sample_root, ledger_path, output = map(Path, (sample_root, ledger_path, output))
    if output.exists():
        raise FileExistsError(f"refusing to overwrite stamp-quality audit: {output}")
    ledger = read_tsv(ledger_path)
    stamped = [row for row in ledger if row.get("stamp", "").strip()]
    if len(stamped) != 333 or len({row["candidate_id"] for row in stamped}) != len(stamped):
        raise ValueError("expected 333 unique hash-verified Gaia star stamps")
    ledger_summary_path = ledger_path.parent / "summary.json"
    ledger_summary = json.loads(ledger_summary_path.read_text(encoding="utf-8"))
    detection_dir = Path(ledger_summary["detection_dir"])
    download_status_path = sample_root / "download_01/download_status.json"
    download_status = json.loads(download_status_path.read_text(encoding="utf-8"))
    download_by_id = {str(row["id"]): row for row in download_status["completed"]}
    grouped: dict[str, list[dict[str, str]]] = {}
    for row in stamped:
        grouped.setdefault(row["cluster_id"], []).append(row)

    results = []
    source_hashes = {str(ledger_path.resolve()): sha256(ledger_path),
                     str(ledger_summary_path.resolve()): sha256(ledger_summary_path),
                     str(download_status_path.resolve()): sha256(download_status_path)}
    mask_section_hashes = {}
    for cluster_id, stars in sorted(grouped.items(), key=lambda item: int(item[0])):
        frame = download_by_id.get(cluster_id)
        if frame is None:
            raise ValueError(f"Gaia star frame missing from verified download status: {cluster_id}")
        image_path = sample_root / "download_01/images" / frame["file"]
        if not image_path.is_file() or image_path.stat().st_size != int(frame["bytes"]):
            raise ValueError(f"HSC image is missing or size changed: {image_path}")
        det_path = detection_dir / f"cluster_{cluster_id}" / "detections.tsv"
        for path in (det_path,):
            source_hashes[str(path.resolve())] = sha256(path)
        detections = _load_detection_table(det_path)
        review = frame["review"]
        mask_hdu = int(review["mask_hdu"])
        section_digest = hashlib.sha256()
        with fits.open(image_path, mode="readonly", memmap=True, lazy_load_hdus=True) as image_hdus:
            raw_mask = image_hdus[mask_hdu]
            for star in stars:
                stamp_path = Path(star["stamp"])
                stamp_digest = sha256(stamp_path)
                if stamp_digest != star["stamp_sha256"]:
                    raise ValueError(f"star stamp checksum mismatch: {stamp_path}")
                source_hashes[str(stamp_path.resolve())] = stamp_digest
                number = int(star["number"])
                if not 1 <= number <= len(detections):
                    raise ValueError(f"detection number is out of range: {star['candidate_id']}")
                detection = detections[number - 1]
                if int(detection["NUMBER"]) != number:
                    raise ValueError("Gaia candidate detection identity changed")
                x, y = float(detection["X_IMAGE"]) - 1.0, float(detection["Y_IMAGE"]) - 1.0
                if not np.allclose([x, y], [float(star["native_x"]), float(star["native_y"])], rtol=0, atol=1e-9):
                    raise ValueError(f"HSC detection center differs from Gaia ledger: {star['candidate_id']}")
                x0, y0 = int(star["x0"]), int(star["y0"])
                x1, y1 = x0 + 401, y0 + 401
                if min(x0, y0) < 0 or x1 > raw_mask.shape[1] or y1 > raw_mask.shape[0]:
                    raise ValueError(f"stamp mask footprint is outside HSC image: {star['candidate_id']}")
                raw = np.asarray(raw_mask.section[y0:y1, x0:x1])
                if raw.shape != (401, 401):
                    raise ValueError(f"HSC mask section has wrong shape: {star['candidate_id']}")
                section_digest.update(
                    f"{star['candidate_id']}:{x0}:{y0}:{x1}:{y1}\n".encode("ascii"))
                section_digest.update(np.ascontiguousarray(raw).tobytes())
                relative_x, relative_y = x - x0, y - y0
                stamp_path = Path(star["stamp"])
                with fits.open(stamp_path, mode="readonly", memmap=False) as star_stamp:
                    if not all(name in star_stamp for name in ("MASK_K2P5", "MASK_K3P5")):
                        raise ValueError(f"expected source masks are missing in {stamp_path}")
                    header = star_stamp[0].header
                    if (int(header.get("CANDID", -1)) != int(star["candidate_id"])
                            or int(header.get("X0", -1)) != x0
                            or int(header.get("Y0", -1)) != y0
                            or not np.allclose(
                                [float(header.get("XCENTER", float("nan"))),
                                 float(header.get("YCENTER", float("nan")))],
                                [relative_x, relative_y], rtol=0, atol=1e-9)):
                        raise ValueError(f"star stamp header geometry mismatch: {stamp_path}")
                    stamp_masks = {
                        name: np.asarray(star_stamp[name].data, dtype=bool)
                        for name in ("MASK_K2P5", "MASK_K3P5")
                    }
                ring_geometry = precompute_rings((401, 401), relative_x, relative_y)
                core_ring = ring_geometry["core_ring"]
                sat_core_pixels = int(np.count_nonzero(((raw & SAT_BIT) != 0) & core_ring))
                bright_core_pixels = int(np.count_nonzero(((raw & BRIGHT_OBJECT_BIT) != 0) & core_ring))
                result = {
                    "candidate_id": star["candidate_id"], "cluster_id": cluster_id,
                    "number": number, "gaia_source_id": star["source_id"],
                    "center_source": "HSC_DETECTION_X_IMAGE_Y_IMAGE",
                    "detection_flags": int(detection["FLAGS"]),
                    "hsc_sat_pixels_core_0_12": sat_core_pixels,
                    "hsc_bright_object_pixels_core_0_12": bright_core_pixels,
                }
                core_sector = ring_geometry["wing_48_64"]["sector"]
                for mask_name in ("MASK_K2P5", "MASK_K3P5"):
                    core_stats = _ring_support_from_geometry(
                        stamp_masks[mask_name], core_ring, core_sector,
                    )
                    prefix = mask_name.lower()
                    result[f"{prefix}_core_0_12_open_fraction"] = core_stats["available_fraction"]
                    result[f"{prefix}_core_0_12_open_sectors"] = core_stats["sectors_at_least_25pct_open"]
                for ring_name, (r0, r1) in RINGS.items():
                    ring = ring_geometry[ring_name]["ring"]
                    sector = ring_geometry[ring_name]["sector"]
                    result[f"hsc_sat_pixels_{ring_name}"] = _flagged_pixels((raw & SAT_BIT) != 0, ring)
                    result[f"hsc_bright_object_pixels_{ring_name}"] = _flagged_pixels((raw & BRIGHT_OBJECT_BIT) != 0, ring)
                    for mask_name in ("MASK_K2P5", "MASK_K3P5"):
                        stats = _ring_support_from_geometry(stamp_masks[mask_name], ring, sector)
                        prefix = f"{mask_name.lower()}_{ring_name}"
                        result[f"{prefix}_open_fraction"] = stats["available_fraction"]
                        result[f"{prefix}_open_sectors"] = stats["sectors_at_least_25pct_open"]
                target_flux = float(detection["FLUX_ISO"])
                result.update(neighbor_summary(detections, number, x, y, target_flux))
                results.append(result)
        mask_section_hashes[str(image_path.resolve())] = section_digest.hexdigest()

    output.mkdir(parents=True, exist_ok=False)
    table_path = output / "star_stamp_quality.tsv"
    fields = list(results[0])
    with table_path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(results)
    summary = {
        "schema": "astrafek-empirical-psf-star-stamp-quality-v1",
        "status": "complete_mask_and_catalog_quality_diagnostic",
        "n_stamps": len(results), "n_target_frames": len(grouped),
        "n_with_saturation_in_core": sum(row["hsc_sat_pixels_core_0_12"] > 0 for row in results),
        "n_with_bright_object_bit_in_core": sum(row["hsc_bright_object_pixels_core_0_12"] > 0 for row in results),
        "n_k2p5_core_open_below_99pct": sum(row["mask_k2p5_core_0_12_open_fraction"] < 0.99 for row in results),
        "n_k3p5_core_open_below_99pct": sum(row["mask_k3p5_core_0_12_open_fraction"] < 0.99 for row in results),
        "median_k2p5_core_open_fraction": float(np.median([
            row["mask_k2p5_core_0_12_open_fraction"] for row in results])),
        "median_k3p5_core_open_fraction": float(np.median([
            row["mask_k3p5_core_0_12_open_fraction"] for row in results])),
        "n_with_neighbor_within_12pix": sum(row["n_neighbors_within_12pix"] > 0 for row in results),
        "n_with_neighbor_within_48pix": sum(row["n_neighbors_within_48pix"] > 0 for row in results),
        "n_with_neighbor_within_96pix": sum(row["n_neighbors_within_96pix"] > 0 for row in results),
        "median_nearest_neighbor_distance_pix": float(np.median([
            row["nearest_neighbor_distance_pix"] for row in results
            if row["nearest_neighbor_distance_pix"] is not None])) if any(
                row["nearest_neighbor_distance_pix"] is not None for row in results) else None,
        "median_k2p5_open_fraction_48_96pix": float(np.median([
            min(row["mask_k2p5_wing_48_64_open_fraction"], row["mask_k2p5_wing_64_96_open_fraction"])
            for row in results])),
        "median_k3p5_open_fraction_48_96pix": float(np.median([
            min(row["mask_k3p5_wing_48_64_open_fraction"], row["mask_k3p5_wing_64_96_open_fraction"])
            for row in results])),
        "source_hashes": source_hashes,
        "hsc_mask_sections_sha256": mask_section_hashes,
        "star_stamp_quality_sha256": sha256(table_path),
        "script_sha256": sha256(Path(__file__)),
        "science_image_pixels_read": False,
        "hsc_requests_made": False,
        "psf_approved": False,
        "science_measurement": False,
        "limitations": [
            "Neighbor counts are HSC detections, not confirmed astrophysical contaminants.",
            "Open-mask fractions are support diagnostics and do not correct flux for missing pixels.",
            "Saturation and bright-object bits are reported separately from source-mask support.",
            "No profile values were recomputed and no candidate was promoted to an approved PSF star.",
        ],
        "completed_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    }
    summary_path = output / "summary.json"
    summary_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2, allow_nan=False) + "\n",
                            encoding="utf-8")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-root", type=Path, default=DEFAULT_SAMPLE)
    parser.add_argument("--ledger", type=Path, default=DEFAULT_LEDGER)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = run(args.sample_root.resolve(), args.ledger.resolve(), args.output.resolve())
    keys = ("status", "n_stamps", "n_target_frames", "n_with_saturation_in_core",
            "n_with_bright_object_bit_in_core", "n_with_neighbor_within_12pix",
            "n_k2p5_core_open_below_99pct", "n_k3p5_core_open_below_99pct",
            "median_k2p5_core_open_fraction", "median_k3p5_core_open_fraction",
            "n_with_neighbor_within_48pix", "n_with_neighbor_within_96pix",
            "median_nearest_neighbor_distance_pix", "science_image_pixels_read",
            "psf_approved", "science_measurement")
    print(json.dumps({key: result[key] for key in keys}, ensure_ascii=False, indent=2))
    print(f"summary_sha256={sha256(args.output / 'summary.json')}")


if __name__ == "__main__":
    main()
