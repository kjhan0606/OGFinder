#!/usr/bin/env python3
"""Audit matched-control coverage at Gaia-star positions without reading pixels.

The audit maps each target-star position to the same angular offset from its
matched control center, then confirms that both existing HSC-I coadds contain
the complete 401-pixel diagnostic footprint. It does not estimate PSF wings.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import platform
from pathlib import Path

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS
import astropy.units as u

from fetch_ancillary_catalogues import sha

BASE = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
STAMP_RADIUS = 200


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def center_offset_pixel(
    target_wcs: WCS,
    control_wcs: WCS,
    target_center: SkyCoord,
    control_center: SkyCoord,
    target_xy: tuple[float, float],
) -> tuple[float, float, float, float]:
    """Map a pixel at target-relative sky offset to the control image.

    Returns control x/y and east/north angular offsets in arcsec.
    """
    target_star = SkyCoord(*target_wcs.pixel_to_world_values(*target_xy), unit="deg", frame="icrs")
    east, north = target_center.spherical_offsets_to(target_star)
    control_sky = control_center.spherical_offsets_by(east, north)
    control_x, control_y = control_wcs.world_to_pixel(control_sky)
    return float(control_x), float(control_y), float(east.to_value(u.arcsec)), float(north.to_value(u.arcsec))


def stamp_inside(shape: tuple[int, int], x: float, y: float, radius: int = STAMP_RADIUS) -> bool:
    height, width = shape
    return (math.isfinite(x) and math.isfinite(y)
            and radius <= x < width - radius and radius <= y < height - radius)


def filters_match(target: dict, control: dict) -> bool:
    return (target.get("band") == control.get("band")
            and target.get("coadd_header_filter") == control.get("coadd_header_filter"))


def _image_header(path: Path, hdu_index: int, expected_shape: tuple[int, int]) -> tuple[WCS, tuple[int, int], str]:
    if not path.is_file():
        raise FileNotFoundError(path)
    with fits.open(path, memmap=True, lazy_load_hdus=True) as hdus:
        if not 0 <= hdu_index < len(hdus):
            raise ValueError(f"invalid science HDU {hdu_index}: {path}")
        hdu = hdus[hdu_index]
        shape = tuple(int(v) for v in hdu.shape)
        if shape != expected_shape:
            raise ValueError(f"image shape differs from download audit for {path}: {shape} != {expected_shape}")
        header = hdu.header.copy()
        primary_filter = str(hdus[0].header.get("FILTER", "")).strip()
    wcs = WCS(header).celestial
    if not wcs.has_celestial or not wcs.pixel_n_dim == 2:
        raise ValueError(f"science image lacks a two-axis celestial WCS: {path}")
    return wcs, shape, primary_filter


def audit(sample_root: Path, ledger_path: Path) -> dict:
    target_dir = sample_root / "download_01"
    control_dir = sample_root / "control_download_01"
    target_manifest_path = target_dir / "download_manifest.json"
    control_manifest_path = control_dir / "download_manifest.json"
    target_status_path = target_dir / "download_status.json"
    control_status_path = control_dir / "download_status.json"
    for path in (target_manifest_path, control_manifest_path, target_status_path, control_status_path, ledger_path):
        if not path.is_file():
            raise FileNotFoundError(path)

    target_doc = json.loads(target_manifest_path.read_text(encoding="utf-8"))
    control_doc = json.loads(control_manifest_path.read_text(encoding="utf-8"))
    target_status = json.loads(target_status_path.read_text(encoding="utf-8"))
    control_status = json.loads(control_status_path.read_text(encoding="utf-8"))
    if target_doc.get("band") != "HSC-I" or control_doc.get("band") != "HSC-I":
        raise ValueError("target/control parent manifests are not HSC-I")

    target_by_id = {str(row["id"]): row for row in target_doc["clusters"]}
    control_by_target: dict[str, list[dict]] = {}
    for row in control_doc["clusters"]:
        control_by_target.setdefault(str(row["target_cluster_id"]), []).append(row)
    target_status_by_id = {str(row["id"]): row for row in target_status["completed"]}
    control_status_by_id = {str(row["id"]): row for row in control_status["completed"]}
    if len(target_by_id) != 800 or len(control_status_by_id) != int(control_status.get("n_expected", -1)):
        raise ValueError("target/control download manifests or completed status have unexpected coverage")

    ledger = read_tsv(ledger_path)
    if not ledger or not {"candidate_id", "cluster_id", "stamp", "stamp_sha256"} <= set(ledger[0]):
        raise ValueError("unexpected Gaia star-ledger schema")
    stamped = [row for row in ledger if row.get("stamp")]
    if not stamped:
        raise ValueError("star ledger contains no complete stamps")

    # Cache only WCS/header objects for the frames in this ledger. No image pixels are read.
    wcs_cache: dict[tuple[str, str], tuple[WCS, tuple[int, int], dict]] = {}
    output = []
    for row in stamped:
        target_id = str(row["cluster_id"])
        if target_id not in target_by_id or len(control_by_target.get(target_id, [])) != 1:
            raise ValueError(f"missing or nonunique target/control pair for {target_id}")
        target = target_by_id[target_id]
        control = control_by_target[target_id][0]
        target_record = target_status_by_id.get(target_id)
        control_record = control_status_by_id.get(str(control["id"]))
        if target_record is None or control_record is None:
            raise ValueError(f"target/control image lacks a completed download record for {target_id}")
        if sha(Path(row["stamp"])) != row["stamp_sha256"]:
            raise ValueError(f"star stamp hash mismatch for candidate {row['candidate_id']}")

        if not filters_match(target, control):
            output.append({
                "candidate_id": row["candidate_id"], "target_cluster_id": target_id,
                "control_id": control["id"], "control_random_object_id": control["random_object_id"],
                "target_filter": target["coadd_header_filter"], "control_filter": control["coadd_header_filter"],
                "target_star_x": None, "target_star_y": None, "control_x": None, "control_y": None,
                "east_offset_arcsec": None, "north_offset_arcsec": None,
                "target_stamp_inside": None, "control_stamp_inside": None,
                "control_margin_x_pix": None, "control_margin_y_pix": None,
                "pair_status": "blocked_filter_mismatch",
            })
            continue

        target_key = ("target", target_id)
        control_key = ("control", target_id)
        if target_key not in wcs_cache:
            target_path = target_dir / "images" / target_record["file"]
            wcs_cache[target_key] = _image_header(
                target_path, int(target_record["review"]["science_hdu"]),
                tuple(target_record["review"]["shape"]),
            )
        if control_key not in wcs_cache:
            control_path = control_dir / "images" / control_record["file"]
            wcs_cache[control_key] = _image_header(
                control_path, int(control_record["review"]["science_hdu"]),
                tuple(control_record["review"]["shape"]),
            )
        target_wcs, target_shape, target_filter = wcs_cache[target_key]
        control_wcs, control_shape, control_filter = wcs_cache[control_key]
        for actual, expected in ((target_filter, target["coadd_header_filter"]),
                                 (control_filter, control["coadd_header_filter"])):
            if actual != str(expected).strip():
                raise ValueError(f"FITS FILTER header disagrees with manifest for target {target_id}")

        # Native coordinates are stored in the original cutout pixel grid; cross-check
        # their local stamp coordinates before using them in the target WCS.
        with fits.open(row["stamp"], memmap=True) as stamp_hdus:
            stamp_header = stamp_hdus[0].header
        x = float(stamp_header["X0"]) + float(stamp_header["XCENTER"])
        y = float(stamp_header["Y0"]) + float(stamp_header["YCENTER"])
        if not np.allclose([x, y], [float(row["native_x"]), float(row["native_y"])], rtol=0, atol=1e-6):
            raise ValueError(f"ledger and stamp native coordinates disagree for candidate {row['candidate_id']}")

        target_center = SkyCoord(float(target["ra_deg"]) * u.deg, float(target["dec_deg"]) * u.deg, frame="icrs")
        control_center = SkyCoord(float(control["ra_deg"]) * u.deg, float(control["dec_deg"]) * u.deg, frame="icrs")
        cx, cy, east_arcsec, north_arcsec = center_offset_pixel(
            target_wcs, control_wcs, target_center, control_center, (x, y),
        )
        target_star_sky = target_wcs.pixel_to_world(x, y)
        tx, ty = target_wcs.world_to_pixel(target_star_sky)
        if not stamp_inside(target_shape, float(tx), float(ty)):
            raise ValueError(f"star stamp footprint no longer fits target image for {row['candidate_id']}")
        output.append({
            "candidate_id": row["candidate_id"], "target_cluster_id": target_id,
            "control_id": control["id"], "control_random_object_id": control["random_object_id"],
            "target_filter": target["coadd_header_filter"], "control_filter": control["coadd_header_filter"],
            "target_star_x": float(tx), "target_star_y": float(ty),
            "control_x": cx, "control_y": cy,
            "east_offset_arcsec": east_arcsec, "north_offset_arcsec": north_arcsec,
            "target_stamp_inside": True, "control_stamp_inside": stamp_inside(control_shape, cx, cy),
            "control_margin_x_pix": min(cx - STAMP_RADIUS, control_shape[1] - STAMP_RADIUS - cx),
            "control_margin_y_pix": min(cy - STAMP_RADIUS, control_shape[0] - STAMP_RADIUS - cy),
            "pair_status": "same_filter_offset_available" if stamp_inside(control_shape, cx, cy)
            else "blocked_control_stamp_edge",
        })

    return {
        "schema": "astrafek-empirical-psf-matched-control-support-v1",
        "status": "complete_header_wcs_preflight",
        "n_stamps": len(stamped), "n_unique_target_frames": len({row["cluster_id"] for row in stamped}),
        "n_same_filter_pairs": sum(row["pair_status"] != "blocked_filter_mismatch" for row in output),
        "n_filter_mismatch": sum(row["pair_status"] == "blocked_filter_mismatch" for row in output),
        "n_control_401px_footprints_inside": sum(row["pair_status"] == "same_filter_offset_available" for row in output),
        "n_control_footprints_outside": sum(row["pair_status"] == "blocked_control_stamp_edge" for row in output),
        "target_manifest_sha256": sha(target_manifest_path), "control_manifest_sha256": sha(control_manifest_path),
        "target_status_sha256": sha(target_status_path), "control_status_sha256": sha(control_status_path),
        "star_ledger_sha256": sha(ledger_path),
        "n_target_wcs_headers": len([key for key in wcs_cache if key[0] == "target"]),
        "n_control_wcs_headers": len([key for key in wcs_cache if key[0] == "control"]),
        "pixel_data_read": False, "hsc_query_submitted": False, "data_downloaded": False,
        "science_measurement": False, "psf_approved": False,
        "interpretation": "Checks only whether matched controls can supply same-offset 401x401 image stamps; no pixels or PSF wings were measured.",
        "rows": output,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-root", type=Path, default=BASE)
    parser.add_argument("--ledger", type=Path, default=BASE / "empirical_psf_gaia_full_sample_2026-10-20_02/star_ledger.tsv")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if not platform.node().lower().startswith("lageunha"):
        parser.error("Run this FITS-header preflight on lageunha CPU")
    if args.output.exists():
        parser.error(f"Refusing to overwrite output directory: {args.output}")
    result = audit(args.sample_root.resolve(), args.ledger.resolve())
    args.output.mkdir(parents=True, exist_ok=False)
    rows_path = args.output / "matched_control_support.tsv"
    rows = result.pop("rows")
    with rows_path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    result["rows_path"] = str(rows_path.resolve())
    result["rows_sha256"] = sha(rows_path)
    summary_path = args.output / "summary.json"
    summary_path.write_text(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({key: result[key] for key in (
        "status", "n_stamps", "n_unique_target_frames", "n_same_filter_pairs",
        "n_filter_mismatch", "n_control_401px_footprints_inside", "n_control_footprints_outside",
        "pixel_data_read", "hsc_query_submitted", "data_downloaded",
    )}, indent=2))
    print(f"summary_sha256={sha(summary_path)}")


if __name__ == "__main__":
    main()
