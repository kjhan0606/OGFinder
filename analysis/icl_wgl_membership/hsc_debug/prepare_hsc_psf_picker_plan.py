#!/usr/bin/env python3
"""Prepare an offline HSC PSF Picker list for existing injection positions.

This script performs no network access. It reconstructs the exact cluster
centers and matched satellite positions used by the saved joint-injection
sample, then writes a coordinate-only PSF Picker request plan.
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
from astropy import units as u
from astropy.cosmology import FlatLambdaCDM
from astropy.io import fits
from astropy.wcs import WCS

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import measure_satellite_outer_profiles as satellite
from download_pilot import validate_fits

DEFAULT_DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
DEFAULT_INJECTIONS = DEFAULT_DATA / "joint_icl_satellite_redetection_2026-10-19_02/summary.json"
DEFAULT_MASKS = DEFAULT_DATA / "prestack_cluster_masks_2026-10-16_01"
DEFAULT_SATELLITE_PREP = DEFAULT_DATA / "satellite_detection_psf_2026-10-18_01"
DEFAULT_OUTPUT = DEFAULT_DATA / "icl_psf_picker_injection_plan_2026-10-19_02"
TEST_CLUSTER_ID = "10"
TEST_FILTER = "i2"
TEST_RA_DEG = 0.339748
TEST_DEC_DEG = 1.030287


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def frame_wcs(image: Path, science_hdu: int):
    with fits.open(image, memmap=True, checksum=True) as hdus:
        header = hdus[science_hdu].header.copy()
        if header.get("INHERIT") is True:
            header.extend(hdus[0].header, unique=True)
    wcs = WCS(header).celestial
    if not wcs.has_celestial:
        raise ValueError(f"image has no celestial WCS: {image}")
    return wcs


def frame_positions(frame: dict, record: dict, mask_record: dict,
                    satellite_prep: Path, h0: float, omega_m: float) -> list[dict]:
    ident = str(record["cluster_id"])
    image = Path(mask_record["image"])
    info = validate_fits(image, frame)
    if info["coadd_header_filter"] != mask_record["primary_filter"]:
        raise ValueError(f"FITS/filter mismatch for cluster {ident}")
    if str(info["coadd_header_filter"]).lower() != str(record["filter"]).lower():
        raise ValueError(f"injection/filter mismatch for cluster {ident}")

    prep = satellite_prep / f"cluster_{ident}"
    detections_path = prep / "detections.tsv"
    members_path = prep / "members.tsv"
    if not detections_path.is_file() or not members_path.is_file():
        raise FileNotFoundError(f"satellite-preparation inputs missing for {ident}")
    detections = satellite.parse_detections(satellite.read_tsv(detections_path))
    members = satellite.read_tsv(members_path)
    matched, _ = satellite.match_members(members, detections, info["pixel_scale_arcsec"])

    center_xy = np.asarray(info["center_pixel_xy"], dtype=float)
    pixel_scale = float(info["pixel_scale_arcsec"])
    cosmology = FlatLambdaCDM(H0=float(h0), Om0=float(omega_m), Tcmb0=0 * u.K)
    kpc_per_pixel = float(cosmology.kpc_proper_per_arcmin(float(frame["z"])).value
                           * pixel_scale / 60.0)
    targets = []
    for item in matched:
        det = detections[item["detection_index"]]
        radius = math.hypot(det["x"] - center_xy[0], det["y"] - center_xy[1]) * kpc_per_pixel
        if satellite.MIN_CLUSTER_RADIUS_KPC <= radius < satellite.MAX_CLUSTER_RADIUS_KPC:
            targets.append((det, radius))

    wcs = frame_wcs(image, info["science_hdu"])
    positions = [{
        "cluster_id": ident,
        "frame_filter": str(record["filter"]).lower(),
        "role": "cluster_center",
        "ra_deg": float(frame["ra_deg"]),
        "dec_deg": float(frame["dec_deg"]),
        "projected_radius_kpc": 0.0,
    }]
    for index, (det, radius) in enumerate(targets):
        sky = wcs.pixel_to_world(float(det["x"]), float(det["y"]))
        ra, dec = float(sky.ra.deg), float(sky.dec.deg)
        if not (math.isfinite(ra) and 0 <= ra < 360 and math.isfinite(dec) and -90 <= dec <= 90):
            raise ValueError(f"invalid member sky position for cluster {ident}")
        positions.append({
            "cluster_id": ident,
            "frame_filter": str(record["filter"]).lower(),
            "role": "satellite_injection_center",
            "target_index": index,
            "ra_deg": ra,
            "dec_deg": dec,
            "projected_radius_kpc": radius,
        })

    expected = int(record["support"]["n_injected_targets_in_50_500_kpc"])
    if len(targets) != expected:
        raise ValueError(f"satellite target count changed for {ident}: {len(targets)} != {expected}")
    return positions


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=DEFAULT_DATA)
    parser.add_argument("--injections", type=Path, default=DEFAULT_INJECTIONS)
    parser.add_argument("--masks", type=Path, default=DEFAULT_MASKS)
    parser.add_argument("--satellite-prep", type=Path, default=DEFAULT_SATELLITE_PREP)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    output = args.output.resolve()
    if output.exists():
        parser.error(f"refusing to overwrite existing output: {output}")

    summary = json.loads(args.injections.read_text(encoding="utf-8"))
    records = summary.get("records", [])
    if len(records) != 57 or summary.get("n_completed") != 57:
        raise ValueError("expected the verified 57-frame joint-injection sample")
    manifest = json.loads((args.data / "sample_manifest.json").read_text(encoding="utf-8"))
    frames = {str(row["id"]): row for row in manifest["clusters"]}
    if len(frames) != len(manifest["clusters"]):
        raise ValueError("sample manifest has duplicate cluster IDs")
    mask_rows = read_tsv(args.masks / "records.tsv")
    mask_records = {str(row["id"]): row for row in mask_rows}
    if len(mask_records) != len(mask_rows):
        raise ValueError("mask records have duplicate cluster IDs")

    positions = []
    for record in records:
        ident = str(record["cluster_id"])
        if ident not in frames or ident not in mask_records:
            raise ValueError(f"missing frame or mask record for cluster {ident}")
        frame = dict(frames[ident])
        frame["image"] = mask_records[ident]["image"]
        positions.extend(frame_positions(frame, record, mask_records[ident],
                                         args.satellite_prep,
                                         manifest.get("h0", 70.0),
                                         manifest.get("omega_m", 0.3)))

    if len(positions) != 745:
        raise ValueError(f"expected 745 total positions (57+688); got {len(positions)}")
    test_rows = [p for p in positions if p["cluster_id"] == TEST_CLUSTER_ID
                 and p["frame_filter"] == TEST_FILTER and p["role"] == "cluster_center"]
    if len(test_rows) != 1:
        raise ValueError("approved cluster-10 i2 PSF test position is not unique in the selected sample")
    test = test_rows[0]
    if not (math.isclose(test["ra_deg"], TEST_RA_DEG, abs_tol=1e-9)
            and math.isclose(test["dec_deg"], TEST_DEC_DEG, abs_tol=1e-9)):
        raise ValueError("cluster-10 PSF test coordinate differs from the approved request")

    # Collapse exact duplicate sky positions while preserving all source roles.
    unique = {}
    for row in positions:
        key = (round(row["ra_deg"], 9), round(row["dec_deg"], 9))
        if key not in unique:
            unique[key] = {**row, "roles": []}
        unique[key]["roles"].append({
            "cluster_id": row["cluster_id"], "frame_filter": row["frame_filter"],
            "role": row["role"], "target_index": row.get("target_index"),
        })
    rows = list(unique.values())
    test_key = (round(TEST_RA_DEG, 9), round(TEST_DEC_DEG, 9))
    if test_key not in unique:
        raise ValueError("test coordinate was lost during duplicate handling")
    remaining = [row for row in rows if (round(row["ra_deg"], 9), round(row["dec_deg"], 9)) != test_key]
    if len(remaining) != 744:
        raise ValueError(f"approved remaining-PSF count is 744; unique remainder is {len(remaining)}")

    output.mkdir(parents=True, exist_ok=False)
    coordinate_path = output / "psf_picker_coordinates.txt"
    with coordinate_path.open("x", encoding="ascii", newline="") as stream:
        stream.write("#? ra dec\n")
        for row in remaining:
            stream.write(f"{row['ra_deg']:.9f} {row['dec_deg']:+.9f}\n")
    mapping_path = output / "position_mapping.json"
    mapping_path.write_text(json.dumps(remaining, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    plan = {
        "schema": "astrafek-hsc-psf-picker-injection-plan-v1",
        "status": "prepared_not_submitted",
        "user_approved_requests": 744,
        "n_original_frame_records": len(records),
        "n_cluster_centers_total": 57,
        "n_satellite_injection_centers_total": 688,
        "n_unique_total_positions": len(rows),
        "n_previously_received_test_positions_excluded": 1,
        "n_requests_remaining": len(remaining),
        "rerun": "pdr3_wide",
        "type": "coadd",
        "filter": "HSC-I",
        "centered": True,
        "coordinate_list": str(coordinate_path),
        "coordinate_list_sha256": sha256(coordinate_path),
        "position_mapping": str(mapping_path),
        "position_mapping_sha256": sha256(mapping_path),
        "injection_summary": str(args.injections.resolve()),
        "injection_summary_sha256": sha256(args.injections),
        "sample_manifest_sha256": sha256(args.data / "sample_manifest.json"),
        "output_name_template": str(output / "psf_{tract}_{patch[0]}_{patch[1]}_{ra:.9f}_{dec:+.9f}_{filter}"),
        "images_requested": False,
        "shape_rows_requested": False,
        "photoz_requested": False,
        "science_measurement": False,
    }
    plan_path = output / "psf_picker_plan.json"
    plan_path.write_text(json.dumps(plan, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({k: plan[k] for k in (
        "status", "n_original_frame_records", "n_cluster_centers_total",
        "n_satellite_injection_centers_total", "n_unique_total_positions",
        "n_requests_remaining", "filter", "rerun", "type",
        "coordinate_list_sha256", "position_mapping_sha256",
    )}, indent=2))
    print(f"plan_sha256={sha256(plan_path)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
