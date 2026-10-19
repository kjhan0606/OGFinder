#!/usr/bin/env python3
"""Audit downloaded HSC PSF Picker FITS files against an approved plan."""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path
import re

import numpy as np
from astropy.io import fits


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan-dir", type=Path, required=True)
    parser.add_argument("--download-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    plan_path = args.plan_dir / "psf_picker_plan.json"
    mapping_path = args.plan_dir / "position_mapping.json"
    coordinate_path = args.plan_dir / "psf_picker_coordinates.txt"
    plan = json.loads(plan_path.read_text(encoding="utf-8"))
    mapping = json.loads(mapping_path.read_text(encoding="utf-8"))
    if plan.get("schema") != "astrafek-hsc-psf-picker-injection-plan-v1":
        raise ValueError("unexpected PSF Picker plan schema")
    if plan.get("status") != "prepared_not_submitted" or plan.get("user_approved_requests") != 744:
        raise ValueError("plan approval scope is not the approved 744-position request")
    if sha256(coordinate_path) != plan.get("coordinate_list_sha256"):
        raise ValueError("coordinate-list checksum differs from the approved plan")
    if sha256(mapping_path) != plan.get("position_mapping_sha256"):
        raise ValueError("position mapping checksum differs from the approved plan")
    if len(mapping) != 744:
        raise ValueError("position mapping does not contain exactly 744 requests")

    expected = {(round(float(row["ra_deg"]), 9), round(float(row["dec_deg"]), 9)): row
                for row in mapping}
    if len(expected) != 744:
        raise ValueError("position mapping contains duplicate coordinates")
    lines = coordinate_path.read_text(encoding="ascii").splitlines()
    if not lines or lines[0].strip() != "#? ra dec":
        raise ValueError("unexpected coordinate-list header")
    coordinates = []
    for line in lines[1:]:
        cols = line.split()
        if len(cols) != 2:
            raise ValueError("coordinate list has a malformed row")
        coordinates.append((round(float(cols[0]), 9), round(float(cols[1]), 9)))
    if len(coordinates) != 744 or set(coordinates) != set(expected):
        raise ValueError("coordinate list and mapping do not contain the same 744 positions")

    file_re = re.compile(r"^psf_(\d+)_(\d+)_(\d+)_(-?\d+\.\d+)_([+-]?\d+\.\d+)_(HSC-I)\.fits$")
    entries = list(args.download_dir.iterdir())
    output_path = args.output.resolve()
    if any(
        not path.is_file()
        or (
            path.suffix.lower() != ".fits"
            and path.suffix.lower() != ".json"
            and path.resolve() != output_path
        )
        for path in entries
    ):
        raise ValueError("download directory contains a non-FITS or non-file entry")
    fits_entries = [path for path in entries if path.suffix.lower() == ".fits"]
    if len(fits_entries) != 744:
        raise ValueError(f"expected 744 FITS outputs, found {len(fits_entries)}")
    found = {}
    records = []
    shapes = Counter()
    sum_min = float("inf")
    sum_max = float("-inf")
    neg_total = 0
    max_centroid_offset = 0.0
    total_bytes = 0
    for path in sorted(fits_entries):
        match = file_re.fullmatch(path.name)
        if not match:
            raise ValueError(f"filename does not preserve tract/patch/position/filter: {path.name}")
        tract, patch_x, patch_y = map(int, match.group(1, 2, 3))
        ra, dec = float(match.group(4)), float(match.group(5))
        filter_name = match.group(6)
        key = (round(ra, 9), round(dec, 9))
        if key not in expected or key in found:
            raise ValueError(f"unexpected or duplicate PSF coordinate in {path.name}")
        with fits.open(path, memmap=False, checksum=True) as hdus:
            if len(hdus) != 1 or hdus[0].data is None or np.ndim(hdus[0].data) != 2:
                raise ValueError(f"unexpected FITS structure: {path.name}")
            image = np.asarray(hdus[0].data, dtype=np.float64)
            if not np.isfinite(image).all():
                raise ValueError(f"non-finite model PSF pixels: {path.name}")
            flux_sum = float(image.sum())
            if not np.isfinite(flux_sum) or abs(flux_sum - 1.0) > 1e-3:
                raise ValueError(f"model PSF is not unit-normalized: {path.name} sum={flux_sum}")
            positive = np.clip(image, 0, None)
            positive_sum = float(positive.sum())
            if positive_sum <= 0:
                raise ValueError(f"model PSF has no positive flux: {path.name}")
            yy, xx = np.indices(image.shape)
            cx = float((positive * xx).sum() / positive_sum)
            cy = float((positive * yy).sum() / positive_sum)
            geometric = ((image.shape[1] - 1) / 2, (image.shape[0] - 1) / 2)
            offset = float(np.hypot(cx - geometric[0], cy - geometric[1]))
            max_centroid_offset = max(max_centroid_offset, offset)
            negative = int(np.count_nonzero(image < 0))
            neg_total += negative
            sum_min, sum_max = min(sum_min, flux_sum), max(sum_max, flux_sum)
            shapes[str(image.shape)] += 1
        byte_count = path.stat().st_size
        total_bytes += byte_count
        found[key] = path.name
        row = expected[key]
        records.append({
            "file": path.name,
            "bytes": byte_count,
            "sha256": sha256(path),
            "tract": tract,
            "patch": [patch_x, patch_y],
            "filter": filter_name,
            "ra_deg_filename": ra,
            "dec_deg_filename": dec,
            "cluster_ids": sorted({str(r["cluster_id"]) for r in row["roles"]}),
            "roles": sorted({str(r["role"]) for r in row["roles"]}),
        })
    if set(found) != set(expected):
        raise ValueError("some approved coordinates have no returned PSF file")
    if max_centroid_offset > 1.0:
        raise ValueError(f"model PSF centroids depart by >1 pixel: {max_centroid_offset}")
    if args.output.exists():
        raise FileExistsError(f"refusing to overwrite audit output: {args.output}")

    result = {
        "schema": "astrafek-hsc-psf-picker-batch-audit-v1",
        "status": "verified",
        "plan": str(plan_path.resolve()),
        "plan_sha256": sha256(plan_path),
        "coordinate_list_sha256": sha256(coordinate_path),
        "position_mapping_sha256": sha256(mapping_path),
        "download_directory": str(args.download_dir.resolve()),
        "n_requested": len(expected),
        "n_returned": len(records),
        "n_matched_coordinates": len(found),
        "filter": "HSC-I",
        "rerun": "pdr3_wide",
        "type": "coadd",
        "total_bytes": total_bytes,
        "shape_counts": dict(shapes),
        "sum_range": [sum_min, sum_max],
        "negative_pixel_total": neg_total,
        "max_positive_flux_centroid_offset_pixels": max_centroid_offset,
        "credentials_saved": False,
        "images_returned": False,
        "shape_rows_returned": False,
        "photoz_rows_returned": False,
        "science_measurement": False,
        "files": records,
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({k: result[k] for k in (
        "status", "n_requested", "n_returned", "total_bytes", "shape_counts",
        "sum_range", "negative_pixel_total", "max_positive_flux_centroid_offset_pixels",
        "credentials_saved", "images_returned", "shape_rows_returned", "photoz_rows_returned",
        "science_measurement",
    )}, indent=2))
    print(f"audit_sha256={sha256(args.output)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
