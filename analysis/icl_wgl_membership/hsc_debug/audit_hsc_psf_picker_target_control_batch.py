#!/usr/bin/env python3
"""Audit target/control-center HSC PDR3 PSF Picker FITS against its plan."""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
import math
from pathlib import Path
import re

import numpy as np
from astropy.io import fits


EXPECTED_PLAN_SCHEMA = "astrafek-hsc-psf-picker-target-control-center-plan-v1"
EXPECTED_AUDIT_SCHEMA = "astrafek-hsc-psf-picker-target-control-center-audit-v1"
EXPECTED_REQUESTS = 1542
FILENAME_RE = re.compile(
    r"^psf_(\d+)_(\d+)_(\d+)_(-?\d+\.\d+)_([+-]?\d+\.\d+)_(HSC-I)\.fits$"
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _coordinate_key(ra: float, dec: float) -> tuple[float, float]:
    return round(float(ra), 9), round(float(dec), 9)


def _load_plan(plan_dir: Path, expected_requests: int):
    plan_path = plan_dir / "psf_picker_plan.json"
    mapping_path = plan_dir / "position_mapping.json"
    coordinate_path = plan_dir / "psf_picker_coordinates.txt"
    plan = json.loads(plan_path.read_text(encoding="utf-8"))
    mapping = json.loads(mapping_path.read_text(encoding="utf-8"))
    if plan.get("schema") != EXPECTED_PLAN_SCHEMA:
        raise ValueError("unexpected target/control PSF plan schema")
    # This offline plan is immutable and therefore remains marked
    # prepared_not_submitted even after a separately approved retrieval.
    if plan.get("status") != "prepared_not_submitted":
        raise ValueError("unexpected status in immutable offline plan")
    if (plan.get("n_requests_remaining") != expected_requests
            or plan.get("filter") != "HSC-I"
            or plan.get("rerun") != "pdr3_wide"
            or plan.get("type") != "coadd"
            or plan.get("centered") is not True):
        raise ValueError("plan does not match the requested PSF batch")
    if sha256(coordinate_path) != plan.get("coordinate_list_sha256"):
        raise ValueError("coordinate-list checksum differs from the plan")
    if sha256(mapping_path) != plan.get("position_mapping_sha256"):
        raise ValueError("position-mapping checksum differs from the plan")
    if not isinstance(mapping, list) or len(mapping) != expected_requests:
        raise ValueError(f"expected {expected_requests} mapped positions")
    expected = {}
    for row in mapping:
        key = _coordinate_key(row["ra_deg"], row["dec_deg"])
        if key in expected:
            raise ValueError(f"duplicate coordinate in mapping: {key}")
        if not isinstance(row.get("roles"), list) or not row["roles"]:
            raise ValueError(f"position has no target/control role: {key}")
        expected[key] = row

    lines = coordinate_path.read_text(encoding="ascii").splitlines()
    if not lines or lines[0].strip() != "#? ra dec":
        raise ValueError("unexpected coordinate-list header")
    coordinates = []
    for line in lines[1:]:
        columns = line.split()
        if len(columns) != 2:
            raise ValueError("malformed coordinate-list row")
        coordinates.append(_coordinate_key(*map(float, columns)))
    if len(coordinates) != expected_requests or len(set(coordinates)) != expected_requests:
        raise ValueError("coordinate list has a wrong count or duplicate rows")
    if set(coordinates) != set(expected):
        raise ValueError("coordinate list and position mapping disagree")
    return plan_path, mapping_path, coordinate_path, plan, expected


def audit_batch(plan_dir: Path, download_dir: Path, output: Path,
                *, expected_requests: int = EXPECTED_REQUESTS) -> dict:
    plan_dir, download_dir = Path(plan_dir), Path(download_dir)
    output = Path(output)
    plan_path, mapping_path, coordinate_path, plan, expected = _load_plan(
        plan_dir, expected_requests
    )
    if output.exists():
        raise FileExistsError(f"refusing to overwrite audit output: {output}")
    if not download_dir.is_dir():
        raise FileNotFoundError(download_dir)
    entries = list(download_dir.iterdir())
    if any(not p.is_file() or p.suffix.lower() not in {".fits", ".json"}
           for p in entries):
        raise ValueError("download directory contains an unexpected entry")
    fits_files = [p for p in entries if p.suffix.lower() == ".fits"]
    if len(fits_files) != expected_requests:
        raise ValueError(f"expected {expected_requests} FITS files, found {len(fits_files)}")

    found = {}
    records = []
    shapes = Counter()
    sums = []
    total_bytes = 0
    negative_pixels = 0
    max_centroid_offset = 0.0
    for path in sorted(fits_files):
        match = FILENAME_RE.fullmatch(path.name)
        if not match:
            raise ValueError(f"filename lacks expected HSC coordinate metadata: {path.name}")
        tract, patch_x, patch_y = map(int, match.group(1, 2, 3))
        ra, dec = map(float, match.group(4, 5))
        key = _coordinate_key(ra, dec)
        if key not in expected or key in found:
            raise ValueError(f"unexpected or duplicate returned coordinate: {path.name}")
        with fits.open(path, memmap=False, checksum=True) as hdus:
            if len(hdus) != 1 or hdus[0].data is None or np.ndim(hdus[0].data) != 2:
                raise ValueError(f"unexpected FITS structure: {path.name}")
            image = np.asarray(hdus[0].data, dtype=np.float64)
        if any(size % 2 != 1 for size in image.shape):
            raise ValueError(f"PSF stamp dimensions must be odd: {path.name}")
        if not np.isfinite(image).all():
            raise ValueError(f"non-finite PSF pixels: {path.name}")
        flux_sum = float(image.sum())
        if not math.isfinite(flux_sum) or abs(flux_sum - 1.0) > 1e-3:
            raise ValueError(f"PSF is not unit-normalized: {path.name}; sum={flux_sum}")
        positive = np.clip(image, 0, None)
        positive_sum = float(positive.sum())
        if positive_sum <= 0:
            raise ValueError(f"PSF has no positive flux: {path.name}")
        yy, xx = np.indices(image.shape)
        cx = float((positive * xx).sum() / positive_sum)
        cy = float((positive * yy).sum() / positive_sum)
        offset = float(np.hypot(cx - (image.shape[1] - 1) / 2,
                                cy - (image.shape[0] - 1) / 2))
        if offset > 1.0:
            raise ValueError(f"PSF centroid offset exceeds 1 pixel: {path.name}; {offset}")

        row = expected[key]
        nbytes = path.stat().st_size
        found[key] = path.name
        shapes[str(image.shape)] += 1
        sums.append(flux_sum)
        total_bytes += nbytes
        negative_pixels += int(np.count_nonzero(image < 0))
        max_centroid_offset = max(max_centroid_offset, offset)
        records.append({
            "file": path.name,
            "bytes": nbytes,
            "sha256": sha256(path),
            "tract": tract,
            "patch": [patch_x, patch_y],
            "filter": "HSC-I",
            "ra_deg_filename": ra,
            "dec_deg_filename": dec,
            "roles": row["roles"],
        })
    if set(found) != set(expected):
        raise ValueError(f"coordinate coverage incomplete: {len(found)}/{expected_requests}")

    result = {
        "schema": EXPECTED_AUDIT_SCHEMA,
        "status": "verified",
        "plan": str(plan_path.resolve()),
        "plan_sha256": sha256(plan_path),
        "coordinate_list_sha256": sha256(coordinate_path),
        "position_mapping_sha256": sha256(mapping_path),
        "download_directory": str(download_dir.resolve()),
        "n_approved": expected_requests,
        "n_returned": len(records),
        "n_matched_coordinates": len(found),
        "rerun": plan["rerun"],
        "type": plan["type"],
        "filter": plan["filter"],
        "centered": plan["centered"],
        "total_bytes": total_bytes,
        "shape_counts": dict(shapes),
        "file_size_range_bytes": [min(r["bytes"] for r in records),
                                  max(r["bytes"] for r in records)],
        "sum_range": [min(sums), max(sums)],
        "negative_pixel_total": negative_pixels,
        "max_positive_flux_centroid_offset_pixels": max_centroid_offset,
        "credentials_saved": False,
        "images_returned": False,
        "shape_rows_returned": False,
        "photoz_rows_returned": False,
        "science_measurement": False,
        "psf_approved": False,
        "files": records,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n",
                      encoding="utf-8")
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan-dir", type=Path, required=True)
    parser.add_argument("--download-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = audit_batch(args.plan_dir, args.download_dir, args.output)
    print(json.dumps({key: result[key] for key in (
        "status", "n_approved", "n_returned", "n_matched_coordinates",
        "total_bytes", "shape_counts", "file_size_range_bytes", "sum_range",
        "negative_pixel_total", "max_positive_flux_centroid_offset_pixels",
        "credentials_saved", "images_returned", "shape_rows_returned",
        "photoz_rows_returned", "science_measurement", "psf_approved",
    )}, indent=2))
    print(f"audit_sha256={sha256(args.output)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
