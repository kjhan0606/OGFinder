#!/usr/bin/env python3
"""Independently audit the approved S23B target cutout download.

Local-only verification of the 1,043-frame HSC-I coadd/bg download. It checks
the exact approved target manifest, checkpoint inventory, per-file checksum,
FITS image/mask/variance layout, WCS registration, filter family, and cutout
footprint. It does not measure ICL or claim science readiness.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from download_pilot import validate_fits
from prepare_s23b_icl_target_image_manifest import sha256


DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
DEFAULT_PLAN = DATA / "s23b_icl_target_image_plan_2026-11-07_01"
DEFAULT_WORK = DEFAULT_PLAN / "target_download_01"
DEFAULT_OUTPUT = DATA / "s23b_icl_target_image_audit_2026-11-07_01"
EXPECTED_MANIFEST_SHA256 = "e487582ba40cd3e58da67704242db4d7219e381fb08f4b33fed47ce4ff268b44"


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def audit(plan_dir: Path, work_dir: Path, output_dir: Path) -> dict:
    plan_dir, work_dir = Path(plan_dir).resolve(), Path(work_dir).resolve()
    output_dir = Path(output_dir).resolve()
    if output_dir.exists():
        raise FileExistsError(f"refusing to overwrite {output_dir}")
    plan_summary_path = plan_dir / "summary.json"
    target_manifest_path = plan_dir / "target_manifest.json"
    status_path = work_dir / "download_status.json"
    download_manifest_path = work_dir / "download_manifest.json"
    for path in (plan_summary_path, target_manifest_path, status_path, download_manifest_path):
        if not path.is_file():
            raise FileNotFoundError(path)
    plan_summary = read_json(plan_summary_path)
    target_manifest = read_json(target_manifest_path)
    status = read_json(status_path)
    downloaded = read_json(download_manifest_path)
    manifest_sha = sha256(target_manifest_path)
    if (plan_summary.get("schema") != "astrafek-s23b-icl-target-image-plan-v1"
            or plan_summary.get("status") != "prepared_after_user_image_approval"
            or plan_summary.get("user_approval_received") is None
            or manifest_sha != EXPECTED_MANIFEST_SHA256
            or manifest_sha != plan_summary.get("outputs_sha256", {}).get("target_manifest.json")):
        raise ValueError("target image plan is not the pinned, approved manifest")
    if (status.get("status") != "complete" or status.get("manifest_sha256") != manifest_sha
            or status.get("n_expected") != 1043
            or status.get("science_measurement") is not False
            or status.get("credentials_saved") is not False):
        raise ValueError("download checkpoint is incomplete or exceeds the approved scope")
    if (downloaded.get("schema") != "hsc-bcg-diffuse-sample-v1"
            or downloaded.get("science_measurement") is not False
            or downloaded.get("bulk_download") is not True
            or len(downloaded.get("clusters", [])) != 1043):
        raise ValueError("completed download manifest has an unexpected sample contract")
    targets = target_manifest.get("clusters", [])
    if len(targets) != 1043 or len({str(row["id"]) for row in targets}) != 1043:
        raise ValueError("approved target manifest does not contain 1,043 unique frames")
    expected = {str(row["id"]): row for row in targets}
    status_rows = status.get("completed", [])
    if len(status_rows) != 1043 or {str(row.get("id")) for row in status_rows} != set(expected):
        raise ValueError("download checkpoint frame inventory differs from target manifest")
    output_rows = {str(row["id"]): row for row in downloaded["clusters"]}
    if set(output_rows) != set(expected):
        raise ValueError("download manifest frame inventory differs from target manifest")

    image_dir = work_dir / "images"
    if not image_dir.is_dir():
        raise FileNotFoundError(image_dir)
    if list(image_dir.glob("*.part")):
        raise ValueError("partial FITS files remain in the completed image directory")
    records, total_bytes = [], 0
    by_id = {str(row["id"]): row for row in status_rows}
    for frame_id in sorted(expected):
        source = expected[frame_id]
        item, output = by_id[frame_id], output_rows[frame_id]
        filename = str(item.get("file", ""))
        if not filename or Path(filename).name != filename or output.get("image") != filename:
            raise ValueError(f"unsafe or inconsistent downloaded filename for {frame_id}")
        path = image_dir / filename
        if not path.is_file() or path.stat().st_size != int(item.get("bytes", -1)):
            raise ValueError(f"missing file or byte count mismatch for {frame_id}")
        digest = sha256(path)
        if digest != item.get("sha256"):
            raise ValueError(f"downloaded FITS SHA-256 mismatch for {frame_id}")
        review = validate_fits(path, source)
        if json.dumps(review, sort_keys=True) != json.dumps(item.get("review"), sort_keys=True):
            raise ValueError(f"FITS/WCS review result differs from download checkpoint for {frame_id}")
        if (output.get("product") != "coadd/bg" or output.get("band") != "HSC-I"
                or output.get("coadd_header_filter") not in {"i", "i2"}
                or output.get("science_hdu") != review["science_hdu"]):
            raise ValueError(f"download manifest has inconsistent product/filter evidence for {frame_id}")
        total_bytes += path.stat().st_size
        records.append({
            "id": frame_id, "file": filename, "bytes": path.stat().st_size,
            "sha256": digest, "primary_filter": review["coadd_header_filter"],
            "shape_y": review["shape"][0], "shape_x": review["shape"][1],
            "pixel_scale_arcsec": review["pixel_scale_arcsec"],
            "center_x": review["center_pixel_xy"][0], "center_y": review["center_pixel_xy"][1],
            "requested_half_width_arcsec": source["requested_half_width_arcsec"],
            "wcs_hdu_alignment_verified": True,
        })
    output_dir.mkdir(parents=True, exist_ok=False)
    rows_path = output_dir / "verified_target_frames.tsv"
    with rows_path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(records[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(records)
    summary = {
        "schema": "astrafek-s23b-icl-target-image-audit-v1",
        "status": "complete_header_wcs_file_integrity_audit",
        "n_expected_frames": len(expected), "n_verified_frames": len(records),
        "total_verified_file_bytes": total_bytes,
        "filter_family_counts": {
            band: sum(row["primary_filter"] == band for row in records)
            for band in ("i", "i2")
        },
        "product": "HSC PDR3 pdr3_wide coadd/bg",
        "planes": ["IMAGE", "MASK", "VARIANCE"],
        "science_measurement": False, "images_ready_for_stacking": False,
        "credentials_saved": False, "network_used_by_audit": False,
        "inputs_sha256": {
            str(path): sha256(path) for path in (
                plan_summary_path, target_manifest_path, status_path, download_manifest_path)
        },
        "output_sha256": {rows_path.name: sha256(rows_path)},
        "limitations": [
            "this audit verifies file integrity, headers, WCS alignment, filter evidence, and requested footprint only",
            "it does not validate pixel-level usable coverage, sky residuals, noise covariance, masking completeness, or PSF wings",
            "centers are S23B CAMIRA cluster centroids; BCG positions and member catalogs remain unavailable in this sample manifest",
            "no ICL stack or satellite association has been measured by this audit",
        ],
    }
    summary_path = output_dir / "summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({key: summary[key] for key in (
        "status", "n_expected_frames", "n_verified_frames", "total_verified_file_bytes",
        "filter_family_counts", "science_measurement", "images_ready_for_stacking")}, indent=2))
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan-dir", type=Path, default=DEFAULT_PLAN)
    parser.add_argument("--work-dir", type=Path, default=DEFAULT_WORK)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    audit(args.plan_dir, args.work_dir, args.output_dir)


if __name__ == "__main__":
    main()
