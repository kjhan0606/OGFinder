#!/usr/bin/env python3
"""Offline integrity and WCS audit for the approved S23B control cutouts."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from download_pilot import validate_fits
from prepare_s23b_icl_target_image_manifest import sha256


DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
DEFAULT_PLAN = DATA / "s23b_icl_control_image_plan_2026-11-07_02"
DEFAULT_WORK = DATA / "s23b_icl_control_image_plan_2026-11-07_01/control_download_01"
DEFAULT_OUTPUT = DATA / "s23b_icl_control_image_audit_2026-11-07_01"
EXPECTED_MANIFEST_SHA256 = "ae9eb406cfc179d9980031bb67f85f8ee71c4580ba59e97713278823db5c06df"
EXPECTED_N = 1041


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def audit(plan_dir: Path, work_dir: Path, output_dir: Path) -> dict:
    plan_dir, work_dir, output_dir = map(Path.resolve, (plan_dir, work_dir, output_dir))
    if output_dir.exists():
        raise FileExistsError(f"refusing to overwrite {output_dir}")
    manifest_path = plan_dir / "control_manifest.json"
    plan_summary_path = plan_dir / "summary.json"
    status_path = work_dir / "download_status.json"
    download_manifest_path = work_dir / "download_manifest.json"
    for path in (manifest_path, plan_summary_path, status_path, download_manifest_path):
        if not path.is_file():
            raise FileNotFoundError(path)

    manifest, plan, status, downloaded = map(read_json,
        (manifest_path, plan_summary_path, status_path, download_manifest_path))
    digest = sha256(manifest_path)
    if (digest != EXPECTED_MANIFEST_SHA256
            or plan.get("approved_manifest_sha256") != digest
            or plan.get("status") not in ("approved_resume_in_progress", "complete")
            or plan.get("science_measurement") is not False
            or manifest.get("schema") != "hsc-bcg-diffuse-sample-v1"
            or manifest.get("controls_ready_for_science") is not False
            or manifest.get("science_measurement") is not False):
        raise ValueError("control plan is not the approved pinned v2 manifest")
    if (status.get("status") != "complete" or status.get("manifest_sha256") != digest
            or status.get("n_expected") != EXPECTED_N
            or status.get("science_measurement") is not False
            or status.get("credentials_saved") is not False):
        raise ValueError("control download checkpoint is incomplete or outside approved scope")
    if (downloaded.get("schema") != "hsc-bcg-diffuse-sample-v1"
            or downloaded.get("science_measurement") is not False
            or len(downloaded.get("clusters", [])) != EXPECTED_N):
        raise ValueError("completed control download manifest has an unexpected contract")

    controls = manifest.get("clusters", [])
    status_rows = status.get("completed", [])
    download_rows = downloaded.get("clusters", [])
    if not (len(controls) == len(status_rows) == len(download_rows) == EXPECTED_N):
        raise ValueError("control frame counts do not match the approved 1,041-frame scope")
    expected = {str(row["id"]): row for row in controls}
    recorded = {str(row.get("id")): row for row in status_rows}
    products = {str(row["id"]): row for row in download_rows}
    if len(expected) != EXPECTED_N or set(expected) != set(recorded) or set(expected) != set(products):
        raise ValueError("control identities differ across manifest, checkpoint, and output manifest")
    targets = [str(row.get("target_cluster_id", "")) for row in controls]
    if len(set(targets)) != EXPECTED_N:
        raise ValueError("controls do not have a unique one-to-one target mapping")
    if list((work_dir / "images").glob("*.part")):
        raise ValueError("partial FITS files remain in the completed control image directory")

    rows = []
    total_bytes = 0
    for ident in sorted(expected):
        frame, item, output = expected[ident], recorded[ident], products[ident]
        if (frame.get("frame_role") != "matched_control"
                or frame.get("bcg_policy") != "none"
                or frame.get("members_radec_deg") != []
                or frame.get("randoms") != []
                or frame.get("science_measurement") is not False
                or float(frame.get("z", -1)) <= 0):
            raise ValueError(f"{ident}: invalid matched-control identity or redshift policy")
        filename = str(item.get("file", ""))
        if not filename or Path(filename).name != filename or output.get("image") != filename:
            raise ValueError(f"{ident}: unsafe or inconsistent FITS filename")
        path = work_dir / "images" / filename
        if not path.is_file() or path.stat().st_size != int(item.get("bytes", -1)):
            raise ValueError(f"{ident}: missing FITS or size mismatch")
        digest_file = sha256(path)
        if digest_file != item.get("sha256"):
            raise ValueError(f"{ident}: FITS SHA-256 mismatch")
        review = validate_fits(path, frame)
        if json.dumps(review, sort_keys=True) != json.dumps(item.get("review"), sort_keys=True):
            raise ValueError(f"{ident}: FITS/WCS review differs from download checkpoint")
        if (output.get("product") != "coadd/bg" or output.get("band") != "HSC-I"
                or output.get("coadd_header_filter") not in {"i", "i2"}
                or output.get("science_hdu") != review["science_hdu"]
                or str(output.get("target_cluster_id")) != str(frame["target_cluster_id"])):
            raise ValueError(f"{ident}: downloaded control manifest fields are inconsistent")
        if ident == "CTRL_S23B_456" and str(frame.get("random_object_id")) != "42740795815504337":
            raise ValueError("replacement control does not match the approved v2 position")
        total_bytes += path.stat().st_size
        rows.append({
            "id": ident, "target_cluster_id": str(frame["target_cluster_id"]),
            "random_object_id": str(frame["random_object_id"]), "file": filename,
            "bytes": path.stat().st_size, "sha256": digest_file,
            "filter": review["coadd_header_filter"], "shape_y": review["shape"][0],
            "shape_x": review["shape"][1], "pixel_scale_arcsec": review["pixel_scale_arcsec"],
            "center_x": review["center_pixel_xy"][0], "center_y": review["center_pixel_xy"][1],
            "wcs_alignment_verified": True,
        })

    output_dir.mkdir(parents=True, exist_ok=False)
    table_path = output_dir / "verified_control_frames.tsv"
    with table_path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    summary = {
        "schema": "astrafek-s23b-icl-control-image-audit-v1",
        "status": "complete_header_wcs_file_integrity_audit",
        "n_expected_frames": EXPECTED_N, "n_verified_frames": len(rows),
        "total_verified_file_bytes": total_bytes,
        "filter_family_counts": {band: sum(row["filter"] == band for row in rows) for band in ("i", "i2")},
        "product": "HSC PDR3 pdr3_wide coadd/bg", "planes": ["IMAGE", "MASK", "VARIANCE"],
        "manifest_sha256": digest, "science_measurement": False,
        "images_ready_for_stacking": False, "credentials_saved": False,
        "network_used_by_audit": False,
        "inputs_sha256": {str(path): sha256(path) for path in
                           (manifest_path, plan_summary_path, status_path, download_manifest_path)},
        "output_sha256": {table_path.name: sha256(table_path)},
        "limitations": [
            "This is a file integrity, FITS layout, filter, and WCS audit; it does not assess pixel usability or sky residuals.",
            "Control positions are matched random sky locations, never BCGs or satellite-member positions.",
            "Redshifts are copied from paired S23B targets only for angular-to-proper scale conversion.",
            "No extinction, PSF, sky-transfer, covariance, or ICL correction is applied.",
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
