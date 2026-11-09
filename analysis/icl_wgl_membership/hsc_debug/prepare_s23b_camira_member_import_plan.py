#!/usr/bin/env python3
"""Prepare an offline, exact-scope intake plan for the S23B CAMIRA member file.

This script reads only the already-approved local target manifest. It does not
contact Google Drive, download the member catalog, or read any member rows.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
DEFAULT_TARGET_MANIFEST = (
    DATA / "s23b_icl_target_image_plan_2026-11-07_01/target_manifest.json"
)
DEFAULT_OUTPUT = DATA / "s23b_camira_member_import_plan_2026-11-07_01"
EXPECTED_TARGET_MANIFEST_SHA256 = (
    "e487582ba40cd3e58da67704242db4d7219e381fb08f4b33fed47ce4ff268b44"
)
EXPECTED_CATALOG_SHA256 = "5af057028feab7976cd3ed389a7246cbe73add8bbb96891101f0f692b1de4d71"
EXPECTED_TARGETS = 1043
SOURCE_FILE_ID = "1DEnjGj9P1ilUooEm2U5_A-qXUrxeiJdE"
SOURCE_NAME = "camira_s23b_wide_sm_member.dat.gz"
SOURCE_BYTES = 10626590
SOURCE_MIME_TYPE = "application/x-gzip"
SOURCE_URL = (
    "https://drive.google.com/file/d/1DEnjGj9P1ilUooEm2U5_A-qXUrxeiJdE/view"
    "?usp=drivesdk"
)
SOURCE_FOLDER_URL = (
    "https://drive.google.com/drive/folders/19dKiNs7Wdq44X3AtYy2opOf8ulHTLWeM"
    "?usp=sharing"
)
SOURCE_README_URL = (
    "https://github.com/oguri/cluster_catalogs/tree/main/hsc_s23b_camira"
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def write_json(path: Path, value: dict) -> None:
    path.write_text(
        json.dumps(value, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def prepare(target_manifest_path: Path, output_dir: Path) -> dict:
    target_manifest_path = Path(target_manifest_path).resolve()
    output_dir = Path(output_dir).resolve()
    if output_dir.exists():
        raise FileExistsError(f"refusing to overwrite {output_dir}")
    target_sha = sha256(target_manifest_path)
    if target_sha != EXPECTED_TARGET_MANIFEST_SHA256:
        raise ValueError("target manifest is not the pinned, approved 1,043-center sample")
    manifest = json.loads(target_manifest_path.read_text(encoding="utf-8"))
    rows = manifest.get("clusters", [])
    ids = [str(row["id"]) for row in rows]
    if (manifest.get("catalog_sha256") != EXPECTED_CATALOG_SHA256
            or len(rows) != EXPECTED_TARGETS
            or len(set(ids)) != EXPECTED_TARGETS
            or manifest.get("member_catalog_sha256") is not None
            or any(row.get("members_radec_deg") for row in rows)):
        raise ValueError("target member catalog state differs from the expected S23B sample")
    target_ids_sha = hashlib.sha256(
        ("\n".join(sorted(ids)) + "\n").encode("utf-8")
    ).hexdigest()
    plan = {
        "schema": "astrafek-s23b-camira-member-import-plan-v1",
        "status": "awaiting_explicit_download_approval",
        "approval_required": True,
        "scope": {
            "files": 1,
            "name": SOURCE_NAME,
            "drive_file_id": SOURCE_FILE_ID,
            "mime_type": SOURCE_MIME_TYPE,
            "compressed_bytes": SOURCE_BYTES,
            "compressed_mib": SOURCE_BYTES / 1024**2,
            "source_url": SOURCE_URL,
            "source_folder_url": SOURCE_FOLDER_URL,
            "source_catalog_readme": SOURCE_README_URL,
            "sample_variant": "S23B CAMIRA Wide with bright-star mask",
            "intended_use": "positions and CAMIRA membership fields for the pinned 1,043-target ICL/satellite sample",
            "excluded": [
                "all other Drive files, including unmasked CAMIRA members and LRG FITS catalogs",
                "HSC images, shape rows, and photo-z rows",
            ],
        },
        "destination_if_approved": str(
            DATA / "s23b_camira_member_catalog_2026-11-07_01" / SOURCE_NAME
        ),
        "target_sample": {
            "n_targets": EXPECTED_TARGETS,
            "target_manifest": str(target_manifest_path),
            "target_manifest_sha256": target_sha,
            "target_catalog_sha256": EXPECTED_CATALOG_SHA256,
            "sorted_target_id_list_sha256": target_ids_sha,
            "member_positions_present_before_import": False,
        },
        "post_receipt_checks": [
            "verify exact compressed byte count and compute SHA-256",
            "verify gzip stream integrity before parsing",
            "record the table header/schema without assuming column meanings",
            "match catalog/cluster IDs to the pinned S23B target list",
            "measure member-coordinate coverage inside the approved image footprints",
        ],
        "source_visibility": "Drive metadata says shared/can-download; connector did not verify public visibility",
        "download_performed": False,
        "network_used_by_this_plan": False,
        "credentials_saved": False,
        "science_measurement": False,
        "ready_for_analysis": False,
    }
    output_dir.mkdir(parents=True, exist_ok=False)
    plan_path = output_dir / "plan.json"
    write_json(plan_path, plan)
    summary = {
        "status": plan["status"],
        "plan_path": str(plan_path),
        "plan_sha256": sha256(plan_path),
        "file_name": SOURCE_NAME,
        "drive_file_id": SOURCE_FILE_ID,
        "compressed_bytes": SOURCE_BYTES,
        "n_targets": EXPECTED_TARGETS,
        "download_performed": False,
        "science_measurement": False,
    }
    write_json(output_dir / "summary.json", summary)
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target-manifest", type=Path, default=DEFAULT_TARGET_MANIFEST)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    prepare(args.target_manifest, args.output_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
