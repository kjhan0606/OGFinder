#!/usr/bin/env python3
"""Validate the approved S23B CAMIRA member file and join target positions.

This is an offline importer. It never contacts HSC or Google Drive. It checks
the approved compressed file against its receipt, validates its eight-column
rows against the pinned S23B cluster catalog, joins the pinned 1,043-target
sample, and tests candidate member positions against each target FITS WCS.

The seventh and eighth source fields are preserved as raw values and are not
given scientific meanings by this script.
"""
from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import statistics
import sys

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS


DATA = Path("/gpfs/kjhan/Astrafek/data")
SAMPLE = DATA / "sample_800_2026-10-13"
IMPORT_DIR = SAMPLE / "s23b_camira_member_catalog_2026-11-07_01"
SOURCE_PART = IMPORT_DIR / "camira_s23b_wide_sm_member.dat.gz.part"
SOURCE = IMPORT_DIR / "camira_s23b_wide_sm_member.dat.gz"
CATALOG = DATA / "catalogue_s23b_camira_2026-10-20_01/camira_s23b_wide_sm.tbl"
TARGET_PLAN = SAMPLE / "s23b_icl_target_image_plan_2026-11-07_01"
TARGET_MANIFEST = TARGET_PLAN / "target_manifest.json"
TARGET_DOWNLOAD = TARGET_PLAN / "target_download_01"
TARGET_STATUS = TARGET_DOWNLOAD / "download_status.json"
TARGET_DOWNLOAD_MANIFEST = TARGET_DOWNLOAD / "download_manifest.json"
TARGET_IMAGES = TARGET_DOWNLOAD / "images"
TARGET_AUDIT = SAMPLE / "s23b_icl_target_image_audit_2026-11-07_01"
TARGET_AUDIT_SUMMARY = TARGET_AUDIT / "summary.json"
TARGET_VERIFIED_TABLE = TARGET_AUDIT / "verified_target_frames.tsv"
PLAN = SAMPLE / "s23b_camira_member_import_plan_2026-11-07_01/plan.json"
OUTPUT = IMPORT_DIR / "import_2026-11-07_02"

EXPECTED_SOURCE_BYTES = 10_626_590
EXPECTED_SOURCE_SHA256 = "66529e02c6ae11ebdccc3fb1d4756b639364f372f9d0f09a721786a19906945a"
EXPECTED_SOURCE_ROWS = 637_146
EXPECTED_SOURCE_CLUSTERS = 12_162
EXPECTED_TARGETS = 1_043
EXPECTED_TARGET_CATALOG_SHA256 = "5af057028feab7976cd3ed389a7246cbe73add8bbb96891101f0f692b1de4d71"
EXPECTED_TARGET_MANIFEST_SHA256 = "e487582ba40cd3e58da67704242db4d7219e381fb08f4b33fed47ce4ff268b44"
EXPECTED_FILE_ID = "1DEnjGj9P1ilUooEm2U5_A-qXUrxeiJdE"
RICHNESS_TOL = 0.00051
REDSHIFT_TOL = 0.00051
MAX_CENTER_MATCH_ARCSEC = 0.10


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_json(path: Path, payload: dict) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("x", encoding="utf-8") as stream:
        json.dump(payload, stream, ensure_ascii=False, indent=2, allow_nan=False)
        stream.write("\n")
        stream.flush()
        os.fsync(stream.fileno())
    temporary.replace(path)


def angular_separation_arcsec(ra1, dec1, ra2, dec2):
    """Vectorized great-circle separation in arcseconds."""
    ra1 = np.deg2rad(np.asarray(ra1, dtype=float))
    dec1 = np.deg2rad(np.asarray(dec1, dtype=float))
    ra2 = math.radians(float(ra2))
    dec2 = math.radians(float(dec2))
    dra = (ra1 - ra2 + math.pi) % (2 * math.pi) - math.pi
    a = (
        np.sin((dec1 - dec2) / 2.0) ** 2
        + math.cos(dec2) * np.cos(dec1) * np.sin(dra / 2.0) ** 2
    )
    return np.rad2deg(2.0 * np.arcsin(np.sqrt(np.clip(a, 0.0, 1.0)))) * 3600.0


def read_catalog(path: Path) -> tuple[dict, dict]:
    by_id = {}
    by_coord = {}
    with path.open(encoding="ascii") as stream:
        for line_number, line in enumerate(stream, 1):
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            fields = stripped.split()
            if len(fields) != 7:
                raise ValueError(f"cluster catalog line {line_number} has {len(fields)} fields")
            row = {
                "catalog_id": fields[0],
                "name": fields[1],
                "ra_deg": float(fields[2]),
                "dec_deg": float(fields[3]),
                "z_cl": float(fields[4]),
                "richness": float(fields[5]),
                "z_spec_bcg": float(fields[6]),
            }
            if not all(math.isfinite(row[k]) for k in (
                "ra_deg", "dec_deg", "z_cl", "richness", "z_spec_bcg"
            )):
                raise ValueError(f"nonfinite cluster catalog value at line {line_number}")
            if row["catalog_id"] in by_id:
                raise ValueError(f"duplicate cluster catalog ID {row['catalog_id']}")
            key = (round(row["ra_deg"], 6), round(row["dec_deg"], 6))
            by_id[row["catalog_id"]] = row
            by_coord.setdefault(key, []).append(row)
    return by_id, by_coord


def load_verified_frames() -> tuple[dict, dict]:
    audit = json.loads(TARGET_AUDIT_SUMMARY.read_text(encoding="utf-8"))
    status = json.loads(TARGET_STATUS.read_text(encoding="utf-8"))
    manifest_sha = sha256(TARGET_MANIFEST)
    status_sha = sha256(TARGET_STATUS)
    download_manifest_sha = sha256(TARGET_DOWNLOAD_MANIFEST)
    verified_sha = sha256(TARGET_VERIFIED_TABLE)
    if audit.get("status") != "complete_header_wcs_file_integrity_audit":
        raise ValueError("target-image audit is not complete")
    if audit.get("n_verified_frames") != EXPECTED_TARGETS:
        raise ValueError("target-image audit frame count differs from the pinned sample")
    expected_inputs = audit.get("inputs_sha256", {})
    for path, found in (
        (TARGET_MANIFEST, manifest_sha),
        (TARGET_STATUS, status_sha),
        (TARGET_DOWNLOAD_MANIFEST, download_manifest_sha),
    ):
        if expected_inputs.get(str(path)) != found:
            raise ValueError(f"target audit input hash changed: {path}")
    if audit.get("output_sha256", {}).get("verified_target_frames.tsv") != verified_sha:
        raise ValueError("verified target-frame table hash changed")
    if status.get("status") != "complete" or len(status.get("completed", [])) != EXPECTED_TARGETS:
        raise ValueError("target download status is not complete")

    rows = {}
    with TARGET_VERIFIED_TABLE.open(newline="", encoding="utf-8") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        for row in reader:
            if row["id"] in rows:
                raise ValueError(f"duplicate verified image ID {row['id']}")
            rows[row["id"]] = row
    completed = {str(item["id"]): item for item in status["completed"]}
    if len(rows) != EXPECTED_TARGETS or set(rows) != set(completed):
        raise ValueError("verified target-frame table and download status IDs differ")
    for cluster_id, frame in completed.items():
        row = rows[cluster_id]
        review = frame["review"]
        if (row["file"] != frame["file"]
                or int(row["bytes"]) != int(frame["bytes"])
                or row["sha256"] != frame["sha256"]
                or row["primary_filter"] != review["coadd_header_filter"]
                or int(row["shape_y"]) != int(review["shape"][0])
                or int(row["shape_x"]) != int(review["shape"][1])):
            raise ValueError(f"target frame table disagrees with download record: {cluster_id}")
        path = TARGET_IMAGES / frame["file"]
        if not path.is_file() or path.stat().st_size != int(frame["bytes"]):
            raise ValueError(f"target image missing or size changed: {cluster_id}")
    return rows, completed


def process(args) -> dict:
    source_input = args.source_part if args.source_part.exists() else args.source
    if source_input == args.source and not args.source.is_file():
        raise FileNotFoundError(args.source)
    if args.source.exists() and args.source_part.exists():
        raise ValueError("both staged and final member files exist; inspect before continuing")
    if source_input.stat().st_size != EXPECTED_SOURCE_BYTES:
        raise ValueError("compressed member file byte count differs from approved plan")
    source_sha = sha256(source_input)
    if source_sha != EXPECTED_SOURCE_SHA256:
        raise ValueError("compressed member file SHA-256 differs from received-file audit")
    for path in (args.catalog, args.target_manifest):
        if not path.is_file():
            raise FileNotFoundError(path)
    catalog_sha = sha256(args.catalog)
    target_manifest_sha = sha256(args.target_manifest)
    if catalog_sha != EXPECTED_TARGET_CATALOG_SHA256:
        raise ValueError("pinned S23B cluster catalog changed")
    if target_manifest_sha != EXPECTED_TARGET_MANIFEST_SHA256:
        raise ValueError("pinned S23B target manifest changed")

    plan = json.loads(args.plan.read_text(encoding="utf-8"))
    if (plan.get("scope", {}).get("drive_file_id") != EXPECTED_FILE_ID
            or plan.get("scope", {}).get("compressed_bytes") != EXPECTED_SOURCE_BYTES
            or plan.get("target_sample", {}).get("n_targets") != EXPECTED_TARGETS):
        raise ValueError("member import plan no longer matches the approved file scope")
    if args.output.exists() or args.output.with_name(args.output.name + ".tmp").exists():
        raise FileExistsError("refusing to overwrite member import output")

    targets_manifest = json.loads(args.target_manifest.read_text(encoding="utf-8"))
    targets = {}
    for target in targets_manifest.get("clusters", []):
        catalog_id = str(target["catalog_id"])
        if catalog_id in targets:
            raise ValueError(f"duplicate target catalog ID {catalog_id}")
        targets[catalog_id] = target
    if len(targets) != EXPECTED_TARGETS:
        raise ValueError("target manifest has an unexpected cluster count")

    catalog_by_id, catalog_by_coord = read_catalog(args.catalog)
    if len(catalog_by_id) != EXPECTED_SOURCE_CLUSTERS:
        raise ValueError("S23B catalog row count differs from member-file row audit")
    for catalog_id, target in targets.items():
        row = catalog_by_id.get(catalog_id)
        if row is None or row["name"] != target["name"]:
            raise ValueError(f"target ID/name not found in pinned catalog: {catalog_id}")
        sep = float(angular_separation_arcsec(
            [target["ra_deg"]], [target["dec_deg"]], row["ra_deg"], row["dec_deg"]
        )[0])
        if (sep > MAX_CENTER_MATCH_ARCSEC
                or abs(float(target["z_cl"]) - row["z_cl"]) > REDSHIFT_TOL
                or abs(float(target["richness"]) - row["richness"]) > 0.0011):
            raise ValueError(f"target manifest properties disagree with catalog: {catalog_id}")

    target_ids = set(targets)
    signature_by_id = {}
    group_counts = {}
    target_rows = {catalog_id: [] for catalog_id in targets}
    n_source_rows = 0
    n_bad_width = 0
    n_nonfinite = 0
    n_bad_member_coordinates = 0
    unmatched_signatures = []

    with gzip.open(source_input, "rt", encoding="ascii", newline="") as stream:
        for source_row, line in enumerate(stream, 1):
            fields = line.split()
            if len(fields) != 8:
                n_bad_width += 1
                continue
            try:
                values = [float(value) for value in fields]
            except ValueError:
                n_nonfinite += 1
                continue
            if not all(math.isfinite(value) for value in values):
                n_nonfinite += 1
                continue
            ra_cluster, dec_cluster, richness, z_cl, ra_member, dec_member = values[:6]
            if not (0.0 <= ra_cluster < 360.0 and -90.0 <= dec_cluster <= 90.0
                    and 0.0 <= ra_member < 360.0 and -90.0 <= dec_member <= 90.0):
                n_bad_member_coordinates += 1
                continue
            key = (round(ra_cluster, 6), round(dec_cluster, 6))
            candidates = catalog_by_coord.get(key, [])
            matches = [
                row for row in candidates
                if abs(row["richness"] - richness) <= RICHNESS_TOL
                and abs(row["z_cl"] - z_cl) <= REDSHIFT_TOL
            ]
            if len(matches) != 1:
                if len(unmatched_signatures) < 8:
                    unmatched_signatures.append({
                        "source_row": source_row,
                        "ra_cluster_deg": ra_cluster,
                        "dec_cluster_deg": dec_cluster,
                        "richness": richness,
                        "z_cl": z_cl,
                        "catalog_candidates": len(candidates),
                        "matching_catalog_rows": len(matches),
                    })
                continue
            row = matches[0]
            catalog_id = row["catalog_id"]
            signature = (fields[0], fields[1], fields[2], fields[3])
            previous = signature_by_id.setdefault(catalog_id, signature)
            if previous != signature:
                raise ValueError(f"inconsistent cluster signature in member rows: {catalog_id}")
            group_counts[catalog_id] = group_counts.get(catalog_id, 0) + 1
            if catalog_id in target_ids:
                target_rows[catalog_id].append({
                    "source_row": source_row,
                    "member_ra_deg": ra_member,
                    "member_dec_deg": dec_member,
                    "field7_raw": fields[6],
                    "field8_raw": fields[7],
                })
            n_source_rows += 1

    if n_bad_width or n_nonfinite or n_bad_member_coordinates or unmatched_signatures:
        raise ValueError(
            "member rows failed structural/catalog validation: "
            f"bad_width={n_bad_width}, nonfinite={n_nonfinite}, "
            f"coordinate_bounds={n_bad_member_coordinates}, "
            f"unmatched_examples={unmatched_signatures}"
        )
    if n_source_rows != EXPECTED_SOURCE_ROWS:
        raise ValueError(f"member row count {n_source_rows} differs from receipt audit")
    if len(signature_by_id) != EXPECTED_SOURCE_CLUSTERS:
        raise ValueError("not every distinct member cluster maps uniquely to the S23B catalog")
    missing_targets = sorted(catalog_id for catalog_id in targets if catalog_id not in group_counts)
    if missing_targets:
        raise ValueError(f"target clusters missing member rows: {missing_targets[:12]}")

    verified_rows, completed = load_verified_frames()
    flattened = []
    counts = []
    all_offsets = []
    total_inside = 0
    total_positions = 0
    for catalog_id, target in targets.items():
        cluster_id = str(target["id"])
        frame = completed.get(cluster_id)
        if frame is None:
            raise ValueError(f"no audited target image for {cluster_id}")
        frame_row = verified_rows[cluster_id]
        image_path = TARGET_IMAGES / frame["file"]
        with fits.open(image_path, mode="readonly", memmap=True, lazy_load_hdus=True) as hdus:
            hdu_index = int(frame["review"]["science_hdu"])
            header = hdus[hdu_index].header
            nx, ny = int(header["NAXIS1"]), int(header["NAXIS2"])
            wcs = WCS(header).celestial
            if not wcs.has_celestial:
                raise ValueError(f"target science HDU has no celestial WCS: {cluster_id}")
            rows = target_rows[catalog_id]
            ra = np.asarray([r["member_ra_deg"] for r in rows], dtype=float)
            dec = np.asarray([r["member_dec_deg"] for r in rows], dtype=float)
            x, y = wcs.world_to_pixel_values(ra, dec)
            inside = (
                np.isfinite(x) & np.isfinite(y)
                & (x >= 0.0) & (x <= nx - 1.0)
                & (y >= 0.0) & (y <= ny - 1.0)
            )
            offsets = angular_separation_arcsec(
                ra, dec, float(target["ra_deg"]), float(target["dec_deg"])
            )
            all_offsets.extend(offsets.tolist())
            total_inside += int(np.count_nonzero(inside))
            total_positions += len(rows)
            counts.append({
                "target_cluster_id": cluster_id,
                "catalog_id": catalog_id,
                "cluster_name": target["name"],
                "n_member_rows": len(rows),
                "n_inside_image_pixel_centers": int(np.count_nonzero(inside)),
                "n_outside_image": int(len(rows) - np.count_nonzero(inside)),
                "image_footprint_fraction": float(np.count_nonzero(inside) / len(rows)),
            })
            for index, member in enumerate(rows):
                flattened.append({
                    "target_cluster_id": cluster_id,
                    "catalog_id": catalog_id,
                    "cluster_name": target["name"],
                    "source_row": member["source_row"],
                    "cluster_ra_deg": float(target["ra_deg"]),
                    "cluster_dec_deg": float(target["dec_deg"]),
                    "z_cl": float(target["z_cl"]),
                    "richness": float(target["richness"]),
                    "member_ra_deg": member["member_ra_deg"],
                    "member_dec_deg": member["member_dec_deg"],
                    "cluster_offset_arcsec": float(offsets[index]),
                    "image_x_pixel": float(x[index]),
                    "image_y_pixel": float(y[index]),
                    "inside_image_pixel_centers": bool(inside[index]),
                    "member_field7_raw": member["field7_raw"],
                    "member_field8_raw": member["field8_raw"],
                    "image_file": frame["file"],
                    "image_sha256_from_pinned_audit": frame["sha256"],
                    "image_file_bytes_from_pinned_audit": int(frame_row["bytes"]),
                })

    if total_positions != 45_634:
        raise ValueError(f"target-member row count differs from structural audit: {total_positions}")
    if not counts or len(counts) != EXPECTED_TARGETS:
        raise ValueError("target/member image-footprint join is incomplete")

    output_tmp = args.output.with_name(args.output.name + ".tmp")
    output_tmp.mkdir(parents=False, exist_ok=False)
    members_path = output_tmp / "target_member_positions.tsv"
    counts_path = output_tmp / "target_member_footprint_counts.tsv"
    member_fields = [
        "target_cluster_id", "catalog_id", "cluster_name", "source_row",
        "cluster_ra_deg", "cluster_dec_deg", "z_cl", "richness",
        "member_ra_deg", "member_dec_deg", "cluster_offset_arcsec",
        "image_x_pixel", "image_y_pixel", "inside_image_pixel_centers",
        "member_field7_raw", "member_field8_raw", "image_file",
        "image_sha256_from_pinned_audit", "image_file_bytes_from_pinned_audit",
    ]
    with members_path.open("x", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=member_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(flattened)
    count_fields = [
        "target_cluster_id", "catalog_id", "cluster_name", "n_member_rows",
        "n_inside_image_pixel_centers", "n_outside_image", "image_footprint_fraction",
    ]
    with counts_path.open("x", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=count_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(sorted(counts, key=lambda row: row["target_cluster_id"]))

    target_counts = sorted(row["n_member_rows"] for row in counts)
    inside_fraction = total_inside / total_positions
    offsets = np.asarray(all_offsets, dtype=float)
    summary = {
        "schema": "astrafek-s23b-camira-member-import-v1",
        "status": "complete_catalog_join_and_image_footprint_audit",
        "science_measurement": False,
        "ready_for_exploratory_satellite_analysis": True,
        "ready_for_science": False,
        "download_scope": {
            "approved_file_id": EXPECTED_FILE_ID,
            "approved_file_name": args.source.name,
            "approved_compressed_bytes": EXPECTED_SOURCE_BYTES,
            "other_files_downloaded": False,
            "credentials_saved": False,
        },
        "member_file": {
            "path": str(args.source),
            "compressed_bytes": EXPECTED_SOURCE_BYTES,
            "sha256": source_sha,
            "gzip_crc_and_stream_validated": True,
            "rows": n_source_rows,
            "distinct_cluster_groups": len(signature_by_id),
            "column_count": 8,
            "header_present": False,
            "interpreted_columns": {
                "1": "cluster RA [deg], matched to pinned cluster catalog",
                "2": "cluster Dec [deg], matched to pinned cluster catalog",
                "3": "cluster richness N_mem, matched within rounding tolerance",
                "4": "cluster photometric redshift z_cl, matched within rounding tolerance",
                "5": "member RA [deg], coordinate range and image-WCS checked",
                "6": "member Dec [deg], coordinate range and image-WCS checked",
                "7": "opaque numeric value preserved; meaning not assigned",
                "8": "opaque numeric value preserved; meaning not assigned",
            },
        },
        "catalog_join": {
            "cluster_catalog": str(args.catalog),
            "cluster_catalog_sha256": catalog_sha,
            "catalog_cluster_rows": len(catalog_by_id),
            "member_cluster_groups_matched": len(signature_by_id),
            "unmatched_or_ambiguous_cluster_groups": 0,
            "target_manifest": str(args.target_manifest),
            "target_manifest_sha256": target_manifest_sha,
            "target_clusters": EXPECTED_TARGETS,
            "target_clusters_with_member_rows": len(counts),
            "target_member_rows": total_positions,
            "target_member_rows_min_median_max": [
                min(target_counts), int(statistics.median(target_counts)), max(target_counts)
            ],
        },
        "image_footprint_audit": {
            "target_image_audit_summary": str(TARGET_AUDIT_SUMMARY),
            "target_image_audit_status": "complete_header_wcs_file_integrity_audit",
            "target_images_checked_by_wcs": len(counts),
            "target_member_positions": total_positions,
            "member_positions_inside_pixel_center_bounds": total_inside,
            "member_positions_outside_image": total_positions - total_inside,
            "overall_inside_fraction": inside_fraction,
            "cluster_offset_arcsec_q50_q95_max": [
                float(np.percentile(offsets, 50)),
                float(np.percentile(offsets, 95)),
                float(np.max(offsets)),
            ],
            "limitations": [
                "A position inside the FITS footprint is not proof of an unmasked or usable image pixel.",
                "This footprint audit does not assess membership purity, detection completeness, photometry, or BCG identity.",
                "The catalog's two headerless trailing fields are preserved without assigning semantics.",
            ],
        },
        "inputs_sha256": {
            str(args.source): source_sha,
            str(args.catalog): catalog_sha,
            str(args.target_manifest): target_manifest_sha,
            str(TARGET_AUDIT_SUMMARY): sha256(TARGET_AUDIT_SUMMARY),
            str(TARGET_VERIFIED_TABLE): sha256(TARGET_VERIFIED_TABLE),
            str(TARGET_STATUS): sha256(TARGET_STATUS),
        },
        "output_sha256": {
            members_path.name: sha256(members_path),
            counts_path.name: sha256(counts_path),
        },
        "import_script_sha256": sha256(Path(__file__).resolve()),
        "created_utc": datetime.now(timezone.utc).isoformat(),
    }
    atomic_json(output_tmp / "summary.json", summary)

    if source_input == args.source_part:
        os.replace(args.source_part, args.source)
    os.replace(output_tmp, args.output)

    plan["status"] = "complete_member_catalog_received_and_validated"
    plan["download_performed"] = True
    plan["ready_for_analysis"] = True
    plan["science_measurement"] = False
    plan["download_receipt"] = {
        "approved": True,
        "approved_scope": "the single S23B CAMIRA Wide star-mask member file",
        "drive_file_id": EXPECTED_FILE_ID,
        "bytes": EXPECTED_SOURCE_BYTES,
        "sha256": source_sha,
        "received_utc": summary["created_utc"],
        "credentials_saved": False,
        "other_files_downloaded": False,
    }
    plan["import_output"] = str(args.output)
    plan["import_script_sha256"] = sha256(Path(__file__).resolve())
    plan["import_summary_sha256"] = sha256(args.output / "summary.json")
    atomic_json(args.plan, plan)
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-part", type=Path, default=SOURCE_PART)
    parser.add_argument("--source", type=Path, default=SOURCE)
    parser.add_argument("--catalog", type=Path, default=CATALOG)
    parser.add_argument("--target-manifest", type=Path, default=TARGET_MANIFEST)
    parser.add_argument("--plan", type=Path, default=PLAN)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    summary = process(args)
    print(json.dumps(summary, ensure_ascii=False, indent=2, allow_nan=False))
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception as error:
        print(f"import_s23b_camira_member_catalog: {type(error).__name__}: {error}", file=sys.stderr)
        sys.exit(1)
