#!/usr/bin/env python3
"""Prepare an unsubmitted, position-specific HSC PSF plan for S23B members.

The plan is built only from the already received target coadds, member table,
and the saved 427-frame source-aware injection-recovery run. It does not
contact HSC or download anything. The exact coordinate list and provenance
hashes are written for explicit user review before any PSF request.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import sys
from collections import Counter, defaultdict

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from download_pilot import validate_fits

DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
RECOVERY = DATA / "s23b_camira_member_outer_injection_recovery_2026-11-07_01"
PROFILE = DATA / "s23b_camira_member_outer_profiles_full_2026-11-07_01"
MEMBER_IMPORT = DATA / "s23b_camira_member_catalog_2026-11-07_01/import_2026-11-07_02"
TARGET_PLAN = DATA / "s23b_icl_target_image_plan_2026-11-07_01"
TARGET_MASKS = DATA / "s23b_icl_mask_support_2026-11-07_01/target"
OLD_PSF_PLAN_DIRS = [
    DATA / "icl_psf_picker_injection_plan_2026-10-19_02",
    DATA / "icl_psf_picker_target_control_centers_plan_2026-10-21_01",
    DATA / "icl_psf_picker_satellite_support224_plan_2026-11-05_01",
]
OLD_PSF_DOWNLOAD_DIRS = [
    DATA / "icl_psf_picker_injection_2026-10-19_01",
    DATA / "icl_psf_picker_target_control_centers_download_2026-10-21_01",
    DATA / "icl_psf_picker_satellite_support224_download_2026-11-06_01",
]
DEFAULT_OUTPUT = DATA / "s23b_camira_member_profile_psf_plan_2026-11-07_02"
EXPECTED_RECOVERY_SUMMARY_SHA256 = "e9d4cced59f0b0f8ea3f189989886f28af74a339cc5a9b3421f8bca187defc77"
EXPECTED_MEMBER_TABLE_SHA256 = "d1eb604922146c4d953ba4cad2f92d4c515fb4d4f10cd553f62b683535bfa6ae"
EXPECTED_PROFILE_BINS_SHA256 = "56ff7feaaf2c2a7dd138c62f919df3d631d7798c9ff44cdc9482564045b963e1"
EXPECTED_MAX_PSF_BYTES = 20160

OLD_PSF_SOURCES = [
    (OLD_PSF_PLAN_DIRS[0], "psf_picker_plan.json", OLD_PSF_DOWNLOAD_DIRS[0],
     "psf_picker_audit_verified_2026-10-20.json", 744),
    (OLD_PSF_PLAN_DIRS[1], "psf_picker_plan.json", OLD_PSF_DOWNLOAD_DIRS[1],
     "psf_picker_audit_reproducible.json", 1542),
    (OLD_PSF_PLAN_DIRS[2], "psf_request_plan.json", OLD_PSF_DOWNLOAD_DIRS[2],
     "psf_picker_audit.json", 2005),
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def angular_sep_arcsec(ra1: float, dec1: float, ra2: float, dec2: float) -> float:
    a1, d1, a2, d2 = map(math.radians, (ra1, dec1, ra2, dec2))
    hav = math.sin((d2 - d1) / 2) ** 2 + math.cos(d1) * math.cos(d2) * math.sin((a2 - a1) / 2) ** 2
    return math.degrees(2 * math.asin(min(1.0, math.sqrt(hav)))) * 3600.0


def exact_psf_positions() -> tuple[set[tuple[float, float]], int, int, dict[str, str]]:
    positions: set[tuple[float, float]] = set()
    n_files = 0
    n_bytes = 0
    audit_hashes = {}
    for plan_dir, plan_name, download_dir, audit_name, expected_count in OLD_PSF_SOURCES:
        mapping_path = plan_dir / "position_mapping.json"
        plan_path = plan_dir / plan_name
        audit_path = download_dir / audit_name
        if not all(path.is_file() for path in (mapping_path, plan_path, audit_path)):
            raise FileNotFoundError("existing PSF bank is missing its plan, mapping, or audit")
        plan = read_json(plan_path)
        audit = read_json(audit_path)
        mapping = read_json(mapping_path)
        files = list(download_dir.glob("*.fits"))
        total_bytes = sum(path.stat().st_size for path in files)
        if (len(mapping) != expected_count
                or len(files) != expected_count
                or audit.get("status") != "verified"
                or audit.get("n_returned") != expected_count
                or audit.get("n_matched_coordinates") != expected_count
                or audit.get("total_bytes") != total_bytes
                or audit.get("position_mapping_sha256") != sha256(mapping_path)
                or plan.get("position_mapping_sha256") != sha256(mapping_path)
                or audit.get("plan_sha256") not in (None, sha256(plan_path))):
            raise ValueError(f"an existing PSF bank failed its plan/file/audit check: {plan_dir}")
        positions.update((round(float(row["ra_deg"]), 9), round(float(row["dec_deg"]), 9))
                         for row in mapping)
        n_files += len(files)
        n_bytes += total_bytes
        audit_hashes[str(audit_path.resolve())] = sha256(audit_path)
    return positions, n_files, n_bytes, audit_hashes


def prepare(output: Path) -> dict:
    output = output.resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite member PSF plan: {output}")

    recovery_summary_path = RECOVERY / "summary.json"
    selected_path = RECOVERY / "selected_members.tsv"
    frame_results_path = RECOVERY / "frame_results.jsonl"
    member_table_path = MEMBER_IMPORT / "target_member_positions.tsv"
    member_summary_path = MEMBER_IMPORT / "summary.json"
    profile_bins_path = PROFILE / "member_profile_bins.tsv"
    target_manifest_path = TARGET_PLAN / "target_manifest.json"
    target_plan_summary_path = TARGET_PLAN / "summary.json"
    mask_records_path = TARGET_MASKS / "records.tsv"
    required = [recovery_summary_path, selected_path, frame_results_path,
                member_table_path, member_summary_path, profile_bins_path,
                target_manifest_path, target_plan_summary_path, mask_records_path]
    for path in required:
        if not path.is_file():
            raise FileNotFoundError(path)

    recovery = read_json(recovery_summary_path)
    if (sha256(recovery_summary_path) != EXPECTED_RECOVERY_SUMMARY_SHA256
            or recovery.get("status") != "complete_exploratory_member_halo_injection_recovery"
            or recovery.get("n_frames_complete") != 427
            or recovery.get("n_frames_failed") != 0
            or recovery.get("n_baseline_profile_replays") != 427
            or recovery.get("science_measurement") is not False
            or recovery.get("images_ready_for_science") is not False):
        raise ValueError("source-aware injection-recovery receipt differs from the reviewed run")
    expected_outputs = recovery.get("output_sha256", {})
    for name, path in (("selected_members.tsv", selected_path),
                       ("frame_results.jsonl", frame_results_path)):
        if expected_outputs.get(name) != sha256(path):
            raise ValueError(f"injection-recovery {name} checksum mismatch")

    member_summary = read_json(member_summary_path)
    if (sha256(member_table_path) != EXPECTED_MEMBER_TABLE_SHA256
            or member_summary.get("status") != "complete_catalog_join_and_image_footprint_audit"
            or member_summary.get("science_measurement") is not False):
        raise ValueError("approved member-position table does not match its imported receipt")
    if sha256(profile_bins_path) != EXPECTED_PROFILE_BINS_SHA256:
        raise ValueError("saved full-profile bin table hash changed")
    target_plan_summary = read_json(target_plan_summary_path)
    if target_plan_summary.get("science_measurement") is not False:
        raise ValueError("target-image plan has an unexpected science state")
    if target_plan_summary.get("outputs_sha256", {}).get("target_manifest.json") != sha256(target_manifest_path):
        raise ValueError("target manifest no longer matches its pinned plan")

    selected = read_tsv(selected_path)
    if len(selected) != 427:
        raise ValueError(f"expected one selected member in each of 427 frames, got {len(selected)}")
    selected_by_frame = {(row["target_id"], row["filter"]): row for row in selected}
    if len(selected_by_frame) != len(selected):
        raise ValueError("duplicate target/filter entries in selected-member table")
    if Counter(row["filter"] for row in selected) != Counter({"i": 62, "i2": 365}):
        raise ValueError("selected-member filter split changed")

    frame_results = [json.loads(line) for line in frame_results_path.read_text(encoding="utf-8").splitlines() if line]
    result_by_frame = {(row["target_id"], row["filter"]): row for row in frame_results}
    if len(frame_results) != 427 or len(result_by_frame) != 427:
        raise ValueError("saved frame receipts do not contain exactly 427 unique frames")
    if any(row.get("status") != "complete" for row in frame_results):
        raise ValueError("at least one source-aware frame receipt is incomplete")

    manifest = read_json(target_manifest_path)
    frame_by_id = {str(row["id"]): row for row in manifest.get("clusters", [])}
    mask_rows = read_tsv(mask_records_path)
    mask_by_id = {row["id"]: row for row in mask_rows}
    if len(frame_by_id) != len(manifest.get("clusters", [])) or len(mask_by_id) != len(mask_rows):
        raise ValueError("duplicate target-frame identifiers in manifest or mask records")

    member_rows = read_tsv(member_table_path)
    member_by_key: dict[tuple[str, str], list[dict[str, str]]] = defaultdict(list)
    for row in member_rows:
        member_by_key[(row["target_cluster_id"], row["source_row"])].append(row)

    input_sha256 = {str(path.resolve()): sha256(path) for path in required}
    if input_sha256[str(member_table_path.resolve())] != EXPECTED_MEMBER_TABLE_SHA256:
        raise ValueError("member input digest differs from expected imported table")

    position_rows = []
    filter_counts: Counter[str] = Counter()
    for key in sorted(selected_by_frame):
        row = selected_by_frame[key]
        result = result_by_frame.get(key)
        if result is None:
            raise ValueError(f"selected member has no complete frame receipt: {key}")
        target_id, filt = key
        frame = frame_by_id.get(target_id)
        mask = mask_by_id.get(target_id)
        if frame is None or mask is None:
            raise ValueError(f"selected member frame missing from manifest/mask table: {target_id}")
        if mask.get("primary_filter", "").lower() != filt:
            raise ValueError(f"selected member and mask filter differ for {target_id}")
        image_path = Path(mask["image"])
        if (not image_path.is_file()
                or image_path.stat().st_size != int(mask["image_bytes"])
                or sha256(image_path) != result.get("image_sha256")
                or result.get("image_sha256") != mask.get("image_sha256")):
            raise ValueError(f"target coadd hash differs from the saved analysis receipt: {target_id}")
        info = validate_fits(image_path, frame)
        actual_filter = str(info["coadd_header_filter"]).lower()
        if actual_filter != filt:
            raise ValueError(f"FITS header filter differs from selected profile filter for {target_id}")

        member_key = (target_id, row["source_row"])
        candidates = member_by_key.get(member_key, [])
        if len(candidates) != 1:
            raise ValueError(f"source row is not unique in imported member table: {member_key}")
        source = candidates[0]
        ra, dec = float(row["member_ra_deg"]), float(row["member_dec_deg"])
        if (abs(float(source["member_ra_deg"]) - ra) > 1e-8
                or abs(float(source["member_dec_deg"]) - dec) > 1e-8
                or source.get("inside_image_pixel_centers") != "True"):
            raise ValueError(f"selected injection coordinate differs from approved member row: {member_key}")
        if not (math.isfinite(ra) and 0 <= ra < 360 and math.isfinite(dec) and -90 <= dec <= 90):
            raise ValueError(f"invalid member coordinate: {member_key}")

        position_rows.append({
            "cluster_id": target_id,
            "frame_filter": filt,
            "request_filter": "HSC-I",
            "role": "selected_camira_member_injection_center",
            "member_index": int(row["member_index"]),
            "source_row": row["source_row"],
            "ra_deg": round(ra, 9),
            "dec_deg": round(dec, 9),
            "cluster_radius_kpc": float(row["cluster_radius_kpc"]),
            "image_path": str(image_path.resolve()),
            "image_sha256": result["image_sha256"],
            "image_primary_filter": actual_filter,
        })
        filter_counts[actual_filter] += 1

    unique_by_coord: dict[tuple[float, float], dict] = {}
    for row in position_rows:
        coord = (round(row["ra_deg"], 9), round(row["dec_deg"], 9))
        if coord in unique_by_coord:
            unique_by_coord[coord]["roles"].append({
                "cluster_id": row["cluster_id"], "frame_filter": row["frame_filter"],
                "member_index": row["member_index"], "source_row": row["source_row"],
            })
        else:
            unique_by_coord[coord] = {**row, "roles": [{
                "cluster_id": row["cluster_id"], "frame_filter": row["frame_filter"],
                "member_index": row["member_index"], "source_row": row["source_row"],
            }]}
    requests = sorted(unique_by_coord.values(), key=lambda row: (row["ra_deg"], row["dec_deg"]))

    previous_psf_positions, previous_psf_files, previous_psf_bytes, previous_psf_audits = exact_psf_positions()
    exact_reuse = [row for row in requests
                   if (round(row["ra_deg"], 9), round(row["dec_deg"], 9)) in previous_psf_positions]
    if exact_reuse:
        raise ValueError("current member positions unexpectedly overlap a previously received PSF position")
    nearest_distances = [
        min(angular_sep_arcsec(row["ra_deg"], row["dec_deg"], ra, dec)
            for ra, dec in previous_psf_positions)
        for row in requests
    ]
    nearest_previous_arcsec = min(nearest_distances)

    output.mkdir(parents=True, exist_ok=False)
    coordinate_path = output / "psf_picker_coordinates.txt"
    with coordinate_path.open("x", encoding="ascii", newline="") as stream:
        stream.write("#? ra dec\n")
        for row in requests:
            stream.write(f"{row['ra_deg']:.9f} {row['dec_deg']:+.9f}\n")
    mapping_path = output / "position_mapping.json"
    mapping_path.write_text(json.dumps(requests, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    estimated_bytes = len(requests) * EXPECTED_MAX_PSF_BYTES
    plan = {
        "schema": "astrafek-s23b-camira-member-profile-psf-plan-v1",
        "status": "prepared_not_submitted",
        "purpose": "position-specific HSC coadd model-PSF cores for the source-aware CAMIRA-member profile injection/recovery",
        "request": {
            "rerun": "pdr3_wide", "type": "coadd", "filter": "HSC-I", "centered": True,
            "product": "position-specific coadd model PSF FITS core",
        },
        "n_member_profile_frames": len(position_rows),
        "n_unique_member_positions": len(requests),
        "n_exact_positions_covered_by_existing_psf_banks": len(exact_reuse),
        "n_existing_psf_files_in_audited_banks": previous_psf_files,
        "existing_psf_bank_total_bytes": previous_psf_bytes,
        "nearest_member_to_existing_psf_position_arcsec": nearest_previous_arcsec,
        "frame_filter_counts": dict(sorted(filter_counts.items())),
        "expected_fits_count": len(requests),
        "observed_max_existing_fits_bytes": EXPECTED_MAX_PSF_BYTES,
        "estimated_bytes_if_each_file_matches_largest_audited_existing_psf": estimated_bytes,
        "estimated_mib": estimated_bytes / (1024 ** 2),
        "coordinate_list": str(coordinate_path),
        "coordinate_list_sha256": sha256(coordinate_path),
        "position_mapping": str(mapping_path),
        "position_mapping_sha256": sha256(mapping_path),
        "analysis_script": str(Path(__file__).resolve()),
        "analysis_script_sha256": sha256(Path(__file__).resolve()),
        "inputs_sha256": input_sha256,
        "prior_psf_mapping_sha256": {
            str(path.resolve()): sha256(path / "position_mapping.json")
            for path in OLD_PSF_PLAN_DIRS
        },
        "prior_psf_audit_sha256": previous_psf_audits,
        "images_requested": False,
        "coadd_images_requested": False,
        "shape_rows_requested": False,
        "photoz_requested": False,
        "science_measurement": False,
        "psf_approved": False,
        "network_submitted": False,
        "requires_explicit_user_approval": True,
        "max_total_bytes_for_review": estimated_bytes,
    }
    plan_path = output / "psf_request_plan.json"
    plan_path.write_text(json.dumps(plan, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    result = {
        "schema": "astrafek-s23b-member-profile-psf-plan-preparation-v1",
        "status": "prepared_not_submitted",
        "plan_path": str(plan_path),
        "plan_sha256": sha256(plan_path),
        "n_member_profile_frames": len(position_rows),
        "n_unique_member_positions": len(requests),
        "n_exact_positions_covered_by_existing_psf_banks": len(exact_reuse),
        "nearest_member_to_existing_psf_position_arcsec": nearest_previous_arcsec,
        "n_member_positions_within_1arcsec_of_existing_psf": sum(distance < 1.0 for distance in nearest_distances),
        "n_member_positions_within_60arcsec_of_existing_psf": sum(distance < 60.0 for distance in nearest_distances),
        "frame_filter_counts": dict(sorted(filter_counts.items())),
        "estimated_bytes_if_each_file_matches_largest_audited_existing_psf": estimated_bytes,
        "max_total_bytes_for_review": estimated_bytes,
        "coordinate_list_sha256": sha256(coordinate_path),
        "position_mapping_sha256": sha256(mapping_path),
        "network_submitted": False,
        "science_measurement": False,
    }
    (output / "preparation_summary.json").write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    print(json.dumps(prepare(args.output), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
