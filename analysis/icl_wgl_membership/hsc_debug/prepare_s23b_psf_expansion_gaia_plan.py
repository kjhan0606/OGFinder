#!/usr/bin/env python3
"""Audit expanded local PSF candidates and prepare an unsubmitted Gaia plan.

This script reads only the already-received S23B target/control coadds and
local preparation products. It writes a content-pinned Gaia DR3 positional
query plan but never submits it.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import Counter
from pathlib import Path

import match_psf_gaia
from match_psf_gaia_batched import GAIA_FIELDS, build_query


DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
EXPANSION = DATA / "s23b_psf_expansion_2026-11-07_01"
PLAN_DIR = EXPANSION / "plan"
TARGET_PREP = EXPANSION / "target_preparation"
CONTROL_PREP = EXPANSION / "control_preparation"
DEFAULT_OUTPUT = EXPANSION / "gaia_plan"
BATCH_SIZE = 250
MAX_MATCH_ROWS = 20_000
MATCH_RADIUS_ARCSEC = 2.0


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, payload: dict) -> None:
    with path.open("x", encoding="utf-8") as stream:
        json.dump(payload, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")


def load_json(path: Path) -> dict:
    with path.open(encoding="utf-8") as stream:
        return json.load(stream)


def frame_candidates(population: str, manifest_path: Path, preparation_dir: Path,
                     prep_summary_path: Path, image_hashes: dict,
                     table_hashes: dict) -> tuple[list[dict], list[dict]]:
    manifest = load_json(manifest_path)
    prep_summary = load_json(prep_summary_path)
    if (prep_summary.get("manifest_sha256") != sha256(manifest_path)
            or prep_summary.get("science_measurement") is not False):
        raise ValueError(f"{population} preparation summary does not pin its manifest")
    if prep_summary.get("status") not in ("complete", "partial"):
        raise ValueError(f"{population} preparation has an unexpected status")
    manifest_frames = {row["id"]: row for row in manifest.get("clusters", [])}
    prepared_records = {row["id"]: row for row in prep_summary.get("records", [])}
    errors = {row["id"]: row for row in prep_summary.get("errors", [])}
    if set(prepared_records) & set(errors):
        raise ValueError(f"{population} frame is both prepared and failed")
    if set(prepared_records) | set(errors) != set(manifest_frames):
        raise ValueError(f"{population} preparation does not account for every selected frame")

    requests = []
    audit_rows = []
    for frame_id, frame in manifest_frames.items():
        prep = prepared_records.get(frame_id)
        if prep is None:
            failure = errors[frame_id]
            audit_rows.append({
                "population": population,
                "frame_id": frame_id,
                "target_cluster_id": frame.get("target_cluster_id", frame_id),
                "filter": frame.get("coadd_header_filter", ""),
                "status": "frame_preparation_failed",
                "error_type": failure.get("error_type", ""),
                "error": failure.get("error", ""),
                "n_detections": 0,
                "n_candidates": 0,
                "image_sha256": frame.get("input_image_sha256", ""),
            })
            continue
        if prep.get("source_image") != frame.get("image"):
            raise ValueError(f"{population} preparation image path changed for {frame_id}")
        if prep.get("center_source_identification_enabled") is not False:
            raise ValueError(f"center was treated as a BCG in {population} frame {frame_id}")
        image_path = Path(frame["image"])
        expected_image_sha = frame.get("input_image_sha256")
        if (not expected_image_sha or image_path.stat().st_size != int(frame["input_image_bytes"])
                or sha256(image_path) != expected_image_sha):
            raise ValueError(f"{population} input image hash differs for {frame_id}")
        image_hashes[str(image_path.resolve())] = expected_image_sha

        frame_dir = preparation_dir / f"cluster_{frame_id}"
        candidates_path = frame_dir / "psf_candidates.tsv"
        detections_path = frame_dir / "detections.tsv"
        prep_record_path = frame_dir / "preparation.json"
        for path in (candidates_path, detections_path, prep_record_path):
            if not path.is_file():
                raise FileNotFoundError(path)
            table_hashes[str(path.resolve())] = sha256(path)
        local_prep = load_json(prep_record_path)
        if (local_prep.get("id") != frame_id
                or local_prep.get("science_measurement") is not False
                or local_prep.get("center_source_identification_enabled") is not False):
            raise ValueError(f"per-frame preparation receipt mismatch for {frame_id}")

        candidate_rows = read_tsv(candidates_path)
        detection_rows = read_tsv(detections_path)
        detections_by_number = {}
        for row in detection_rows:
            number = int(row["NUMBER"])
            if number in detections_by_number:
                raise ValueError(f"duplicate detection number in {frame_id}: {number}")
            detections_by_number[number] = row
        if len(candidate_rows) != int(prep["n_psf_candidates"]):
            raise ValueError(f"candidate table row count differs from summary for {frame_id}")
        if len(candidate_rows) > 25:
            raise ValueError(f"candidate cap exceeded for {frame_id}")
        seen_numbers = set()
        filt = str(frame.get("coadd_header_filter", "")).lower()
        if filt not in ("i", "i2"):
            raise ValueError(f"unrecognized primary HSC-I filter for {frame_id}")
        for candidate in candidate_rows:
            number = int(candidate["NUMBER"])
            if number in seen_numbers:
                raise ValueError(f"duplicate PSF candidate number in {frame_id}: {number}")
            seen_numbers.add(number)
            detection = detections_by_number.get(number)
            if detection is None:
                raise ValueError(f"PSF candidate lacks matching detection row in {frame_id}")
            if candidate.get("CLASSIFICATION") != "provisional_compact_candidate":
                raise ValueError(f"unexpected morphology class in {frame_id}")
            ra, dec = float(detection["RA"]), float(detection["DEC"])
            if not (math.isfinite(ra) and 0 <= ra < 360 and math.isfinite(dec) and -90 <= dec <= 90):
                raise ValueError(f"invalid candidate coordinate in {frame_id} / {number}")
            target_id = frame.get("target_cluster_id", frame_id)
            request = {
                "local_index": -1,
                "population": population,
                "cluster_id": frame_id,
                "target_cluster_id": target_id,
                "control_id": frame_id if population == "control" else frame.get("paired_control_id", ""),
                "number": number,
                "ra_deg": ra,
                "dec_deg": dec,
                "classification": candidate["CLASSIFICATION"],
                "primary_filter": filt,
                "fwhm_moment_pix": float(candidate["FWHM_MOMENT_PIX"]),
                "snr_nominal": float(candidate["SNR_NOMINAL"]),
                "wing_valid_fraction": float(candidate["WING_VALID_FRACTION"]),
                "image_sha256": expected_image_sha,
                "preparation_candidate_table_sha256": table_hashes[str(candidates_path.resolve())],
            }
            requests.append(request)
        audit_rows.append({
            "population": population,
            "frame_id": frame_id,
            "target_cluster_id": frame.get("target_cluster_id", frame_id),
            "filter": filt,
            "status": "prepared",
            "error_type": "",
            "error": "",
            "n_detections": int(prep["n_detections"]),
            "n_candidates": len(candidate_rows),
            "image_sha256": expected_image_sha,
        })
    return requests, audit_rows


def prepare(output: Path) -> dict:
    output = output.resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite Gaia expansion plan: {output}")
    inputs = [
        PLAN_DIR / "summary.json", PLAN_DIR / "selection.tsv",
        PLAN_DIR / "target_manifest.json", PLAN_DIR / "control_manifest.json",
        TARGET_PREP / "summary.json", CONTROL_PREP / "summary.json",
    ]
    for path in inputs:
        if not path.is_file():
            raise FileNotFoundError(path)
    plan_summary = load_json(PLAN_DIR / "summary.json")
    if (plan_summary.get("n_pairs_selected") != 200
            or plan_summary.get("n_pairs_by_filter") != {"i": 100, "i2": 100}
            or plan_summary.get("external_queries_submitted") is not False):
        raise ValueError("local matched-pair selection plan differs from reviewed scope")
    if sha256(PLAN_DIR / "target_manifest.json") != plan_summary.get("target_manifest_sha256"):
        raise ValueError("target expansion manifest hash differs from its plan")
    if sha256(PLAN_DIR / "control_manifest.json") != plan_summary.get("control_manifest_sha256"):
        raise ValueError("control expansion manifest hash differs from its plan")
    selection_rows = read_tsv(PLAN_DIR / "selection.tsv")
    if len(selection_rows) != 200:
        raise ValueError("matched-pair selection table must contain 200 pairs")

    image_hashes: dict[str, str] = {}
    table_hashes: dict[str, str] = {}
    target_requests, target_audit = frame_candidates(
        "target", PLAN_DIR / "target_manifest.json", TARGET_PREP,
        TARGET_PREP / "summary.json", image_hashes, table_hashes,
    )
    control_requests, control_audit = frame_candidates(
        "control", PLAN_DIR / "control_manifest.json", CONTROL_PREP,
        CONTROL_PREP / "summary.json", image_hashes, table_hashes,
    )
    requests = target_requests + control_requests
    if not requests:
        raise ValueError("no morphology candidates available")
    request_ids = set()
    for index, request in enumerate(requests):
        key = (request["population"], request["cluster_id"], request["number"])
        if key in request_ids:
            raise ValueError(f"duplicate candidate identity: {key}")
        request_ids.add(key)
        request["local_index"] = index

    by_filter_population = {
        f"{population}_{filt}": sum(
            row["population"] == population and row["primary_filter"] == filt for row in requests
        )
        for population in ("target", "control") for filt in ("i", "i2")
    }
    frames_by_filter_population = {
        f"{population}_{filt}": len({
            row["cluster_id"] for row in requests
            if row["population"] == population and row["primary_filter"] == filt
        })
        for population in ("target", "control") for filt in ("i", "i2")
    }
    existing_target_status_path = DATA / "s23b_psf_gaia_match_2026-11-07_01/status.json"
    existing_target_response_path = DATA / "s23b_psf_gaia_match_2026-11-07_01/gaia_matches.csv"
    existing_control_status_path = DATA / "s23b_psf_control_gaia_match_2026-11-07_01/status.json"
    existing_control_response_path = DATA / "s23b_psf_control_gaia_match_2026-11-07_01/gaia_matches.csv"
    for path in (existing_target_status_path, existing_target_response_path,
                 existing_control_status_path, existing_control_response_path):
        if not path.is_file():
            raise FileNotFoundError(path)
    target_status = load_json(existing_target_status_path)
    control_status = load_json(existing_control_status_path)
    if (target_status.get("status") != "complete" or control_status.get("status") != "complete"
            or sha256(existing_target_response_path) != target_status.get("response_sha256")
            or sha256(existing_control_response_path) != control_status.get("response_sha256")):
        raise ValueError("existing approved target/control Gaia query receipts are not complete")
    baseline_candidates = int(target_status["n_candidates"]) + int(control_status["n_candidates"])
    baseline_matches = int(target_status["n_match_rows"]) + int(control_status["n_match_rows"])
    expected_rows = round(len(requests) * baseline_matches / baseline_candidates)
    baseline_response_bytes = sum(
        (DATA / f"{name}/batch_{batch:03d}/gaia_cone_union.csv").stat().st_size
        for name, status in (
            ("s23b_psf_gaia_match_2026-11-07_01", target_status),
            ("s23b_psf_control_gaia_match_2026-11-07_01", control_status),
        )
        for batch in range(int(status["completed_batches"]))
    )
    expected_bytes = round(len(requests) * baseline_response_bytes / baseline_candidates)

    output.mkdir(parents=True, exist_ok=False)
    requests_path = output / "requests.json"
    write_json(requests_path, {"candidates": requests})
    positions_path = output / "candidate_positions.tsv"
    positions = [{
        "local_index": row["local_index"],
        "population": row["population"],
        "target_cluster_id": row["target_cluster_id"],
        "control_id": row["control_id"],
        "frame_id": row["cluster_id"],
        "filter": row["primary_filter"],
        "number": row["number"],
        "ra_deg": f"{row['ra_deg']:.12f}",
        "dec_deg": f"{row['dec_deg']:.12f}",
        "classification": row["classification"],
    } for row in requests]
    write_tsv(positions_path, positions, list(positions[0]))
    candidate_audit_path = output / "candidate_frame_audit.tsv"
    all_audit_rows = target_audit + control_audit
    write_tsv(candidate_audit_path, all_audit_rows, list(all_audit_rows[0]))

    query_dir = output / "query_batches"
    query_dir.mkdir()
    query_batches = []
    for batch_number, offset in enumerate(range(0, len(requests), BATCH_SIZE)):
        batch_rows = requests[offset:offset + BATCH_SIZE]
        query_path = query_dir / f"batch_{batch_number:03d}.adql"
        query_text = build_query(batch_rows) + "\n"
        with query_path.open("x", encoding="utf-8") as stream:
            stream.write(query_text)
        query_batches.append({
            "batch": batch_number,
            "first_candidate_id": offset,
            "n_candidates": len(batch_rows),
            "query_file": str(query_path.resolve()),
            "query_sha256": sha256(query_path),
        })

    query_manifest = json.dumps(query_batches, sort_keys=True, separators=(",", ":")).encode()
    input_hashes = {str(path.resolve()): sha256(path) for path in inputs}
    input_hashes.update({str(path.resolve()): sha256(path) for path in (
        existing_target_status_path, existing_target_response_path,
        existing_control_status_path, existing_control_response_path,
    )})
    plan = {
        "schema": "astrafek-s23b-matched-psf-expansion-gaia-plan-v1",
        "status": "prepared_not_submitted",
        "external_query": "Gaia DR3 TAP positional cross-match only",
        "endpoint": match_psf_gaia.ENDPOINT,
        "request_parameters": {
            "REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv", "MAXREC": str(MAX_MATCH_ROWS),
        },
        "match_radius_arcsec": MATCH_RADIUS_ARCSEC,
        "n_matched_target_control_pairs": 200,
        "n_candidate_positions": len(requests),
        "n_target_candidate_positions": len(target_requests),
        "n_control_candidate_positions": len(control_requests),
        "candidate_positions_by_population_filter": by_filter_population,
        "frames_with_candidates_by_population_filter": frames_by_filter_population,
        "candidate_frames_selected_by_filter": plan_summary["n_pairs_by_filter"],
        "batch_size_recommended": BATCH_SIZE,
        "expected_batches": len(query_batches),
        "expected_match_rows_rough_estimate": expected_rows,
        "estimated_response_bytes_rough": expected_bytes,
        "maxrec_per_batch": MAX_MATCH_ROWS,
        "max_match_rows_before_truncation_guard": MAX_MATCH_ROWS,
        "gaia_columns_queried": GAIA_FIELDS,
        "columns_returned": ["candidate_id", *GAIA_FIELDS, "separation_arcsec"],
        "query_method": (
            "Synchronous Gaia TAP cone-union ADQL; each batch embeds the approved candidate ICRS "
            "positions. Candidate indices and exact Gaia/HSC separations are associated locally."
        ),
        "images_requested": False,
        "hsc_shape_rows_requested": False,
        "photoz_requested": False,
        "cas_query_requested": False,
        "science_measurement": False,
        "psf_approved": False,
        "user_approval_required_before_submission": True,
        "candidate_list": str(requests_path.resolve()),
        "candidate_list_sha256": sha256(requests_path),
        "candidate_positions_tsv": str(positions_path.resolve()),
        "candidate_positions_tsv_sha256": sha256(positions_path),
        "candidate_frame_audit_tsv": str(candidate_audit_path.resolve()),
        "candidate_frame_audit_sha256": sha256(candidate_audit_path),
        "query_batches": query_batches,
        "query_batch_manifest_sha256": hashlib.sha256(query_manifest).hexdigest(),
        "input_sha256": input_hashes,
        "candidate_preparation_tables_sha256": table_hashes,
        "selected_hsc_image_sha256": image_hashes,
        "estimate_basis": (
            f"Scaled from two completed, separately approved Gaia DR3 positional samples: "
            f"{baseline_matches} matched rows from {baseline_candidates} candidate cones; "
            f"raw union-response bytes across those requests={baseline_response_bytes}."
        ),
        "notes": [
            "The current plan covers only the 200 newly selected, seeing-stratified matched pairs.",
            "Existing approvals for the 492 target and 458 control pilot positions do not cover this expansion.",
            "The exact candidate list, 2-arcsec search cones, endpoint, Gaia fields, and all SQL batch hashes are pinned locally.",
            "The remote request contains candidate coordinates only; local HSC morphology and member-coordinate data are not sent.",
            "No HSC images, shape rows, photo-z values, or CAS data are requested.",
            "Gaia non-matches are not evidence that a compact HSC candidate is non-stellar.",
        ],
    }
    plan_path = output / "plan.json"
    write_json(plan_path, plan)
    summary = {
        "schema": "astrafek-s23b-matched-psf-expansion-gaia-summary-v1",
        "status": "prepared_not_submitted",
        "plan_path": str(plan_path.resolve()),
        "plan_sha256": sha256(plan_path),
        "n_candidate_positions": len(requests),
        "n_target_candidate_positions": len(target_requests),
        "n_control_candidate_positions": len(control_requests),
        "candidate_positions_by_population_filter": by_filter_population,
        "frames_with_candidates_by_population_filter": frames_by_filter_population,
        "n_candidate_frames_with_failed_preparation": len([row for row in all_audit_rows if row["status"] == "frame_preparation_failed"]),
        "expected_batches": len(query_batches),
        "estimated_match_rows_rough": expected_rows,
        "estimated_response_bytes_rough": expected_bytes,
        "images_requested": False,
        "hsc_shape_rows_requested": False,
        "photoz_requested": False,
        "cas_query_requested": False,
        "science_measurement": False,
        "psf_approved": False,
        "query_submitted": False,
    }
    write_json(output / "summary.json", summary)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    result = prepare(args.output)
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
