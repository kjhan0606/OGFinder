#!/usr/bin/env python3
"""Prepare an exact, unsubmitted Gaia plan for matched-control PSF candidates.

The builder reads the completed local 49-pair control-candidate audit and
creates a hash-pinned Gaia DR3 positional-query plan. It never submits the
query. A user must separately approve the exact resulting plan SHA-256.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import match_psf_gaia
from match_psf_gaia_batched import GAIA_FIELDS, build_query


DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
AUDIT_DIR = DATA / "s23b_psf_control_candidate_audit_2026-11-07_01"
CONTROL_PILOT = DATA / "s23b_psf_control_candidates_2026-11-07_01"
TARGET_RUN = DATA / "s23b_psf_gaia_match_2026-11-07_01"
TARGET_REQUESTS = DATA / "s23b_psf_gaia_plan_2026-11-07_01/requests.json"
DEFAULT_OUTPUT = DATA / "s23b_psf_control_gaia_plan_2026-11-07_03"
BATCH_SIZE = 250
MAX_MATCH_ROWS = 20_000


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def write_json(path: Path, payload: dict) -> None:
    with path.open("x", encoding="utf-8") as stream:
        json.dump(payload, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def prepare(output: Path) -> dict:
    if output.exists():
        raise FileExistsError(f"refusing to overwrite Gaia plan directory: {output}")
    audit_summary_path = AUDIT_DIR / "summary.json"
    candidates_path = AUDIT_DIR / "control_candidate_positions.tsv"
    target_run_status_path = TARGET_RUN / "status.json"
    target_response_path = TARGET_RUN / "gaia_matches.csv"
    for path in (audit_summary_path, candidates_path, target_run_status_path,
                 target_response_path, TARGET_REQUESTS):
        if not path.is_file():
            raise FileNotFoundError(path)

    audit = json.loads(audit_summary_path.read_text(encoding="utf-8"))
    target_status = json.loads(target_run_status_path.read_text(encoding="utf-8"))
    if (audit.get("status") != "complete_local_diagnostic_only"
            or audit.get("n_matched_target_control_pairs") != 49
            or audit.get("gaia_query_submitted") is not False
            or audit.get("science_measurement") is not False):
        raise ValueError("matched-control candidate audit is incomplete or exceeds scope")
    if (sha256(candidates_path) != audit["outputs_sha256"][candidates_path.name]
            or sha256(target_response_path) != target_status.get("response_sha256")
            or target_status.get("status") != "complete"):
        raise ValueError("candidate audit or approved target Gaia receipt failed hash validation")

    sources = read_tsv(candidates_path)
    if not sources:
        raise ValueError("no control morphology candidates are available for a Gaia plan")
    requests = []
    candidate_tables = {}
    seen = set()
    for row in sources:
        key = (row["control_id"], int(row["number"]))
        if key in seen:
            raise ValueError(f"duplicate control candidate identity: {key}")
        seen.add(key)
        ra, dec = float(row["ra_deg"]), float(row["dec_deg"])
        if not (math.isfinite(ra) and 0 <= ra < 360 and math.isfinite(dec) and -90 <= dec <= 90):
            raise ValueError(f"invalid control candidate coordinates: {key}")
        if row["classification"] != "provisional_compact_candidate":
            raise ValueError(f"unexpected candidate class at {key}")
        requests.append({
            "local_index": len(requests),
            "cluster_id": row["control_id"],
            "target_cluster_id": row["target_cluster_id"],
            "number": int(row["number"]),
            "ra_deg": ra,
            "dec_deg": dec,
            "classification": row["classification"],
            "primary_filter": row["filter"],
            "fwhm_moment_pix": float(row["fwhm_moment_pix"]),
            "snr_nominal": float(row["snr_nominal"]),
            "wing_valid_fraction": float(row["wing_valid_fraction"]),
        })

    if len(requests) > 49 * 25:
        raise ValueError("control candidate count exceeds the pilot runner's 25-per-frame cap")
    target_requests = json.loads(TARGET_REQUESTS.read_text(encoding="utf-8"))["candidates"]
    target_request_n = len(target_requests)
    if target_status.get("n_candidates") != target_request_n:
        raise ValueError("target Gaia response and approved target request counts differ")
    expected_match_rows = round(len(requests) * target_status["n_match_rows"] / target_request_n)
    expected_response_bytes = round(
        expected_match_rows * target_response_path.stat().st_size / max(1, target_status["n_match_rows"])
    )

    output.mkdir(parents=True, exist_ok=False)
    requests_doc = {"candidates": requests}
    requests_path = output / "requests.json"
    write_json(requests_path, requests_doc)
    positions_path = output / "candidate_positions.tsv"
    position_rows = [{
        "local_index": row["local_index"], "target_cluster_id": row["target_cluster_id"],
        "control_id": row["cluster_id"], "filter": row["primary_filter"],
        "number": row["number"], "ra_deg": row["ra_deg"], "dec_deg": row["dec_deg"],
        "classification": row["classification"],
    } for row in requests]
    write_tsv(positions_path, position_rows, list(position_rows[0]))
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

    for row in requests:
        control_id = row["cluster_id"]
        if control_id in candidate_tables:
            continue
        table_path = CONTROL_PILOT / "preparation" / f"cluster_{control_id}" / "psf_candidates.tsv"
        if not table_path.is_file():
            raise FileNotFoundError(table_path)
        candidate_tables[control_id] = sha256(table_path)

    query_batch_manifest = json.dumps(query_batches, sort_keys=True, separators=(",", ":")).encode()
    plan = {
        "schema": "astrafek-s23b-control-psf-gaia-plan-v3",
        "status": "prepared_not_submitted",
        "external_query": "Gaia DR3 TAP positional cross-match only",
        "endpoint": match_psf_gaia.ENDPOINT,
        "request_parameters": {
            "REQUEST": "doQuery",
            "LANG": "ADQL",
            "FORMAT": "csv",
            "MAXREC": str(MAX_MATCH_ROWS),
        },
        "match_radius_arcsec": 2.0,
        "n_control_frames": len({row["cluster_id"] for row in requests}),
        "candidate_positions_by_filter": {
            band: sum(row["primary_filter"] == band for row in requests)
            for band in ("i", "i2")
        },
        "control_frames_by_filter": {
            band: len({row["cluster_id"] for row in requests if row["primary_filter"] == band})
            for band in ("i", "i2")
        },
        "n_matched_target_control_pairs": 49,
        "n_new_candidate_positions": len(requests),
        "n_candidate_positions": len(requests),
        "batch_size_recommended": BATCH_SIZE,
        "expected_batches": math.ceil(len(requests) / BATCH_SIZE),
        "expected_match_rows_rough_estimate": expected_match_rows,
        "estimated_response_bytes_rough": expected_response_bytes,
        "maxrec_per_batch": MAX_MATCH_ROWS,
        "max_match_rows_before_truncation_guard": MAX_MATCH_ROWS,
        "columns_returned": ["candidate_id", *GAIA_FIELDS, "separation_arcsec"],
        "gaia_columns_queried": GAIA_FIELDS,
        "query_method": (
            "Synchronous Gaia TAP cone-union ADQL; each batch embeds only its approved ICRS "
            "candidate coordinates. Candidate IDs and exact separations are associated locally."
        ),
        "images_requested": False,
        "hsc_shape_rows_requested": False,
        "photoz_requested": False,
        "science_measurement": False,
        "psf_approved": False,
        "candidate_list": str(requests_path.resolve()),
        "candidate_list_sha256": sha256(requests_path),
        "candidate_positions_tsv": str(positions_path.resolve()),
        "candidate_positions_tsv_sha256": sha256(positions_path),
        "query_batches": query_batches,
        "query_batch_manifest_sha256": hashlib.sha256(query_batch_manifest).hexdigest(),
        "input_sha256": {
            str(audit_summary_path.resolve()): sha256(audit_summary_path),
            str(candidates_path.resolve()): sha256(candidates_path),
            str(target_run_status_path.resolve()): sha256(target_run_status_path),
            str(target_response_path.resolve()): sha256(target_response_path),
            str(TARGET_REQUESTS.resolve()): sha256(TARGET_REQUESTS),
        },
        "candidate_tables_sha256": candidate_tables,
        "estimate_basis": (
            f"Scaled from the completed, separately approved 492-position S23B target query "
            f"({target_status['n_match_rows']} rows; {target_response_path.stat().st_size} response bytes)."
        ),
        "notes": [
            "These are morphology-only positions in matched random-sky control images.",
            "The existing approval covered only 492 target positions, not these control positions.",
            "This plan is not submitted; explicit approval of this exact plan is required.",
            "The remote ADQL requests only the listed Gaia DR3 columns and the approved 2-arcsec cones.",
            "HSC candidate metadata stays local; the TAP query contains positions, not HSC measurements.",
            "No images, HSC shape rows, photo-z values, or CAS data are requested.",
        ],
    }
    plan_path = output / "plan.json"
    write_json(plan_path, plan)
    summary = {
        "schema": "astrafek-s23b-control-psf-gaia-plan-summary-v3",
        "status": "prepared_not_submitted",
        "plan_path": str(plan_path.resolve()),
        "plan_sha256": sha256(plan_path),
        "candidate_list_sha256": sha256(requests_path),
        "candidate_positions_tsv_sha256": sha256(positions_path),
        "n_candidate_positions": len(requests),
        "n_control_frames": len({row["cluster_id"] for row in requests}),
        "expected_batches": math.ceil(len(requests) / BATCH_SIZE),
        "estimated_match_rows_rough": expected_match_rows,
        "estimated_response_bytes_rough": expected_response_bytes,
        "images_requested": False,
        "hsc_shape_rows_requested": False,
        "photoz_requested": False,
        "science_measurement": False,
        "query_submitted": False,
    }
    write_json(output / "summary.json", summary)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    summary = prepare(args.output.resolve())
    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
