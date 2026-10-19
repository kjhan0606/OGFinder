#!/usr/bin/env python3
"""Independently verify an approved, batched Gaia DR3 candidate cross-match."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

from match_psf_gaia_batched import GAIA_FIELDS, associate, build_query, sha256, write_json


def read_csv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream)
        return reader.fieldnames or [], list(reader)


def verify(plan_dir: Path, work: Path, expected_plan_sha256: str,
           output: Path) -> dict:
    plan_path = plan_dir / "plan.json"
    request_path = plan_dir / "candidates.json"
    if sha256(plan_path) != expected_plan_sha256:
        raise ValueError("plan checksum differs from the user's approved plan")
    plan = json.loads(plan_path.read_text(encoding="utf-8"))
    targets = json.loads(request_path.read_text(encoding="utf-8"))["candidates"]
    status = json.loads((work / "status.json").read_text(encoding="utf-8"))
    if (plan.get("status") != "prepared_not_submitted"
            or status.get("status") != "complete"
            or status.get("plan_sha256") != expected_plan_sha256
            or status.get("request_sha256") != sha256(request_path)
            or status.get("n_candidates") != len(targets)
            or status.get("psf_approved") is not False
            or status.get("science_measurement") is not False):
        raise ValueError("plan, request, or execution status does not match")
    final_columns = ["candidate_id", *GAIA_FIELDS, "separation_arcsec"]
    fields, final_rows = read_csv(work / "gaia_matches.csv")
    if fields != final_columns or len(final_rows) != status.get("n_match_rows"):
        raise ValueError("final CSV schema or row count differs from execution status")
    if sha256(work / "gaia_matches.csv") != status.get("response_sha256"):
        raise ValueError("final CSV checksum differs from execution status")

    expected_batches = (len(targets) + status["batch_size"] - 1) // status["batch_size"]
    if len(status.get("batches", [])) != expected_batches:
        raise ValueError("batch manifest is incomplete")
    batch_matches = []
    for index, record in enumerate(status["batches"]):
        first = index * status["batch_size"]
        batch_targets = targets[first:first + status["batch_size"]]
        batch_dir = work / f"batch_{index:03d}"
        batch_status = json.loads((batch_dir / "status.json").read_text(encoding="utf-8"))
        if (record.get("status") != "complete" or batch_status != record
                or record.get("first_candidate_id") != first
                or record.get("n_candidates") != len(batch_targets)):
            raise ValueError(f"incomplete or misaligned batch {index}")
        query_path = batch_dir / "query.adql"
        raw_path = batch_dir / "gaia_cone_union.csv"
        if (query_path.read_text(encoding="utf-8") != build_query(batch_targets) + "\n"
                or sha256(query_path) != record.get("query_sha256")
                or sha256(raw_path) != record.get("response_sha256")):
            raise ValueError(f"query coordinates or raw-response checksum mismatch in batch {index}")
        raw_fields, sources = read_csv(raw_path)
        if raw_fields != GAIA_FIELDS or len(sources) != record.get("n_gaia_union_rows"):
            raise ValueError(f"raw Gaia response schema or count mismatch in batch {index}")
        recalculated = associate(batch_targets, sources, first)
        pair_fields, saved_pairs = read_csv(batch_dir / "matched_pairs.csv")
        if pair_fields != final_columns or len(saved_pairs) != len(recalculated):
            raise ValueError(f"local match-pair output mismatch in batch {index}")
        saved_keys = {(r["candidate_id"], r["source_id"]): float(r["separation_arcsec"])
                      for r in saved_pairs}
        recalc_keys = {(r["candidate_id"], r["source_id"]): float(r["separation_arcsec"])
                       for r in recalculated}
        if saved_keys.keys() != recalc_keys.keys() or any(
                abs(saved_keys[key] - recalc_keys[key]) > 1e-8 for key in saved_keys):
            raise ValueError(f"candidate associations differ from raw Gaia rows in batch {index}")
        batch_matches.extend(recalculated)

    final_keys = {(r["candidate_id"], r["source_id"]): float(r["separation_arcsec"])
                  for r in final_rows}
    calculated_keys = {(r["candidate_id"], r["source_id"]): float(r["separation_arcsec"])
                       for r in batch_matches}
    if final_keys.keys() != calculated_keys.keys():
        raise ValueError("assembled output pairs differ from per-batch results")
    if len(final_keys) != len(final_rows):
        raise ValueError("duplicate candidate–Gaia source pairs")
    if any(not 0 <= separation <= 2.00001 for separation in final_keys.values()):
        raise ValueError("a Gaia match lies outside the approved 2-arcsec radius")

    evidence_path = work / "evidence_summary.json"
    evidence = json.loads(evidence_path.read_text(encoding="utf-8"))
    if (evidence.get("n_candidates") != len(targets)
            or evidence.get("science_measurement") is not False
            or evidence.get("psf_approved") is not False
            or sum(evidence["status_counts"].values()) != len(targets)):
        raise ValueError("evidence summary is inconsistent with the approved candidate set")
    if output.exists():
        raise FileExistsError(f"refusing to overwrite verification report: {output}")
    separations = list(final_keys.values())
    result = {
        "schema": "astrafek-psf-gaia-batched-audit-v1",
        "status": "verified",
        "approved_plan_sha256": expected_plan_sha256,
        "candidate_list_sha256": sha256(request_path),
        "response_csv_sha256": sha256(work / "gaia_matches.csv"),
        "evidence_summary_sha256": sha256(evidence_path),
        "n_candidates": len(targets),
        "n_batches": len(status["batches"]),
        "n_match_rows": len(final_rows),
        "n_unique_gaia_sources": len({row["source_id"] for row in final_rows}),
        "n_matched_candidates": len({row["candidate_id"] for row in final_rows}),
        "separation_arcsec_min": min(separations) if separations else None,
        "separation_arcsec_max": max(separations) if separations else None,
        "evidence_status_counts": evidence["status_counts"],
        "images_requested": False,
        "hsc_shape_rows_requested": False,
        "photoz_requested": False,
        "science_measurement": False,
        "psf_approved": False,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    write_json(output, result)
    result["verification_report_sha256"] = sha256(output)
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan-dir", type=Path, required=True)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--approved-plan-sha256", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(verify(args.plan_dir, args.work,
                            args.approved_plan_sha256, args.output), indent=2))


if __name__ == "__main__":
    main()
