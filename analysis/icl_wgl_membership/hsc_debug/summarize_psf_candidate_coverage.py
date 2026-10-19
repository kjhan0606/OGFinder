#!/usr/bin/env python3
"""Summarize the audited 800-frame PSF-core candidate coverage.

These are morphology-selected compact candidates, not confirmed stars or
science-approved PSFs. The script verifies each candidate table against the
independent artifact audit and reports only core-candidate availability.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import statistics
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def summarize(prepared: Path, *, expected_frames: int = 800) -> dict:
    prepared = Path(prepared).resolve()
    summary = json.loads((prepared / "summary.json").read_text(encoding="utf-8"))
    audit_path = prepared / "artifact_audit.json"
    audit = json.loads(audit_path.read_text(encoding="utf-8"))
    if (summary.get("status") != "complete" or summary.get("science_measurement") is not False
            or audit.get("passed") is not True or audit.get("science_measurement") is not False
            or audit.get("mask_ready_for_science") is not False
            or audit.get("n_frames") != expected_frames):
        raise ValueError(f"expected a complete, independently audited, non-science {expected_frames}-frame preparation")
    table_hashes = audit.get("table_sha256", {})
    records = audit.get("records", [])
    if len(records) != expected_frames or len({str(row.get("id")) for row in records}) != expected_frames:
        raise ValueError("audited frame records are missing or duplicated")

    by_filter = {"i": [], "i2": []}
    classifications = {}
    total_candidates = 0
    verified_tables = 0
    for record in records:
        filt = record.get("filter")
        if filt not in by_filter:
            raise ValueError(f"unexpected coadd filter {filt}")
        path = prepared / f"cluster_{record['id']}" / "psf_candidates.tsv"
        expected_sha = table_hashes.get(str(path))
        if expected_sha is None or sha256(path) != expected_sha:
            raise ValueError(f"PSF candidate table is missing or changed: {path}")
        with path.open(encoding="utf-8", newline="") as stream:
            rows = list(csv.DictReader(stream, delimiter="\t"))
        if len(rows) != int(record["n_psf_candidates"]):
            raise ValueError(f"candidate count differs from independent audit: {path}")
        fwhm_values = []
        for row in rows:
            classification = row.get("CLASSIFICATION")
            classifications[classification] = classifications.get(classification, 0) + 1
            if classification != "provisional_compact_candidate":
                raise ValueError("candidate table contains an unexpected classification")
            fwhm = float(row["FWHM_MOMENT_PIX"])
            snr = float(row["SNR_NOMINAL"])
            wing = float(row["WING_VALID_FRACTION"])
            if not math.isfinite(fwhm) or fwhm <= 0 or not math.isfinite(snr) or snr < 50:
                raise ValueError("candidate core measurements are invalid")
            if not math.isfinite(wing) or not 0 <= wing <= 1:
                raise ValueError("candidate stamp coverage is invalid")
            fwhm_values.append(fwhm)
        by_filter[filt].append({
            "id": str(record["id"]), "n_candidates": len(rows),
            "median_fwhm_moment_pixels": statistics.median(fwhm_values) if fwhm_values else None,
        })
        total_candidates += len(rows)
        verified_tables += 1

    result_by_filter = {}
    for filt, frames in by_filter.items():
        available = [row["median_fwhm_moment_pixels"] for row in frames
                     if row["median_fwhm_moment_pixels"] is not None]
        result_by_filter[filt] = {
            "n_frames": len(frames),
            "n_frames_with_at_least_4_candidates": sum(row["n_candidates"] >= 4 for row in frames),
            "n_frames_with_no_candidates": sum(row["n_candidates"] == 0 for row in frames),
            "n_candidate_rows": sum(row["n_candidates"] for row in frames),
            "median_of_frame_median_fwhm_pixels": statistics.median(available) if available else None,
            "range_of_frame_median_fwhm_pixels": [min(available), max(available)] if available else None,
        }
    return {
        "schema": "astrafek-psf-candidate-coverage-summary-v1",
        "status": "complete_descriptive",
        "science_measurement": False,
        "psf_approved": False,
        "psf_wings_validated": False,
        "prepared_directory": str(prepared),
        "preparation_summary_sha256": sha256(prepared / "summary.json"),
        "independent_artifact_audit_sha256": sha256(audit_path),
        "n_frames_verified": verified_tables,
        "n_candidate_rows": total_candidates,
        "classification_counts": classifications,
        "by_filter": result_by_filter,
        "core_stamp_side_pixels": 101,
        "core_candidate_interpretation": "morphology-only compact candidates; not externally confirmed stars",
        "limitations": [
            "Candidates come from each target/coadd image and can include compact galaxies.",
            "The 101x101-pixel stamps extend only 50 pixels from the center; they do not validate extended PSF wings relevant to ICL outskirts.",
            "No per-position epoch-matched PSF model, recentering, background subtraction, or flux normalization is certified.",
            "Candidate coverage does not activate a PSF correction or science measurement.",
        ],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prepared", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    output = args.output.resolve()
    if output.exists():
        parser.error(f"refusing to overwrite {output}")
    result = summarize(args.prepared)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False) + "\n",
                      encoding="utf-8")
    print(json.dumps({"status": result["status"], "science_measurement": False,
                      "n_frames_verified": result["n_frames_verified"],
                      "n_candidate_rows": result["n_candidate_rows"],
                      "by_filter": result["by_filter"], "output": str(output)},
                     ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
