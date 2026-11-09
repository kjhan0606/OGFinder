#!/usr/bin/env python3
"""Compare morphology-only PSF candidates in matched S23B target/control fields.

This offline audit prepares a hash-pinned control-candidate table for a
possible, separately approved Gaia query. It performs no external query or
download and makes no science measurement.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from statistics import median


DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
TARGET_PILOT = DATA / "s23b_detection_psf_pilot_2026-11-07_01"
CONTROL_PILOT = DATA / "s23b_psf_control_candidates_2026-11-07_01"
DEFAULT_OUTPUT = DATA / "s23b_psf_control_candidate_audit_2026-11-07_01"


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


def audit(output: Path) -> dict:
    if output.exists():
        raise FileExistsError(f"refusing to overwrite candidate audit: {output}")
    target_manifest_path = TARGET_PILOT / "manifest.json"
    target_plan_path = TARGET_PILOT / "plan.json"
    target_summary_path = TARGET_PILOT / "preparation/summary.json"
    control_pilot_manifest_path = CONTROL_PILOT / "pilot_manifest.json"
    control_pilot_summary_path = CONTROL_PILOT / "summary.json"
    control_preparation_summary_path = CONTROL_PILOT / "preparation/summary.json"
    inputs = (
        target_manifest_path, target_plan_path, target_summary_path,
        control_pilot_manifest_path, control_pilot_summary_path,
        control_preparation_summary_path,
    )
    for path in inputs:
        if not path.is_file():
            raise FileNotFoundError(path)

    target_manifest = json.loads(target_manifest_path.read_text(encoding="utf-8"))
    target_plan = json.loads(target_plan_path.read_text(encoding="utf-8"))
    target_summary = json.loads(target_summary_path.read_text(encoding="utf-8"))
    control_manifest = json.loads(control_pilot_manifest_path.read_text(encoding="utf-8"))
    control_plan = json.loads(control_pilot_summary_path.read_text(encoding="utf-8"))
    control_summary = json.loads(control_preparation_summary_path.read_text(encoding="utf-8"))

    if target_summary.get("n_prepared") != 49 or target_summary.get("science_measurement") is not False:
        raise ValueError("target pilot must have 49 prepared diagnostic-only frames")
    if target_summary.get("manifest_sha256") != sha256(target_manifest_path):
        raise ValueError("target preparation summary does not pin the target manifest")
    if control_plan.get("pilot_manifest_sha256") != sha256(control_pilot_manifest_path):
        raise ValueError("control pilot summary does not pin its prepared manifest")
    if (control_summary.get("status") != "complete" or control_summary.get("n_requested") != 49
            or control_summary.get("n_prepared") != 49 or control_summary.get("errors")):
        raise ValueError("matched-control preparation did not complete all 49 frames")
    if control_summary.get("manifest_sha256") != sha256(control_pilot_manifest_path):
        raise ValueError("control candidate preparation did not use the pinned manifest")
    if control_summary.get("science_measurement") is not False:
        raise ValueError("control preparation unexpectedly promotes science")

    target_selection = {row["id"]: row for row in target_plan["selection"]}
    target_frames = {row["id"]: row for row in target_manifest["clusters"]}
    control_by_target = {row["target_cluster_id"]: row for row in control_manifest["clusters"]}
    target_prepared = {row["id"]: row for row in target_summary["records"]}
    control_prepared = {row["id"]: row for row in control_summary["records"]}
    if len(target_prepared) != 49 or len(control_prepared) != 49:
        raise ValueError("target/control preparation record count is not 49")

    pair_rows = []
    candidate_rows = []
    table_hashes = {}
    target_candidate_total = 0
    control_candidate_total = 0
    for target_id in target_prepared:
        target = target_frames[target_id]
        selection = target_selection[target_id]
        control = control_by_target[target_id]
        control_id = control["id"]
        if control_id not in control_prepared:
            raise ValueError(f"prepared control is missing: {control_id}")
        target_dir = TARGET_PILOT / "preparation" / f"cluster_{target_id}"
        control_dir = CONTROL_PILOT / "preparation" / f"cluster_{control_id}"
        target_candidates_path = target_dir / "psf_candidates.tsv"
        target_detections_path = target_dir / "detections.tsv"
        control_candidates_path = control_dir / "psf_candidates.tsv"
        control_detections_path = control_dir / "detections.tsv"
        for path in (target_candidates_path, target_detections_path,
                     control_candidates_path, control_detections_path):
            if not path.is_file():
                raise FileNotFoundError(path)
            table_hashes[str(path.resolve())] = sha256(path)

        target_candidates = read_tsv(target_candidates_path)
        target_detections = {row["NUMBER"]: row for row in read_tsv(target_detections_path)}
        control_candidates = read_tsv(control_candidates_path)
        control_detections = {row["NUMBER"]: row for row in read_tsv(control_detections_path)}
        expected_filter = selection["primary_filter"]
        if control.get("primary_filter") != expected_filter:
            raise ValueError(f"target/control filter mismatch for {target_id}")
        if (len(target_candidates) != target_prepared[target_id]["n_psf_candidates"]
                or len(control_candidates) != control_prepared[control_id]["n_psf_candidates"]):
            raise ValueError(f"candidate table count mismatch for pair {target_id}/{control_id}")

        def numeric(rows: list[dict[str, str]], key: str) -> list[float]:
            values = [float(row[key]) for row in rows if row.get(key) not in (None, "")]
            if not all(math.isfinite(value) for value in values):
                raise ValueError(f"non-finite {key} in {target_id}/{control_id}")
            return values

        target_fwhm = numeric(target_candidates, "FWHM_MOMENT_PIX")
        control_fwhm = numeric(control_candidates, "FWHM_MOMENT_PIX")
        target_snr = numeric(target_candidates, "SNR_NOMINAL")
        control_snr = numeric(control_candidates, "SNR_NOMINAL")
        pair_rows.append({
            "target_cluster_id": target_id,
            "control_id": control_id,
            "filter": expected_filter,
            "target_n_detections": target_prepared[target_id]["n_detections"],
            "control_n_detections": control_prepared[control_id]["n_detections"],
            "target_n_candidates": len(target_candidates),
            "control_n_candidates": len(control_candidates),
            "candidate_count_difference_control_minus_target": len(control_candidates) - len(target_candidates),
            "target_fwhm_locus_pix": target_prepared[target_id].get("psf_fwhm_locus_pix"),
            "control_fwhm_locus_pix": control_prepared[control_id].get("psf_fwhm_locus_pix"),
            "target_candidate_fwhm_median_pix": median(target_fwhm) if target_fwhm else None,
            "control_candidate_fwhm_median_pix": median(control_fwhm) if control_fwhm else None,
            "target_candidate_snr_median": median(target_snr) if target_snr else None,
            "control_candidate_snr_median": median(control_snr) if control_snr else None,
            "control_match_score": control["match_score"],
            "target_seeing_arcsec": control["target_seeing_arcsec"],
            "control_seeing_arcsec": control["candidate_seeing_arcsec"],
            "center_note": "target CAMIRA centroid / control random sky position; not a BCG label",
        })
        target_candidate_total += len(target_candidates)
        control_candidate_total += len(control_candidates)

        for candidate in control_candidates:
            number = candidate["NUMBER"]
            detection = control_detections.get(number)
            if detection is None:
                raise ValueError(f"control PSF candidate lacks detection row: {control_id}/{number}")
            if candidate["CLASSIFICATION"] != "provisional_compact_candidate":
                raise ValueError(f"unexpected control candidate class: {control_id}/{number}")
            candidate_rows.append({
                "local_index": len(candidate_rows),
                "target_cluster_id": target_id,
                "control_id": control_id,
                "filter": expected_filter,
                "number": number,
                "ra_deg": detection["RA"],
                "dec_deg": detection["DEC"],
                "fwhm_moment_pix": candidate["FWHM_MOMENT_PIX"],
                "snr_nominal": candidate["SNR_NOMINAL"],
                "wing_valid_fraction": candidate["WING_VALID_FRACTION"],
                "classification": candidate["CLASSIFICATION"],
            })

    output.mkdir(parents=True, exist_ok=False)
    pair_path = output / "paired_candidate_counts.tsv"
    candidate_path = output / "control_candidate_positions.tsv"
    write_tsv(pair_path, pair_rows, list(pair_rows[0]))
    write_tsv(candidate_path, candidate_rows, list(candidate_rows[0]))
    by_filter = {}
    for band in ("i", "i2"):
        selected = [row for row in pair_rows if row["filter"] == band]
        by_filter[band] = {
            "n_pairs": len(selected),
            "target_candidate_total": sum(row["target_n_candidates"] for row in selected),
            "control_candidate_total": sum(row["control_n_candidates"] for row in selected),
            "target_frames_with_at_least_four": sum(row["target_n_candidates"] >= 4 for row in selected),
            "control_frames_with_at_least_four": sum(row["control_n_candidates"] >= 4 for row in selected),
            "median_target_candidate_count": median(row["target_n_candidates"] for row in selected),
            "median_control_candidate_count": median(row["control_n_candidates"] for row in selected),
        }
    summary = {
        "schema": "astrafek-s23b-control-psf-candidate-audit-v1",
        "status": "complete_local_diagnostic_only",
        "n_matched_target_control_pairs": len(pair_rows),
        "n_target_morphology_candidates": target_candidate_total,
        "n_control_morphology_candidates": control_candidate_total,
        "n_control_candidate_positions_for_possible_gaia_plan": len(candidate_rows),
        "by_filter": by_filter,
        "input_sha256": {str(path.resolve()): sha256(path) for path in inputs},
        "candidate_table_sha256": table_hashes,
        "outputs_sha256": {
            pair_path.name: sha256(pair_path),
            candidate_path.name: sha256(candidate_path),
        },
        "science_measurement": False,
        "gaia_query_submitted": False,
        "notes": [
            "The target/control pairs are the 49 exact matches from the existing stratified target pilot.",
            "Candidate positions are morphology-only; Gaia confirmation has not been requested for controls.",
            "The target approval covered only 492 target positions and does not authorize these control positions.",
            "Central-position diagnostics are not BCG identifications; controls are random sky positions.",
            "Counts capped at 25 candidates per field by the shared engineering runner.",
        ],
    }
    write_json(output / "summary.json", summary)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    summary = audit(args.output.resolve())
    print(json.dumps({key: summary[key] for key in (
        "status", "n_matched_target_control_pairs", "n_target_morphology_candidates",
        "n_control_morphology_candidates", "by_filter", "science_measurement",
        "gaia_query_submitted",
    )}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
