#!/usr/bin/env python3
"""Audit the approved 4,299-position S23B Gaia expansion locally.

This read-only audit reproduces each batch's candidate/source associations
from its saved Gaia cone-union response, verifies the approved query and all
response hashes, and writes diagnostic evidence tables. It performs no
network access and does not approve PSFs or promote science measurements.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import sys
from collections import Counter, defaultdict
from pathlib import Path

from match_psf_gaia_batched import GAIA_FIELDS, associate, build_query


DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
EXPANSION = DATA / "s23b_psf_expansion_2026-11-07_01"
PLAN_DIR = EXPANSION / "gaia_plan"
RUN_DIR = EXPANSION / "gaia_match_2026-11-07_01"
TARGET_MANIFEST = EXPANSION / "plan/target_manifest.json"
CONTROL_MANIFEST = EXPANSION / "plan/control_manifest.json"
TARGET_PREPARATION = EXPANSION / "target_preparation/summary.json"
CONTROL_PREPARATION = EXPANSION / "control_preparation/summary.json"
APPROVED_PLAN_SHA256 = "c00eda5f896f8160022a7e4198b3f8ce69803f0e30adc8f560ce04c8beba0e7d"
EXPECTED_COLUMNS = ["candidate_id", *GAIA_FIELDS, "separation_arcsec"]
DEFAULT_OUTPUT = EXPANSION / "gaia_audit_2026-11-07_01"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


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


def number(value: str | None) -> float | None:
    if value in (None, ""):
        return None
    result = float(value)
    return result if math.isfinite(result) else None


def flag(value: str | None) -> bool | None:
    if value in (None, ""):
        return None
    lowered = value.strip().lower()
    if lowered in ("true", "1", "t"):
        return True
    if lowered in ("false", "0", "f"):
        return False
    raise ValueError(f"unexpected Gaia boolean value: {value!r}")


def audit(output: Path) -> dict:
    if output.exists():
        raise FileExistsError(f"refusing to overwrite Gaia audit directory: {output}")

    plan_path = PLAN_DIR / "plan.json"
    requests_path = PLAN_DIR / "requests.json"
    position_path = PLAN_DIR / "candidate_positions.tsv"
    status_path = RUN_DIR / "status.json"
    response_path = RUN_DIR / "gaia_matches.csv"
    candidate_snapshot = RUN_DIR / "candidates.json"
    for path in (plan_path, requests_path, position_path, status_path,
                 response_path, candidate_snapshot, TARGET_MANIFEST,
                 CONTROL_MANIFEST, TARGET_PREPARATION, CONTROL_PREPARATION):
        if not path.is_file():
            raise FileNotFoundError(path)

    plan = read_json(plan_path)
    request_doc = read_json(requests_path)
    requests = request_doc.get("candidates")
    status = read_json(status_path)
    if sha256(plan_path) != APPROVED_PLAN_SHA256:
        raise ValueError("Gaia plan SHA-256 differs from the user's approved plan")
    if (plan.get("schema") != "astrafek-s23b-matched-psf-expansion-gaia-plan-v1"
            or plan.get("status") != "prepared_not_submitted"
            or not isinstance(requests, list)
            or len(requests) != 4299
            or plan.get("n_candidate_positions") != len(requests)
            or plan.get("candidate_list_sha256") != sha256(requests_path)
            or plan.get("candidate_positions_tsv_sha256") != sha256(position_path)
            or plan.get("match_radius_arcsec") != 2.0
            or plan.get("columns_returned") != EXPECTED_COLUMNS
            or plan.get("expected_batches") != 18
            or plan.get("batch_size_recommended") != 250):
        raise ValueError("approved plan/request linkage or query definition failed")
    for key in ("images_requested", "hsc_shape_rows_requested", "photoz_requested",
                "cas_query_requested", "science_measurement", "psf_approved"):
        if plan.get(key) is not False:
            raise ValueError(f"approved Gaia plan has unexpected scope/state for {key}")
    if json.loads(candidate_snapshot.read_text(encoding="utf-8")) != requests:
        raise ValueError("submitted candidate snapshot differs from the approved request")
    if (status.get("status") != "complete"
            or status.get("plan_sha256") != sha256(plan_path)
            or status.get("request_sha256") != sha256(requests_path)
            or status.get("response_sha256") != sha256(response_path)
            or status.get("n_candidates") != len(requests)
            or status.get("match_radius_arcsec") != 2.0
            or status.get("columns") != EXPECTED_COLUMNS
            or status.get("completed_batches") != 18
            or len(status.get("batches", [])) != 18):
        raise ValueError("completed run receipt does not match the approved plan")
    for key in ("images_requested", "hsc_shape_rows_requested", "photoz_requested",
                "science_measurement", "psf_approved"):
        if status.get(key) is not False:
            raise ValueError(f"Gaia run has unexpected scope/state for {key}")

    target_manifest = read_json(TARGET_MANIFEST)
    control_manifest = read_json(CONTROL_MANIFEST)
    target_preparation = read_json(TARGET_PREPARATION)
    control_preparation = read_json(CONTROL_PREPARATION)
    for prep, manifest_path, population, expected_n in (
        (target_preparation, TARGET_MANIFEST, "target", 200),
        (control_preparation, CONTROL_MANIFEST, "control", 200),
    ):
        if (prep.get("n_requested") != expected_n
                or prep.get("manifest_sha256") != sha256(manifest_path)
                or prep.get("science_measurement") is not False):
            raise ValueError(f"{population} candidate preparation no longer matches its manifest")
    manifest_by_population = {
        "target": {row["id"]: row for row in target_manifest["clusters"]},
        "control": {row["id"]: row for row in control_manifest["clusters"]},
    }
    prep_by_population = {
        "target": {row["id"]: row for row in target_preparation["records"]},
        "control": {row["id"]: row for row in control_preparation["records"]},
    }
    prep_table_hashes = plan.get("candidate_preparation_tables_sha256", {})
    candidate_tables: dict[tuple[str, str], dict[int, dict[str, str]]] = {}
    detection_tables: dict[tuple[str, str], dict[int, dict[str, str]]] = {}

    local_indices = [int(row["local_index"]) for row in requests]
    if sorted(local_indices) != list(range(len(requests))):
        raise ValueError("request local_index values are incomplete or duplicated")
    request_by_index = {int(row["local_index"]): row for row in requests}
    positions = read_tsv(position_path)
    if len(positions) != len(requests):
        raise ValueError("candidate positions table count differs from Gaia request count")
    position_by_index = {int(row["local_index"]): row for row in positions}
    if len(position_by_index) != len(positions):
        raise ValueError("candidate positions table contains duplicate local indices")

    for index, request in request_by_index.items():
        population = request.get("population")
        cluster_id = request.get("cluster_id")
        if population not in ("target", "control") or cluster_id not in manifest_by_population[population]:
            raise ValueError(f"candidate {index} refers to an unknown expansion frame")
        frame = manifest_by_population[population][cluster_id]
        prep_record = prep_by_population[population].get(cluster_id)
        if prep_record is None:
            raise ValueError(f"candidate {index} frame is absent from its preparation summary")
        pair_id = frame.get("paired_control_id") if population == "target" else frame.get("paired_target_cluster_id")
        if (request.get("image_sha256") != frame.get("input_image_sha256")
                or request.get("primary_filter") != frame.get("coadd_header_filter")
                or (population == "target" and request.get("control_id") != pair_id)
                or (population == "control" and request.get("target_cluster_id") != pair_id)):
            raise ValueError(f"candidate {index} no longer maps to its paired audited frame")
        table_path = EXPANSION / f"{population}_preparation/cluster_{cluster_id}/psf_candidates.tsv"
        detections_path = EXPANSION / f"{population}_preparation/cluster_{cluster_id}/detections.tsv"
        expected_table_hash = prep_table_hashes.get(str(table_path))
        if (expected_table_hash is None or sha256(table_path) != expected_table_hash
                or request.get("preparation_candidate_table_sha256") != expected_table_hash):
            raise ValueError(f"candidate table hash differs for {population}/{cluster_id}")
        expected_detections_hash = prep_table_hashes.get(str(detections_path))
        if (expected_detections_hash is None
                or sha256(detections_path) != expected_detections_hash):
            raise ValueError(f"detection table hash differs for {population}/{cluster_id}")
        table_key = (population, cluster_id)
        if table_key not in candidate_tables:
            rows = read_tsv(table_path)
            candidate_tables[table_key] = {int(row["NUMBER"]): row for row in rows}
            detection_rows = read_tsv(detections_path)
            detection_tables[table_key] = {int(row["NUMBER"]): row for row in detection_rows}
        det = candidate_tables[table_key].get(int(request["number"]))
        source = detection_tables[table_key].get(int(request["number"]))
        if (det is None or source is None
                or det.get("CLASSIFICATION") != "provisional_compact_candidate"
                or abs(float(source["RA"]) - float(request["ra_deg"])) > 1e-9
                or abs(float(source["DEC"]) - float(request["dec_deg"])) > 1e-9):
            raise ValueError(f"candidate {index} no longer maps to its morphology detection")
        position = position_by_index[index]
        if (position.get("population") != population
                or position.get("frame_id") != cluster_id
                or int(position.get("number", -1)) != int(request["number"])
                or abs(float(position["ra_deg"]) - float(request["ra_deg"])) > 1e-9
                or abs(float(position["dec_deg"]) - float(request["dec_deg"])) > 1e-9):
            raise ValueError(f"candidate {index} differs between request and pinned position table")

    all_matches: list[dict[str, str]] = []
    batch_audit_rows = []
    planned_batches = plan.get("query_batches")
    if not isinstance(planned_batches, list) or len(planned_batches) != 18:
        raise ValueError("approved plan does not contain its exact 18-batch SQL manifest")
    for batch_number in range(18):
        offset = batch_number * 250
        batch_targets = requests[offset:offset + 250]
        planned = planned_batches[batch_number]
        receipt = status["batches"][batch_number]
        batch_dir = RUN_DIR / f"batch_{batch_number:03d}"
        query_path = batch_dir / "query.adql"
        union_path = batch_dir / "gaia_cone_union.csv"
        pairs_path = batch_dir / "matched_pairs.csv"
        batch_status_path = batch_dir / "status.json"
        for path in (query_path, union_path, pairs_path, batch_status_path):
            if not path.is_file():
                raise FileNotFoundError(path)
        batch_status = read_json(batch_status_path)
        expected_query = build_query(batch_targets) + "\n"
        query_hash = sha256(query_path)
        if (planned.get("batch") != batch_number
                or planned.get("first_candidate_id") != offset
                or planned.get("n_candidates") != len(batch_targets)
                or query_path.read_text(encoding="utf-8") != expected_query
                or query_hash != planned.get("query_sha256")
                or query_hash != receipt.get("query_sha256")
                or batch_status.get("status") != "complete"
                or batch_status.get("query_sha256") != query_hash
                or batch_status.get("response_sha256") != sha256(union_path)
                or receipt.get("response_sha256") != sha256(union_path)
                or receipt.get("status") != "complete"):
            raise ValueError(f"batch {batch_number} query or response hash failed")
        with union_path.open(encoding="utf-8-sig", newline="") as stream:
            reader = csv.DictReader(stream)
            if not set(GAIA_FIELDS).issubset(reader.fieldnames or []):
                raise ValueError(f"batch {batch_number} has unexpected Gaia source columns")
            sources = list(reader)
        source_ids = [row["source_id"] for row in sources]
        if (len(sources) != batch_status.get("n_gaia_union_rows")
                or len(source_ids) != len(set(source_ids))):
            raise ValueError(f"batch {batch_number} Gaia union table is incomplete or duplicated")
        expected_matches = associate(batch_targets, sources, offset)
        with pairs_path.open(encoding="utf-8", newline="") as stream:
            stored_matches = list(csv.DictReader(stream))
        if stored_matches != expected_matches:
            raise ValueError(f"batch {batch_number} stored associations do not reproduce from Gaia rows")
        if len(expected_matches) != receipt.get("n_match_rows"):
            raise ValueError(f"batch {batch_number} match row count differs from its receipt")
        all_matches.extend(expected_matches)
        batch_audit_rows.append({
            "batch": batch_number,
            "first_candidate_id": offset,
            "n_candidates": len(batch_targets),
            "n_gaia_union_rows": len(sources),
            "n_match_rows": len(expected_matches),
            "query_sha256": query_hash,
            "raw_response_sha256": sha256(union_path),
            "match_table_sha256": sha256(pairs_path),
            "status": "verified",
        })

    with response_path.open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames != EXPECTED_COLUMNS:
            raise ValueError(f"aggregate response columns differ: {reader.fieldnames}")
        response_rows = list(reader)
    if response_rows != all_matches or len(response_rows) != status.get("n_match_rows"):
        raise ValueError("aggregate Gaia response differs from the verified batch associations")

    matches_by_index: dict[int, list[dict[str, str]]] = {index: [] for index in request_by_index}
    seen_pairs: set[tuple[int, str]] = set()
    max_separation = 0.0
    for row in response_rows:
        index = int(row["candidate_id"])
        if index not in request_by_index:
            raise ValueError(f"aggregate response contains unknown candidate index {index}")
        pair = (index, row["source_id"])
        if pair in seen_pairs:
            raise ValueError(f"duplicate candidate/source association: {pair}")
        seen_pairs.add(pair)
        separation = number(row["separation_arcsec"])
        if separation is None or separation > 2.00001:
            raise ValueError("a Gaia match lies outside the approved 2-arcsec cone")
        max_separation = max(max_separation, separation)
        matches_by_index[index].append(row)

    evidence_rows = []
    levels: Counter[str] = Counter()
    by_population_filter: dict[tuple[str, str], Counter[str]] = defaultdict(Counter)
    strict_frames: dict[tuple[str, str], set[str]] = defaultdict(set)
    for index in range(len(requests)):
        request = request_by_index[index]
        matches = sorted(matches_by_index[index],
                         key=lambda row: (float(row["separation_arcsec"]), row["source_id"]))
        population = request["population"]
        band = request["primary_filter"]
        key = (population, band)
        detail: dict[str, object] = {}
        if not matches:
            level = "no_gaia_match"
        else:
            best = matches[0]
            parallax_snr = number(best["parallax_over_error"])
            pmra, pmdec = number(best["pmra"]), number(best["pmdec"])
            pmra_error, pmdec_error = number(best["pmra_error"]), number(best["pmdec_error"])
            ruwe = number(best["ruwe"])
            star_probability = number(best["classprob_dsc_combmod_star"])
            duplicated = flag(best["duplicated_source"])
            qso = flag(best["in_qso_candidates"])
            galaxy = flag(best["in_galaxy_candidates"])
            pm_terms = []
            if pmra is not None and pmra_error not in (None, 0.0):
                pm_terms.append((pmra / pmra_error) ** 2)
            if pmdec is not None and pmdec_error not in (None, 0.0):
                pm_terms.append((pmdec / pmdec_error) ** 2)
            pm_significance = math.sqrt(sum(pm_terms)) if pm_terms else None
            pm_total = math.hypot(pmra, pmdec) if pmra is not None and pmdec is not None else None
            clean_fit = (ruwe is not None and ruwe < 1.4 and duplicated is False
                         and qso is not True and galaxy is not True)
            significant_motion = ((parallax_snr is not None and parallax_snr >= 3.0)
                                  or (pm_significance is not None and pm_significance >= 3.0))
            astrometric_evidence = len(matches) == 1 and clean_fit and significant_motion
            strict = (astrometric_evidence and star_probability is not None
                      and star_probability >= 0.5)
            if len(matches) > 1:
                level = "ambiguous_cone"
            elif strict:
                level = "strict_candidate"
            elif astrometric_evidence:
                level = "clean_astrometric_evidence"
            else:
                level = "matched_no_clean_astrometry"
            pm_shift = {years: (pm_total * years / 1000.0 if pm_total is not None else None)
                        for years in (2, 5, 10)}
            detail = {
                "best_source_id": best["source_id"],
                "best_separation_arcsec": float(best["separation_arcsec"]),
                "gaia_ref_epoch": number(best["ref_epoch"]),
                "parallax_over_error": parallax_snr,
                "pmra_mas_per_year": pmra,
                "pmdec_mas_per_year": pmdec,
                "pm_total_mas_per_year": pm_total,
                "pm_significance": pm_significance,
                "pm_shift_arcsec_at_2yr": pm_shift[2],
                "pm_shift_arcsec_at_5yr": pm_shift[5],
                "pm_shift_arcsec_at_10yr": pm_shift[10],
                "ruwe": ruwe,
                "duplicated_source": duplicated,
                "in_qso_candidates": qso,
                "in_galaxy_candidates": galaxy,
                "star_probability": star_probability,
                "phot_g_mean_mag": number(best["phot_g_mean_mag"]),
                "phot_variable_flag": best["phot_variable_flag"],
                "astrometric_star_evidence": bool(astrometric_evidence),
                "strict_candidate": bool(strict),
            }
            if strict:
                strict_frames[key].add(request["cluster_id"])
        levels[level] += 1
        by_population_filter[key][level] += 1
        evidence_rows.append({
            "local_index": index,
            "population": population,
            "target_cluster_id": request["target_cluster_id"],
            "cluster_id": request["cluster_id"],
            "control_id": request["control_id"],
            "primary_filter": band,
            "number": request["number"],
            "ra_deg": request["ra_deg"],
            "dec_deg": request["dec_deg"],
            "morphology_classification": request["classification"],
            "fwhm_moment_pix": request["fwhm_moment_pix"],
            "snr_nominal": request["snr_nominal"],
            "n_gaia_matches_within_2arcsec": len(matches),
            "evidence_level": level,
            **detail,
        })

    output.mkdir(parents=True, exist_ok=False)
    evidence_path = output / "candidate_evidence.tsv"
    batch_audit_path = output / "batch_audit.tsv"
    evidence_fields = list(evidence_rows[0])
    write_tsv(evidence_path, evidence_rows, evidence_fields)
    write_tsv(batch_audit_path, batch_audit_rows, list(batch_audit_rows[0]))
    summary = {
        "schema": "astrafek-s23b-psf-expansion-gaia-audit-v1",
        "status": "complete_local_diagnostic_only",
        "approved_plan_sha256": sha256(plan_path),
        "request_sha256": sha256(requests_path),
        "aggregate_response_sha256": sha256(response_path),
        "n_candidate_positions": len(requests),
        "n_target_candidate_positions": sum(row["population"] == "target" for row in requests),
        "n_control_candidate_positions": sum(row["population"] == "control" for row in requests),
        "n_completed_batches": len(batch_audit_rows),
        "n_gaia_match_rows": len(response_rows),
        "n_positions_with_gaia_match": sum(bool(rows) for rows in matches_by_index.values()),
        "max_separation_arcsec": max_separation,
        "evidence_levels": dict(levels),
        "by_population_filter": {
            f"{population}_{band}": {
                "candidate_positions": sum(
                    row["population"] == population and row["primary_filter"] == band
                    for row in requests
                ),
                **dict(counts),
                "strict_candidate_frames": len(strict_frames[(population, band)]),
            }
            for (population, band), counts in sorted(by_population_filter.items())
        },
        "astrometric_evidence_definition": {
            "cone_source_count": "exactly one Gaia DR3 source within 2 arcsec",
            "clean_fit": "RUWE < 1.4; duplicated_source=false; not flagged as Gaia QSO or galaxy candidate",
            "significant_motion": "parallax_over_error >= 3 OR two-component proper-motion significance >= 3",
            "strict_candidate_additional_condition": "Gaia combined star-class probability >= 0.5",
            "interpretation": "diagnostic thresholds only; not a calibrated star or PSF classifier",
        },
        "scope": {
            "images_requested": False,
            "hsc_shape_rows_requested": False,
            "photoz_requested": False,
            "cas_query_requested": False,
            "science_measurement": False,
            "psf_approved": False,
        },
        "inputs_sha256": {
            str(path.resolve()): sha256(path)
            for path in (plan_path, requests_path, position_path, status_path,
                         response_path, candidate_snapshot, TARGET_MANIFEST,
                         CONTROL_MANIFEST, TARGET_PREPARATION, CONTROL_PREPARATION)
        },
        "output_sha256": {
            evidence_path.name: sha256(evidence_path),
            batch_audit_path.name: sha256(batch_audit_path),
        },
        "notes": [
            "Gaia non-match is not evidence that an HSC candidate is non-stellar.",
            "A Gaia match does not by itself establish an isolated point source or PSF suitability.",
            "The 448 strict candidates are diagnostic inputs only; PSF approval remains false.",
            "The audit checks the approved query scope and reproduces every positional association from saved Gaia responses.",
            "HSC candidate positions are retained for image stamps; Gaia positions are not propagated to an unknown effective coadd epoch.",
        ],
        "science_measurement": False,
        "psf_approved": False,
    }
    write_json(output / "summary.json", summary)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    summary = audit(args.output.resolve())
    print(json.dumps({key: summary[key] for key in (
        "status", "approved_plan_sha256", "n_candidate_positions",
        "n_gaia_match_rows", "n_positions_with_gaia_match", "evidence_levels",
        "by_population_filter", "science_measurement", "psf_approved",
    )}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
