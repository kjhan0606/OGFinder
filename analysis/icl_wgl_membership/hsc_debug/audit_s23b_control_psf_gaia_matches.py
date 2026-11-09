#!/usr/bin/env python3
"""Audit a separately approved Gaia match for S23B control-field candidates.

This is a local-only response audit. It verifies the exact approved plan,
request, and Gaia response, then assigns conservative diagnostic evidence
levels. It does not approve a PSF or make a science measurement.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import Counter
from pathlib import Path

from match_psf_gaia_batched import GAIA_FIELDS


DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
PLAN_DIR = DATA / "s23b_psf_control_gaia_plan_2026-11-07_01"
DEFAULT_RUN = DATA / "s23b_psf_control_gaia_match_2026-11-07_01"
DEFAULT_OUTPUT = DATA / "s23b_psf_control_gaia_audit_2026-11-07_01"
EXPECTED_COLUMNS = ["candidate_id", *GAIA_FIELDS, "separation_arcsec"]


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


def audit(plan: Path, requests_path: Path, run: Path, output: Path) -> dict:
    if output.exists():
        raise FileExistsError(f"refusing to overwrite Gaia audit directory: {output}")
    status_path = run / "status.json"
    response_path = run / "gaia_matches.csv"
    for path in (plan, requests_path, status_path, response_path):
        if not path.is_file():
            raise FileNotFoundError(path)

    plan_doc = json.loads(plan.read_text(encoding="utf-8"))
    status = json.loads(status_path.read_text(encoding="utf-8"))
    requests_doc = json.loads(requests_path.read_text(encoding="utf-8"))
    candidates = requests_doc.get("candidates")
    if not isinstance(candidates, list) or not candidates:
        raise ValueError("approved request candidate list is empty or invalid")
    if (plan_doc.get("status") != "prepared_not_submitted"
            or status.get("status") != "complete"
            or status.get("plan_sha256") != sha256(plan)
            or status.get("request_sha256") != sha256(requests_path)
            or status.get("response_sha256") != sha256(response_path)
            or status.get("n_candidates") != len(candidates)
            or status.get("match_radius_arcsec") != 2.0
            or plan_doc.get("candidate_list_sha256") != sha256(requests_path)
            or plan_doc.get("n_new_candidate_positions") != len(candidates)
            or plan_doc.get("columns_returned") != EXPECTED_COLUMNS):
        raise ValueError("plan/request/response linkage or exact Gaia scope validation failed")
    for key in ("images_requested", "hsc_shape_rows_requested", "photoz_requested",
                "science_measurement", "psf_approved"):
        if plan_doc.get(key) is not False or status.get(key) is not False:
            raise ValueError(f"unexpected scope/state for {key}")

    with response_path.open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames != EXPECTED_COLUMNS:
            raise ValueError(f"unexpected Gaia response columns: {reader.fieldnames}")
        response_rows = list(reader)
    if len(response_rows) != status.get("n_match_rows"):
        raise ValueError("Gaia response row count differs from run status")

    request_by_index = {int(row["local_index"]): row for row in candidates}
    if sorted(request_by_index) != list(range(len(candidates))):
        raise ValueError("candidate local_index values are incomplete or duplicated")
    matches_by_index: dict[int, list[dict]] = {index: [] for index in range(len(candidates))}
    for row in response_rows:
        index = int(row["candidate_id"])
        if index not in matches_by_index:
            raise ValueError(f"response has an unknown candidate ID: {index}")
        separation = number(row["separation_arcsec"])
        if separation is None or separation > 2.00001:
            raise ValueError(f"Gaia response match falls outside the 2-arcsec cone: {index}")
        item = dict(row)
        item["_separation"] = separation
        item["_parallax_snr"] = number(row["parallax_over_error"])
        item["_pmra"] = number(row["pmra"])
        item["_pmra_error"] = number(row["pmra_error"])
        item["_pmdec"] = number(row["pmdec"])
        item["_pmdec_error"] = number(row["pmdec_error"])
        item["_ruwe"] = number(row["ruwe"])
        item["_star_probability"] = number(row["classprob_dsc_combmod_star"])
        item["_duplicated"] = flag(row["duplicated_source"])
        item["_qso"] = flag(row["in_qso_candidates"])
        item["_galaxy"] = flag(row["in_galaxy_candidates"])
        pm_terms = []
        if item["_pmra"] is not None and item["_pmra_error"] not in (None, 0.0):
            pm_terms.append((item["_pmra"] / item["_pmra_error"]) ** 2)
        if item["_pmdec"] is not None and item["_pmdec_error"] not in (None, 0.0):
            pm_terms.append((item["_pmdec"] / item["_pmdec_error"]) ** 2)
        item["_pm_significance"] = math.sqrt(sum(pm_terms)) if pm_terms else None
        matches_by_index[index].append(item)

    evidence_rows = []
    levels: Counter[str] = Counter()
    by_filter: dict[str, Counter[str]] = {}
    for index, request in enumerate(candidates):
        matches = sorted(matches_by_index[index], key=lambda row: (row["_separation"], row["source_id"]))
        band = request["primary_filter"]
        counts = by_filter.setdefault(band, Counter())
        counts["candidate_positions"] += 1
        if not matches:
            level = "no_gaia_match"
            best = {}
        else:
            counts["matched_positions"] += 1
            best = matches[0]
            clean_fit = (
                best["_ruwe"] is not None and best["_ruwe"] < 1.4
                and best["_duplicated"] is False
                and best["_qso"] is not True
                and best["_galaxy"] is not True
            )
            significant_motion = (
                (best["_parallax_snr"] is not None and best["_parallax_snr"] >= 3.0)
                or (best["_pm_significance"] is not None and best["_pm_significance"] >= 3.0)
            )
            strict = bool(
                len(matches) == 1 and clean_fit and significant_motion
                and best["_star_probability"] is not None
                and best["_star_probability"] >= 0.5
            )
            if len(matches) > 1:
                level = "ambiguous_cone"
            elif strict:
                level = "strict_candidate"
            elif clean_fit and significant_motion:
                level = "clean_astrometric_evidence"
            else:
                level = "matched_no_clean_astrometry"
        levels[level] += 1
        counts[level] += 1
        row = {
            "local_index": index,
            "target_cluster_id": request["target_cluster_id"],
            "control_id": request["cluster_id"],
            "number": request["number"],
            "primary_filter": band,
            "ra_deg": request["ra_deg"],
            "dec_deg": request["dec_deg"],
            "n_gaia_matches_within_2arcsec": len(matches),
            "evidence_level": level,
            "best_source_id": best.get("source_id"),
            "best_separation_arcsec": best.get("separation_arcsec"),
            "parallax_over_error": best.get("parallax_over_error"),
            "pm_total_mas_per_year": (
                math.hypot(float(best["pmra"]), float(best["pmdec"]))
                if best and number(best.get("pmra")) is not None and number(best.get("pmdec")) is not None
                else None
            ),
            "pm_significance": best.get("_pm_significance"),
            "ruwe": best.get("ruwe"),
            "duplicated_source": best.get("duplicated_source"),
            "in_qso_candidates": best.get("in_qso_candidates"),
            "in_galaxy_candidates": best.get("in_galaxy_candidates"),
            "star_probability": best.get("classprob_dsc_combmod_star"),
            "phot_g_mean_mag": best.get("phot_g_mean_mag"),
            "phot_variable_flag": best.get("phot_variable_flag"),
            "strict_candidate": level == "strict_candidate",
        }
        evidence_rows.append(row)

    output.mkdir(parents=True, exist_ok=False)
    evidence_path = output / "control_candidate_evidence.tsv"
    write_tsv(evidence_path, evidence_rows, list(evidence_rows[0]))
    summary = {
        "schema": "astrafek-s23b-control-psf-gaia-audit-v1",
        "status": "complete_local_diagnostic_only",
        "n_candidate_positions": len(candidates),
        "n_gaia_match_rows": len(response_rows),
        "evidence_levels": dict(levels),
        "by_filter": {band: dict(counts) for band, counts in by_filter.items()},
        "astrometric_evidence_definition": {
            "cone_source_count": "exactly one Gaia source within 2 arcsec",
            "clean_fit": "RUWE < 1.4; duplicated_source=false; not flagged as Gaia QSO or galaxy candidate",
            "significant_motion": "parallax_over_error >= 3 OR two-component proper-motion significance >= 3",
            "strict_candidate_additional_condition": "Gaia combined star-class probability >= 0.5",
            "status": "diagnostic thresholds only; not a calibrated star/PSF classifier",
        },
        "inputs_sha256": {
            str(path.resolve()): sha256(path)
            for path in (plan, requests_path, status_path, response_path)
        },
        "output_sha256": {evidence_path.name: sha256(evidence_path)},
        "science_measurement": False,
        "psf_approved": False,
        "notes": [
            "A Gaia non-match is not evidence that an HSC candidate is non-stellar.",
            "A Gaia match alone does not prove an isolated point source or approve a PSF.",
            "The only Gaia query audited here is the exact separately approved control-candidate plan.",
            "The effective HSC coadd epoch must be reviewed before propagating Gaia positions.",
        ],
    }
    write_json(output / "summary.json", summary)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", type=Path, required=True)
    parser.add_argument("--requests", type=Path, required=True)
    parser.add_argument("--run", type=Path, default=DEFAULT_RUN)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    summary = audit(args.plan.resolve(), args.requests.resolve(), args.run.resolve(), args.output.resolve())
    print(json.dumps({key: summary[key] for key in (
        "status", "n_candidate_positions", "n_gaia_match_rows",
        "evidence_levels", "by_filter", "science_measurement", "psf_approved",
    )}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
