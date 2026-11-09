#!/usr/bin/env python3
"""Audit an approved Gaia DR3 match for S23B morphology-only PSF candidates.

This is a local, read-only audit of saved Gaia response rows, candidate tables,
and HSC FITS headers. It performs no network access and does not approve a PSF
or promote a science measurement.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

from astropy.io import fits


DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
PILOT = DATA / "s23b_detection_psf_pilot_2026-11-07_01"
PLAN_DIR = DATA / "s23b_psf_gaia_plan_2026-11-07_01"
RUN_DIR = DATA / "s23b_psf_gaia_match_2026-11-07_01"
DEFAULT_OUTPUT = DATA / "s23b_psf_gaia_audit_2026-11-07_01"

GAIA_FIELDS = [
    "source_id", "ra", "dec", "ref_epoch", "parallax_over_error",
    "pmra", "pmra_error", "pmdec", "pmdec_error", "ruwe",
    "duplicated_source", "phot_g_mean_mag", "phot_variable_flag",
    "in_qso_candidates", "in_galaxy_candidates",
    "classprob_dsc_combmod_star",
]
MATCH_FIELDS = ["candidate_id", *GAIA_FIELDS, "separation_arcsec"]
TIME_KEY_TOKENS = ("DATE", "MJD", "EPOCH")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def number(value: str | None) -> float | None:
    if value in (None, ""):
        return None
    parsed = float(value)
    return parsed if math.isfinite(parsed) else None


def flag(value: str | None) -> bool | None:
    if value in (None, ""):
        return None
    lowered = value.strip().lower()
    if lowered in ("true", "1", "t"):
        return True
    if lowered in ("false", "0", "f"):
        return False
    raise ValueError(f"unexpected boolean field: {value!r}")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def write_json(path: Path, payload: dict) -> None:
    with path.open("x", encoding="utf-8") as stream:
        json.dump(payload, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")


def audit(output: Path) -> dict:
    if output.exists():
        raise FileExistsError(f"refusing to overwrite audit directory: {output}")

    plan_path = PLAN_DIR / "plan.json"
    requests_path = PLAN_DIR / "requests.json"
    run_status_path = RUN_DIR / "status.json"
    response_path = RUN_DIR / "gaia_matches.csv"
    pilot_plan_path = PILOT / "plan.json"
    pilot_manifest_path = PILOT / "manifest.json"
    prep_summary_path = PILOT / "preparation/summary.json"

    for path in (
        plan_path, requests_path, run_status_path, response_path,
        pilot_plan_path, pilot_manifest_path, prep_summary_path,
    ):
        if not path.is_file():
            raise FileNotFoundError(path)

    plan = json.loads(plan_path.read_text(encoding="utf-8"))
    requests_doc = json.loads(requests_path.read_text(encoding="utf-8"))
    requests = requests_doc.get("candidates")
    run_status = json.loads(run_status_path.read_text(encoding="utf-8"))
    pilot_plan = json.loads(pilot_plan_path.read_text(encoding="utf-8"))
    pilot_manifest = json.loads(pilot_manifest_path.read_text(encoding="utf-8"))
    prep_summary = json.loads(prep_summary_path.read_text(encoding="utf-8"))

    if not isinstance(requests, list) or len(requests) != 492:
        raise ValueError("the hash-pinned Gaia request must contain 492 candidate positions")
    if plan.get("status") != "prepared_not_submitted":
        raise ValueError("saved plan status is not the approved, immutable prepared plan")
    if run_status.get("status") != "complete":
        raise ValueError("Gaia query did not complete")
    if run_status.get("plan_sha256") != sha256(plan_path):
        raise ValueError("run status does not pin the current plan SHA-256")
    if run_status.get("request_sha256") != sha256(requests_path):
        raise ValueError("run status does not pin the current request SHA-256")
    if run_status.get("response_sha256") != sha256(response_path):
        raise ValueError("run status does not pin the current Gaia response SHA-256")
    if run_status.get("n_candidates") != len(requests):
        raise ValueError("Gaia candidate count differs from the approved request")
    if run_status.get("match_radius_arcsec") != 2.0:
        raise ValueError("Gaia run radius differs from the approved plan")
    if any(run_status.get(key) is not False for key in (
        "images_requested", "hsc_shape_rows_requested", "photoz_requested",
        "science_measurement", "psf_approved",
    )):
        raise ValueError("the Gaia run flags exceed the approved query scope")
    if prep_summary.get("science_measurement") is not False:
        raise ValueError("candidate-preparation run unexpectedly promotes science")

    local_indices = [int(row["local_index"]) for row in requests]
    if sorted(local_indices) != list(range(len(requests))):
        raise ValueError("request local_index values are incomplete or duplicated")
    request_by_index = {int(row["local_index"]): row for row in requests}
    selection_by_cluster = {
        row["id"]: row["primary_filter"] for row in pilot_plan["selection"]
    }
    if len(selection_by_cluster) != 50:
        raise ValueError("pilot selection no longer contains 50 frames")

    prep_records = {row["id"]: row for row in prep_summary.get("records", [])}
    gaia_rows: list[dict[str, str]] = []
    with response_path.open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream)
        if reader.fieldnames != MATCH_FIELDS:
            raise ValueError(f"unexpected Gaia response columns: {reader.fieldnames}")
        gaia_rows.extend(reader)
    if len(gaia_rows) != run_status.get("n_match_rows"):
        raise ValueError("Gaia response row count differs from completed run status")

    matches_by_index: dict[int, list[dict]] = {i: [] for i in range(len(requests))}
    source_to_candidate_indices: dict[str, set[int]] = {}
    for row in gaia_rows:
        index = int(row["candidate_id"])
        if index not in request_by_index:
            raise ValueError(f"response has unknown local candidate index {index}")
        separation = number(row["separation_arcsec"])
        if separation is None or separation > 2.00001:
            raise ValueError("a returned Gaia match falls outside the approved 2-arcsec cone")
        enriched = dict(row)
        enriched["_separation"] = separation
        enriched["_parallax_over_error"] = number(row["parallax_over_error"])
        enriched["_pmra"] = number(row["pmra"])
        enriched["_pmra_error"] = number(row["pmra_error"])
        enriched["_pmdec"] = number(row["pmdec"])
        enriched["_pmdec_error"] = number(row["pmdec_error"])
        enriched["_ruwe"] = number(row["ruwe"])
        enriched["_star_probability"] = number(row["classprob_dsc_combmod_star"])
        enriched["_duplicated"] = flag(row["duplicated_source"])
        enriched["_qso"] = flag(row["in_qso_candidates"])
        enriched["_galaxy"] = flag(row["in_galaxy_candidates"])
        pm_significance_terms = []
        if enriched["_pmra"] is not None and enriched["_pmra_error"] not in (None, 0.0):
            pm_significance_terms.append((enriched["_pmra"] / enriched["_pmra_error"]) ** 2)
        if enriched["_pmdec"] is not None and enriched["_pmdec_error"] not in (None, 0.0):
            pm_significance_terms.append((enriched["_pmdec"] / enriched["_pmdec_error"]) ** 2)
        enriched["_pm_significance"] = math.sqrt(sum(pm_significance_terms)) if pm_significance_terms else None
        enriched["_pm_total"] = (
            math.hypot(enriched["_pmra"], enriched["_pmdec"])
            if enriched["_pmra"] is not None and enriched["_pmdec"] is not None
            else None
        )
        matches_by_index[index].append(enriched)
        source_to_candidate_indices.setdefault(row["source_id"], set()).add(index)

    for rows in matches_by_index.values():
        rows.sort(key=lambda row: (row["_separation"], row["source_id"]))

    candidate_rows = []
    evidence_levels: dict[str, int] = {
        "no_gaia_match": 0,
        "ambiguous_cone": 0,
        "matched_no_clean_astrometry": 0,
        "clean_astrometric_evidence": 0,
        "strict_candidate": 0,
    }
    duplicate_source_match_rows = sum(
        len(indices) > 1 for indices in source_to_candidate_indices.values()
    )
    match_counts_by_filter: dict[str, dict[str, int]] = {}

    for index in range(len(requests)):
        request = request_by_index[index]
        cluster_id = request["cluster_id"]
        band = selection_by_cluster.get(cluster_id)
        if band not in ("i", "i2"):
            raise ValueError(f"candidate references an unselected frame: {cluster_id}")
        band_counts = match_counts_by_filter.setdefault(band, {"candidate_positions": 0, "matched_positions": 0, "match_rows": 0})
        band_counts["candidate_positions"] += 1
        matches = matches_by_index[index]
        band_counts["match_rows"] += len(matches)
        if matches:
            band_counts["matched_positions"] += 1
            best = matches[0]
            ppar = best["_parallax_over_error"]
            pm_sig = best["_pm_significance"]
            pm_total = best["_pm_total"]
            clean_fit = (
                best["_ruwe"] is not None and best["_ruwe"] < 1.4
                and best["_duplicated"] is False
                and best["_qso"] is not True
                and best["_galaxy"] is not True
            )
            significant_motion = (
                (ppar is not None and ppar >= 3.0)
                or (pm_sig is not None and pm_sig >= 3.0)
            )
            astrometric_evidence = bool(len(matches) == 1 and clean_fit and significant_motion)
            strict = bool(
                astrometric_evidence
                and best["_star_probability"] is not None
                and best["_star_probability"] >= 0.5
            )
            if len(matches) > 1:
                level = "ambiguous_cone"
            elif strict:
                level = "strict_candidate"
            elif astrometric_evidence:
                level = "clean_astrometric_evidence"
            else:
                level = "matched_no_clean_astrometry"
            shifts = {
                str(years): (pm_total * years / 1000.0 if pm_total is not None else None)
                for years in (2, 5, 10)
            }
            best_record = {
                "best_source_id": best["source_id"],
                "best_separation_arcsec": best["_separation"],
                "parallax_over_error": ppar,
                "pm_total_mas_per_year": pm_total,
                "pm_significance": pm_sig,
                "pm_shift_arcsec_at_2yr": shifts["2"],
                "pm_shift_arcsec_at_5yr": shifts["5"],
                "pm_shift_arcsec_at_10yr": shifts["10"],
                "ruwe": best["_ruwe"],
                "duplicated_source": best["_duplicated"],
                "in_qso_candidates": best["_qso"],
                "in_galaxy_candidates": best["_galaxy"],
                "star_probability": best["_star_probability"],
                "phot_g_mean_mag": number(best["phot_g_mean_mag"]),
                "phot_variable_flag": best["phot_variable_flag"],
                "astrometric_star_evidence": astrometric_evidence,
                "strict_candidate": strict,
            }
        else:
            level = "no_gaia_match"
            best_record = {
                "best_source_id": None,
                "best_separation_arcsec": None,
                "parallax_over_error": None,
                "pm_total_mas_per_year": None,
                "pm_significance": None,
                "pm_shift_arcsec_at_2yr": None,
                "pm_shift_arcsec_at_5yr": None,
                "pm_shift_arcsec_at_10yr": None,
                "ruwe": None,
                "duplicated_source": None,
                "in_qso_candidates": None,
                "in_galaxy_candidates": None,
                "star_probability": None,
                "phot_g_mean_mag": None,
                "phot_variable_flag": None,
                "astrometric_star_evidence": False,
                "strict_candidate": False,
            }
        evidence_levels[level] += 1
        candidate_rows.append({
            "local_index": index,
            "cluster_id": cluster_id,
            "number": request["number"],
            "primary_filter": band,
            "ra_deg": request["ra_deg"],
            "dec_deg": request["dec_deg"],
            "fwhm_moment_pix": request["fwhm_moment_pix"],
            "snr_nominal": request["snr_nominal"],
            "n_gaia_matches_within_2arcsec": len(matches),
            "evidence_level": level,
            **best_record,
        })

    output.mkdir(parents=True, exist_ok=False)
    evidence_path = output / "candidate_evidence.tsv"
    evidence_fields = list(candidate_rows[0])
    with evidence_path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=evidence_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(candidate_rows)

    frame_headers = []
    time_key_count: dict[str, int] = {}
    frame_by_id = {row["id"]: row for row in pilot_manifest["clusters"]}
    for cluster_id in sorted(selection_by_cluster, key=lambda value: int(value.split("_")[1])):
        frame = frame_by_id[cluster_id]
        fits_path = Path(frame["image"])
        hdu_index = int(frame["science_hdu"])
        temporal_values = {}
        with fits.open(fits_path, memmap=True, do_not_scale_image_data=True, lazy_load_hdus=True) as hdus:
            headers = [hdus[0].header, hdus[hdu_index].header]
            for header in headers:
                for key in header.keys():
                    upper = str(key).upper()
                    if any(token in upper for token in TIME_KEY_TOKENS) and key not in temporal_values:
                        value = header.get(key)
                        if value is not None and isinstance(value, (str, int, float)):
                            temporal_values[str(key)] = value
        for key in temporal_values:
            time_key_count[key] = time_key_count.get(key, 0) + 1
        frame_headers.append({"cluster_id": cluster_id, "temporal_header_values": temporal_values})

    time_audit_path = output / "coadd_epoch_header_audit.json"
    write_json(time_audit_path, {
        "schema": "astrafek-s23b-coadd-temporal-header-audit-v1",
        "n_frames": len(frame_headers),
        "n_frames_with_any_date_mjd_or_epoch_key": sum(bool(row["temporal_header_values"]) for row in frame_headers),
        "keyword_frame_counts": time_key_count,
        "frames": frame_headers,
        "note": "Header-only inspection; no FITS science pixels read. A missing epoch does not identify the actual coadd observation dates.",
    })

    summary = {
        "schema": "astrafek-s23b-psf-gaia-audit-v1",
        "status": "complete",
        "science_measurement": False,
        "psf_approved": False,
        "n_candidate_positions": len(requests),
        "n_gaia_match_rows": len(gaia_rows),
        "n_positions_with_at_least_one_match": sum(bool(rows) for rows in matches_by_index.values()),
        "n_positions_with_no_match": sum(not rows for rows in matches_by_index.values()),
        "n_positions_with_multiple_gaia_sources": sum(len(rows) > 1 for rows in matches_by_index.values()),
        "n_unique_gaia_source_ids": len(source_to_candidate_indices),
        "n_gaia_sources_associated_with_multiple_hsc_candidates": duplicate_source_match_rows,
        "match_rate_of_all_candidate_positions": sum(bool(rows) for rows in matches_by_index.values()) / len(requests),
        "evidence_levels": evidence_levels,
        "by_filter": match_counts_by_filter,
        "astrometric_evidence_definition": {
            "cone_source_count": "exactly one Gaia source within 2 arcsec",
            "clean_fit": "RUWE < 1.4; duplicated_source=false; not flagged as Gaia QSO or galaxy candidate",
            "significant_motion": "parallax_over_error >= 3 OR two-component proper-motion significance >= 3",
            "strict_candidate_additional_condition": "Gaia combined star-class probability >= 0.5",
            "status": "diagnostic thresholds only; not a calibrated star/PSF classifier",
        },
        "proper_motion_epoch_sensitivity": {
            "reference_epoch_year": 2016.0,
            "time_offsets_tested_years": [2, 5, 10],
            "n_matches_with_gt_0p5_arcsec_shift_at_5yr": sum(
                (number(row["pm_shift_arcsec_at_5yr"]) or 0.0) > 0.5 for row in candidate_rows
            ),
            "n_matches_with_gt_2arcsec_shift_at_5yr": sum(
                (number(row["pm_shift_arcsec_at_5yr"]) or 0.0) > 2.0 for row in candidate_rows
            ),
            "interpretation": "scenario-only; HSC coadd observation epoch is not assumed",
        },
        "inputs_sha256": {
            str(plan_path): sha256(plan_path),
            str(requests_path): sha256(requests_path),
            str(run_status_path): sha256(run_status_path),
            str(response_path): sha256(response_path),
            str(pilot_plan_path): sha256(pilot_plan_path),
            str(pilot_manifest_path): sha256(pilot_manifest_path),
            str(prep_summary_path): sha256(prep_summary_path),
        },
        "output_sha256": {
            "candidate_evidence.tsv": sha256(evidence_path),
            "coadd_epoch_header_audit.json": sha256(time_audit_path),
        },
        "notes": [
            "A Gaia non-match is not evidence that an HSC compact candidate is not a star; Gaia is shallower than HSC.",
            "A Gaia match does not alone prove that the HSC source is an isolated point source.",
            "No S23B candidate PSF is approved; target/control PSF matching and extended wings remain uncalibrated.",
        ],
    }
    write_json(output / "summary.json", summary)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    summary = audit(args.output.resolve())
    print(json.dumps({
        key: summary[key] for key in (
            "status", "n_candidate_positions", "n_gaia_match_rows",
            "n_positions_with_at_least_one_match", "n_positions_with_no_match",
            "n_positions_with_multiple_gaia_sources", "evidence_levels",
            "by_filter", "psf_approved", "science_measurement",
        )
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
