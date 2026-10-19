#!/usr/bin/env python3
"""Prepare a Gaia DR3 cross-match plan for unqueried HSC PSF candidates.

This is an offline preparation step. It reads existing HSC detection tables
and FITS WCS headers, excludes candidate detections already present in the
approved Gaia match request, and writes a coordinate-only candidate list for
review. It performs no network access.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
from pathlib import Path

from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS

import match_psf_gaia
import summarize_psf_candidate_coverage


DEFAULT_DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
DEFAULT_PREPARED = DEFAULT_DATA / "satellite_detection_psf_2026-10-18_01"
DEFAULT_DOWNLOAD_STATUS = DEFAULT_DATA / "download_01/download_status.json"
DEFAULT_PREVIOUS_CANDIDATES = DEFAULT_DATA / "psf_cds_2026-10-13_01/evidence/candidates.json"
DEFAULT_PREVIOUS_STATUS = DEFAULT_DATA / "psf_cds_2026-10-13_01/evidence/status.json"
DEFAULT_PREVIOUS_MATCHES = DEFAULT_DATA / "psf_cds_2026-10-13_01/evidence/gaia_matches.csv"
DEFAULT_COVERAGE = DEFAULT_DATA / "psf_candidate_coverage_2026-10-19_01/summary.json"
DEFAULT_OUTPUT = DEFAULT_DATA / "psf_gaia_full_sample_plan_2026-10-20_01"
MAX_MATCH_ROWS = 20_000


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def image_wcs(image: Path, science_hdu: int) -> WCS:
    with fits.open(image, memmap=True) as hdus:
        header = hdus[science_hdu].header.copy()
        if header.get("INHERIT") is True:
            header.extend(hdus[0].header, unique=True)
    wcs = WCS(header).celestial
    if not wcs.has_celestial:
        raise ValueError(f"image has no celestial WCS: {image}")
    return wcs


def prepare(prepared: Path, download_status_path: Path,
            previous_candidates_path: Path, previous_status_path: Path,
            previous_matches_path: Path, coverage_path: Path,
            output: Path) -> dict:
    if output.exists():
        raise FileExistsError(f"refusing to overwrite plan directory: {output}")

    coverage = summarize_psf_candidate_coverage.summarize(prepared)
    saved_coverage = json.loads(coverage_path.read_text(encoding="utf-8"))
    coverage_keys = (
        "status", "n_frames_verified", "n_candidate_rows", "by_filter",
        "preparation_summary_sha256", "independent_artifact_audit_sha256",
    )
    if any(saved_coverage.get(key) != coverage.get(key) for key in coverage_keys):
        raise ValueError("saved candidate-coverage summary does not match the audited tables")

    previous_status = json.loads(previous_status_path.read_text(encoding="utf-8"))
    previous = json.loads(previous_candidates_path.read_text(encoding="utf-8"))
    if previous_status.get("status") != "complete":
        raise ValueError("existing Gaia candidate cross-match is not complete")
    with previous_matches_path.open(encoding="utf-8-sig", newline="") as stream:
        previous_rows = list(csv.DictReader(stream))
    if (sha256(previous_matches_path) != previous_status.get("response_sha256")
            or len(previous) != previous_status.get("n_candidates")
            or len(previous_rows) != previous_status.get("n_match_rows")):
        raise ValueError("existing Gaia candidate cross-match failed its saved checksum/count checks")
    old_keys = {(str(row["cluster_id"]), int(row["number"])) for row in previous}
    if len(old_keys) != len(previous):
        raise ValueError("existing Gaia candidate list has duplicate cluster/detection keys")

    prep_summary_path = prepared / "summary.json"
    prep_summary = json.loads(prep_summary_path.read_text(encoding="utf-8"))
    if (prep_summary.get("status") != "complete"
            or prep_summary.get("n_prepared") != 800
            or prep_summary.get("science_measurement") is not False):
        raise ValueError("expected the audited 800-frame provisional candidate preparation")
    download_status = json.loads(download_status_path.read_text(encoding="utf-8"))
    downloaded = {str(row["id"]): row for row in download_status.get("completed", [])}
    records = {str(row["id"]): row for row in prep_summary.get("records", [])}
    if len(downloaded) != 800 or len(records) != 800 or set(downloaded) != set(records):
        raise ValueError("prepared candidate and downloaded-image manifests do not cover the same 800 frames")

    candidates = []
    seen_keys = set()
    candidate_frames = set()
    all_count = 0
    for cluster_id in sorted(records, key=lambda value: int(value)):
        record = records[cluster_id]
        file_record = downloaded[cluster_id]
        image = Path(record["source_image"]).resolve()
        expected_image = (DEFAULT_DATA / "download_01/images" / file_record["file"]).resolve()
        if image != expected_image or not image.is_file():
            raise ValueError(f"prepared source image differs from the download manifest for {cluster_id}")
        science_hdu = int(file_record["review"]["science_hdu"])
        wcs = image_wcs(image, science_hdu)
        table_path = prepared / f"cluster_{cluster_id}" / "psf_candidates.tsv"
        with table_path.open(encoding="utf-8", newline="") as stream:
            rows = list(csv.DictReader(stream, delimiter="\t"))
        if len(rows) != int(record["n_psf_candidates"]):
            raise ValueError(f"candidate count differs from preparation summary: {table_path}")
        for row in rows:
            if row.get("CLASSIFICATION") != "provisional_compact_candidate":
                raise ValueError(f"unexpected candidate classification: {table_path}")
            number = int(row["NUMBER"])
            key = (cluster_id, number)
            if key in seen_keys:
                raise ValueError(f"duplicate candidate identity: {key}")
            seen_keys.add(key)
            all_count += 1
            if key in old_keys:
                continue
            candidate_frames.add(cluster_id)
            x, y = float(row["X_IMAGE"]) - 1.0, float(row["Y_IMAGE"]) - 1.0
            sky = wcs.pixel_to_world(x, y)
            ra, dec = float(sky.ra.deg), float(sky.dec.deg)
            if not (math.isfinite(ra) and 0 <= ra < 360 and math.isfinite(dec) and -90 <= dec <= 90):
                raise ValueError(f"invalid WCS coordinate for candidate {key}")
            candidates.append({
                "cluster_id": cluster_id,
                "number": number,
                "ra_deg": ra,
                "dec_deg": dec,
                "classification": "provisional_compact_candidate",
                "fwhm_moment_pix": float(row["FWHM_MOMENT_PIX"]),
                "snr_nominal": float(row["SNR_NOMINAL"]),
                "wing_valid_fraction": float(row["WING_VALID_FRACTION"]),
                "required_evidence": [
                    "Gaia DR3 point-source cross-match",
                    "astrometric quality and epoch review",
                    "isolation and extended-wing image review",
                ],
            })

    if all_count != 8078 or len(old_keys - seen_keys) != 0:
        raise ValueError(f"candidate set changed: found {all_count}, prior unmatched={len(old_keys-seen_keys)}")
    if len(candidates) != all_count - len(old_keys):
        raise ValueError("new candidate list did not exclude exactly the prior Gaia query inputs")
    previous_frames = {str(row["cluster_id"]) for row in previous}
    votable = io.BytesIO()
    Table(dict(candidate_id=list(range(len(candidates))),
               ra=[row["ra_deg"] for row in candidates],
               dec=[row["dec_deg"] for row in candidates])).write(votable, format="votable")
    votable_bytes = votable.getvalue()
    estimated_rows = round(len(candidates) * len(previous_rows) / len(old_keys))
    estimated_result_bytes = round(estimated_rows * previous_matches_path.stat().st_size / len(previous_rows))
    output.mkdir(parents=True, exist_ok=False)
    candidate_path = output / "candidates.json"
    candidate_path.write_text(json.dumps({"candidates": candidates}, indent=2) + "\n", encoding="utf-8")
    query_path = output / "query.adql"
    query_path.write_text(match_psf_gaia.QUERY + "\n", encoding="utf-8")
    plan = {
        "schema": "astrafek-psf-gaia-candidate-plan-v1",
        "status": "prepared_not_submitted",
        "external_query": "ESA Gaia DR3 TAP, positional cross-match only",
        "match_radius_arcsec": 2.0,
        "n_hsc_frames": 800,
        "n_provisional_candidates": all_count,
        "n_candidates_in_previous_gaia_query": len(old_keys),
        "n_new_candidate_positions": len(candidates),
        "n_frames_with_new_candidates": len(candidate_frames),
        "n_frames_in_previous_gaia_query": len(previous_frames),
        "n_frames_covered_after_merge": len(candidate_frames | previous_frames),
        "previous_gaia_match_rows": len(previous_rows),
        "expected_match_rows_rough_estimate": estimated_rows,
        "estimated_response_bytes_rough": estimated_result_bytes,
        "coordinate_votable_upload_bytes": len(votable_bytes),
        "max_match_rows_before_truncation_guard": MAX_MATCH_ROWS,
        "columns_returned": [
            "candidate_id", "source_id", "ra", "dec", "ref_epoch",
            "parallax_over_error", "pmra", "pmra_error", "pmdec", "pmdec_error",
            "ruwe", "duplicated_source", "phot_g_mean_mag", "phot_variable_flag",
            "in_qso_candidates", "in_galaxy_candidates", "classprob_dsc_combmod_star",
            "separation_arcsec",
        ],
        "images_requested": False,
        "hsc_shape_rows_requested": False,
        "photoz_requested": False,
        "science_measurement": False,
        "previous_gaia_candidate_list_sha256": sha256(previous_candidates_path),
        "previous_gaia_status_sha256": sha256(previous_status_path),
        "previous_gaia_matches_sha256": sha256(previous_matches_path),
        "candidate_coverage_summary_sha256": sha256(coverage_path),
        "preparation_summary_sha256": sha256(prep_summary_path),
        "download_status_sha256": sha256(download_status_path),
        "candidate_list": str(candidate_path.resolve()),
        "candidate_list_sha256": sha256(candidate_path),
        "query": str(query_path.resolve()),
        "query_sha256": sha256(query_path),
        "input_hsc_images_sha256": {
            cluster_id: downloaded[cluster_id]["sha256"] for cluster_id in sorted(downloaded, key=int)
        },
        "candidate_tables_sha256": {
            cluster_id: sha256(prepared / f"cluster_{cluster_id}" / "psf_candidates.tsv")
            for cluster_id in sorted(records, key=int)
        },
        "psf_approved": False,
    }
    (output / "plan.json").write_text(json.dumps(plan, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return plan


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prepared", type=Path, default=DEFAULT_PREPARED)
    parser.add_argument("--download-status", type=Path, default=DEFAULT_DOWNLOAD_STATUS)
    parser.add_argument("--previous-candidates", type=Path, default=DEFAULT_PREVIOUS_CANDIDATES)
    parser.add_argument("--previous-status", type=Path, default=DEFAULT_PREVIOUS_STATUS)
    parser.add_argument("--previous-matches", type=Path, default=DEFAULT_PREVIOUS_MATCHES)
    parser.add_argument("--coverage-summary", type=Path, default=DEFAULT_COVERAGE)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    result = prepare(args.prepared.resolve(), args.download_status.resolve(),
                     args.previous_candidates.resolve(), args.previous_status.resolve(),
                     args.previous_matches.resolve(), args.coverage_summary.resolve(),
                     args.output.resolve())
    print(json.dumps({k: result[k] for k in (
        "status", "n_hsc_frames", "n_provisional_candidates",
        "n_candidates_in_previous_gaia_query", "n_new_candidate_positions",
        "previous_gaia_match_rows", "expected_match_rows_rough_estimate",
        "match_radius_arcsec", "candidate_list_sha256", "query_sha256",
    )}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
