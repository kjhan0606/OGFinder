#!/usr/bin/env python3
"""Prepare HSC random-point metadata queries for S23B blank-field controls.

Offline only. The plan is limited to random positions and scalar image-quality
flags for the 1,043 S23B candidates passing the nine-point coverage screen.
It requests no images, source detections, shapes, or photo-z.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path


DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
DEFAULT_COVERAGE_PLAN = DATA / "icl_s23b_coverage_plan_2026-11-06_04"
DEFAULT_COVERAGE_AUDIT = DATA / "icl_s23b_coverage_audit_2026-11-06_02"
DEFAULT_OUTPUT = DATA / "s23b_icl_control_candidate_plan_2026-11-06_01"
EXPECTED_COVERAGE_PLAN_SHA256 = "62c28446f553bf45a940ea0126af1521d4ee0c1356c12cc73f42712df482c0c7"
EXPECTED_COVERAGE_AUDIT_SHA256 = "9ad4bbdd84881e63e22c6e3dec21b87a5ade71f2978b1844d523bccd9d4030d2"
BATCH_SIZE_CANDIDATES = 10
MAX_RANDOM_PER_CANDIDATE = 200
SEARCH_RADIUS_ARCSEC = 1800.0
OUTPUT_COLUMNS = [
    "candidate_id", "random_object_id", "ra", "dec", "skymap_id", "tract", "patch", "patch_s",
    "isprimary", "detect_ispatchinner", "detect_istractinner", "i_pix_variance", "i_sky_mean",
    "i_sky_std", "i_inputcount_value", "i_inputcount_flag", "i_inputcount_flag_noinputs",
    "i_inputcount_flag_badcentroid", "i_pixelflags_offimage", "i_pixelflags_edge",
    "i_pixelflags_bad", "i_pixelflags_interpolatedcenter", "i_pixelflags_saturatedcenter",
    "i_pixelflags_crcenter", "i_pixelflags_suspectcenter", "i_pixelflags_bright_objectcenter",
    "icalexp_exists", "imag_psf_depth", "iseeing",
]


def sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def sql_for_candidate(candidate: dict[str, str], seed: str) -> str:
    ident = candidate["candidate_id"].replace("'", "''")
    ra, dec = float(candidate["ra_deg"]), float(candidate["dec_deg"])
    random_columns = [
        "r.object_id AS random_object_id", "r.ra", "r.dec", "r.skymap_id", "r.tract", "r.patch", "r.patch_s",
        "r.isprimary", "r.detect_ispatchinner", "r.detect_istractinner", "r.i_pix_variance", "r.i_sky_mean",
        "r.i_sky_std", "r.i_inputcount_value", "r.i_inputcount_flag", "r.i_inputcount_flag_noinputs",
        "r.i_inputcount_flag_badcentroid", "r.i_pixelflags_offimage", "r.i_pixelflags_edge", "r.i_pixelflags_bad",
        "r.i_pixelflags_interpolatedcenter", "r.i_pixelflags_saturatedcenter", "r.i_pixelflags_crcenter",
        "r.i_pixelflags_suspectcenter", "r.i_pixelflags_bright_objectcenter", "q.icalexp_exists",
        "q.imag_psf_depth", "q.iseeing",
    ]
    return (
        "SELECT * FROM (\n"
        f"  SELECT '{ident}' AS candidate_id, {', '.join(random_columns)}\n"
        "  FROM pdr3_wide.random AS r\n"
        "  LEFT JOIN pdr3_wide.patch_qa AS q USING (skymap_id)\n"
        f"  WHERE conesearch(r.coord, {ra:.10f}, {dec:.10f}, {SEARCH_RADIUS_ARCSEC:.1f})\n"
        "    AND r.i_pixelflags_bright_objectcenter IS FALSE\n"
        "    AND r.i_pixelflags_edge IS FALSE\n"
        f"  ORDER BY md5(r.object_id::text || ':{seed}')\n"
        f"  LIMIT {MAX_RANDOM_PER_CANDIDATE}\n"
        ") AS random_candidates"
    )


def batch_sql(candidates: list[dict[str, str]], batch_index: int) -> str:
    seed = f"S23B-control-20261106-{batch_index:03d}"
    return "\nUNION ALL\n".join(f"({sql_for_candidate(row, seed)})" for row in candidates) + ";\n"


def prepare(coverage_plan_dir: Path, audit_dir: Path, output_dir: Path) -> dict:
    coverage_plan_dir, audit_dir, output_dir = map(Path.resolve, (coverage_plan_dir, audit_dir, output_dir))
    if output_dir.exists():
        raise FileExistsError(f"refusing to overwrite {output_dir}")
    coverage_plan_path = coverage_plan_dir / "coverage_query_plan.json"
    audit_summary_path = audit_dir / "summary.json"
    candidate_path = coverage_plan_dir / "candidate_centers.tsv"
    coverage_table_path = audit_dir / "per_candidate_coverage.tsv"
    for path in (coverage_plan_path, audit_summary_path, candidate_path, coverage_table_path):
        if not path.is_file():
            raise FileNotFoundError(path)
    if sha(coverage_plan_path) != EXPECTED_COVERAGE_PLAN_SHA256:
        raise ValueError("S23B footprint plan differs from the completed approved plan")
    audit = json.loads(audit_summary_path.read_text(encoding="utf-8"))
    if EXPECTED_COVERAGE_AUDIT_SHA256 and sha(audit_summary_path) != EXPECTED_COVERAGE_AUDIT_SHA256:
        raise ValueError("S23B coverage-audit summary checksum differs from the reviewed audit")
    if (audit.get("schema") != "astrafek-s23b-icl-coverage-metadata-audit-v1"
            or audit.get("status") != "complete_exploratory_audit"
            or audit.get("plan_sha256") != EXPECTED_COVERAGE_PLAN_SHA256
            or audit.get("n_candidates_all_9_probes_covered") != 1043
            or audit.get("n_candidates_all_9_probes_with_imeas_exists") != 1043
            or audit.get("n_candidates_all_9_probes_with_idet_bkgd_exists") != 1043
            or audit.get("science_measurement") is not False):
        raise ValueError("S23B full-coverage candidate set changed")
    centers = read_tsv(candidate_path)
    coverage = read_tsv(coverage_table_path)
    by_id = {row["candidate_id"]: row for row in centers}
    if len(by_id) != 1333 or len(centers) != 1333 or len(coverage) != 1333:
        raise ValueError("S23B candidate inventory changed")
    eligible = [row for row in coverage
                if row["all_nine_probes_covered"] == "True"
                and row["all_nine_have_image_product_flag"] == "True"
                and row["all_nine_have_background_product_flag"] == "True"]
    if len(eligible) != 1043:
        raise ValueError("expected exactly 1,043 nine-point screened candidates")
    eligible.sort(key=lambda row: row["candidate_id"])
    output_dir.mkdir(parents=True, exist_ok=False)
    candidate_output = output_dir / "eligible_s23b_candidates.tsv"
    fields = list(centers[0]) + [
        "covered_probe_count", "image_product_probe_count", "background_product_probe_count",
        "min_imag_psf_depth", "median_imag_psf_depth", "median_seeing_arcsec",
    ]
    with candidate_output.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for quality in eligible:
            writer.writerow({**by_id[quality["candidate_id"]], **{
                key: quality[key] for key in fields if key not in by_id[quality["candidate_id"]]
            }})
    batches = []
    for start in range(0, len(eligible), BATCH_SIZE_CANDIDATES):
        chunk = eligible[start:start + BATCH_SIZE_CANDIDATES]
        name = f"random_control_candidates_{start // BATCH_SIZE_CANDIDATES:03d}.sql"
        sql_path = output_dir / name
        sql_path.write_text(batch_sql([by_id[row["candidate_id"]] for row in chunk], start // BATCH_SIZE_CANDIDATES), encoding="utf-8")
        batches.append({
            "name": sql_path.stem, "path": str(sql_path), "sha256": sha(sql_path),
            "n_candidates": len(chunk), "max_rows": len(chunk) * MAX_RANDOM_PER_CANDIDATE,
            "sql_bytes": sql_path.stat().st_size,
        })
    inputs = [coverage_plan_path, audit_summary_path, candidate_path, coverage_table_path]
    summary = {
        "schema": "astrafek-s23b-icl-control-candidate-query-plan-v1",
        "status": "prepared_pending_user_approval",
        "purpose": "retrieve candidate HSC random-point positions and limited scalar I-band QA metadata for matched-control screening",
        "release": "pdr3", "table": "pdr3_wide.random plus pdr3_wide.patch_qa",
        "n_target_candidates": len(eligible), "search_radius_arcsec": SEARCH_RADIUS_ARCSEC,
        "max_random_rows_per_candidate": MAX_RANDOM_PER_CANDIDATE,
        "maximum_returned_rows": len(eligible) * MAX_RANDOM_PER_CANDIDATE,
        "n_sql_batches": len(batches), "batch_size_candidates": BATCH_SIZE_CANDIDATES,
        "output_columns": OUTPUT_COLUMNS,
        "shape_or_photoz_rows_requested": False, "images_requested": False,
        "source_detections_requested": False, "credentials_saved": False,
        "cas_query_submitted": False, "science_measurement": False,
        "user_approval_required_before_submission": True,
        "approved_parent_metadata_plan_sha256": EXPECTED_COVERAGE_PLAN_SHA256,
        "input_sha256": {str(path): sha(path) for path in inputs},
        "candidate_table": str(candidate_output),
        "candidate_table_sha256": sha(candidate_output),
        "batches": batches,
        "estimated_csv_bytes": 80_000_000,
        "estimate_basis": "conservative estimate from existing approved random-table CSV row sizes, scaled to 208,600 rows and increased for joined QA columns",
        "data_contract": "Random-point coordinates/IDs and I-band flags, local pixel-variance/sky/input-count scalars, and patch depth/seeing/exposure availability only; no detected-source photometry, shape catalog, photo-z, or images.",
        "limitations": [
            "returned random points are candidates, not validated matched controls",
            "control matching, SFD extinction, CAMIRA contamination exclusion, angular-footprint coverage, and spatial dependence are evaluated locally after retrieval",
            "the 200-position cap may leave some target candidates without usable controls",
            "the 1-Mpc cluster-exclusion scale is a declared screening parameter, not a mass-derived radius",
        ],
    }
    summary_path = output_dir / "control_candidate_query_plan.json"
    summary["output_sha256"] = {str(candidate_output): sha(candidate_output), **{b["path"]: b["sha256"] for b in batches}}
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({k: summary[k] for k in (
        "status", "n_target_candidates", "maximum_returned_rows", "n_sql_batches",
        "output_columns", "estimated_csv_bytes", "science_measurement")}, indent=2))
    print("plan_sha256", sha(summary_path))
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--coverage-plan-dir", type=Path, default=DEFAULT_COVERAGE_PLAN)
    parser.add_argument("--coverage-audit-dir", type=Path, default=DEFAULT_COVERAGE_AUDIT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    prepare(args.coverage_plan_dir, args.coverage_audit_dir, args.output_dir)


if __name__ == "__main__":
    main()
