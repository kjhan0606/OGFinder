#!/usr/bin/env python3
"""Audit whether S23B-flagged control pairs enter the ICL--satellite sample.

This reads only saved JSON summaries and TSV tables. It performs no HSC query,
download, FITS access, image-pixel measurement, or association re-fit.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path


ROOT = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
DEFAULT_STACK = ROOT / "icl_2d_s23b_control_exclusion_2026-10-22_01"
DEFAULT_ICL = ROOT / "icl_annulus_support_2026-10-18_01"
DEFAULT_ASSOCIATION = ROOT / "satellite_outer_icl_association_perbin_2026-10-18_03"
DEFAULT_SATELLITE = ROOT / "satellite_outer_profile_perbin_full_2026-10-18_01"
DEFAULT_OUTPUT = ROOT / "icl_s23b_association_overlap_2026-10-22_01/summary.json"


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


def build_records(audit_matches: list[dict], stack_summary: dict,
                  icl_rows: list[dict[str, str]], association_rows: list[dict[str, str]],
                  satellite_rows: list[dict[str, str]]) -> list[dict]:
    effective_excluded = {str(row["target_id"]): row
                          for row in stack_summary.get("s23b_candidate_control_exclusions", [])}
    preexcluded = {str(row["target_id"]): row
                   for row in stack_summary.get("candidate_controls_already_filter_excluded", [])}
    if len(effective_excluded) != 6 or len(preexcluded) != 1:
        raise ValueError("stack summary must distinguish six new exclusions and one existing filter exclusion")
    icl_by_key = {}
    for row in icl_rows:
        if row.get("radial_bin_kpc") != "200-300":
            continue
        key = (str(row.get("id", "")), str(row.get("filter", "")))
        if key in icl_by_key:
            raise ValueError(f"duplicate ICL 200-300 kpc support row: {key}")
        icl_by_key[key] = row
    association_by_key = {}
    for row in association_rows:
        key = (str(row.get("cluster_id", "")), str(row.get("filter", "")))
        if key in association_by_key:
            raise ValueError(f"duplicate ICL--satellite association row: {key}")
        association_by_key[key] = row
    satellite_count_by_key = {}
    for row in satellite_rows:
        if row.get("radial_bin_kpc") != "50-75":
            continue
        key = (str(row.get("cluster_id", "")), str(row.get("filter", "")))
        satellite_count_by_key[key] = satellite_count_by_key.get(key, 0) + 1

    records, seen_targets, seen_candidates = [], set(), set()
    for match in audit_matches:
        target_id = str(match["target_cluster_id"])
        control_id = str(match["frame_id"])
        band = str(match["filter"])
        if target_id in seen_targets:
            raise ValueError(f"duplicate candidate-control overlap for target {target_id}")
        seen_targets.add(target_id)
        candidate_ids = [str(row["candidate_id"]) for row in match["candidates"]]
        if not candidate_ids or any(candidate in seen_candidates for candidate in candidate_ids):
            raise ValueError("candidate IDs are missing or repeated across audited control fields")
        seen_candidates.update(candidate_ids)
        if target_id in effective_excluded:
            exclusion = effective_excluded[target_id]
            if exclusion["control_id"] != control_id or exclusion["filter"] != band:
                raise ValueError("effective stack exclusion differs from the WCS candidate audit")
            baseline_eligible = True
        elif target_id in preexcluded:
            exclusion = preexcluded[target_id]
            if exclusion["control_id"] != control_id or exclusion["filter"] != band:
                raise ValueError("pre-existing filter exclusion differs from the WCS candidate audit")
            baseline_eligible = False
        else:
            raise ValueError(f"WCS overlap is absent from both stack exclusion sets: {target_id}")
        icl_row = icl_by_key.get((target_id, band))
        if icl_row is not None and str(icl_row.get("control_id")) != control_id:
            raise ValueError("ICL support uses a different control than the candidate WCS audit")
        association_row = association_by_key.get((target_id, band))
        satellite_count = satellite_count_by_key.get((target_id, band), 0)
        records.append({
            "target_id": target_id, "control_id": control_id, "filter": band,
            "candidate_ids": candidate_ids,
            "baseline_pair_eligible": baseline_eligible,
            "excluded_from_alternate_stack": target_id in effective_excluded,
            "has_icl_200_300_support": icl_row is not None,
            "has_icl_satellite_association_row": association_row is not None,
            "n_satellite_profiles_50_75": satellite_count,
            "has_satellite_50_75_profiles": satellite_count > 0,
            "icl_delta_200_300_njy_arcsec2": (
                float(icl_row["delta_sb_njy_arcsec2"]) if icl_row is not None else None),
            "satellite_delta_50_75_njy_arcsec2": (
                float(association_row["satellite_delta_50_75"])
                if association_row is not None else None),
        })
    if len(records) != 7 or len(seen_candidates) != 7:
        raise ValueError("expected exactly seven unique audited S23B candidate overlaps")
    return records


def run(stack_dir: Path, icl_dir: Path, association_dir: Path,
        satellite_dir: Path) -> dict:
    stack_dir, icl_dir, association_dir, satellite_dir = map(
        Path.resolve, (stack_dir, icl_dir, association_dir, satellite_dir))
    stack_summary_path = stack_dir / "summary.json"
    stack = read_json(stack_summary_path)
    if (stack.get("schema") != "astrafek-hsc-icl-target-control-2d-stack-v1"
            or stack.get("status") != "complete_exploratory"
            or stack.get("science_measurement") is not False
            or stack.get("n_pairs_stacked") != 777):
        raise ValueError("candidate-exclusion stack is incomplete or outside its reviewed scope")
    audit_ref = stack.get("candidate_control_audit") or {}
    audit_path = Path(audit_ref.get("path", ""))
    if (not audit_path.is_file() or audit_ref.get("schema")
            != "astrafek-s23b-candidate-existing-cutout-support-v1"
            or sha256(audit_path) != audit_ref.get("sha256")
            or stack.get("input_sha256", {}).get(str(audit_path)) != audit_ref.get("sha256")):
        raise ValueError("candidate WCS audit path or hash is invalid")
    audit = read_json(audit_path)
    if (audit.get("status") != "header_wcs_audit_complete"
            or audit.get("science_measurement") is not False
            or audit.get("image_pixels_read") is not False
            or audit.get("n_candidates") != 1333
            or audit.get("candidate_table_sha256")
            != "1b59ca9941849fa7970c43c6e0374db6587716fc4ae5d71b11f84a245ebe6b2c"
            or audit.get("target", {}).get("n_unique_candidates_inside") != 0):
        raise ValueError("candidate WCS audit scope is invalid")

    icl_summary_path = icl_dir / "summary.json"
    icl_summary = read_json(icl_summary_path)
    icl_table = icl_dir / "annulus_pair_profiles.tsv"
    if (icl_summary.get("schema") != "hsc-icl-per-annulus-support-v1"
            or icl_summary.get("status") != "complete"
            or icl_summary.get("science_gates", {}).get("science_measurement") is not False
            or icl_summary.get("output_sha256", {}).get(str(icl_table)) != sha256(icl_table)):
        raise ValueError("ICL support table failed its saved provenance check")

    assoc_summary_path = association_dir / "summary.json"
    assoc_summary = read_json(assoc_summary_path)
    assoc_table = association_dir / "matched_cluster_values.tsv"
    if (assoc_summary.get("schema") != "hsc-satellite-outer-icl-correlation-v1"
            or assoc_summary.get("status") != "complete"
            or assoc_summary.get("science_measurement") is not False
            or assoc_summary.get("outputs_sha256", {}).get(assoc_table.name) != sha256(assoc_table)):
        raise ValueError("ICL--satellite table failed its saved provenance check")
    satellite_summary_path = satellite_dir / "summary.json"
    satellite_summary = read_json(satellite_summary_path)
    satellite_table = satellite_dir / "satellite_profiles.tsv"
    if (satellite_summary.get("status") != "complete"
            or satellite_summary.get("science_measurement") is not False
            or satellite_summary.get("output_sha256", {}).get(str(satellite_table))
            != sha256(satellite_table)):
        raise ValueError("satellite 50--75 kpc table failed its saved provenance check")
    association_rows = read_tsv(assoc_table)
    if len(association_rows) != 94:
        raise ValueError("the hash-pinned ICL--satellite association table must contain 94 rows")
    records = build_records(audit["matched_controls"]["matches"], stack,
                            read_tsv(icl_table), association_rows,
                            read_tsv(satellite_table))
    return {
        "schema": "astrafek-s23b-icl-satellite-support-overlap-v1",
        "status": "verified_descriptive_audit",
        "science_measurement": False,
        "analysis_script": str(Path(__file__).resolve()),
        "analysis_script_sha256": sha256(Path(__file__).resolve()),
        "n_candidate_control_overlaps": len(records),
        "n_baseline_filter_eligible_pairs": sum(r["baseline_pair_eligible"] for r in records),
        "n_icl_200_300_support_rows": sum(r["has_icl_200_300_support"] for r in records),
        "n_icl_satellite_association_rows": sum(r["has_icl_satellite_association_row"] for r in records),
        "n_with_satellite_50_75_profiles": sum(r["has_satellite_50_75_profiles"] for r in records),
        "inputs_sha256": {
            str(stack_summary_path): sha256(stack_summary_path),
            str(audit_path): sha256(audit_path),
            str(icl_summary_path): sha256(icl_summary_path),
            str(icl_table): sha256(icl_table),
            str(assoc_summary_path): sha256(assoc_summary_path),
            str(assoc_table): sha256(assoc_table),
            str(satellite_summary_path): sha256(satellite_summary_path),
            str(satellite_table): sha256(satellite_table),
        },
        "matches": records,
        "limitations": [
            "This checks whether control-field WCS candidate flags overlap saved support rows; it does not infer flux contamination.",
            "The association table has 94 selected cluster/filter rows and is exploratory; this audit does not refit or resample the association.",
            "No image pixels, HSC services, new data, or new measurements were used.",
        ],
    }


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stack", type=Path, default=DEFAULT_STACK)
    parser.add_argument("--icl", type=Path, default=DEFAULT_ICL)
    parser.add_argument("--association", type=Path, default=DEFAULT_ASSOCIATION)
    parser.add_argument("--satellite", type=Path, default=DEFAULT_SATELLITE)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args(argv)
    output = args.output.resolve()
    if output.exists():
        parser.error(f"refusing to overwrite {output}")
    result = run(args.stack, args.icl, args.association, args.satellite)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n",
                      encoding="utf-8")
    print(json.dumps({k: result[k] for k in ("status", "n_candidate_control_overlaps",
                                             "n_icl_200_300_support_rows",
                                             "n_with_satellite_50_75_profiles",
                                             "n_icl_satellite_association_rows", "inputs_sha256")},
                     indent=2))
    print(f"output={output}")
    print(f"output_sha256={sha256(output)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
