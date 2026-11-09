#!/usr/bin/env python3
"""Compare defensible S23B CAMIRA redshift/spectroscopy cuts offline.

The audit deduplicates S23B candidates against the full local S20A CAMIRA
catalog at 1 arcsec, excludes the current 800 lenses, applies a 6-arcmin
richness-priority separation, and counts newly occupied 5-degree cells. It
does not select a download sample, query HSC, or read image pixels.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree

from audit_s23b_camira_catalog import (
    EXPECTED_S23B_BLOB_SHA1,
    EXPECTED_S23B_BYTES,
    greedy_additions,
    read_catalog,
    sha256,
    sky_cell,
    unit_vectors,
)


DATA = Path("/gpfs/kjhan/Astrafek/data")
S23B = DATA / "catalogue_s23b_camira_2026-10-20_01/camira_s23b_wide_sm.tbl"
S20A = DATA / "catalog/camira_s20a_wide_sm.tbl"
MANIFEST = DATA / "sample_800_2026-10-13/sample_manifest.json"
BASELINE_AUDIT = DATA / "catalogue_s23b_camira_audit_2026-10-20_02"
BASELINE_SUMMARY_SHA256 = "f98e9bdb0b7f5cfaf217132cafa6c8a0cf9187a0ec448cc9ee1c50826c800e20"
BASELINE_CANDIDATE_SHA256 = "1b59ca9941849fa7970c43c6e0374db6587716fc4ae5d71b11f84a245ebe6b2c"
EXPECTED_S23B_SHA256 = "5af057028feab7976cd3ed389a7246cbe73add8bbb96891101f0f692b1de4d71"
EXPECTED_S20A_SHA256 = "d929985d7f6ee38c4435fc12099e9df1e5dbfd9ce9603cea319d5704129d22aa"
EXPECTED_MANIFEST_SHA256 = "95388814062e95f6e8429e2543a5432433515c9d9607680d2e873488d2758401"
CURRENT_MATCH_ARCSEC = 1.0
MIN_SEPARATION_DEG = 0.1
SCENARIOS = (
    ("baseline_spec", 0.30, 0.55, True),
    ("broader_spec", 0.25, 0.65, True),
    ("wide_spec", 0.20, 0.70, True),
    ("baseline_cluster_z_only", 0.30, 0.55, False),
    ("broader_cluster_z_only", 0.25, 0.65, False),
    ("wide_cluster_z_only", 0.20, 0.70, False),
)
DEFAULT_OUTPUT = DATA / "sample_800_2026-10-13/icl_s23b_camira_selection_scenarios_2026-11-05_02"


def center_rows(manifest_path: Path) -> list[dict]:
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    rows = [{
        "id": str(row["id"]),
        "name": str(row.get("name", row["id"])),
        "ra_deg": float(row["ra_deg"]),
        "dec_deg": float(row["dec_deg"]),
        "z": float(row["z_cl"]),
        "richness": float(row["richness"]),
        "z_bcg": float(row["z_spec_bcg"]),
    } for row in manifest.get("clusters", [])]
    if len(rows) != 800 or len({row["id"] for row in rows}) != 800:
        raise ValueError("sample manifest must contain exactly 800 unique lenses")
    return rows


def match_mask(candidates: list[dict], reference: list[dict], radius_arcsec: float) -> np.ndarray:
    radius_chord = 2 * math.sin(math.radians(radius_arcsec / 3600) / 2)
    tree = cKDTree(unit_vectors(reference))
    return np.asarray([
        bool(hits) for hits in tree.query_ball_point(unit_vectors(candidates), radius_chord)
    ], dtype=bool)


def eligible(row: dict, z_min: float, z_max: float, require_spec: bool) -> bool:
    if not z_min <= float(row["z"]) < z_max:
        return False
    if not require_spec:
        return True
    z_bcg = float(row["z_bcg"])
    return z_bcg > 0 and abs(z_bcg - float(row["z"])) / (1 + float(row["z"])) < 0.01


def distribution(rows: list[dict], key: str) -> dict:
    values = np.asarray([float(row[key]) for row in rows], dtype=float)
    if not values.size:
        return {"n": 0}
    return {
        "n": int(values.size), "min": float(values.min()),
        "median": float(np.median(values)), "p90": float(np.quantile(values, 0.9)),
        "max": float(values.max()),
    }


def run(output: Path) -> dict:
    if output.exists():
        raise FileExistsError(f"refusing to overwrite {output}")
    for path in (S23B, S20A, MANIFEST):
        if not path.is_file():
            raise FileNotFoundError(path)
    if (sha256(S23B) != EXPECTED_S23B_SHA256 or S23B.stat().st_size != EXPECTED_S23B_BYTES
            or sha256(S20A) != EXPECTED_S20A_SHA256
            or sha256(MANIFEST) != EXPECTED_MANIFEST_SHA256):
        raise ValueError("one or more pinned CAMIRA/sample input hashes changed")
    manifest_blob_sha1 = hashlib.sha1(
        f"blob {S23B.stat().st_size}\0".encode("ascii") + S23B.read_bytes()
    ).hexdigest()
    if manifest_blob_sha1 != EXPECTED_S23B_BLOB_SHA1:
        raise ValueError("S23B Git blob hash differs from the approved catalog")

    current = center_rows(MANIFEST)
    s23b = read_catalog(S23B)
    s20a = read_catalog(S20A)
    s20a_match = match_mask(s23b, s20a, CURRENT_MATCH_ARCSEC)
    current_match = match_mask(s23b, current, CURRENT_MATCH_ARCSEC)
    current_cells = {sky_cell(row) for row in current}
    scenarios = {}
    candidate_tables: dict[str, list[dict]] = {}

    for name, z_min, z_max, require_spec in SCENARIOS:
        before_cross_match = [row for row in s23b if eligible(row, z_min, z_max, require_spec)]
        before_ids = {row["id"] for row in before_cross_match}
        deduplicated = [row for index, row in enumerate(s23b)
                        if eligible(row, z_min, z_max, require_spec)
                        and not s20a_match[index]]
        additions = greedy_additions(current, deduplicated)
        additions_by_cell = {sky_cell(row) for row in additions}
        candidate_tables[name] = additions
        scenarios[name] = {
            "redshift_range": [z_min, z_max],
            "requires_consistent_bcg_spectroscopic_redshift": require_spec,
            "eligible_s23b_rows_before_cross_release_dedup": len(before_cross_match),
            "eligible_rows_within_1arcsec_of_s20a": int(sum(
                eligible(row, z_min, z_max, require_spec) and bool(s20a_match[i])
                for i, row in enumerate(s23b))),
            "eligible_rows_matching_current_800_within_1arcsec": int(sum(
                eligible(row, z_min, z_max, require_spec) and bool(current_match[i])
                for i, row in enumerate(s23b))),
            "eligible_after_s20a_dedup_before_current_sample_separation": len(deduplicated),
            "greedy_additions_after_6arcmin_separation": len(additions),
            "new_5deg_cells_outside_current_800": len(additions_by_cell - current_cells),
            "total_5deg_cells_after_additions": len(current_cells | additions_by_cell),
            "z_distribution": distribution(additions, "z"),
            "richness_distribution": distribution(additions, "richness"),
            "candidates_are_not_image_supported": True,
        }

    baseline_summary = BASELINE_AUDIT / "summary.json"
    baseline_table = BASELINE_AUDIT / "greedy_candidate_additions.tsv"
    if (sha256(baseline_summary) != BASELINE_SUMMARY_SHA256
            or sha256(baseline_table) != BASELINE_CANDIDATE_SHA256):
        raise ValueError("baseline S23B candidate audit differs from its pinned review")
    baseline_ids = set()
    with baseline_table.open(encoding="utf-8", newline="") as stream:
        baseline_ids = {row["id"] for row in csv.DictReader(stream, delimiter="\t")}
    recomputed_baseline_ids = {row["id"] for row in candidate_tables["baseline_spec"]}
    if baseline_ids != recomputed_baseline_ids:
        raise ValueError("recomputed baseline candidate set differs from the approved audit")

    output.mkdir(parents=True, exist_ok=False)
    tables = {}
    for name, rows in candidate_tables.items():
        table_path = output / f"{name}_candidates.tsv"
        fields = ["id", "name", "ra_deg", "dec_deg", "z", "richness", "z_bcg"]
        with table_path.open("x", encoding="utf-8", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
        tables[name] = {"path": str(table_path), "rows": len(rows), "sha256": sha256(table_path)}
    result = {
        "schema": "astrafek-s23b-camira-selection-scenario-audit-v1",
        "status": "complete_offline_catalog_scenario_audit",
        "science_measurement": False,
        "hsc_query_submitted": False,
        "coordinates_exported_to_external_service": False,
        "images_or_shapes_downloaded": False,
        "image_pixels_read": False,
        "catalog_rows": {"s23b": len(s23b), "s20a": len(s20a), "current_sample": len(current)},
        "cross_release_match_radius_arcsec": CURRENT_MATCH_ARCSEC,
        "minimum_selected_center_separation_arcmin": MIN_SEPARATION_DEG * 60,
        "current_parent_5deg_cells": len(current_cells),
        "scenarios": scenarios,
        "candidate_tables": tables,
        "baseline_reproduction": {
            "baseline_summary_path": str(baseline_summary),
            "baseline_summary_sha256": sha256(baseline_summary),
            "baseline_candidate_table_sha256": sha256(baseline_table),
            "exact_candidate_ids_reproduced": True,
            "n_candidate_ids": len(baseline_ids),
        },
        "inputs": {str(path): sha256(path) for path in (S23B, S20A, MANIFEST)},
        "limitations": [
            "A wider redshift range changes the cluster population and requires a prespecified science question and lens-redshift calibration.",
            "Cluster-redshift-only scenarios omit the BCG spectroscopic consistency cut and require additional photo-z/lensing calibration.",
            "Catalog rows and new 5-degree cells do not establish image, mask, depth, PSF, source completeness, or matched-control support.",
            "The 6-arcmin greedy selector is a deterministic footprint upper-bound procedure, not a globally optimal independent set.",
            "5-degree cells upper-bound connected covariance blocks under the current design and are not proven independent.",
        ],
    }
    summary_path = output / "summary.json"
    summary_path.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    result["summary_path"] = str(summary_path)
    result["summary_sha256"] = sha256(summary_path)
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    result = run(args.output.resolve())
    print(json.dumps({
        "status": result["status"],
        "science_measurement": False,
        "scenarios": {name: {
            "additions": row["greedy_additions_after_6arcmin_separation"],
            "new_cells": row["new_5deg_cells_outside_current_800"],
            "total_cells": row["total_5deg_cells_after_additions"],
        } for name, row in result["scenarios"].items()},
        "baseline_reproduction": result["baseline_reproduction"],
        "output": result["summary_path"],
        "summary_sha256": result["summary_sha256"],
    }, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
