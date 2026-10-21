#!/usr/bin/env python3
"""Estimate a consistent, cross-release CAMIRA expansion upper bound offline.

The audit combines already downloaded S20A and S23B catalogs under the current
0.30 <= z_cl < 0.55 plus consistent BCG spectroscopy cut. It excludes current
targets, then greedily enforces six-arcminute separation across releases. It
does not query CAS, download data, or claim image support or scientific power.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import sys

import numpy as np
from scipy.spatial import cKDTree

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from audit_camira_sample_expansion import angular_sep_deg, greedy_additions, read_catalog, sky_cell

DATA = Path("/gpfs/kjhan/Astrafek/data")
DEFAULT_S20A = DATA / "catalog/camira_s20a_wide_sm.tbl"
DEFAULT_S23B = DATA / "catalogue_s23b_camira_2026-10-20_01/camira_s23b_wide_sm.tbl"
DEFAULT_MANIFEST = DATA / "sample_800_2026-10-13/sample_manifest.json"
DEFAULT_POWER = DATA / "sample_800_2026-10-13/icl_sample_size_plan_2026-10-21_01/sample_size_plan.json"
DEFAULT_OUTPUT = DATA / "sample_800_2026-10-13/icl_cross_release_expansion_2026-10-21_01"
EXPECTED_HASHES = {
    "s20a": "d929985d7f6ee38c4435fc12099e9df1e5dbfd9ce9603cea319d5704129d22aa",
    "s23b": "5af057028feab7976cd3ed389a7246cbe73add8bbb96891101f0f692b1de4d71",
    "manifest": "95388814062e95f6e8429e2543a5432433515c9d9607680d2e873488d2758401",
    "power": "df3d6d06adecaa310ab94df94bbebd1441262087520c16385927f285846f0b88",
}
Z_MIN, Z_MAX, SPEC_FRACTION_MAX = 0.30, 0.55, 0.01
MIN_SEPARATION_DEG = 0.1
CURRENT_MATCH_ARCSEC = 1.0


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def unit_vectors(rows):
    ra = np.deg2rad([float(row["ra_deg"]) for row in rows])
    dec = np.deg2rad([float(row["dec_deg"]) for row in rows])
    return np.column_stack((np.cos(dec) * np.cos(ra), np.cos(dec) * np.sin(ra), np.sin(dec)))


def angular_match_count(candidates, reference, radius_arcsec):
    if not candidates or not reference:
        return 0
    tree = cKDTree(unit_vectors(reference))
    radius = 2.0 * math.sin(math.radians(radius_arcsec / 3600.0) / 2.0)
    return sum(bool(hits) for hits in tree.query_ball_point(unit_vectors(candidates), radius))


def eligible_rows(rows, release):
    out = []
    for row in rows:
        z, z_bcg = float(row["z"]), float(row["z_bcg"])
        if not (Z_MIN <= z < Z_MAX and z_bcg > 0):
            continue
        if abs(z_bcg - z) / (1.0 + z) >= SPEC_FRACTION_MAX:
            continue
        out.append({
            **row,
            "release": release,
            "id": f"{release}:{row['id']}",
            "z": z,
            "z_bcg": z_bcg,
            "richness": float(row["richness"]),
        })
    return out


def evaluate(current, s20a_rows, s23b_rows, power_plan):
    current_ids = {str(row["id"]) for row in current}
    if len(current) != 800 or len(current_ids) != 800:
        raise ValueError("Expected 800 unique current CAMIRA targets")
    candidates_by_release = {
        "S20A": eligible_rows(s20a_rows, "S20A"),
        "S23B": eligible_rows(s23b_rows, "S23B"),
    }
    all_eligible = candidates_by_release["S20A"] + candidates_by_release["S23B"]
    current_matches = angular_match_count(all_eligible, current, CURRENT_MATCH_ARCSEC)
    current_vectors = unit_vectors(current)
    current_tree = cKDTree(current_vectors)
    match_chord = 2.0 * math.sin(math.radians(CURRENT_MATCH_ARCSEC / 3600.0) / 2.0)
    candidate_vectors = unit_vectors(all_eligible)
    valid = [row for row, vector in zip(all_eligible, candidate_vectors)
             if not current_tree.query_ball_point(vector, match_chord)]

    cross_pairs = 0
    left, right = candidates_by_release["S20A"], candidates_by_release["S23B"]
    if left and right:
        right_tree = cKDTree(unit_vectors(right))
        cross_chord = 2.0 * math.sin(math.radians(CURRENT_MATCH_ARCSEC / 3600.0) / 2.0)
        cross_pairs = sum(len(hits) for hits in right_tree.query_ball_point(unit_vectors(left), cross_chord))

    # The archived greedy selector uses richness then ID; release-prefixed IDs
    # make the cross-release tie-break deterministic (S20A before S23B).
    additions = greedy_additions(current, valid, MIN_SEPARATION_DEG)
    selected_by_release = {
        release: sum(row["release"] == release for row in additions)
        for release in ("S20A", "S23B")
    }
    current_cells = {sky_cell(row) for row in current}
    added_cells = {sky_cell(row) for row in additions} - current_cells

    parent_n = int(power_plan["input"]["current_parent_lenses"])
    current_support = {k: int(v) for k, v in power_plan["input"]["current_joint_support_rows"].items()}
    yields = {k: float(v) for k, v in power_plan["input"]["observed_support_yield"].items()}
    projected = {band: current_support[band] + len(additions) * yields[band] for band in yields}
    scenarios = []
    for effect in (0.20, 0.293):
        for deff in (1.0, 1.5):
            match = next((row for row in power_plan["scenarios"]
                          if row["planning_partial_r"] == effect and row["target_power"] == 0.80
                          and row["multiplicity_scenario"] == "bonferroni_168"
                          and row["design_effect"] == deff), None)
            if match is None:
                raise ValueError(f"Missing sample-size scenario for r={effect}, design effect={deff}")
            required = int(match["inflated_usable_rows"])
            scenarios.append({
                "assumed_partial_r": effect,
                "design_effect": deff,
                "required_i2_rows": required,
                "projected_i2_rows_if_current_yield_transfers": projected["i2"],
                "fraction_of_required": projected["i2"] / required,
                "meets_scenario": projected["i2"] >= required,
            })
    return {
        "schema": "astrafek-icl-cross-release-expansion-audit-v1",
        "status": "offline_catalog_feasibility_only",
        "science_measurement": False,
        "cas_query_submitted": False,
        "images_or_shapes_downloaded": False,
        "selection": {
            "z_cl": [Z_MIN, Z_MAX],
            "require_consistent_bcg_spectroscopic_redshift": True,
            "relative_redshift_difference_limit": SPEC_FRACTION_MAX,
            "exclude_current_targets_within_arcsec": CURRENT_MATCH_ARCSEC,
            "greedy_minimum_separation_arcmin": MIN_SEPARATION_DEG * 60,
            "greedy_order": "richness descending, then release-prefixed catalog ID",
        },
        "catalog_counts": {
            "current_targets": len(current),
            "s20a_rows": len(s20a_rows),
            "s23b_rows": len(s23b_rows),
            "s20a_eligible": len(candidates_by_release["S20A"]),
            "s23b_eligible": len(candidates_by_release["S23B"]),
            "eligible_rows_within_1arcsec_of_current_targets": current_matches,
            "s20a_s23b_eligible_cross_pairs_within_1arcsec": cross_pairs,
            "union_rows_after_current_target_exclusion": len(valid),
            "greedy_additions": len(additions),
            "selected_from_each_release": selected_by_release,
            "new_5deg_ra_dec_cells": len(added_cells),
            "total_5deg_ra_dec_cells": len(current_cells | added_cells),
            "catalog_only_parent_upper_bound": parent_n + len(additions),
        },
        "support_projection": {
            "observed_yields_assumed": yields,
            "expected_joint_support_rows_if_yield_transfers": {
                band: len(additions) * yields[band] for band in yields
            },
            "projected_total_support_rows_if_yield_transfers": projected,
            "sample_size_scenarios": scenarios,
        },
        "limitations": [
            "catalogue rows are not evidence of HSC-I image, mask, depth, PSF, or matched-control support",
            "six-arcminute separation and 5-degree cells do not establish independent covariance blocks",
            "S20A/S23B catalogue entries can differ; cross-release near-duplicates are naturally removed by greedy separation but are not a physical identity proof",
            "current i/i2 joint-support yield may not transfer to the selected union",
            "power values are exploratory assumptions, not measured effect sizes or detection claims",
            "PSF-wing, sky-transfer, satellite-completeness, and coadd covariance gates remain unresolved",
        ],
    }, additions


def run(s20a_path: Path, s23b_path: Path, manifest_path: Path, power_path: Path):
    for key, path in (("s20a", s20a_path), ("s23b", s23b_path),
                      ("manifest", manifest_path), ("power", power_path)):
        if not path.is_file():
            raise FileNotFoundError(path)
        if sha256(path) != EXPECTED_HASHES[key]:
            raise ValueError(f"Reviewed {key} input hash mismatch")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    power_plan = json.loads(power_path.read_text(encoding="utf-8"))
    if len(manifest.get("clusters", [])) != 800 or manifest.get("catalog_sha256") != EXPECTED_HASHES["s20a"]:
        raise ValueError("Current sample manifest does not match the pinned S20A catalog")
    current = [{
        "id": str(row["id"]), "ra_deg": float(row["ra_deg"]), "dec_deg": float(row["dec_deg"]),
        "z": float(row["z_cl"]), "richness": float(row["richness"]), "z_bcg": float(row["z_spec_bcg"]),
    } for row in manifest["clusters"]]
    result, additions = evaluate(current, read_catalog(s20a_path), read_catalog(s23b_path), power_plan)
    result["inputs"] = {str(path.resolve()): sha256(path)
                        for path in (s20a_path, s23b_path, manifest_path, power_path)}
    return result, additions


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--s20a", type=Path, default=DEFAULT_S20A)
    parser.add_argument("--s23b", type=Path, default=DEFAULT_S23B)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--power-plan", type=Path, default=DEFAULT_POWER)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args(argv)
    if args.output_dir.exists():
        parser.error(f"refusing to overwrite {args.output_dir}")
    result, additions = run(args.s20a, args.s23b, args.manifest, args.power_plan)
    args.output_dir.mkdir(parents=True, exist_ok=False)
    candidate_path = args.output_dir / "greedy_cross_release_candidates.tsv"
    with candidate_path.open("x", encoding="utf-8", newline="") as stream:
        fields = ["id", "name", "release", "ra_deg", "dec_deg", "z", "z_bcg", "richness"]
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for row in additions:
            writer.writerow({key: row.get(key) for key in fields})
    result["candidate_table"] = {"path": str(candidate_path.resolve()), "sha256": sha256(candidate_path)}
    summary_path = args.output_dir / "summary.json"
    summary_path.write_text(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False) + "\n",
                            encoding="utf-8")
    print(json.dumps({
        "status": result["status"],
        "greedy_additions": result["catalog_counts"]["greedy_additions"],
        "selected_from_each_release": result["catalog_counts"]["selected_from_each_release"],
        "new_5deg_cells": result["catalog_counts"]["new_5deg_ra_dec_cells"],
        "projected_joint_support_rows": result["support_projection"]["projected_total_support_rows_if_yield_transfers"],
        "output_dir": str(args.output_dir.resolve()),
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
