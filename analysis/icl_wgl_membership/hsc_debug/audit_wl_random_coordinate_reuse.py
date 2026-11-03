#!/usr/bin/env python3
"""Compare strict global and lens-scoped reuse rules for random controls.

This is an offline sensitivity audit of the completed geometry table. It does
not query HSC, alter the geometry screen, select random centers, or measure
weak lensing. The relaxed candidate pool allows a sky coordinate to be paired
with more than one matched target lens, while still requiring no CAMIRA
3-proper-Mpc aperture overlap for each pair and no duplicate within a lens.
"""
from __future__ import annotations

import argparse
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path

from fetch_ancillary_catalogues import sha


def summarize_rows(rows):
    """Summarize strict and lens-scoped reuse policies from audit rows."""
    strict_lenses = Counter()
    relaxed_lenses = Counter()
    relaxed_positions = defaultdict(set)
    seen_lens_positions = set()
    total = 0
    strict_total = 0
    relaxed_total = 0
    duplicate_within_lens = 0

    for row in rows:
        total += 1
        lens_id = str(row["lens_id"])
        ra, dec = str(row["ra_deg"]), str(row["dec_deg"])
        key = (ra, dec)
        lens_key = (lens_id, ra, dec)
        overlaps = int(row["camira_3mpc_overlap_count"])
        strict = row["geometry_candidate"] in (True, "True", "true", "1")
        if strict:
            strict_total += 1
            strict_lenses[lens_id] += 1
        if overlaps != 0:
            continue
        if lens_key in seen_lens_positions:
            duplicate_within_lens += 1
            continue
        seen_lens_positions.add(lens_key)
        relaxed_total += 1
        relaxed_lenses[lens_id] += 1
        relaxed_positions[key].add(lens_id)

    shared = {position: lens_ids for position, lens_ids in relaxed_positions.items()
              if len(lens_ids) > 1}
    n_lenses = len(set(strict_lenses) | set(relaxed_lenses))
    return {
        "n_coordinate_rows": total,
        "strict_global_unique": {
            "n_candidate_rows": strict_total,
            "n_lenses_with_candidates": sum(v > 0 for v in strict_lenses.values()),
            "n_lenses_at_least_5": sum(v >= 5 for v in strict_lenses.values()),
            "n_lenses_at_least_10": sum(v >= 10 for v in strict_lenses.values()),
            "per_lens_counts": dict(strict_lenses),
        },
        "lens_scoped_reuse": {
            "n_candidate_rows": relaxed_total,
            "n_lenses_with_candidates": sum(v > 0 for v in relaxed_lenses.values()),
            "n_lenses_at_least_5": sum(v >= 5 for v in relaxed_lenses.values()),
            "n_lenses_at_least_10": sum(v >= 10 for v in relaxed_lenses.values()),
            "n_unique_sky_positions": len(relaxed_positions),
            "n_positions_shared_across_lenses": len(shared),
            "n_rows_at_shared_positions": sum(len(v) for v in shared.values()),
            "n_excess_reuse_rows": sum(len(v) - 1 for v in shared.values()),
            "max_lenses_per_position": max(map(len, relaxed_positions.values()), default=0),
            "n_duplicate_rows_within_lens_excluded": duplicate_within_lens,
            "per_lens_counts": dict(relaxed_lenses),
        },
        "n_distinct_lenses": n_lenses,
    }


def audit(geometry_dir: Path, output: Path):
    geometry_dir, output = Path(geometry_dir).resolve(), Path(output).resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite {output}")
    summary_path = geometry_dir / "summary.json"
    table_path = geometry_dir / "candidate_geometry_audit.tsv"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    if (summary.get("schema") != "astrafek-wl-random-coordinate-geometry-audit-v1"
            or summary.get("status") != "complete_geometry_screen_only"
            or summary.get("science_measurement") is not False
            or summary.get("random_controls_selected") is not False):
        raise ValueError("geometry screen is incomplete or has unexpected science state")
    if sha(table_path) != summary.get("candidate_geometry_audit_sha256"):
        raise ValueError("candidate geometry table checksum mismatch")

    with table_path.open(encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    result = summarize_rows(rows)
    if result["n_coordinate_rows"] != summary.get("n_coordinates"):
        raise ValueError("geometry row count differs from summary")
    strict = result["strict_global_unique"]
    if strict["n_candidate_rows"] != summary.get("n_geometry_candidates"):
        raise ValueError("strict candidate count differs from summary")
    result.update({
        "schema": "astrafek-wl-random-coordinate-reuse-audit-v1",
        "status": "offline_reuse_sensitivity_only",
        "geometry_summary_sha256": sha(summary_path),
        "geometry_table_sha256": sha(table_path),
        "coordinate_plan_sha256": summary["coordinate_plan_sha256"],
        "policies": {
            "strict_global_unique": "retain the geometry_candidate flag from the original screen",
            "lens_scoped_reuse": "require zero CAMIRA overlap for this lens-coordinate pair and allow reuse across different lens IDs; duplicate pairs within a lens are excluded",
        },
        "random_controls_selected": False,
        "hsc_query_submitted": False,
        "science_ready": False,
        "science_measurement": False,
    })
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
                      encoding="utf-8")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--geometry", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = audit(args.geometry, args.output)
    concise = {key: result[key] for key in (
        "status", "n_coordinate_rows", "strict_global_unique", "lens_scoped_reuse",
        "hsc_query_submitted", "science_measurement")}
    for policy in ("strict_global_unique", "lens_scoped_reuse"):
        concise[policy] = {key: value for key, value in concise[policy].items()
                           if key != "per_lens_counts"}
    print(json.dumps(concise, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
