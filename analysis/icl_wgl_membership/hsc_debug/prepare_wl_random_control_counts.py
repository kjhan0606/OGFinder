#!/usr/bin/env python3
"""Prepare count-only availability queries for HSC PDR3 random positions.

This program writes local SQL and a provenance manifest; it never contacts HSC
CAS and never retrieves random coordinates. Each query returns one integer
candidate count for each of the 800 approved lens centers.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from submit_wl_shape_counts import load_plan


DEFAULT_PREFLIGHT = Path(
    "/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/"
    "wl_shape_preflight_2026-10-18_01"
)
DEFAULT_OUTPUT = Path(
    "/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/"
    "wl_random_point_count_plan_2026-10-19_01"
)
BATCH_SIZE = 25
OUTER_RADIUS_ARCSEC = 1800.0
ADJUST_DENSITY_MAX = 0.01
TARGET_CANDIDATES_PER_LENS = 200


def _sql_string(value):
    text = str(value)
    if not text or any(char not in "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_-" for char in text):
        raise ValueError("lens_id contains unsafe characters")
    return "'" + text + "'"


def _sql_number(value):
    number = float(value)
    if not math.isfinite(number):
        raise ValueError("lens coordinates and radii must be finite")
    return format(number, ".12g")


def random_count_query(lens):
    """Count sparse eligible randoms in the non-overlapping local annulus."""
    lens_id = _sql_string(lens["lens_id"])
    ra = _sql_number(lens["ra_deg"])
    dec = _sql_number(lens["dec_deg"])
    aperture_radius = float(lens["radius_arcsec"])
    if not math.isfinite(aperture_radius) or aperture_radius <= 0:
        raise ValueError("lens aperture radius must be finite and positive")
    # Two equal 3 Mpc apertures do not overlap if their centers are at least
    # twice the lens aperture radius apart. The 30 arcmin outer limit matches
    # the existing control-candidate planning radius.
    inner = 2.0 * aperture_radius
    if inner >= OUTER_RADIUS_ARCSEC:
        raise ValueError("non-overlapping control annulus is empty for a lens")
    inner_text = _sql_number(inner)
    outer_text = _sql_number(OUTER_RADIUS_ARCSEC)
    density_text = _sql_number(ADJUST_DENSITY_MAX)
    return (
        "SELECT " + lens_id + " AS lens_id, COUNT(*) AS n_random_candidates\n"
        "FROM pdr3_wide.random AS r\n"
        f"WHERE coneSearch(r.coord, {ra}, {dec}, {outer_text})\n"
        f"  AND NOT coneSearch(r.coord, {ra}, {dec}, {inner_text})\n"
        f"  AND r.adjust_density < {density_text}\n"
        "  AND r.isprimary IS TRUE\n"
        "  AND r.patch_inner IS TRUE\n"
        "  AND r.tract_inner IS TRUE\n"
        "  AND r.i_pixelflags_edge IS FALSE\n"
        "  AND r.i_pixelflags_bright_objectcenter IS FALSE"
    )


def prepare(preflight, output):
    preflight = Path(preflight).resolve()
    output = Path(output).resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite existing plan directory: {output}")
    source_plan = load_plan(preflight)
    lenses = source_plan["lenses"]
    if len(lenses) != 800:
        raise ValueError("expected exactly 800 approved lens centers")
    if len({str(item["lens_id"]) for item in lenses}) != len(lenses):
        raise ValueError("lens centers must have unique IDs")
    chunks = [lenses[start:start + BATCH_SIZE]
              for start in range(0, len(lenses), BATCH_SIZE)]
    queries = ["\nUNION ALL\n".join(random_count_query(lens) for lens in chunk) + "\n"
               for chunk in chunks]
    output.mkdir(parents=True, exist_ok=False)
    batches = []
    for index, (chunk, sql) in enumerate(zip(chunks, queries)):
        path = output / f"random_counts_{index:03d}.sql"
        path.write_text(sql, encoding="utf-8")
        batches.append({
            "path": str(path),
            "n_lenses": len(chunk),
            "first_lens_id": str(chunk[0]["lens_id"]),
            "last_lens_id": str(chunk[-1]["lens_id"]),
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        })

    manifest = {
        "schema": "hsc-pdr3-random-control-count-plan-v1",
        "status": "prepared_not_submitted",
        "release": "pdr3",
        "source_preflight": str(preflight),
        "source_preflight_sha256": hashlib.sha256(
            (preflight / "preflight_manifest.json").read_bytes()
        ).hexdigest(),
        "random_table": "pdr3_wide.random",
        "n_lenses": len(lenses),
        "n_query_batches": len(batches),
        "batch_size": BATCH_SIZE,
        "candidate_query": {
            "outer_radius_arcsec": OUTER_RADIUS_ARCSEC,
            "inner_radius": "twice each lens's angular 3 Mpc aperture",
            "random_density_thinning": "adjust_density < 0.01 (nominally 1 point/arcmin^2 from the documented 100/arcmin^2)",
            "quality_filters": [
                "isprimary IS TRUE", "patch_inner IS TRUE", "tract_inner IS TRUE",
                "i_pixelflags_edge IS FALSE", "i_pixelflags_bright_objectcenter IS FALSE",
            ],
            "candidate_target_per_lens": TARGET_CANDIDATES_PER_LENS,
        },
        "output_columns": ["lens_id", "n_random_candidates"],
        "random_coordinates_returned": False,
        "shape_rows_returned": False,
        "source_rows_returned": False,
        "network_used": False,
        "query_submitted": False,
        "execution_requires_explicit_user_approval": True,
        "later_gates": [
            "available count is only a candidate-pool screen, not random-control selection",
            "after coordinate export, exclude all CAMIRA cluster apertures and audit footprint coverage",
            "match local HSC-I depth, PSF, sky, extinction, and mask completeness",
            "audit random-center uniqueness, overlapping apertures, and spatial-block covariance",
            "choose final random count from measured covariance, not the preliminary target of 200",
        ],
        "batches": batches,
    }
    (output / "random_control_count_plan.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8",
    )
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--preflight", type=Path, default=DEFAULT_PREFLIGHT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    manifest = prepare(args.preflight, args.output)
    print(json.dumps({
        "status": manifest["status"],
        "network_used": False,
        "query_submitted": False,
        "lenses": manifest["n_lenses"],
        "batches": manifest["n_query_batches"],
        "output_columns": manifest["output_columns"],
        "output": str(Path(args.output).resolve()),
    }, ensure_ascii=False))


if __name__ == "__main__":
    main()
