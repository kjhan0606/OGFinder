#!/usr/bin/env python3
"""Select a reproducible, balanced random-center subset for a count-only gate.

This is an offline planner. It consumes the hashed CAMIRA geometry screen,
assigns each center the redshift of its matched CAMIRA lens, and writes a
manifest for the existing random shape-count planner. It does not contact CAS
or retrieve images, source rows, shapes, or redshifts.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

from fetch_ancillary_catalogues import sha
from screen_wl_random_control_coordinates import load_preflight_cosmology
from submit_wl_random_control_coordinates import validate_coordinate_plan


def _integer(value, field):
    try:
        parsed = int(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"invalid {field}") from exc
    if parsed < 0:
        raise ValueError(f"negative {field}")
    return parsed


def deterministic_select(candidates, redshifts, per_lens, seed):
    """Choose up to N geometry-clean coordinates per lens, without replacement."""
    if not isinstance(per_lens, int) or per_lens < 1:
        raise ValueError("per_lens must be a positive integer")
    if not str(seed):
        raise ValueError("seed must be non-empty")
    groups = {}
    seen_positions = set()
    for source in candidates:
        lens_id = str(source["lens_id"])
        if lens_id not in redshifts:
            raise ValueError(f"unknown matched lens {lens_id}")
        ra, dec = float(source["ra_deg"]), float(source["dec_deg"])
        if not (0 <= ra < 360 and -90 <= dec <= 90):
            raise ValueError("invalid candidate coordinates")
        position = (ra, dec)
        if position in seen_positions:
            raise ValueError("geometry-clean candidate pool contains duplicate coordinates")
        seen_positions.add(position)
        rank_material = f"{seed}\0{lens_id}\0{ra:.12f}\0{dec:.12f}".encode()
        rank = hashlib.sha256(rank_material).hexdigest()
        groups.setdefault(lens_id, []).append((rank, ra, dec))

    chosen, per_lens_counts = [], {}
    for lens_id in sorted(redshifts):
        ranked = sorted(groups.get(lens_id, []))
        take = ranked[:per_lens]
        per_lens_counts[lens_id] = {"available": len(ranked), "selected": len(take)}
        for rank, ra, dec in take:
            center_material = f"{lens_id}\0{ra:.12f}\0{dec:.12f}".encode()
            center_id = "RC-" + hashlib.sha256(center_material).hexdigest()[:20]
            chosen.append({
                "random_center_id": center_id,
                "matched_lens_id": lens_id,
                "ra_deg": ra,
                "dec_deg": dec,
                "z_lens": float(redshifts[lens_id]),
                "geometry_rank_sha256": rank,
            })
    return chosen, per_lens_counts


def prepare(geometry_dir: Path, output: Path, *, per_lens: int, seed: str):
    geometry_dir, output = Path(geometry_dir).resolve(), Path(output).resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite {output}")
    summary_path = geometry_dir / "summary.json"
    candidate_path = geometry_dir / "candidate_geometry_audit.tsv"
    per_lens_path = geometry_dir / "per_lens_candidate_counts.tsv"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    if (summary.get("schema") != "astrafek-wl-random-coordinate-geometry-audit-v1"
            or summary.get("status") != "complete_geometry_screen_only"
            or summary.get("science_measurement") is not False
            or summary.get("random_controls_selected") is not False):
        raise ValueError("geometry screen is incomplete or has unexpected science state")
    if (sha(candidate_path) != summary.get("candidate_geometry_audit_sha256")
            or sha(per_lens_path) != summary.get("per_lens_candidate_counts_sha256")):
        raise ValueError("geometry candidate table checksum mismatch")

    coordinate_plan = Path(summary["coordinate_plan"]).resolve()
    coord_manifest, _, lenses, lens_ids, expected_counts = validate_coordinate_plan(coordinate_plan)
    if (sha(coordinate_plan / "random_control_coordinate_plan.json") != summary.get("coordinate_plan_sha256")
            or sha(coordinate_plan / "cas_coordinates/coordinate_summary.json") != summary.get("coordinate_summary_sha256")):
        raise ValueError("geometry audit does not match the approved coordinate export")
    coord_status = json.loads((coordinate_plan / "cas_coordinates/coordinate_status.json").read_text())
    coord_summary = json.loads((coordinate_plan / "cas_coordinates/coordinate_summary.json").read_text())
    if coord_status.get("status") != "complete" or coord_summary.get("status") != "complete":
        raise ValueError("coordinate export is not complete")
    if (coord_status.get("credentials_saved") is not False
            or any(coord_status.get(key) is not False for key in (
                "source_rows_returned", "shape_rows_returned", "photoz_rows_returned",
                "images_returned", "science_measurement"))
            or coord_summary.get("output_columns") != ["lens_id", "ra_deg", "dec_deg"]
            or coord_summary.get("coordinate_plan_sha256") != sha(coordinate_plan / "random_control_coordinate_plan.json")
            or coord_summary.get("n_coordinate_rows") != sum(expected_counts.values())
            or coord_status.get("total_coordinate_rows") != coord_summary.get("n_coordinate_rows")):
        raise ValueError("coordinate export scope or totals are invalid")
    if summary.get("n_lenses") != 800 or summary.get("n_coordinates") != sum(expected_counts.values()):
        raise ValueError("geometry screen totals disagree with approved coordinate plan")
    if len(lens_ids) != 800 or len(lenses) != 800:
        raise ValueError("expected exactly 800 matched CAMIRA lenses")
    redshifts = {str(row["lens_id"]): float(row["z_lens"]) for row in lenses}

    raw_counts = {}
    with per_lens_path.open(encoding="utf-8", newline="") as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            lens_id = str(row["lens_id"])
            if lens_id in raw_counts:
                raise ValueError(f"duplicate per-lens geometry summary: {lens_id}")
            raw_counts[lens_id] = {
                "available": _integer(row["n_geometry_candidates"], "geometry candidate count"),
                "n_candidates": _integer(row["n_candidates"], "candidate count"),
            }
    if set(raw_counts) - set(redshifts):
        raise ValueError("per-lens geometry table contains an unknown lens")
    for lens_id in redshifts:
        raw = raw_counts.get(lens_id, {"available": 0, "n_candidates": 0})
        if raw["n_candidates"] != expected_counts[lens_id]:
            raise ValueError(f"per-lens coordinate count mismatch for {lens_id}")

    candidates, seen_geometry = [], {lens_id: 0 for lens_id in redshifts}
    with candidate_path.open(encoding="utf-8", newline="") as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            if row["geometry_candidate"] != "True":
                continue
            if row.get("science_ready") != "False":
                raise ValueError("geometry candidate must remain explicitly non-science-ready")
            lens_id = str(row["lens_id"])
            if lens_id not in redshifts:
                raise ValueError(f"geometry candidate references unknown lens {lens_id}")
            candidates.append({"lens_id": lens_id, "ra_deg": row["ra_deg"], "dec_deg": row["dec_deg"]})
            seen_geometry[lens_id] += 1
    for lens_id in redshifts:
        if seen_geometry[lens_id] != raw_counts.get(lens_id, {}).get("available", 0):
            raise ValueError(f"per-lens geometry count mismatch for {lens_id}")
    if len(candidates) != summary.get("n_geometry_candidates"):
        raise ValueError("geometry candidate total mismatch")

    centers, counts = deterministic_select(candidates, redshifts, per_lens, seed)
    if not centers:
        raise ValueError("geometry screen yielded no selectable centers")
    count_manifest = json.loads(
        (Path(coord_manifest["count_plan"]) / "random_control_count_plan.json").read_text(encoding="utf-8")
    )
    h0, omega_m = load_preflight_cosmology(Path(count_manifest["source_preflight"]))
    center_manifest = {
        "schema": "astrafek-wl-random-centers-v1",
        "status": "selected_for_shape_count_preflight",
        "h0_km_s_mpc": h0,
        "omega_m": omega_m,
        "science_measurement": False,
        "selection": {
            "method": "lowest SHA-256 ranks per matched lens; coordinates selected without replacement",
            "seed": str(seed), "requested_per_lens": per_lens,
            "n_matched_lenses": len(redshifts), "n_random_centers": len(centers),
            "per_lens_counts": counts,
        },
        "provenance": {
            "coordinate_plan": str(coordinate_plan),
            "coordinate_plan_sha256": summary["coordinate_plan_sha256"],
            "coordinate_summary_sha256": summary["coordinate_summary_sha256"],
            "geometry_summary": str(summary_path), "geometry_summary_sha256": sha(summary_path),
            "candidate_geometry_audit_sha256": sha(candidate_path),
            "per_lens_candidate_counts_sha256": sha(per_lens_path),
        },
        "full_aperture_footprint_verified": False,
        "depth_psf_sky_matched": False,
        "random_center_overlap_audited": False,
        "spatial_covariance_validated": False,
        "execution_requires_explicit_user_approval": True,
        "centers": centers,
    }
    output.mkdir(parents=True, exist_ok=False)
    manifest_path = output / "random_centers.json"
    manifest_path.write_text(json.dumps(center_manifest, ensure_ascii=False, indent=2, allow_nan=False) + "\n")
    result = {
        "schema": "astrafek-wl-random-center-selection-audit-v1",
        "status": "selected_for_shape_count_preflight",
        "n_geometry_candidates": len(candidates), "n_random_centers": len(centers),
        "n_matched_lenses_with_selected_centers": sum(v["selected"] > 0 for v in counts.values()),
        "requested_per_lens": per_lens, "seed": str(seed),
        "manifest": str(manifest_path), "manifest_sha256": sha(manifest_path),
        "science_measurement": False, "science_ready": False,
        "hsc_query_submitted": False, "hsc_data_returned": False,
    }
    (output / "selection_audit.json").write_text(json.dumps(result, indent=2) + "\n")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--geometry", type=Path, required=True)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--per-lens", type=int, required=True,
                        help="maximum selected centers per target lens; must be chosen explicitly")
    parser.add_argument("--seed", required=True, help="recorded deterministic selection seed")
    args = parser.parse_args()
    result = prepare(args.geometry, args.work, per_lens=args.per_lens, seed=args.seed)
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
