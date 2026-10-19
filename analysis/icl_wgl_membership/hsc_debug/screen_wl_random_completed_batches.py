#!/usr/bin/env python3
"""Geometry-screen only the verified HSC random-coordinate batches already returned.

This is an explicitly partial pipeline check while the approved coordinate
export is still running. It never promotes these early batches to a final
random-center sample and requests no additional catalog rows.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
from pathlib import Path
import platform
import statistics
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from fetch_ancillary_catalogues import catalogue_reader, sha
from prepare_control_candidates import CAMIRA_SHA
from screen_wl_random_control_coordinates import load_preflight_cosmology, screen_geometry
from submit_wl_random_control_coordinates import (
    DEFAULT_PLAN, load_count_inputs, read_coordinates,
)

DEFAULT_CATALOG = Path("/gpfs/kjhan/Astrafek/data/catalog/camira_s20a_wide_sm.tbl")


def load_verified_completed_batches(coordinate_dir: Path, batches: list[dict],
                                    lens_ids: list[str], counts: dict[str, int],
                                    batch_size: int, checkpoint: dict,
                                    plan_sha256: str) -> tuple[list[dict], list[dict]]:
    """Load complete, hash-matched batches and skip all not-yet-complete batches."""
    if checkpoint.get("plan_sha256") != plan_sha256:
        raise ValueError("coordinate checkpoint belongs to a different plan")
    if checkpoint.get("status") not in ("running", "stopped_on_error", "interrupted", "complete"):
        raise ValueError("coordinate checkpoint is not an active or terminal export")
    coordinates, batch_records = [], []
    jobs = checkpoint.get("jobs", {})
    for index, batch in enumerate(batches):
        name = f"random_coordinates_{index:03d}"
        item = jobs.get(name, {})
        if item.get("status") != "complete":
            continue
        path = coordinate_dir / f"{name}.csv"
        if not path.is_file() or sha(path) != item.get("sha256"):
            raise ValueError(f"completed coordinate batch is missing or changed: {name}")
        start = index * batch_size
        ids = lens_ids[start:start + int(batch["n_lenses"])]
        expected = {lens_id: counts[lens_id] for lens_id in ids}
        rows = read_coordinates(path, expected)
        if rows != int(batch["expected_rows"]) or rows != int(item.get("rows", -1)):
            raise ValueError(f"completed coordinate batch row count mismatch: {name}")
        batch_rows = []
        with path.open(encoding="utf-8-sig", newline="") as stream:
            batch_rows.extend(dict(row) for row in catalogue_reader(stream))
        if len(batch_rows) != rows:
            raise ValueError(f"parsed coordinate batch row count mismatch: {name}")
        coordinates.extend({"lens_id": str(row["lens_id"]),
                            "ra_deg": float(row["ra_deg"]),
                            "dec_deg": float(row["dec_deg"])}
                           for row in batch_rows)
        batch_records.append({"name": name, "job_id": str(item["job_id"]),
                              "lens_ids": ids,
                              "n_lenses": len(ids), "n_coordinates": rows,
                              "sha256": item["sha256"]})
    if not batch_records:
        raise ValueError("no verified completed coordinate batches are available")
    return coordinates, batch_records


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", type=Path, default=DEFAULT_PLAN)
    parser.add_argument("--catalog", type=Path, default=DEFAULT_CATALOG)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--chunk-size", type=int, default=10000)
    args = parser.parse_args()
    if (not platform.node().lower().startswith("syntax")
            and not os.environ.get("SLURM_JOB_ID")):
        parser.error("Run this geometry screen on syntax or in a Slurm job")
    if args.work.exists():
        parser.error(f"Refusing to overwrite existing output: {args.work}")
    if args.chunk_size < 1:
        parser.error("chunk size must be positive")

    plan_dir = args.plan.resolve()
    plan_path = plan_dir / "random_control_coordinate_plan.json"
    coordinate_plan = json.loads(plan_path.read_text(encoding="utf-8"))
    count_plan_dir = Path(coordinate_plan["count_plan"])
    manifest, lenses, lens_ids, counts, count_plan_sha = load_count_inputs(count_plan_dir)
    if count_plan_sha != coordinate_plan.get("count_plan_sha256"):
        parser.error("coordinate plan references a changed count-only plan")
    batches = coordinate_plan["batches"]
    coordinate_dir = plan_dir / "cas_coordinates"
    checkpoint_path = coordinate_dir / "coordinate_status.json"
    checkpoint = json.loads(checkpoint_path.read_text(encoding="utf-8"))
    coordinates, completed_batches = load_verified_completed_batches(
        coordinate_dir, batches, lens_ids, counts, int(manifest["batch_size"]),
        checkpoint, sha(plan_path))
    selected_ids = {lens_id for batch in completed_batches for lens_id in batch["lens_ids"]}
    selected_lenses = [{"lens_id": str(row["lens_id"]), "z_lens": float(row["z_lens"])}
                       for row in lenses if str(row["lens_id"]) in selected_ids]
    if len(selected_lenses) != sum(batch["n_lenses"] for batch in completed_batches):
        raise ValueError("completed batches do not map to the expected unique lenses")

    catalog = args.catalog.resolve()
    digest = sha(catalog)
    if digest != CAMIRA_SHA:
        raise ValueError("CAMIRA source checksum changed")
    camira = np.loadtxt(catalog, usecols=(0, 2, 3, 4), ndmin=2)
    if camira.shape[1] != 4 or len(np.unique(camira[:, 0])) != len(camira):
        raise ValueError("CAMIRA catalogue schema or IDs are invalid")
    if (np.any(camira[:, 1] < 0) or np.any(camira[:, 1] >= 360)
            or np.any(np.abs(camira[:, 2]) > 90)):
        raise ValueError("CAMIRA catalogue has invalid sky coordinates")
    count_plan_path = count_plan_dir / "random_control_count_plan.json"
    count_manifest = json.loads(count_plan_path.read_text(encoding="utf-8"))
    preflight = Path(count_manifest["source_preflight"])
    h0, omega_m = load_preflight_cosmology(preflight)
    audited = screen_geometry(coordinates, selected_lenses, camira,
                              h0, omega_m, args.chunk_size)

    per_lens = {str(row["lens_id"]): {"n_coordinates": 0,
                                      "n_geometry_candidates": 0}
                for row in selected_lenses}
    for row in audited:
        stats = per_lens.setdefault(row["lens_id"], {"n_coordinates": 0,
                                                       "n_geometry_candidates": 0})
        stats["n_coordinates"] += 1
        stats["n_geometry_candidates"] += int(row["geometry_candidate"])
    args.work.mkdir(parents=True, exist_ok=False)
    geometry_path = args.work / "partial_candidate_geometry_audit.tsv"
    fields = list(audited[0])
    with geometry_path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(audited)
    per_lens_path = args.work / "partial_per_lens_candidate_counts.tsv"
    with per_lens_path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, ["lens_id", "n_coordinates", "n_geometry_candidates"],
                                delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for lens_id in sorted(per_lens, key=lambda value: (int(value) if value.isdigit() else value)):
            writer.writerow({"lens_id": lens_id, **per_lens[lens_id]})
    coordinate_plan_sha = sha(plan_path)
    summary = {
        "schema": "astrafek-wl-random-coordinate-partial-geometry-audit-v1",
        "status": "complete_partial_geometry_pilot",
        "coordinate_export_terminal": False,
        "coordinate_plan_sha256": coordinate_plan_sha,
        "count_plan_sha256": count_plan_sha,
        "completed_batches": completed_batches,
        "n_completed_batches": len(completed_batches),
        "n_lenses": len(selected_lenses),
        "n_coordinates": len(audited),
        "n_camira": len(camira),
        "n_geometry_candidates": int(sum(bool(row["geometry_candidate"]) for row in audited)),
        "targets_with_geometry_candidates": sum(v["n_geometry_candidates"] > 0 for v in per_lens.values()),
        "targets_without_geometry_candidates": sum(v["n_geometry_candidates"] == 0 for v in per_lens.values()),
        "per_lens_geometry_candidate_counts": {
            "minimum": min(v["n_geometry_candidates"] for v in per_lens.values()),
            "median": statistics.median(v["n_geometry_candidates"] for v in per_lens.values()),
            "maximum": max(v["n_geometry_candidates"] for v in per_lens.values()),
        },
        "n_rows_overlapping_camira_apertures": sum(
            int(row["camira_3mpc_overlap_count"]) > 0 for row in audited),
        "n_rows_with_duplicate_coordinates": sum(
            int(row["coordinate_multiplicity"]) > 1 for row in audited),
        "n_rows_both_overlap_and_duplicate": sum(
            int(row["camira_3mpc_overlap_count"]) > 0
            and int(row["coordinate_multiplicity"]) > 1 for row in audited),
        "camira_catalog_sha256": digest,
        "full_aperture_footprint_verified": False,
        "hsc_i_depth_psf_sky_matched": False,
        "random_controls_selected": False,
        "science_ready": False,
        "science_measurement": False,
        "limitations": [
            "Only currently completed early export batches were screened; this is not representative of all 800 target lenses.",
            "The geometry screen only removes duplicate coordinates and candidates whose 3 proper-Mpc aperture overlaps any CAMIRA cluster aperture.",
            "Full-aperture image coverage, masks, depth, PSF, sky, extinction, control count, and spatial covariance remain untested.",
            "No random centers were selected and no shape, photo-z, image, or source rows were requested.",
        ],
        "outputs_sha256": {geometry_path.name: sha(geometry_path),
                           per_lens_path.name: sha(per_lens_path)},
        "partial_screen_script_sha256": sha(Path(__file__).resolve()),
        "geometry_screen_script_sha256": sha(Path(screen_geometry.__code__.co_filename).resolve()),
    }
    summary_path = args.work / "summary.json"
    summary_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n",
                            encoding="utf-8")
    print(json.dumps({key: summary[key] for key in (
        "status", "n_completed_batches", "n_lenses", "n_coordinates",
        "n_geometry_candidates", "science_ready", "science_measurement")}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
