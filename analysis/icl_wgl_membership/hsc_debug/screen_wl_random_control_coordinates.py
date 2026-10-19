#!/usr/bin/env python3
"""Offline CAMIRA-overlap screen for an approved HSC random-coordinate export.

This is a geometry-only candidate audit. It does not select final random lenses,
validate full-aperture image coverage, match depth/PSF, or measure lensing.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import csv
import json
import math
from pathlib import Path
import platform
import sys

import numpy as np
from astropy.cosmology import FlatLambdaCDM
from scipy.spatial import cKDTree

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from fetch_ancillary_catalogues import catalogue_reader, sha
from prepare_control_candidates import CAMIRA_SHA
from submit_wl_random_control_coordinates import (
    DEFAULT_PLAN, read_coordinates, validate_coordinate_plan,
)

DEFAULT_CATALOG = Path("/gpfs/kjhan/Astrafek/data/catalog/camira_s20a_wide_sm.tbl")


def unit_vectors(ra_deg, dec_deg):
    ra = np.deg2rad(np.asarray(ra_deg, dtype=float))
    dec = np.deg2rad(np.asarray(dec_deg, dtype=float))
    return np.column_stack((np.cos(dec) * np.cos(ra),
                            np.cos(dec) * np.sin(ra), np.sin(dec)))


def angular_separation_deg(a, b):
    """Great-circle separation in degrees, stable across RA=0."""
    delta = np.linalg.norm(np.asarray(a) - np.asarray(b), axis=-1)
    return np.rad2deg(2.0 * np.arcsin(np.clip(delta / 2.0, 0.0, 1.0)))


def screen_geometry(coordinates, lenses, camira, h0, omega_m, chunk_size=10000):
    """Annotate candidate rows for full-CAMIRA 3-Mpc aperture overlap.

    `camira` is an N x 4 array of (cluster_id, ra_deg, dec_deg, z_cluster).
    Each control aperture uses the redshift of its associated target lens. An
    overlap is recorded when the two angular 3-Mpc radii intersect or touch.
    Duplicate coordinates are counted by exact numeric RA/Dec equality; the
    CAS random table returns the same stored coordinate for a repeated point.
    """
    if chunk_size < 1:
        raise ValueError("chunk size must be positive")
    if len(camira) == 0 or not np.all(np.isfinite(camira)) or np.any(camira[:, 3] <= 0):
        raise ValueError("CAMIRA catalogue must have finite coordinates and positive redshifts")
    lens_by_id = {str(row["lens_id"]): row for row in lenses}
    if len(lens_by_id) != len(lenses) or not lens_by_id:
        raise ValueError("lens IDs must be non-empty and unique")
    cosmo = FlatLambdaCDM(H0=float(h0), Om0=float(omega_m), Tcmb0=0)
    camira_vectors = unit_vectors(camira[:, 1], camira[:, 2])
    camira_tree = cKDTree(camira_vectors)
    camira_radius_deg = np.rad2deg(3.0 / cosmo.angular_diameter_distance(camira[:, 3]).value)
    lens_radius_deg = {}
    for lens_id, row in lens_by_id.items():
        z = float(row["z_lens"])
        if not math.isfinite(z) or z <= 0:
            raise ValueError(f"invalid target redshift for lens {lens_id}")
        lens_radius_deg[lens_id] = math.degrees(3.0 / float(cosmo.angular_diameter_distance(z).value))

    coordinate_keys = Counter(
        (float(row["ra_deg"]), float(row["dec_deg"]))
        for row in coordinates
    )
    output = []
    max_radius_deg = max(lens_radius_deg.values()) + float(camira_radius_deg.max())
    max_chord = 2.0 * math.sin(math.radians(max_radius_deg) / 2.0)
    for start in range(0, len(coordinates), chunk_size):
        block = coordinates[start:start + chunk_size]
        points = unit_vectors([r["ra_deg"] for r in block], [r["dec_deg"] for r in block])
        nearby_lists = camira_tree.query_ball_point(points, max_chord)
        for row, point, nearby in zip(block, points, nearby_lists):
            lens_id = str(row["lens_id"])
            if lens_id not in lens_by_id:
                raise ValueError(f"coordinate row references unknown lens {lens_id}")
            limit_deg = lens_radius_deg[lens_id]
            if nearby:
                ids = np.asarray(nearby, dtype=int)
                sep = angular_separation_deg(point, camira_vectors[ids])
                overlaps = int(np.count_nonzero(sep <= limit_deg + camira_radius_deg[ids]))
            else:
                overlaps = 0
            key = (float(row["ra_deg"]), float(row["dec_deg"]))
            multiplicity = coordinate_keys[key]
            output.append({
                "lens_id": lens_id,
                "ra_deg": float(row["ra_deg"]),
                "dec_deg": float(row["dec_deg"]),
                "target_3mpc_radius_arcsec": limit_deg * 3600.0,
                "camira_3mpc_overlap_count": overlaps,
                "coordinate_multiplicity": multiplicity,
                "geometry_candidate": overlaps == 0 and multiplicity == 1,
                "science_ready": False,
            })
    return output


def load_coordinates(plan_dir: Path):
    manifest, batches, lenses, lens_ids, counts = validate_coordinate_plan(plan_dir)
    result_dir = plan_dir / "cas_coordinates"
    status_path = result_dir / "coordinate_status.json"
    summary_path = result_dir / "coordinate_summary.json"
    if not status_path.is_file() or not summary_path.is_file():
        raise ValueError("coordinate export has no terminal checkpoint and summary")
    status = json.loads(status_path.read_text(encoding="utf-8"))
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    plan_sha = sha(plan_dir / "random_control_coordinate_plan.json")
    if status.get("status") != "complete" or status.get("plan_sha256") != plan_sha:
        raise ValueError("coordinate export is not complete for this plan")
    if summary.get("status") != "complete" or summary.get("coordinate_plan_sha256") != plan_sha:
        raise ValueError("coordinate summary is incomplete or belongs to another plan")
    scope_flags = ("source_rows_returned", "shape_rows_returned", "photoz_rows_returned",
                   "images_returned", "science_measurement")
    if status.get("credentials_saved") is not False or any(status.get(key) is not False for key in scope_flags):
        raise ValueError("coordinate checkpoint reports data or state outside the approved coordinate scope")
    if summary.get("output_columns") != ["lens_id", "ra_deg", "dec_deg"] or any(
            summary.get(key) is not False for key in scope_flags):
        raise ValueError("coordinate summary reports fields outside the approved coordinate scope")
    if summary.get("n_lenses") != 800 or summary.get("n_coordinate_rows") != sum(counts.values()):
        raise ValueError("coordinate summary does not match approved lens counts")
    if status.get("total_coordinate_rows") != summary.get("n_coordinate_rows"):
        raise ValueError("coordinate checkpoint and summary totals disagree")
    rows = []
    for index, batch in enumerate(batches):
        name = Path(batch["path"]).stem
        item = status.get("jobs", {}).get(name, {})
        path = result_dir / f"{name}.csv"
        if item.get("status") != "complete" or not path.is_file() or sha(path) != item.get("sha256"):
            raise ValueError(f"coordinate result missing or changed: {name}")
        start = index * int(manifest["batch_size"])
        expected = {lens_id: counts[lens_id] for lens_id in lens_ids[start:start + int(manifest["batch_size"])]}
        read_coordinates(path, expected)
        with path.open(encoding="utf-8-sig", newline="") as stream:
            rows.extend(dict(row) for row in catalogue_reader(stream))
    if len(rows) != summary["n_coordinate_rows"]:
        raise ValueError("loaded coordinate total disagrees with summary")
    normalized = []
    for row in rows:
        ra, dec = float(row["ra_deg"]), float(row["dec_deg"])
        if not math.isfinite(ra) or not 0 <= ra < 360 or not math.isfinite(dec) or not -90 <= dec <= 90:
            raise ValueError("invalid coordinate in saved result")
        normalized.append({"lens_id": str(row["lens_id"]), "ra_deg": ra, "dec_deg": dec})
    lenses = [{"lens_id": str(row["lens_id"]), "z_lens": float(row["z_lens"])} for row in lenses]
    return manifest, summary, normalized, lenses


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", type=Path, default=DEFAULT_PLAN)
    parser.add_argument("--catalog", type=Path, default=DEFAULT_CATALOG)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--chunk-size", type=int, default=10000)
    args = parser.parse_args()
    if not platform.node().lower().startswith("lageunha"):
        parser.error("Run the 780k-row geometry screen on lageunha CPU")
    if args.work.exists():
        parser.error(f"Refusing to overwrite existing output: {args.work}")
    if args.chunk_size < 1:
        parser.error("chunk size must be positive")

    plan_dir = args.plan.resolve()
    manifest, coordinate_summary, coordinates, lenses = load_coordinates(plan_dir)
    catalog = args.catalog.resolve()
    digest = sha(catalog)
    if digest != CAMIRA_SHA:
        raise ValueError("CAMIRA source checksum changed")
    values = np.loadtxt(catalog, usecols=(0, 2, 3, 4), ndmin=2)
    if values.shape[1] != 4 or len(np.unique(values[:, 0])) != len(values):
        raise ValueError("CAMIRA catalogue schema or IDs are invalid")
    if np.any(values[:, 1] < 0) or np.any(values[:, 1] >= 360) or np.any(np.abs(values[:, 2]) > 90):
        raise ValueError("CAMIRA catalogue has invalid sky coordinates")
    preflight = Path(json.loads((Path(manifest["count_plan"]) / "random_control_count_plan.json").read_text())["source_preflight"])
    preflight_doc = json.loads((preflight / "preflight_manifest.json").read_text())
    audited = screen_geometry(coordinates, lenses, values,
                              preflight_doc["h0"], preflight_doc["omega_m"], args.chunk_size)

    per_lens = defaultdict(lambda: {"n_candidates": 0, "n_geometry_candidates": 0})
    for row in audited:
        stats = per_lens[row["lens_id"]]
        stats["n_candidates"] += 1
        stats["n_geometry_candidates"] += int(row["geometry_candidate"])
    args.work.mkdir(parents=True, exist_ok=False)
    audit_path = args.work / "candidate_geometry_audit.tsv"
    fields = list(audited[0]) if audited else ["lens_id", "ra_deg", "dec_deg", "target_3mpc_radius_arcsec",
                                               "camira_3mpc_overlap_count", "coordinate_multiplicity",
                                               "geometry_candidate", "science_ready"]
    with audit_path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(audited)
    per_lens_path = args.work / "per_lens_candidate_counts.tsv"
    with per_lens_path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, ["lens_id", "n_candidates", "n_geometry_candidates"], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for lens_id in sorted(per_lens):
            writer.writerow({"lens_id": lens_id, **per_lens[lens_id]})
    summary = {
        "schema": "astrafek-wl-random-coordinate-geometry-audit-v1",
        "status": "complete_geometry_screen_only",
        "coordinate_plan": str(plan_dir),
        "coordinate_plan_sha256": sha(plan_dir / "random_control_coordinate_plan.json"),
        "coordinate_summary_sha256": sha(plan_dir / "cas_coordinates/coordinate_summary.json"),
        "camira_catalog": str(catalog), "camira_catalog_sha256": digest,
        "n_lenses": len(lenses), "n_coordinates": len(audited), "n_camira": len(values),
        "n_geometry_candidates": sum(row["geometry_candidate"] for row in audited),
        "targets_with_geometry_candidates": sum(v["n_geometry_candidates"] > 0 for v in per_lens.values()),
        "overlap_rule": "control 3 proper-Mpc angular aperture plus each CAMIRA 3 proper-Mpc angular aperture, using their respective redshifts",
        "duplicate_rule": "exact numeric RA/Dec equality across all lens-associated candidate rows; only unique rows remain geometry candidates",
        "full_aperture_footprint_verified": False,
        "hsc_i_depth_psf_sky_matched": False,
        "random_controls_selected": False,
        "science_ready": False,
        "science_measurement": False,
        "next_gates": [
            "verify full-aperture image coverage and i-band masks",
            "match HSC-I depth, PSF, and sky to each target using validated image metadata",
            "choose random-center count using measured spatial-block covariance",
            "request explicit approval before exporting HSC shape rows around selected random centers",
        ],
    }
    summary["screen_script_sha256"] = sha(Path(__file__))
    summary["candidate_geometry_audit_sha256"] = sha(audit_path)
    summary["per_lens_candidate_counts_sha256"] = sha(per_lens_path)
    (args.work / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({k: summary[k] for k in ("status", "n_lenses", "n_coordinates", "n_geometry_candidates", "targets_with_geometry_candidates", "science_ready")}, indent=2))


if __name__ == "__main__":
    main()
