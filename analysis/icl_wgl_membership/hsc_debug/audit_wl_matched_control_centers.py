#!/usr/bin/env python3
"""Test whether saved ICL matched-control centers can support WL random subtraction.

This offline geometry audit uses no HSC query and reads no source-shape rows.
It checks the previously downloaded image-control centers against the full
CAMIRA catalogue and against one another at the paired target redshift.
"""
from __future__ import annotations

import argparse
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

from fetch_ancillary_catalogues import sha
from prepare_control_candidates import CAMIRA_SHA
from screen_wl_random_control_coordinates import screen_geometry


def control_overlap_counts(controls, h0: float, omega_m: float):
    """Count intersecting 3-proper-Mpc apertures among control centers."""
    cosmo = FlatLambdaCDM(H0=h0, Om0=omega_m, Tcmb0=0)
    vectors = np.column_stack((
        np.cos(np.deg2rad([r["dec_deg"] for r in controls]))
        * np.cos(np.deg2rad([r["ra_deg"] for r in controls])),
        np.cos(np.deg2rad([r["dec_deg"] for r in controls]))
        * np.sin(np.deg2rad([r["ra_deg"] for r in controls])),
        np.sin(np.deg2rad([r["dec_deg"] for r in controls])),
    ))
    radii_deg = np.rad2deg(np.array([
        3.0 / cosmo.angular_diameter_distance(r["z_lens"]).value
        for r in controls
    ]))
    tree = cKDTree(vectors)
    counts = np.zeros(len(controls), dtype=int)
    max_chord = 2 * math.sin(math.radians(2 * float(radii_deg.max())) / 2)
    for i, point in enumerate(vectors):
        neighbors = tree.query_ball_point(point, max_chord)
        for j in neighbors:
            if j <= i:
                continue
            chord = float(np.linalg.norm(point - vectors[j]))
            separation_deg = math.degrees(2 * math.asin(min(1.0, chord / 2)))
            if separation_deg <= radii_deg[i] + radii_deg[j]:
                counts[i] += 1
                counts[j] += 1
    return counts, radii_deg


def audit(controls_manifest: Path, target_manifest: Path, camira_catalog: Path):
    controls_doc = json.loads(controls_manifest.read_text(encoding="utf-8"))
    targets_doc = json.loads(target_manifest.read_text(encoding="utf-8"))
    if controls_doc.get("purpose") != "matched_control_image_validation":
        raise ValueError("control manifest is not the saved ICL matched-control sample")
    controls = controls_doc.get("clusters", [])
    targets = targets_doc.get("clusters", [])
    if not controls or not targets:
        raise ValueError("target or control manifest has no centers")
    target_by_id = {str(r["id"]): r for r in targets}
    if len(target_by_id) != len(targets):
        raise ValueError("target manifest contains duplicate IDs")
    ids = [str(r.get("target_cluster_id", "")) for r in controls]
    object_ids = [str(r.get("random_object_id", "")) for r in controls]
    coordinates = [(float(r["ra_deg"]), float(r["dec_deg"])) for r in controls]
    if len(set(ids)) != len(ids) or len(set(object_ids)) != len(object_ids):
        raise ValueError("control target IDs or HSC random object IDs are duplicated")
    if len(set(coordinates)) != len(coordinates):
        raise ValueError("control sky coordinates are duplicated")

    normalized = []
    for row in controls:
        target_id = str(row["target_cluster_id"])
        if target_id not in target_by_id:
            raise ValueError(f"control references unknown target {target_id}")
        target = target_by_id[target_id]
        z = float(row["z"])
        if not math.isfinite(z) or z <= 0 or abs(z - float(target["z"])) > 1e-8:
            raise ValueError(f"paired target redshift mismatch for {target_id}")
        ra, dec = float(row["ra_deg"]), float(row["dec_deg"])
        if not math.isfinite(ra) or not 0 <= ra < 360 or not math.isfinite(dec) or not -90 <= dec <= 90:
            raise ValueError(f"invalid control coordinate for {target_id}")
        normalized.append({
            "lens_id": target_id,
            "ra_deg": ra,
            "dec_deg": dec,
            "z_lens": z,
            "random_object_id": str(row["random_object_id"]),
        })

    camira_sha = sha(camira_catalog)
    if camira_sha != CAMIRA_SHA:
        raise ValueError("CAMIRA catalogue checksum differs from the audited S20A input")
    camira = np.loadtxt(camira_catalog, usecols=(0, 2, 3, 4), ndmin=2)
    if camira.shape[1] != 4 or not np.all(np.isfinite(camira)) or np.any(camira[:, 3] <= 0):
        raise ValueError("CAMIRA catalogue has invalid coordinate/redshift fields")

    overlap_rows = screen_geometry(normalized, [
        {"lens_id": r["lens_id"], "z_lens": r["z_lens"]} for r in normalized
    ], camira, float(controls_doc["h0"]), float(controls_doc["omega_m"]))
    camira_clean_indices = [
        i for i, row in enumerate(overlap_rows)
        if row["camira_3mpc_overlap_count"] == 0
    ]
    clean_controls = [normalized[i] for i in camira_clean_indices]
    clean_overlap_counts, _ = control_overlap_counts(
        clean_controls, float(controls_doc["h0"]), float(controls_doc["omega_m"])
    )
    clean_overlap_by_index = dict(zip(camira_clean_indices, clean_overlap_counts))
    # Keep radii for the full table for reference, but use the CAMIRA-clean pool
    # when testing random-aperture self-overlap.
    cosmo = FlatLambdaCDM(H0=float(controls_doc["h0"]), Om0=float(controls_doc["omega_m"]), Tcmb0=0)
    radii_deg = np.rad2deg(np.array([
        3.0 / cosmo.angular_diameter_distance(r["z_lens"]).value
        for r in normalized
    ]))
    rows = []
    for i, (control, overlap) in enumerate(zip(normalized, overlap_rows)):
        camira_clean = overlap["camira_3mpc_overlap_count"] == 0
        n_clean_neighbors = int(clean_overlap_by_index.get(i, 0))
        rows.append({
            **control,
            "target_3mpc_radius_arcsec": 3600.0 * float(radii_deg[i]),
            "n_overlapping_full_camira_3mpc_apertures": int(overlap["camira_3mpc_overlap_count"]),
            "n_overlapping_camira_clean_control_3mpc_apertures": n_clean_neighbors,
            "geometry_candidate": bool(
                camira_clean and n_clean_neighbors == 0
            ),
        })
    summary = {
        "schema": "astrafek-wl-icl-control-center-geometry-audit-v1",
        "status": "complete_geometry_only",
        "controls_manifest": str(controls_manifest.resolve()),
        "controls_manifest_sha256": sha(controls_manifest),
        "target_manifest": str(target_manifest.resolve()),
        "target_manifest_sha256": sha(target_manifest),
        "camira_catalog": str(camira_catalog.resolve()),
        "camira_catalog_sha256": camira_sha,
        "n_target_lenses": len(targets),
        "n_control_centers": len(normalized),
        "n_controls_with_full_camira_3mpc_overlap": int(sum(
            r["n_overlapping_full_camira_3mpc_apertures"] > 0 for r in rows
        )),
        "n_camira_clean_controls": len(camira_clean_indices),
        "n_camira_clean_controls_with_control_3mpc_overlap": int(sum(
            r["n_overlapping_camira_clean_control_3mpc_apertures"] > 0 for r in rows
        )),
        "n_geometry_candidates": int(sum(r["geometry_candidate"] for r in rows)),
        "exclusion_rule": "control aperture and any supplied S20A CAMIRA 3 proper-Mpc aperture may not intersect",
        "control_independence_rule": "control 3 proper-Mpc apertures may not intersect each other",
        "script_sha256": sha(Path(__file__).resolve()),
        "source_shape_rows_read": False,
        "hsc_query_submitted": False,
        "science_measurement": False,
        "random_controls_selected": False,
        "full_aperture_depth_psf_sky_match_verified": False,
        "limitations": [
            "The CAMIRA veto covers the supplied S20A catalogue, not every real structure.",
            "The original control sample was selected for ICL image matching, not as a calibrated WL random sample.",
            "Passing aperture geometry alone does not validate shape-catalog coverage, source selection, mask area, or spatial covariance.",
        ],
        "controls": rows,
    }
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--controls", type=Path, required=True, help="control_download_01/download_manifest.json")
    parser.add_argument("--targets", type=Path, required=True, help="download_01/download_manifest.json")
    parser.add_argument("--catalog", type=Path, default=Path("/gpfs/kjhan/Astrafek/data/catalog/camira_s20a_wide_sm.tbl"))
    parser.add_argument("--output", type=Path, required=True, help="new output directory")
    args = parser.parse_args()
    if not platform.node().lower().startswith("lageunha"):
        parser.error("Run this CPU audit on lageunha")
    if args.output.exists():
        parser.error(f"Refusing to overwrite output: {args.output}")
    result = audit(args.controls.resolve(), args.targets.resolve(), args.catalog.resolve())
    args.output.mkdir(parents=True, exist_ok=False)
    tsv = args.output / "control_geometry.tsv"
    with tsv.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=[
            "lens_id", "random_object_id", "ra_deg", "dec_deg", "z_lens",
            "target_3mpc_radius_arcsec", "n_overlapping_full_camira_3mpc_apertures",
            "n_overlapping_camira_clean_control_3mpc_apertures", "geometry_candidate",
        ], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows([{key: row[key] for key in writer.fieldnames} for row in result["controls"]])
    result["control_table"] = str(tsv.resolve())
    result["control_table_sha256"] = sha(tsv)
    summary_path = args.output / "summary.json"
    summary_path.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: result[key] for key in (
        "status", "n_target_lenses", "n_control_centers",
        "n_controls_with_full_camira_3mpc_overlap",
        "n_camira_clean_controls", "n_camira_clean_controls_with_control_3mpc_overlap",
        "n_geometry_candidates",
        "source_shape_rows_read", "hsc_query_submitted", "science_measurement",
    )}, indent=2))
    print(f"summary_sha256={sha(summary_path)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
