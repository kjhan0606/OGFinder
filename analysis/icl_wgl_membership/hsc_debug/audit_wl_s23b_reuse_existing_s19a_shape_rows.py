#!/usr/bin/env python3
"""Audit whether the existing S19A shape export can support S23B apertures.

This is a read-only local analysis. It verifies the already downloaded
800-lens S19A source-row export, deduplicates source IDs, and measures source
and geometric coverage for the 427 S23B outer-profile targets. It makes no
network requests and exports no new source rows.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from astropy.cosmology import FlatLambdaCDM
from scipy.spatial import cKDTree
from scipy.stats import qmc


DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
PREFLIGHT_DIR = DATA / "wl_shape_preflight_2026-10-18_01"
PREFLIGHT = PREFLIGHT_DIR / "preflight_manifest.json"
ROW_PLAN_DIR = DATA / "wl_shape_rows_plan_2026-10-18_04"
ROW_PLAN = ROW_PLAN_DIR / "row_export_plan.json"
ROW_EXPORT_DIR = DATA / "wl_shape_rows_export_2026-10-18_03"
ROW_STATUS = ROW_EXPORT_DIR / "row_status.json"
ROW_SUMMARY = ROW_EXPORT_DIR / "row_summary.json"
S23B_MANIFEST = DATA / "s23b_icl_target_image_plan_2026-11-07_01/target_manifest.json"
S23B_RELATION = DATA / "s23b_camira_member_outer_icl_relation_2026-11-07_01/cluster_outer_relation_inputs.tsv"
DEFAULT_OUTPUT = DATA / "wl_s23b_existing_export_overlap_2026-11-07_01"

EXPECTED_ROW_PLAN_SHA256 = "3fc2f11105f956c5f4c940225c7b268bb70b22bcf17aaf0c62f3306a6a527557"
EXPECTED_PREFLIGHT_SHA256 = "a33932d2300368be14cf8286dbbb494c8ecf92df6424288f54b2552904be9dca"
EXPECTED_ROW_SUMMARY_SHA256 = "919e7d0bdc3bcbcf81a5336c1ba020a70e1ba23d8e42558bdf19703a9e2dc898"
EXPECTED_ROW_STATUS_SHA256 = "033e85eb98fdfb05226245fb7c7df9fd2b82fb7f14d20430c996532570026670"
EXPECTED_S23B_MANIFEST_SHA256 = "e487582ba40cd3e58da67704242db4d7219e381fb08f4b33fed47ce4ff268b44"
EXPECTED_S23B_RELATION_SHA256 = "b44c375ed003a460e289589472b2fddbe8281a079b336fb408c9509028b28674"

OUTPUT_COLUMNS = [
    "lens_id", "lens_ra_deg", "lens_dec_deg", "lens_z", "object_id",
    "source_ra_deg", "source_dec_deg", "e1", "e2", "sigma_e", "rms_e",
    "weight", "m", "c1", "c2", "mag_a10", "resolution", "psf_shape11",
    "psf_shape22", "psf_shape12", "b_mode_mask_raw", "source_zbin",
]
SOURCE_FIELDS = OUTPUT_COLUMNS[4:]
NUMERIC_SOURCE_FIELDS = [
    field for field in SOURCE_FIELDS
    if field not in ("object_id", "b_mode_mask_raw", "source_zbin")
]
OUTPUT_PER_TARGET_FIELDS = [
    "target_id", "filter", "ra_deg", "dec_deg", "z_cl",
    "aperture_radius_arcsec", "n_intersecting_old_cones",
    "nearest_old_cone_center_arcsec", "single_old_cone_contains_full_aperture",
    "old_cone_union_area_fraction_256", "old_cone_union_area_fraction_1024",
    "old_cone_union_area_fraction_4096", "n_existing_unique_sources_total",
    "n_shape_quality", "n_source_selection", "n_source_selection_zbin0",
    "n_source_selection_zbin1", "n_source_selection_zbin2",
    "n_source_selection_zbin3", "n_source_selection_zbin4",
    "n_source_selection_zbin1_4", "n_source_selection_zbin1_4_bmode_raw0",
    "n_source_selection_zbin1_4_bmode_raw1",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def unit_vector(ra_deg: float, dec_deg: float) -> np.ndarray:
    ra, dec = math.radians(ra_deg), math.radians(dec_deg)
    return np.array([math.cos(dec) * math.cos(ra),
                     math.cos(dec) * math.sin(ra),
                     math.sin(dec)], dtype=float)


def unit_vectors(ra_deg: np.ndarray, dec_deg: np.ndarray) -> np.ndarray:
    ra = np.deg2rad(ra_deg)
    dec = np.deg2rad(dec_deg)
    cos_dec = np.cos(dec)
    return np.column_stack((cos_dec * np.cos(ra),
                            cos_dec * np.sin(ra),
                            np.sin(dec)))


def separation_arcsec_vec(ra1, dec1, ra2, dec2) -> np.ndarray:
    ra1 = np.deg2rad(np.asarray(ra1, dtype=float))
    dec1 = np.deg2rad(np.asarray(dec1, dtype=float))
    ra2 = np.deg2rad(np.asarray(ra2, dtype=float))
    dec2 = np.deg2rad(np.asarray(dec2, dtype=float))
    dra = (ra2 - ra1 + np.pi) % (2.0 * np.pi) - np.pi
    ddec = dec2 - dec1
    hav = np.sin(ddec / 2.0) ** 2 + np.cos(dec1) * np.cos(dec2) * np.sin(dra / 2.0) ** 2
    hav = np.clip(hav, 0.0, 1.0)
    return 2.0 * np.arcsin(np.sqrt(hav)) * 206264.80624709636


def disk_sample(center: np.ndarray, radius_rad: float, uv: np.ndarray) -> np.ndarray:
    ra = math.atan2(center[1], center[0])
    dec = math.asin(float(center[2]))
    east = np.array([-math.sin(ra), math.cos(ra), 0.0])
    north = np.array([-math.sin(dec) * math.cos(ra),
                      -math.sin(dec) * math.sin(ra), math.cos(dec)])
    theta = np.arccos(np.clip(1.0 - uv[:, 0] * (1.0 - math.cos(radius_rad)), -1.0, 1.0))
    phi = 2.0 * np.pi * uv[:, 1]
    tangent = np.cos(phi)[:, None] * east + np.sin(phi)[:, None] * north
    points = (np.cos(theta)[:, None] * center[None, :]
              + np.sin(theta)[:, None] * tangent)
    return points / np.linalg.norm(points, axis=1)[:, None]


def require_sha(path: Path, expected: str, label: str) -> str:
    if not path.is_file():
        raise FileNotFoundError(path)
    actual = sha256(path)
    if actual != expected:
        raise ValueError(f"{label} SHA-256 changed: {actual}")
    return actual


def verify_export_and_load_sources():
    preflight_sha = require_sha(PREFLIGHT, EXPECTED_PREFLIGHT_SHA256, "S19A preflight")
    row_plan_sha = require_sha(ROW_PLAN, EXPECTED_ROW_PLAN_SHA256, "S19A row plan")
    row_summary_sha = require_sha(ROW_SUMMARY, EXPECTED_ROW_SUMMARY_SHA256, "S19A row summary")
    row_status_sha = require_sha(ROW_STATUS, EXPECTED_ROW_STATUS_SHA256, "S19A row status")
    preflight = read_json(PREFLIGHT)
    plan = read_json(ROW_PLAN)
    state = read_json(ROW_STATUS)
    summary = read_json(ROW_SUMMARY)
    if (preflight.get("status") != "prepared_not_submitted"
            or preflight.get("n_lenses") != 800
            or preflight.get("radius_mpc_proper") != 3.0
            or plan.get("status") != "prepared_not_submitted"
            or plan.get("n_lenses") != 800
            or plan.get("n_query_batches") != 32
            or plan.get("photoz_pdf_rows_included") is not False
            or state.get("status") != "complete"
            or state.get("plan_sha256") != row_plan_sha
            or state.get("source_rows_downloaded") is not True
            or state.get("science_measurement") is not False
            or summary.get("status") != "complete"
            or summary.get("n_lenses") != 800
            or summary.get("n_pair_rows") != 3069544
            or summary.get("columns") != OUTPUT_COLUMNS
            or summary.get("photoz_pdf_rows_included") is not False
            or summary.get("science_measurement") is not False):
        raise ValueError("S19A export plan/status/summary does not match the audited release")
    old_lenses = {str(row["lens_id"]): row for row in preflight["lenses"]}
    if len(old_lenses) != 800 or len(plan["batches"]) != 32 or len(state["jobs"]) != 32:
        raise ValueError("S19A lens or batch count is incomplete")
    for field, expected in (("lens_ra_deg", "ra_deg"),
                            ("lens_dec_deg", "dec_deg"),
                            ("lens_z", "z_lens")):
        if any(not math.isfinite(float(lens[expected])) for lens in old_lenses.values()):
            raise ValueError(f"invalid old lens geometry in {expected}")

    sources: dict[str, tuple] = {}
    batch_audit = []
    total_rows = 0
    max_cone_excess_arcsec = -math.inf
    source_duplicates = 0
    source_conflicts = 0
    for index, batch in enumerate(plan["batches"]):
        job_key = f"shape_rows_{index:03d}"
        job = state["jobs"].get(job_key)
        expected_name = f"shape_rows_{index:03d}.csv"
        path = ROW_EXPORT_DIR / expected_name
        if (job is None or job.get("status") != "complete"
                or Path(batch["path"]).stem != f"shape_rows_{index:03d}"
                or not path.is_file()):
            raise ValueError(f"missing or incomplete source-row batch {job_key}")
        file_sha = sha256(path)
        if (file_sha != job.get("sha256")
                or path.stat().st_size != int(job.get("bytes", -1))):
            raise ValueError(f"source-row batch receipt failed for {job_key}")
        expected_batch_lenses = {
            str(row["lens_id"])
            for row in preflight["lenses"][index * 25:(index + 1) * 25]
        }
        if len(expected_batch_lenses) != int(batch.get("n_lenses", -1)):
            raise ValueError(f"preflight batch lens membership differs for {job_key}")
        lens_ids_seen = set()
        pair_ids = set()
        source_ra_rows, source_dec_rows = [], []
        lens_ra_rows, lens_dec_rows, lens_radius_rows = [], [], []
        rows_in_batch = 0
        with path.open(encoding="utf-8-sig", newline="") as stream:
            reader = csv.reader(stream)
            raw_header = next(reader, None)
            if raw_header is None:
                raise ValueError(f"empty source-row file: {path}")
            header = [value.lstrip("#").strip() for value in raw_header]
            if header != OUTPUT_COLUMNS:
                raise ValueError(f"unexpected CSV header in {path}: {header}")
            for values in reader:
                if len(values) != len(OUTPUT_COLUMNS):
                    raise ValueError(f"malformed CSV row in {path}")
                row = dict(zip(OUTPUT_COLUMNS, values))
                lens_id = row["lens_id"]
                object_id = row["object_id"]
                if (lens_id not in expected_batch_lenses or not object_id
                        or (lens_id, object_id) in pair_ids):
                    raise ValueError(f"unexpected or repeated lens-source pair in {path}")
                pair_ids.add((lens_id, object_id))
                lens_ids_seen.add(lens_id)
                lens = old_lenses[lens_id]
                for field, source_key in (("lens_ra_deg", "ra_deg"),
                                          ("lens_dec_deg", "dec_deg"),
                                          ("lens_z", "z_lens")):
                    if abs(float(row[field]) - float(lens[source_key])) > 1e-9:
                        raise ValueError(f"changed lens metadata for {lens_id} in {path}")
                numeric = {field: float(row[field]) for field in NUMERIC_SOURCE_FIELDS}
                if not all(math.isfinite(value) for value in numeric.values()):
                    raise ValueError(f"non-finite shape value in {path}")
                source_ra, source_dec = numeric["source_ra_deg"], numeric["source_dec_deg"]
                if not (0.0 <= source_ra < 360.0 and -90.0 <= source_dec <= 90.0):
                    raise ValueError(f"invalid source coordinates in {path}")
                mask = int(row["b_mode_mask_raw"])
                zbin = int(row["source_zbin"])
                if str(mask) != row["b_mode_mask_raw"] or mask not in (0, 1):
                    raise ValueError(f"invalid raw B-mode mask in {path}")
                if str(zbin) != row["source_zbin"] or zbin < 0:
                    raise ValueError(f"invalid source redshift bin in {path}")
                lens_radius = float(lens["radius_arcsec"])
                lens_ra_rows.append(float(lens["ra_deg"]))
                lens_dec_rows.append(float(lens["dec_deg"]))
                lens_radius_rows.append(lens_radius)
                source_ra_rows.append(source_ra)
                source_dec_rows.append(source_dec)

                digest_payload = "\x1f".join(row[field] for field in SOURCE_FIELDS).encode("utf-8")
                source_digest = hashlib.blake2b(digest_payload, digest_size=16).digest()
                shape_quality = (
                    0.0 <= numeric["sigma_e"] <= 0.4
                    and math.hypot(numeric["e1"], numeric["e2"]) < 2.0
                )
                source_selection = (
                    shape_quality and numeric["mag_a10"] <= 25.5
                    and numeric["resolution"] >= 0.3
                )
                source_value = (
                    source_digest, source_ra, source_dec, shape_quality,
                    source_selection, zbin, mask,
                )
                previous = sources.get(object_id)
                if previous is None:
                    sources[object_id] = source_value
                else:
                    source_duplicates += 1
                    if previous != source_value:
                        source_conflicts += 1
                        raise ValueError(f"source-level fields conflict for object_id {object_id}")
                rows_in_batch += 1
        if rows_in_batch != int(job.get("rows", -1)):
            raise ValueError(f"source-row count differs from receipt for {job_key}")
        if rows_in_batch != int(batch.get("expected_pair_rows_from_count_only", -2)):
            raise ValueError(f"source-row count differs from count-only plan for {job_key}")
        sep = separation_arcsec_vec(
            np.asarray(lens_ra_rows), np.asarray(lens_dec_rows),
            np.asarray(source_ra_rows), np.asarray(source_dec_rows),
        )
        cone_excess = sep - np.asarray(lens_radius_rows)
        batch_max_excess = float(np.max(cone_excess)) if len(cone_excess) else None
        if batch_max_excess is not None and batch_max_excess > 0.05:
            raise ValueError(f"source rows fall outside their old query cones in {job_key}")
        if batch_max_excess is not None:
            max_cone_excess_arcsec = max(max_cone_excess_arcsec, batch_max_excess)
        total_rows += rows_in_batch
        batch_audit.append({
            "batch": index, "file": path.name, "rows": rows_in_batch,
            "bytes": path.stat().st_size, "sha256": file_sha,
            "n_lens_ids_with_rows": len(lens_ids_seen),
            "max_cone_edge_excess_arcsec": batch_max_excess,
        })
    if (total_rows != 3069544 or total_rows != int(state.get("total_pair_rows", -1))
            or total_rows != int(summary.get("count_summary_sum", -2))
            or source_conflicts != 0):
        raise ValueError("S19A aggregate row or source-consistency audit failed")
    return {
        "preflight": preflight, "old_lenses": old_lenses, "sources": sources,
        "batch_audit": batch_audit, "total_rows": total_rows,
        "unique_sources": len(sources), "repeated_source_pair_rows": source_duplicates,
        "max_cone_edge_excess_arcsec": max_cone_excess_arcsec,
        "input_hashes": {
            str(PREFLIGHT): preflight_sha, str(ROW_PLAN): row_plan_sha,
            str(ROW_STATUS): row_status_sha, str(ROW_SUMMARY): row_summary_sha,
        },
    }


def load_s23b_targets():
    manifest_sha = require_sha(S23B_MANIFEST, EXPECTED_S23B_MANIFEST_SHA256, "S23B target manifest")
    relation_sha = require_sha(S23B_RELATION, EXPECTED_S23B_RELATION_SHA256, "S23B relation inputs")
    manifest = read_json(S23B_MANIFEST)
    cluster_by_id = {str(row["id"]): row for row in manifest["clusters"]}
    if len(cluster_by_id) != 1043:
        raise ValueError("S23B target manifest does not contain 1,043 unique CAMIRA centers")
    targets = []
    with S23B_RELATION.open(encoding="utf-8", newline="") as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            target_id = str(row["target_id"])
            cluster = cluster_by_id.get(target_id)
            if cluster is None:
                raise ValueError(f"unknown S23B target in relation inputs: {target_id}")
            for field, key in (("ra_deg", "ra_deg"), ("dec_deg", "dec_deg"), ("z_cl", "z_cl")):
                if abs(float(row[field]) - float(cluster[key])) > 1e-9:
                    raise ValueError(f"S23B target geometry differs for {target_id}")
            if row["science_measurement"].strip().lower() != "false":
                raise ValueError("S23B relation inputs unexpectedly promote a science measurement")
            targets.append({
                "target_id": target_id, "filter": row["filter"],
                "ra_deg": float(row["ra_deg"]), "dec_deg": float(row["dec_deg"]),
                "z_cl": float(row["z_cl"]),
            })
    if (len(targets) != 427 or len({row["target_id"] for row in targets}) != 427
            or sum(row["filter"] == "i" for row in targets) != 62
            or sum(row["filter"] == "i2" for row in targets) != 365):
        raise ValueError("S23B relation inputs do not contain the pinned 427 targets")
    return targets, {
        str(S23B_MANIFEST): manifest_sha,
        str(S23B_RELATION): relation_sha,
    }


def classify_target(target, old_tree, old_vectors, old_radii_arcsec,
                    source_tree, source_values, cosmology, uv):
    da_mpc = float(cosmology.angular_diameter_distance(target["z_cl"]).value)
    radius_rad = 3.0 / da_mpc
    radius_arcsec = math.degrees(radius_rad) * 3600.0
    center = unit_vector(target["ra_deg"], target["dec_deg"])
    max_old_radius_rad = math.radians(float(np.max(old_radii_arcsec)) / 3600.0)
    possible = old_tree.query_ball_point(
        center, 2.0 * math.sin((radius_rad + max_old_radius_rad) / 2.0)
    )
    if possible:
        old_candidate_vecs = old_vectors[np.asarray(possible, dtype=int)]
        candidate_sep_rad = 2.0 * np.arcsin(np.clip(
            np.linalg.norm(old_candidate_vecs - center[None, :], axis=1) / 2.0,
            0.0, 1.0,
        ))
        candidate_sep_arcsec = np.rad2deg(candidate_sep_rad) * 3600.0
        candidate_radii_arcsec = old_radii_arcsec[np.asarray(possible, dtype=int)]
        intersects = candidate_sep_arcsec <= radius_arcsec + candidate_radii_arcsec
        overlap_indices = np.asarray(possible, dtype=int)[intersects]
        overlap_sep_arcsec = candidate_sep_arcsec[intersects]
        overlap_radii_arcsec = candidate_radii_arcsec[intersects]
        single_contains = bool(np.any(
            overlap_sep_arcsec + radius_arcsec <= overlap_radii_arcsec + 1e-8
        ))
    else:
        overlap_indices = np.array([], dtype=int)
        overlap_sep_arcsec = np.array([], dtype=float)
        overlap_radii_arcsec = np.array([], dtype=float)
        single_contains = False

    disk_points = disk_sample(center, radius_rad, uv)
    covered = np.zeros(len(disk_points), dtype=bool)
    for old_index in overlap_indices:
        dot = np.clip(disk_points @ old_vectors[old_index], -1.0, 1.0)
        separations_arcsec = np.arccos(dot) * 206264.80624709636
        covered |= separations_arcsec <= old_radii_arcsec[old_index]
        if np.all(covered):
            break

    sample_counts = {
        str(count): float(np.mean(covered[:count]))
        for count in (256, 1024, 4096)
    }
    source_radius_chord = 2.0 * math.sin(radius_rad / 2.0)
    source_indices = source_tree.query_ball_point(center, source_radius_chord)
    selected = [source_values[index] for index in source_indices]
    n_raw = len(selected)
    n_shape = sum(row[3] for row in selected)
    n_source = sum(row[4] for row in selected)
    by_zbin = {
        zbin: sum(row[4] and row[5] == zbin for row in selected)
        for zbin in range(5)
    }
    zbin_1_4 = sum(by_zbin[index] for index in range(1, 5))
    bmode0 = sum(row[4] and 1 <= row[5] <= 4 and row[6] == 0 for row in selected)
    bmode1 = sum(row[4] and 1 <= row[5] <= 4 and row[6] == 1 for row in selected)
    if overlap_indices.size:
        nearest_sep = float(np.min(overlap_sep_arcsec))
    else:
        nearest_idx = int(old_tree.query(center, k=1)[1])
        nearest_sep = float(np.rad2deg(
            2.0 * math.asin(min(1.0, np.linalg.norm(old_vectors[nearest_idx] - center) / 2.0))
        ) * 3600.0)
    return {
        **target,
        "aperture_radius_arcsec": radius_arcsec,
        "n_intersecting_old_cones": int(len(overlap_indices)),
        "nearest_old_cone_center_arcsec": nearest_sep,
        "single_old_cone_contains_full_aperture": single_contains,
        "old_cone_union_area_fraction_256": sample_counts["256"],
        "old_cone_union_area_fraction_1024": sample_counts["1024"],
        "old_cone_union_area_fraction_4096": sample_counts["4096"],
        "n_existing_unique_sources_total": n_raw,
        "n_shape_quality": n_shape,
        "n_source_selection": n_source,
        "n_source_selection_zbin0": by_zbin[0],
        "n_source_selection_zbin1": by_zbin[1],
        "n_source_selection_zbin2": by_zbin[2],
        "n_source_selection_zbin3": by_zbin[3],
        "n_source_selection_zbin4": by_zbin[4],
        "n_source_selection_zbin1_4": zbin_1_4,
        "n_source_selection_zbin1_4_bmode_raw0": bmode0,
        "n_source_selection_zbin1_4_bmode_raw1": bmode1,
    }


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def run(output: Path, n_area_samples: int) -> dict:
    if output.exists():
        raise FileExistsError(f"refusing to overwrite existing local audit: {output}")
    if n_area_samples != 4096:
        raise ValueError("this pinned audit uses exactly 4,096 Sobol area samples per target")
    export = verify_export_and_load_sources()
    targets, target_hashes = load_s23b_targets()
    old_lens_rows = list(export["old_lenses"].values())
    old_ra = np.asarray([float(row["ra_deg"]) for row in old_lens_rows])
    old_dec = np.asarray([float(row["dec_deg"]) for row in old_lens_rows])
    old_radii_arcsec = np.asarray([float(row["radius_arcsec"]) for row in old_lens_rows])
    old_vectors = unit_vectors(old_ra, old_dec)
    old_tree = cKDTree(old_vectors)
    source_values = list(export["sources"].values())
    del export["sources"]
    source_vectors = unit_vectors(
        np.asarray([row[1] for row in source_values], dtype=float),
        np.asarray([row[2] for row in source_values], dtype=float),
    )
    source_tree = cKDTree(source_vectors)
    cosmology = FlatLambdaCDM(H0=70.0, Om0=0.3)
    uv = qmc.Sobol(d=2, scramble=False).random_base2(m=12)
    per_target = [
        classify_target(target, old_tree, old_vectors, old_radii_arcsec,
                        source_tree, source_values, cosmology, uv)
        for target in targets
    ]
    fractions = np.asarray([
        row["old_cone_union_area_fraction_4096"] for row in per_target
    ], dtype=float)
    nearest_60_arcsec = sum(row["nearest_old_cone_center_arcsec"] <= 60.0 for row in per_target)
    summary = {
        "schema": "astrafek-s19a-shape-reuse-s23b-local-audit-v1",
        "status": "complete_local_only",
        "network_used": False,
        "data_written_to_remote_service": False,
        "existing_export_rows": export["total_rows"],
        "existing_export_old_lens_count": len(export["old_lenses"]),
        "existing_export_unique_source_ids": export["unique_sources"],
        "repeated_source_pair_rows": export["repeated_source_pair_rows"],
        "source_metadata_conflicts": 0,
        "old_query_cone_membership_max_edge_excess_arcsec": export["max_cone_edge_excess_arcsec"],
        "s23b_target_count": len(per_target),
        "s23b_filter_counts": {
            "i": sum(row["filter"] == "i" for row in per_target),
            "i2": sum(row["filter"] == "i2" for row in per_target),
        },
        "target_cones_with_any_reused_sources": sum(
            row["n_existing_unique_sources_total"] > 0 for row in per_target
        ),
        "target_cones_intersecting_old_query_cones": sum(
            row["n_intersecting_old_cones"] > 0 for row in per_target
        ),
        "target_cones_with_any_old_center_within_60_arcsec": nearest_60_arcsec,
        "target_cones_fully_contained_by_one_old_cone": sum(
            row["single_old_cone_contains_full_aperture"] for row in per_target
        ),
        "estimated_old_cone_union_area_coverage_fraction": {
            "min": float(np.min(fractions)),
            "p10": float(np.quantile(fractions, 0.1)),
            "median": float(np.median(fractions)),
            "p90": float(np.quantile(fractions, 0.9)),
            "max": float(np.max(fractions)),
            "targets_at_least_95_percent": int(np.count_nonzero(fractions >= 0.95)),
            "targets_at_least_99_percent": int(np.count_nonzero(fractions >= 0.99)),
        },
        "selection_counts_are_diagnostic_only": True,
        "source_selection": "0<=sigma_e<=0.4, hypot(e1,e2)<2, mag_a10<=25.5, resolution>=0.3",
        "b_mode_raw_polarity_decoded": False,
        "lensing_efficiency_or_nz_applied": False,
        "photoz_pdf_rows_included": False,
        "science_measurement": False,
        "interpretation": (
            "Existing source rows can only be reused where their 800 old 3-proper-Mpc "
            "query cones cover the S23B aperture. The area fractions are deterministic "
            "Sobol estimates, not exact sky-union areas. This audit does not create a "
            "complete S23B lens-source export or validate a DeltaSigma measurement."
        ),
        "area_estimator": {
            "method": "uniform solid-angle Sobol samples in each 3 proper-Mpc target aperture",
            "samples_per_target": 4096,
            "reported_prefixes": [256, 1024, 4096],
            "old_cone_membership": "exact great-circle distance against each candidate 3 proper-Mpc old query cone",
        },
        "input_sha256": {**export["input_hashes"], **target_hashes},
        "per_target_tsv": "per_target_coverage.tsv",
        "batch_audit_tsv": "batch_integrity.tsv",
    }
    output.mkdir(parents=True, exist_ok=False)
    write_tsv(output / "per_target_coverage.tsv", per_target, OUTPUT_PER_TARGET_FIELDS)
    write_tsv(output / "batch_integrity.tsv", export["batch_audit"], [
        "batch", "file", "rows", "bytes", "sha256",
        "n_lens_ids_with_rows", "max_cone_edge_excess_arcsec",
    ])
    summary["per_target_tsv_sha256"] = sha256(output / "per_target_coverage.tsv")
    summary["batch_audit_tsv_sha256"] = sha256(output / "batch_integrity.tsv")
    with (output / "summary.json").open("x", encoding="utf-8") as stream:
        json.dump(summary, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--area-samples", type=int, default=4096)
    args = parser.parse_args()
    result = run(args.output.resolve(), args.area_samples)
    print(json.dumps({
        "status": result["status"],
        "network_used": result["network_used"],
        "old_pair_rows_verified": result["existing_export_rows"],
        "old_unique_sources": result["existing_export_unique_source_ids"],
        "s23b_targets": result["s23b_target_count"],
        "target_cones_with_reused_sources": result["target_cones_with_any_reused_sources"],
        "target_cones_fully_contained_by_one_old_cone": result["target_cones_fully_contained_by_one_old_cone"],
        "area_coverage_median": result["estimated_old_cone_union_area_coverage_fraction"]["median"],
        "output": str(args.output.resolve()),
    }, ensure_ascii=False))


if __name__ == "__main__":
    main()
