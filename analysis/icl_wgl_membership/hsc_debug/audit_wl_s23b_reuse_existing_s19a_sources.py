#!/usr/bin/env python3
"""Audit the union of already-downloaded S19A lens and random shape rows.

The audit tests whether existing approved S19A source rows cover the 427 S23B
member-profile apertures. It makes no network requests and writes only
per-target counts, geometric coverage summaries, and provenance receipts.
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
from astropy.cosmology import FlatLambdaCDM
from scipy.spatial import cKDTree
from scipy.stats import qmc


HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from audit_wl_s23b_reuse_existing_s19a_shape_rows import (  # noqa: E402
    DATA, EXPECTED_S23B_MANIFEST_SHA256, EXPECTED_S23B_RELATION_SHA256,
    EXPECTED_ROW_PLAN_SHA256, EXPECTED_PREFLIGHT_SHA256,
    EXPECTED_ROW_STATUS_SHA256, EXPECTED_ROW_SUMMARY_SHA256,
    NUMERIC_SOURCE_FIELDS, OUTPUT_COLUMNS as LENS_COLUMNS,
    disk_sample, read_json, require_sha, separation_arcsec_vec,
    unit_vector, unit_vectors, verify_export_and_load_sources,
)
from submit_wl_random_shape_rows import validate_plan  # noqa: E402


RANDOM_PLAN_DIR = DATA / "wl_random_shape_rows_lens_reuse_2026-11-03_01"
RANDOM_PLAN = RANDOM_PLAN_DIR / "random_shape_row_plan.json"
RANDOM_RESULTS = RANDOM_PLAN_DIR / "results"
RANDOM_STATUS = RANDOM_RESULTS / "row_status.json"
RANDOM_SUMMARY = RANDOM_RESULTS / "row_summary.json"
RANDOM_PLAN_SHA256 = "fc883444d2b23663d5f6034d773b76db948add29f1cda0c23384d17cfed2996a"
RANDOM_STATUS_SHA256 = "43e63f84d207a5fd4eccde3ee5cf133bcab4026486dc5b716ab406d45a26e2e9"
RANDOM_SUMMARY_SHA256 = "256e324bf4a8f3ffa31f3493743588b11ecc46ea6a87824380d74b6fe0a380a0"
OLD_AUDIT_DIR = DATA / "wl_s23b_existing_export_overlap_2026-11-07_01"
OLD_AUDIT_SUMMARY_SHA256 = "61026a5d82d47303724f576922415c451589f2516f39337adee3abe021c2c54a"
OLD_AUDIT_TABLE_SHA256 = "8dd03433ddb5603ac19cba03554d53e284b16a2551aba449818f832fccc0282e"
S23B_MANIFEST = DATA / "s23b_icl_target_image_plan_2026-11-07_01/target_manifest.json"
S23B_RELATION = DATA / "s23b_camira_member_outer_icl_relation_2026-11-07_01/cluster_outer_relation_inputs.tsv"
DEFAULT_OUTPUT = DATA / "wl_s23b_existing_source_union_2026-11-07_01"

RANDOM_COLUMNS = [
    "random_center_id", "matched_lens_id", "random_ra_deg", "random_dec_deg",
    "lens_z", "object_id", "source_ra_deg", "source_dec_deg", "e1", "e2",
    "sigma_e", "rms_e", "weight", "m", "c1", "c2", "mag_a10",
    "resolution", "psf_shape11", "psf_shape22", "psf_shape12",
    "b_mode_mask_raw", "source_zbin",
]
SOURCE_FIELDS = LENS_COLUMNS[4:]
RANDOM_SOURCE_FIELDS = [
    "object_id", "source_ra_deg", "source_dec_deg", "e1", "e2", "sigma_e",
    "rms_e", "weight", "m", "c1", "c2", "mag_a10", "resolution",
    "psf_shape11", "psf_shape22", "psf_shape12", "b_mode_mask_raw", "source_zbin",
]
PER_TARGET_FIELDS = [
    "target_id", "filter", "ra_deg", "dec_deg", "z_cl", "aperture_radius_arcsec",
    "n_old_lens_cones_intersecting", "n_old_random_cones_intersecting",
    "single_old_cone_contains_aperture", "single_lens_cone_contains_aperture",
    "single_random_cone_contains_aperture",
    "old_lens_union_area_fraction_256", "old_lens_union_area_fraction_1024",
    "old_lens_union_area_fraction_4096",
    "all_query_union_area_fraction_256", "all_query_union_area_fraction_1024",
    "all_query_union_area_fraction_4096",
    "n_lens_export_unique_sources", "n_random_export_unique_sources",
    "n_both_export_unique_sources", "n_union_unique_sources",
    "n_union_shape_quality", "n_union_source_selection",
    "n_union_source_selection_zbin1_4",
    "n_union_source_selection_zbin1_4_bmode_raw0",
    "n_union_source_selection_zbin1_4_bmode_raw1",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def load_targets() -> list[dict]:
    require_sha(S23B_MANIFEST, EXPECTED_S23B_MANIFEST_SHA256, "S23B target manifest")
    require_sha(S23B_RELATION, EXPECTED_S23B_RELATION_SHA256, "S23B relation inputs")
    manifest = read_json(S23B_MANIFEST)
    cluster_by_id = {str(row["id"]): row for row in manifest["clusters"]}
    targets = []
    with S23B_RELATION.open(encoding="utf-8", newline="") as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            target_id = str(row["target_id"])
            cluster = cluster_by_id.get(target_id)
            if cluster is None:
                raise ValueError(f"unknown S23B profile target {target_id}")
            for field in ("ra_deg", "dec_deg", "z_cl"):
                if abs(float(row[field]) - float(cluster[field])) > 1e-9:
                    raise ValueError(f"changed S23B coordinates/redshift for {target_id}")
            if row["science_measurement"].strip().lower() != "false":
                raise ValueError("S23B input unexpectedly promotes a science measurement")
            targets.append({
                "target_id": target_id, "filter": row["filter"],
                "ra_deg": float(row["ra_deg"]), "dec_deg": float(row["dec_deg"]),
                "z_cl": float(row["z_cl"]),
            })
    if len(targets) != 427 or len({row["target_id"] for row in targets}) != 427:
        raise ValueError("the pinned S23B target list must contain 427 unique centers")
    return targets


def source_value(row: dict[str, str]) -> list:
    numeric = {field: float(row[field]) for field in NUMERIC_SOURCE_FIELDS}
    if not all(math.isfinite(value) for value in numeric.values()):
        raise ValueError("non-finite source field in random-center export")
    ra, dec = numeric["source_ra_deg"], numeric["source_dec_deg"]
    if not (0.0 <= ra < 360.0 and -90.0 <= dec <= 90.0):
        raise ValueError("invalid source coordinates in random-center export")
    mask = int(row["b_mode_mask_raw"])
    zbin = int(row["source_zbin"])
    if str(mask) != row["b_mode_mask_raw"] or mask not in (0, 1):
        raise ValueError("raw B-mode mask is not binary")
    if str(zbin) != row["source_zbin"] or zbin < 0:
        raise ValueError("source tomographic bin is not a non-negative integer")
    payload = "\x1f".join(row[field] for field in RANDOM_SOURCE_FIELDS).encode("utf-8")
    digest = hashlib.blake2b(payload, digest_size=16).digest()
    shape_quality = (
        0.0 <= numeric["sigma_e"] <= 0.4
        and math.hypot(numeric["e1"], numeric["e2"]) < 2.0
    )
    selected = (
        shape_quality and numeric["mag_a10"] <= 25.5
        and numeric["resolution"] >= 0.3
    )
    return [digest, ra, dec, shape_quality, selected, zbin, mask]


def verify_random_export_and_extend_sources(old_sources: dict[str, tuple]) -> dict:
    require_sha(RANDOM_PLAN, RANDOM_PLAN_SHA256, "random shape-row plan")
    require_sha(RANDOM_STATUS, RANDOM_STATUS_SHA256, "random shape-row status")
    require_sha(RANDOM_SUMMARY, RANDOM_SUMMARY_SHA256, "random shape-row summary")
    plan, centers, batches, counts, count_plan_sha, row_plan_sha = validate_plan(RANDOM_PLAN_DIR)
    state = read_json(RANDOM_STATUS)
    summary = read_json(RANDOM_SUMMARY)
    if (row_plan_sha != RANDOM_PLAN_SHA256
            or plan.get("status") != "prepared_not_submitted"
            or plan.get("n_random_centers") != 3100
            or plan.get("n_matched_lenses") != 620
            or plan.get("expected_random_source_pair_rows") != 9010324
            or plan.get("output_columns") != RANDOM_COLUMNS
            or state.get("status") != "complete"
            or state.get("plan_sha256") != row_plan_sha
            or state.get("source_shape_rows_downloaded") is not True
            or summary.get("status") != "complete"
            or summary.get("plan_sha256") != row_plan_sha
            or summary.get("n_source_pairs") != 9010324
            or summary.get("n_random_centers") != len(centers)
            or summary.get("science_measurement") is not False):
        raise ValueError("random shape-row export plan or receipt differs from its pinned scope")
    for flag in ("credentials_saved", "photoz_pdf_rows_included", "images_returned",
                 "science_measurement"):
        if state.get(flag) is not False:
            raise ValueError(f"random shape-row run has unexpected scope flag {flag}")

    center_by_id = {str(row["random_center_id"]): row for row in centers}
    if len(center_by_id) != 3100:
        raise ValueError("random-center plan does not contain 3,100 unique centers")
    jobs = state.get("jobs", {})
    if len(jobs) != 124 or len(batches) != 124:
        raise ValueError("random shape-row export does not contain all 124 batches")

    for object_id, value in list(old_sources.items()):
        old_sources[object_id] = [*value, True, False]

    total_rows = 0
    total_bytes = 0
    random_unique_ids = 0
    random_only_new_ids = 0
    random_ids_overlapping_old = 0
    random_repeat_rows = 0
    conflict_count = 0
    max_cone_excess = -math.inf
    batch_records = []
    batch_size = int(plan["batch_size"])

    for index, batch in enumerate(batches):
        job_name = f"random_shape_rows_{index:03d}"
        job = jobs.get(job_name)
        path = RANDOM_RESULTS / f"{job_name}.csv"
        sql_path = Path(batch["path"])
        if (job is None or job.get("status") != "complete"
                or not path.is_file()
                or sql_path.parent != RANDOM_PLAN_DIR
                or sql_path.name != f"random_shape_rows_{index:03d}.sql"):
            raise ValueError(f"missing or incomplete approved random batch {job_name}")
        actual_sha = sha256(path)
        actual_bytes = path.stat().st_size
        if (actual_sha != job.get("sha256")
                or actual_bytes != int(job.get("bytes", -1))):
            raise ValueError(f"random source-row receipt failed for {job_name}")

        expected_centers = centers[index * batch_size:(index + 1) * batch_size]
        center_ids = {str(row["random_center_id"]) for row in expected_centers}
        if (len(expected_centers) != int(batch.get("n_random_centers", -1))
                or int(batch.get("expected_pair_rows_from_count_only", -2))
                != sum(int(counts[row["random_center_id"]]["n_shape_sources"])
                       for row in expected_centers)):
            raise ValueError(f"random batch count plan differs for {job_name}")
        seen_pairs = set()
        coordinates_by_center: dict[str, tuple[list, list]] = {}
        rows_in_batch = 0
        with path.open(encoding="utf-8-sig", newline="") as stream:
            reader = csv.reader(stream)
            raw_header = next(reader, None)
            header = [value.lstrip("#").strip() for value in (raw_header or [])]
            if header != RANDOM_COLUMNS:
                raise ValueError(f"random-row CSV schema changed in {job_name}")
            for values in reader:
                if len(values) != len(RANDOM_COLUMNS):
                    raise ValueError(f"malformed random source row in {job_name}")
                row = dict(zip(RANDOM_COLUMNS, values))
                center_id = str(row["random_center_id"])
                center = center_by_id.get(center_id)
                if center_id not in center_ids or center is None:
                    raise ValueError(f"unrequested random center in {job_name}")
                if (str(row["matched_lens_id"]) != str(center["matched_lens_id"])
                        or abs(float(row["random_ra_deg"]) - float(center["ra_deg"])) > 1e-9
                        or abs(float(row["random_dec_deg"]) - float(center["dec_deg"])) > 1e-9
                        or abs(float(row["lens_z"]) - float(center["z_lens"])) > 1e-10):
                    raise ValueError(f"random-center/lens metadata changed in {job_name}")
                object_id = str(row["object_id"])
                if not object_id or (center_id, object_id) in seen_pairs:
                    raise ValueError(f"empty or duplicate center-source pair in {job_name}")
                seen_pairs.add((center_id, object_id))
                if row["b_mode_mask_raw"] not in ("0", "1"):
                    raise ValueError(f"invalid B-mode field in {job_name}")
                if not row["source_zbin"].isdigit():
                    raise ValueError(f"invalid redshift bin in {job_name}")
                value = source_value(row)
                previous = old_sources.get(object_id)
                if previous is None:
                    old_sources[object_id] = [*value, False, True]
                    random_unique_ids += 1
                    random_only_new_ids += 1
                else:
                    if previous[:7] != value:
                        conflict_count += 1
                        raise ValueError(f"source fields conflict for object_id {object_id}")
                    if not previous[8]:
                        random_unique_ids += 1
                        if previous[7]:
                            random_ids_overlapping_old += 1
                    else:
                        random_repeat_rows += 1
                    previous[8] = True

                ra_list, dec_list = coordinates_by_center.setdefault(center_id, ([], []))
                ra_list.append(value[1])
                dec_list.append(value[2])
                rows_in_batch += 1
        if rows_in_batch != int(job.get("rows", -1)):
            raise ValueError(f"random source-row count differs from receipt in {job_name}")
        if rows_in_batch != int(batch.get("expected_pair_rows_from_count_only", -2)):
            raise ValueError(f"random source-row count differs from count preflight in {job_name}")

        for center_id, (ra_values, dec_values) in coordinates_by_center.items():
            center = center_by_id[center_id]
            separations = separation_arcsec_vec(
                float(center["ra_deg"]), float(center["dec_deg"]),
                np.asarray(ra_values), np.asarray(dec_values),
            )
            excess = separations - float(center["radius_arcsec"])
            batch_center_excess = float(np.max(excess))
            if batch_center_excess > 0.05:
                raise ValueError(f"source rows lie outside random query cone {center_id}")
            max_cone_excess = max(max_cone_excess, batch_center_excess)

        total_rows += rows_in_batch
        total_bytes += actual_bytes
        batch_records.append({
            "batch": index, "file": path.name, "rows": rows_in_batch,
            "bytes": actual_bytes, "sha256": actual_sha,
            "centers_with_rows": len(coordinates_by_center),
        })
        if (index + 1) % 20 == 0:
            print(f"verified {index + 1}/124 random export batches", flush=True)

    if (total_rows != 9010324 or total_rows != int(summary["n_source_pairs"])
            or total_bytes != int(summary["n_csv_bytes"])
            or conflict_count != 0):
        raise ValueError("random-row aggregate or source-field consistency check failed")
    return {
        "plan": plan, "centers": centers, "row_plan_sha256": row_plan_sha,
        "count_plan_sha256": count_plan_sha, "batch_records": batch_records,
        "total_rows": total_rows, "total_bytes": total_bytes,
        "random_unique_ids": random_unique_ids,
        "random_only_new_ids": random_only_new_ids,
        "random_ids_overlapping_old": random_ids_overlapping_old,
        "random_repeat_rows": random_repeat_rows,
        "source_metadata_conflicts": conflict_count,
        "max_random_cone_edge_excess_arcsec": max_cone_excess,
    }


def area_union(center: np.ndarray, radius_rad: float, uv: np.ndarray,
               tree: cKDTree, cone_vectors: np.ndarray, cone_radii_arcsec: np.ndarray,
               is_random: np.ndarray) -> dict:
    max_radius_rad = math.radians(float(np.max(cone_radii_arcsec)) / 3600.0)
    possible = tree.query_ball_point(
        center, 2.0 * math.sin((radius_rad + max_radius_rad) / 2.0)
    )
    if possible:
        possible = np.asarray(possible, dtype=int)
        vecs = cone_vectors[possible]
        sep_arcsec = np.rad2deg(2.0 * np.arcsin(np.clip(
            np.linalg.norm(vecs - center[None, :], axis=1) / 2.0, 0.0, 1.0
        ))) * 3600.0
        intersects = sep_arcsec <= math.degrees(radius_rad) * 3600.0 + cone_radii_arcsec[possible]
        overlapping = possible[intersects]
        overlapping_sep = sep_arcsec[intersects]
        radii = cone_radii_arcsec[overlapping]
        random_flags = is_random[overlapping]
        target_radius_arcsec = math.degrees(radius_rad) * 3600.0
        single_lens = bool(np.any((~random_flags) & (overlapping_sep + target_radius_arcsec <= radii + 1e-8)))
        single_random = bool(np.any(random_flags & (overlapping_sep + target_radius_arcsec <= radii + 1e-8)))
        single_any = single_lens or single_random
    else:
        overlapping = np.array([], dtype=int)
        overlapping_sep = np.array([], dtype=float)
        random_flags = np.array([], dtype=bool)
        single_lens = single_random = single_any = False
    points = disk_sample(center, radius_rad, uv)
    covered_lens = np.zeros(len(points), dtype=bool)
    covered_all = np.zeros(len(points), dtype=bool)
    for cone_index in overlapping:
        if np.all(covered_all) and np.all(covered_lens):
            break
        distances = np.arccos(np.clip(points @ cone_vectors[cone_index], -1.0, 1.0)) * 206264.80624709636
        in_cone = distances <= cone_radii_arcsec[cone_index]
        covered_all |= in_cone
        if not is_random[cone_index]:
            covered_lens |= in_cone
    return {
        "n_lens_intersections": int(np.count_nonzero(~random_flags)),
        "n_random_intersections": int(np.count_nonzero(random_flags)),
        "single_any": single_any, "single_lens": single_lens,
        "single_random": single_random,
        "lens_fraction": {str(n): float(np.mean(covered_lens[:n])) for n in (256, 1024, 4096)},
        "all_fraction": {str(n): float(np.mean(covered_all[:n])) for n in (256, 1024, 4096)},
    }


def summarize_target(target: dict, source_tree: cKDTree, source_values: list,
                     cosmology, cone_tree: cKDTree, cone_vectors: np.ndarray,
                     cone_radii_arcsec: np.ndarray, is_random_cone: np.ndarray,
                     uv: np.ndarray, old_target_counts: dict) -> dict:
    da_mpc = float(cosmology.angular_diameter_distance(target["z_cl"]).value)
    radius_rad = 3.0 / da_mpc
    radius_arcsec = math.degrees(radius_rad) * 3600.0
    center = unit_vector(target["ra_deg"], target["dec_deg"])
    geometry = area_union(center, radius_rad, uv, cone_tree, cone_vectors,
                          cone_radii_arcsec, is_random_cone)
    indices = source_tree.query_ball_point(center, 2.0 * math.sin(radius_rad / 2.0))
    rows = [source_values[index] for index in indices]
    old_count = sum(row[6] for row in rows)
    random_count = sum(row[7] for row in rows)
    both_count = sum(row[6] and row[7] for row in rows)
    if old_count != int(old_target_counts[target["target_id"]]):
        raise ValueError(f"union audit does not reproduce old-only source count for {target['target_id']}")
    selected = sum(row[3] for row in rows)
    zbin_1_4 = [row for row in rows if row[3] and 1 <= row[4] <= 4]
    return {
        **target,
        "aperture_radius_arcsec": radius_arcsec,
        "n_old_lens_cones_intersecting": geometry["n_lens_intersections"],
        "n_old_random_cones_intersecting": geometry["n_random_intersections"],
        "single_old_cone_contains_aperture": geometry["single_any"],
        "single_lens_cone_contains_aperture": geometry["single_lens"],
        "single_random_cone_contains_aperture": geometry["single_random"],
        "old_lens_union_area_fraction_256": geometry["lens_fraction"]["256"],
        "old_lens_union_area_fraction_1024": geometry["lens_fraction"]["1024"],
        "old_lens_union_area_fraction_4096": geometry["lens_fraction"]["4096"],
        "all_query_union_area_fraction_256": geometry["all_fraction"]["256"],
        "all_query_union_area_fraction_1024": geometry["all_fraction"]["1024"],
        "all_query_union_area_fraction_4096": geometry["all_fraction"]["4096"],
        "n_lens_export_unique_sources": old_count,
        "n_random_export_unique_sources": random_count,
        "n_both_export_unique_sources": both_count,
        "n_union_unique_sources": len(rows),
        "n_union_shape_quality": sum(row[2] for row in rows),
        "n_union_source_selection": selected,
        "n_union_source_selection_zbin1_4": len(zbin_1_4),
        "n_union_source_selection_zbin1_4_bmode_raw0": sum(row[5] == 0 for row in zbin_1_4),
        "n_union_source_selection_zbin1_4_bmode_raw1": sum(row[5] == 1 for row in zbin_1_4),
    }


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def run(output: Path, samples: int) -> dict:
    if output.exists():
        raise FileExistsError(f"refusing to overwrite existing S23B source-union audit: {output}")
    if samples != 4096:
        raise ValueError("the pinned geometric audit uses exactly 4,096 Sobol samples")

    old = verify_export_and_load_sources()
    old_summary_path = OLD_AUDIT_DIR / "summary.json"
    old_table_path = OLD_AUDIT_DIR / "per_target_coverage.tsv"
    require_sha(old_summary_path, OLD_AUDIT_SUMMARY_SHA256, "previous S19A overlap summary")
    require_sha(old_table_path, OLD_AUDIT_TABLE_SHA256, "previous S19A per-target table")
    previous_targets = {}
    with old_table_path.open(encoding="utf-8", newline="") as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            previous_targets[row["target_id"]] = int(row["n_existing_unique_sources_total"])

    random_audit = verify_random_export_and_extend_sources(old["sources"])
    targets = load_targets()

    # Convert the identity audit map to compact source records before building
    # the position tree; no raw source IDs are written to the result tables.
    source_values = [
        (row[1], row[2], row[3], row[4], row[5], row[6], row[7], row[8])
        for row in old["sources"].values()
    ]
    source_vectors = unit_vectors(
        np.asarray([row[0] for row in source_values], dtype=float),
        np.asarray([row[1] for row in source_values], dtype=float),
    )
    source_tree = cKDTree(source_vectors)
    del old["sources"], source_vectors

    old_lenses = list(old["old_lenses"].values())
    random_centers = random_audit["centers"]
    cone_ra = np.asarray(
        [float(row["ra_deg"]) for row in old_lenses]
        + [float(row["ra_deg"]) for row in random_centers], dtype=float
    )
    cone_dec = np.asarray(
        [float(row["dec_deg"]) for row in old_lenses]
        + [float(row["dec_deg"]) for row in random_centers], dtype=float
    )
    cone_radii_arcsec = np.asarray(
        [float(row["radius_arcsec"]) for row in old_lenses]
        + [float(row["radius_arcsec"]) for row in random_centers], dtype=float
    )
    is_random_cone = np.zeros(len(cone_ra), dtype=bool)
    is_random_cone[len(old_lenses):] = True
    cone_vectors = unit_vectors(cone_ra, cone_dec)
    cone_tree = cKDTree(cone_vectors)
    cosmology = FlatLambdaCDM(H0=70.0, Om0=0.3, Tcmb0=0)
    uv = qmc.Sobol(d=2, scramble=False).random_base2(m=12)

    per_target = [
        summarize_target(target, source_tree, source_values, cosmology,
                         cone_tree, cone_vectors, cone_radii_arcsec,
                         is_random_cone, uv, previous_targets)
        for target in targets
    ]
    old_fraction = np.asarray([row["old_lens_union_area_fraction_4096"] for row in per_target])
    all_fraction = np.asarray([row["all_query_union_area_fraction_4096"] for row in per_target])
    total_unique = int(len(source_values))
    summary = {
        "schema": "astrafek-s19a-existing-source-union-s23b-audit-v1",
        "status": "complete_local_only",
        "network_used": False,
        "data_written_to_remote_service": False,
        "old_lens_query_pair_rows": old["total_rows"],
        "old_lens_query_unique_sources": old["unique_sources"],
        "old_random_query_pair_rows": random_audit["total_rows"],
        "old_random_query_csv_bytes": random_audit["total_bytes"],
        "old_random_query_unique_sources": random_audit["random_unique_ids"],
        "old_random_sources_new_to_lens_export": random_audit["random_only_new_ids"],
        "old_random_sources_overlapping_lens_export": random_audit["random_ids_overlapping_old"],
        "all_query_union_unique_source_ids": total_unique,
        "source_metadata_conflicts": random_audit["source_metadata_conflicts"],
        "random_query_max_cone_edge_excess_arcsec": random_audit["max_random_cone_edge_excess_arcsec"],
        "s23b_target_count": len(per_target),
        "s23b_filter_counts": {
            "i": sum(row["filter"] == "i" for row in per_target),
            "i2": sum(row["filter"] == "i2" for row in per_target),
        },
        "targets_with_any_old_lens_sources": sum(row["n_lens_export_unique_sources"] > 0 for row in per_target),
        "targets_with_any_old_random_sources": sum(row["n_random_export_unique_sources"] > 0 for row in per_target),
        "targets_with_any_union_sources": sum(row["n_union_unique_sources"] > 0 for row in per_target),
        "targets_with_no_union_sources": sum(row["n_union_unique_sources"] == 0 for row in per_target),
        "targets_with_one_old_cone_containing_full_aperture": sum(
            row["single_old_cone_contains_aperture"] for row in per_target
        ),
        "targets_with_one_old_lens_cone_containing_full_aperture": sum(
            row["single_lens_cone_contains_aperture"] for row in per_target
        ),
        "targets_with_one_old_random_cone_containing_full_aperture": sum(
            row["single_random_cone_contains_aperture"] for row in per_target
        ),
        "old_lens_union_area_fraction": {
            "median": float(np.median(old_fraction)),
            "p90": float(np.quantile(old_fraction, 0.9)),
            "max": float(np.max(old_fraction)),
        },
        "all_query_union_area_fraction": {
            "median": float(np.median(all_fraction)),
            "p90": float(np.quantile(all_fraction, 0.9)),
            "max": float(np.max(all_fraction)),
            "targets_at_least_95_percent": int(np.count_nonzero(all_fraction >= 0.95)),
            "targets_at_least_99_percent": int(np.count_nonzero(all_fraction >= 0.99)),
            "targets_at_least_99p9_percent": int(np.count_nonzero(all_fraction >= 0.999)),
        },
        "area_estimator": {
            "method": "uniform solid-angle Sobol samples in each S23B 3 proper-Mpc aperture",
            "samples_per_target": samples,
            "reported_prefixes": [256, 1024, 4096],
            "query_union": "800 old S19A lens cones plus 3,100 approved S19A random-center cones",
        },
        "selection_counts_are_diagnostic_only": True,
        "source_selection": "0<=sigma_e<=0.4, hypot(e1,e2)<2, mag_a10<=25.5, resolution>=0.3",
        "b_mode_raw_polarity_decoded": False,
        "lensing_efficiency_or_nz_applied": False,
        "photoz_pdf_rows_included": False,
        "science_measurement": False,
        "interpretation": (
            "Existing source rows from both approved S19A exports were deduplicated by object_id "
            "and searched locally inside each S23B aperture. Geometric query-cone coverage is an "
            "upper bound on reusable catalog coverage; it does not establish uniform HSC depth, "
            "matched random support, or a calibrated DeltaSigma result."
        ),
        "input_sha256": {
            "old_lens_row_plan": EXPECTED_ROW_PLAN_SHA256,
            "old_lens_row_status": EXPECTED_ROW_STATUS_SHA256,
            "old_lens_row_summary": EXPECTED_ROW_SUMMARY_SHA256,
            "previous_lens_only_audit_summary": OLD_AUDIT_SUMMARY_SHA256,
            "previous_lens_only_audit_table": OLD_AUDIT_TABLE_SHA256,
            "random_row_plan": random_audit["row_plan_sha256"],
            "random_row_status": RANDOM_STATUS_SHA256,
            "random_row_summary": RANDOM_SUMMARY_SHA256,
            "random_count_plan": random_audit["count_plan_sha256"],
            "s23b_target_manifest": EXPECTED_S23B_MANIFEST_SHA256,
            "s23b_relation_inputs": EXPECTED_S23B_RELATION_SHA256,
        },
        "per_target_tsv": "per_target_source_union.tsv",
        "batch_audit_tsv": "random_batch_integrity.tsv",
    }
    output.mkdir(parents=True, exist_ok=False)
    write_tsv(output / "per_target_source_union.tsv", per_target, PER_TARGET_FIELDS)
    write_tsv(output / "random_batch_integrity.tsv", random_audit["batch_records"], [
        "batch", "file", "rows", "bytes", "sha256", "centers_with_rows",
    ])
    summary["per_target_tsv_sha256"] = sha256(output / "per_target_source_union.tsv")
    summary["random_batch_audit_tsv_sha256"] = sha256(output / "random_batch_integrity.tsv")
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
        "status": result["status"], "network_used": result["network_used"],
        "old_lens_rows": result["old_lens_query_pair_rows"],
        "old_random_rows": result["old_random_query_pair_rows"],
        "random_unique_sources": result["old_random_query_unique_sources"],
        "s23b_targets": result["s23b_target_count"],
        "targets_with_any_union_sources": result["targets_with_any_union_sources"],
        "targets_with_one_old_cone_containing_full_aperture": result[
            "targets_with_one_old_cone_containing_full_aperture"
        ],
        "all_query_area_coverage": result["all_query_union_area_fraction"],
        "output": str(args.output.resolve()),
    }, ensure_ascii=False))


if __name__ == "__main__":
    main()
