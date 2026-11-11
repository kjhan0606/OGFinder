#!/usr/bin/env python3
"""Prepare, but never submit, an S16A count-only plan for five WL controls/lens.

The plan screens the saved, provisional five-random-centers-per-lens set for
S16A catalog coverage. It writes SQL and provenance locally; it does not
contact HSC CAS or return source-level data.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path

from astropy.cosmology import FlatLambdaCDM

from prepare_s23b_s16a_shape_count_plan import (
    EXPECTED_TARGET_SHA256,
    TARGET_MANIFEST,
    canonical_sha,
    sql_literal,
    sql_number,
)

DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
SELECTION_DIR = DATA / "s23b_s16a_expansion_wl_random_centers_2026-11-11_02"
SELECTION_SUMMARY = SELECTION_DIR / "summary.json"
RANDOM_CENTERS = SELECTION_DIR / "random_centers.tsv"
TARGET_SELECTION = SELECTION_DIR / "target_selection.tsv"
OUTPUT = DATA / "wl_s23b_s16a_expansion_random_count_plan_2026-11-11_01"

EXPECTED_SELECTION_SUMMARY_SHA256 = "197dd377f54dd645a0acdaa3059a50c5c6a6a8d87e494faa75f60970981bce33"
EXPECTED_RANDOM_CENTERS_SHA256 = "d63bc851903882e8ad21bc85de522ac2722a49d7b8533f7fc1ea8335fe9cac24"
EXPECTED_TARGET_SELECTION_SHA256 = "904744b0d6efd0c7ded8570fe1bc941b9e97367aaaae5282df5448151da5a2f6"
EXPECTED_CANDIDATE_AUDIT_SHA256 = "ac9cb72ca769b19ace23e2b90de3f6dd5676c9f3813a525e8f8319f6fae6cb14"
EXPECTED_TARGETS = 114
EXPECTED_RANDOM_CENTERS = 570
CENTERS_PER_LENS = 5
RADIUS_MPC = 3.0
BATCH_SIZE = 25
QUERY_TABLES = [
    "s16a_wide.meas2",
    "s16a_wide.weaklensing_hsm_regauss",
    "s16a_wide.meas",
]
OUTPUT_COLUMNS = ["random_center_id", "n_shape_sources"]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with Path(path).open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def count_query(center: dict) -> str:
    return (
        f"SELECT {sql_literal(center['random_center_id'])} AS random_center_id, "
        "COUNT(DISTINCT b.object_id) AS n_shape_sources\n"
        "FROM s16a_wide.meas2 AS a\n"
        "INNER JOIN s16a_wide.weaklensing_hsm_regauss AS b USING (object_id)\n"
        "INNER JOIN s16a_wide.meas AS c USING (object_id)\n"
        f"WHERE coneSearch(c.icoord, {sql_number(center['ra_deg'])}, "
        f"{sql_number(center['dec_deg'])}, {sql_number(center['radius_arcsec'])})"
    )


def prepare(output: Path = OUTPUT) -> dict:
    output = output.resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite S16A random-count plan: {output}")
    for path, expected, label in (
        (SELECTION_SUMMARY, EXPECTED_SELECTION_SUMMARY_SHA256, "selection summary"),
        (RANDOM_CENTERS, EXPECTED_RANDOM_CENTERS_SHA256, "random-center table"),
        (TARGET_SELECTION, EXPECTED_TARGET_SELECTION_SHA256, "target-selection ledger"),
        (TARGET_MANIFEST, EXPECTED_TARGET_SHA256, "S23B target manifest"),
    ):
        if not path.is_file() or sha256(path) != expected:
            raise ValueError(f"pinned {label} is missing or changed")

    summary = json.loads(SELECTION_SUMMARY.read_text(encoding="utf-8"))
    if (summary.get("schema") != "astrafek-s23b-s16a-expansion-wl-random-centers-v2"
            or summary.get("status") != "complete_local_provisional_selection"
            or summary.get("science_measurement") is not False
            or summary.get("network_used") is not False
            or summary.get("cas_query_submitted") is not False
            or summary.get("source_rows_downloaded") is not False
            or summary.get("images_read") is not False
            or summary.get("n_candidate_targets") != 116
            or summary.get("n_targets_with_five_globally_unique_random_centers") != EXPECTED_TARGETS
            or summary.get("n_random_centers") != EXPECTED_RANDOM_CENTERS
            or summary.get("controls_ready_for_science") is not False
            or summary.get("random_center_table_sha256") != EXPECTED_RANDOM_CENTERS_SHA256
            or summary.get("target_selection_table_sha256") != EXPECTED_TARGET_SELECTION_SHA256
            or summary.get("candidate_audit_csv_sha256") != EXPECTED_CANDIDATE_AUDIT_SHA256):
        raise ValueError("saved random-center selection has an unexpected scope or status")

    manifest = json.loads(TARGET_MANIFEST.read_text(encoding="utf-8"))
    if (manifest.get("schema") != "hsc-bcg-diffuse-sample-v1"
            or manifest.get("science_measurement") is not False):
        raise ValueError("unexpected S23B target-manifest schema or status")
    targets = {str(row["id"]): row for row in manifest.get("clusters", [])}
    if len(targets) != 1043:
        raise ValueError("pinned target manifest does not contain 1,043 unique clusters")

    selection_rows = read_tsv(TARGET_SELECTION)
    selection_by_lens = {row["matched_lens_id"]: row for row in selection_rows}
    if len(selection_rows) != 116 or len(selection_by_lens) != 116:
        raise ValueError("target-selection ledger must contain 116 unique candidate lenses")
    selected_ids = {
        lens_id for lens_id, row in selection_by_lens.items()
        if row["selection_status"] == "selected"
    }
    if len(selected_ids) != EXPECTED_TARGETS:
        raise ValueError("target-selection ledger does not select exactly 114 lenses")
    for lens_id, row in selection_by_lens.items():
        if lens_id not in targets or str(targets[lens_id].get("sample_role")) != "target":
            raise ValueError(f"selection ledger refers to an unknown target: {lens_id}")
        if (int(row["n_selected_random_centers"]) != (CENTERS_PER_LENS if lens_id in selected_ids else 0)
                or int(row["n_eligible_random_candidates"]) < int(row["n_selected_random_centers"])):
            raise ValueError(f"invalid selected/eligible random-center counts: {lens_id}")

    h0 = float(manifest["h0"])
    omega_m = float(manifest["omega_m"])
    cosmology = FlatLambdaCDM(H0=h0, Om0=omega_m, Tcmb0=0)
    centers = read_tsv(RANDOM_CENTERS)
    if len(centers) != EXPECTED_RANDOM_CENTERS:
        raise ValueError("random-center table does not contain exactly 570 rows")
    required_columns = [
        "random_center_id", "matched_lens_id", "matched_catalog_id", "icl_filter",
        "ra_deg", "dec_deg", "z_lens", "target_s16a_raw_count", "random_object_id",
        "control_slot", "offset_from_target_mpc", "offset_from_target_arcmin",
        "match_score", "tract", "patch", "skymap_id", "delta_depth_mag",
        "seeing_ratio", "delta_ebv_sfd", "cluster_overlap_2mpc",
    ]
    if list(centers[0]) != required_columns:
        raise ValueError("random-center table columns changed")

    center_ids: set[str] = set()
    random_object_ids: set[str] = set()
    sky_positions: set[tuple[float, float]] = set()
    per_lens: dict[str, list[dict]] = {lens_id: [] for lens_id in selected_ids}
    for row in centers:
        center_id = str(row["random_center_id"])
        lens_id = str(row["matched_lens_id"])
        if lens_id not in selected_ids or center_id in center_ids:
            raise ValueError("random center has an unknown lens or duplicate center ID")
        target = targets[lens_id]
        ra, dec, z = float(row["ra_deg"]), float(row["dec_deg"]), float(row["z_lens"])
        offset_mpc = float(row["offset_from_target_mpc"])
        target_count = int(row["target_s16a_raw_count"])
        match_score = float(row["match_score"])
        if (not (math.isfinite(ra) and 0 <= ra < 360
                 and math.isfinite(dec) and -90 <= dec <= 90)
                or not math.isclose(z, float(target["z_cl"]), rel_tol=0, abs_tol=1e-12)
                or str(row["matched_catalog_id"]) != str(target["catalog_id"])
                or str(row["icl_filter"]) != str(selection_by_lens[lens_id]["icl_filter"])
                or target_count < 1000
                or target_count != int(selection_by_lens[lens_id]["target_s16a_raw_count"])
                or int(row["cluster_overlap_2mpc"]) != 0
                or not math.isfinite(offset_mpc) or not offset_mpc > 6.0
                or not math.isfinite(match_score)):
            raise ValueError(f"random-center metadata or candidate cuts changed: {center_id}")
        random_object_id = str(row["random_object_id"])
        if random_object_id in random_object_ids:
            raise ValueError("candidate random-object ID was reused")
        position = (round(ra, 12), round(dec, 12))
        if position in sky_positions:
            raise ValueError("random-center sky position was reused")
        expected_center_id = "RC-" + hashlib.sha256(
            f"{lens_id}\0{random_object_id}\0{ra:.12f}\0{dec:.12f}".encode()
        ).hexdigest()[:20]
        if center_id != expected_center_id:
            raise ValueError(f"random-center ID does not match its pinned identity: {center_id}")

        da_mpc = float(cosmology.angular_diameter_distance(z).value)
        radius_arcsec = math.degrees(RADIUS_MPC / da_mpc) * 3600.0
        center = {
            "random_center_id": center_id,
            "matched_lens_id": lens_id,
            "icl_filter": str(row["icl_filter"]),
            "ra_deg": ra,
            "dec_deg": dec,
            "z_lens": z,
            "radius_mpc_proper": RADIUS_MPC,
            "radius_arcsec": radius_arcsec,
            "control_slot": int(row["control_slot"]),
            "random_object_id": random_object_id,
        }
        center_ids.add(center_id)
        random_object_ids.add(random_object_id)
        sky_positions.add(position)
        per_lens[lens_id].append(center)

    for lens_id, lens_centers in per_lens.items():
        if (len(lens_centers) != CENTERS_PER_LENS
                or sorted(row["control_slot"] for row in lens_centers) != [1, 2, 3, 4, 5]):
            raise ValueError(f"lens does not have exactly five unique random centers: {lens_id}")
    centers = [center for lens_id in sorted(per_lens) for center in per_lens[lens_id]]
    centers.sort(key=lambda row: (row["matched_lens_id"], row["control_slot"]))

    batches = []
    batch_sql = []
    for index, start in enumerate(range(0, len(centers), BATCH_SIZE)):
        chunk = centers[start:start + BATCH_SIZE]
        sql = "\nUNION ALL\n".join(count_query(center) for center in chunk) + ";\n"
        name = f"s16a_wl_random_counts_{index:03d}.sql"
        batch_sql.append((name, sql))
        batches.append({
            "name": name,
            "n_centers": len(chunk),
            "first_random_center_id": chunk[0]["random_center_id"],
            "last_random_center_id": chunk[-1]["random_center_id"],
            "sql_bytes": len(sql.encode("utf-8")),
            "sha256": hashlib.sha256(sql.encode("utf-8")).hexdigest(),
        })

    header_bytes = len(("\t".join(OUTPUT_COLUMNS) + "\n").encode("utf-8"))
    max_result_bytes = header_bytes + sum(
        len(center["random_center_id"].encode("utf-8")) + 1 + 19 + 1
        for center in centers
    )
    plan = {
        "schema": "astrafek-s23b-s16a-wl-random-count-plan-v1",
        "status": "prepared_not_submitted",
        "catalog_release_documentation": "HSC S16A shape catalog; CAS namespace s16a_wide",
        "cas_release_version": "pdr1",
        "selection_summary": str(SELECTION_SUMMARY),
        "selection_summary_sha256": sha256(SELECTION_SUMMARY),
        "random_center_table": str(RANDOM_CENTERS),
        "random_center_table_sha256": sha256(RANDOM_CENTERS),
        "target_selection_ledger": str(TARGET_SELECTION),
        "target_selection_ledger_sha256": sha256(TARGET_SELECTION),
        "source_target_manifest": str(TARGET_MANIFEST),
        "source_target_manifest_sha256": sha256(TARGET_MANIFEST),
        "candidate_audit_csv_sha256": summary["candidate_audit_csv_sha256"],
        "target_scope": "114 of the 116 existing ICL-paired S16A high-raw-count expansion candidates; two i-band lenses lack five eligible random centers",
        "random_center_selection": "Five globally unique saved candidate positions per lens; center flags and metadata match pass; CAMIRA overlap within 2 proper-Mpc is zero; center is >6 proper-Mpc from its matched lens; candidates sorted by match score and random-object ID before deterministic assignment",
        "n_candidate_lenses": 116,
        "n_selected_lenses": len(per_lens),
        "n_excluded_lenses": 116 - len(per_lens),
        "n_random_centers": len(centers),
        "centers_per_lens": CENTERS_PER_LENS,
        "radius_mpc_proper": RADIUS_MPC,
        "cosmology": {"H0_km_s_mpc": h0, "Omega_m": omega_m},
        "n_query_batches": len(batches),
        "batch_size": BATCH_SIZE,
        "query_tables": QUERY_TABLES,
        "spatial_function": "coneSearch(c.icoord, ra_deg, dec_deg, radius_arcsec)",
        "output_columns": OUTPUT_COLUMNS,
        "query_semantics": "COUNT(DISTINCT S16A weaklensing_hsm_regauss object_id) in each matched random-center 3 proper-Mpc cone at its lens redshift, before source-z, shape-quality, B-mode, mask, radial, and estimator cuts",
        "estimated_result_rows": len(centers),
        "estimated_result_bytes_upper_bound": max_result_bytes,
        "estimated_result_scope": "570 random-center IDs and integer catalog counts only; no source IDs, coordinates returned, shapes, photo-z values, images, or PSF values",
        "coordinate_transfer": "If submitted, the 570 saved random-center coordinates and derived cone radii are sent to HSC CAS as query predicates; the result contains only the center ID and one integer count per center",
        "network_used": False,
        "query_submitted": False,
        "source_rows_downloaded": False,
        "images_returned": False,
        "photoz_returned": False,
        "science_measurement": False,
        "controls_ready_for_science": False,
        "execution_requires_explicit_user_approval": True,
        "interpretation_limit": "Raw catalog availability only; does not establish usable source support or a weak-lensing measurement. These provisional random centers have not passed image/PSF/mask QA.",
        "release_caveat": "S16A and S19A shape rows use different processing and calibration; do not combine them without separate validation.",
        "batches": batches,
        "centers": centers,
    }
    plan["query_plan_sha256"] = canonical_sha(plan)

    output.mkdir(parents=True, exist_ok=False)
    for name, sql in batch_sql:
        (output / name).write_text(sql, encoding="utf-8")
    plan_path = output / "shape_count_plan.json"
    plan_path.write_text(json.dumps(plan, ensure_ascii=False, indent=2, allow_nan=False) + "\n",
                         encoding="utf-8")
    (output / "PLAN_SHA256").write_text(
        plan["query_plan_sha256"] + "  shape_count_plan.json\n", encoding="ascii"
    )
    return plan


def main() -> None:
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    plan = prepare(args.output)
    print(json.dumps({
        "status": plan["status"],
        "network_used": plan["network_used"],
        "query_submitted": plan["query_submitted"],
        "candidate_lenses": plan["n_candidate_lenses"],
        "selected_lenses": plan["n_selected_lenses"],
        "random_centers": plan["n_random_centers"],
        "centers_per_lens": plan["centers_per_lens"],
        "batches": plan["n_query_batches"],
        "sql_bytes": sum(batch["sql_bytes"] for batch in plan["batches"]),
        "estimated_result_rows": plan["estimated_result_rows"],
        "estimated_result_bytes_upper_bound": plan["estimated_result_bytes_upper_bound"],
        "query_plan_sha256": plan["query_plan_sha256"],
        "output": str(args.output.resolve()),
    }, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
