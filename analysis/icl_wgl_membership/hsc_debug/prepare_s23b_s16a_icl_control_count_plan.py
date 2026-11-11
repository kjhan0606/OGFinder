#!/usr/bin/env python3
"""Prepare, but never submit, an S16A count-only query for saved ICL controls.

This plan screens the one existing ICL control coordinate associated with each
of 116 S16A high-count expansion candidates. It is a local preparation step;
it does not contact HSC CAS or return source-level data.
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
AUDIT_DIR = DATA / "s23b_s16a_shape_coverage_selection_audit_2026-11-11_04"
CONTROL_TABLE = AUDIT_DIR / "s16a_ge1000_outside_s19a_icl_control_centers.tsv"
AUDIT_SUMMARY = AUDIT_DIR / "summary.json"
OUTPUT = DATA / "wl_s23b_s16a_icl_control_count_plan_2026-11-11_01"

EXPECTED_AUDIT_SUMMARY_SHA256 = "c08602e1b45492a92389dabd074e5f854050e8e14c2cd55b4683b5ae6fc0b398"
EXPECTED_CONTROL_TABLE_SHA256 = "81c702166871e75570918f770ffad44702cebe9a29920c80d3d73d2263c21706"
RADIUS_MPC = 3.0
BATCH_SIZE = 25
EXPECTED_CANDIDATES = 126
EXPECTED_CONTROLS = 116


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def control_count_query(center: dict) -> str:
    return (
        f"SELECT {sql_literal(center['control_id'])} AS control_id, "
        "COUNT(DISTINCT b.object_id) AS n_shape_sources\n"
        "FROM s16a_wide.meas2 AS a\n"
        "INNER JOIN s16a_wide.weaklensing_hsm_regauss AS b USING (object_id)\n"
        "INNER JOIN s16a_wide.meas AS c USING (object_id)\n"
        f"WHERE coneSearch(c.icoord, {sql_number(center['ra_deg'])}, "
        f"{sql_number(center['dec_deg'])}, {sql_number(center['radius_arcsec'])})"
    )


def prepare() -> dict:
    if OUTPUT.exists():
        raise FileExistsError(f"refusing to overwrite control-count plan: {OUTPUT}")
    if sha256(CONTROL_TABLE) != EXPECTED_CONTROL_TABLE_SHA256:
        raise ValueError("pinned candidate ICL-control table is missing or changed")
    if sha256(AUDIT_SUMMARY) != EXPECTED_AUDIT_SUMMARY_SHA256:
        raise ValueError("pinned local selection-audit summary is missing or changed")
    if sha256(TARGET_MANIFEST) != EXPECTED_TARGET_SHA256:
        raise ValueError("pinned S23B target manifest is missing or changed")

    audit = json.loads(AUDIT_SUMMARY.read_text(encoding="utf-8"))
    candidate_scope = audit.get("expansion_candidate_icl_control_pairs", {})
    if (audit.get("status") != "complete_local_descriptive_audit"
            or audit.get("science_measurement") is not False
            or audit.get("network_used") is not False
            or audit.get("cas_query_submitted") is not False
            or audit.get("source_rows_downloaded") is not False
            or candidate_scope.get("n_candidate_targets") != EXPECTED_CANDIDATES
            or candidate_scope.get("n_targets_with_saved_icl_control_pair") != EXPECTED_CONTROLS
            or candidate_scope.get("n_unique_control_centers") != EXPECTED_CONTROLS):
        raise ValueError("local audit does not match the reviewed 126-target/116-control scope")

    manifest = json.loads(TARGET_MANIFEST.read_text(encoding="utf-8"))
    if (manifest.get("schema") != "hsc-bcg-diffuse-sample-v1"
            or manifest.get("science_measurement") is not False):
        raise ValueError("unexpected pinned S23B target-manifest schema or status")
    target_map = {str(row["id"]): row for row in manifest.get("clusters", [])}
    if len(target_map) != 1043:
        raise ValueError("pinned target manifest does not contain 1,043 unique targets")

    cosmology = FlatLambdaCDM(
        H0=float(manifest["h0"]), Om0=float(manifest["omega_m"]), Tcmb0=0
    )
    with CONTROL_TABLE.open(encoding="utf-8-sig", newline="") as stream:
        controls = list(csv.DictReader(stream, delimiter="\t"))
    required = [
        "target_id", "control_id", "icl_filter", "target_ra_deg", "target_dec_deg",
        "control_ra_deg", "control_dec_deg", "z_cl", "richness",
        "s16a_target_count_3pmpc", "delta_ebv_sfd_target_minus_control",
        "control_s16a_count_status",
    ]
    if not controls or list(controls[0]) != required or len(controls) != EXPECTED_CONTROLS:
        raise ValueError("candidate ICL-control table columns or row count changed")

    centers = []
    target_ids, control_ids = set(), set()
    for row in controls:
        target_id, control_id = row["target_id"], row["control_id"]
        if (target_id not in target_map or target_id in target_ids
                or not control_id.startswith("CTRL_S23B_") or control_id in control_ids):
            raise ValueError("candidate target/control IDs are missing, repeated, or unrecognized")
        if row["control_s16a_count_status"] != "not_queried":
            raise ValueError("control centers already have a count status; refusing a duplicate plan")
        target = target_map[target_id]
        try:
            z = float(row["z_cl"])
            ra, dec = float(row["control_ra_deg"]), float(row["control_dec_deg"])
            target_count = int(row["s16a_target_count_3pmpc"])
            delta_ebv = float(row["delta_ebv_sfd_target_minus_control"])
        except (TypeError, ValueError) as error:
            raise ValueError(f"invalid coordinate or count for {control_id}") from error
        if (target_count < 1000
                or not math.isfinite(z) or z <= 0
                or not math.isfinite(ra) or not 0 <= ra < 360
                or not math.isfinite(dec) or not -90 <= dec <= 90
                or not math.isfinite(delta_ebv)
                or z != float(target["z_cl"])
                or float(row["target_ra_deg"]) != float(target["ra_deg"])
                or float(row["target_dec_deg"]) != float(target["dec_deg"])
                or float(row["richness"]) != float(target["richness"])):
            raise ValueError(f"invalid or mismatched matched-center metadata for {control_id}")
        da_mpc = float(cosmology.angular_diameter_distance(z).value)
        radius_arcsec = math.degrees(RADIUS_MPC / da_mpc) * 3600.0
        centers.append({
            "control_id": control_id,
            "matched_lens_id": target_id,
            "icl_filter": row["icl_filter"],
            "ra_deg": ra,
            "dec_deg": dec,
            "z_lens": z,
            "radius_mpc_proper": RADIUS_MPC,
            "radius_arcsec": radius_arcsec,
            "matched_target_s16a_raw_count": target_count,
            "delta_ebv_sfd_target_minus_control": delta_ebv,
        })
        target_ids.add(target_id)
        control_ids.add(control_id)

    if len(target_ids) != EXPECTED_CONTROLS or len(control_ids) != EXPECTED_CONTROLS:
        raise ValueError("control plan does not contain 116 unique matched target/control IDs")

    batches = []
    all_sql = []
    for index, start in enumerate(range(0, len(centers), BATCH_SIZE)):
        chunk = centers[start:start + BATCH_SIZE]
        sql = "\nUNION ALL\n".join(control_count_query(center) for center in chunk) + ";\n"
        name = f"s16a_icl_control_counts_{index:03d}.sql"
        all_sql.append((name, sql, chunk))
        batches.append({
            "name": name,
            "n_centers": len(chunk),
            "first_control_id": chunk[0]["control_id"],
            "last_control_id": chunk[-1]["control_id"],
            "sql_bytes": len(sql.encode("utf-8")),
            "sha256": hashlib.sha256(sql.encode("utf-8")).hexdigest(),
        })

    header_bytes = len(b"control_id\tn_shape_sources\n")
    max_result_bytes = header_bytes + sum(
        len(center["control_id"].encode("utf-8")) + 1 + 19 + 1 for center in centers
    )
    plan = {
        "schema": "astrafek-s23b-s16a-icl-control-count-plan-v1",
        "status": "prepared_not_submitted",
        "catalog_release_documentation": "HSC S16A shape catalog; CAS namespace s16a_wide",
        "cas_release_version": "pdr1",
        "candidate_control_table": str(CONTROL_TABLE),
        "candidate_control_table_sha256": sha256(CONTROL_TABLE),
        "selection_audit_summary": str(AUDIT_SUMMARY),
        "selection_audit_summary_sha256": sha256(AUDIT_SUMMARY),
        "source_target_manifest": str(TARGET_MANIFEST),
        "source_target_manifest_sha256": sha256(TARGET_MANIFEST),
        "target_scope": "one previously saved ICL control center for each S16A raw-count >= 1,000 target outside the S19A 427-center count audit, where a saved ICL pair exists",
        "n_candidate_targets": EXPECTED_CANDIDATES,
        "n_controls": len(centers),
        "radius_mpc_proper": RADIUS_MPC,
        "cosmology": {"H0_km_s_mpc": float(manifest["h0"]), "Omega_m": float(manifest["omega_m"])},
        "n_query_batches": len(batches),
        "batch_size": BATCH_SIZE,
        "query_tables": [
            "s16a_wide.meas2", "s16a_wide.weaklensing_hsm_regauss", "s16a_wide.meas"],
        "spatial_function": "coneSearch(c.icoord, ra_deg, dec_deg, radius_arcsec)",
        "output_columns": ["control_id", "n_shape_sources"],
        "query_semantics": "COUNT(DISTINCT S16A weaklensing_hsm_regauss object_id) in each matched ICL-control 3 proper-Mpc cone, using its target cluster redshift, before source-z, shape-quality, B-mode, mask, radial, and estimator cuts",
        "estimated_result_rows": len(centers),
        "estimated_result_bytes_upper_bound": max_result_bytes,
        "estimated_result_scope": "116 control IDs and integer source counts only; no source identifiers, coordinates, shapes, photo-z values, images, or PSF values",
        "network_used": False,
        "query_submitted": False,
        "source_rows_downloaded": False,
        "images_returned": False,
        "photoz_returned": False,
        "science_measurement": False,
        "execution_requires_explicit_user_approval": True,
        "control_interpretation": "This one-to-one ICL-image control screen does not replace the five-random-centers-per-lens weak-lensing control design and is not sufficient for a lensing stack.",
        "later_gates": [
            "compare target/control raw catalog-count coverage without treating it as usable source support",
            "design and audit the five-random-centers-per-lens S16A weak-lensing control sample",
            "apply source-z, shape-quality, mask, radial, and estimator support cuts",
            "prepare an exact S16A source-row manifest and obtain separate approval before any source-row transfer",
            "keep S16A and S19A source rows separate pending release/calibration validation",
        ],
        "batches": batches,
        "centers": centers,
    }
    plan["query_plan_sha256"] = canonical_sha(plan)

    OUTPUT.mkdir(parents=True, exist_ok=False)
    for name, sql, _ in all_sql:
        (OUTPUT / name).write_text(sql, encoding="utf-8")
    plan_path = OUTPUT / "shape_count_plan.json"
    plan_path.write_text(json.dumps(plan, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    (OUTPUT / "PLAN_SHA256").write_text(
        plan["query_plan_sha256"] + "  shape_count_plan.json\n", encoding="ascii")
    return plan


def main() -> None:
    plan = prepare()
    print(json.dumps({
        "status": plan["status"],
        "network_used": plan["network_used"],
        "query_submitted": plan["query_submitted"],
        "n_candidate_targets": plan["n_candidate_targets"],
        "n_controls": plan["n_controls"],
        "n_query_batches": plan["n_query_batches"],
        "sql_bytes": sum(batch["sql_bytes"] for batch in plan["batches"]),
        "estimated_result_rows": plan["estimated_result_rows"],
        "estimated_result_bytes_upper_bound": plan["estimated_result_bytes_upper_bound"],
        "query_plan_sha256": plan["query_plan_sha256"],
        "output": str(OUTPUT),
    }, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
