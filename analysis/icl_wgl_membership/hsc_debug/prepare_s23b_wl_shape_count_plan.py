#!/usr/bin/env python3
"""Prepare, but never submit, S19A count-only queries for S23B profile targets."""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

from astropy.cosmology import FlatLambdaCDM


DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
TARGET_MANIFEST = DATA / "s23b_icl_target_image_plan_2026-11-07_01/target_manifest.json"
PROFILE_INPUTS = DATA / "s23b_camira_member_outer_icl_relation_2026-11-07_01/cluster_outer_relation_inputs.tsv"
MEMBER_FILE = DATA / "s23b_camira_member_catalog_2026-11-07_01/camira_s23b_wide_sm_member.dat.gz"
OUTPUT = DATA / "wl_s23b_shape_count_plan_2026-11-07_01"
BATCH_SIZE = 25
RADIUS_MPC = 3.0
H0 = 70.0
OMEGA_M = 0.3


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def sql_literal(value: str) -> str:
    return "'" + str(value).replace("'", "''") + "'"


def sql_number(value: float) -> str:
    value = float(value)
    if not math.isfinite(value):
        raise ValueError("query coordinates and radii must be finite")
    return format(value, ".12g")


def main() -> None:
    if OUTPUT.exists():
        raise FileExistsError(f"refusing to overwrite query plan: {OUTPUT}")
    for path in (TARGET_MANIFEST, PROFILE_INPUTS, MEMBER_FILE):
        if not path.is_file():
            raise FileNotFoundError(path)

    target_doc = json.loads(TARGET_MANIFEST.read_text(encoding="utf-8"))
    target_rows = target_doc.get("clusters", [])
    target_by_catalog_id = {str(row["catalog_id"]): row for row in target_rows
                            if row.get("sample_role") == "target"}
    catalog_ids: set[str] = set()
    with PROFILE_INPUTS.open(encoding="utf-8", newline="") as stream:
        header = stream.readline().rstrip("\n").split("\t")
        if not {"target_id", "filter"}.issubset(header):
            raise ValueError("member-profile input table has an unexpected schema")
        for line in stream:
            if not line.strip():
                continue
            values = dict(zip(header, line.rstrip("\n").split("\t")))
            catalog_id = str(values["target_id"]).removeprefix("S23B_")
            if catalog_id not in target_by_catalog_id:
                raise ValueError(f"profile target is absent from target manifest: {catalog_id}")
            catalog_ids.add(catalog_id)
    if len(catalog_ids) != 427:
        raise ValueError(f"expected 427 unique profile targets; found {len(catalog_ids)}")

    cosmology = FlatLambdaCDM(H0=H0, Om0=OMEGA_M, Tcmb0=0)
    centers = []
    for catalog_id in sorted(catalog_ids, key=lambda value: (int(value), value)):
        row = target_by_catalog_id[catalog_id]
        z = float(row["z_cl"])
        da_mpc = float(cosmology.angular_diameter_distance(z).value)
        radius_arcsec = math.degrees(RADIUS_MPC / da_mpc) * 3600.0
        centers.append({
            "lens_id": f"S23B_{catalog_id}",
            "catalog_id": catalog_id,
            "ra_deg": float(row["ra_deg"]),
            "dec_deg": float(row["dec_deg"]),
            "z_lens": z,
            "radius_mpc_proper": RADIUS_MPC,
            "radius_arcsec": radius_arcsec,
        })

    OUTPUT.mkdir(parents=True, exist_ok=False)
    batches = []
    for index, start in enumerate(range(0, len(centers), BATCH_SIZE)):
        chunk = centers[start:start + BATCH_SIZE]
        statements = []
        for center in chunk:
            statements.append(
                f"SELECT {sql_literal(center['lens_id'])} AS lens_id, "
                "COUNT(DISTINCT b.object_id) AS n_shape_sources\n"
                "FROM s19a_wide.meas2 a\n"
                "INNER JOIN s19a_wide.weaklensing_hsm_regauss b USING (object_id)\n"
                "INNER JOIN s19a_wide.meas c USING (object_id)\n"
                f"WHERE coneSearch(c.i_coord, {sql_number(center['ra_deg'])}, "
                f"{sql_number(center['dec_deg'])}, {sql_number(center['radius_arcsec'])})"
            )
        path = OUTPUT / f"shape_counts_{index:03d}.sql"
        path.write_text("\nUNION ALL\n".join(statements) + ";\n", encoding="utf-8")
        batches.append({
            "name": path.name,
            "n_lenses": len(chunk),
            "first_lens_id": chunk[0]["lens_id"],
            "last_lens_id": chunk[-1]["lens_id"],
            "sql_bytes": path.stat().st_size,
            "sha256": sha256(path),
        })

    plan = {
        "schema": "astrafek-s23b-s19a-wl-shape-count-plan-v1",
        "status": "prepared_not_submitted",
        "release": "HSC PDR3 CAS; S19A Wide weak-lensing shape catalog",
        "source_target_manifest": str(TARGET_MANIFEST),
        "source_target_manifest_sha256": sha256(TARGET_MANIFEST),
        "source_member_profile_inputs": str(PROFILE_INPUTS),
        "source_member_profile_inputs_sha256": sha256(PROFILE_INPUTS),
        "approved_member_catalog_file": str(MEMBER_FILE),
        "approved_member_catalog_sha256": sha256(MEMBER_FILE),
        "target_scope": "427 unique S23B CAMIRA clusters that appear in the saved member outer-profile relation; lens_id is explicitly namespaced as S23B_<catalog_id>",
        "n_lenses": len(centers),
        "radius_mpc_proper": RADIUS_MPC,
        "cosmology": {"H0_km_s_mpc": H0, "Omega_m": OMEGA_M},
        "n_query_batches": len(batches),
        "batch_size": BATCH_SIZE,
        "query_tables": ["s19a_wide.meas2", "s19a_wide.weaklensing_hsm_regauss", "s19a_wide.meas"],
        "output_columns": ["lens_id", "n_shape_sources"],
        "query_semantics": "COUNT(DISTINCT object_id) inside a 3 proper-Mpc lens cone, before shape-quality, source-redshift, tomography, B-mode, or radial cuts",
        "estimated_result_rows": len(centers),
        "estimated_result_scope": "427 aggregate integer counts; no source coordinates, shapes, photo-z values, images, PSF, or variance data",
        "network_used": False,
        "query_submitted": False,
        "source_rows_downloaded": False,
        "images_returned": False,
        "photoz_returned": False,
        "science_measurement": False,
        "execution_requires_explicit_user_approval": True,
        "batches": batches,
        "centers": centers,
    }
    canonical = json.dumps(plan, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode("utf-8")
    plan["query_plan_sha256"] = hashlib.sha256(canonical).hexdigest()
    plan_path = OUTPUT / "shape_count_plan.json"
    plan_path.write_text(json.dumps(plan, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    (OUTPUT / "PLAN_SHA256").write_text(plan["query_plan_sha256"] + "  shape_count_plan.json\n", encoding="ascii")
    print(json.dumps({
        "status": plan["status"], "query_submitted": False,
        "network_used": False, "lenses": len(centers), "batches": len(batches),
        "sql_bytes": sum(batch["sql_bytes"] for batch in batches),
        "result_rows": plan["estimated_result_rows"],
        "query_plan_sha256": plan["query_plan_sha256"],
        "output": str(OUTPUT),
    }, ensure_ascii=False))


if __name__ == "__main__":
    main()
