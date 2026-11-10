#!/usr/bin/env python3
"""Prepare, but never submit, a count-only S16A shape-coverage plan.

The plan tests whether the existing 1,043 S23B ICL targets overlap the public
S16A weak-lensing shape catalog. It writes aggregate-count SQL locally and
never contacts HSC CAS or returns object-level data.
"""
from __future__ import annotations

import hashlib
import json
import math
from pathlib import Path

from astropy.cosmology import FlatLambdaCDM


DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
TARGET_MANIFEST = DATA / "s23b_icl_target_image_plan_2026-11-07_01/target_manifest.json"
OUTPUT = DATA / "wl_s23b_s16a_shape_count_plan_2026-11-11_01"
EXPECTED_TARGETS = 1043
EXPECTED_TARGET_SHA256 = "e487582ba40cd3e58da67704242db4d7219e381fb08f4b33fed47ce4ff268b44"
RADIUS_MPC = 3.0
BATCH_SIZE = 25


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def sql_literal(value: str) -> str:
    value = str(value)
    if not value or any(ord(char) < 32 for char in value):
        raise ValueError("IDs must be non-empty printable strings")
    return "'" + value.replace("'", "''") + "'"


def sql_number(value: float) -> str:
    value = float(value)
    if not math.isfinite(value):
        raise ValueError("query coordinates and radius must be finite")
    return format(value, ".12g")


def count_query(center: dict) -> str:
    return (
        f"SELECT {sql_literal(center['lens_id'])} AS lens_id, "
        "COUNT(DISTINCT b.object_id) AS n_shape_sources\n"
        "FROM s16a_wide.meas2 AS a\n"
        "INNER JOIN s16a_wide.weaklensing_hsm_regauss AS b USING (object_id)\n"
        "INNER JOIN s16a_wide.meas AS c USING (object_id)\n"
        f"WHERE coneSearch(c.icoord, {sql_number(center['ra_deg'])}, "
        f"{sql_number(center['dec_deg'])}, {sql_number(center['radius_arcsec'])})"
    )


def canonical_sha(plan: dict) -> str:
    payload = dict(plan)
    payload.pop("query_plan_sha256", None)
    encoded = json.dumps(payload, sort_keys=True, separators=(",", ":"),
                         ensure_ascii=False, allow_nan=False).encode("utf-8")
    return hashlib.sha256(encoded).hexdigest()


def prepare() -> dict:
    if OUTPUT.exists():
        raise FileExistsError(f"refusing to overwrite S16A count plan: {OUTPUT}")
    if not TARGET_MANIFEST.is_file() or sha256(TARGET_MANIFEST) != EXPECTED_TARGET_SHA256:
        raise ValueError("pinned 1,043-target S23B manifest is missing or changed")
    document = json.loads(TARGET_MANIFEST.read_text(encoding="utf-8"))
    rows = document.get("clusters", [])
    if (document.get("schema") != "hsc-bcg-diffuse-sample-v1"
            or document.get("science_measurement") is not False
            or len(rows) != EXPECTED_TARGETS):
        raise ValueError("unexpected S23B target manifest or scientific status")
    ids = [str(row.get("id", "")) for row in rows]
    if (any(not value.startswith("S23B_") for value in ids)
            or len(set(ids)) != EXPECTED_TARGETS
            or any(row.get("sample_role") != "target" for row in rows)):
        raise ValueError("target manifest IDs must be unique S23B cluster centers")

    h0 = float(document["h0"])
    omega_m = float(document["omega_m"])
    if not math.isfinite(h0) or h0 <= 0 or not math.isfinite(omega_m) or not 0 < omega_m < 1:
        raise ValueError("invalid cosmology in the pinned target manifest")
    cosmology = FlatLambdaCDM(H0=h0, Om0=omega_m, Tcmb0=0)
    centers = []
    for row in rows:
        ra = float(row["ra_deg"])
        dec = float(row["dec_deg"])
        redshift = float(row["z_cl"])
        if not (math.isfinite(ra) and 0 <= ra < 360
                and math.isfinite(dec) and -90 <= dec <= 90
                and math.isfinite(redshift) and redshift > 0):
            raise ValueError(f"invalid S23B center or redshift: {row.get('id')}")
        da_mpc = float(cosmology.angular_diameter_distance(redshift).value)
        radius_arcsec = math.degrees(RADIUS_MPC / da_mpc) * 3600.0
        centers.append({
            "lens_id": str(row["id"]),
            "ra_deg": ra,
            "dec_deg": dec,
            "z_lens": redshift,
            "radius_mpc_proper": RADIUS_MPC,
            "radius_arcsec": radius_arcsec,
        })

    OUTPUT.mkdir(parents=True, exist_ok=False)
    batches = []
    for index, start in enumerate(range(0, len(centers), BATCH_SIZE)):
        chunk = centers[start:start + BATCH_SIZE]
        sql = "\nUNION ALL\n".join(count_query(center) for center in chunk) + ";\n"
        path = OUTPUT / f"s16a_shape_counts_{index:03d}.sql"
        path.write_text(sql, encoding="utf-8")
        batches.append({
            "name": path.name,
            "n_centers": len(chunk),
            "first_lens_id": chunk[0]["lens_id"],
            "last_lens_id": chunk[-1]["lens_id"],
            "sql_bytes": path.stat().st_size,
            "sha256": sha256(path),
        })

    plan = {
        "schema": "astrafek-s23b-s16a-wl-shape-count-plan-v1",
        "status": "prepared_not_submitted",
        "catalog_release_documentation": "HSC S16A shape catalog; CAS table namespace s16a_wide",
        "cas_release_version": "pdr1",
        "source_target_manifest": str(TARGET_MANIFEST),
        "source_target_manifest_sha256": sha256(TARGET_MANIFEST),
        "target_scope": "all 1,043 S23B CAMIRA centers in the existing HSC ICL target manifest",
        "n_centers": len(centers),
        "radius_mpc_proper": RADIUS_MPC,
        "cosmology": {"H0_km_s_mpc": h0, "Omega_m": omega_m},
        "n_query_batches": len(batches),
        "batch_size": BATCH_SIZE,
        "query_tables": [
            "s16a_wide.meas2", "s16a_wide.weaklensing_hsm_regauss", "s16a_wide.meas"],
        "spatial_function": "coneSearch(c.icoord, ra_deg, dec_deg, radius_arcsec)",
        "output_columns": ["lens_id", "n_shape_sources"],
        "query_semantics": "COUNT(DISTINCT S16A weaklensing_hsm_regauss object_id) in each 3 proper-Mpc cone before Astrafek source-z, radial, B-mode, and estimator cuts",
        "estimated_result_rows": len(centers),
        "estimated_result_scope": "1,043 aggregate integer shape-catalog counts; no source identifiers, positions, shapes, photo-z values, images, or PSF values",
        "network_used": False,
        "query_submitted": False,
        "source_rows_downloaded": False,
        "images_returned": False,
        "photoz_returned": False,
        "science_measurement": False,
        "execution_requires_explicit_user_approval": True,
        "catalog_notes": [
            "S16A shape catalog and S19A shape rows use different processing releases and calibrations; do not concatenate their rows before a separate cross-release and calibration audit",
            "S16A raw counts are a coverage ceiling only; they do not establish usable source counts after tomography, masking, and estimator cuts",
        ],
        "batches": batches,
        "centers": centers,
    }
    plan["query_plan_sha256"] = canonical_sha(plan)
    plan_path = OUTPUT / "shape_count_plan.json"
    plan_path.write_text(json.dumps(plan, ensure_ascii=False, indent=2,
                                    allow_nan=False) + "\n", encoding="utf-8")
    (OUTPUT / "PLAN_SHA256").write_text(
        plan["query_plan_sha256"] + "  shape_count_plan.json\n", encoding="ascii")
    return plan


def main() -> None:
    plan = prepare()
    print(json.dumps({
        "status": plan["status"],
        "network_used": plan["network_used"],
        "query_submitted": plan["query_submitted"],
        "n_centers": plan["n_centers"],
        "n_query_batches": plan["n_query_batches"],
        "sql_bytes": sum(batch["sql_bytes"] for batch in plan["batches"]),
        "estimated_result_rows": plan["estimated_result_rows"],
        "query_plan_sha256": plan["query_plan_sha256"],
        "output": str(OUTPUT),
    }, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
