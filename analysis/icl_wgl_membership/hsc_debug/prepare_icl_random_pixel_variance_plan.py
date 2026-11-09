#!/usr/bin/env python3
"""Prepare a bounded HSC PDR3 random-pixel variance export plan for ICL.

This runner prepares SQL only. Submission is a separate operation requiring
explicit approval of the exact generated plan SHA-256.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

from astropy.cosmology import FlatLambdaCDM


ROOT = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
TARGETS = ROOT / "sample_manifest.json"
CONTROLS = ROOT / "control_download_01/download_manifest.json"
DEFAULT_OUTPUT = ROOT / "icl_random_pixel_variance_plan_2026-11-06_01"
BATCH_SIZE = 100
THINNING = 0.2  # HSC adjust_density < 0.2 gives a nominal ~20 points/arcmin^2.
INNER_KPC = 200.0
OUTER_KPC = 300.0
COSMOLOGY = FlatLambdaCDM(H0=70.0, Om0=0.3, Tcmb0=2.725)
OUTPUT_COLUMNS = [
    "center_id", "sample_role", "target_cluster_id", "center_ra_deg",
    "center_dec_deg", "cluster_z_for_scale", "random_object_id", "ra", "dec",
    "skymap_id", "tract", "patch", "patch_s", "adjust_density",
    "i_pix_variance", "i_sky_mean", "i_sky_std", "i_inputcount_value",
    "i_pixelflags_interpolatedcenter", "i_pixelflags_saturatedcenter",
    "i_pixelflags_crcenter", "i_pixelflags_suspectcenter",
    "i_pixelflags_bright_objectcenter", "i_pixelflags_edge",
    "i_sdssshape_psf_shape11", "i_sdssshape_psf_shape22",
    "i_sdssshape_psf_shape12",
]
NO_DATA_FIELDS = ["coadd image pixels", "shape catalog rows", "photo-z/PDF", "galaxy object catalog rows", "credentials"]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def quote(value: str) -> str:
    return "'" + value.replace("'", "''") + "'"


def centers_from_manifests() -> list[dict]:
    target_doc = json.loads(TARGETS.read_text(encoding="utf-8"))
    control_doc = json.loads(CONTROLS.read_text(encoding="utf-8"))
    targets = {str(row["id"]): row for row in target_doc["clusters"]}
    if len(targets) != 800:
        raise ValueError("expected exactly 800 target clusters")
    controls_by_target = {}
    for row in control_doc["clusters"]:
        target_id = str(row["target_cluster_id"])
        if target_id not in targets or target_id in controls_by_target:
            raise ValueError("control manifest has an unknown or duplicate target mapping")
        controls_by_target[target_id] = row
    if len(controls_by_target) != 798:
        raise ValueError("expected 798 matched control centers")

    centers = []
    for target_id, target in targets.items():
        centers.append({
            "center_id": f"target:{target_id}", "role": "target",
            "target_cluster_id": target_id, "ra_deg": float(target["ra_deg"]),
            "dec_deg": float(target["dec_deg"]), "z": float(target["z"]),
        })
        control = controls_by_target.get(target_id)
        if control is not None:
            centers.append({
                "center_id": f"control:{control['id']}", "role": "control",
                "target_cluster_id": target_id, "ra_deg": float(control["ra_deg"]),
                "dec_deg": float(control["dec_deg"]), "z": float(target["z"]),
            })
    if len(centers) != 1598:
        raise ValueError(f"expected 1,598 target/control centers; found {len(centers)}")
    for row in centers:
        kpc_per_arcsec = COSMOLOGY.kpc_proper_per_arcmin(row["z"]).to_value("kpc/arcsec")
        row["inner_arcsec"] = INNER_KPC / kpc_per_arcsec
        row["outer_arcsec"] = OUTER_KPC / kpc_per_arcsec
    return centers


def sql_for_batch(rows: list[dict]) -> str:
    values = ",\n        ".join(
        "(" + ", ".join((quote(row["center_id"]), quote(row["role"]),
            quote(row["target_cluster_id"]), f"{row['ra_deg']:.10f}",
            f"{row['dec_deg']:.10f}", f"{row['z']:.8f}",
            f"{row['inner_arcsec']:.8f}", f"{row['outer_arcsec']:.8f}")) + ")"
        for row in rows
    )
    return f"""-- Candidate HSC PDR3 random-point variance sample for 200--300 kpc ICL annuli.
-- Do not submit without approval of this exact plan SHA-256.
WITH centers(center_id, sample_role, target_cluster_id, ra0, dec0, cluster_z,
             theta_inner_arcsec, theta_outer_arcsec) AS (
    VALUES
        {values}
)
SELECT c.center_id, c.sample_role, c.target_cluster_id,
       c.ra0 AS center_ra_deg, c.dec0 AS center_dec_deg,
       c.cluster_z AS cluster_z_for_scale,
       r.object_id AS random_object_id, r.ra, r.dec,
       r.skymap_id, r.tract, r.patch, r.patch_s, r.adjust_density,
       r.i_pix_variance, r.i_sky_mean, r.i_sky_std, r.i_inputcount_value,
       r.i_pixelflags_interpolatedcenter, r.i_pixelflags_saturatedcenter,
       r.i_pixelflags_crcenter, r.i_pixelflags_suspectcenter,
       r.i_pixelflags_bright_objectcenter, r.i_pixelflags_edge,
       r.i_sdssshape_psf_shape11, r.i_sdssshape_psf_shape22,
       r.i_sdssshape_psf_shape12
FROM centers AS c
JOIN pdr3_wide.random AS r
  ON conesearch(r.coord, c.ra0, c.dec0, c.theta_outer_arcsec)
WHERE r.adjust_density < {THINNING:.2f}
  AND r.i_pixelflags_bright_objectcenter IS FALSE
  AND r.i_pixelflags_edge IS FALSE
  AND 2.0 * asin(sqrt(
        power(sin(radians(r.dec - c.dec0) / 2.0), 2.0)
        + cos(radians(c.dec0)) * cos(radians(r.dec))
          * power(sin(radians(r.ra - c.ra0) / 2.0), 2.0)
      )) * 206264.80624709636 BETWEEN c.theta_inner_arcsec AND c.theta_outer_arcsec
ORDER BY c.center_id, r.object_id;
"""


def build(output_dir: Path) -> dict:
    output_dir = Path(output_dir).resolve()
    if output_dir.exists():
        raise FileExistsError(f"refusing to overwrite {output_dir}")
    centers = centers_from_manifests()
    batches = []
    expected_rows = 0.0
    for start in range(0, len(centers), BATCH_SIZE):
        index = start // BATCH_SIZE
        group = centers[start:start + BATCH_SIZE]
        path = output_dir / f"icl_random_variance_{index:03d}.sql"
        sql = sql_for_batch(group)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(sql, encoding="utf-8")
        # Nominal HSC random density is 100/arcmin^2; adjust_density < 0.2
        # retains approximately 20/arcmin^2. This is only a size estimate.
        area_arcmin2 = [math.pi * (r["outer_arcsec"]**2 - r["inner_arcsec"]**2) / 3600.0
                        for r in group]
        batch_expected = THINNING * 100.0 * sum(area_arcmin2)
        expected_rows += batch_expected
        batches.append({
            "name": path.stem, "path": str(path), "sha256": sha256(path),
            "centers": len(group), "first_center_id": group[0]["center_id"],
            "last_center_id": group[-1]["center_id"], "sql_bytes": path.stat().st_size,
            "center_ids": [row["center_id"] for row in group],
            "expected_rows_from_nominal_random_density": round(batch_expected),
        })
    if len(batches) != 16:
        raise ValueError("expected sixteen 100-center SQL batches")
    input_paths = (TARGETS, CONTROLS)
    plan = {
        "schema": "astrafek-icl-random-pixel-variance-plan-v1",
        "status": "prepared_not_submitted",
        "release": "HSC PDR3 Wide",
        "table": "pdr3_wide.random",
        "purpose": "sample HSC random-point variance, sky, input-count, local PSF and mask metadata in paired 200--300 proper-kpc target/control ICL annuli",
        "n_target_centers": 800, "n_control_centers": 798, "n_centers": len(centers),
        "n_batches": len(batches), "batch_size_centers": BATCH_SIZE,
        "annulus_proper_kpc": [INNER_KPC, OUTER_KPC],
        "cosmology": {"H0_km_s_Mpc": 70.0, "Om0": 0.3, "redshift": "paired target cluster z"},
        "random_density_thinning": {"predicate": "adjust_density < 0.2", "nominal_points_per_arcmin2": 20},
        "estimated_rows": round(expected_rows),
        "estimated_csv_size": f"approximately {round(expected_rows * 220 / 1_000_000)} MB at 220 bytes/row; actual yield depends on masks and HSC random density",
        "output_columns": OUTPUT_COLUMNS,
        "not_requested": NO_DATA_FIELDS,
        "images_returned": False, "shape_rows_returned": False,
        "photoz_rows_returned": False, "object_catalog_rows_returned": False,
        "network_used": False, "query_submitted": False,
        "credentials_saved": False, "science_measurement": False,
        "input_sha256": {str(path): sha256(path) for path in input_paths},
        "centers": centers,
        "batches": batches,
        "sql_aggregate_sha256": hashlib.sha256("\n".join(row["sha256"] for row in batches).encode()).hexdigest(),
        "execution_requires_explicit_user_approval": True,
        "limitations": [
            "The HSC random-table pix_variance is a random-point variance observable, not a complete per-pixel coadd weight map.",
            "The 200--300 kpc annulus is the ICL aperture; this design does not provide adequate random-point density for 50--75 kpc satellite noise correction.",
            "The nominal 100 random points/arcmin^2 density is approximate; returned row count is data-dependent.",
            "This export alone cannot establish sky subtraction, correlated-noise covariance, source-Poisson correction, or a science measurement.",
        ],
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    plan_path = output_dir / "plan.json"
    plan_path.write_text(json.dumps(plan, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return {"plan": str(plan_path), "plan_sha256": sha256(plan_path),
            "sql_aggregate_sha256": plan["sql_aggregate_sha256"],
            "batches": len(batches), "centers": len(centers),
            "estimated_rows": round(expected_rows), "sql_bytes": sum(r["sql_bytes"] for r in batches),
            "network_used": False, "query_submitted": False,
            "images": False, "shape_rows": False, "photoz": False,
            "science_measurement": False}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    print(json.dumps(build(args.output), ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
