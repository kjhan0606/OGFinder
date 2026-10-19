#!/usr/bin/env python3
"""Prepare a count-only HSC S19A support plan for selected random lenses.

This writes SQL and provenance locally. It never contacts HSC CAS or retrieves
source rows. Submission is a separate user-approved operation.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

from astropy.cosmology import FlatLambdaCDM


OUTPUT_COLUMNS = ["random_center_id", "matched_lens_id", "n_shape_sources"]
DEFAULT_BATCH_SIZE = 25


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def sql_literal(value):
    value = str(value)
    if not value or any(ord(char) < 32 for char in value):
        raise ValueError("IDs must be non-empty printable strings")
    return "'" + value.replace("'", "''") + "'"


def sql_number(value):
    number = float(value)
    if not math.isfinite(number):
        raise ValueError("center coordinates and redshifts must be finite")
    return format(number, ".12g")


def shape_count_query(center):
    center_id = sql_literal(center["random_center_id"])
    matched_id = sql_literal(center["matched_lens_id"])
    ra = sql_number(center["ra_deg"])
    dec = sql_number(center["dec_deg"])
    radius = sql_number(center["radius_arcsec"])
    return (
        f"SELECT {center_id} AS random_center_id, {matched_id} AS matched_lens_id, "
        "COUNT(DISTINCT b.object_id) AS n_shape_sources\n"
        "FROM s19a_wide.meas2 a\n"
        "INNER JOIN s19a_wide.weaklensing_hsm_regauss b USING (object_id)\n"
        "INNER JOIN s19a_wide.meas c USING (object_id)\n"
        f"WHERE coneSearch(c.i_coord, {ra}, {dec}, {radius})"
    )


def normalize_centers(document, *, radius_mpc=3.0, batch_size=DEFAULT_BATCH_SIZE):
    if document.get("schema") != "astrafek-wl-random-centers-v1":
        raise ValueError("unexpected random-center manifest schema")
    if document.get("status") != "selected_for_shape_count_preflight":
        raise ValueError("centers must be explicitly selected before shape-count planning")
    if document.get("science_measurement") is not False:
        raise ValueError("random-center manifest cannot claim a science measurement")
    h0, omega_m = float(document["h0_km_s_mpc"]), float(document["omega_m"])
    if not math.isfinite(h0) or h0 <= 0 or not math.isfinite(omega_m) or not 0 < omega_m < 1:
        raise ValueError("invalid flat-LambdaCDM parameters")
    centers = document.get("centers")
    if (not isinstance(centers, list) or not centers or batch_size < 1
            or not math.isfinite(float(radius_mpc)) or radius_mpc <= 0):
        raise ValueError("centers, batch size, and physical radius must be positive and finite")
    ids = [str(row.get("random_center_id", "")) for row in centers]
    if any(not value for value in ids) or len(set(ids)) != len(ids):
        raise ValueError("random_center_id values must be non-empty and unique")
    lens_ids = {str(row.get("matched_lens_id", "")) for row in centers}
    if "" in lens_ids:
        raise ValueError("every random center must reference its matched lens")
    cosmo = FlatLambdaCDM(H0=h0, Om0=omega_m, Tcmb0=0)
    normalized = []
    matched_redshift = {}
    for row in centers:
        ra, dec, redshift = (float(row[name]) for name in ("ra_deg", "dec_deg", "z_lens"))
        if not math.isfinite(ra) or not 0 <= ra < 360 or not math.isfinite(dec) or not -90 <= dec <= 90:
            raise ValueError("invalid random-center sky coordinate")
        if not math.isfinite(redshift) or redshift <= 0:
            raise ValueError("matched lens redshift must be finite and positive")
        matched_id = str(row["matched_lens_id"])
        if matched_id in matched_redshift and not math.isclose(
                redshift, matched_redshift[matched_id], rel_tol=0, abs_tol=1e-12):
            raise ValueError("random centers matched to one lens must use the same lens redshift")
        matched_redshift[matched_id] = redshift
        da = float(cosmo.angular_diameter_distance(redshift).value)
        radius_arcsec = math.degrees(radius_mpc / da) * 3600.0
        normalized.append({
            "random_center_id": str(row["random_center_id"]),
            "matched_lens_id": str(row["matched_lens_id"]),
            "ra_deg": ra, "dec_deg": dec, "z_lens": redshift,
            "radius_mpc_proper": float(radius_mpc),
            "radius_arcsec": radius_arcsec,
        })
    return normalized, h0, omega_m, len(lens_ids)


def prepare(centers_path: Path, output: Path, *, radius_mpc=3.0, batch_size=DEFAULT_BATCH_SIZE):
    centers_path, output = Path(centers_path).resolve(), Path(output).resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite existing plan: {output}")
    document = json.loads(centers_path.read_text(encoding="utf-8"))
    centers, h0, omega_m, n_matched_lenses = normalize_centers(
        document, radius_mpc=radius_mpc, batch_size=batch_size,
    )
    output.mkdir(parents=True, exist_ok=False)
    batches = []
    for index, start in enumerate(range(0, len(centers), batch_size)):
        chunk = centers[start:start + batch_size]
        sql = "\nUNION ALL\n".join(shape_count_query(center) for center in chunk) + "\n"
        path = output / f"random_shape_counts_{index:03d}.sql"
        path.write_text(sql, encoding="utf-8")
        batches.append({
            "path": str(path), "n_random_centers": len(chunk),
            "first_random_center_id": chunk[0]["random_center_id"],
            "last_random_center_id": chunk[-1]["random_center_id"],
            "sha256": sha(path),
        })
    manifest = {
        "schema": "hsc-s19a-random-shape-count-plan-v1",
        "status": "prepared_not_submitted",
        "source_centers": str(centers_path),
        "source_centers_sha256": sha(centers_path),
        "release": "pdr3",
        "query_tables": ["s19a_wide.meas2", "s19a_wide.weaklensing_hsm_regauss", "s19a_wide.meas"],
        "n_random_centers": len(centers), "n_matched_lenses": n_matched_lenses,
        "n_query_batches": len(batches), "batch_size": int(batch_size),
        "radius_mpc_proper": float(radius_mpc), "cosmology": {"h0_km_s_mpc": h0, "omega_m": omega_m},
        "output_columns": OUTPUT_COLUMNS,
        "query_semantics": "raw distinct S19A shape-object count inside each random-center proper-radius cone; before shape, redshift, B-mode, and radial cuts",
        "source_shape_rows_returned": False,
        "photoz_rows_returned": False,
        "images_returned": False,
        "network_used": False,
        "query_submitted": False,
        "science_measurement": False,
        "execution_requires_explicit_user_approval": True,
        "batches": batches,
        "later_gates": [
            "compare per-center counts with matched cluster-lens support",
            "choose a balanced random-center count per matched lens",
            "estimate row and transfer volume from count results",
            "prepare exact source-row manifest and obtain separate explicit approval before export",
        ],
    }
    (output / "random_shape_count_plan.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8",
    )
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--centers", type=Path, required=True,
                        help="selected random-center manifest; planning only")
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--radius-mpc", type=float, default=3.0)
    parser.add_argument("--batch-size", type=int, default=DEFAULT_BATCH_SIZE)
    args = parser.parse_args()
    manifest = prepare(args.centers, args.work, radius_mpc=args.radius_mpc,
                       batch_size=args.batch_size)
    print(json.dumps({
        "status": manifest["status"], "network_used": False,
        "query_submitted": False, "random_centers": manifest["n_random_centers"],
        "matched_lenses": manifest["n_matched_lenses"],
        "batches": manifest["n_query_batches"], "output_columns": manifest["output_columns"],
        "work": str(args.work.resolve()),
    }, ensure_ascii=False))


if __name__ == "__main__":
    main()
