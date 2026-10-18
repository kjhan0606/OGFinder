#!/usr/bin/env python3
"""Prepare count-only HSC S19A shape-catalog cone queries for the 800 lenses.

This writes SQL and a provenance manifest locally. It makes no network request
and does not retrieve shape rows or image data.
"""
import argparse
import hashlib
import json
import math
from pathlib import Path

from astropy.cosmology import FlatLambdaCDM


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--radius-mpc", type=float, default=3.0)
    parser.add_argument("--batch-size", type=int, default=25)
    args = parser.parse_args()
    source = args.manifest.resolve()
    work = args.work.resolve()
    if work.exists():
        parser.error(f"Refusing to overwrite existing output: {work}")
    if args.radius_mpc <= 0 or args.batch_size < 1:
        parser.error("radius and batch size must be positive")

    doc = json.loads(source.read_text(encoding="utf-8"))
    frames = doc.get("clusters", [])
    if len(frames) != 800 or len({str(row["id"]) for row in frames}) != 800:
        parser.error("Expected exactly 800 unique lens rows")
    cosmo = FlatLambdaCDM(H0=float(doc["h0"]), Om0=float(doc["omega_m"]), Tcmb0=0)
    requests = []
    for row in frames:
        da_mpc = float(cosmo.angular_diameter_distance(float(row["z"])).value)
        radius_arcsec = math.degrees(float(args.radius_mpc) / da_mpc) * 3600.0
        requests.append({
            "lens_id": str(row["id"]), "ra_deg": float(row["ra_deg"]),
            "dec_deg": float(row["dec_deg"]), "z_lens": float(row["z"]),
            "radius_mpc_proper": float(args.radius_mpc), "radius_arcsec": radius_arcsec,
        })

    work.mkdir(parents=True, exist_ok=False)
    files = []
    for start in range(0, len(requests), args.batch_size):
        batch = requests[start:start + args.batch_size]
        statements = []
        for item in batch:
            lens_id = item["lens_id"].replace("'", "''")
            ra, dec, radius = item["ra_deg"], item["dec_deg"], item["radius_arcsec"]
            statements.append(
                "SELECT '" + lens_id + "' AS lens_id, COUNT(DISTINCT b.object_id) AS n_shape_sources\n"
                "FROM s19a_wide.meas2 a\n"
                "INNER JOIN s19a_wide.weaklensing_hsm_regauss b USING (object_id)\n"
                "INNER JOIN s19a_wide.meas c USING (object_id)\n"
                f"WHERE coneSearch(c.i_coord, {ra:.8f}, {dec:.8f}, {radius:.4f})"
            )
        path = work / f"shape_counts_{start // args.batch_size:03d}.sql"
        path.write_text("\nUNION ALL\n".join(statements) + ";\n", encoding="utf-8")
        files.append({"path": str(path), "first_lens_id": batch[0]["lens_id"],
                      "last_lens_id": batch[-1]["lens_id"], "n_lenses": len(batch),
                      "sha256": sha256(path)})

    manifest = {
        "schema": "hsc-s19a-shape-count-preflight-v1",
        "status": "prepared_not_submitted",
        "source_manifest": str(source), "source_sha256": sha256(source),
        "n_lenses": len(requests), "radius_mpc_proper": float(args.radius_mpc),
        "radius_arcsec_range": [min(x["radius_arcsec"] for x in requests),
                                 max(x["radius_arcsec"] for x in requests)],
        "batch_size": args.batch_size, "n_sql_batches": len(files),
        "sql_files": files, "lenses": requests,
        "query_tables": ["s19a_wide.meas2", "s19a_wide.weaklensing_hsm_regauss", "s19a_wide.meas"],
        "query_function": "coneSearch(c.i_coord, ra_deg, dec_deg, radius_arcsec)",
        "query_output": "one count-only row per lens; no source rows are returned",
        "data_transfer_started": False, "credentials_requested": False,
        "science_measurement": False,
        "notes": [
            "The output counts raw S19A shape-catalog rows inside each lens cone, before source-redshift selection or lensing cuts.",
            "Repeated source object_id values across different lenses are expected and must be retained as distinct lens-source pairs in a later analysis.",
            "The 3 Mpc radius is a planning aperture for support counts, not the final science radial cut.",
            "S19A tables are documented for the PDR3 CAS release; the query has not been submitted or schema-validated against CAS.",
        ],
    }
    out = work / "preflight_manifest.json"
    out.write_text(json.dumps(manifest, ensure_ascii=False, indent=2, allow_nan=False) + "\n",
                   encoding="utf-8")
    print(json.dumps({"status": manifest["status"], "n_lenses": len(requests),
                      "n_sql_batches": len(files), "output": str(work)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
