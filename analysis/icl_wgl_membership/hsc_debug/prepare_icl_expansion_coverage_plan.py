#!/usr/bin/env python3
"""Prepare offline PDR3 coadd-footprint metadata queries for S20A CAMIRA candidates.

This script selects the archived 0.20 <= z_cl < 0.70 cluster-redshift pool,
verifies its local CAMIRA member positions, and writes HSC-I mosaic metadata
SQL only. It performs no network request and downloads no images or catalogs.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
REPOSITORY_ROOT = HERE.parents[2]
if str(REPOSITORY_ROOT) not in sys.path:
    sys.path.insert(0, str(REPOSITORY_ROOT))
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from audit_camira_sample_expansion import (
    DEFAULT_CATALOG, DEFAULT_MANIFEST, MIN_SEP_DEG,
    angular_sep_deg, greedy_additions, read_catalog, sky_cell,
)
from ogfkit.wlcluster import ARCSEC, Cosmology


DATA = Path("/gpfs/kjhan/Astrafek/data")
DEFAULT_MEMBER_CATALOG = DATA / "catalog/camira_s20a_wide_sm_member.dat.gz"
DEFAULT_EXPANSION_AUDIT = DATA / "sample_800_2026-10-13/camira_expansion_pool_2026-10-19_02/summary.json"
DEFAULT_OUTPUT = DATA / "sample_800_2026-10-13/icl_expansion_coverage_plan_2026-10-21_01"
EXPECTED_EXPANSION_AUDIT_SHA256 = "27ada4545891397f0da0bd8ed4f410a32cda0dfeea197b088a3dcc204523c2ed"
EXPECTED_MEMBER_CATALOG_SHA256 = "bf0c1fec5e11c6378ea69e3a65e49e9b32beecc83e5e68000bfa997e4ef2059e"
EXPECTED_CATALOG_EXPANSIONS = 2596
EXPECTED_NEW_5DEG_CELLS = 6
BOUNDARY_KPC = 1000.0
BOUNDARY_BEARINGS_DEG = tuple(range(0, 360, 45))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def destination(ra_deg: float, dec_deg: float, distance_rad: float,
                bearing_deg: float) -> tuple[float, float]:
    """Great-circle destination; bearing is clockwise from celestial north."""
    lon1, lat1, bearing = map(math.radians, (ra_deg, dec_deg, bearing_deg))
    sin_lat2 = (math.sin(lat1) * math.cos(distance_rad)
                + math.cos(lat1) * math.sin(distance_rad) * math.cos(bearing))
    lat2 = math.asin(max(-1.0, min(1.0, sin_lat2)))
    lon2 = lon1 + math.atan2(
        math.sin(bearing) * math.sin(distance_rad) * math.cos(lat1),
        math.cos(distance_rad) - math.sin(lat1) * math.sin(lat2),
    )
    return math.degrees(lon2) % 360.0, math.degrees(lat2)


def make_probes(row: dict, cosmology: Cosmology) -> list[dict]:
    da_mpc = float(cosmology.angular_diameter(float(row["z"])))
    if not math.isfinite(da_mpc) or da_mpc <= 0:
        raise ValueError(f"invalid angular-diameter distance for CAMIRA ID {row['id']}")
    radius_rad = (BOUNDARY_KPC / 1000.0) / da_mpc
    points = [{"candidate_id": str(row["id"]), "probe_id": "center",
               "ra_deg": float(row["ra_deg"]), "dec_deg": float(row["dec_deg"]),
               "radius_kpc": 0.0}]
    for bearing in BOUNDARY_BEARINGS_DEG:
        ra, dec = destination(float(row["ra_deg"]), float(row["dec_deg"]),
                              radius_rad, bearing)
        points.append({"candidate_id": str(row["id"]),
                       "probe_id": f"r1Mpc_pa{bearing:03d}",
                       "ra_deg": ra, "dec_deg": dec,
                       "radius_kpc": BOUNDARY_KPC})
    return points


def select_candidates(catalog: list[dict], current: list[dict]) -> list[dict]:
    current_ids = {str(row["id"]) for row in current}
    if len(current_ids) != 800:
        raise ValueError("current sample must contain 800 unique CAMIRA IDs")
    if not current_ids.issubset({str(row["id"]) for row in catalog}):
        raise ValueError("current sample includes IDs absent from the local S20A catalog")
    eligible = [row for row in catalog
                if str(row["id"]) not in current_ids and 0.20 <= row["z"] < 0.70]
    selected = greedy_additions(current, eligible, MIN_SEP_DEG)
    if len(selected) != EXPECTED_CATALOG_EXPANSIONS:
        raise ValueError(f"expected {EXPECTED_CATALOG_EXPANSIONS} additions, got {len(selected)}")
    cells = {sky_cell(row) for row in current}
    new_cells = {sky_cell(row) for row in selected} - cells
    if len(new_cells) != EXPECTED_NEW_5DEG_CELLS:
        raise ValueError(f"expected {EXPECTED_NEW_5DEG_CELLS} new 5-degree cells, got {len(new_cells)}")
    return selected


def attach_members(candidates: list[dict], member_path: Path) -> dict[str, list[tuple[float, float]]]:
    by_cluster_key = {}
    rows_by_id = {str(row["id"]): row for row in candidates}
    for row in candidates:
        key = (round(float(row["ra_deg"]), 5), round(float(row["dec_deg"]), 5),
               round(float(row["richness"]), 3), round(float(row["z"]), 3))
        if key in by_cluster_key:
            raise ValueError("two candidates share the rounded member-catalog identity")
        by_cluster_key[key] = str(row["id"])
    members = {ident: [] for ident in rows_by_id}
    with gzip.open(member_path, "rt", encoding="utf-8") as stream:
        for line_no, line in enumerate(stream, 1):
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            values = [float(value) for value in line.split()]
            if len(values) != 8:
                raise ValueError(f"member row {line_no} has {len(values)} fields, expected 8")
            # Several CAMIRA releases can place distinct rows at the same
            # rounded center; include richness and z_cl to select the matching
            # cluster instead of treating coordinate coincidence as identity.
            ident = by_cluster_key.get((round(values[0], 5), round(values[1], 5),
                                        round(values[2], 3), round(values[3], 3)))
            if ident is None:
                continue
            cluster = rows_by_id[ident]
            # Member file columns 3–4 are N_mem then z_cl, not the order used
            # in the cluster table.
            if abs(values[2] - cluster["richness"]) > 0.001 or abs(values[3] - cluster["z"]) > 0.001:
                raise ValueError(f"CAMIRA member row does not match cluster {ident}")
            ra, dec = values[4], values[5]
            if not (math.isfinite(ra) and 0 <= ra < 360 and math.isfinite(dec) and -90 <= dec <= 90):
                raise ValueError(f"invalid member coordinate for CAMIRA ID {ident}")
            members[ident].append((ra, dec))
    missing = [ident for ident, values in members.items() if not values]
    if missing:
        raise ValueError(f"candidate clusters lack local member rows: {missing[:10]}")
    return members


def render_query(probes: list[dict]) -> str:
    values = ",\n".join(
        "    ('%s', '%s', %.9f, %.9f)" % (
            p["candidate_id"].replace("'", "''"), p["probe_id"].replace("'", "''"),
            p["ra_deg"], p["dec_deg"],
        ) for p in probes
    )
    return (
        "WITH probes(candidate_id, probe_id, ra_deg, dec_deg) AS (\n"
        "  VALUES\n" + values + "\n)\n"
        "SELECT p.candidate_id, p.probe_id, m.skymap_id, m.tract, m.patch, m.band,\n"
        "       m.naxis1, m.naxis2, m.seeing, m.zeropt,\n"
        "       q.imag_psf_depth, q.iseeing, q.imeas_exists, q.idet_bkgd_exists\n"
        "FROM probes AS p\n"
        "JOIN pdr3_wide.mosaic AS m\n"
        "JOIN public.skymap AS s USING (skymap_id)\n"
        "LEFT JOIN pdr3_wide.patch_qa AS q USING (skymap_id)\n"
        "WHERE m.band = 'HSC-I'\n"
        "  AND patch_contains(s.patch_area, s.wcs, p.ra_deg, p.dec_deg)\n"
        "ORDER BY p.candidate_id, p.probe_id, m.skymap_id;\n"
    )


def prepare(catalog_path: Path, manifest_path: Path, member_path: Path,
            expansion_summary_path: Path, output_dir: Path, batch_size: int = 50) -> dict:
    if output_dir.exists():
        raise FileExistsError(f"refusing to overwrite {output_dir}")
    if batch_size <= 0:
        raise ValueError("batch size must be positive")
    for path in (catalog_path, manifest_path, member_path, expansion_summary_path):
        if not path.is_file():
            raise FileNotFoundError(path)
    if sha256(member_path) != EXPECTED_MEMBER_CATALOG_SHA256:
        raise ValueError("local CAMIRA member catalog hash differs from the reviewed file")
    if sha256(expansion_summary_path) != EXPECTED_EXPANSION_AUDIT_SHA256:
        raise ValueError("S20A expansion summary differs from the reviewed result")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    catalog_hash = sha256(catalog_path)
    if manifest.get("catalog_sha256") != catalog_hash:
        raise ValueError("parent S20A catalog hash does not match sample manifest")
    expansion = json.loads(expansion_summary_path.read_text(encoding="utf-8"))
    scenario = expansion.get("selection_scenarios", {}).get("wide_z_cluster_redshift", {})
    if scenario.get("greedy_additions_after_6arcmin_separation") != EXPECTED_CATALOG_EXPANSIONS:
        raise ValueError("archived expansion count is inconsistent")

    catalog = read_catalog(catalog_path)
    current = [{
        "id": str(row["id"]), "name": row.get("name", str(row["id"])),
        "ra_deg": float(row["ra_deg"]), "dec_deg": float(row["dec_deg"]),
        "z": float(row["z_cl"]), "richness": float(row["richness"]),
        "z_bcg": float(row["z_spec_bcg"]),
    } for row in manifest.get("clusters", [])]
    candidates = select_candidates(catalog, current)
    member_positions = attach_members(candidates, member_path)
    cosmology = Cosmology(H0=float(manifest["h0"]), Om=float(manifest["omega_m"]))
    probes = []
    probes_by_id = {}
    for candidate in candidates:
        candidate_probes = make_probes(candidate, cosmology)
        probes.extend(candidate_probes)
        probes_by_id[str(candidate["id"])] = candidate_probes

    output_dir.mkdir(parents=True, exist_ok=False)
    candidate_path = output_dir / "candidate_centers.tsv"
    with candidate_path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=[
            "id", "name", "ra_deg", "dec_deg", "z_cl", "richness", "z_spec_bcg",
            "n_camira_member_positions", "ra_dec_cell_5deg",
        ], delimiter="\t")
        writer.writeheader()
        for row in candidates:
            ident = str(row["id"])
            writer.writerow({
                "id": ident, "name": row["name"], "ra_deg": f"{row['ra_deg']:.8f}",
                "dec_deg": f"{row['dec_deg']:.8f}", "z_cl": f"{row['z']:.5f}",
                "richness": f"{row['richness']:.3f}", "z_spec_bcg": f"{row['z_bcg']:.5f}",
                "n_camira_member_positions": len(member_positions[ident]),
                "ra_dec_cell_5deg": "%d,%d" % sky_cell(row),
            })
    member_path_out = output_dir / "camira_member_positions.tsv"
    with member_path_out.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream, delimiter="\t")
        writer.writerow(["candidate_id", "member_ra_deg", "member_dec_deg"])
        for ident in sorted(member_positions, key=int):
            for ra, dec in member_positions[ident]:
                writer.writerow([ident, f"{ra:.8f}", f"{dec:.8f}"])
    probe_path = output_dir / "coverage_probes.tsv"
    with probe_path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=[
            "candidate_id", "probe_id", "ra_deg", "dec_deg", "radius_kpc",
        ], delimiter="\t")
        writer.writeheader()
        for row in probes:
            writer.writerow({**row, "ra_deg": f"{row['ra_deg']:.9f}",
                             "dec_deg": f"{row['dec_deg']:.9f}",
                             "radius_kpc": f"{row['radius_kpc']:.1f}"})

    batch_records = []
    for start in range(0, len(candidates), batch_size):
        chunk = candidates[start:start + batch_size]
        chunk_probes = [probe for row in chunk for probe in probes_by_id[str(row["id"])]]
        index = start // batch_size
        path = output_dir / f"coverage_mosaic_{index:03d}.sql"
        path.write_text(render_query(chunk_probes), encoding="utf-8")
        batch_records.append({
            "path": str(path.resolve()), "n_candidates": len(chunk),
            "n_probe_points": len(chunk_probes), "expected_rows": None,
            "sha256": sha256(path),
        })

    output_files = [candidate_path, member_path_out, probe_path] + [Path(r["path"]) for r in batch_records]
    result = {
        "schema": "astrafek-icl-pdr3-coverage-query-plan-v1",
        "status": "prepared_not_submitted",
        "science_measurement": False,
        "cas_query_submitted": False,
        "images_downloaded": False,
        "shapes_downloaded": False,
        "photoz_rows_requested": False,
        "selection": {
            "source": "local CAMIRA S20A member catalog",
            "redshift_rule": "0.20 <= z_cl < 0.70; cluster photo-z only",
            "exclude_current_800": True,
            "minimum_center_separation_arcmin": 6.0,
            "order": "richness descending, then catalog ID ascending",
            "candidate_count": len(candidates),
            "new_5deg_ra_dec_cells": EXPECTED_NEW_5DEG_CELLS,
            "member_positions_verified": sum(map(len, member_positions.values())),
        },
        "coverage_test": {
            "catalog_release": "HSC PDR3 Wide",
            "table": "pdr3_wide.mosaic joined to public.skymap and pdr3_wide.patch_qa",
            "band": "HSC-I (PDR3 metadata merges i/i2 into HSC-I)",
            "geometry": "patch_contains at each candidate center plus eight great-circle points at 1,000 proper kpc",
            "output_columns": ["candidate_id", "probe_id", "skymap_id", "tract", "patch", "band", "naxis1", "naxis2", "seeing", "zeropt", "imag_psf_depth", "iseeing", "imeas_exists", "idet_bkgd_exists"],
            "images_or_catalog_object_rows": False,
            "limitation": "nine probes are a coverage screen, not a proof that every pixel in the 1-Mpc disk is covered or usable",
        },
        "inputs": {
            "catalog": {"path": str(catalog_path.resolve()), "sha256": catalog_hash},
            "sample_manifest": {"path": str(manifest_path.resolve()), "sha256": sha256(manifest_path)},
            "member_catalog": {"path": str(member_path.resolve()), "sha256": sha256(member_path)},
            "prior_expansion_audit": {"path": str(expansion_summary_path.resolve()), "sha256": sha256(expansion_summary_path)},
        },
        "n_candidates": len(candidates), "n_probes": len(probes),
        "batch_size_candidates": batch_size, "n_sql_batches": len(batch_records),
        "batches": batch_records,
        "post_query_gates": [
            "require center and all eight 1-Mpc boundary probes to match an HSC-I mosaic patch",
            "evaluate patch depth/seeing and full aperture support with actual WCS and mask after image authorization",
            "design and verify matched-control positions before image retrieval",
            "recompute connected spatial blocks and sample selection before stacking",
        ],
        "output_sha256": {str(path.resolve()): sha256(path) for path in output_files},
    }
    (output_dir / "coverage_query_plan.json").write_text(
        json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    return result


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--catalog", type=Path, default=DEFAULT_CATALOG)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--members", type=Path, default=DEFAULT_MEMBER_CATALOG)
    parser.add_argument("--expansion-summary", type=Path, default=DEFAULT_EXPANSION_AUDIT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--batch-size", type=int, default=50)
    args = parser.parse_args(argv)
    result = prepare(args.catalog, args.manifest, args.members,
                     args.expansion_summary, args.output_dir, args.batch_size)
    print(json.dumps({
        "status": result["status"], "n_candidates": result["n_candidates"],
        "n_probes": result["n_probes"], "n_sql_batches": result["n_sql_batches"],
        "member_positions_verified": result["selection"]["member_positions_verified"],
        "output_dir": str(args.output_dir.resolve()),
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
