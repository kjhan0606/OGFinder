#!/usr/bin/env python3
"""Prepare, but never submit, HSC-I coverage metadata SQL for S23B candidates.

This plan is specific to the 1,333 S23B additions relevant to the current
CAMIRA selection. It reuses the already-reviewed PDR3 metadata query shape and
does not contact HSC, retrieve metadata, or download science data.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from prepare_icl_expansion_coverage_plan import make_probes, render_query
from ogfkit.wlcluster import Cosmology


DATA = Path("/gpfs/kjhan/Astrafek/data")
SAMPLE = DATA / "sample_800_2026-10-13"
DEFAULT_CANDIDATES = DATA / "catalogue_s23b_camira_audit_2026-10-20_02/greedy_candidate_additions.tsv"
DEFAULT_PARENT = SAMPLE / "sample_manifest.json"
DEFAULT_LEDGER = SAMPLE / "icl_s23b_spatial_candidate_ledger_2026-10-22_01/summary.json"
DEFAULT_PRIOR_PLAN = SAMPLE / "icl_expansion_coverage_plan_2026-10-21_01/coverage_query_plan.json"
DEFAULT_OUTPUT = SAMPLE / "icl_s23b_coverage_plan_2026-10-22_01"

EXPECTED = {
    "candidates": "1b59ca9941849fa7970c43c6e0374db6587716fc4ae5d71b11f84a245ebe6b2c",
    "parent": "95388814062e95f6e8429e2543a5432433515c9d9607680d2e873488d2758401",
    "ledger": "f6ac867881d8184d02a580e9a18540b7513a90f530e51e6cdf04205db535e8d1",
    "prior_plan": "9da220c05f1f24c409a602fcc4ddffc8cbd5bc33f459d335cc8587f1cd235006",
}
BATCH_SIZE = 50
PROBES_PER_CANDIDATE = 9
OUTPUT_COLUMNS = [
    "candidate_id", "probe_id", "skymap_id", "tract", "patch", "band",
    "naxis1", "naxis2", "seeing", "zeropt", "imag_psf_depth", "iseeing",
    "imeas_exists", "idet_bkgd_exists",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_candidates(path: Path) -> list[dict]:
    with path.open(encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    required = {"id", "name", "ra_deg", "dec_deg", "z", "richness"}
    if not rows or not required.issubset(rows[0]):
        raise ValueError("S23B candidate table lacks required columns")
    result, seen = [], set()
    for row in rows:
        source_id = str(row["id"])
        ra, dec, z = map(float, (row["ra_deg"], row["dec_deg"], row["z"]))
        richness = float(row["richness"])
        if (not source_id.isdigit() or source_id in seen
                or not all(map(math.isfinite, (ra, dec, z, richness)))
                or not (0 <= ra < 360 and -90 <= dec <= 90 and z > 0 and richness > 0)):
            raise ValueError("invalid or duplicate S23B candidate row")
        seen.add(source_id)
        result.append({
            "id": f"S23B:{source_id}", "source_id": source_id,
            "name": row["name"], "ra_deg": ra, "dec_deg": dec,
            "z": z, "richness": richness,
        })
    if len(result) != 1333:
        raise ValueError(f"expected 1,333 reviewed S23B candidates, got {len(result)}")
    return result


def build_probes(candidates: list[dict], cosmology: Cosmology) -> list[dict]:
    probes = []
    for candidate in candidates:
        probes.extend(make_probes(candidate, cosmology))
    if len(probes) != len(candidates) * PROBES_PER_CANDIDATE:
        raise ValueError("probe generation did not produce center + 8 boundary points")
    keys = {(row["candidate_id"], row["probe_id"]) for row in probes}
    if len(keys) != len(probes):
        raise ValueError("duplicate candidate/probe key")
    return probes


def prepare(candidate_path: Path, parent_path: Path, ledger_path: Path,
            prior_plan_path: Path, output_dir: Path) -> dict:
    output_dir = output_dir.resolve()
    if output_dir.exists():
        raise FileExistsError(f"refusing to overwrite {output_dir}")
    paths = {
        "candidates": candidate_path.resolve(), "parent": parent_path.resolve(),
        "ledger": ledger_path.resolve(), "prior_plan": prior_plan_path.resolve(),
    }
    for key, path in paths.items():
        if not path.is_file():
            raise FileNotFoundError(path)
        if sha256(path) != EXPECTED[key]:
            raise ValueError(f"reviewed {key} input checksum mismatch")
    parent = json.loads(parent_path.read_text(encoding="utf-8"))
    ledger = json.loads(ledger_path.read_text(encoding="utf-8"))
    prior = json.loads(prior_plan_path.read_text(encoding="utf-8"))
    if (parent.get("schema") != "hsc-bcg-diffuse-sample-v1"
            or parent.get("science_measurement") is not False
            or len(parent.get("clusters", [])) != 800):
        raise ValueError("unexpected current parent manifest")
    if (ledger.get("schema") != "astrafek-s23b-spatial-candidate-ledger-v1"
            or ledger.get("status") != "offline_coordinate_and_footprint_ledger_complete"
            or ledger.get("statistics", {}).get("candidate_count") != 1333
            or ledger.get("science_measurement") is not False):
        raise ValueError("unexpected S23B spatial ledger")
    if (prior.get("schema") != "astrafek-icl-pdr3-coverage-query-plan-v1"
            or prior.get("status") != "prepared_not_submitted"
            or prior.get("cas_query_submitted") is not False
            or prior.get("n_candidates") != 2596 or prior.get("n_probes") != 23364):
        raise ValueError("existing approved S20A coverage plan changed or was submitted")

    candidates = read_candidates(candidate_path)
    cosmology = Cosmology(H0=float(parent["h0"]), Om=float(parent["omega_m"]))
    probes = build_probes(candidates, cosmology)
    if len(probes) != 11997:
        raise ValueError("expected exactly 11,997 metadata probes")

    output_dir.mkdir(parents=True, exist_ok=False)
    centers_path = output_dir / "candidate_centers.tsv"
    with centers_path.open("x", encoding="utf-8", newline="") as stream:
        fields = ["candidate_id", "catalog_id", "name", "ra_deg", "dec_deg", "z_cl", "richness"]
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for row in candidates:
            writer.writerow({
                "candidate_id": row["id"], "catalog_id": row["source_id"], "name": row["name"],
                "ra_deg": f"{row['ra_deg']:.8f}", "dec_deg": f"{row['dec_deg']:.8f}",
                "z_cl": f"{row['z']:.5f}", "richness": f"{row['richness']:.5f}",
            })
    probes_path = output_dir / "coverage_probes.tsv"
    with probes_path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=[
            "candidate_id", "probe_id", "ra_deg", "dec_deg", "radius_kpc",
        ], delimiter="\t")
        writer.writeheader()
        for row in probes:
            writer.writerow({**row, "ra_deg": f"{row['ra_deg']:.9f}",
                             "dec_deg": f"{row['dec_deg']:.9f}",
                             "radius_kpc": f"{row['radius_kpc']:.1f}"})

    by_candidate = {}
    for probe in probes:
        by_candidate.setdefault(probe["candidate_id"], []).append(probe)
    batches = []
    for start in range(0, len(candidates), BATCH_SIZE):
        chunk = candidates[start:start + BATCH_SIZE]
        sql_path = output_dir / f"coverage_mosaic_{start // BATCH_SIZE:03d}.sql"
        sql_path.write_text(render_query([
            probe for candidate in chunk for probe in by_candidate[candidate["id"]]
        ]), encoding="utf-8")
        batches.append({
            "path": str(sql_path), "n_candidates": len(chunk),
            "n_probe_points": sum(len(by_candidate[row["id"]]) for row in chunk),
            "expected_rows": None, "sha256": sha256(sql_path),
        })
    output_files = [centers_path, probes_path] + [Path(row["path"]) for row in batches]
    result = {
        "schema": "astrafek-s23b-pdr3-coverage-query-plan-v1",
        "status": "prepared_pending_user_approval",
        "science_measurement": False,
        "user_approval_required_before_submission": True,
        "cas_query_submitted": False,
        "data_downloaded": False,
        "images_downloaded": False,
        "shape_or_photoz_rows_requested": False,
        "selection": {
            "source": "pinned S23B CAMIRA Wide candidate additions",
            "candidate_count": len(candidates),
            "catalog_id_namespace": "S23B-prefixed; never joined to S20A by numeric ID",
            "current_association_new_target_cells_i2": 67,
            "candidate_centers_in_existing_target_cutouts": 0,
            "candidate_centers_in_existing_control_cutouts": 7,
        },
        "coverage_test": {
            "catalog_release": "HSC PDR3 Wide",
            "table": "pdr3_wide.mosaic joined to public.skymap and pdr3_wide.patch_qa",
            "band": "HSC-I (PDR3 metadata merges i/i2 into HSC-I)",
            "geometry": "patch_contains at each candidate center plus eight great-circle points at 1,000 proper kpc",
            "output_columns": OUTPUT_COLUMNS,
            "images_or_catalog_object_rows": False,
            "estimated_max_metadata_rows": len(probes) * 4,
            "limitation": "nine points screen footprint only; they do not prove full-disk coverage, masks, or matched-control availability",
        },
        "prior_approved_plan": {
            "path": str(prior_plan_path.resolve()),
            "sha256": sha256(prior_plan_path),
            "status": prior.get("status"),
            "n_candidates": prior.get("n_candidates"),
            "n_probes": prior.get("n_probes"),
            "cas_query_submitted": prior.get("cas_query_submitted"),
            "note": "This 2,596-candidate S20A plan is a distinct approved scope; it does not substitute for these S23B positions.",
        },
        "inputs_sha256": {str(path): sha256(path) for path in paths.values()},
        "n_candidates": len(candidates), "n_probes": len(probes),
        "batch_size_candidates": BATCH_SIZE, "n_sql_batches": len(batches),
        "batches": batches,
        "post_query_gates": [
            "verify candidate centers and all eight 1-Mpc boundary probes against HSC-I mosaic coverage",
            "screen depth, seeing, image-product flags, and actual WCS/masks before image request",
            "design and verify uncontaminated matched controls for each eligible target",
            "recompute usable connected spatial blocks before stacking",
        ],
        "output_sha256": {str(path.resolve()): sha256(path) for path in output_files},
    }
    plan_path = output_dir / "coverage_query_plan.json"
    plan_path.write_text(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False) + "\n",
                         encoding="utf-8")
    return result


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidates", type=Path, default=DEFAULT_CANDIDATES)
    parser.add_argument("--parent", type=Path, default=DEFAULT_PARENT)
    parser.add_argument("--ledger", type=Path, default=DEFAULT_LEDGER)
    parser.add_argument("--prior-plan", type=Path, default=DEFAULT_PRIOR_PLAN)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args(argv)
    result = prepare(args.candidates, args.parent, args.ledger, args.prior_plan, args.output_dir)
    print(json.dumps({
        "status": result["status"], "n_candidates": result["n_candidates"],
        "n_probes": result["n_probes"], "n_sql_batches": result["n_sql_batches"],
        "estimated_max_metadata_rows": result["coverage_test"]["estimated_max_metadata_rows"],
        "cas_query_submitted": result["cas_query_submitted"],
        "output_dir": str(args.output_dir.resolve()),
    }, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
