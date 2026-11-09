#!/usr/bin/env python3
"""Audit whether the local S23B candidate pool adds HSC tract components.

The audit uses only the already saved S23B candidate-center ledger, the HSC
patch-QA catalogue, and the existing 800-lens manifest. It requests no data.
Tract-center components approximate field connectivity and do not define
official HSC survey fields.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
PROJECT_ANALYSIS = HERE.parents[2]
if str(PROJECT_ANALYSIS) not in sys.path:
    sys.path.insert(0, str(PROJECT_ANALYSIS))

from fetch_ancillary_catalogues import sha
from run_wl_hsc_tract_region_influence import (
    nearest_region_labels,
    tract_centers_from_patch_qa,
    tract_region_labels,
)


ROOT = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
DEFAULT_CANDIDATES = ROOT / "icl_s23b_coverage_plan_2026-10-22_01/candidate_centers.tsv"
DEFAULT_PATCH_QA = ROOT / "catalogues_2026-10-13_01/patch_qa.csv"
DEFAULT_DISCOVERY_MANIFEST = ROOT / "download_01/download_manifest.json"
DEFAULT_OUTPUT = ROOT / "icl_s23b_tract_region_support_2026-11-06_01"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with Path(path).open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def run(args):
    output = args.output.resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite output directory: {output}")
    candidates = read_tsv(args.candidates)
    discovery = json.loads(args.discovery_manifest.read_text(encoding="utf-8"))
    discovery_rows = discovery.get("clusters", [])
    if len(candidates) != 1333 or len(discovery_rows) != 800:
        raise ValueError("candidate/discovery samples differ from the pinned local ledgers")
    required = {"candidate_id", "ra_deg", "dec_deg"}
    if not candidates or not required.issubset(candidates[0]):
        raise ValueError("S23B candidate-center ledger lacks identifiers or coordinates")
    candidate_ids = [str(row["candidate_id"]) for row in candidates]
    lens_ids = [str(row["id"]) for row in discovery_rows]
    if len(set(candidate_ids)) != len(candidate_ids) or len(set(lens_ids)) != len(lens_ids):
        raise ValueError("candidate or discovery ledger contains duplicate IDs")
    candidate_ra = np.asarray([float(row["ra_deg"]) for row in candidates])
    candidate_dec = np.asarray([float(row["dec_deg"]) for row in candidates])
    lens_ra = np.asarray([float(row["ra_deg"]) for row in discovery_rows])
    lens_dec = np.asarray([float(row["dec_deg"]) for row in discovery_rows])
    if not all(np.isfinite(x).all() for x in (candidate_ra, candidate_dec, lens_ra, lens_dec)):
        raise ValueError("candidate or discovery coordinates contain non-finite values")

    tract_keys, tract_xyz = tract_centers_from_patch_qa(args.patch_qa)
    tract_labels = tract_region_labels(tract_xyz, args.adjacency_deg)
    candidate_labels = nearest_region_labels(candidate_ra, candidate_dec, tract_xyz, tract_labels)
    discovery_labels = nearest_region_labels(lens_ra, lens_dec, tract_xyz, tract_labels)
    component_rows = []
    for component in sorted(set(map(int, tract_labels.tolist()))):
        member_tracts = [tract_keys[index] for index, label in enumerate(tract_labels)
                         if int(label) == component]
        component_rows.append({
            "component_id": component,
            "n_tracts": len(member_tracts),
            "tract_ids_sha256": hashlib.sha256("\n".join(member_tracts).encode()).hexdigest(),
            "discovery_lenses": int(np.count_nonzero(discovery_labels == component)),
            "s23b_candidates": int(np.count_nonzero(candidate_labels == component)),
            "new_to_discovery": int(np.count_nonzero(discovery_labels == component)) == 0
                                 and int(np.count_nonzero(candidate_labels == component)) > 0,
        })
    if sum(row["discovery_lenses"] for row in component_rows) != 800:
        raise ValueError("tract component labels failed to cover all discovery lenses")
    if sum(row["s23b_candidates"] for row in component_rows) != 1333:
        raise ValueError("tract component labels failed to cover all candidates")

    new_candidates = sum(row["s23b_candidates"] for row in component_rows if row["new_to_discovery"])
    result = {
        "schema": "astrafek-s23b-candidate-tract-region-support-v1",
        "status": "complete_offline_feasibility_audit",
        "science_measurement": False,
        "network_used": False,
        "hsc_query_submitted": False,
        "data_downloaded": False,
        "image_pixels_read": False,
        "method": "assign candidate and discovery lens centers to nearest spherical-mean patch-QA tract center, then group tract centers into connected components",
        "tract_adjacency_deg": float(args.adjacency_deg),
        "n_tracts": len(tract_keys),
        "n_connected_components": len(component_rows),
        "n_discovery_lenses": len(discovery_rows),
        "n_s23b_catalog_candidates": len(candidates),
        "components": component_rows,
        "candidate_centers_in_components_already_occupied": len(candidates) - new_candidates,
        "candidate_centers_in_components_unoccupied_by_discovery": new_candidates,
        "fraction_in_existing_components": float((len(candidates) - new_candidates) / len(candidates)),
        "interpretation": "the catalog pool adds candidates within the existing connected HSC footprint; the small unoccupied component population cannot supply an independent validation sample, and catalog candidates are not images or measured clusters",
        "limitations": [
            "nearest tract-center assignment approximates the footprint and does not use tract polygons",
            "connected components do not establish independent covariance regions",
            "catalog candidate positions do not establish image availability, depth, PSF, mask, ICL support, or satellite support",
        ],
        "inputs": {
            "candidate_centers_sha256": sha(args.candidates),
            "patch_qa_sha256": sha(args.patch_qa),
            "discovery_manifest_sha256": sha(args.discovery_manifest),
            "runner_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        },
    }
    output.mkdir(parents=True, exist_ok=False)
    with (output / "summary.json").open("x", encoding="utf-8") as stream:
        json.dump(result, stream, ensure_ascii=False, indent=2, allow_nan=False)
        stream.write("\n")
    print(json.dumps({
        "status": result["status"], "science_measurement": False,
        "tract_components": len(component_rows), "discovery_lenses": 800,
        "s23b_candidates": 1333, "candidates_in_new_component": new_candidates,
        "output": str(output / "summary.json"),
    }, ensure_ascii=False))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidates", type=Path, default=DEFAULT_CANDIDATES)
    parser.add_argument("--patch-qa", type=Path, default=DEFAULT_PATCH_QA)
    parser.add_argument("--discovery-manifest", type=Path, default=DEFAULT_DISCOVERY_MANIFEST)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--adjacency-deg", type=float, default=2.0)
    run(parser.parse_args())


if __name__ == "__main__":
    main()
