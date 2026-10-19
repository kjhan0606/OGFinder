#!/usr/bin/env python3
"""Estimate additional separated CAMIRA clusters from the saved parent catalog.

Offline count/planning diagnostic only: no CAS request, coordinate export,
image download, or science measurement is performed.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path


DEFAULT_CATALOG = Path("/gpfs/kjhan/Astrafek/data/catalog/camira_s20a_wide_sm.tbl")
DEFAULT_MANIFEST = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/sample_manifest.json")
DEFAULT_OUTPUT = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/camira_expansion_pool_2026-10-19_01/summary.json")
MIN_SEP_DEG = 0.1
SCENARIOS = (
    ("current_selection", 0.30, 0.55, True),
    ("broader_z_with_bcg_spec", 0.25, 0.65, True),
    ("wide_z_with_bcg_spec", 0.20, 0.70, True),
    ("current_z_cluster_redshift", 0.30, 0.55, False),
    ("broader_z_cluster_redshift", 0.25, 0.65, False),
    ("wide_z_cluster_redshift", 0.20, 0.70, False),
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def angular_sep_deg(a: dict, b: dict) -> float:
    ra1, dec1 = math.radians(float(a["ra_deg"])), math.radians(float(a["dec_deg"]))
    ra2, dec2 = math.radians(float(b["ra_deg"])), math.radians(float(b["dec_deg"]))
    cos_sep = (math.sin(dec1) * math.sin(dec2)
               + math.cos(dec1) * math.cos(dec2) * math.cos(ra1 - ra2))
    return math.degrees(math.acos(min(1.0, max(-1.0, cos_sep))))


def greedy_additions(current: list[dict], candidates: list[dict], minimum_sep_deg: float) -> list[dict]:
    """Apply the archived richness-priority non-overlap rule to new candidates."""
    if not math.isfinite(minimum_sep_deg) or minimum_sep_deg < 0:
        raise ValueError("minimum separation must be finite and non-negative")
    occupied_ids = {str(row["id"]) for row in current}
    chosen = list(current)
    additions = []
    ordered = sorted(candidates, key=lambda row: (-float(row["richness"]), str(row["id"])))
    for row in ordered:
        ident = str(row["id"])
        if ident in occupied_ids:
            continue
        if all(angular_sep_deg(row, other) > minimum_sep_deg for other in chosen):
            chosen.append(row)
            additions.append(row)
            occupied_ids.add(ident)
    return additions


def sky_cell(row: dict, width_deg: float = 5.0) -> tuple[int, int]:
    ra, dec = float(row["ra_deg"]), float(row["dec_deg"])
    if not (math.isfinite(ra) and 0 <= ra < 360 and math.isfinite(dec) and -90 <= dec <= 90):
        raise ValueError("cluster coordinate is outside valid sky bounds")
    return math.floor(ra / width_deg), math.floor((dec + 90.0) / width_deg)


def read_catalog(path: Path) -> list[dict]:
    rows = []
    with Path(path).open(encoding="utf-8") as stream:
        for line_no, line in enumerate(stream, 1):
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            fields = line.split()
            if len(fields) != 7:
                raise ValueError(f"CAMIRA row {line_no} has {len(fields)} fields, expected 7")
            ident, name = fields[:2]
            ra, dec, z, richness, z_bcg = map(float, fields[2:])
            if not (math.isfinite(ra) and 0 <= ra < 360 and math.isfinite(dec)
                    and -90 <= dec <= 90 and math.isfinite(z) and z >= 0
                    and math.isfinite(richness) and math.isfinite(z_bcg)):
                raise ValueError(f"CAMIRA row {line_no} contains invalid values")
            rows.append({"id": ident, "name": name, "ra_deg": ra, "dec_deg": dec,
                         "z": z, "richness": richness, "z_bcg": z_bcg})
    if len({row["id"] for row in rows}) != len(rows):
        raise ValueError("CAMIRA catalog contains duplicate cluster IDs")
    return rows


def evaluate(catalog: list[dict], current: list[dict], minimum_sep_deg: float = MIN_SEP_DEG) -> dict:
    current_ids = {str(row["id"]) for row in current}
    if len(current_ids) != 800:
        raise ValueError("the current manifest must contain 800 unique cluster IDs")
    if not current_ids.issubset({str(row["id"]) for row in catalog}):
        raise ValueError("some current cluster IDs are absent from the parent catalog")
    current_cells = {sky_cell(row) for row in current}
    scenarios = {}
    for name, z_min, z_max, require_bcg_spec in SCENARIOS:
        eligible = []
        for row in catalog:
            if str(row["id"]) in current_ids or not z_min <= row["z"] < z_max:
                continue
            has_spec = row["z_bcg"] > 0
            if require_bcg_spec:
                if not has_spec:
                    continue
                if abs(row["z_bcg"] - row["z"]) / (1 + row["z"]) >= 0.01:
                    continue
            eligible.append(row)
        additions = greedy_additions(current, eligible, minimum_sep_deg)
        new_cells = {sky_cell(row) for row in additions} - current_cells
        scenarios[name] = {
            "redshift_range": [z_min, z_max],
            "requires_consistent_bcg_spectroscopic_redshift": require_bcg_spec,
            "eligible_additional_catalog_rows_before_separation": len(eligible),
            "greedy_additions_after_6arcmin_separation": len(additions),
            "additional_5deg_cells_occupied": len(new_cells),
            "total_5deg_cells_occupied": len(current_cells | {sky_cell(row) for row in additions}),
            "interpretation": "catalog pool estimate only; not a selected sample or download request",
        }
    return {
        "schema": "camira-offline-expansion-pool-audit-v1",
        "status": "complete_catalog_pool_diagnostic",
        "science_measurement": False,
        "cas_query_submitted": False,
        "coordinates_exported": False,
        "images_downloaded": False,
        "current_sample_n": len(current),
        "catalog_rows": len(catalog),
        "minimum_center_separation_arcmin": minimum_sep_deg * 60,
        "current_occupied_5deg_cells": len(current_cells),
        "selection_scenarios": scenarios,
        "limitations": [
            "spectroscopic-redshift scenarios reproduce the existing cut; cluster-redshift scenarios relax it and need separate photo-z and lensing calibration",
            "greedy richness ordering estimates a pool, not a globally optimal maximum independent set",
            "a newly occupied 5-degree coordinate cell does not guarantee an independent target/control covariance block",
            "member positions, image availability, masking, depth, PSF, and redshift-selection effects were not tested for new candidates",
        ],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--catalog", type=Path, default=DEFAULT_CATALOG)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--minimum-separation-arcmin", type=float, default=6.0)
    args = parser.parse_args()
    output = args.output.resolve()
    if output.exists():
        parser.error(f"refusing to overwrite {output}")
    manifest = json.loads(args.manifest.read_text(encoding="utf-8"))
    catalog_digest = sha256(args.catalog)
    if catalog_digest != manifest.get("catalog_sha256"):
        raise ValueError("local CAMIRA parent catalog does not match the 800-sample manifest")
    catalog = read_catalog(args.catalog)
    current = [{
        "id": str(row["id"]), "name": row.get("name", str(row["id"])),
        "ra_deg": float(row["ra_deg"]), "dec_deg": float(row["dec_deg"]),
        "z": float(row["z_cl"]), "richness": float(row["richness"]),
        "z_bcg": float(row["z_spec_bcg"]),
    } for row in manifest.get("clusters", [])]
    result = evaluate(catalog, current, args.minimum_separation_arcmin / 60.0)
    result["inputs"] = {
        "catalog_path": str(args.catalog.resolve()),
        "catalog_sha256": catalog_digest,
        "sample_manifest_path": str(args.manifest.resolve()),
        "sample_manifest_sha256": sha256(args.manifest),
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False) + "\n")
    print(json.dumps({
        "status": result["status"],
        "current_sample": result["current_sample_n"],
        "parent_catalog_rows": result["catalog_rows"],
        "current_occupied_5deg_cells": result["current_occupied_5deg_cells"],
        "scenarios": {name: {k: v[k] for k in (
            "eligible_additional_catalog_rows_before_separation",
            "greedy_additions_after_6arcmin_separation",
            "additional_5deg_cells_occupied",
        )} for name, v in result["selection_scenarios"].items()},
        "science_measurement": False,
        "output": str(output),
        "output_sha256": sha256(output),
    }, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
