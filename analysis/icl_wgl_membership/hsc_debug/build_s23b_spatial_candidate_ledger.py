#!/usr/bin/env python3
"""Build an offline, coordinate-keyed S23B candidate/ICL footprint ledger.

The ledger compares an already-audited S23B candidate list with the current
800-target manifest, the 94-row ICL--satellite analysis sample, and the saved
header/WCS-only candidate audit. It does not query HSC, download data, inspect
image pixels, select a final science sample, or make a science measurement.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

from prestack_core import spatial_labels


DATA = Path("/gpfs/kjhan/Astrafek/data")
SAMPLE = DATA / "sample_800_2026-10-13"
DEFAULT_CANDIDATES = DATA / "catalogue_s23b_camira_audit_2026-10-20_02/greedy_candidate_additions.tsv"
DEFAULT_PARENT = SAMPLE / "sample_manifest.json"
DEFAULT_ASSOCIATION = SAMPLE / "satellite_outer_icl_association_perbin_2026-10-18_03"
DEFAULT_WCS_AUDIT = SAMPLE / "icl_s23b_existing_cutout_support_2026-10-22_01/summary.json"
DEFAULT_OUTPUT = SAMPLE / "icl_s23b_spatial_candidate_ledger_2026-10-22_01"

EXPECTED_SHA256 = {
    "candidates": "1b59ca9941849fa7970c43c6e0374db6587716fc4ae5d71b11f84a245ebe6b2c",
    "parent": "95388814062e95f6e8429e2543a5432433515c9d9607680d2e873488d2758401",
    "association_summary": "df50d282910e1aea60bd37f190ad3488940a825c786c9eff8d25890abed82ea4",
    "association_rows": "9a6befa87a81257f8faf1d34fa7c019a91d64484fbfa64136bd93496dfc7b58c",
    "wcs_audit": "071377813fcfbf15f298a5d48f3e7a3f0b74f60899d88c46eafb7d5bb58ffebb",
}
ASSOCIATION_COUNTS = {"i": 26, "i2": 68}
MATCH_RADIUS_ARCSEC = 1.0


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def read_candidates(path: Path) -> list[dict]:
    rows = read_tsv(path)
    required = {"id", "name", "ra_deg", "dec_deg", "z", "richness"}
    if not rows or not required.issubset(rows[0]):
        raise ValueError("S23B candidate table is missing required columns")
    seen = set()
    result = []
    for row in rows:
        ident = str(row["id"])
        ra, dec = float(row["ra_deg"]), float(row["dec_deg"])
        z, richness = float(row["z"]), float(row["richness"])
        if (not ident or ident in seen or not all(map(math.isfinite, (ra, dec, z, richness)))
                or not 0 <= ra < 360 or not -90 <= dec <= 90):
            raise ValueError("invalid or duplicate S23B candidate")
        seen.add(ident)
        result.append({"id": ident, "name": row["name"], "ra_deg": ra,
                       "dec_deg": dec, "z": z, "richness": richness})
    if len(result) != 1333:
        raise ValueError(f"expected 1,333 reviewed candidates, found {len(result)}")
    return result


def sky_cell(ra: float, dec: float) -> str:
    # Must match prestack_core.spatial_labels, which defines current 5-degree
    # rectangular RA/Dec bins used in the covariance diagnostics.
    return str(spatial_labels([ra], [dec], width=5.)[0])


def angular_sep_arcsec(ra1: float, dec1: float, ra2: float, dec2: float) -> float:
    a1, d1, a2, d2 = map(math.radians, (ra1, dec1, ra2, dec2))
    sd = math.sin((d2-d1)/2.)
    sr = math.sin((a2-a1)/2.)
    h = sd*sd + math.cos(d1)*math.cos(d2)*sr*sr
    return math.degrees(2.*math.asin(min(1., math.sqrt(max(0., h))))) * 3600.


def build_ledger(candidates: list[dict], parent: list[dict],
                 association: list[dict], control_hits: set[str]) -> tuple[list[dict], dict]:
    parent_by_id = {}
    parent_cells = set()
    for row in parent:
        ident = str(row["id"])
        if ident in parent_by_id:
            raise ValueError("duplicate parent target ID")
        record = {"id": ident, "ra_deg": float(row["ra_deg"]), "dec_deg": float(row["dec_deg"])}
        parent_by_id[ident] = record
        parent_cells.add(sky_cell(record["ra_deg"], record["dec_deg"]))
    if len(parent_by_id) != 800:
        raise ValueError(f"expected 800 current targets, found {len(parent_by_id)}")

    association_by_filter: dict[str, set[str]] = {"i": set(), "i2": set()}
    association_keys = set()
    for row in association:
        ident, band = str(row["cluster_id"]), str(row["filter"])
        if band not in association_by_filter or (ident, band) in association_keys:
            raise ValueError("unexpected filter or duplicate ICL--satellite row")
        if ident not in parent_by_id:
            raise ValueError(f"association target {ident} is absent from the parent manifest")
        association_keys.add((ident, band))
        p = parent_by_id[ident]
        association_by_filter[band].add(sky_cell(p["ra_deg"], p["dec_deg"]))
    counts = {band: sum(key[1] == band for key in association_keys) for band in association_by_filter}
    if counts != ASSOCIATION_COUNTS:
        raise ValueError(f"unexpected association row counts: {counts}")

    rows = []
    candidate_cells = set()
    for candidate in candidates:
        cell = sky_cell(candidate["ra_deg"], candidate["dec_deg"])
        candidate_cells.add(cell)
        nearest = min(angular_sep_arcsec(candidate["ra_deg"], candidate["dec_deg"],
                                         p["ra_deg"], p["dec_deg"])
                      for p in parent_by_id.values())
        rows.append({
            "candidate_key": f"S23B:{candidate['id']}",
            "catalog_id": candidate["id"],
            "name": candidate["name"],
            "ra_deg": f"{candidate['ra_deg']:.8f}",
            "dec_deg": f"{candidate['dec_deg']:.8f}",
            "z": f"{candidate['z']:.5f}",
            "richness": f"{candidate['richness']:.5f}",
            "spatial_cell_5deg": cell,
            "nearest_parent_sep_arcsec": f"{nearest:.6f}",
            "coordinate_match_parent_within_1arcsec": nearest <= MATCH_RADIUS_ARCSEC,
            "cell_in_current_800_targets": cell in parent_cells,
            "cell_in_current_i_association": cell in association_by_filter["i"],
            "cell_in_current_i2_association": cell in association_by_filter["i2"],
            "center_inside_existing_target_cutout": False,
            "center_inside_existing_control_cutout": candidate["id"] in control_hits,
            "new_target_coadd_support": "unverified_not_downloaded",
            "mask_depth_psf_support": "unverified",
            "satellite_member_support": "unverified",
            "weak_lensing_shape_support": "unverified",
            "matched_control_pair": "unknown",
        })

    stats = {
        "candidate_count": len(candidates),
        "candidate_unique_5deg_cells": len(candidate_cells),
        "parent_target_count": len(parent_by_id),
        "parent_unique_5deg_cells": len(parent_cells),
        "candidate_cells_overlapping_parent": len(candidate_cells & parent_cells),
        "candidate_cells_new_vs_parent": len(candidate_cells - parent_cells),
        "parent_coordinate_matches_within_1arcsec": sum(
            row["coordinate_match_parent_within_1arcsec"] for row in rows),
        "association_rows": counts,
        "association_unique_target_cells": {
            band: len(cells) for band, cells in association_by_filter.items()},
        "candidate_cells_overlapping_association": {
            band: len(candidate_cells & cells) for band, cells in association_by_filter.items()},
        "candidate_cells_new_vs_association": {
            band: len(candidate_cells - cells) for band, cells in association_by_filter.items()},
        "candidate_centers_in_existing_target_cutouts": 0,
        "candidate_centers_in_existing_control_cutouts": len(control_hits),
        "control_overlap_candidate_ids": sorted(control_hits),
        "candidate_cell_count_is_only_a_target_footprint_upper_bound": True,
    }
    return rows, stats


def run(candidate_path: Path, parent_path: Path, association_dir: Path,
        wcs_audit_path: Path) -> tuple[list[dict], dict]:
    association_summary_path = association_dir / "summary.json"
    association_path = association_dir / "matched_cluster_values.tsv"
    inputs = {
        "candidates": candidate_path.resolve(), "parent": parent_path.resolve(),
        "association_summary": association_summary_path.resolve(),
        "association_rows": association_path.resolve(), "wcs_audit": wcs_audit_path.resolve(),
    }
    for key, path in inputs.items():
        if not path.is_file():
            raise FileNotFoundError(path)
        if sha256(path) != EXPECTED_SHA256[key]:
            raise ValueError(f"reviewed {key} input hash mismatch")
    parent_json = json.loads(parent_path.read_text(encoding="utf-8"))
    assoc_summary = json.loads(association_summary_path.read_text(encoding="utf-8"))
    wcs_audit = json.loads(wcs_audit_path.read_text(encoding="utf-8"))
    if parent_json.get("schema") != "hsc-bcg-diffuse-sample-v1" or parent_json.get("science_measurement") is not False:
        raise ValueError("unexpected current parent manifest")
    if (assoc_summary.get("schema") != "hsc-satellite-outer-icl-correlation-v1"
            or assoc_summary.get("status") != "complete"
            or assoc_summary.get("science_measurement") is not False
            or assoc_summary.get("outputs_sha256", {}).get("matched_cluster_values.tsv") != EXPECTED_SHA256["association_rows"]):
        raise ValueError("unexpected ICL--satellite association summary")
    if (wcs_audit.get("schema") != "astrafek-s23b-candidate-existing-cutout-support-v1"
            or wcs_audit.get("status") != "header_wcs_audit_complete"
            or wcs_audit.get("science_measurement") is not False
            or wcs_audit.get("hsc_query_submitted") is not False
            or wcs_audit.get("data_downloaded") is not False
            or wcs_audit.get("image_pixels_read") is not False
            or wcs_audit.get("candidate_table_sha256") != EXPECTED_SHA256["candidates"]
            or wcs_audit.get("target", {}).get("n_unique_candidates_inside") != 0):
        raise ValueError("unexpected candidate WCS audit or expanded scope")
    control_hits = {
        str(candidate["candidate_id"])
        for match in wcs_audit.get("matched_controls", {}).get("matches", [])
        for candidate in match.get("candidates", [])
    }
    if len(control_hits) != wcs_audit.get("matched_controls", {}).get("n_unique_candidates_inside"):
        raise ValueError("matched-control candidate-center count is inconsistent")
    candidates = read_candidates(candidate_path)
    parent = parent_json.get("clusters", [])
    association = read_tsv(association_path)
    rows, stats = build_ledger(candidates, parent, association, control_hits)
    result = {
        "schema": "astrafek-s23b-spatial-candidate-ledger-v1",
        "status": "offline_coordinate_and_footprint_ledger_complete",
        "science_measurement": False,
        "hsc_query_submitted": False,
        "data_downloaded": False,
        "image_pixels_read": False,
        "coordinate_join_key": "sky coordinates; catalog IDs are release-prefixed and never joined across releases",
        "spatial_cell_definition": "prestack_core.spatial_labels: floor(RA/5), floor((Dec+90)/5)",
        "statistics": stats,
        "inputs_sha256": {str(path): sha256(path) for path in inputs.values()},
        "limitations": [
            "The S23B candidate table is a catalog-only pool and does not establish coadd, mask, depth, PSF, member, shape, or matched-control support.",
            "Only candidate centers inside existing target/control image WCS footprints were checked; a control overlap is a possible contaminant flag, not a flux measurement.",
            "New target 5-degree cells are a target-footprint upper bound on added connected pair blocks, not a count of usable or independent blocks; controls can merge components.",
            "The 5-degree rectangular cells are a coarse covariance grouping and do not prove statistical independence.",
            "No final candidate selection, HSC query, data transfer, pixel measurement, or science inference was performed.",
        ],
    }
    return rows, result


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidates", type=Path, default=DEFAULT_CANDIDATES)
    parser.add_argument("--parent", type=Path, default=DEFAULT_PARENT)
    parser.add_argument("--association-dir", type=Path, default=DEFAULT_ASSOCIATION)
    parser.add_argument("--wcs-audit", type=Path, default=DEFAULT_WCS_AUDIT)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args(argv)
    out = args.output_dir.resolve()
    if out.exists():
        parser.error(f"refusing to overwrite {out}")
    rows, result = run(args.candidates, args.parent, args.association_dir, args.wcs_audit)
    out.mkdir(parents=True, exist_ok=False)
    ledger_path = out / "candidate_ledger.tsv"
    with ledger_path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    result["ledger"] = {"path": str(ledger_path), "rows": len(rows), "sha256": sha256(ledger_path)}
    result["analysis_script"] = str(Path(__file__).resolve())
    result["analysis_script_sha256"] = sha256(Path(__file__).resolve())
    summary_path = out / "summary.json"
    summary_path.write_text(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({"status": result["status"], **result["statistics"],
                      "ledger": str(ledger_path), "summary": str(summary_path)},
                     ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
