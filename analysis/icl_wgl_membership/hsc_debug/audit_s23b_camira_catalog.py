#!/usr/bin/env python3
"""Audit the approved S23B CAMIRA table against the local S20A/800 sample.

Catalog-only feasibility diagnostic: no query, image, shape, or photo-z fetch.
All reported counts are candidate-pool checks, not scientific measurements.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import platform
from pathlib import Path

import numpy as np
from scipy.spatial import cKDTree


ROOT = Path("/gpfs/kjhan/Astrafek/data")
S23B = ROOT / "catalogue_s23b_camira_2026-10-20_01/camira_s23b_wide_sm.tbl"
S20A = ROOT / "catalog/camira_s20a_wide_sm.tbl"
SAMPLE = ROOT / "sample_800_2026-10-13/sample_manifest.json"
DOWNLOAD_STATUS = ROOT / "sample_800_2026-10-13/download_01/download_status.json"
OUTPUT = ROOT / "catalogue_s23b_camira_audit_2026-10-20_02"
MIN_SEP_DEG = 0.1
MATCH_RADII_ARCSEC = (1.0, 60.0)
EXPECTED_S23B_BYTES = 827_339
EXPECTED_S23B_BLOB_SHA1 = "066c2c8beca8f3a65186db7511fddca52d298908"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_catalog(path: Path) -> list[dict]:
    rows = []
    with Path(path).open(encoding="utf-8") as stream:
        for line_no, line in enumerate(stream, 1):
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            fields = line.split()
            if len(fields) != 7:
                raise ValueError(f"{path}:{line_no}: expected 7 columns, found {len(fields)}")
            ident, name = fields[:2]
            ra, dec, z, richness, z_bcg = map(float, fields[2:])
            values = (ra, dec, z, richness, z_bcg)
            if not all(math.isfinite(value) for value in values):
                raise ValueError(f"{path}:{line_no}: non-finite value")
            if not (0 <= ra < 360 and -90 <= dec <= 90 and z >= 0 and richness >= 0):
                raise ValueError(f"{path}:{line_no}: value outside valid range")
            rows.append({"id": ident, "name": name, "ra_deg": ra, "dec_deg": dec,
                         "z": z, "richness": richness, "z_bcg": z_bcg})
    if len({row["id"] for row in rows}) != len(rows):
        raise ValueError(f"{path}: duplicate IDs")
    return rows


def unit_vectors(rows: list[dict]) -> np.ndarray:
    ra = np.deg2rad([row["ra_deg"] for row in rows])
    dec = np.deg2rad([row["dec_deg"] for row in rows])
    return np.column_stack((np.cos(dec) * np.cos(ra),
                            np.cos(dec) * np.sin(ra), np.sin(dec)))


def nearest_matches(left: list[dict], right: list[dict]) -> tuple[np.ndarray, np.ndarray]:
    tree = cKDTree(unit_vectors(right))
    chord, index = tree.query(unit_vectors(left), k=1)
    angle_deg = np.rad2deg(2 * np.arcsin(np.clip(chord / 2, 0, 1)))
    return angle_deg * 3600, index


def sky_cell(row: dict) -> tuple[int, int]:
    return math.floor(row["ra_deg"] / 5), math.floor((row["dec_deg"] + 90) / 5)


def greedy_additions(current: list[dict], candidates: list[dict]) -> list[dict]:
    vectors = unit_vectors(current)
    chosen = list(current)
    chosen_vectors = [vector for vector in vectors]
    additions = []
    cos_limit = math.cos(math.radians(MIN_SEP_DEG))
    for row in sorted(candidates, key=lambda item: (-item["richness"], item["id"])):
        vector = unit_vectors([row])[0]
        if all(float(vector @ other) < cos_limit for other in chosen_vectors):
            chosen.append(row)
            chosen_vectors.append(vector)
            additions.append(row)
    return additions


def run(s23b_path: Path, s20a_path: Path, sample_path: Path,
        download_status_path: Path, output: Path) -> dict:
    if not platform.node().lower().startswith("lageunha"):
        raise RuntimeError("Run the catalog audit on lageunha CPU")
    output = output.resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite {output}")
    for path in (s23b_path, s20a_path, sample_path, download_status_path):
        if not path.is_file():
            raise FileNotFoundError(path)
    acquisition_path = s23b_path.with_suffix(s23b_path.suffix + ".manifest.json")
    acquisition = json.loads(acquisition_path.read_text(encoding="utf-8"))
    blob_sha1 = hashlib.sha1(f"blob {s23b_path.stat().st_size}\0".encode("ascii")
                             + s23b_path.read_bytes()).hexdigest()
    if (s23b_path.stat().st_size != EXPECTED_S23B_BYTES
            or blob_sha1 != EXPECTED_S23B_BLOB_SHA1
            or acquisition.get("git_blob_sha1") != EXPECTED_S23B_BLOB_SHA1
            or acquisition.get("bytes") != EXPECTED_S23B_BYTES):
        raise ValueError("S23B acquisition manifest/byte/blob verification failed")
    sample_manifest = json.loads(sample_path.read_text(encoding="utf-8"))
    download_status = json.loads(download_status_path.read_text(encoding="utf-8"))
    sample = [{"id": str(row["id"]), "name": row.get("name", str(row["id"])),
               "ra_deg": float(row["ra_deg"]), "dec_deg": float(row["dec_deg"]),
               "z": float(row["z_cl"]), "richness": float(row["richness"]),
               "z_bcg": float(row["z_spec_bcg"])}
              for row in sample_manifest.get("clusters", [])]
    if len(sample) != 800 or len({row["id"] for row in sample}) != 800:
        raise ValueError("expected the existing manifest to contain 800 unique lenses")
    downloaded = {str(row["id"]) for row in download_status.get("completed", [])}
    if (download_status.get("status") != "complete" or len(downloaded) != 800
            or downloaded != {row["id"] for row in sample}):
        raise ValueError("existing image inventory is not a complete 800-target match")

    s23b, s20a = read_catalog(s23b_path), read_catalog(s20a_path)
    sep_s20a, idx_s20a = nearest_matches(s23b, s20a)
    sep_sample, idx_sample = nearest_matches(s23b, sample)
    duplicates = {}
    for radius in MATCH_RADII_ARCSEC:
        duplicates[f"within_{radius:g}_arcsec"] = int(np.count_nonzero(sep_s20a <= radius))
    sample_matches = {f"within_{radius:g}_arcsec": int(np.count_nonzero(sep_sample <= radius))
                      for radius in MATCH_RADII_ARCSEC}

    s20a_match = sep_s20a <= MATCH_RADII_ARCSEC[0]
    eligible = [row for i, row in enumerate(s23b)
                if not s20a_match[i]
                and 0.30 <= row["z"] < 0.55
                and row["z_bcg"] > 0
                and abs(row["z_bcg"] - row["z"]) / (1 + row["z"]) < 0.01]
    additions = greedy_additions(sample, eligible)
    sample_cells = {sky_cell(row) for row in sample}
    addition_cells = {sky_cell(row) for row in additions}

    # Link only to the already validated 800-target image sample by sky position.
    # A miss says nothing about other PDR3 tiles that are not inventoried here.
    close_to_sample = sep_sample <= MATCH_RADII_ARCSEC[0]
    supported_pool = [row for i, row in enumerate(s23b)
                      if close_to_sample[i] and 0.30 <= row["z"] < 0.55
                      and row["z_bcg"] > 0
                      and abs(row["z_bcg"] - row["z"]) / (1 + row["z"]) < 0.01]
    result = {
        "schema": "astrafek-s23b-camira-catalog-audit-v1",
        "status": "complete_catalog_only_feasibility_audit",
        "science_measurement": False,
        "query_submitted": False,
        "images_shape_rows_photoz_downloaded": False,
        "selection": {"z_cl": [0.30, 0.55], "bcg_spec_required": True,
                      "frac_abs_zbcg_minus_zcl_over_1_plus_zcl_lt": 0.01,
                      "minimum_center_separation_arcmin": 6.0,
                      "greedy_order": "descending richness, then catalog ID"},
        "catalogs": {
            "s23b": {"path": str(s23b_path.resolve()), "rows": len(s23b),
                     "sha256": sha256(s23b_path)},
            "s20a": {"path": str(s20a_path.resolve()), "rows": len(s20a),
                     "sha256": sha256(s20a_path)},
            "current_800_manifest": {"path": str(sample_path.resolve()),
                                      "sha256": sha256(sample_path), "rows": len(sample)},
            "current_800_image_inventory": {"path": str(download_status_path.resolve()),
                                             "sha256": sha256(download_status_path),
                                             "status": download_status["status"],
                                             "downloaded": len(downloaded)},
        },
        "nearest_s20a_sky_match_counts": duplicates,
        "nearest_existing_800_target_match_counts": sample_matches,
        "candidate_pool": {
            "eligible_new_after_s20a_1arcsec_dedup": len(eligible),
            "greedy_additions_after_existing_800_and_6arcmin_separation": len(additions),
            "new_5deg_ra_dec_cells": len(addition_cells - sample_cells),
            "total_5deg_ra_dec_cells_after_additions": len(sample_cells | addition_cells),
            "candidate_rows_matching_existing_800_image_positions_and_selection": len(supported_pool),
            "interpretation": "catalog-only upper bound; not verified full-aperture image support or an independent covariance-block count",
        },
        "image_release_limit": "Existing image manifests cover the current 800 targets only. A position match identifies an existing target; misses do not prove PDR3 unavailability. New candidates need a separate exact image-support audit and download approval.",
        "covariance_limit": "5-degree coordinate cells are only an upper bound. The connected target/control block count cannot be recomputed until matched controls are assigned and supported.",
    }
    output.mkdir(parents=True, exist_ok=False)
    (output / "summary.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n",
                                          encoding="utf-8")
    with (output / "greedy_candidate_additions.tsv").open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=["id", "name", "ra_deg", "dec_deg", "z",
                                                    "richness", "z_bcg", "ra_dec_cell_5deg"],
                                delimiter="\t")
        writer.writeheader()
        for row in additions:
            writer.writerow({**row, "ra_dec_cell_5deg": "%d,%d" % sky_cell(row)})
    result["outputs"] = {name: str(output / name) for name in
                         ("summary.json", "greedy_candidate_additions.tsv")}
    # Rewrite summary with its output paths now that all output paths are fixed.
    (output / "summary.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n",
                                          encoding="utf-8")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--s23b", type=Path, default=S23B)
    parser.add_argument("--s20a", type=Path, default=S20A)
    parser.add_argument("--sample", type=Path, default=SAMPLE)
    parser.add_argument("--download-status", type=Path, default=DOWNLOAD_STATUS)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    result = run(args.s23b, args.s20a, args.sample, args.download_status, args.output)
    print(json.dumps({"status": result["status"], "catalog_rows": result["catalogs"]["s23b"]["rows"],
                      "s20a_matches": result["nearest_s20a_sky_match_counts"],
                      "sample_matches": result["nearest_existing_800_target_match_counts"],
                      "candidate_pool": result["candidate_pool"], "outputs": result["outputs"],
                      "science_measurement": False}, indent=2))


if __name__ == "__main__":
    main()
