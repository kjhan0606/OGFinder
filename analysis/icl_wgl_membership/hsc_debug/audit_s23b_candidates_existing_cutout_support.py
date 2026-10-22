#!/usr/bin/env python3
"""Audit whether extra CAMIRA candidates fall inside already downloaded cutouts.

Reads only local manifests and FITS headers/WCS. It does not read image pixels,
contact HSC, submit queries, or download files. A candidate inside a random
control cutout flags a possible control-field cluster contaminant; it is not
automatically accepted as a new science lens.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

from astropy.io import fits
from astropy.wcs import WCS


DATA_ROOT = Path("/gpfs/kjhan/Astrafek/data")
SAMPLE_ROOT = DATA_ROOT / "sample_800_2026-10-13"
DEFAULT_CANDIDATE_TABLE = DATA_ROOT / "catalogue_s23b_camira_audit_2026-10-20_02/greedy_candidate_additions.tsv"
DEFAULT_OUTPUT = SAMPLE_ROOT / "icl_s23b_existing_cutout_support_2026-10-22_01/summary.json"
EXPECTED_CANDIDATE_SHA256 = "1b59ca9941849fa7970c43c6e0374db6587716fc4ae5d71b11f84a245ebe6b2c"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_candidates(path: Path) -> list[dict]:
    if sha256(path) != EXPECTED_CANDIDATE_SHA256:
        raise ValueError("candidate list differs from the reviewed S23B table")
    with path.open(encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        required = {"id", "ra_deg", "dec_deg", "z", "richness"}
        if not reader.fieldnames or not required.issubset(reader.fieldnames):
            raise ValueError("candidate table omits required CAMIRA columns")
        candidates = []
        seen = set()
        for row in reader:
            ident = str(row["id"])
            ra, dec = float(row["ra_deg"]), float(row["dec_deg"])
            z, richness = float(row["z"]), float(row["richness"])
            if (not ident or ident in seen or not all(map(math.isfinite, (ra, dec, z, richness)))
                    or not 0 <= ra < 360 or not -90 <= dec <= 90):
                raise ValueError("candidate table contains duplicate IDs or invalid coordinates")
            seen.add(ident)
            candidates.append({"id": ident, "ra_deg": ra, "dec_deg": dec,
                               "z": z, "richness": richness})
    if len(candidates) != 1333:
        raise ValueError("expected the pinned 1,333-candidate S23B expansion table")
    return candidates


def image_records(manifest_path: Path, status_path: Path, image_dir: Path,
                  role: str) -> tuple[list[dict], dict[str, str]]:
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    status = json.loads(status_path.read_text(encoding="utf-8"))
    if role == "target":
        rows = manifest.get("clusters", [])
        manifest_key = "id"
    elif role == "matched_control":
        rows = manifest.get("clusters", [])
        manifest_key = "id"
    else:
        raise ValueError("role must be target or matched_control")
    completed = {str(row["id"]): row for row in status.get("completed", [])}
    if len(completed) != len(status.get("completed", [])):
        raise ValueError(f"duplicate completed download IDs for {role}")
    if len(rows) != len(completed) or len({str(row[manifest_key]) for row in rows}) != len(rows):
        raise ValueError(f"manifest/status image count or IDs disagree for {role}")
    records = []
    for row in rows:
        ident = str(row[manifest_key])
        saved = completed.get(ident)
        if saved is None or saved.get("file") != row.get("image"):
            raise ValueError(f"manifest/download status mapping mismatch for {role}/{ident}")
        image_path = image_dir / str(saved["file"])
        if not image_path.is_file() or image_path.stat().st_size != int(saved["bytes"]):
            raise ValueError(f"downloaded FITS missing or byte size changed: {image_path}")
        records.append({**row, "frame_id": ident, "path": image_path,
                        "download_sha256": str(saved["sha256"]),
                        "download_bytes": int(saved["bytes"])})
    return records, {"manifest": sha256(manifest_path), "download_status": sha256(status_path)}


def image_wcs(hdul) -> tuple[WCS, tuple[int, int], str]:
    primary = hdul[0].header
    image_index = None
    for index, hdu in enumerate(hdul):
        if str(hdu.header.get("EXTTYPE", "")).strip().upper() == "IMAGE":
            if image_index is not None:
                raise ValueError("duplicate IMAGE HDUs")
            image_index = index
    if image_index is None:
        raise ValueError("FITS has no IMAGE HDU")
    header = hdul[image_index].header.copy()
    if header.get("INHERIT") is True:
        header.extend(primary, unique=True)
    if len(hdul[image_index].shape) != 2:
        raise ValueError("IMAGE HDU is not two-dimensional")
    wcs = WCS(header).celestial
    if not wcs.has_celestial:
        raise ValueError("IMAGE HDU has no celestial WCS")
    filt = str(primary.get("FILTER", "")).strip().lower()
    if filt not in ("i", "i2"):
        raise ValueError(f"unexpected HSC-I family filter: {filt}")
    return wcs, tuple(hdul[image_index].shape), filt


def contained_candidate_rows(wcs: WCS, shape: tuple[int, int],
                             candidates: list[dict]) -> list[dict]:
    ra = [row["ra_deg"] for row in candidates]
    dec = [row["dec_deg"] for row in candidates]
    x, y = wcs.world_to_pixel_values(ra, dec)
    height, width = shape
    rows = []
    for candidate, px, py in zip(candidates, x, y):
        if math.isfinite(float(px)) and math.isfinite(float(py)) and 0 <= px < width and 0 <= py < height:
            rows.append({"candidate_id": candidate["id"],
                         "x_pixel": float(px), "y_pixel": float(py),
                         "z": candidate["z"], "richness": candidate["richness"]})
    return rows


def inspect_role(records: list[dict], candidates: list[dict], role: str) -> tuple[list[dict], dict]:
    matches = []
    filters = {"i": 0, "i2": 0}
    for index, record in enumerate(records, 1):
        with fits.open(record["path"], memmap=True, lazy_load_hdus=True) as hdul:
            wcs, shape, filter_name = image_wcs(hdul)
            matches_here = contained_candidate_rows(wcs, shape, candidates)
        filters[filter_name] += 1
        if matches_here:
            matches.append({
                "role": role,
                "frame_id": record["frame_id"],
                "target_cluster_id": (str(record.get("target_cluster_id"))
                                       if role == "matched_control" else record["frame_id"]),
                "filter": filter_name,
                "image": record["path"].name,
                "download_sha256_at_acquisition": record["download_sha256"],
                "shape": list(shape),
                "n_candidate_centers_inside": len(matches_here),
                "candidates": matches_here,
            })
        if index % 200 == 0:
            print(f"{role}: inspected {index}/{len(records)} FITS headers", flush=True)
    return matches, {"n_images": len(records), "filters": filters,
                     "n_images_with_candidate_centers": len(matches),
                     "n_unique_candidates_inside": len({c["candidate_id"] for m in matches for c in m["candidates"]})}


def run(candidate_table: Path, sample_root: Path) -> dict:
    candidate_table, sample_root = candidate_table.resolve(), sample_root.resolve()
    candidates = read_candidates(candidate_table)
    target_root = sample_root / "download_01"
    control_root = sample_root / "control_download_01"
    target_manifest = target_root / "download_manifest.json"
    target_status = target_root / "download_status.json"
    control_manifest = control_root / "download_manifest.json"
    control_status = control_root / "download_status.json"
    targets, target_hashes = image_records(target_manifest, target_status,
                                           target_root / "images", "target")
    controls, control_hashes = image_records(control_manifest, control_status,
                                              control_root / "images", "matched_control")
    target_matches, target_summary = inspect_role(targets, candidates, "target")
    control_matches, control_summary = inspect_role(controls, candidates, "matched_control")
    return {
        "schema": "astrafek-s23b-candidate-existing-cutout-support-v1",
        "status": "header_wcs_audit_complete",
        "science_measurement": False,
        "hsc_query_submitted": False,
        "data_downloaded": False,
        "image_pixels_read": False,
        "n_candidates": len(candidates),
        "candidate_table_sha256": sha256(candidate_table),
        "target": {**target_summary, "inputs_sha256": target_hashes,
                   "matches": target_matches},
        "matched_controls": {**control_summary, "inputs_sha256": control_hashes,
                             "matches": control_matches},
        "limitations": [
            "Containment means only that the candidate center falls inside an existing image WCS; it does not establish mask, depth, PSF, member, or aperture support.",
            "S23B candidates are a restricted catalog subset; no match does not prove the image contains no other cluster.",
            "Candidate/control overlaps flag possible random-control contamination and are not automatic new science lenses.",
            "Download-time image SHA-256 values are recorded from the existing status files; full FITS files were not rehashed in this header-only audit.",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidate-table", type=Path, default=DEFAULT_CANDIDATE_TABLE)
    parser.add_argument("--sample-root", type=Path, default=SAMPLE_ROOT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(f"refusing to overwrite {args.output}")
    result = run(args.candidate_table, args.sample_root)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n",
                           encoding="utf-8")
    print(json.dumps({"status": result["status"], "n_candidates": result["n_candidates"],
                      "target": {k: result["target"][k] for k in ("n_images", "n_images_with_candidate_centers", "n_unique_candidates_inside")},
                      "matched_controls": {k: result["matched_controls"][k] for k in ("n_images", "n_images_with_candidate_centers", "n_unique_candidates_inside")},
                      "output": str(args.output.resolve())}, indent=2))


if __name__ == "__main__":
    main()
