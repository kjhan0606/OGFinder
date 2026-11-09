#!/usr/bin/env python3
"""Audit CAMIRA member positions available for forced satellite profiles.

This offline inventory tests whether the existing image/preparation sample can
extend beyond SEP-detected CAMIRA members. It reads catalogues and manifests,
not science image pixels, and performs no HSC queries or downloads.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
from astropy import units as u
from astropy.cosmology import FlatLambdaCDM
from astropy.coordinates import SkyCoord

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import measure_satellite_outer_profiles as core

DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
PREP = DATA / "satellite_detection_psf_2026-10-18_01"
MANIFEST = DATA / "sample_manifest.json"
PERBIN = DATA / "satellite_outer_profile_perbin_full_2026-10-18_01"
DEFAULT_OUTPUT = DATA / "satellite_forced_position_support_2026-11-06_03"
RMIN_KPC = 50.0
RMAX_KPC = 500.0
FORCED_MIN_SEP_ARCSEC = 3.0


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def write_tsv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty table: {path}")
    with path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t",
                                lineterminator="\n", extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def run(output: Path) -> dict:
    if output.exists():
        raise FileExistsError(f"refusing to overwrite {output}")
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    cosmology = FlatLambdaCDM(H0=float(manifest.get("h0", 70.0)),
                              Om0=float(manifest.get("omega_m", 0.3)), Tcmb0=0 * u.K)
    perbin_summary_path = PERBIN / "summary.json"
    perbin_summary = json.loads(perbin_summary_path.read_text(encoding="utf-8"))
    diagnostics_path = PERBIN / "candidate_diagnostics.tsv"
    diagnostics = core.read_tsv(diagnostics_path)
    if len(manifest.get("clusters", [])) != 800 or perbin_summary.get("n_frames_completed") != 800:
        raise ValueError("expected the complete existing 800-target image sample")
    if perbin_summary.get("n_candidate_diagnostics") != len(diagnostics):
        raise ValueError("per-bin summary does not match candidate diagnostic table")
    expected_matched = {(str(row["cluster_id"]), str(row["member_index"]))
                        for row in diagnostics}
    if len(expected_matched) != len(diagnostics):
        raise ValueError("duplicate cluster/member keys in existing candidate diagnostics")

    frame_summary = {str(row["cluster_id"]): row
                     for row in perbin_summary.get("frame_records", [])}
    if len(frame_summary) != 800:
        raise ValueError(f"expected 800 existing frame records; got {len(frame_summary)}")
    output_rows = []
    counts = {}
    input_hashes = {str(MANIFEST): sha256(MANIFEST),
                    str(perbin_summary_path): sha256(perbin_summary_path),
                    str(diagnostics_path): sha256(diagnostics_path)}
    for frame in manifest["clusters"]:
        ident = str(frame["id"])
        prep_dir = PREP / f"cluster_{ident}"
        member_path = prep_dir / "members.tsv"
        detection_path = prep_dir / "detections.tsv"
        preparation_path = prep_dir / "preparation.json"
        for path in (member_path, detection_path, preparation_path):
            if not path.is_file():
                raise FileNotFoundError(path)
        input_hashes[str(member_path)] = sha256(member_path)
        input_hashes[str(detection_path)] = sha256(detection_path)
        input_hashes[str(preparation_path)] = sha256(preparation_path)
        preparation = core.read_json(preparation_path)
        if preparation.get("science_measurement") is not False:
            raise ValueError(f"unexpected preparation gate for cluster {ident}")
        detections = core.parse_detections(core.read_tsv(detection_path))
        members = core.read_tsv(member_path)
        unique_matches, _ = core.match_members(members, detections, 0.168)
        matches_by_index = {int(row["member_index"]): row for row in unique_matches}
        matched_ids = set(matches_by_index)
        in_frame = [row for row in members if row["IN_FRAME"] == "True"]
        if not in_frame:
            continue
        member_sky = SkyCoord(
            ra=np.asarray([float(row["RA"]) for row in in_frame]) * u.deg,
            dec=np.asarray([float(row["DEC"]) for row in in_frame]) * u.deg,
        )
        cluster_sky = SkyCoord(ra=float(frame["ra_deg"]) * u.deg,
                               dec=float(frame["dec_deg"]) * u.deg)
        radii_kpc = (member_sky.separation(cluster_sky).to(u.arcmin)
                     * cosmology.kpc_proper_per_arcmin(float(frame["z"]))).value
        diagnostic_radii = {
            str(row["member_index"]): float(row["cluster_radius_kpc"])
            for row in diagnostics if str(row["cluster_id"]) == ident
        }
        frame_filter = frame_summary.get(ident, {}).get("filter")
        if frame_filter not in ("i", "i2"):
            raise ValueError(f"missing/invalid filter in existing frame audit for {ident}")
        per_frame = {"n_positions_in_frame": len(in_frame), "n_outer_positions": 0,
                     "n_unique_sep_matches": 0, "n_forced_clean_candidates": 0,
                     "n_near_or_ambiguous_unmatched": 0, "n_missing_sep_unmatched": 0,
                     "n_match_candidates_with_any_profile": 0}
        for member_row, member_radius in zip(in_frame, radii_kpc):
            member_index = int(member_row["INDEX"])
            if member_index in matches_by_index:
                radius = diagnostic_radii.get(str(member_index), float(member_radius))
            else:
                radius = float(member_radius)
            if not (RMIN_KPC <= radius < RMAX_KPC):
                continue
            separation = float(member_row["DETECTION_SEPARATION_ARCSEC"])
            matched = member_index in matched_ids
            per_frame["n_outer_positions"] += 1
            if matched:
                per_frame["n_unique_sep_matches"] += 1
                if (ident, str(member_index)) in expected_matched:
                    per_frame["n_match_candidates_with_any_profile"] += 1
                category = "unique_sep_match_le_1arcsec"
            elif np.isfinite(separation) and separation > FORCED_MIN_SEP_ARCSEC:
                per_frame["n_forced_clean_candidates"] += 1
                category = "no_sep_detection_within_3arcsec"
            elif np.isfinite(separation):
                per_frame["n_near_or_ambiguous_unmatched"] += 1
                category = "unmatched_sep_detection_within_3arcsec"
            else:
                per_frame["n_missing_sep_unmatched"] += 1
                category = "missing_nearest_sep"
            output_rows.append({
                "cluster_id": ident,
                "filter": frame_filter,
                "redshift": float(frame["z"]),
                "richness": float(frame["richness"]),
                "member_index": member_index,
                "cluster_radius_kpc": float(radius),
                "nearest_sep_detection_arcsec": separation,
                "unique_sep_match": matched,
                "existing_profile_candidate": (ident, str(member_index)) in expected_matched,
                "category": category,
                "science_measurement": False,
            })
        counts[ident] = per_frame

    output_rows.sort(key=lambda row: (row["filter"], int(row["cluster_id"]), row["member_index"]))
    if len(expected_matched) != 9468:
        raise ValueError(f"unexpected existing matched candidate count: {len(expected_matched)}")
    if sum(row["existing_profile_candidate"] for row in output_rows) != 9468:
        raise ValueError("existing matched candidates do not reproduce all per-bin diagnostics")
    totals = {key: sum(value[key] for value in counts.values())
              for key in next(iter(counts.values()))}
    by_filter = {}
    for band in ("i", "i2"):
        selected = [row for row in output_rows if row["filter"] == band]
        by_filter[band] = {
            "n_positions_50_500kpc": len(selected),
            "n_unique_sep_matches_by_catalog_radius": sum(row["category"] == "unique_sep_match_le_1arcsec" for row in selected),
            "n_forced_clean_candidates": sum(row["category"] == "no_sep_detection_within_3arcsec" for row in selected),
            "n_near_or_ambiguous_unmatched": sum(row["category"] == "unmatched_sep_detection_within_3arcsec" for row in selected),
            "n_missing_sep_unmatched": sum(row["category"] == "missing_nearest_sep" for row in selected),
        }
    result = {
        "schema": "hsc-camira-forced-satellite-position-support-v1",
        "status": "complete_exploratory",
        "science_measurement": False,
        "sample": "existing 800 HSC-I target cutouts; catalog/preparation metadata only",
        "radial_range_proper_kpc": [RMIN_KPC, RMAX_KPC],
        "forced_candidate_definition": f"IN_FRAME CAMIRA position with no unique SEP match within 1 arcsec and nearest saved detection farther than {FORCED_MIN_SEP_ARCSEC:g} arcsec",
        "n_target_frames": 800,
        "n_outer_member_positions": totals["n_outer_positions"],
        "n_unique_sep_matches_by_catalog_radius": totals["n_unique_sep_matches"],
        "n_existing_unique_matches": len(expected_matched),
        "n_existing_matched_candidates_with_profile_attempt": totals["n_match_candidates_with_any_profile"],
        "n_forced_clean_candidates": totals["n_forced_clean_candidates"],
        "n_near_or_ambiguous_unmatched": totals["n_near_or_ambiguous_unmatched"],
        "n_missing_sep_unmatched": totals["n_missing_sep_unmatched"],
        "by_filter": by_filter,
        "limitations": [
            "This audit counts candidates; it does not measure forced profiles or demonstrate usable annular coverage.",
            "A CAMIRA member position without a nearby SEP detection can still be a false/uncertain member or a blend.",
            "The current PSF, BCG center, source masks, sky transfer, and coadd-noise gates remain unresolved.",
            "The radial range and separation categories are exploratory and do not certify a science sample.",
        ],
        "input_hashes": {
            "sample_manifest": sha256(MANIFEST),
            "perbin_summary": sha256(perbin_summary_path),
            "candidate_diagnostics": sha256(diagnostics_path),
            "all_member_detection_preparation_files": hashlib.sha256(
                "\n".join(f"{path}:{digest}" for path, digest in sorted(input_hashes.items())).encode()
            ).hexdigest(),
        },
    }
    output.mkdir(parents=True, exist_ok=False)
    table_path = output / "member_position_support.tsv"
    write_tsv(table_path, output_rows)
    result["output_sha256"] = {"member_position_support.tsv": sha256(table_path)}
    result["script_sha256"] = sha256(Path(__file__).resolve())
    (output / "summary.json").write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n",
                                         encoding="utf-8")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    result = run(args.output)
    print(json.dumps({key: result[key] for key in (
        "status", "output", "n_outer_member_positions", "n_unique_sep_matches_by_catalog_radius",
        "n_existing_unique_matches",
        "n_forced_clean_candidates", "n_near_or_ambiguous_unmatched", "by_filter",
        "science_measurement") if key in result}, indent=2))


if __name__ == "__main__":
    main()
