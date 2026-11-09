#!/usr/bin/env python3
"""Audit positional overlap between the existing WL stack and S23B CAMIRA.

This is a local-only provenance and sky-match audit. It reads existing
manifests and count summaries; it does not query HSC CAS, download data, or
read science-image pixels.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
from collections import Counter
from pathlib import Path


DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
WL_SAMPLE = DATA / "sample_manifest.json"
WL_PREFLIGHT = DATA / "wl_shape_preflight_2026-10-18_01"
WL_RANDOM_COUNTS = DATA / "wl_random_shape_count_results_lens_reuse_2026-11-03_01"
WL_RANDOM_CENTERS = DATA / "wl_random_centers_lens_reuse_2026-11-03_01/random_centers.json"
WL_RANDOM_STACK = DATA / "wl_random_shape_stack_lens_reuse_2026-11-04_01/summary.json"
S23B_TARGETS = DATA / "s23b_icl_target_image_plan_2026-11-07_01/target_manifest.json"
S23B_PROFILE_INPUTS = DATA / "s23b_camira_member_outer_icl_relation_2026-11-07_01/cluster_outer_relation_inputs.tsv"
S23B_MEMBER_FILE = DATA / "s23b_camira_member_catalog_2026-11-07_01/camira_s23b_wide_sm_member.dat.gz"
OUT = DATA / "wl_s23b_camira_overlap_audit_2026-11-07_01"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def read_count_csvs(directory: Path) -> dict[str, int]:
    values: dict[str, int] = {}
    for path in sorted(directory.glob("*.csv")):
        with path.open(encoding="utf-8-sig", newline="") as stream:
            first = stream.readline().strip()
            if not first.startswith("# "):
                raise ValueError(f"unexpected count CSV header: {path}")
            fieldnames = first[2:].split(",")
            for row in csv.DictReader(stream, fieldnames=fieldnames):
                key, value = str(row[fieldnames[0]]), int(row[fieldnames[-1]])
                if key in values:
                    raise ValueError(f"duplicate count key {key} in {directory}")
                values[key] = value
    return values


def angular_sep_arcsec(ra1: float, dec1: float, ra2: float, dec2: float) -> float:
    """Great-circle distance, stable at small angles and across RA=0."""
    r1, d1, r2, d2 = map(math.radians, (ra1, dec1, ra2, dec2))
    sd = math.sin((d2 - d1) / 2.0)
    sr = math.sin((r2 - r1) / 2.0)
    hav = sd * sd + math.cos(d1) * math.cos(d2) * sr * sr
    return math.degrees(2.0 * math.asin(math.sqrt(min(1.0, max(0.0, hav))))) * 3600.0


def read_profile_rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def run() -> dict:
    if OUT.exists():
        raise FileExistsError(f"refusing to overwrite existing audit: {OUT}")
    required = (
        WL_SAMPLE,
        WL_PREFLIGHT / "preflight_manifest.json",
        WL_PREFLIGHT / "cas_counts/count_summary.json",
        WL_RANDOM_COUNTS / "count_summary.json",
        WL_RANDOM_CENTERS,
        WL_RANDOM_STACK,
        S23B_TARGETS,
        S23B_PROFILE_INPUTS,
        S23B_MEMBER_FILE,
    )
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise FileNotFoundError("required local audit inputs are missing: " + ", ".join(missing))

    wl_manifest = read_json(WL_SAMPLE)
    preflight = read_json(WL_PREFLIGHT / "preflight_manifest.json")
    random_stack = read_json(WL_RANDOM_STACK)
    random_centers = read_json(WL_RANDOM_CENTERS)
    camira_manifest = read_json(S23B_TARGETS)
    wl_lenses = preflight["lenses"]
    camira_targets = [row for row in camira_manifest["clusters"] if row.get("sample_role") == "target"]
    if len(wl_lenses) != 800 or len(camira_targets) != 1043:
        raise ValueError("unexpected WL/S23B target counts; refusing an unreviewed sample change")
    wl_by_id = {str(row["lens_id"]): row for row in wl_lenses}
    if len(wl_by_id) != len(wl_lenses):
        raise ValueError("WL lens IDs are not unique")
    sample_by_id = {str(row["id"]): row for row in wl_manifest["clusters"]}
    if set(sample_by_id) != set(wl_by_id):
        raise ValueError("WL preflight and original sample manifest ID sets differ")
    for lens_id, lens in wl_by_id.items():
        source = sample_by_id[lens_id]
        if (float(lens["ra_deg"]) != float(source["ra_deg"])
                or float(lens["dec_deg"]) != float(source["dec_deg"])
                or float(lens["z_lens"]) != float(source["z"])):
            raise ValueError(f"WL coordinates/redshift disagree across manifests for lens {lens_id}")

    lens_counts = read_count_csvs(WL_PREFLIGHT / "cas_counts")
    random_counts = read_count_csvs(WL_RANDOM_COUNTS)
    matched_ids = {str(value) for value in random_stack["input"]["matched_lens_ids"]}
    center_to_lens = {str(row["random_center_id"]): str(row["matched_lens_id"])
                      for row in random_centers["centers"]}
    random_support = {
        center_to_lens[center_id]
        for center_id, count in random_counts.items()
        if int(count) > 0 and center_id in center_to_lens
    }
    common_support = matched_ids & {key for key, count in lens_counts.items() if count > 0} & random_support
    common_digest = hashlib.sha256("\n".join(sorted(common_support)).encode("utf-8")).hexdigest()
    expected_common = random_stack["input"]["matched_sample"]["common_lens_ids_sha256"]
    if len(common_support) != 435 or common_digest != expected_common:
        raise ValueError("reconstructed common WL support does not match the saved stack provenance")

    profile_rows = read_profile_rows(S23B_PROFILE_INPUTS)
    profile_catalog_ids = {row["target_id"].removeprefix("S23B_") for row in profile_rows}
    if len(profile_rows) != 427 or len(profile_catalog_ids) != 427:
        raise ValueError("unexpected profile relation target count")
    if not profile_catalog_ids.issubset({str(row["catalog_id"]) for row in camira_targets}):
        raise ValueError("profile relation contains targets absent from the S23B manifest")

    camira_by_catalog_id = {str(row["catalog_id"]): row for row in camira_targets}
    if len(camira_by_catalog_id) != len(camira_targets):
        raise ValueError("S23B catalog IDs are not unique")
    exact_id_overlaps = set(wl_by_id) & set(camira_by_catalog_id)
    exact_profile_id_overlaps = exact_id_overlaps & profile_catalog_ids

    rows = []
    for lens_id, lens in sorted(wl_by_id.items(), key=lambda item: (int(item[0]) if item[0].isdigit() else 10**12, item[0])):
        nearest = min(
            camira_targets,
            key=lambda candidate: angular_sep_arcsec(
                float(lens["ra_deg"]), float(lens["dec_deg"]),
                float(candidate["ra_deg"]), float(candidate["dec_deg"]),
            ),
        )
        nearest_sep = angular_sep_arcsec(
            float(lens["ra_deg"]), float(lens["dec_deg"]),
            float(nearest["ra_deg"]), float(nearest["dec_deg"]),
        )
        exact = camira_by_catalog_id.get(lens_id)
        exact_sep = (angular_sep_arcsec(
            float(lens["ra_deg"]), float(lens["dec_deg"]),
            float(exact["ra_deg"]), float(exact["dec_deg"]),
        ) if exact else None)
        rows.append({
            "wl_lens_id": lens_id,
            "wl_ra_deg": float(lens["ra_deg"]),
            "wl_dec_deg": float(lens["dec_deg"]),
            "wl_z": float(lens["z_lens"]),
            "wl_random_matched_620": lens_id in matched_ids,
            "wl_common_shape_random_support_435": lens_id in common_support,
            "nearest_s23b_catalog_id": str(nearest["catalog_id"]),
            "nearest_s23b_name": str(nearest["name"]),
            "nearest_s23b_ra_deg": float(nearest["ra_deg"]),
            "nearest_s23b_dec_deg": float(nearest["dec_deg"]),
            "nearest_s23b_z": float(nearest["z_cl"]),
            "nearest_sep_arcsec": nearest_sep,
            "nearest_abs_dz": abs(float(lens["z_lens"]) - float(nearest["z_cl"])),
            "nearest_has_member_outer_profile": str(nearest["catalog_id"]) in profile_catalog_ids,
            "numeric_catalog_id_collision": exact is not None,
            "numeric_id_candidate_sep_arcsec": exact_sep,
            "numeric_id_candidate_abs_dz": (abs(float(lens["z_lens"]) - float(exact["z_cl"])) if exact else None),
        })

    OUT.mkdir(parents=True, exist_ok=False)
    table_path = OUT / "wl_lenses_nearest_s23b_target.tsv"
    with table_path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    separations = [row["nearest_sep_arcsec"] for row in rows]
    sorted_separations = sorted(separations)
    def quantile(p: float) -> float:
        position = (len(sorted_separations) - 1) * p
        low, high = math.floor(position), math.ceil(position)
        frac = position - low
        return sorted_separations[low] * (1.0 - frac) + sorted_separations[high] * frac

    thresholds = [1, 3, 10, 30, 60]
    profile_separations = [row["nearest_sep_arcsec"] for row in rows if row["nearest_has_member_outer_profile"]]
    common_rows = [row for row in rows if row["wl_common_shape_random_support_435"]]
    summary = {
        "schema": "astrafek-wl-s23b-camira-positional-overlap-audit-v1",
        "status": "complete_exploratory",
        "science_measurement": False,
        "network_used": False,
        "image_pixels_read": False,
        "source_rows_downloaded": False,
        "inputs": {
            str(path): sha256(path) for path in required
        },
        "sample_counts": {
            "existing_wl_lenses": len(wl_lenses),
            "wl_lenses_with_random_centers": len(matched_ids),
            "wl_lenses_in_verified_common_shape_random_support": len(common_support),
            "s23b_camira_targets": len(camira_targets),
            "s23b_targets_with_saved_member_outer_profile_relation": len(profile_catalog_ids),
            "member_outer_profile_relation_rows_by_filter": dict(Counter(row["filter"] for row in profile_rows)),
        },
        "positional_match": {
            "method": "nearest great-circle S23B CAMIRA target to each existing WL lens center",
            "separation_thresholds_arcsec": {str(threshold): sum(value <= threshold for value in separations) for threshold in thresholds},
            "n_matches_within_60_arcsec_in_wl_random_matched_620": sum(row["nearest_sep_arcsec"] <= 60 for row in rows if row["wl_random_matched_620"]),
            "n_matches_within_60_arcsec_in_common_support_435": sum(row["nearest_sep_arcsec"] <= 60 for row in common_rows),
            "n_matches_within_60_arcsec_among_member_profile_targets": sum(
                min(angular_sep_arcsec(float(w["wl_ra_deg"]), float(w["wl_dec_deg"]),
                                       float(camira_by_catalog_id[cid]["ra_deg"]), float(camira_by_catalog_id[cid]["dec_deg"]))
                    for cid in profile_catalog_ids) <= 60
                for w in rows
            ),
            "nearest_separation_arcsec_quantiles": {
                "min": sorted_separations[0], "q25": quantile(0.25), "median": quantile(0.5),
                "q75": quantile(0.75), "q90": quantile(0.9), "q99": quantile(0.99), "max": sorted_separations[-1],
            },
            "common_support_nearest_separation_min_arcsec": min(row["nearest_sep_arcsec"] for row in common_rows),
            "member_profile_nearest_separation_min_arcsec": min(profile_separations),
        },
        "numeric_id_diagnostic_only": {
            "n_exact_numeric_id_overlaps_all_wl_lenses": len(exact_id_overlaps),
            "n_exact_numeric_id_overlaps_with_member_profile_targets": len(exact_profile_id_overlaps),
            "n_numeric_id_overlaps_within_60_arcsec": sum(
                row["numeric_catalog_id_collision"] and row["numeric_id_candidate_sep_arcsec"] <= 60
                for row in rows
            ),
            "interpretation": "The IDs come from different catalog products. Numeric equality is not a cross-match; use sky position and redshift.",
        },
        "reconstructed_wl_common_support": {
            "n_lenses": len(common_support), "sorted_id_sha256": common_digest,
            "saved_expected_sha256": expected_common,
        },
        "conclusion": "No positional overlap exists within 60 arcsec between the current 800-lens WL sample and the S23B CAMIRA targets; the current WL–ICL–member outer-profile three-way relation cannot be measured by joining these two samples.",
        "next_step": "For a same-sample three-way analysis, either prepare a new HSC S19A shape-count plan for the S23B member-profile targets plus separately validated matched random centers, or verify and use authoritative member positions for the existing WL sample. Any new CAS query or row export requires its own explicit approval.",
        "outputs": {"nearest_target_table": str(table_path)},
    }
    summary_path = OUT / "summary.json"
    summary_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    return summary


if __name__ == "__main__":
    result = run()
    print(json.dumps({
        "status": result["status"],
        "science_measurement": result["science_measurement"],
        "wl_lenses": result["sample_counts"]["existing_wl_lenses"],
        "wl_common_support": result["sample_counts"]["wl_lenses_in_verified_common_shape_random_support"],
        "s23b_targets": result["sample_counts"]["s23b_camira_targets"],
        "nearest_sep_min_arcsec": result["positional_match"]["nearest_separation_arcsec_quantiles"]["min"],
        "n_matches_60_arcsec": result["positional_match"]["separation_thresholds_arcsec"]["60"],
        "numeric_id_overlaps_not_validated": result["numeric_id_diagnostic_only"]["n_exact_numeric_id_overlaps_all_wl_lenses"],
        "output": result["outputs"]["nearest_target_table"],
    }, ensure_ascii=False))
