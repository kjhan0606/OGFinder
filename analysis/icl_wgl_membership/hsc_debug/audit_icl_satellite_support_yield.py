#!/usr/bin/env python3
"""Quantify the ICL/satellite radial support intersection in the 800-parent set.

This is an offline support audit. It reads saved Astrafek tables only; it does
not request HSC rows, inspect image pixels, or relax any measurement cut.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path


DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
ICL_PAIRS = DATA / "icl_annulus_support_2026-10-18_01/pair_metrics_k2p5_c40_r50_300.tsv"
ICL_AUDIT = DATA / "icl_annulus_support_2026-10-18_01/summary.json"
DEFAULT_SATELLITE_DIR = DATA / "satellite_outer_profile_perbin_full_2026-10-18_01"
BINS = ("5-10", "10-15", "15-20", "20-30", "30-40", "40-50", "50-75")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty table: {path}")
    with path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t",
                                extrasaction="raise", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def run(output: Path, satellite_dir: Path) -> dict:
    satellite_clusters = satellite_dir / "cluster_profiles.tsv"
    satellite_candidates = satellite_dir / "candidate_diagnostics.tsv"
    satellite_audit_path = satellite_dir / "summary.json"
    inputs = (ICL_PAIRS, satellite_clusters, satellite_candidates, ICL_AUDIT, satellite_audit_path)
    for path in inputs:
        if not path.is_file():
            raise FileNotFoundError(path)
    if output.exists():
        raise FileExistsError(f"refusing to overwrite {output}")

    icl_audit = json.loads(ICL_AUDIT.read_text(encoding="utf-8"))
    satellite_audit = json.loads(satellite_audit_path.read_text(encoding="utf-8"))
    if (icl_audit.get("schema") != "hsc-icl-per-annulus-support-v1"
            or icl_audit.get("status") != "complete"
            or icl_audit.get("n_200_300_pairs") != 224
            or satellite_audit.get("schema") != "hsc-satellite-outer-profile-per-annulus-v1"
            or satellite_audit.get("status") != "complete"
            or satellite_audit.get("science_measurement") is not False):
        raise ValueError("input support summaries are incomplete or use unexpected schemas")

    icl_rows = read_tsv(ICL_PAIRS)
    satellite_rows = read_tsv(satellite_clusters)
    candidate_rows = read_tsv(satellite_candidates)
    icl_by_key = {(r["id"], r["filter"]): r for r in icl_rows}
    sat_by_key = {(r["cluster_id"], r["filter"], r["radial_bin_kpc"]): r for r in satellite_rows}
    if len(icl_rows) != 224 or len(icl_by_key) != len(icl_rows):
        raise ValueError("expected 224 unique ICL-supported cluster/filter pairs")
    if len(sat_by_key) != len(satellite_rows) or not set(r["radial_bin_kpc"] for r in satellite_rows) <= set(BINS):
        raise ValueError("duplicate satellite profile key or unknown radial bin")
    if len(candidate_rows) != satellite_audit.get("n_candidate_diagnostics"):
        raise ValueError("satellite candidate diagnostic rows differ from their saved summary")

    radial_rows = []
    by_filter = defaultdict(dict)
    for radial_bin in BINS:
        for band in ("i", "i2"):
            satellite_keys = {(ident, flt) for ident, flt, rb in sat_by_key if flt == band and rb == radial_bin}
            icl_keys = {key for key in icl_by_key if key[1] == band}
            common = satellite_keys & icl_keys
            center_eligible = sum(icl_by_key[key]["center_masked_k2p5"] == "False" for key in common)
            row = {
                "satellite_radial_bin_kpc": radial_bin,
                "filter": band,
                "satellite_profile_cluster_filter_pairs": len(satellite_keys),
                "icl_200_300_cluster_filter_pairs": len(icl_keys),
                "common_cluster_filter_pairs": len(common),
                "common_center_eligible_pairs": center_eligible,
                "fraction_of_icl_support_with_satellite_profile": len(common) / len(icl_keys) if icl_keys else None,
                "science_measurement": False,
            }
            radial_rows.append(row)
            by_filter[band][radial_bin] = row

    status_counts = defaultdict(Counter)
    candidate_groups = Counter()
    for row in candidate_rows:
        status = row["status"]
        band = row["filter"]
        cohort = row["cohort"]
        status_counts[(band, cohort)][status] += 1
        candidate_groups[(band, cohort)] += 1
    allowed = {"no_usable_annulus", "retained", "insufficient_bin_controls", "invalid_target_shape"}
    if any(set(counts) - allowed for counts in status_counts.values()):
        raise ValueError("unexpected satellite-candidate status in source diagnostics")
    status_rows = []
    for (band, cohort), counts in sorted(status_counts.items()):
        total = candidate_groups[(band, cohort)]
        for status in sorted(allowed):
            count = counts.get(status, 0)
            status_rows.append({
                "filter": band, "cohort": cohort, "candidate_status": status,
                "candidate_count": count, "candidates_in_filter_cohort": total,
                "fraction": count / total if total else None,
                "science_measurement": False,
            })

    output.mkdir(parents=True, exist_ok=False)
    radial_path = output / "radial_support_intersection.tsv"
    status_path = output / "satellite_candidate_status.tsv"
    write_tsv(radial_path, radial_rows)
    write_tsv(status_path, status_rows)
    summary = {
        "schema": "astrafek-icl-satellite-support-yield-v1",
        "status": "complete_offline_diagnostic",
        "science_measurement": False,
        "hsc_network_used": False,
        "image_pixels_read": False,
        "data_requested": False,
        "parent_sample_size": 800,
        "icl_200_300_supported_pairs": len(icl_rows),
        "satellite_candidate_diagnostics": len(candidate_rows),
        "satellite_candidate_status_counts": dict(Counter(r["status"] for r in candidate_rows)),
        "satellite_cluster_profiles_by_radial_bin_filter": {
            band: {rb: sum(1 for r in satellite_rows if r["filter"] == band and r["radial_bin_kpc"] == rb)
                   for rb in BINS}
            for band in ("i", "i2")
        },
        "support_intersection_by_radial_bin_filter": by_filter,
        "interpretation": "support-yield audit only; varying annuli and filters change the selected sample; no association estimator or measurement threshold was changed",
        "input_sha256": {str(path): sha256(path) for path in inputs},
        "output_sha256": {str(radial_path): sha256(radial_path), str(status_path): sha256(status_path)},
    }
    summary_path = output / "summary.json"
    summary_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    summary["summary_sha256"] = sha256(summary_path)
    print(json.dumps({
        "status": summary["status"], "science_measurement": False,
        "candidate_status_counts": summary["satellite_candidate_status_counts"],
        "intersection_at_50_75": {band: by_filter[band]["50-75"] for band in ("i", "i2")},
        "output": str(output.resolve()),
    }, ensure_ascii=False, indent=2))
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--satellite-dir", type=Path, default=DEFAULT_SATELLITE_DIR)
    args = parser.parse_args()
    run(args.output.resolve(), args.satellite_dir.resolve())


if __name__ == "__main__":
    main()
