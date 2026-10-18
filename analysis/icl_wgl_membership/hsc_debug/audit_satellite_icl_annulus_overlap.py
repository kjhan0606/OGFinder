#!/usr/bin/env python3
"""Quantify geometric overlap of satellite outskirts with the ICL annulus.

This is a support-geometry audit only. It does not model source masks, pixel
noise, PSF wings, or flux covariance and never upgrades the exploratory stacks
to science measurements.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path


DEFAULT_SAMPLE = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
DEFAULT_ICL = DEFAULT_SAMPLE / "prestack_flux_profiles_2026-10-18_02/pair_metrics_k2p5_c40_r50_300.tsv"
DEFAULT_SATELLITE = DEFAULT_SAMPLE / "satellite_outer_profile_perbin_full_2026-10-18_01/satellite_profiles.tsv"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def disk_intersection_area(r1: float, r2: float, separation: float) -> float:
    """Area shared by two disks, including containment and disjoint cases."""
    if min(r1, r2, separation) < 0 or r1 == 0 or r2 == 0:
        raise ValueError("radii and separation must be non-negative; radii must be positive")
    if separation >= r1 + r2:
        return 0.0
    if separation <= abs(r1 - r2):
        return math.pi * min(r1, r2) ** 2
    c1 = (separation * separation + r1 * r1 - r2 * r2) / (2 * separation * r1)
    c2 = (separation * separation + r2 * r2 - r1 * r1) / (2 * separation * r2)
    radicand = max(0.0, (-separation + r1 + r2) * (separation + r1 - r2)
                   * (separation - r1 + r2) * (separation + r1 + r2))
    return r1 * r1 * math.acos(max(-1.0, min(1.0, c1))) \
        + r2 * r2 * math.acos(max(-1.0, min(1.0, c2))) - 0.5 * math.sqrt(radicand)


def annulus_intersection_area(center_separation: float,
                              icl_inner: float = 200.0,
                              icl_outer: float = 300.0,
                              satellite_inner: float = 50.0,
                              satellite_outer: float = 75.0) -> float:
    """Exact flat-sky overlap area of two circular annuli (kpc^2)."""
    if not (0 <= icl_inner < icl_outer and 0 <= satellite_inner < satellite_outer):
        raise ValueError("annulus radii must be ordered and non-negative")
    d = center_separation
    overlap = (
        disk_intersection_area(icl_outer, satellite_outer, d)
        - disk_intersection_area(icl_inner, satellite_outer, d)
        - disk_intersection_area(icl_outer, satellite_inner, d)
        + disk_intersection_area(icl_inner, satellite_inner, d)
    )
    # Inclusion-exclusion may produce tiny negative roundoff at tangencies.
    return max(0.0, overlap)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def analyze(icl_path: Path, satellite_path: Path,
            icl_inner: float = 200.0, icl_outer: float = 300.0,
            sat_inner: float = 50.0, sat_outer: float = 75.0) -> dict:
    icl_rows = read_tsv(icl_path)
    satellite_rows = read_tsv(satellite_path)
    required_icl = {"id", "filter", "coverage_200_300"}
    required_sat = {"cluster_id", "filter", "member_index", "cluster_radius_kpc",
                    "radial_bin_kpc", "target_coverage", "control_coverage"}
    if not icl_rows or not required_icl.issubset(icl_rows[0]):
        raise ValueError(f"ICL table lacks columns: {sorted(required_icl)}")
    if not satellite_rows or not required_sat.issubset(satellite_rows[0]):
        raise ValueError(f"satellite table lacks columns: {sorted(required_sat)}")

    icl_pairs = set()
    for row in icl_rows:
        if float(row["coverage_200_300"]) >= 0.40:
            icl_pairs.add((row["id"], row["filter"]))

    unique = {}
    for row in satellite_rows:
        if row["radial_bin_kpc"] != f"{sat_inner:g}-{sat_outer:g}":
            continue
        key = (row["cluster_id"], row["filter"], row["member_index"])
        if key in unique:
            raise ValueError(f"duplicate satellite profile row: {key}")
        unique[key] = row

    records = []
    for (cluster_id, filt, member_index), row in sorted(unique.items()):
        if (cluster_id, filt) not in icl_pairs:
            continue
        radius = float(row["cluster_radius_kpc"])
        overlap = annulus_intersection_area(radius, icl_inner, icl_outer,
                                            sat_inner, sat_outer)
        sat_area = math.pi * (sat_outer ** 2 - sat_inner ** 2)
        icl_area = math.pi * (icl_outer ** 2 - icl_inner ** 2)
        records.append({
            "cluster_id": cluster_id,
            "filter": filt,
            "member_index": member_index,
            "cluster_radius_kpc": radius,
            "overlap_area_kpc2": overlap,
            "satellite_annulus_overlap_fraction": overlap / sat_area,
            "icl_annulus_overlap_fraction": overlap / icl_area,
            "target_coverage": float(row["target_coverage"]),
            "control_coverage": float(row["control_coverage"]),
        })

    by_filter = {}
    for filt in sorted({r["filter"] for r in records}):
        subset = [r for r in records if r["filter"] == filt]
        fracs = sorted(r["satellite_annulus_overlap_fraction"] for r in subset)
        radii = [r["cluster_radius_kpc"] for r in subset]
        def quantile(values: list[float], q: float) -> float:
            if len(values) == 1:
                return values[0]
            pos = (len(values) - 1) * q
            lo = int(math.floor(pos)); hi = int(math.ceil(pos))
            return values[lo] + (values[hi] - values[lo]) * (pos - lo)
        by_filter[filt] = {
            "n_satellite_profiles": len(subset),
            "n_clusters": len({r["cluster_id"] for r in subset}),
            "n_with_geometric_overlap": sum(x > 0 for x in fracs),
            "n_with_satellite_area_overlap_ge_10pct": sum(x >= 0.10 for x in fracs),
            "median_cluster_radius_kpc": quantile(sorted(radii), 0.5),
            "satellite_annulus_overlap_fraction_q50_q90": [
                quantile(fracs, 0.5), quantile(fracs, 0.9)],
            "median_icl_annulus_overlap_fraction": quantile(
                sorted(r["icl_annulus_overlap_fraction"] for r in subset), 0.5),
        }
    return {
        "schema": "icl-satellite-annulus-geometric-overlap-v1",
        "status": "complete",
        "science_measurement": False,
        "flux_covariance_available": False,
        "geometry": {
            "icl_annulus_kpc": [icl_inner, icl_outer],
            "satellite_outer_annulus_kpc": [sat_inner, sat_outer],
            "overlap_definition": "exact flat-sky intersection of circular annuli using each member's cluster-centric projected radius",
        },
        "limitations": [
            "This is an unmasked support-geometry upper bound; it does not measure shared valid pixels or flux covariance.",
            "ICL and satellite estimators use different source masks and controls; no mask overlap is inferred.",
            "PSF wings, sky transfer, correlated noise, and flux allocation remain untested.",
            "Kron 2.5 and 40% ICL coverage are exploratory settings, not science-approved cuts.",
        ],
        "n_icl_cluster_filter_pairs": len(icl_pairs),
        "n_satellite_outer_profiles_total": len(unique),
        "n_matched_satellite_profiles": len(records),
        "inputs": {
            "icl_pair_metrics": {"path": str(icl_path), "sha256": sha256(icl_path)},
            "satellite_profiles": {"path": str(satellite_path), "sha256": sha256(satellite_path)},
        },
        "by_filter": by_filter,
        "records": records,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--icl", type=Path, default=DEFAULT_ICL)
    parser.add_argument("--satellite", type=Path, default=DEFAULT_SATELLITE)
    parser.add_argument("--output", type=Path, required=True,
                        help="new JSON output path; existing files are never overwritten")
    args = parser.parse_args()
    if args.output.exists():
        parser.error(f"refusing to overwrite existing output: {args.output}")
    result = analyze(args.icl, args.satellite)
    result["script_sha256"] = sha256(Path(__file__))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n",
                           encoding="utf-8")
    print(json.dumps({k: result[k] for k in
                      ("status", "n_icl_cluster_filter_pairs", "n_satellite_outer_profiles_total",
                       "n_matched_satellite_profiles", "by_filter", "science_measurement",
                       "flux_covariance_available")}, indent=2))


if __name__ == "__main__":
    main()
