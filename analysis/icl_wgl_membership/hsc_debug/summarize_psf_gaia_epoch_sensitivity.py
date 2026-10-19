#!/usr/bin/env python3
"""Audit Gaia/HSC candidate-match stability under unknown HSC epochs."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import time

MAS_TO_RAD = math.pi / (180.0 * 3600.0 * 1000.0)
ARCSEC_PER_RAD = 180.0 * 3600.0 / math.pi
DEFAULT_DT_YEARS = [-10, -5, -2, 0, 2, 5, 10]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def propagate_tangent_motion(ra_deg: float, dec_deg: float,
                             pmra_mas_yr: float, pmdec_mas_yr: float,
                             dt_years: float) -> tuple[float, float]:
    """Move a unit direction along the sphere with constant tangent proper motion."""
    alpha, delta = math.radians(ra_deg), math.radians(dec_deg)
    n = (math.cos(delta) * math.cos(alpha),
         math.cos(delta) * math.sin(alpha), math.sin(delta))
    east = (-math.sin(alpha), math.cos(alpha), 0.0)
    north = (-math.sin(delta) * math.cos(alpha),
             -math.sin(delta) * math.sin(alpha), math.cos(delta))
    mu_east = pmra_mas_yr * MAS_TO_RAD
    mu_north = pmdec_mas_yr * MAS_TO_RAD
    speed = math.hypot(mu_east, mu_north)
    if speed == 0.0 or dt_years == 0.0:
        return ra_deg % 360.0, dec_deg
    tangent = tuple(mu_east * east[i] + mu_north * north[i] for i in range(3))
    theta = speed * dt_years
    moved = tuple(math.cos(theta) * n[i] + (math.sin(theta) / speed) * tangent[i]
                  for i in range(3))
    ra = math.degrees(math.atan2(moved[1], moved[0])) % 360.0
    dec = math.degrees(math.asin(max(-1.0, min(1.0, moved[2]))))
    return ra, dec


def separation_arcsec(ra1: float, dec1: float, ra2: float, dec2: float) -> float:
    alpha1, delta1, alpha2, delta2 = map(math.radians, (ra1, dec1, ra2, dec2))
    dra = (alpha2 - alpha1 + math.pi) % (2.0 * math.pi) - math.pi
    ddec = delta2 - delta1
    hav = math.sin(ddec / 2.0) ** 2 + math.cos(delta1) * math.cos(delta2) * math.sin(dra / 2.0) ** 2
    return 2.0 * math.asin(math.sqrt(min(1.0, max(0.0, hav)))) * ARCSEC_PER_RAD


def quantiles(values: list[float]) -> dict[str, float | None]:
    if not values:
        return {"q16": None, "median": None, "q84": None, "maximum": None}
    ordered = sorted(values)

    def percentile(p: float) -> float:
        position = (len(ordered) - 1) * p
        low = int(position)
        high = min(low + 1, len(ordered) - 1)
        return ordered[low] + (ordered[high] - ordered[low]) * (position - low)

    return {"q16": percentile(0.16), "median": percentile(0.5),
            "q84": percentile(0.84), "maximum": ordered[-1]}


def run(candidates_path: Path, matches_path: Path, output: Path,
        dt_years: list[float]) -> dict:
    if output.exists():
        raise FileExistsError(f"refusing to overwrite Gaia epoch sensitivity: {output}")
    candidates = read_tsv(candidates_path)
    if len(candidates) != 480:
        raise ValueError(f"expected the approved 480-candidate set, found {len(candidates)}")
    by_candidate = {row["candidate_id"]: row for row in candidates}
    if len(by_candidate) != len(candidates):
        raise ValueError("duplicate candidate IDs")
    selected_matches = {}
    for row in read_csv(matches_path):
        candidate_id = str(row["candidate_id"])
        if candidate_id not in by_candidate:
            continue
        if candidate_id in selected_matches:
            raise ValueError(f"multiple Gaia matches for candidate {candidate_id}")
        selected_matches[candidate_id] = row
    if set(selected_matches) != set(by_candidate):
        raise ValueError("Gaia match rows do not cover the full selected candidate set")

    scenario_years = sorted(set(float(value) for value in dt_years))
    if 0.0 not in scenario_years or not scenario_years:
        raise ValueError("epoch scenarios must include dt=0")
    output.mkdir(parents=True)
    result = []
    for candidate_id, candidate in by_candidate.items():
        gaia = selected_matches[candidate_id]
        if str(gaia["source_id"]) != str(candidate["gaia_source_id"]):
            raise ValueError(f"Gaia source ID changed for candidate {candidate_id}")
        numeric = ("ra", "dec", "ref_epoch", "pmra", "pmdec", "separation_arcsec")
        values = {name: float(gaia[name]) for name in numeric}
        hsc_ra, hsc_dec = float(candidate["ra_deg"]), float(candidate["dec_deg"])
        recomputed = separation_arcsec(values["ra"], values["dec"], hsc_ra, hsc_dec)
        if abs(recomputed - values["separation_arcsec"]) > 2e-6:
            raise ValueError(f"Gaia/HSC separation mismatch for candidate {candidate_id}")
        row = {
            "candidate_id": candidate_id,
            "cluster_id": candidate["cluster_id"],
            "gaia_source_id": str(gaia["source_id"]),
            "gaia_g_mag": candidate["gaia_g_mag"],
            "gaia_ref_epoch": values["ref_epoch"],
            "hsc_candidate_ra_deg": hsc_ra,
            "hsc_candidate_dec_deg": hsc_dec,
            "gaia_ra_ref_deg": values["ra"],
            "gaia_dec_ref_deg": values["dec"],
            "gaia_hsc_separation_at_ref_arcsec": recomputed,
            "pmra_mas_yr": values["pmra"],
            "pmdec_mas_yr": values["pmdec"],
            "pm_total_mas_yr": math.hypot(values["pmra"], values["pmdec"]),
        }
        for dt in scenario_years:
            moved_ra, moved_dec = propagate_tangent_motion(
                values["ra"], values["dec"], values["pmra"], values["pmdec"], dt)
            row[f"gaia_hsc_separation_dt_{dt:g}yr_arcsec"] = separation_arcsec(
                moved_ra, moved_dec, hsc_ra, hsc_dec)
        result.append(row)

    fields = list(result[0])
    table_path = output / "candidate_epoch_sensitivity.tsv"
    with table_path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(result)
    scenario_summary = {}
    for dt in scenario_years:
        key = f"gaia_hsc_separation_dt_{dt:g}yr_arcsec"
        separations = [float(row[key]) for row in result]
        scenario_summary[f"dt_{dt:g}_years"] = {
            "separation_arcsec": quantiles(separations),
            "within_2_arcsec": sum(value <= 2.0 for value in separations),
        }
    proper_motions = [float(row["pm_total_mas_yr"]) for row in result]
    summary = {
        "status": "complete_scenario_sensitivity",
        "n_candidates": len(result),
        "gaia_reference_epoch_range": [
            min(float(row["gaia_ref_epoch"]) for row in result),
            max(float(row["gaia_ref_epoch"]) for row in result),
        ],
        "delta_epoch_scenarios_years": scenario_years,
        "proper_motion_mas_yr": quantiles(proper_motions),
        "scenario_results": scenario_summary,
        "inputs_sha256": {str(candidates_path.resolve()): sha256(candidates_path),
                           str(matches_path.resolve()): sha256(matches_path)},
        "candidate_table_sha256": sha256(table_path),
        "script_sha256": sha256(Path(__file__)),
        "method": "constant Gaia proper motion propagated as a great-circle tangent vector; HSC candidate position is held fixed; no coadd epoch is assumed",
        "limitations": [
            "Scenario offsets are not actual HSC exposure/coadd epochs.",
            "Perspective acceleration, parallax factors, binary acceleration, and HSC centroid errors are not modeled.",
            "A match-radius pass does not establish stellar identity, isolation, or PSF-wing validity.",
        ],
        "science_measurement": False,
        "psf_approved": False,
        "completed_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    }
    summary_path = output / "summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidates", type=Path, required=True)
    parser.add_argument("--matches", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--dt-years", type=float, nargs="+", default=DEFAULT_DT_YEARS)
    args = parser.parse_args()
    print(json.dumps(run(args.candidates, args.matches, args.output, args.dt_years),
                     indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
