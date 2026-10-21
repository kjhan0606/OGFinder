#!/usr/bin/env python3
"""Measure radial availability of geometry-clean random coordinates in the partial pilot."""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
import statistics

import numpy as np

from fetch_ancillary_catalogues import sha
from submit_wl_random_control_coordinates import load_count_inputs

RADIAL_OFFSET_EDGES_ARCMIN = (0.0, 5.0, 10.0, 15.0, 20.0, 25.0, 30.0)
OUTER_QUERY_ARCMIN = 30.0
DEFAULT_PILOT = Path(
    "/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/"
    "wl_random_point_geometry_partial_2026-10-20_04"
)


def separation_arcmin(ra1_deg: float, dec1_deg: float,
                      ra2_deg: float, dec2_deg: float) -> float:
    """Great-circle angular separation stable at RA wrap and high declination."""
    values = (ra1_deg, dec1_deg, ra2_deg, dec2_deg)
    if not np.all(np.isfinite(values)) or not (-90 <= dec1_deg <= 90 and -90 <= dec2_deg <= 90):
        raise ValueError("invalid sky coordinate")
    a1, d1, a2, d2 = map(math.radians, values)
    sd = math.sin((d2 - d1) / 2) ** 2
    sr = math.sin((a2 - a1) / 2) ** 2
    hav = sd + math.cos(d1) * math.cos(d2) * sr
    return math.degrees(2 * math.asin(math.sqrt(min(1.0, max(0.0, hav))))) * 60.0


def spherical_annulus_area_arcmin2(inner_arcmin: float, outer_arcmin: float) -> float:
    if not np.isfinite([inner_arcmin, outer_arcmin]).all() or inner_arcmin < 0 or outer_arcmin <= inner_arcmin:
        raise ValueError("annulus radii must be finite and increasing")
    inner = math.radians(inner_arcmin / 60.0)
    outer = math.radians(outer_arcmin / 60.0)
    return 2 * math.pi * (math.cos(inner) - math.cos(outer)) * (180.0 / math.pi * 60.0) ** 2


def write_tsv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f"refusing to write an empty table: {path}")
    with path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t", extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def audit(pilot_dir: Path, output_dir: Path) -> dict:
    pilot_dir, output_dir = Path(pilot_dir).resolve(), Path(output_dir).resolve()
    if output_dir.exists():
        raise FileExistsError(f"refusing to overwrite radial-yield audit: {output_dir}")
    source_summary_path = pilot_dir / "summary.json"
    candidate_path = pilot_dir / "partial_candidate_geometry_audit.tsv"
    per_lens_path = pilot_dir / "partial_per_lens_candidate_counts.tsv"
    summary = json.loads(source_summary_path.read_text(encoding="utf-8"))
    if (summary.get("schema") != "astrafek-wl-random-coordinate-partial-geometry-audit-v1"
            or summary.get("coordinate_export_terminal") is not False
            or summary.get("science_measurement") is not False):
        raise ValueError("input is not the expected incomplete, non-science coordinate pilot")
    for name, path in ((candidate_path.name, candidate_path), (per_lens_path.name, per_lens_path)):
        if sha(path) != summary.get("outputs_sha256", {}).get(name):
            raise ValueError(f"partial geometry output hash mismatch: {name}")

    plan_dir = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_random_point_coordinate_plan_2026-10-19_01")
    plan_path = plan_dir / "random_control_coordinate_plan.json"
    plan = json.loads(plan_path.read_text(encoding="utf-8"))
    if sha(plan_path) != summary.get("coordinate_plan_sha256"):
        raise ValueError("pilot coordinates do not match the pinned random-coordinate plan")
    count_plan = Path(plan["count_plan"])
    _count_manifest, lenses, lens_ids, counts, count_plan_sha = load_count_inputs(count_plan)
    if count_plan_sha != summary.get("count_plan_sha256"):
        raise ValueError("pilot references a changed count-only plan")
    completed = summary.get("completed_batches", [])
    expected_ids = [str(lens_id) for item in completed for lens_id in item.get("lens_ids", [])]
    lens_by_id = {str(row["lens_id"]): row for row in lenses}
    if (len(expected_ids) != summary.get("n_lenses")
            or len(set(expected_ids)) != len(expected_ids)
            or not set(expected_ids) <= set(lens_by_id)):
        raise ValueError("completed pilot batch lens mapping is inconsistent")

    radial_rows = []
    per_lens_rows = []
    raw_by_bin = {i: 0 for i in range(len(RADIAL_OFFSET_EDGES_ARCMIN) - 1)}
    clean_by_bin = {i: 0 for i in range(len(RADIAL_OFFSET_EDGES_ARCMIN) - 1)}
    lens_presence = {i: set() for i in range(len(RADIAL_OFFSET_EDGES_ARCMIN) - 1)}
    lens_clean_presence = {i: set() for i in range(len(RADIAL_OFFSET_EDGES_ARCMIN) - 1)}
    per_lens_counts = {lens_id: {i: [0, 0] for i in raw_by_bin} for lens_id in expected_ids}
    with candidate_path.open(encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        required = {"lens_id", "ra_deg", "dec_deg", "target_3mpc_radius_arcsec",
                    "camira_3mpc_overlap_count", "coordinate_multiplicity", "geometry_candidate"}
        if not reader.fieldnames or not required <= set(reader.fieldnames):
            raise ValueError("partial candidate audit schema changed")
        n_rows = 0
        for row in reader:
            n_rows += 1
            lens_id = str(row["lens_id"])
            if lens_id not in per_lens_counts:
                raise ValueError(f"candidate row belongs to an uncompleted lens: {lens_id}")
            lens = lens_by_id[lens_id]
            sep = separation_arcmin(float(lens["ra_deg"]), float(lens["dec_deg"]),
                                    float(row["ra_deg"]), float(row["dec_deg"]))
            inner = 2.0 * float(row["target_3mpc_radius_arcsec"]) / 60.0
            offset = sep - inner
            if offset < -1e-7 or sep > OUTER_QUERY_ARCMIN + 1e-7:
                raise ValueError("coordinate lies outside the approved 2×3-Mpc–30-arcmin query annulus")
            offset = max(0.0, offset)
            bin_index = next((i for i, (lo, hi) in enumerate(zip(
                RADIAL_OFFSET_EDGES_ARCMIN[:-1], RADIAL_OFFSET_EDGES_ARCMIN[1:]))
                if lo <= offset < hi), None)
            if bin_index is None:
                raise ValueError("coordinate falls outside radial-yield bins")
            clean = row["geometry_candidate"].lower() == "true"
            if clean and (int(row["camira_3mpc_overlap_count"]) != 0
                          or int(row["coordinate_multiplicity"]) != 1):
                raise ValueError("geometry candidate contradicts its overlap/duplicate fields")
            raw_by_bin[bin_index] += 1
            clean_by_bin[bin_index] += int(clean)
            lens_presence[bin_index].add(lens_id)
            if clean:
                lens_clean_presence[bin_index].add(lens_id)
            per_lens_counts[lens_id][bin_index][0] += 1
            per_lens_counts[lens_id][bin_index][1] += int(clean)
    if n_rows != summary.get("n_coordinates"):
        raise ValueError("geometry audit row count differs from its summary")
    if sum(raw_by_bin.values()) != n_rows or sum(clean_by_bin.values()) != summary.get("n_geometry_candidates"):
        raise ValueError("radial bins do not reproduce pilot totals")

    for i, (lo, hi) in enumerate(zip(RADIAL_OFFSET_EDGES_ARCMIN[:-1], RADIAL_OFFSET_EDGES_ARCMIN[1:])):
        area_rows = []
        lens_area_by_id = {}
        for lens_id in expected_ids:
            inner3mpc_arcmin = float(lens_by_id[lens_id]["radius_arcsec"]) / 60.0
            shell_inner = 2.0 * inner3mpc_arcmin + lo
            shell_outer = min(2.0 * inner3mpc_arcmin + hi, OUTER_QUERY_ARCMIN)
            lens_area = (spherical_annulus_area_arcmin2(shell_inner, shell_outer)
                         if shell_outer > shell_inner else 0.0)
            lens_area_by_id[lens_id] = lens_area
            area_rows.append(lens_area)
        available_ids = [lens_id for lens_id in expected_ids if lens_area_by_id[lens_id] > 0]
        clean_counts = np.asarray([per_lens_counts[lens_id][i][1] for lens_id in available_ids], dtype=float)
        area = float(sum(area_rows))
        raw = raw_by_bin[i]
        clean = clean_by_bin[i]
        radial_rows.append({
            "offset_inner_arcmin": lo, "offset_outer_arcmin": hi,
            "raw_coordinate_rows": raw, "geometry_clean_rows": clean,
            "geometry_survival_fraction": clean / raw if raw else None,
            "lenses_with_raw_rows": len(lens_presence[i]),
            "lenses_with_geometry_clean_rows": len(lens_clean_presence[i]),
            "lenses_with_annulus_area": len(available_ids),
            "lenses_without_annulus_area": len(expected_ids) - len(available_ids),
            "median_clean_rows_per_lens": float(np.median(clean_counts)) if len(clean_counts) else None,
            "q10_clean_rows_per_lens": float(np.quantile(clean_counts, 0.10)) if len(clean_counts) else None,
            "q90_clean_rows_per_lens": float(np.quantile(clean_counts, 0.90)) if len(clean_counts) else None,
            "lenses_with_zero_clean_rows": int(np.count_nonzero(clean_counts == 0)) if len(clean_counts) else None,
            "ideal_spherical_annulus_area_arcmin2": area,
            "geometry_clean_density_per_ideal_area_arcmin2": clean / area if area > 0 else None,
            "actual_HSC_footprint_correction": "not_available_from_candidate_points_only",
        })
        for index, lens_id in enumerate(expected_ids):
            n_raw, n_clean = per_lens_counts[lens_id][i]
            inner3mpc_arcmin = float(lens_by_id[lens_id]["radius_arcsec"]) / 60.0
            shell_inner = 2.0 * inner3mpc_arcmin + lo
            shell_outer = min(2.0 * inner3mpc_arcmin + hi, OUTER_QUERY_ARCMIN)
            lens_area = (spherical_annulus_area_arcmin2(shell_inner, shell_outer)
                         if shell_outer > shell_inner else 0.0)
            per_lens_rows.append({
                "lens_id": lens_id, "offset_inner_arcmin": lo, "offset_outer_arcmin": hi,
                "raw_coordinate_rows": n_raw, "geometry_clean_rows": n_clean,
                "ideal_spherical_annulus_area_arcmin2": lens_area,
                "query_annulus_available": lens_area > 0,
                "geometry_clean_density_per_ideal_area_arcmin2": n_clean / lens_area if lens_area > 0 else None,
            })

    expansion_area_ratios = []
    for lens_id in expected_ids:
        inner3mpc_arcmin = float(lens_by_id[lens_id]["radius_arcsec"]) / 60.0
        existing_area = spherical_annulus_area_arcmin2(2.0 * inner3mpc_arcmin, OUTER_QUERY_ARCMIN)
        expanded_area = spherical_annulus_area_arcmin2(OUTER_QUERY_ARCMIN, 60.0)
        expansion_area_ratios.append(expanded_area / existing_area)
    output_dir.mkdir(parents=True, exist_ok=False)
    write_tsv(output_dir / "radial_yield.tsv", radial_rows)
    write_tsv(output_dir / "per_lens_radial_yield.tsv", per_lens_rows)
    result = {
        "schema": "astrafek-wl-random-coordinate-radial-yield-v1",
        "status": "complete_partial_radial_yield_audit",
        "coordinate_export_terminal": False, "n_completed_batches": len(completed),
        "n_lenses": len(expected_ids), "n_coordinates": n_rows,
        "n_geometry_candidates": int(sum(clean_by_bin.values())),
        "radial_offset_definition": "great-circle separation from lens center minus twice the angular proper-3-Mpc radius",
        "radial_offset_edges_arcmin": list(RADIAL_OFFSET_EDGES_ARCMIN),
        "outer_query_radius_arcmin": OUTER_QUERY_ARCMIN,
        "available_query_annulus_width_arcmin_by_lens": {
            "minimum": min(OUTER_QUERY_ARCMIN - 2.0 * float(lens_by_id[lens_id]["radius_arcsec"]) / 60.0
                           for lens_id in expected_ids),
            "median": statistics.median(OUTER_QUERY_ARCMIN - 2.0 * float(lens_by_id[lens_id]["radius_arcsec"]) / 60.0
                                        for lens_id in expected_ids),
            "maximum": max(OUTER_QUERY_ARCMIN - 2.0 * float(lens_by_id[lens_id]["radius_arcsec"]) / 60.0
                           for lens_id in expected_ids),
        },
        "ideal_30_60_arcmin_area_divided_by_existing_2R3Mpc_30_arcmin_area": {
            "minimum": min(expansion_area_ratios), "median": statistics.median(expansion_area_ratios),
            "maximum": max(expansion_area_ratios),
        },
        "camira_catalog_sha256": summary.get("camira_catalog_sha256"),
        "input_sha256": {str(path.resolve()): sha(path) for path in (
            source_summary_path, candidate_path, per_lens_path, plan_path,
            count_plan / "random_control_count_plan.json",
        )},
        "outputs_sha256": {
            "radial_yield.tsv": sha(output_dir / "radial_yield.tsv"),
            "per_lens_radial_yield.tsv": sha(output_dir / "per_lens_radial_yield.tsv"),
        },
        "science_ready": False, "science_measurement": False,
        "limitations": [
            "Only the first two completed coordinate batches (50 plan-order lenses) are represented; this is not a random or representative sample of the 800 lenses.",
            "The ideal spherical annulus area is not corrected for actual HSC footprint, depth, mask, or source-selection completeness.",
            "Geometry-clean means only no 3-proper-Mpc CAMIRA aperture overlap and no exact coordinate duplicate in this pilot.",
            "No controls were selected; no source shapes, photo-z values, or images were requested.",
        ],
    }
    (output_dir / "summary.json").write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pilot", type=Path, default=DEFAULT_PILOT)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = audit(args.pilot, args.output)
    print(json.dumps({key: result[key] for key in (
        "status", "n_completed_batches", "n_lenses", "n_coordinates", "n_geometry_candidates",
        "science_ready", "science_measurement",
    )}, indent=2))


if __name__ == "__main__":
    main()
