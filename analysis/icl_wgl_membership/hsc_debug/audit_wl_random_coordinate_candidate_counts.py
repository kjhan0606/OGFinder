#!/usr/bin/env python3
"""Audit approved aggregate HSC random-point counts for all 800 lens centers.

This reads the already returned count-only CSV batches and the prepared
coordinate-only SQL plan. It makes no CAS request and reads no coordinates,
source rows, shapes, photo-z values, or images.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import statistics

import numpy as np
from scipy.stats import spearmanr

from submit_wl_random_control_coordinates import validate_coordinate_plan


Z_EDGES = (0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.70)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def annulus_area_arcmin2(inner_radius_arcsec: float, outer_radius_arcsec: float) -> float:
    inner, outer = float(inner_radius_arcsec), float(outer_radius_arcsec)
    if not math.isfinite(inner) or not math.isfinite(outer) or inner < 0 or outer <= inner:
        raise ValueError("annulus radii must be finite with 0 <= inner < outer")
    return math.pi * (outer * outer - inner * inner) / (60.0 * 60.0)


def _quantiles(values: list[float]) -> dict[str, float]:
    array = np.asarray(values, dtype=float)
    return {key: float(np.quantile(array, q)) for key, q in (
        ("q10", 0.10), ("q25", 0.25), ("q50", 0.50), ("q75", 0.75), ("q90", 0.90),
    )}


def run(coordinate_plan_dir: Path, output_dir: Path) -> dict:
    coordinate_plan_dir, output_dir = Path(coordinate_plan_dir).resolve(), Path(output_dir).resolve()
    if output_dir.exists():
        raise FileExistsError(f"refusing to overwrite candidate-count audit: {output_dir}")
    manifest, batches, lenses, lens_ids, counts = validate_coordinate_plan(coordinate_plan_dir)
    if len(lenses) != 800 or len(lens_ids) != 800 or len(set(lens_ids)) != 800:
        raise ValueError("expected exactly 800 unique approved lens centers")
    if set(counts) != set(lens_ids) or any(int(value) < 0 for value in counts.values()):
        raise ValueError("count-only results do not match all approved lenses")

    count_plan_dir = Path(manifest["count_plan"])
    count_manifest_path = count_plan_dir / "random_control_count_plan.json"
    count_status_path = count_plan_dir / "cas_counts/count_status.json"
    count_summary_path = count_plan_dir / "cas_counts/count_summary.json"
    count_manifest = json.loads(count_manifest_path.read_text(encoding="utf-8"))
    count_summary = json.loads(count_summary_path.read_text(encoding="utf-8"))
    count_query = count_manifest.get("candidate_query", {})
    target_candidates = int(count_query.get("candidate_target_per_lens", 200))
    outer_radius = float(count_query.get("outer_radius_arcsec", 1800.0))
    if count_summary.get("status") != "complete" or count_summary.get("n_lenses") != 800:
        raise ValueError("count-only CAS summary is incomplete")
    if count_summary.get("sum_candidate_counts") != sum(counts.values()):
        raise ValueError("per-lens count results do not match aggregate count summary")

    rows = []
    for lens_id, lens in zip(lens_ids, lenses):
        if str(lens["lens_id"]) != str(lens_id):
            raise ValueError("lens order differs between coordinate plan and count results")
        z_lens = float(lens["z_lens"])
        aperture = float(lens["radius_arcsec"])
        inner_radius = 2.0 * aperture
        area = annulus_area_arcmin2(inner_radius, outer_radius)
        n_candidates = int(counts[lens_id])
        rows.append({
            "lens_id": lens_id, "ra_deg": float(lens["ra_deg"]), "dec_deg": float(lens["dec_deg"]),
            "z_lens": z_lens, "target_aperture_radius_arcsec": aperture,
            "random_center_exclusion_radius_arcsec": inner_radius,
            "query_outer_radius_arcsec": outer_radius, "available_annulus_area_arcmin2": area,
            "n_random_candidates": n_candidates,
            "candidate_density_per_arcmin2": n_candidates / area,
            "planned_target_candidates": target_candidates,
        })
    values = [row["n_random_candidates"] for row in rows]
    expected_raw_rows = int(manifest["candidate_pair_rows_expected"])
    if sum(values) != expected_raw_rows:
        raise ValueError("per-lens counts do not match coordinate-plan expected rows")
    areas = [row["available_annulus_area_arcmin2"] for row in rows]
    density = [row["candidate_density_per_arcmin2"] for row in rows]
    redshift = [row["z_lens"] for row in rows]
    z_bins = []
    for low, high in zip(Z_EDGES[:-1], Z_EDGES[1:]):
        subset = [row for row in rows if low <= row["z_lens"] < high]
        if subset:
            z_bins.append({
                "z_min": low, "z_max_exclusive": high, "n_lenses": len(subset),
                "candidate_count_median": float(statistics.median(row["n_random_candidates"] for row in subset)),
                "candidate_density_median_per_arcmin2": float(statistics.median(
                    row["candidate_density_per_arcmin2"] for row in subset)),
            })

    output_dir.mkdir(parents=True, exist_ok=False)
    table_path = output_dir / "random_candidate_count_by_lens.tsv"
    with table_path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    rho_count = float(spearmanr(redshift, values).statistic)
    rho_density = float(spearmanr(redshift, density).statistic)
    input_hashes = {
        str((coordinate_plan_dir / "random_control_coordinate_plan.json")): sha256(
            coordinate_plan_dir / "random_control_coordinate_plan.json"),
        str(count_manifest_path): sha256(count_manifest_path),
        str(count_status_path): sha256(count_status_path),
        str(count_summary_path): sha256(count_summary_path),
    }
    for batch in count_manifest["batches"]:
        path = Path(batch["path"])
        input_hashes[str(path)] = sha256(path)
    for batch in batches:
        path = Path(batch["path"])
        input_hashes[str(path)] = sha256(path)

    summary = {
        "schema": "astrafek-wl-random-coordinate-candidate-count-audit-v1",
        "status": "complete_aggregate_count_audit",
        "n_lenses": len(rows), "n_count_rows": len(rows),
        "sum_candidate_counts": sum(values), "planned_target_candidates_per_lens": target_candidates,
        "lenses_meeting_planned_candidate_count": sum(value >= target_candidates for value in values),
        "lenses_below_planned_candidate_count": sum(value < target_candidates for value in values),
        "zero_candidate_lenses": sum(value == 0 for value in values),
        "candidate_count_min": min(values), "candidate_count_max": max(values),
        "candidate_count_mean": float(statistics.mean(values)),
        "candidate_count_median": float(statistics.median(values)),
        "candidate_count_quantiles": _quantiles([float(value) for value in values]),
        "available_annulus_area_arcmin2_min": min(areas),
        "available_annulus_area_arcmin2_median": float(statistics.median(areas)),
        "available_annulus_area_arcmin2_max": max(areas),
        "candidate_density_per_arcmin2_min": min(density),
        "candidate_density_per_arcmin2_median": float(statistics.median(density)),
        "candidate_density_per_arcmin2_mean": float(statistics.mean(density)),
        "candidate_density_per_arcmin2_max": max(density),
        "spearman_candidate_count_vs_lens_redshift": rho_count,
        "spearman_area_normalized_candidate_density_vs_lens_redshift": rho_density,
        "lens_redshift_bins": z_bins,
        "candidate_query": count_query,
        "coordinate_plan_sha256": input_hashes[str(coordinate_plan_dir / "random_control_coordinate_plan.json")],
        "count_plan_sha256": manifest["count_plan_sha256"],
        "input_sha256": input_hashes,
        "table_sha256": sha256(table_path),
        "audit_script_sha256": sha256(Path(__file__)),
        "random_coordinates_read": False, "source_rows_read": False, "shape_rows_read": False,
        "photoz_rows_read": False, "images_read": False, "cas_requests_made": False,
        "science_measurement": False,
        "limitations": [
            "Counts are from the already approved count-only query; no random positions were read.",
            "Area normalization assumes a spherical annulus with the plan's 2x 3-Mpc exclusion and 30-arcmin outer radius.",
            "This is availability before the full all-cluster aperture-overlap, duplicate, image-support, and shape-quality screens.",
            "Spearman coefficients describe this fixed 800-lens sample and are not inferential significance tests.",
        ],
    }
    summary_path = output_dir / "summary.json"
    summary_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--coordinate-plan", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = run(args.coordinate_plan, args.output)
    print(json.dumps({key: result[key] for key in (
        "status", "n_lenses", "sum_candidate_counts", "candidate_count_median",
        "candidate_count_min", "candidate_count_max", "lenses_meeting_planned_candidate_count",
        "candidate_density_per_arcmin2_median", "spearman_candidate_count_vs_lens_redshift",
        "spearman_area_normalized_candidate_density_vs_lens_redshift", "science_measurement",
    )}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
