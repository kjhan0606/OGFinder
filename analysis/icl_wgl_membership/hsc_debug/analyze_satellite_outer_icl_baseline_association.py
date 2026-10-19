#!/usr/bin/env python3
"""Compare satellite outer-light association with raw and outer-subtracted ICL.

The ICL response is measured for the same cluster/control pairs in the 200--300
kpc annulus, first as measured and then after subtracting each pair's 400--500
kpc target-minus-control baseline. Satellite light is the equal-satellite
weighted 50--75 kpc profile. Filters are analyzed independently. Results remain
exploratory because the outer annulus is not certified blank sky and neither
measurement has complete PSF, sky, mask, or image-covariance calibration.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

from evaluate_prestack_pair_support import pair_block_labels


DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
FILTERS = ("i", "i2")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def partial_pearson(x: np.ndarray, y: np.ndarray,
                    z: np.ndarray, richness: np.ndarray) -> float:
    covariates = np.column_stack((z, richness))
    scales = np.std(covariates, axis=0)
    if (not np.isfinite(covariates).all() or np.any(scales <= 0)
            or np.std(x) <= 0 or np.std(y) <= 0):
        return math.nan
    design = np.column_stack((np.ones(len(x)),
                              (covariates - np.mean(covariates, axis=0)) / scales))
    xr = x - design @ np.linalg.lstsq(design, x, rcond=None)[0]
    yr = y - design @ np.linalg.lstsq(design, y, rcond=None)[0]
    if np.std(xr) <= 0 or np.std(yr) <= 0:
        return math.nan
    return float(np.corrcoef(xr, yr)[0, 1])


def metrics(rows: list[dict]) -> dict:
    satellite = np.asarray([r["satellite"] for r in rows], dtype=float)
    redshift = np.asarray([r["redshift"] for r in rows], dtype=float)
    richness = np.asarray([r["richness"] for r in rows], dtype=float)
    out = {}
    for label in ("icl_raw", "icl_outer_subtracted"):
        icl = np.asarray([r[label] for r in rows], dtype=float)
        if np.std(satellite) == 0 or np.std(icl) == 0:
            out[label] = {"spearman_rho": None, "partial_pearson_r": None}
            continue
        rho = float(spearmanr(icl, satellite).statistic)
        partial = partial_pearson(icl, satellite, redshift, richness)
        out[label] = {
            "spearman_rho": rho,
            "partial_pearson_r_control_z_richness": partial if math.isfinite(partial) else None,
            "mean_icl_delta_njy_arcsec2": float(np.mean(icl)),
            "mean_satellite_delta_njy_arcsec2": float(np.mean(satellite)),
        }
    return out


def block_bootstrap(rows: list[dict], replicas: int, seed: int) -> dict:
    target = [{"ra_deg": r["ra_deg"], "dec_deg": r["dec_deg"]} for r in rows]
    control = [{"ra_deg": r["control_ra_deg"], "dec_deg": r["control_dec_deg"]} for r in rows]
    _, labels, n_target, n_control, n_blocks, largest = pair_block_labels(
        target, control, np.ones(len(rows), dtype=bool))
    if n_blocks < 3:
        raise ValueError(f"at least three connected target/control blocks required, found {n_blocks}")
    members = [np.flatnonzero(labels == label) for label in np.unique(labels)]
    rng = np.random.default_rng(seed)
    names = ("icl_raw", "icl_outer_subtracted")
    quantities = ("spearman_rho", "partial_pearson_r_control_z_richness")
    draws = {name: {metric: [] for metric in quantities} for name in names}
    draws["corrected_minus_raw"] = {metric: [] for metric in quantities}
    for _ in range(replicas):
        chosen = rng.integers(0, len(members), size=len(members))
        sample = [rows[int(index)] for block in chosen for index in members[int(block)]]
        result = metrics(sample)
        for name in names:
            for metric in draws[name]:
                value = result[name][metric]
                if value is not None and math.isfinite(value):
                    draws[name][metric].append(value)
        for metric in quantities:
            raw, corrected = result["icl_raw"][metric], result["icl_outer_subtracted"][metric]
            if raw is not None and corrected is not None and math.isfinite(raw) and math.isfinite(corrected):
                draws["corrected_minus_raw"][metric].append(corrected - raw)

    intervals = {}
    for name, quantities in draws.items():
        intervals[name] = {}
        for quantity, values in quantities.items():
            values = np.asarray(values, dtype=float)
            if len(values) < 0.8 * replicas:
                raise ValueError(f"too few valid bootstrap draws for {name}/{quantity}: {len(values)}")
            intervals[name][quantity] = {
                "q025": float(np.quantile(values, 0.025)),
                "q975": float(np.quantile(values, 0.975)),
                "valid_replicates": int(len(values)),
            }
    return {
        "n_target_cells": int(n_target), "n_control_cells": int(n_control),
        "n_connected_blocks": int(n_blocks), "largest_block_pairs": int(largest),
        "intervals_95": intervals,
    }


def join_rows(icl_rows: list[dict[str, str]], outer_rows: list[dict[str, str]],
              satellite_rows: list[dict[str, str]], clusters: dict[str, dict],
              controls: dict[str, dict]) -> dict[str, list[dict]]:
    icl_by_key = {}
    for row in icl_rows:
        key = (str(row["target_id"]), row["filter"])
        if key in icl_by_key:
            raise ValueError(f"duplicate ICL pair row for {key}")
        icl_by_key[key] = row
    sat_by_key = {}
    for row in satellite_rows:
        if (row["cohort"] != "full_sample" or row["radial_bin_kpc"] != "50-75"):
            continue
        key = (str(row["cluster_id"]), row["filter"])
        if key in sat_by_key:
            raise ValueError(f"duplicate satellite profile for {key}")
        sat_by_key[key] = row

    joined = {band: [] for band in FILTERS}
    seen = set()
    for row in outer_rows:
        if (float(row["r_min_kpc"]) != 200.0 or float(row["r_max_kpc"]) != 300.0):
            continue
        band, target_id, control_id = row["filter"], str(row["target_id"]), str(row["control_id"])
        key = (target_id, band)
        if band not in joined:
            raise ValueError(f"unexpected filter {band}")
        if key in seen:
            raise ValueError(f"duplicate outer-subtracted ICL row for {key}")
        seen.add(key)
        icl = icl_by_key.get(key)
        satellite = sat_by_key.get(key)
        cluster = clusters.get(target_id)
        control_meta = controls.get(target_id)
        if icl is None or satellite is None or cluster is None or control_meta is None:
            continue
        if str(icl["control_id"]) != control_id or str(control_meta["control_id"]) != control_id:
            raise ValueError(f"target/control mapping mismatch for {key}")
        if (abs(float(icl["z"]) - float(cluster["z"])) > 1e-7
                or abs(float(icl["ra_deg"]) - float(cluster["ra_deg"])) > 1e-7
                or abs(float(icl["dec_deg"]) - float(cluster["dec_deg"])) > 1e-7):
            raise ValueError(f"target metadata mismatch for {key}")
        inner = float(row["inner_delta_njy_arcsec2"])
        outer = float(row["outer_400_500_delta_njy_arcsec2"])
        corrected = float(row["outer_subtracted_delta_njy_arcsec2"])
        inner_coverage = float(row["inner_coverage"])
        outer_coverage = float(row["outer_coverage"])
        if inner_coverage < 0.05 or outer_coverage < 0.05:
            raise ValueError(f"outer-baseline row is below its 5% support cut for {key}")
        if not math.isclose(corrected, inner - outer, rel_tol=1e-9, abs_tol=1e-9):
            raise ValueError(f"outer-baseline subtraction mismatch for {key}")
        joined[band].append({
            "cluster_id": target_id, "filter": band, "control_id": control_id,
            "icl_raw": inner, "outer_baseline": outer,
            "icl_outer_subtracted": corrected,
            "satellite": float(satellite["mean_delta_sb_njy_arcsec2"]),
            "n_satellites": int(satellite["n_satellites"]),
            "redshift": float(cluster["z"]), "richness": float(cluster["richness"]),
            "ra_deg": float(cluster["ra_deg"]), "dec_deg": float(cluster["dec_deg"]),
            "control_ra_deg": float(control_meta["ra_deg"]),
            "control_dec_deg": float(control_meta["dec_deg"]),
            "inner_coverage": inner_coverage, "outer_coverage": outer_coverage,
        })
    return joined


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=DATA)
    parser.add_argument("--stack-dir", type=Path, default=DATA / "icl_2d_matched_stack_2026-10-19_02")
    parser.add_argument("--outer-dir", type=Path, default=DATA / "icl_outer_annulus_baseline_2026-10-20_02")
    parser.add_argument("--satellite-dir", type=Path, default=DATA / "satellite_outer_profile_full_2026-10-18_02")
    parser.add_argument("--manifest", type=Path, default=DATA / "sample_manifest.json")
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--replicas", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=20261020)
    args = parser.parse_args()
    if args.replicas < 100:
        parser.error("replicas must be at least 100")
    if args.work.exists():
        parser.error(f"refusing to overwrite existing output: {args.work}")

    stack_summary_path = args.stack_dir / "summary.json"
    stack_path = args.stack_dir / "paired_radial_profiles.tsv"
    outer_summary_path = args.outer_dir / "summary.json"
    outer_path = args.outer_dir / "paired_outer_subtracted_profiles.tsv"
    sat_summary_path = args.satellite_dir / "summary.json"
    sat_path = args.satellite_dir / "cluster_profiles.tsv"
    extraction_path = args.data / "prestack_extraction_2026-10-15_01/extraction_summary.json"
    required = (stack_summary_path, stack_path, outer_summary_path, outer_path,
                sat_summary_path, sat_path, args.manifest, extraction_path)
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        parser.error("missing input(s): " + ", ".join(missing))

    stack_summary, outer_summary = read_json(stack_summary_path), read_json(outer_summary_path)
    satellite_summary = read_json(sat_summary_path)
    manifest, extraction = read_json(args.manifest), read_json(extraction_path)
    if stack_summary.get("status") != "complete_exploratory" or stack_summary.get("science_measurement") is not False:
        parser.error("ICL stack is not complete exploratory output")
    if outer_summary.get("status") != "complete_exploratory_sensitivity" or outer_summary.get("science_measurement") is not False:
        parser.error("outer-baseline audit is not complete exploratory output")
    if satellite_summary.get("status") != "complete" or satellite_summary.get("science_measurement") is not False:
        parser.error("satellite profile is incomplete or science-promoted")
    if manifest.get("schema") != "hsc-bcg-diffuse-sample-v1" or manifest.get("science_measurement") is not False:
        parser.error("unexpected sample manifest or science promotion")
    if extraction.get("status") != "complete" or extraction.get("science_measurement") is not False:
        parser.error("matched-control extraction is incomplete or science-promoted")
    expected_stack_inputs = stack_summary.get("input_sha256", {})
    if expected_stack_inputs.get(str(args.manifest.resolve())) != sha256(args.manifest):
        parser.error("sample manifest checksum differs from the ICL stack record")
    if outer_summary.get("inputs", {}).get("matched_control_extraction_sha256") != sha256(extraction_path):
        parser.error("matched-control extraction checksum differs from the ICL baseline record")
    if stack_summary.get("outputs", {}).get("paired_profiles_sha256") != sha256(stack_path):
        parser.error("ICL paired-profile checksum mismatch")
    if outer_summary.get("inputs", {}).get("paired_profiles_sha256") != sha256(stack_path):
        parser.error("outer-baseline input ICL profile checksum mismatch")
    if outer_summary.get("output_sha256") != sha256(outer_path):
        parser.error("outer-baseline output checksum mismatch")
    if satellite_summary.get("output_sha256", {}).get(str(sat_path)) != sha256(sat_path):
        parser.error("satellite cluster-profile checksum mismatch")

    clusters = {str(row["id"]): row for row in manifest["clusters"]}
    controls = {str(row["target_cluster_id"]): row for row in extraction["records"]}
    if len(clusters) != len(manifest["clusters"]) or len(controls) != len(extraction["records"]):
        parser.error("duplicate cluster or control metadata rows")
    joined = join_rows(read_tsv(stack_path), read_tsv(outer_path),
                       read_tsv(sat_path), clusters, controls)

    results = []
    detail_rows = []
    for band_index, band in enumerate(FILTERS):
        rows = joined[band]
        if len(rows) < 12:
            parser.error(f"insufficient matched clusters for {band}: {len(rows)}")
        point = metrics(rows)
        boot = block_bootstrap(rows, args.replicas, args.seed + band_index * 10000)
        results.append({"filter": band, "n_clusters": len(rows), "point_estimates": point,
                        "corrected_minus_raw": {
                            metric: (point["icl_outer_subtracted"][metric] - point["icl_raw"][metric]
                                     if point["icl_outer_subtracted"][metric] is not None
                                     and point["icl_raw"][metric] is not None else None)
                            for metric in ("spearman_rho", "partial_pearson_r_control_z_richness")},
                        "block_bootstrap": boot})
        detail_rows.extend(rows)

    args.work.mkdir(parents=True, exist_ok=False)
    fields = ["cluster_id", "filter", "control_id", "icl_raw", "outer_baseline",
              "icl_outer_subtracted", "satellite", "n_satellites", "redshift", "richness",
              "ra_deg", "dec_deg", "control_ra_deg", "control_dec_deg",
              "inner_coverage", "outer_coverage"]
    with (args.work / "matched_cluster_values.tsv").open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", extrasaction="raise")
        writer.writeheader()
        writer.writerows(detail_rows)
    summary = {
        "schema": "astrafek-satellite-outer-icl-baseline-association-v1",
        "status": "complete_exploratory", "science_measurement": False,
        "method": {
            "icl_raw": "same matched target/control pair, 200-300 kpc delta from the ICL profile",
            "icl_outer_subtracted": "raw 200-300 kpc pair delta minus that pair's 400-500 kpc target-control delta",
            "satellite": "full-sample equal-satellite-weighted cluster mean at 50-75 kpc",
            "support": "inherited 5% ICL inner/outer common-area cut; satellite outer-bin support is independent",
            "join": "exact cluster ID, exact filter, and exact matched-control ID",
            "filters": "i and i2 analyzed independently",
            "covariates": "redshift and richness linear residuals for partial Pearson diagnostic",
            "uncertainty": "95% percentile intervals from connected 5-degree target/control block bootstrap",
            "multiplicity": "descriptive comparison; no confirmatory p-values or multiplicity claims",
        },
        "limitations": [
            "the 400-500 kpc annulus is not certified blank sky; subtracting it may remove real extended cluster light",
            "satellite profiles use independent radial-bin support and no member mass or CAMIRA weight",
            "automatic centers and masks, PSF wings, upstream sky transfer, extinction, passband response, and full coadd pixel covariance remain unvalidated",
            "block-bootstrap intervals represent only the chosen connected target/control sky blocks",
            "these selected-sample associations do not establish a physical ICL-satellite relation",
        ],
        "bootstrap": {"replicas": args.replicas, "seed": args.seed},
        "inputs": {
            "icl_stack_summary": {"path": str(stack_summary_path.resolve()), "sha256": sha256(stack_summary_path)},
            "icl_paired_profiles": {"path": str(stack_path.resolve()), "sha256": sha256(stack_path)},
            "outer_baseline_summary": {"path": str(outer_summary_path.resolve()), "sha256": sha256(outer_summary_path)},
            "outer_baseline_profiles": {"path": str(outer_path.resolve()), "sha256": sha256(outer_path)},
            "satellite_summary": {"path": str(sat_summary_path.resolve()), "sha256": sha256(sat_summary_path)},
            "satellite_cluster_profiles": {"path": str(sat_path.resolve()), "sha256": sha256(sat_path)},
            "sample_manifest": {"path": str(args.manifest.resolve()), "sha256": sha256(args.manifest)},
            "matched_control_extraction": {"path": str(extraction_path.resolve()), "sha256": sha256(extraction_path)},
        },
        "results": results,
        "outputs_sha256": {"matched_cluster_values.tsv": sha256(args.work / "matched_cluster_values.tsv")},
        "script_sha256": sha256(Path(__file__).resolve()),
    }
    with (args.work / "summary.json").open("x", encoding="utf-8") as stream:
        json.dump(summary, stream, ensure_ascii=False, indent=2, allow_nan=False)
        stream.write("\n")
    print(json.dumps({"status": summary["status"], "science_measurement": False,
                      "n_clusters": {r["filter"]: r["n_clusters"] for r in results},
                      "results": results, "output": str(args.work)},
                     ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
