#!/usr/bin/env python3
"""Test the fixed 94-row ICL--satellite correlation's area-support sensitivity.

The association sample is not reselected. This diagnostic compares the existing
redshift/richness-adjusted partial Pearson correlation with versions that also
condition on annulus area-support measurements. Spatial-block intervals are
descriptive and do not correct the selected sample or establish causality.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from audit_satellite_icl_region_fixed_effects import (
    DEFAULT_ASSOC, DEFAULT_EXTRACTION, DEFAULT_ICL,
)
from audit_satellite_icl_spatial_blocks import coordinate_labels, load_verified_rows, read_tsv

ROOT = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
DEFAULT_BASELINE = ROOT / "satellite_icl_spatial_blocks_2026-10-20_02/summary.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def join_coverage(rows: list[dict], annulus_rows: list[dict]) -> list[dict]:
    coverage = {}
    for row in annulus_rows:
        if row.get("radial_bin_kpc") != "200-300":
            continue
        key = (str(row["id"]), str(row["filter"]))
        if key in coverage:
            raise ValueError(f"duplicate ICL annulus coverage row: {key}")
        target = float(row["target_coverage"])
        control = float(row["control_coverage"])
        common = float(row["common_coverage"])
        if (not all(math.isfinite(v) and 0 <= v <= 1 for v in (target, control, common))
                or not math.isclose(common, min(target, control), rel_tol=0, abs_tol=1e-12)):
            raise ValueError(f"invalid or inconsistent area support for {key}")
        coverage[key] = {"icl_target_coverage": target,
                         "icl_control_coverage": control,
                         "icl_common_coverage": common}

    output, seen = [], set()
    for row in rows:
        key = (str(row["cluster_id"]), str(row["filter"]))
        if key in seen:
            raise ValueError(f"duplicate association row: {key}")
        seen.add(key)
        if key not in coverage:
            raise ValueError(f"association row lacks ICL area support: {key}")
        output.append({**row, **coverage[key]})
    return output


def partial_correlation(rows: list[dict], covariates: tuple[str, ...]) -> dict:
    rows = list(rows)
    if len(rows) < 5:
        raise ValueError("partial correlation requires at least five rows")
    x = np.asarray([float(row["icl_delta_200_300"]) for row in rows])
    y = np.asarray([float(row["satellite_delta_50_75"]) for row in rows])
    vectors = [np.asarray([float(row[name]) for row in rows]) for name in covariates]
    if not all(np.isfinite(vector).all() for vector in (x, y, *vectors)):
        raise ValueError("partial-correlation inputs must be finite")
    design = np.column_stack((np.ones(len(rows)), *vectors))
    rank = int(np.linalg.matrix_rank(design))
    if rank != design.shape[1]:
        raise ValueError("nuisance design is rank deficient")
    residual_dof = len(rows) - rank
    if residual_dof < 3:
        raise ValueError("fewer than three residual degrees of freedom remain")
    rx = x - design @ np.linalg.lstsq(design, x, rcond=None)[0]
    ry = y - design @ np.linalg.lstsq(design, y, rcond=None)[0]
    denominator = float(np.linalg.norm(rx) * np.linalg.norm(ry))
    if not math.isfinite(denominator) or denominator <= 0:
        raise ValueError("partial correlation is undefined for a constant residual")
    return {"partial_r": float(rx @ ry / denominator),
            "n_rows": len(rows), "design_rank": rank,
            "residual_degrees_of_freedom": residual_dof}


MODELS = {
    "redshift_richness": ("redshift", "richness"),
    "redshift_richness_common_coverage": (
        "redshift", "richness", "icl_common_coverage"),
    "redshift_richness_target_control_coverage": (
        "redshift", "richness", "icl_target_coverage", "icl_control_coverage"),
}


def block_bootstrap(rows: list[dict], labels: np.ndarray, replicas: int,
                    seed: int) -> dict:
    if len(labels) != len(rows):
        raise ValueError("one spatial block label is required per association row")
    blocks = np.unique(labels)
    if len(blocks) < 6:
        raise ValueError("at least six connected spatial blocks are required")
    rng = np.random.default_rng(seed)
    indices_by_block = {block: np.flatnonzero(labels == block) for block in blocks}
    draws = {name: [] for name in MODELS}
    deltas = {name: [] for name in MODELS if name != "redshift_richness"}
    invalid = 0
    for _ in range(replicas):
        chosen = rng.choice(blocks, size=len(blocks), replace=True)
        indices = np.concatenate([indices_by_block[block] for block in chosen])
        sample = [rows[int(index)] for index in indices]
        try:
            results = {name: partial_correlation(sample, covariates)
                       for name, covariates in MODELS.items()}
        except ValueError:
            invalid += 1
            continue
        for name, result in results.items():
            draws[name].append(result["partial_r"])
        base = results["redshift_richness"]["partial_r"]
        for name in deltas:
            deltas[name].append(results[name]["partial_r"] - base)
    if not draws["redshift_richness"]:
        raise ValueError("no valid spatial-block bootstrap replicas")
    quantiles = {}
    for name, values in draws.items():
        quantiles[name] = np.percentile(values, [2.5, 50, 97.5]).tolist()
    delta_quantiles = {name: np.percentile(values, [2.5, 50, 97.5]).tolist()
                       for name, values in deltas.items() if values}
    return {"block_count": int(len(blocks)), "replicas_requested": replicas,
            "replicas_valid": len(draws["redshift_richness"]),
            "replicas_invalid": invalid, "seed": seed,
            "partial_r_percentiles_2_5_50_97_5": quantiles,
            "paired_delta_r_percentiles_2_5_50_97_5": delta_quantiles}


def analyze(association_dir: Path, icl_dir: Path, extraction_summary: Path,
            baseline_path: Path, output_path: Path, replicas: int = 10000,
            seed: int = 20261020) -> dict:
    if replicas < 1000:
        raise ValueError("at least 1000 spatial-block replicas are required")
    output_path = Path(output_path).resolve()
    if output_path.exists():
        raise FileExistsError(f"refusing to overwrite existing output: {output_path}")
    rows, provenance = load_verified_rows(association_dir, icl_dir, extraction_summary)
    icl_dir = Path(icl_dir).resolve()
    annulus_path = icl_dir / "annulus_pair_profiles.tsv"
    rows = join_coverage(rows, read_tsv(annulus_path))
    if len(rows) != 94:
        raise ValueError(f"expected the fixed 94-row association sample; found {len(rows)}")
    baseline_path = Path(baseline_path).resolve()
    baseline = json.loads(baseline_path.read_text(encoding="utf-8"))
    if (baseline.get("schema") != "hsc-satellite-icl-spatial-block-bootstrap-v1"
            or baseline.get("science_measurement") is not False):
        raise ValueError("unexpected or science-promoted archived baseline")
    for key, value in provenance.items():
        if baseline.get("provenance", {}).get(key) != value:
            raise ValueError(f"baseline input provenance mismatch: {key}")
    archived = {row["filter"]: float(row["partial_spatial_block_bootstrap"]["partial_r_point"])
                for row in baseline.get("results", [])
                if row.get("selection") == "all_supported_pairs"}
    if set(archived) != {"i", "i2"}:
        raise ValueError("archived baseline does not contain both filters")

    by_filter = {}
    for band in ("i2", "i"):
        selected = [row for row in rows if row["filter"] == band]
        expected_n = 68 if band == "i2" else 26
        if len(selected) != expected_n:
            raise ValueError(f"fixed {band} support changed: {len(selected)} != {expected_n}")
        point = {name: partial_correlation(selected, covariates)
                 for name, covariates in MODELS.items()}
        if not math.isclose(point["redshift_richness"]["partial_r"], archived[band],
                            rel_tol=0, abs_tol=1e-12):
            raise ValueError(f"baseline correlation does not reproduce for {band}")
        labels = coordinate_labels(selected)
        bootstrap = block_bootstrap(selected, labels, replicas, seed + (1 if band == "i" else 0))
        by_filter[band] = {
            "n_pairs": len(selected),
            "models": {name: {**point[name], "covariates": list(MODELS[name]),
                              "delta_r_from_redshift_richness": (
                                  point[name]["partial_r"] - point["redshift_richness"]["partial_r"])}
                        for name in MODELS},
            "spatial_block_bootstrap": bootstrap,
        }

    result = {
        "schema": "astrafek-satellite-icl-coverage-adjustment-v1",
        "status": "complete_exploratory_sensitivity",
        "science_measurement": False,
        "causal_effect_estimated": False,
        "sample": "fixed hash-verified 94 cluster/filter association rows; no reselection",
        "method": {
            "response_pair": "ICL 200-300 kpc target-minus-control versus satellite 50-75 kpc target-minus-control",
            "statistic": "Pearson correlation between residuals after linear nuisance regression",
            "coverage_models": "common=min(target,control), or separate target and control annulus coverage; no inverse-probability weights",
            "uncertainty": "paired 5-degree connected target/control spatial-block bootstrap; percentile spread is descriptive, not a calibrated confidence interval",
            "interpretation": "measurement-support sensitivity only; area support may be associated with selection and signal, so adjustment is not causal correction",
        },
        "by_filter": by_filter,
        "provenance": provenance,
        "inputs": {
            "icl_annulus_table": str(annulus_path),
            "icl_annulus_table_sha256": sha256(annulus_path),
            "baseline_summary": str(baseline_path),
            "baseline_summary_sha256": sha256(baseline_path),
        },
        "limitations": [
            "the 94 rows are a selected intersection of ICL and satellite support, not the full cluster population",
            "coverage adjustment cannot recover rows excluded by prior support cuts or establish population representativeness",
            "the i filter has only 26 rows and the connected-block intervals can be unstable",
            "the analysis omits unresolved PSF wings, sky transfer, complete image covariance, and random-lens subtraction",
            "no HSC request or download was made",
        ],
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("x", encoding="utf-8") as stream:
        json.dump(result, stream, indent=2, ensure_ascii=False, allow_nan=False)
        stream.write("\n")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--association", type=Path, default=DEFAULT_ASSOC)
    parser.add_argument("--icl", type=Path, default=DEFAULT_ICL)
    parser.add_argument("--extraction-summary", type=Path, default=DEFAULT_EXTRACTION)
    parser.add_argument("--baseline", type=Path, default=DEFAULT_BASELINE)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--replicas", type=int, default=10000)
    parser.add_argument("--seed", type=int, default=20261020)
    args = parser.parse_args()
    result = analyze(args.association, args.icl, args.extraction_summary,
                     args.baseline, args.output, args.replicas, args.seed)
    print(json.dumps({"status": result["status"], "science_measurement": result["science_measurement"],
                      "by_filter": result["by_filter"]}, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
