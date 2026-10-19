#!/usr/bin/env python3
"""Test whether the ICL--satellite relation survives HSC region fixed effects."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from audit_satellite_icl_spatial_blocks import load_verified_rows
from run_wl_hsc_tract_region_influence import (
    nearest_region_labels,
    tract_centers_from_patch_qa,
    tract_region_labels,
)

ROOT = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
DEFAULT_ASSOC = ROOT / "satellite_outer_icl_association_perbin_2026-10-18_03"
DEFAULT_ICL = ROOT / "icl_annulus_support_2026-10-18_01"
DEFAULT_EXTRACTION = ROOT / "prestack_extraction_2026-10-15_01/extraction_summary.json"
DEFAULT_BASELINE = ROOT / "satellite_icl_spatial_blocks_2026-10-20_02/summary.json"
DEFAULT_PATCH_QA = ROOT / "catalogues_2026-10-13_01/patch_qa.csv"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def partial_correlation(rows, *, region_effects=()):
    """Partial correlation after nuisance regression and optional field dummies."""
    rows = list(rows)
    if len(rows) < 5:
        raise ValueError("partial correlation requires at least five rows")
    x = np.asarray([float(row["icl_delta_200_300"]) for row in rows])
    y = np.asarray([float(row["satellite_delta_50_75"]) for row in rows])
    z = np.asarray([float(row["redshift"]) for row in rows])
    richness = np.asarray([float(row["richness"]) for row in rows])
    if not all(np.isfinite(a).all() for a in (x, y, z, richness)):
        raise ValueError("partial-correlation inputs must be finite")
    columns = [np.ones(len(rows)), z, richness]
    terms = ["intercept", "redshift", "richness"]
    redundant_region_terms = []
    for effect in region_effects:
        if effect not in ("target_region", "control_region"):
            raise ValueError(f"unsupported region effect: {effect}")
        categories = sorted({int(row[effect]) for row in rows})
        for category in categories[1:]:
            column = np.asarray([int(row[effect]) == category for row in rows], dtype=float)
            name = f"{effect}[{category}]"
            if any(np.array_equal(column, existing) for existing in columns):
                redundant_region_terms.append(name)
                continue
            columns.append(column)
            terms.append(name)
    design = np.column_stack(columns)
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
        raise ValueError("partial correlation is undefined for constant residuals")
    return {
        "partial_r": float(rx @ ry / denominator),
        "n_rows": len(rows),
        "design_rank": rank,
        "design_columns": terms,
        "redundant_region_terms": redundant_region_terms,
        "residual_degrees_of_freedom": residual_dof,
    }


def assign_pair_regions(rows, patch_qa: Path, adjacency_deg=2.0):
    tract_keys, centers = tract_centers_from_patch_qa(patch_qa)
    tract_regions = tract_region_labels(centers, adjacency_deg)
    target = nearest_region_labels(
        [float(row["target_ra_deg"]) for row in rows],
        [float(row["target_dec_deg"]) for row in rows], centers, tract_regions,
    )
    control = nearest_region_labels(
        [float(row["control_ra_deg"]) for row in rows],
        [float(row["control_dec_deg"]) for row in rows], centers, tract_regions,
    )
    prepared = [{**row, "target_region": int(t), "control_region": int(c)}
                for row, t, c in zip(rows, target, control)]
    return tract_keys, centers, tract_regions, prepared


def analyze(association, icl, extraction, baseline_path, patch_qa,
            output_path, adjacency_deg=2.0):
    output_path = Path(output_path).resolve()
    if output_path.exists():
        raise FileExistsError(f"refusing to overwrite existing output: {output_path}")
    baseline_path = Path(baseline_path).resolve()
    baseline = json.loads(baseline_path.read_text(encoding="utf-8"))
    if baseline.get("schema") != "hsc-satellite-icl-spatial-block-bootstrap-v1" or baseline.get("science_measurement") is not False:
        raise ValueError("unexpected 94-row exploratory baseline")
    rows, provenance = load_verified_rows(association, icl, extraction)
    if len(rows) != 94:
        raise ValueError(f"expected the fixed 94-row sample; found {len(rows)}")
    for key, value in provenance.items():
        if baseline.get("provenance", {}).get(key) != value:
            raise ValueError(f"baseline provenance mismatch for {key}")
    tract_keys, centers, tract_regions, prepared = assign_pair_regions(
        rows, Path(patch_qa), adjacency_deg,
    )
    baseline_r = {
        row["filter"]: float(row["partial_spatial_block_bootstrap"]["partial_r_point"])
        for row in baseline.get("results", []) if row.get("selection") == "all_supported_pairs"
    }
    if set(baseline_r) != {"i", "i2"}:
        raise ValueError("baseline lacks both all-supported filter values")

    effects = {
        "redshift_richness_only": (),
        "target_region_only": ("target_region",),
        "control_region_only": ("control_region",),
        "target_and_control_region": ("target_region", "control_region"),
    }
    filter_results = []
    for band in ("i2", "i"):
        selected = [row for row in prepared if row["filter"] == band]
        if len(selected) != (68 if band == "i2" else 26):
            raise ValueError(f"fixed {band} sample count changed")
        results = []
        for name, region_effects in effects.items():
            result = partial_correlation(selected, region_effects=region_effects)
            if name == "redshift_richness_only":
                if not math.isclose(result["partial_r"], baseline_r[band], rel_tol=0, abs_tol=1e-12):
                    raise ValueError(f"baseline partial correlation mismatch for {band}")
                result["baseline_reproduction"] = "exact_match_within_1e-12"
            else:
                result["baseline_reproduction"] = None
            result["model"] = name
            result["delta_r_from_base"] = (result["partial_r"] - baseline_r[band]
                                           if name != "redshift_richness_only" else 0.0)
            results.append(result)
        filter_results.append({"filter": band, "n_pairs": len(selected), "models": results})

    threshold_sensitivity = []
    for threshold in (1.6, 1.8, 2.0, 2.2, 2.5, 3.0):
        regions = tract_region_labels(centers, threshold)
        target = nearest_region_labels(
            [float(row["target_ra_deg"]) for row in rows],
            [float(row["target_dec_deg"]) for row in rows], centers, regions,
        )
        control = nearest_region_labels(
            [float(row["control_ra_deg"]) for row in rows],
            [float(row["control_dec_deg"]) for row in rows], centers, regions,
        )
        counts = {}
        for row, t, c in zip(rows, target, control):
            counts.setdefault(row["filter"], {})
            for region in set((int(t), int(c))):
                key = str(region)
                counts[row["filter"]][key] = counts[row["filter"]].get(key, 0) + 1
        threshold_sensitivity.append({
            "adjacency_deg": threshold,
            "n_tract_components": int(len(np.unique(regions))),
            "pairs_touching_by_filter_and_region": counts,
        })
    expected = threshold_sensitivity[0]
    if any(row["n_tract_components"] != expected["n_tract_components"]
           or row["pairs_touching_by_filter_and_region"] != expected["pairs_touching_by_filter_and_region"]
           for row in threshold_sensitivity):
        raise ValueError("target/control endpoint partition changes over the threshold grid")

    result = {
        "schema": "astrafek-satellite-icl-region-fixed-effects-v1",
        "status": "complete_exploratory_sensitivity",
        "science_measurement": False,
        "covariance_estimated": False,
        "method": {
            "sample": "fixed 94 validated cluster/filter rows, no reselection",
            "statistic": "partial correlation of 200-300 kpc ICL and 50-75 kpc satellite light",
            "base_nuisance_terms": ["redshift", "richness"],
            "additional_models": ["target-region fixed effects", "control-region fixed effects", "both endpoints' region fixed effects"],
            "collinearity": "identical target/control region indicator columns are included once and recorded as redundant",
            "grouping": "nearest HSC patch-QA tract center; connected components at the configured adjacency threshold",
            "interpretation": "confounding sensitivity, not a confirmatory model or detection test",
        },
        "adjacency_deg": float(adjacency_deg),
        "n_tracts": len(tract_keys),
        "n_tract_components": int(len(np.unique(tract_regions))),
        "filter_results": filter_results,
        "partition_threshold_sensitivity": threshold_sensitivity,
        "provenance": provenance,
        "limitations": [
            "the fixed 94-row sample is selected by prior support/extraction gates",
            "approximate tract-center regions may not match exact survey polygons",
            "fixed effects do not model all sky, PSF, mask, covariance, or source-selection systematics",
            "the small sample, especially i, limits within-region interpretation",
            "no random-lens subtraction or new HSC request was made",
        ],
        "inputs": {
            "baseline": str(baseline_path),
            "baseline_sha256": sha256(baseline_path),
            "patch_qa": str(Path(patch_qa).resolve()),
            "patch_qa_sha256": sha256(Path(patch_qa)),
        },
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("x", encoding="utf-8") as stream:
        json.dump(result, stream, indent=2, ensure_ascii=False, allow_nan=False)
        stream.write("\n")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--association", type=Path, default=DEFAULT_ASSOC)
    parser.add_argument("--icl", type=Path, default=DEFAULT_ICL)
    parser.add_argument("--extraction-summary", type=Path, default=DEFAULT_EXTRACTION)
    parser.add_argument("--baseline", type=Path, default=DEFAULT_BASELINE)
    parser.add_argument("--patch-qa", type=Path, default=DEFAULT_PATCH_QA)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--adjacency-deg", type=float, default=2.0)
    args = parser.parse_args()
    result = analyze(args.association, args.icl, args.extraction_summary,
                     args.baseline, args.patch_qa, args.output,
                     args.adjacency_deg)
    print(json.dumps({
        "status": result["status"], "science_measurement": False,
        "filter_results": result["filter_results"],
    }, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
