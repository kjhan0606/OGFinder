#!/usr/bin/env python3
"""Check ICL--PSF core dependence after approximate HSC field-region effects."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from run_wl_hsc_tract_region_influence import (
    nearest_region_labels,
    tract_centers_from_patch_qa,
    tract_region_labels,
)

ROOT = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
DEFAULT_ASSOCIATION = ROOT / "icl_psf_core_dependence_2026-10-20_02/summary.json"
DEFAULT_PATCH_QA = ROOT / "catalogues_2026-10-13_01/patch_qa.csv"
EXPECTED_ASSOCIATION_SHA256 = "2fbbf6f89c16eea38153f19002b6e14b8b63d28b72d73c6fabd7ab000af5a805"
BASE_MODELS = {
    "unadjusted": (),
    "redshift_richness_common_coverage": ("redshift", "richness", "common_coverage"),
    "redshift_richness_separate_coverage": (
        "redshift", "richness", "target_coverage", "control_coverage"),
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def partial_correlation(rows: list[dict], covariates: tuple[str, ...],
                        region_effect: bool = False) -> dict:
    """Residualize ICL and PSF size against numeric nuisance terms and regions."""
    if len(rows) < 5:
        raise ValueError("partial correlation requires at least five rows")
    x = np.asarray([float(row["icl_delta"]) for row in rows])
    y = np.asarray([float(row["psf_fwhm_major_arcsec"]) for row in rows])
    columns = [np.ones(len(rows))]
    names = ["intercept"]
    for covariate in covariates:
        values = np.asarray([float(row[covariate]) for row in rows])
        if not np.isfinite(values).all():
            raise ValueError(f"non-finite covariate: {covariate}")
        columns.append(values)
        names.append(covariate)
    if region_effect:
        categories = sorted({int(row["target_region"]) for row in rows})
        for category in categories[1:]:
            columns.append(np.asarray(
                [int(row["target_region"]) == category for row in rows], dtype=float))
            names.append(f"tract_region[{category}]")
    if not np.isfinite(x).all() or not np.isfinite(y).all():
        raise ValueError("correlation inputs must be finite")
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
        "n_frames": len(rows),
        "design_rank": rank,
        "design_columns": names,
        "residual_degrees_of_freedom": residual_dof,
    }


def assign_regions(rows: list[dict], patch_qa: Path, adjacency_deg: float) -> tuple[list[dict], dict]:
    tract_keys, centers = tract_centers_from_patch_qa(patch_qa)
    tract_regions = tract_region_labels(centers, adjacency_deg)
    target_regions = nearest_region_labels(
        [float(row["target_ra_deg"]) for row in rows],
        [float(row["target_dec_deg"]) for row in rows], centers, tract_regions,
    )
    control_regions = nearest_region_labels(
        [float(row["control_ra_deg"]) for row in rows],
        [float(row["control_dec_deg"]) for row in rows], centers, tract_regions,
    )
    prepared = [{**row, "target_region": int(target), "control_region": int(control)}
                for row, target, control in zip(rows, target_regions, control_regions)]
    return prepared, {
        "n_tracts": len(tract_keys),
        "n_connected_regions": int(len(np.unique(tract_regions))),
        "target_control_same_region_rows": int(np.count_nonzero(target_regions == control_regions)),
        "target_control_different_region_rows": int(np.count_nonzero(target_regions != control_regions)),
    }


def analyze(association_path: Path, patch_qa: Path, output_path: Path,
            adjacency_deg: float = 2.0) -> dict:
    association_path, patch_qa = Path(association_path).resolve(), Path(patch_qa).resolve()
    output_path = Path(output_path).resolve()
    if output_path.exists():
        raise FileExistsError(f"refusing to overwrite existing output: {output_path}")
    association_digest = sha256(association_path)
    if association_digest != EXPECTED_ASSOCIATION_SHA256:
        raise ValueError("PSF--ICL association summary differs from its audited checksum")
    association = json.loads(association_path.read_text(encoding="utf-8"))
    if (association.get("schema") != "astrafek-icl-psf-core-dependence-v1"
            or association.get("status") != "complete_exploratory_sensitivity"
            or association.get("science_measurement") is not False
            or association.get("psf_wings_validated") is not False):
        raise ValueError("input is incomplete or marked as a science measurement")
    rows = association.get("records", [])
    if len(rows) != 57 or {row.get("filter") for row in rows} != {"i", "i2"}:
        raise ValueError("expected the fixed 57-frame i/i2 support sample")
    prepared, region_summary = assign_regions(rows, patch_qa, adjacency_deg)
    if region_summary["target_control_different_region_rows"] != 0:
        raise ValueError("paired endpoints fall in different approximate tract regions")

    by_filter = {}
    for band, expected_n in (("i", 13), ("i2", 44)):
        selected = [row for row in prepared if row["filter"] == band]
        if len(selected) != expected_n:
            raise ValueError(f"fixed {band} sample count changed")
        n_regions = len({row["target_region"] for row in selected})
        models = {}
        for model_name, covariates in BASE_MODELS.items():
            base = partial_correlation(selected, covariates)
            adjusted = partial_correlation(selected, covariates, region_effect=True)
            models[model_name] = {
                "base": base,
                "plus_tract_region_fixed_effects": adjusted,
                "delta_r_after_region_effects": adjusted["partial_r"] - base["partial_r"],
            }
        by_filter[band] = {"n_frames": len(selected), "n_regions_touched": n_regions,
                           "region_counts": {str(region): sum(
                               row["target_region"] == region for row in selected)
                               for region in sorted({row["target_region"] for row in selected})},
                           "models": models}

    threshold_sensitivity = []
    for threshold in (1.6, 1.8, 2.0, 2.2, 2.5, 3.0):
        labeled, summary = assign_regions(rows, patch_qa, threshold)
        if summary["target_control_different_region_rows"]:
            raise ValueError("target/control region pairing changes over threshold grid")
        estimate = {}
        for band in ("i", "i2"):
            subset = [row for row in labeled if row["filter"] == band]
            estimate[band] = partial_correlation(
                subset, BASE_MODELS["redshift_richness_common_coverage"],
                region_effect=True,
            )
        threshold_sensitivity.append({"adjacency_deg": threshold,
                                      "n_connected_regions": summary["n_connected_regions"],
                                      "by_filter": estimate})

    result = {
        "schema": "astrafek-icl-psf-region-fixed-effects-v1",
        "status": "complete_exploratory_sensitivity",
        "science_measurement": False,
        "psf_wings_validated": False,
        "method": {
            "response_pair": "ICL 200-300 kpc target-minus-control surface brightness and local center-PSF moment FWHM",
            "base_adjustments": "same numeric covariates as the existing PSF-core audit",
            "region_adjustment": "nearest patch-QA tract center; connected tract components at 2 degree adjacency",
            "paired_region_handling": "target and matched control map to the same approximate component in all 57 rows; one component indicator set is used",
            "uncertainty": "no p-values or calibrated confidence intervals; fixed-effect point estimates only",
            "interpretation": "sensitivity to broad field-region mean differences, not full survey-region or causal adjustment",
        },
        "adjacency_deg": adjacency_deg,
        "region_summary": region_summary,
        "by_filter": by_filter,
        "threshold_sensitivity_common_coverage": threshold_sensitivity,
        "provenance": {"association_summary": str(association_path),
                       "association_summary_sha256": association_digest,
                       "patch_qa": str(patch_qa), "patch_qa_sha256": sha256(patch_qa)},
        "limitations": [
            "57 rows are a support-selected subset (13 i and 44 i2), not the full 800-cluster sample",
            "nearest tract-center connected components approximate field regions and are not exact survey footprints",
            "i has only 13 rows, so region-adjusted models retain few residual degrees of freedom",
            "region fixed effects do not model PSF wings, sky, masks, noise, or source-selection covariance",
            "no random-lens subtraction, new HSC request, or science measurement was made",
        ],
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("x", encoding="utf-8") as stream:
        json.dump(result, stream, indent=2, ensure_ascii=False, allow_nan=False)
        stream.write("\n")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--association", type=Path, default=DEFAULT_ASSOCIATION)
    parser.add_argument("--patch-qa", type=Path, default=DEFAULT_PATCH_QA)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--adjacency-deg", type=float, default=2.0)
    args = parser.parse_args()
    result = analyze(args.association, args.patch_qa, args.output, args.adjacency_deg)
    print(json.dumps({"status": result["status"], "science_measurement": False,
                      "region_summary": result["region_summary"],
                      "by_filter": result["by_filter"],
                      "threshold_sensitivity_common_coverage": result["threshold_sensitivity_common_coverage"]},
                     ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
