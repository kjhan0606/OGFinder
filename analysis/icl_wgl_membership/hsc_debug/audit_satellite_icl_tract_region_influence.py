#!/usr/bin/env python3
"""Test ICL--satellite correlation sensitivity to connected HSC tract regions."""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from audit_satellite_icl_perbin_block_scale import partial_pearson
from audit_satellite_icl_spatial_blocks import load_verified_rows, spearman_rho
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
DEFAULT_MANIFEST = ROOT / "sample_manifest.json"
DEFAULT_PATCH_QA = ROOT / "catalogues_2026-10-13_01/patch_qa.csv"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def prepare_rows(rows, centers_xyz, tract_regions):
    """Attach nearest-tract connected-region labels to both pair endpoints."""
    target_regions = nearest_region_labels(
        [float(row["target_ra_deg"]) for row in rows],
        [float(row["target_dec_deg"]) for row in rows], centers_xyz, tract_regions,
    )
    control_regions = nearest_region_labels(
        [float(row["control_ra_deg"]) for row in rows],
        [float(row["control_dec_deg"]) for row in rows], centers_xyz, tract_regions,
    )
    prepared = []
    for row, target_region, control_region in zip(rows, target_regions, control_regions):
        prepared.append({
            **row,
            "satellite_delta": row["satellite_delta_50_75"],
            "target_region": int(target_region),
            "control_region": int(control_region),
        })
    return prepared


def endpoint_touch_counts(rows):
    """Count per filter how many pairs touch each region at either endpoint."""
    result = {}
    for band in ("i", "i2"):
        selected = [row for row in rows if row["filter"] == band]
        region_ids = sorted({int(row[key]) for row in selected
                             for key in ("target_region", "control_region")})
        result[band] = {
            str(region): sum(
                int(row["target_region"]) == region
                or int(row["control_region"]) == region
                for row in selected
            )
            for region in region_ids
        }
    return result


def leave_one_region_out(rows):
    """Recompute partial Pearson and Spearman after paired endpoint deletion."""
    if not rows:
        raise ValueError("fixed association sample is empty")
    full_partial = partial_pearson(rows)
    full_spearman = spearman_rho(rows, "satellite_delta")
    results = []
    regions = sorted({int(row[key]) for row in rows
                      for key in ("target_region", "control_region")})
    for region in regions:
        retained = [row for row in rows
                    if int(row["target_region"]) != region
                    and int(row["control_region"]) != region]
        omitted = [row for row in rows
                   if int(row["target_region"]) == region
                   or int(row["control_region"]) == region]
        record = {
            "region_id": region,
            "n_pairs_omitted": len(omitted),
            "n_pairs_retained": len(retained),
            "omitted_cluster_ids": sorted({str(row["cluster_id"]) for row in omitted}),
            "partial_r_without_region": None,
            "delta_partial_r": None,
            "partial_sign_flip": None,
            "spearman_rho_without_region": None,
            "delta_spearman_rho": None,
            "spearman_sign_flip": None,
            "status": "undefined_after_omission",
        }
        try:
            loo_partial = partial_pearson(retained)
            loo_spearman = spearman_rho(retained, "satellite_delta")
        except (ValueError, np.linalg.LinAlgError):
            results.append(record)
            continue
        record.update({
            "partial_r_without_region": loo_partial,
            "delta_partial_r": loo_partial - full_partial,
            "partial_sign_flip": bool(np.sign(loo_partial) != np.sign(full_partial)),
            "spearman_rho_without_region": loo_spearman,
            "delta_spearman_rho": loo_spearman - full_spearman,
            "spearman_sign_flip": bool(np.sign(loo_spearman) != np.sign(full_spearman)),
            "status": "complete_exploratory",
        })
        results.append(record)
    return {
        "n_pairs": len(rows),
        "full_partial_r": full_partial,
        "full_spearman_rho": full_spearman,
        "leave_one_region_out": results,
    }


def analyze(association, icl, extraction, baseline_path, manifest_path,
            patch_qa_path, output_path, adjacency_deg=2.0):
    output_path = Path(output_path).resolve()
    if output_path.exists():
        raise FileExistsError(f"refusing to overwrite existing output: {output_path}")
    baseline_path = Path(baseline_path).resolve()
    baseline = json.loads(baseline_path.read_text(encoding="utf-8"))
    if baseline.get("schema") != "hsc-satellite-icl-spatial-block-bootstrap-v1" or baseline.get("science_measurement") is not False:
        raise ValueError("unexpected verified 94-pair baseline")
    rows, provenance = load_verified_rows(association, icl, extraction)
    if len(rows) != 94:
        raise ValueError(f"expected fixed 94-row sample; found {len(rows)}")
    for key, value in provenance.items():
        if baseline.get("provenance", {}).get(key) != value:
            raise ValueError(f"baseline provenance mismatch for {key}")

    manifest = json.loads(Path(manifest_path).read_text(encoding="utf-8"))
    lens_rows = manifest.get("clusters", [])
    if len(lens_rows) != 800:
        raise ValueError("expected fixed 800-target sample manifest")
    lens_by_id = {str(row["id"]): row for row in lens_rows}
    if len(lens_by_id) != 800:
        raise ValueError("sample manifest contains duplicate lens IDs")
    if any(str(row["cluster_id"]) not in lens_by_id for row in rows):
        raise ValueError("association row contains a lens absent from the manifest")

    tract_keys, centers = tract_centers_from_patch_qa(Path(patch_qa_path))
    regions = tract_region_labels(centers, adjacency_deg)
    prepared = prepare_rows(rows, centers, regions)
    endpoint_partition_sensitivity = []
    for threshold in (1.6, 1.8, 2.0, 2.2, 2.5, 3.0):
        threshold_regions = tract_region_labels(centers, threshold)
        threshold_rows = prepare_rows(rows, centers, threshold_regions)
        endpoint_partition_sensitivity.append({
            "adjacency_deg": threshold,
            "n_tract_components": int(len(np.unique(threshold_regions))),
            "pairs_touching_region_by_filter": endpoint_touch_counts(threshold_rows),
        })
    reference_partition = endpoint_partition_sensitivity[0]
    if any(item["n_tract_components"] != reference_partition["n_tract_components"]
           or item["pairs_touching_region_by_filter"] != reference_partition["pairs_touching_region_by_filter"]
           for item in endpoint_partition_sensitivity):
        raise ValueError("target/control endpoint partition changes over the adjacency threshold grid")

    baseline_r = {
        row["filter"]: float(row["partial_spatial_block_bootstrap"]["partial_r_point"])
        for row in baseline.get("results", []) if row.get("selection") == "all_supported_pairs"
    }
    if set(baseline_r) != {"i", "i2"}:
        raise ValueError("baseline does not contain both all-supported filters")
    filter_results = []
    baseline_reproduction = {}
    for band in ("i2", "i"):
        selected = [row for row in prepared if row["filter"] == band]
        if len(selected) not in (26, 68):
            raise ValueError(f"unexpected fixed {band} pair count: {len(selected)}")
        full = leave_one_region_out(selected)
        # The full-sample statistic must reproduce the existing 94-pair baseline.
        if not math.isclose(full["full_partial_r"], baseline_r[band], rel_tol=0, abs_tol=1e-12):
            raise ValueError(f"full-sample partial r differs from baseline for {band}")
        baseline_reproduction[band] = {
            "partial_r": full["full_partial_r"],
            "baseline_partial_r": baseline_r[band],
            "status": "exact_match_within_1e-12",
        }
        # Counts refer to matched rows that touch each region at either endpoint.
        full["pairs_touching_region"] = endpoint_touch_counts(selected)[band]
        filter_results.append({"filter": band, **full})

    manifest_regions = nearest_region_labels(
        [row["ra_deg"] for row in lens_rows], [row["dec_deg"] for row in lens_rows],
        centers, regions,
    )
    _, manifest_counts = np.unique(manifest_regions, return_counts=True)
    threshold_sensitivity = []
    for threshold in (1.6, 1.8, 2.0, 2.2, 2.5, 3.0):
        threshold_regions = tract_region_labels(centers, threshold)
        assigned = nearest_region_labels(
            [row["ra_deg"] for row in lens_rows], [row["dec_deg"] for row in lens_rows],
            centers, threshold_regions,
        )
        _, counts = np.unique(assigned, return_counts=True)
        threshold_sensitivity.append({
            "adjacency_deg": threshold,
            "n_tract_components": int(len(np.unique(threshold_regions))),
            "occupied_manifest_lens_counts": sorted(map(int, counts), reverse=True),
        })
    if any(row["n_tract_components"] != threshold_sensitivity[0]["n_tract_components"]
           or row["occupied_manifest_lens_counts"] != threshold_sensitivity[0]["occupied_manifest_lens_counts"]
           for row in threshold_sensitivity):
        raise ValueError("lens-region partition changes over the adjacency threshold grid")

    result = {
        "schema": "astrafek-satellite-icl-tract-region-influence-v1",
        "status": "complete_exploratory_influence",
        "science_measurement": False,
        "covariance_estimated": False,
        "method": {
            "sample": "fixed 94 validated cluster/filter pairs, split by filter",
            "statistic": "partial Pearson r controlling redshift and richness; unadjusted Spearman companion",
            "spatial_deletion": "omit a pair if its target or matched-control coordinate is assigned to the deleted connected tract component",
            "tract_adjacency_deg": float(adjacency_deg),
            "interpretation": "leave-one-region influence only; three occupied regions cannot estimate a stable covariance",
        },
        "n_tracts": len(tract_keys),
        "n_tract_components": int(len(np.unique(regions))),
        "manifest_lens_counts_by_region": {
            str(region): int(count) for region, count in zip(*np.unique(manifest_regions, return_counts=True))
        },
        "pair_endpoint_partition_threshold_sensitivity": endpoint_partition_sensitivity,
        "partition_threshold_sensitivity": threshold_sensitivity,
        "baseline_reproduction": baseline_reproduction,
        "filter_results": filter_results,
        "provenance": provenance,
        "limitations": [
            "tract-center graph is approximate and is not a tract-polygon or official HSC field partition",
            "only three regions contain targets, so deletions are not confidence intervals or a detection test",
            "the fixed 94-pair sample was selected by prior support and extraction gates",
            "mask, sky, PSF, membership completeness, and image covariance limitations remain",
            "no random-lens subtraction, no new HSC query, and no new download were performed",
        ],
        "inputs": {
            "association_dir": str(Path(association).resolve()),
            "icl_dir": str(Path(icl).resolve()),
            "extraction_summary": str(Path(extraction).resolve()),
            "baseline": str(baseline_path),
            "baseline_sha256": sha256(baseline_path),
            "lens_manifest": str(Path(manifest_path).resolve()),
            "lens_manifest_sha256": sha256(Path(manifest_path)),
            "patch_qa": str(Path(patch_qa_path).resolve()),
            "patch_qa_sha256": sha256(Path(patch_qa_path)),
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
    parser.add_argument("--lens-manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--patch-qa", type=Path, default=DEFAULT_PATCH_QA)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--adjacency-deg", type=float, default=2.0)
    args = parser.parse_args()
    result = analyze(args.association, args.icl, args.extraction_summary,
                     args.baseline, args.lens_manifest, args.patch_qa,
                     args.output, args.adjacency_deg)
    print(json.dumps({
        "status": result["status"], "science_measurement": False,
        "baseline_reproduction": result["baseline_reproduction"],
        "filter_results": [{
            "filter": row["filter"], "n_pairs": row["n_pairs"],
            "full_partial_r": row["full_partial_r"],
            "leave_one_region_out": row["leave_one_region_out"],
        } for row in result["filter_results"]],
    }, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
