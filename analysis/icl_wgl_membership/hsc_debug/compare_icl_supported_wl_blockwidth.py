#!/usr/bin/env python3
"""Compare ICL-supported WL point estimates and covariance across sky blocks."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.stats import chi2


DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
STACK_5 = DATA / "icl_supported_wl_subset_stack_2026-11-05_01/summary.json"
STACK_2P5 = DATA / "icl_supported_wl_subset_stack_2026-11-05_02p5/summary.json"
ASSESSMENT_5 = DATA / "icl_supported_wl_subset_assessment_2026-11-05_01/assessment.json"
ASSESSMENT_2P5 = DATA / "icl_supported_wl_subset_assessment_2026-11-05_02p5/assessment.json"
OUTPUT = DATA / "icl_supported_wl_blockwidth_comparison_2026-11-05_01"


def sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def null_summary(stack_filter: dict, field: str, covariance_indices: np.ndarray) -> dict:
    jk = stack_filter["paired_endpoint_jackknife"]
    values = np.asarray([row[field] for row in jk["profile"]], dtype=float)
    covariance = np.asarray(jk["covariance"], dtype=float)
    sub = covariance[np.ix_(covariance_indices, covariance_indices)]
    eig = np.linalg.eigvalsh(sub)
    scale = max(float(np.max(np.abs(eig))), 1.0)
    rank = int(np.count_nonzero(eig > scale * 1e-10))
    q = float(values @ np.linalg.pinv(sub, rcond=1e-10) @ values)
    return {"chi_square": q, "rank": rank, "p_value": float(chi2.sf(q, rank))}


def run(stack_5_path: Path, stack_2p5_path: Path,
        assessment_5_path: Path, assessment_2p5_path: Path, output: Path) -> dict:
    if output.exists():
        raise FileExistsError(f"refusing to overwrite block-width comparison: {output}")
    stack_5, stack_2p5 = (json.loads(path.read_text(encoding="utf-8"))
                          for path in (stack_5_path, stack_2p5_path))
    assessment_5, assessment_2p5 = (json.loads(path.read_text(encoding="utf-8"))
                                    for path in (assessment_5_path, assessment_2p5_path))
    if any(doc.get("science_measurement") is not False
           for doc in (stack_5, stack_2p5, assessment_5, assessment_2p5)):
        raise ValueError("one input has an unexpected science gate")
    if stack_5["radial_edges_mpc"] != stack_2p5["radial_edges_mpc"]:
        raise ValueError("radial bins differ between the block-width runs")

    by_filter = {}
    for filt in ("i", "i2"):
        a, b = stack_5["by_filter"][filt], stack_2p5["by_filter"][filt]
        if (a["lens_ids"] != b["lens_ids"]
                or a["n_target_shape_pairs"] != b["n_target_shape_pairs"]
                or a["n_random_shape_pairs"] != b["n_random_shape_pairs"]):
            raise ValueError(f"{filt} sample or pair counts changed across block widths")
        prof_a = a["paired_endpoint_jackknife"]["profile"]
        prof_b = b["paired_endpoint_jackknife"]["profile"]
        if len(prof_a) != len(prof_b):
            raise ValueError(f"{filt} radial profile length changed")
        profile = []
        for index, (r5, r25) in enumerate(zip(prof_a, prof_b)):
            for key in ("delta_sigma_random_subtracted_msun_mpc2",
                        "cross_delta_sigma_random_subtracted_msun_mpc2"):
                if not np.isclose(float(r5[key]), float(r25[key]), rtol=1e-10, atol=1e-3):
                    raise ValueError(f"{filt} point estimate changed at bin {index}: {key}")
            profile.append({
                "bin": index,
                "r_min_mpc": r5["r_min_mpc"],
                "r_max_mpc": r5["r_max_mpc"],
                "delta_sigma": r5["delta_sigma_random_subtracted_msun_mpc2"],
                "delta_sigma_error_5deg": r5["delta_sigma_jackknife_error_msun_mpc2"],
                "delta_sigma_error_2p5deg": r25["delta_sigma_jackknife_error_msun_mpc2"],
                "delta_sigma_error_ratio_2p5_over_5": (
                    r25["delta_sigma_jackknife_error_msun_mpc2"] /
                    r5["delta_sigma_jackknife_error_msun_mpc2"]
                ),
                "cross_delta_sigma": r5["cross_delta_sigma_random_subtracted_msun_mpc2"],
                "cross_error_5deg": r5["cross_jackknife_error_msun_mpc2"],
                "cross_error_2p5deg": r25["cross_jackknife_error_msun_mpc2"],
            })
        tangential_idx = np.arange(0, 14, 2)
        cross_idx = np.arange(1, 14, 2)
        nulls = {
            "5deg": {
                "tangential": null_summary(a, "delta_sigma_random_subtracted_msun_mpc2", tangential_idx),
                "cross": null_summary(a, "cross_delta_sigma_random_subtracted_msun_mpc2", cross_idx),
            },
            "2p5deg": {
                "tangential": null_summary(b, "delta_sigma_random_subtracted_msun_mpc2", tangential_idx),
                "cross": null_summary(b, "cross_delta_sigma_random_subtracted_msun_mpc2", cross_idx),
            },
        }
        ratios = [row["delta_sigma_error_ratio_2p5_over_5"] for row in profile]
        cross_bad = [nulls[width]["cross"]["p_value"] < 0.05 for width in ("5deg", "2p5deg")]
        by_filter[filt] = {
            "n_lenses": a["n_icl_profile_supported_lenses"],
            "n_target_pairs": a["n_target_shape_pairs"],
            "n_random_pairs": a["n_random_shape_pairs"],
            "n_blocks": {"5deg": a["paired_endpoint_jackknife"]["n_occupied_blocks"],
                         "2p5deg": b["paired_endpoint_jackknife"]["n_occupied_blocks"]},
            "point_estimates_identical": True,
            "tangential_error_ratio_2p5_over_5": {
                "min": float(min(ratios)), "median": float(np.median(ratios)),
                "max": float(max(ratios)),
            },
            "null_tests": nulls,
            "cross_null_rejected_at_both_widths": bool(all(cross_bad)),
            "profile": profile,
        }
    result = {
        "schema": "astrafek-icl-supported-wl-blockwidth-comparison-v1",
        "status": "complete_exploratory",
        "science_measurement": False,
        "block_widths_deg": [5.0, 2.5],
        "by_filter": by_filter,
        "input_hashes": {
            "stack_5deg_sha256": sha(stack_5_path),
            "stack_2p5deg_sha256": sha(stack_2p5_path),
            "assessment_5deg_sha256": sha(assessment_5_path),
            "assessment_2p5deg_sha256": sha(assessment_2p5_path),
        },
        "interpretation": {
            "i": "Tangential point estimates are unchanged, but the cross null is rejected under both covariance partitions. Treat this as a systematic/covariance failure warning; do not interpret the tangential profile as a detection.",
            "i2": "Tangential point estimates are unchanged and the cross null is consistent under both partitions. This is a stronger exploratory lensing result, still conditional on unresolved survey-region, random-aperture, photo-z, boost, and reduced-shear systematics.",
        },
        "science_gates_remaining": [
            "validated contiguous HSC jackknife regions and block-width stability beyond this two-width check",
            "full-aperture random-control matching in depth, masks, PSF, and sky residuals",
            "photo-z dilution and calibration uncertainty propagation",
            "source-density boost with mask and obscuration correction",
            "reduced-shear correction and coadd pixel-covariance validation",
        ],
    }
    output.mkdir(parents=True)
    path = output / "comparison.json"
    path.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": result["status"], "science_measurement": False,
        "by_filter": {
            filt: {
                "n_lenses": item["n_lenses"],
                "tangential_p_5deg": item["null_tests"]["5deg"]["tangential"]["p_value"],
                "tangential_p_2p5deg": item["null_tests"]["2p5deg"]["tangential"]["p_value"],
                "cross_p_5deg": item["null_tests"]["5deg"]["cross"]["p_value"],
                "cross_p_2p5deg": item["null_tests"]["2p5deg"]["cross"]["p_value"],
                "median_error_ratio": item["tangential_error_ratio_2p5_over_5"]["median"],
            } for filt, item in by_filter.items()
        },
        "output": str(path),
    }, ensure_ascii=False))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stack-5deg", type=Path, default=STACK_5)
    parser.add_argument("--stack-2p5deg", type=Path, default=STACK_2P5)
    parser.add_argument("--assessment-5deg", type=Path, default=ASSESSMENT_5)
    parser.add_argument("--assessment-2p5deg", type=Path, default=ASSESSMENT_2P5)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    run(args.stack_5deg, args.stack_2p5deg, args.assessment_5deg,
        args.assessment_2p5deg, args.output)
