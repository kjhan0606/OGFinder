#!/usr/bin/env python3
"""Measure radial-cut sensitivity of saved ICL-split profiles.

For each approximate block width, re-evaluates the conditional tangential and
cross zero-profile Wald diagnostics after successively removing inner radial
bins. This is a post-hoc diagnostic motivated by the source/random cross
decomposition; its nominal p-values are not calibrated detection statistics.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.stats import chi2


DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
DEFAULT_INPUTS = {
    "i": {
        2.5: DATA / "wl_random_subtracted_icl_split_i_2026-11-06_02/summary.json",
        5.0: DATA / "wl_random_subtracted_icl_split_i_2026-11-06_01/summary.json",
        10.0: DATA / "wl_random_subtracted_icl_split_i_2026-11-06_03/summary.json",
    },
    "i2": {
        2.5: DATA / "wl_random_subtracted_icl_split_2026-11-05_03/summary.json",
        5.0: DATA / "wl_random_subtracted_icl_split_2026-11-05_02/summary.json",
    },
}
DEFAULT_OUTPUT = DATA / "wl_icl_split_radial_cut_nulls_2026-11-06_02"
RADIAL_CUTS_MPC = (0.1, 0.2, 0.4, 0.7, 1.0)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def sliced_null(profile: list[dict], covariance: np.ndarray,
                radial_cut_mpc: float, component: str) -> dict:
    if component not in {"tangential", "cross"}:
        raise ValueError(component)
    bins = [i for i, row in enumerate(profile)
            if float(row["r_min_mpc"]) >= radial_cut_mpc]
    if not bins:
        raise ValueError(f"radial cut {radial_cut_mpc} leaves no bins")
    offset = 0 if component == "tangential" else 1
    index = np.asarray([2 * i + offset for i in bins], dtype=int)
    field = ("high_minus_low_delta_sigma_msun_mpc2" if component == "tangential"
             else "high_minus_low_cross_msun_mpc2")
    values = np.asarray([float(profile[i][field]) for i in bins])
    sub_cov = covariance[np.ix_(index, index)]
    eigenvalues = np.linalg.eigvalsh(sub_cov)
    scale = max(float(np.max(np.abs(eigenvalues))), 1.0)
    rank = int(np.count_nonzero(eigenvalues > scale * 1e-10))
    if rank == 0 or not np.all(np.isfinite(values)):
        raise ValueError("invalid profile or singular covariance subset")
    statistic = float(values @ np.linalg.pinv(sub_cov, rcond=1e-10) @ values)
    return {
        "component": component,
        "radial_cut_mpc": float(radial_cut_mpc),
        "included_bins": [
            {"r_min_mpc": float(profile[i]["r_min_mpc"]),
             "r_max_mpc": float(profile[i]["r_max_mpc"]),
             "high_minus_low_msun_mpc2": float(values[j])}
            for j, i in enumerate(bins)
        ],
        "n_bins": len(bins),
        "covariance_rank": rank,
        "chi2": statistic,
        "nominal_p_conditional_on_saved_covariance": float(chi2.sf(statistic, rank)),
        "interpretation": "post-hoc conditional covariance diagnostic; not a calibrated detection significance",
    }


def run(output: Path = DEFAULT_OUTPUT) -> dict:
    if output.exists():
        raise FileExistsError(f"refusing to overwrite radial-cut assessment: {output}")
    by_filter = {}
    for filter_name, inputs in DEFAULT_INPUTS.items():
        by_width = {}
        for width, path in inputs.items():
            doc = json.loads(path.read_text(encoding="utf-8"))
            if (doc.get("status") != "complete_exploratory"
                    or doc.get("science_measurement") is not False
                    or doc.get("filter") != filter_name
                    or float(doc.get("block_width_deg")) != width):
                raise ValueError(f"unexpected {filter_name}-band saved profile at {width} degrees")
            profile = doc["profile"]
            covariance = np.asarray(doc["high_minus_low_covariance"], dtype=float)
            if len(profile) != 7 or covariance.shape != (14, 14):
                raise ValueError("saved split profile/covariance shape changed")
            by_width[str(width)] = {
                "n_matched_pairs": int(doc["n_matched_lens_pairs"]),
                "n_shared_delete_one_blocks": int(doc["n_shared_delete_one_blocks"]),
                "summary_sha256": sha256(path),
                "radial_cuts": {
                    str(cut): {
                        comp: sliced_null(profile, covariance, cut, comp)
                        for comp in ("tangential", "cross")
                    }
                    for cut in RADIAL_CUTS_MPC
                },
            }
        by_filter[filter_name] = by_width

    result = {
        "schema": "astrafek-wl-icl-split-radial-cut-null-audit-v2",
        "status": "complete_posthoc_diagnostic",
        "science_measurement": False,
        "network_used": False,
        "images_read": False,
        "shape_rows_read": False,
        "method": "slice the saved 14-dimensional high-minus-low paired endpoint covariance to radial subsets and recompute conditional zero-profile Wald diagnostics",
        "radial_cuts_mpc": list(RADIAL_CUTS_MPC),
        "filter_results": by_filter,
        "limitations": [
            "the 1 Mpc cut was examined after seeing the inner cross-component structure",
            "the 10-pair sample and 9–24 approximate blocks do not provide a calibrated field covariance",
            "nominal p-values are conditional diagnostics, not detection significances",
            "no photo-z, shear, depth, mask, PSF, or source-Poisson uncertainty is propagated",
        ],
    }
    output.mkdir(parents=True)
    report = output / "summary.json"
    result["analysis_script"] = str(Path(__file__).resolve())
    result["analysis_script_sha256"] = sha256(Path(__file__).resolve())
    report.write_text(json.dumps(result, indent=2, ensure_ascii=False,
                                 allow_nan=False) + "\n", encoding="utf-8")
    result["summary_sha256"] = sha256(report)
    print(json.dumps({
        "status": result["status"],
        "science_measurement": False,
        "conditional_p_by_filter_width_cut": {
            filter_name: {
                width: {
                    cut: {comp: round(test["nominal_p_conditional_on_saved_covariance"], 5)
                          for comp, test in comps.items()}
                    for cut, comps in doc["radial_cuts"].items()
                }
                for width, doc in widths.items()
            }
            for filter_name, widths in by_filter.items()
        },
        "summary": str(report),
        "summary_sha256": result["summary_sha256"],
    }, ensure_ascii=False))
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    run(args.output)


if __name__ == "__main__":
    main()
