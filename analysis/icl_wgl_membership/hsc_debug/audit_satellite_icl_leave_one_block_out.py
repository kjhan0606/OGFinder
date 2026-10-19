#!/usr/bin/env python3
"""Measure how much one connected sky block can influence ICL--satellite r."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from audit_satellite_icl_perbin_block_scale import (
    DEFAULT_ASSOC, DEFAULT_BASELINE, DEFAULT_EXTRACTION, DEFAULT_ICL,
    partial_pearson,
)
from audit_satellite_icl_perbin_block_scale import sha256 as sha
from audit_satellite_icl_spatial_blocks import load_verified_rows
from audit_satellite_icl_block_scale_sensitivity import connected_pair_blocks


def write_tsv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError("refusing to write an empty leave-one-block table")
    with path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t",
                                extrasaction="raise", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def analyze(association: Path, icl: Path, extraction: Path, baseline: Path,
            output: Path, widths=(2.5, 5.0, 10.0)) -> dict:
    output = Path(output).resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite existing output: {output}")
    widths = tuple(float(w) for w in widths)
    if not widths or len(set(widths)) != len(widths) or any(
        not math.isfinite(w) or w <= 0 or w > 180 for w in widths
    ):
        raise ValueError("block widths must be unique and lie in (0, 180]")
    baseline = Path(baseline).resolve()
    base = json.loads(baseline.read_text(encoding="utf-8"))
    if base.get("schema") != "hsc-satellite-icl-spatial-block-bootstrap-v1" or base.get("science_measurement") is not False:
        raise ValueError("unexpected 94-row spatial-block baseline")
    rows, provenance = load_verified_rows(association, icl, extraction)
    if len(rows) != 94:
        raise ValueError(f"expected the fixed 94-row sample; found {len(rows)}")
    for key, value in provenance.items():
        if base.get("provenance", {}).get(key) != value:
            raise ValueError(f"baseline provenance mismatch for {key}")
    baseline_r = {
        row["filter"]: float(row["partial_spatial_block_bootstrap"]["partial_r_point"])
        for row in base.get("results", [])
        if row.get("selection") == "all_supported_pairs"
    }
    if set(baseline_r) != {"i", "i2"}:
        raise ValueError("baseline lacks the two all-supported filter estimates")

    output_rows = []
    summaries = []
    baseline_reproduction = {}
    bands = ("i2", "i")
    for width in widths:
        for band in bands:
            selected = [dict(row, z=row["redshift"],
                             satellite_delta=row["satellite_delta_50_75"])
                        for row in rows if row["filter"] == band]
            labels = connected_pair_blocks(selected, width)
            blocks = np.unique(labels)
            full_r = partial_pearson(selected)
            if band not in baseline_reproduction:
                matches_baseline = math.isclose(full_r, baseline_r[band],
                                                rel_tol=0, abs_tol=1e-12)
                if not matches_baseline:
                    raise ValueError(f"full-sample r differs from verified baseline for {band}")
                baseline_reproduction[band] = {
                    "partial_r": full_r,
                    "baseline_partial_r": baseline_r[band],
                    "status": "exact_match_within_1e-12",
                }
            by_block = {str(block): [row for row, label in zip(selected, labels)
                                     if label == block]
                        for block in blocks}
            influences = []
            for block in sorted(by_block):
                omitted = by_block[block]
                retained = [row for row, label in zip(selected, labels)
                            if str(label) != block]
                record = {
                    "block_width_deg": width,
                    "filter": band,
                    "omitted_block": block,
                    "omitted_cluster_ids": ";".join(sorted({str(row["cluster_id"])
                                                                for row in omitted})),
                    "n_omitted_pairs": len(omitted),
                    "n_retained_pairs": len(retained),
                    "full_partial_r": full_r,
                    "partial_r_without_block": None,
                    "delta_r": None,
                    "status": "undefined_after_omission",
                }
                try:
                    loo_r = partial_pearson(retained)
                except (ValueError, np.linalg.LinAlgError):
                    pass
                else:
                    delta = loo_r - full_r
                    record.update(partial_r_without_block=loo_r, delta_r=delta,
                                  status="complete_exploratory")
                    influences.append((block, len(omitted), loo_r, delta))
                output_rows.append(record)
            if not influences:
                raise ValueError(f"all leave-one-block estimates undefined for {band}/{width}")
            loo_values = [item[2] for item in influences]
            largest = max(influences, key=lambda item: abs(item[3]))
            summaries.append({
                "block_width_deg": width,
                "filter": band,
                "n_pairs": len(selected),
                "n_blocks": len(blocks),
                "largest_block_pairs": max(len(group) for group in by_block.values()),
                "full_partial_r": full_r,
                "leave_one_block_r_min": min(loo_values),
                "leave_one_block_r_max": max(loo_values),
                "sign_reversals": int(sum(np.sign(value) != np.sign(full_r)
                                           for value in loo_values)),
                "valid_leave_one_out": len(influences),
                "undefined_leave_one_out": len(blocks) - len(influences),
                "max_abs_delta_r": abs(largest[3]),
                "most_influential_block": largest[0],
                "most_influential_cluster_ids": ";".join(sorted({
                    str(row["cluster_id"]) for row in by_block[largest[0]]})),
                "most_influential_block_pairs": largest[1],
                "partial_r_without_most_influential_block": largest[2],
            })

    output.mkdir(parents=True, exist_ok=False)
    table = output / "leave_one_block_out.tsv"
    write_tsv(table, output_rows)
    script = Path(__file__).resolve()
    result = {
        "schema": "hsc-satellite-icl-leave-one-block-out-v1",
        "status": "complete_exploratory",
        "science_measurement": False,
        "method": {
            "sample": "fixed 94 per-annulus cluster/filter rows, split by filter",
            "statistic": "partial Pearson r of ICL 200-300 kpc and satellite 50-75 kpc, controlling redshift and richness",
            "spatial_groups": "connected target/control RA-Dec cells using the archived block helper",
            "leave_one_out": "omit each connected block once and refit nuisance residuals",
            "block_widths_deg": list(widths),
            "interpretation": "influence diagnostic; ranges are not confidence intervals",
        },
        "inputs": {
            "association": str(Path(association).resolve()),
            "icl": str(Path(icl).resolve()),
            "extraction_summary": str(Path(extraction).resolve()),
            "baseline": str(baseline),
            "baseline_sha256": sha(baseline),
            "verified_provenance": provenance,
        },
        "full_sample_baseline_reproduction": baseline_reproduction,
        "summary_rows": summaries,
        "results_tsv": str(table),
        "results_tsv_sha256": sha(table),
        "analysis_script": str(script),
        "analysis_script_sha256": sha(script),
        "limitations": [
            "The fixed sample is selected by common-support and extraction-quality gates.",
            "Rectangular RA/Dec block widths are post hoc sensitivity choices and do not capture all long-range covariance.",
            "Leave-one-block ranges are not sampling intervals or significance tests; there is no multiplicity correction.",
            "Sky transfer, extended PSF wings, membership completeness, extinction, and full coadd covariance remain unresolved.",
        ],
    }
    summary = output / "summary.json"
    summary.write_text(json.dumps(result, ensure_ascii=False, indent=2,
                                  allow_nan=False) + "\n", encoding="utf-8")
    result["summary_sha256"] = sha(summary)
    print(json.dumps({"status": result["status"], "science_measurement": False,
                      "summary_rows": summaries,
                      "leave_one_block_out_sha256": result["results_tsv_sha256"],
                      "summary_sha256": result["summary_sha256"]},
                     ensure_ascii=False, indent=2))
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--association", type=Path, default=DEFAULT_ASSOC)
    parser.add_argument("--icl", type=Path, default=DEFAULT_ICL)
    parser.add_argument("--extraction-summary", type=Path, default=DEFAULT_EXTRACTION)
    parser.add_argument("--baseline", type=Path, default=DEFAULT_BASELINE)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--widths-deg", type=float, nargs="+", default=[2.5, 5.0, 10.0])
    args = parser.parse_args()
    analyze(args.association, args.icl, args.extraction_summary, args.baseline,
            args.output, args.widths_deg)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
