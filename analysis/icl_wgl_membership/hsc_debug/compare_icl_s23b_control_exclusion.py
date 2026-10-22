#!/usr/bin/env python3
"""Compare verified baseline and S23B-control-excluded ICL radial stacks.

This is a paired-sample sensitivity description, not a significance test:
the alternate stack is a subset of the baseline and their covariance is not
estimated here.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def read_profile(path: Path) -> dict[tuple[str, float, float], dict]:
    result = {}
    with path.open(encoding="utf-8", newline="") as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            key = (row["filter"], float(row["r_min_kpc"]), float(row["r_max_kpc"]))
            if key in result:
                raise ValueError(f"duplicate radial profile bin {key}")
            result[key] = row
    return result


def area_weighted_mean(rows: list[dict], column: str,
                       r_min: float, r_max: float) -> float:
    selected = [row for row in rows
                if float(row["r_min_kpc"]) >= r_min
                and float(row["r_max_kpc"]) <= r_max]
    if not selected:
        raise ValueError(f"no radial bins within {r_min}-{r_max} kpc")
    ordered = sorted(selected, key=lambda row: float(row["r_min_kpc"]))
    edges = [(float(row["r_min_kpc"]), float(row["r_max_kpc"])) for row in ordered]
    if edges[0][0] != r_min or edges[-1][1] != r_max or any(
            edges[i][1] != edges[i + 1][0] for i in range(len(edges) - 1)):
        raise ValueError(f"radial bins do not exactly tile {r_min}-{r_max} kpc")
    values, weights = [], []
    for row, (lo, hi) in zip(ordered, edges):
        value = float(row[column])
        if not math.isfinite(value):
            raise ValueError(f"non-finite {column} in radial range {r_min}-{r_max}")
        values.append(value)
        weights.append(hi * hi - lo * lo)
    return sum(v * w for v, w in zip(values, weights)) / sum(weights)


def verified_summary(work: Path, verification: Path) -> tuple[dict, dict, dict]:
    summary_path = work / "summary.json"
    summary = read_json(summary_path)
    verification_doc = read_json(verification)
    if (summary.get("schema") != "astrafek-hsc-icl-target-control-2d-stack-v1"
            or summary.get("status") != "complete_exploratory"
            or summary.get("science_measurement") is not False
            or verification_doc.get("status") != "verified"
            or verification_doc.get("science_measurement") is not False
            or verification_doc.get("summary_sha256") != sha256(summary_path)
            or verification_doc.get("mask_factor") != 2.5):
        raise ValueError(f"stack or independent verification is invalid: {work}")
    radial_path = Path(summary["outputs"]["stacked_profiles"])
    if sha256(radial_path) != summary["outputs"]["stacked_profiles_sha256"]:
        raise ValueError(f"radial profile checksum mismatch: {radial_path}")
    return summary, verification_doc, read_profile(radial_path)


def compare(baseline_work: Path, baseline_verification: Path,
            alternate_work: Path, alternate_verification: Path) -> dict:
    base, base_v, base_rows = verified_summary(baseline_work, baseline_verification)
    alt, alt_v, alt_rows = verified_summary(alternate_work, alternate_verification)
    if (base_v.get("mask_factor") != alt_v.get("mask_factor")
            or base_v.get("mask_factor") != 2.5
            or base["grid"] != alt["grid"]
            or base.get("n_pairs_stacked") != 783
            or alt.get("n_pairs_stacked") != 777
            or alt.get("n_pairs_excluded_s23b_candidate_control") != 6
            or alt.get("n_candidate_controls_already_filter_excluded") != 1
            or alt.get("candidate_control_audit") is None):
        raise ValueError("baseline/alternate pair counts, masks, grid, or audit provenance differ from the expected sensitivity")
    expected_counts = dict(base["n_by_filter"])
    for exclusion in alt.get("s23b_candidate_control_exclusions", []):
        expected_counts[exclusion["filter"]] -= 1
    if alt.get("n_by_filter") != expected_counts:
        raise ValueError("alternate per-filter pair counts do not match the candidate-control exclusions")
    if base_v.get("n_eligible_pairs") != 783 or alt_v.get("n_eligible_pairs") != 777:
        raise ValueError("independent verification pair counts are inconsistent")
    if set(base_rows) != set(alt_rows) or len(base_rows) != 40:
        raise ValueError("baseline and alternate radial bins do not match")

    per_bin, interval_summary = [], []
    for band in ("i", "i2"):
        band_rows_base, band_rows_alt = [], []
        for (row_band, rmin, rmax), base_row in sorted(base_rows.items()):
            if row_band != band:
                continue
            alt_row = alt_rows[(row_band, rmin, rmax)]
            b = float(base_row["mean_delta_njy_arcsec2"])
            a = float(alt_row["mean_delta_njy_arcsec2"])
            change = a - b
            per_bin.append({"filter": band, "r_min_kpc": rmin, "r_max_kpc": rmax,
                            "baseline_delta_njy_arcsec2": b,
                            "alternate_delta_njy_arcsec2": a,
                            "change_njy_arcsec2": change,
                            "baseline_n_pairs": int(base_row["n_pairs_delta"]),
                            "alternate_n_pairs": int(alt_row["n_pairs_delta"])})
            band_rows_base.append(base_row)
            band_rows_alt.append(alt_row)
        for lo, hi in ((50.0, 200.0), (200.0, 300.0)):
            bmean = area_weighted_mean(band_rows_base, "mean_delta_njy_arcsec2", lo, hi)
            amean = area_weighted_mean(band_rows_alt, "mean_delta_njy_arcsec2", lo, hi)
            diff = amean - bmean
            interval_summary.append({"filter": band, "r_min_kpc": lo, "r_max_kpc": hi,
                                     "area_weighted_baseline_delta_njy_arcsec2": bmean,
                                     "area_weighted_alternate_delta_njy_arcsec2": amean,
                                     "change_njy_arcsec2": diff,
                                     "relative_change_percent": (100.0 * diff / bmean
                                         if abs(bmean) > 1e-12 else None)})
    return {
        "schema": "astrafek-icl-s23b-control-exclusion-comparison-v1",
        "status": "verified_exploratory_sensitivity",
        "science_measurement": False,
        "baseline": {"work": str(baseline_work.resolve()),
                     "summary_sha256": sha256(baseline_work / "summary.json"),
                     "verification_sha256": sha256(baseline_verification),
                     "n_pairs": base["n_pairs_stacked"], "n_by_filter": base["n_by_filter"]},
        "alternate": {"work": str(alternate_work.resolve()),
                      "summary_sha256": sha256(alternate_work / "summary.json"),
                      "verification_sha256": sha256(alternate_verification),
                      "n_pairs": alt["n_pairs_stacked"], "n_by_filter": alt["n_by_filter"],
                      "candidate_control_audit": alt["candidate_control_audit"]},
        "intervals": interval_summary,
        "per_bin": per_bin,
        "interpretation_limits": [
            "The alternate stack is a subset of the baseline; this comparison is descriptive and does not estimate the covariance of their difference.",
            "Area weighting combines annular surface-brightness means within each requested radial interval.",
            "Candidate centers in control images flag possible cluster-light contamination only; they do not establish detected diffuse-light bias.",
            "The result remains exploratory because extinction, sky-transfer, PSF-wing, correlated-noise, and full spatial-covariance corrections are incomplete.",
        ],
    }


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline-work", type=Path, required=True)
    parser.add_argument("--baseline-verification", type=Path, required=True)
    parser.add_argument("--alternate-work", type=Path, required=True)
    parser.add_argument("--alternate-verification", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args(argv)
    output = args.output.resolve()
    if output.exists():
        parser.error(f"refusing to overwrite {output}")
    result = compare(args.baseline_work.resolve(), args.baseline_verification.resolve(),
                     args.alternate_work.resolve(), args.alternate_verification.resolve())
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n",
                      encoding="utf-8")
    print(json.dumps({"status": result["status"], "intervals": result["intervals"],
                      "output": str(output)}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
