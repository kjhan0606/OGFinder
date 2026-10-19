#!/usr/bin/env python3
"""Test ICL profiles after a paired 400--500 kpc outer-annulus baseline."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np

from evaluate_prestack_pair_support import pair_block_labels


ROOT = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
DEFAULT_STACK = ROOT / "icl_2d_matched_stack_2026-10-19_02"
DEFAULT_EXTRACTION = ROOT / "prestack_extraction_2026-10-15_01/extraction_summary.json"
ANNULI = ((50.0, 100.0), (100.0, 200.0), (200.0, 300.0))
OUTER = (400.0, 500.0)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def weighted_annulus(profile: list[dict], r_min: float, r_max: float) -> tuple[float | None, int, float]:
    """Combine saved sub-annuli by valid-pixel count; return mean, pixels, coverage."""
    selected = [row for row in profile
                if float(row["r_min_kpc"]) >= r_min
                and float(row["r_max_kpc"]) <= r_max]
    if not selected:
        raise ValueError(f"no saved radial bins span {r_min:g}--{r_max:g} kpc")
    geometric = sum(int(row["geometric_pixels"]) for row in selected)
    valid = sum(int(row["n_valid_pixels"]) for row in selected)
    values = [(float(row["delta_mean_njy_arcsec2"]), int(row["n_valid_pixels"]))
              for row in selected
              if row["delta_mean_njy_arcsec2"] is not None
              and int(row["n_valid_pixels"]) > 0]
    if valid == 0:
        return None, 0, 0.0
    if sum(weight for _, weight in values) != valid:
        raise ValueError("radial value count and valid-pixel count disagree")
    mean = sum(value * weight for value, weight in values) / valid
    return float(mean), valid, valid / geometric


def percentile_interval(values: np.ndarray) -> list[float | None]:
    finite = np.asarray(values, dtype=float)
    finite = finite[np.isfinite(finite)]
    if len(finite) == 0:
        return [None, None]
    return [float(np.quantile(finite, 0.025)), float(np.quantile(finite, 0.975))]


def render_report(result: dict) -> str:
    lines = [
        "# ICL outer-annulus baseline sensitivity — 2026-10-20",
        "",
        "This is an exploratory comparison of each matched target-control pair's",
        "inner annular mean with the same pair's 400-500 kpc mean. The outer",
        "annulus is not assumed to be blank sky: subtracting it measures inner",
        "excess relative to that region and may remove real extended ICL.",
        "",
        "All profiles require at least a 5% common-pixel area fraction in the",
        "inner and outer annuli. Pair means are weighted equally. Intervals are",
        "95% percentile ranges from resampling connected 5-degree target/control",
        "blocks; they are descriptive and omit residual-sky, PSF-wing, sample",
        "selection, and full coadd-covariance uncertainties.",
        "",
        "| Filter | Annulus (kpc) | Pairs / blocks | Raw ΔSB | 400–500 kpc | Outer-subtracted |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for band in ("i", "i2"):
        for name, row in result["by_filter"][band].items():
            lines.append(
                f"| {band} | {name} | {row['n_pairs']} / {row['n_connected_blocks']} | "
                f"{row['raw_mean_njy_arcsec2']:.3f} "
                f"[{row['raw_block_bootstrap_95_interval'][0]:.3f}, {row['raw_block_bootstrap_95_interval'][1]:.3f}] | "
                f"{row['outer_400_500_mean_njy_arcsec2']:.3f} "
                f"[{row['outer_400_500_block_bootstrap_95_interval'][0]:.3f}, {row['outer_400_500_block_bootstrap_95_interval'][1]:.3f}] | "
                f"{row['outer_subtracted_mean_njy_arcsec2']:.3f} "
                f"[{row['outer_subtracted_block_bootstrap_95_interval'][0]:.3f}, {row['outer_subtracted_block_bootstrap_95_interval'][1]:.3f}] |"
            )
    lines.extend([
        "",
        "The i-band outer baseline is positive in this support-selected analysis,",
        "while the i2 estimate is near zero. This quantifies the known",
        "filter-dependent large-radius offset but does not determine whether it",
        "comes from residual sky or genuine extended light. No measurement is",
        "corrected or promoted to a detection (`science_measurement=false`).",
        "",
        f"Replicates: {result['bootstrap']['replicas']}; seed base: {result['bootstrap']['seed']}",
        f"Input profile SHA-256: `{result['inputs']['paired_profiles_sha256']}`",
        f"Output table SHA-256: `{result['output_sha256']}`",
    ])
    return "\n".join(lines) + "\n"


def block_bootstrap(records: list[dict], replicas: int, seed: int) -> dict:
    if replicas < 100:
        raise ValueError("at least 100 connected-block bootstrap replicates are required")
    target_rows = [{"ra_deg": float(row["ra_deg"]), "dec_deg": float(row["dec_deg"])}
                   for row in records]
    control_rows = [{"ra_deg": float(row["control_ra_deg"]),
                     "dec_deg": float(row["control_dec_deg"])} for row in records]
    _, labels, n_target_cells, n_control_cells, n_blocks, largest = pair_block_labels(
        target_rows, control_rows, np.ones(len(records), dtype=bool))
    members = [np.flatnonzero(labels == label) for label in np.unique(labels)]
    if n_blocks < 2:
        raise ValueError("fewer than two connected target/control blocks")
    raw = np.asarray([row["inner_delta"] for row in records], dtype=float)
    baseline = np.asarray([row["outer_delta"] for row in records], dtype=float)
    corrected = raw - baseline
    rng = np.random.default_rng(seed)
    draws = np.empty((replicas, 3), dtype=float)
    for i in range(replicas):
        selected = rng.integers(0, len(members), size=len(members))
        index = np.concatenate([members[int(j)] for j in selected])
        draws[i] = (np.mean(raw[index]), np.mean(baseline[index]), np.mean(corrected[index]))
    return {
        "n_pairs": len(records),
        "n_target_cells": n_target_cells,
        "n_control_cells": n_control_cells,
        "n_connected_blocks": n_blocks,
        "largest_block_pairs": largest,
        "raw_mean_njy_arcsec2": float(np.mean(raw)),
        "raw_block_bootstrap_95_interval": percentile_interval(draws[:, 0]),
        "outer_400_500_mean_njy_arcsec2": float(np.mean(baseline)),
        "outer_400_500_block_bootstrap_95_interval": percentile_interval(draws[:, 1]),
        "outer_subtracted_mean_njy_arcsec2": float(np.mean(corrected)),
        "outer_subtracted_block_bootstrap_95_interval": percentile_interval(draws[:, 2]),
        "paired_outer_subtraction_shift_njy_arcsec2": float(np.mean(corrected - raw)),
        "median_inner_coverage": float(np.median([row["inner_coverage"] for row in records])),
        "median_outer_coverage": float(np.median([row["outer_coverage"] for row in records])),
        "bootstrap_replicates": replicas,
        "bootstrap_seed": seed,
    }


def analyze(stack_dir: Path, extraction_summary_path: Path, output_dir: Path,
            min_coverage: float = 0.05, replicas: int = 5000, seed: int = 20261020) -> dict:
    if output_dir.exists():
        raise FileExistsError(f"refusing to overwrite {output_dir}")
    summary_path = stack_dir / "summary.json"
    profile_path = stack_dir / "paired_radial_profiles.tsv"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    if summary.get("status") != "complete_exploratory" or summary.get("science_measurement") is not False:
        raise ValueError("input stack must be complete and explicitly exploratory")
    outputs = summary["outputs"]
    if (outputs.get("paired_profiles") != str(profile_path)
            or outputs.get("paired_profiles_sha256") != sha256(profile_path)):
        raise ValueError("paired profile path or checksum differs from the verified stack summary")
    extraction = json.loads(extraction_summary_path.read_text(encoding="utf-8"))
    if extraction.get("status") != "complete" or extraction.get("science_measurement") is not False:
        raise ValueError("matched-control extraction summary is not complete exploratory output")
    control_by_target = {str(row["target_cluster_id"]): row for row in extraction["records"]}
    if len(control_by_target) != len(extraction["records"]):
        raise ValueError("matched-control extraction has duplicate target IDs")

    rows_by_filter: dict[str, list[dict]] = {"i": [], "i2": []}
    with profile_path.open(encoding="utf-8", newline="") as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            if row["filter"] not in rows_by_filter:
                continue
            target_id = str(row["target_id"])
            control = control_by_target.get(target_id)
            if control is None or str(control["control_id"]) != str(row["control_id"]):
                raise ValueError(f"target/control identity mismatch for {target_id}")
            bins = json.loads(row["profile_json"])
            outer, outer_n, outer_coverage = weighted_annulus(bins, *OUTER)
            if outer is None or outer_coverage < min_coverage:
                continue
            for annulus_index, (lo, hi) in enumerate(ANNULI):
                inner, inner_n, inner_coverage = weighted_annulus(bins, lo, hi)
                if inner is None or inner_coverage < min_coverage:
                    continue
                rows_by_filter[row["filter"]].append({
                    "target_id": target_id,
                    "control_id": str(row["control_id"]),
                    "ra_deg": float(row["ra_deg"]),
                    "dec_deg": float(row["dec_deg"]),
                    "control_ra_deg": float(control["ra_deg"]),
                    "control_dec_deg": float(control["dec_deg"]),
                    "inner_delta": inner,
                    "outer_delta": outer,
                    "inner_coverage": inner_coverage,
                    "outer_coverage": outer_coverage,
                    "inner_pixels": inner_n,
                    "outer_pixels": outer_n,
                    "annulus_index": annulus_index,
                })

    metrics = {}
    row_table = []
    for band, all_rows in rows_by_filter.items():
        metrics[band] = {}
        for index, (lo, hi) in enumerate(ANNULI):
            records = [row for row in all_rows if row["annulus_index"] == index]
            if len(records) < 8:
                raise ValueError(f"insufficient paired support for {band} {lo:g}--{hi:g} kpc")
            result = block_bootstrap(records, replicas, seed + (0 if band == "i" else 100) + index)
            result["r_min_kpc"] = lo
            result["r_max_kpc"] = hi
            metrics[band][f"{lo:g}-{hi:g}"] = result
            for row in records:
                row_table.append({
                    "filter": band,
                    "r_min_kpc": lo,
                    "r_max_kpc": hi,
                    "target_id": row["target_id"],
                    "control_id": row["control_id"],
                    "inner_delta_njy_arcsec2": row["inner_delta"],
                    "outer_400_500_delta_njy_arcsec2": row["outer_delta"],
                    "outer_subtracted_delta_njy_arcsec2": row["inner_delta"] - row["outer_delta"],
                    "inner_coverage": row["inner_coverage"],
                    "outer_coverage": row["outer_coverage"],
                    "inner_valid_pixels": row["inner_pixels"],
                    "outer_valid_pixels": row["outer_pixels"],
                })

    output_dir.mkdir(parents=True, exist_ok=False)
    rows_path = output_dir / "paired_outer_subtracted_profiles.tsv"
    fields = list(row_table[0])
    with rows_path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", extrasaction="raise")
        writer.writeheader()
        writer.writerows(row_table)
    result = {
        "schema": "astrafek-icl-outer-annulus-baseline-audit-v1",
        "status": "complete_exploratory_sensitivity",
        "science_measurement": False,
        "observational_validation_passed": False,
        "method": {
            "primary_comparison": "per-pair 400-500 kpc weighted mean subtracted from that same pair's inner annular mean",
            "radial_bin_weight": "number of valid common-support pixels in each saved 25 kpc bin",
            "pair_weight": "equal across pairs after coverage threshold",
            "minimum_common_area_fraction_per_inner_and_outer_annulus": min_coverage,
            "outer_annulus_kpc": list(OUTER),
            "inner_annuli_kpc": [list(bounds) for bounds in ANNULI],
            "uncertainty": "95% percentile interval from resampling connected 5-degree target/control sky blocks",
        },
        "by_filter": metrics,
        "bootstrap": {"replicas": replicas, "seed": seed,
                      "method": "connected 5-degree target/control block resampling"},
        "inputs": {
            "stack_summary": str(summary_path.resolve()),
            "stack_summary_sha256": sha256(summary_path),
            "paired_profiles": str(profile_path.resolve()),
            "paired_profiles_sha256": sha256(profile_path),
            "matched_control_extraction": str(extraction_summary_path.resolve()),
            "matched_control_extraction_sha256": sha256(extraction_summary_path),
        },
        "output": str(rows_path.resolve()),
        "output_sha256": sha256(rows_path),
        "limitations": [
            "The 400-500 kpc annulus may contain genuine extended ICL; subtracting it measures inner excess relative to that annulus, not total ICL flux.",
            "This is a shared-support profile sensitivity with a 5% annular coverage floor, not an alternative science estimator.",
            "Connected 5-degree block intervals do not include PSF-wing, residual-sky, selection, or full coadd-covariance systematics.",
            "The analysis uses existing stacks only; it makes no HSC request and retains science_measurement=false.",
        ],
    }
    result["script_sha256"] = sha256(Path(__file__).resolve())
    report_path = output_dir / "assessment.md"
    report_path.write_text(render_report(result), encoding="utf-8")
    result["assessment_markdown"] = str(report_path.resolve())
    result["assessment_markdown_sha256"] = sha256(report_path)
    (output_dir / "summary.json").write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stack", type=Path, default=DEFAULT_STACK)
    parser.add_argument("--extraction-summary", type=Path, default=DEFAULT_EXTRACTION)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--min-coverage", type=float, default=0.05)
    parser.add_argument("--replicas", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=20261020)
    args = parser.parse_args()
    result = analyze(args.stack.resolve(), args.extraction_summary.resolve(),
                     args.output_dir.resolve(), args.min_coverage, args.replicas, args.seed)
    print(json.dumps({"status": result["status"], "science_measurement": result["science_measurement"],
                      "by_filter": result["by_filter"], "output": str(args.output_dir.resolve())},
                     ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
