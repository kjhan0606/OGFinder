#!/usr/bin/env python3
"""Compare the approved target/control Gaia diagnostic samples locally.

The comparison uses the completed 492-position target match, the separately
approved 458-position control match, and their 49 matched image pairs. It does
not make another network request or approve a PSF/science measurement.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path

import numpy as np


DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
TARGET_AUDIT = DATA / "s23b_psf_gaia_audit_2026-11-07_01"
CONTROL_AUDIT = DATA / "s23b_psf_control_gaia_audit_2026-11-07_01"
CONTROL_PILOT = DATA / "s23b_psf_control_candidates_2026-11-07_01"
TARGET_PLAN = DATA / "s23b_psf_gaia_plan_2026-11-07_01"
CONTROL_PLAN = DATA / "s23b_psf_control_gaia_plan_2026-11-07_03"
DEFAULT_OUTPUT = DATA / "s23b_psf_target_control_gaia_comparison_2026-11-07_01"
BOOTSTRAP_DRAWS = 10_000
SEED = 20261107


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, payload: dict) -> None:
    with path.open("x", encoding="utf-8") as stream:
        json.dump(payload, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")


def bool_value(value: str) -> bool:
    normalized = value.strip().lower()
    if normalized in ("true", "1", "t"):
        return True
    if normalized in ("false", "0", "f"):
        return False
    raise ValueError(f"invalid boolean value: {value!r}")


def bootstrap_interval(values: list[float], rng: np.random.Generator) -> dict:
    if not values:
        return {"estimate_mean": None, "estimate_median": None, "ci95_mean": None,
                "ci95_median": None}
    data = np.asarray(values, dtype=float)
    indices = rng.integers(0, len(data), size=(BOOTSTRAP_DRAWS, len(data)))
    draws = data[indices]
    mean_draws = np.mean(draws, axis=1)
    median_draws = np.median(draws, axis=1)
    return {
        "estimate_mean": float(np.mean(data)),
        "estimate_median": float(np.median(data)),
        "ci95_mean": [float(x) for x in np.quantile(mean_draws, [0.025, 0.975])],
        "ci95_median": [float(x) for x in np.quantile(median_draws, [0.025, 0.975])],
    }


def summarize_band(rows: list[dict], rng: np.random.Generator) -> dict:
    if not rows:
        return {"n_matched_image_pairs": 0}
    target_counts = [int(row["target_morphology_candidates"]) for row in rows]
    control_counts = [int(row["control_morphology_candidates"]) for row in rows]
    target_strict = [int(row["target_gaia_strict_diagnostic_candidates"]) for row in rows]
    control_strict = [int(row["control_gaia_strict_diagnostic_candidates"]) for row in rows]
    count_diff = [c - t for t, c in zip(target_counts, control_counts)]
    strict_diff = [c - t for t, c in zip(target_strict, control_strict)]
    both_supported = [
        row for row in rows
        if int(row["target_morphology_candidates"]) >= 4
        and int(row["control_morphology_candidates"]) >= 4
    ]
    target_fwhm = [float(row["target_strict_median_fwhm_pix"]) for row in rows
                   if row["target_strict_median_fwhm_pix"] != ""]
    control_fwhm = [float(row["control_strict_median_fwhm_pix"]) for row in rows
                    if row["control_strict_median_fwhm_pix"] != ""]
    return {
        "n_matched_image_pairs": len(rows),
        "n_pairs_with_at_least_four_morphology_candidates_on_both_sides": len(both_supported),
        "morphology_candidates": {
            "target_total": sum(target_counts),
            "control_total": sum(control_counts),
            "target_median_per_pair": float(np.median(target_counts)),
            "control_median_per_pair": float(np.median(control_counts)),
            "control_minus_target_paired_count_difference": bootstrap_interval(count_diff, rng),
            "frames_with_at_least_four": {
                "target": sum(value >= 4 for value in target_counts),
                "control": sum(value >= 4 for value in control_counts),
            },
        },
        "gaia_diagnostic_strict_candidates": {
            "target_total": sum(target_strict),
            "control_total": sum(control_strict),
            "target_median_per_pair": float(np.median(target_strict)),
            "control_median_per_pair": float(np.median(control_strict)),
            "control_minus_target_paired_count_difference": bootstrap_interval(strict_diff, rng),
            "pairs_by_nonzero_strict_candidate_support": {
                "both": sum(t > 0 and c > 0 for t, c in zip(target_strict, control_strict)),
                "target_only": sum(t > 0 and c == 0 for t, c in zip(target_strict, control_strict)),
                "control_only": sum(t == 0 and c > 0 for t, c in zip(target_strict, control_strict)),
                "neither": sum(t == 0 and c == 0 for t, c in zip(target_strict, control_strict)),
            },
        },
        "strict_candidate_fwhm_median_pix": {
            "target": float(np.median(target_fwhm)) if target_fwhm else None,
            "control": float(np.median(control_fwhm)) if control_fwhm else None,
            "n_target": len(target_fwhm),
            "n_control": len(control_fwhm),
        },
        "strict_candidate_yield_among_morphology_candidates": {
            "target": sum(target_strict) / sum(target_counts) if sum(target_counts) else None,
            "control": sum(control_strict) / sum(control_counts) if sum(control_counts) else None,
        },
        "candidate_count_difference_definition": "control minus target, paired by matched target/control image",
        "bootstrap": {
            "unit": "matched image pair",
            "draws": BOOTSTRAP_DRAWS,
            "seed": SEED,
            "interval": "percentile 95% interval; descriptive with this small pilot",
        },
    }


def analyze(output: Path) -> dict:
    if output.exists():
        raise FileExistsError(f"refusing to overwrite comparison output: {output}")
    target_summary_path = TARGET_AUDIT / "summary.json"
    target_evidence_path = TARGET_AUDIT / "candidate_evidence.tsv"
    control_summary_path = CONTROL_AUDIT / "summary.json"
    control_evidence_path = CONTROL_AUDIT / "control_candidate_evidence.tsv"
    selection_path = CONTROL_PILOT / "selection.tsv"
    target_requests_path = TARGET_PLAN / "requests.json"
    control_requests_path = CONTROL_PLAN / "requests.json"
    required = [target_summary_path, target_evidence_path, control_summary_path,
                control_evidence_path, selection_path, target_requests_path,
                control_requests_path, TARGET_PLAN / "plan.json",
                CONTROL_PLAN / "plan.json",
                CONTROL_PILOT / "summary.json"]
    for path in required:
        if not path.is_file():
            raise FileNotFoundError(path)

    target_summary = read_json(target_summary_path)
    control_summary = read_json(control_summary_path)
    if (target_summary.get("status") != "complete"
            or target_summary.get("n_candidate_positions") != 492
            or target_summary.get("science_measurement") is not False
            or target_summary.get("psf_approved") is not False):
        raise ValueError("approved target Gaia audit is incomplete or has unexpected state")
    if (control_summary.get("status") != "complete_local_diagnostic_only"
            or control_summary.get("n_candidate_positions") != 458
            or control_summary.get("science_measurement") is not False
            or control_summary.get("psf_approved") is not False):
        raise ValueError("separately approved control Gaia audit is incomplete or has unexpected state")

    for audit_summary in (target_summary, control_summary):
        for input_name, expected_hash in audit_summary.get("inputs_sha256", {}).items():
            input_path = Path(input_name)
            if not input_path.is_file() or sha256(input_path) != expected_hash:
                raise ValueError(f"Gaia audit input hash mismatch: {input_path}")
    if (target_summary.get("output_sha256", {}).get(target_evidence_path.name) != sha256(target_evidence_path)
            or control_summary.get("output_sha256", {}).get(control_evidence_path.name) != sha256(control_evidence_path)):
        raise ValueError("Gaia candidate-evidence output hash differs from its audit summary")

    pairs = read_tsv(selection_path)
    target_rows = read_tsv(target_evidence_path)
    control_rows = read_tsv(control_evidence_path)
    target_requests = read_json(target_requests_path)["candidates"]
    control_requests = read_json(control_requests_path)["candidates"]
    if len(pairs) != 49 or len(target_rows) != 492 or len(control_rows) != 458:
        raise ValueError("unexpected target/control pair or candidate table size")
    if len(target_requests) != len(target_rows) or len(control_requests) != len(control_rows):
        raise ValueError("request and Gaia audit row counts differ")

    target_by_cluster: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in target_rows:
        target_by_cluster[row["cluster_id"]].append(row)
    control_by_id: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in control_rows:
        control_by_id[row["control_id"]].append(row)
    target_requests_by_index = {int(row["local_index"]): row for row in target_requests}
    control_requests_by_index = {int(row["local_index"]): row for row in control_requests}
    if len(target_requests_by_index) != len(target_requests) or len(control_requests_by_index) != len(control_requests):
        raise ValueError("request candidate indices are duplicated")

    target_fwhm_by_key = {
        (row["cluster_id"], int(row["number"])): float(row["fwhm_moment_pix"])
        for row in target_requests
    }
    control_fwhm_by_key = {
        (row["cluster_id"], int(row["number"])): float(row["fwhm_moment_pix"])
        for row in control_requests
    }
    if len(target_fwhm_by_key) != len(target_requests) or len(control_fwhm_by_key) != len(control_requests):
        raise ValueError("duplicate candidate identity in a request list")
    control_pilot_summary = read_json(CONTROL_PILOT / "summary.json")
    if (sha256(selection_path) != control_pilot_summary.get("selection_sha256")
            or sha256(TARGET_PLAN / "plan.json") not in target_summary.get("inputs_sha256", {}).values()
            or sha256(target_requests_path) not in target_summary.get("inputs_sha256", {}).values()
            or sha256(CONTROL_PLAN / "plan.json") not in control_summary.get("inputs_sha256", {}).values()
            or sha256(control_requests_path) not in control_summary.get("inputs_sha256", {}).values()):
        raise ValueError("paired selection or Gaia plan/request is not linked to its audit receipts")

    pair_rows = []
    seen_targets, seen_controls = set(), set()
    for pair in pairs:
        target_id = pair["target_cluster_id"]
        control_id = pair["control_id"]
        band = pair["target_primary_filter"]
        if (target_id in seen_targets or control_id in seen_controls
                or pair["control_primary_filter"] != band):
            raise ValueError("duplicate or filter-mismatched image pair")
        seen_targets.add(target_id)
        seen_controls.add(control_id)
        target_evidence = target_by_cluster.get(target_id, [])
        control_evidence = control_by_id.get(control_id, [])
        target_strict_rows = [row for row in target_evidence if bool_value(row["strict_candidate"])]
        control_strict_rows = [row for row in control_evidence if bool_value(row["strict_candidate"])]
        target_count = len(target_evidence)
        control_count = len(control_evidence)
        target_strict = len(target_strict_rows)
        control_strict = len(control_strict_rows)
        target_fwhm = [target_fwhm_by_key[(target_id, int(row["number"]))]
                       for row in target_strict_rows]
        control_fwhm = [control_fwhm_by_key[(control_id, int(row["number"]))]
                        for row in control_strict_rows]
        pair_rows.append({
            "target_cluster_id": target_id,
            "control_id": control_id,
            "filter": band,
            "target_morphology_candidates": target_count,
            "control_morphology_candidates": control_count,
            "target_gaia_match_positions": sum(int(row["n_gaia_matches_within_2arcsec"]) > 0 for row in target_evidence),
            "control_gaia_match_positions": sum(int(row["n_gaia_matches_within_2arcsec"]) > 0 for row in control_evidence),
            "target_gaia_strict_diagnostic_candidates": target_strict,
            "control_gaia_strict_diagnostic_candidates": control_strict,
            "control_minus_target_strict_count": control_strict - target_strict,
            "target_strict_median_fwhm_pix": f"{float(np.median(target_fwhm)):.8g}" if target_fwhm else "",
            "control_strict_median_fwhm_pix": f"{float(np.median(control_fwhm)):.8g}" if control_fwhm else "",
            "target_seeing_arcsec": pair["target_seeing_arcsec"],
            "control_seeing_arcsec": pair["control_seeing_arcsec"],
            "control_minus_target_seeing_arcsec": (
                float(pair["control_seeing_arcsec"]) - float(pair["target_seeing_arcsec"])
            ),
        })

    if len(seen_targets) != 49 or len(seen_controls) != 49:
        raise ValueError("matched image selection does not contain 49 unique pairs")
    output.mkdir(parents=True, exist_ok=False)
    pair_path = output / "matched_pair_gaia_candidate_support.tsv"
    write_tsv(pair_path, pair_rows, list(pair_rows[0]))

    rng = np.random.default_rng(SEED)
    by_filter = {
        band: summarize_band([row for row in pair_rows if row["filter"] == band], rng)
        for band in ("i", "i2")
    }
    input_paths = [*required]
    summary = {
        "schema": "astrafek-s23b-matched-psf-gaia-candidate-comparison-v1",
        "status": "complete_local_paired_diagnostic_only",
        "target_candidates": len(target_rows),
        "control_candidates": len(control_rows),
        "n_matched_image_pairs": len(pair_rows),
        "n_gaia_target_match_rows": target_summary.get("n_gaia_match_rows"),
        "n_gaia_control_match_rows": control_summary.get("n_gaia_match_rows"),
        "target_evidence_levels": target_summary.get("evidence_levels"),
        "control_evidence_levels": control_summary.get("evidence_levels"),
        "by_filter": by_filter,
        "interpretation_limits": [
            "Strict candidates are threshold-based Gaia astrometric diagnostics, not a calibrated stellar or PSF classifier.",
            "Gaia non-matches do not establish that HSC candidates are non-stellar.",
            "The effective HSC coadd epoch is not known; no Gaia position propagation was applied.",
            "Target and random-sky control fields differ in cluster light and source crowding; morphology counts are not stellar surface densities.",
            "This paired comparison does not measure PSF wings, approve a PSF, or correct the ICL stack.",
        ],
        "inputs_sha256": {str(path.resolve()): sha256(path) for path in input_paths},
        "output_sha256": {pair_path.name: sha256(pair_path)},
        "science_measurement": False,
        "psf_approved": False,
        "images_downloaded": False,
        "additional_external_queries_submitted": False,
    }
    summary_path = output / "summary.json"
    write_json(summary_path, summary)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    result = analyze(args.output.resolve())
    print(json.dumps({
        "status": result["status"],
        "n_matched_image_pairs": result["n_matched_image_pairs"],
        "target_candidates": result["target_candidates"],
        "control_candidates": result["control_candidates"],
        "target_evidence_levels": result["target_evidence_levels"],
        "control_evidence_levels": result["control_evidence_levels"],
        "by_filter": result["by_filter"],
        "science_measurement": result["science_measurement"],
        "psf_approved": result["psf_approved"],
        "additional_external_queries_submitted": result["additional_external_queries_submitted"],
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
