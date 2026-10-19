#!/usr/bin/env python3
"""Block-width sensitivity for the full 94-row per-annulus ICL association."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from audit_satellite_icl_block_scale_sensitivity import connected_pair_blocks
from audit_satellite_icl_spatial_blocks import (
    DEFAULT_ASSOC, DEFAULT_EXTRACTION, DEFAULT_ICL, load_verified_rows,
)

ROOT = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
DEFAULT_BASELINE = ROOT / "satellite_icl_spatial_blocks_2026-10-20_02/summary.json"
DEFAULT_OUTPUT = ROOT / "satellite_icl_perbin_block_scale_2026-10-19_01"
WIDTHS = (2.5, 5.0, 10.0)
BANDS = ("i2", "i")
SELECTIONS = ("all_supported_pairs", "center_eligible")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def partial_pearson(rows: list[dict]) -> float:
    x = np.asarray([float(row["icl_delta_200_300"]) for row in rows])
    y = np.asarray([float(row["satellite_delta"]) for row in rows])
    z = np.asarray([float(row["redshift"]) for row in rows])
    richness = np.asarray([float(row["richness"]) for row in rows])
    design = np.column_stack((np.ones(len(rows)), z, richness))
    if len(rows) < 5 or np.linalg.matrix_rank(design) < design.shape[1]:
        raise ValueError("partial correlation requires >=5 rows and full-rank covariates")
    rx = x - design @ np.linalg.lstsq(design, x, rcond=None)[0]
    ry = y - design @ np.linalg.lstsq(design, y, rcond=None)[0]
    denom = float(np.linalg.norm(rx) * np.linalg.norm(ry))
    if not math.isfinite(denom) or denom <= 0:
        raise ValueError("partial correlation is undefined")
    return float(rx @ ry / denom)


def bootstrap(rows: list[dict], width_deg: float, replicas: int, seed: int,
              min_blocks: int = 6) -> dict:
    labels = connected_pair_blocks(rows, width_deg)
    blocks = np.unique(labels)
    point = partial_pearson(rows)
    counts = np.bincount(np.searchsorted(blocks, labels))
    base = {
        "n_pairs": len(rows), "n_spatial_blocks": len(blocks),
        "largest_block_pairs": int(counts.max()), "partial_r": point,
    }
    if len(blocks) < min_blocks:
        return {**base, "status": "insufficient_spatial_blocks",
                "replicas_requested": replicas, "replicas_valid": 0,
                "replicas_invalid": 0, "ci95": None}
    groups = {block: np.flatnonzero(labels == block) for block in blocks}
    rng = np.random.default_rng(seed)
    values = []
    invalid = 0
    for _ in range(replicas):
        draw = rng.choice(blocks, size=len(blocks), replace=True)
        indices = np.concatenate([groups[block] for block in draw])
        sample = [rows[int(index)] for index in indices]
        try:
            values.append(partial_pearson(sample))
        except (ValueError, np.linalg.LinAlgError):
            invalid += 1
    return {
        **base, "status": "complete_exploratory" if values else "no_valid_replicates",
        "replicas_requested": replicas, "replicas_valid": len(values),
        "replicas_invalid": invalid,
        "ci95": np.percentile(values, [2.5, 97.5]).tolist() if values else None,
    }


def compare_baseline(metrics: dict, baseline: dict) -> dict:
    comparisons = 0
    reference = {(row["filter"], row["selection"]): row
                 for row in baseline.get("results", [])}
    for band in BANDS:
        for selection in SELECTIONS:
            actual = metrics[band][selection]
            row = reference.get((band, selection))
            if row is None:
                raise ValueError(f"missing archived 5-degree result for {band}/{selection}")
            expected = row["partial_spatial_block_bootstrap"]
            count_pairs = (actual["n_pairs"], expected["n_pairs"])
            count_blocks = (actual["n_spatial_blocks"], expected["n_blocks"])
            if count_pairs[0] != count_pairs[1] or count_blocks[0] != count_blocks[1]:
                raise ValueError(f"5-degree support mismatch for {band}/{selection}")
            comparisons += 2
            for left, right, label in (
                (actual["partial_r"], expected["partial_r_point"], "point_r"),
                (actual["replicas_valid"], expected["replicas_valid"], "valid_replicas"),
                (actual["replicas_invalid"], expected["replicas_invalid"], "invalid_replicas"),
            ):
                if label in ("valid_replicas", "invalid_replicas"):
                    same = int(left) == int(right)
                else:
                    same = math.isclose(float(left), float(right), rel_tol=0, abs_tol=1e-12)
                if not same:
                    raise ValueError(f"5-degree mismatch for {band}/{selection}/{label}")
                comparisons += 1
            actual_ci, expected_ci = actual["ci95"], expected["ci95_partial_r"]
            if len(actual_ci) != len(expected_ci) or any(
                not math.isclose(float(a), float(b), rel_tol=0, abs_tol=1e-12)
                for a, b in zip(actual_ci, expected_ci)
            ):
                raise ValueError(f"5-degree interval mismatch for {band}/{selection}")
            comparisons += len(actual_ci)
    return {"status": "exact_match_within_1e-12", "numeric_comparisons": comparisons}


def analyze(association: Path, icl: Path, extraction: Path, baseline_path: Path,
            output: Path, widths=WIDTHS, replicas: int = 5000,
            seed: int = 20261020, min_blocks: int = 6) -> dict:
    output = Path(output).resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite existing output: {output}")
    if replicas < 1000 or min_blocks < 2:
        raise ValueError("replicas must be >=1000 and minimum blocks >=2")
    widths = tuple(float(width) for width in widths)
    if not widths or len(set(widths)) != len(widths) or any(
        not math.isfinite(width) or not 0 < width <= 180 for width in widths
    ):
        raise ValueError("block widths must be unique and lie in (0, 180]")
    rows, provenance = load_verified_rows(association, icl, extraction)
    if len(rows) != 94:
        raise ValueError(f"expected the fixed 94-row per-annulus sample, found {len(rows)}")
    baseline_path = Path(baseline_path).resolve()
    baseline = json.loads(baseline_path.read_text(encoding="utf-8"))
    if (baseline.get("schema") != "hsc-satellite-icl-spatial-block-bootstrap-v1"
            or baseline.get("status") != "complete_exploratory"
            or baseline.get("science_measurement") is not False):
        raise ValueError("unexpected archived spatial-block baseline")
    for key, value in provenance.items():
        if baseline.get("provenance", {}).get(key) != value:
            raise ValueError(f"baseline provenance mismatch for {key}")

    metrics_by_width = {}
    flat_rows = []
    for width in widths:
        metrics = {}
        for band_index, band in enumerate(BANDS):
            metrics[band] = {}
            band_rows = [row for row in rows if row["filter"] == band]
            for selection_index, selection in enumerate(SELECTIONS):
                selected = band_rows if selection == "all_supported_pairs" else [
                    row for row in band_rows if str(row["center_eligible"]).lower() == "true"
                ]
                normalized = [{
                    **row,
                    "z": row["redshift"],
                    "satellite_delta": row["satellite_delta_50_75"],
                } for row in selected]
                # Match the archived per-group seeds exactly for the 5-degree gate.
                group_seed = seed + 2 * band_index + selection_index
                result = bootstrap(normalized, width, replicas, group_seed, min_blocks)
                metrics[band][selection] = result
                flat_rows.append({"block_width_deg": width, "filter": band,
                                  "selection": selection, **result, "seed": group_seed})
        metrics_by_width[str(width)] = metrics
    baseline_check = None
    if 5.0 in widths:
        baseline_check = compare_baseline(metrics_by_width["5.0"], baseline)

    output.mkdir(parents=True, exist_ok=False)
    result_path = output / "perbin_block_scale_sensitivity.tsv"
    fields = list(dict.fromkeys(key for row in flat_rows for key in row))
    with result_path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t",
                                extrasaction="raise", lineterminator="\n")
        writer.writeheader()
        writer.writerows(flat_rows)
    result = {
        "schema": "hsc-satellite-icl-perbin-block-scale-sensitivity-v1",
        "status": "complete_exploratory", "science_measurement": False,
        "method": {
            "fixed_input_sample": "all 94 per-annulus matched cluster/filter rows; no reselection",
            "primary_filter": "i2", "cell_widths_deg": list(widths),
            "bootstrap": "connected target/control RA/Dec cells; resample spatial components with replacement",
            "replicas_per_group": replicas, "min_spatial_blocks": min_blocks,
            "seed_rule": "20261020 + 2*band_index + selection_index; same seed at each cell width",
            "multiplicity_adjustment": "none; descriptive sensitivity only",
        },
        "input": {"association": str(Path(association).resolve()),
                  "icl": str(Path(icl).resolve()),
                  "extraction_summary": str(Path(extraction).resolve()),
                  "verified_provenance": provenance,
                  "baseline_summary": str(baseline_path),
                  "baseline_summary_sha256": sha256(baseline_path)},
        "metrics_by_width_deg": metrics_by_width,
        "five_degree_baseline_reproduction": baseline_check,
        "analysis_script_sha256": sha256(Path(__file__).resolve()),
        "results_tsv": str(result_path), "results_tsv_sha256": sha256(result_path),
        "limitations": [
            "Block widths are post hoc sensitivity choices; no width is selected from the observed association.",
            "RA/Dec rectangular cells are not equal-area spherical regions and do not model all long-range covariance.",
            "The selected support sample and its mask, PSF, sky, extinction, membership, and coadd-noise uncertainties remain.",
            "Intervals are not confirmatory significance tests and no multiplicity correction is applied.",
        ],
    }
    summary_out = output / "summary.json"
    summary_out.write_text(json.dumps(result, ensure_ascii=False, indent=2,
                                      allow_nan=False) + "\n", encoding="utf-8")
    result["summary_sha256"] = sha256(summary_out)
    print(json.dumps({"status": result["status"], "science_measurement": False,
                      "five_degree_baseline_reproduction": baseline_check,
                      "i2": {width: metrics_by_width[width]["i2"]
                             for width in metrics_by_width},
                      "results_tsv_sha256": result["results_tsv_sha256"],
                      "summary_sha256": result["summary_sha256"]},
                     ensure_ascii=False, indent=2))
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--association", type=Path, default=DEFAULT_ASSOC)
    parser.add_argument("--icl", type=Path, default=DEFAULT_ICL)
    parser.add_argument("--extraction-summary", type=Path, default=DEFAULT_EXTRACTION)
    parser.add_argument("--baseline", type=Path, default=DEFAULT_BASELINE)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--widths-deg", type=float, nargs="+", default=list(WIDTHS))
    parser.add_argument("--replicas", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=20261020)
    parser.add_argument("--min-blocks", type=int, default=6)
    args = parser.parse_args()
    analyze(args.association, args.icl, args.extraction_summary, args.baseline,
            args.output, args.widths_deg, args.replicas, args.seed, args.min_blocks)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
