#!/usr/bin/env python3
"""Test how ICL--satellite uncertainty changes with the spatial block scale.

This is a pre-existing-sample sensitivity audit. It changes only the width of
the RA/Dec cells used to connect target/control pairs into bootstrap blocks;
it does not choose a scale from the observed association or make a science
measurement.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np

from analyze_weighted_satellite_icl_association import partial_pearson

ROOT = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
DEFAULT_INPUT = ROOT / "satellite_icl_geometric_strata_2026-10-18_01"
DEFAULT_OUTPUT = ROOT / "satellite_icl_block_scale_sensitivity_2026-10-19_01"
DEFAULT_WIDTHS = (2.5, 5.0, 10.0)
CATEGORIES = ("all", "no_geometric_overlap", "overlap_possible")
BANDS = ("i", "i2")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with Path(path).open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def output_fields(rows: list[dict]) -> list[str]:
    return list(dict.fromkeys(key for row in rows for key in row))


def _cell(row: dict[str, str], ra_key: str, dec_key: str,
          width_deg: float) -> tuple[int, int]:
    ra, dec = float(row[ra_key]), float(row[dec_key])
    if not (math.isfinite(ra) and 0 <= ra < 360):
        raise ValueError("RA must be finite and lie in [0, 360)")
    if not (math.isfinite(dec) and -90 <= dec <= 90):
        raise ValueError("Dec must be finite and lie in [-90, 90]")
    return math.floor(ra / width_deg), math.floor((dec + 90.0) / width_deg)


def connected_pair_blocks(rows: list[dict[str, str]], width_deg: float) -> np.ndarray:
    """Connect rows sharing a target or control cell, then take components."""
    width = float(width_deg)
    if not math.isfinite(width) or not 0 < width <= 180:
        raise ValueError("block width must lie in (0, 180]")
    parent = list(range(len(rows)))

    def find(index: int) -> int:
        while parent[index] != index:
            parent[index] = parent[parent[index]]
            index = parent[index]
        return index

    def union(left: int, right: int) -> None:
        root_left, root_right = find(left), find(right)
        if root_left != root_right:
            parent[max(root_left, root_right)] = min(root_left, root_right)

    owner: dict[tuple[str, tuple[int, int]], int] = {}
    for index, row in enumerate(rows):
        for role, ra_key, dec_key in (
            ("target", "target_ra_deg", "target_dec_deg"),
            ("control", "control_ra_deg", "control_dec_deg"),
        ):
            key = (role, _cell(row, ra_key, dec_key, width))
            if key in owner:
                union(index, owner[key])
            else:
                owner[key] = index
    roots = sorted({find(index) for index in range(len(rows))})
    root_labels = {root: label for label, root in enumerate(roots)}
    return np.asarray([f"pair_block_{root_labels[find(i)]}" for i in range(len(rows))])


def _correlation(rows: list[dict[str, str]], key: str) -> float:
    x = np.asarray([float(row["icl_delta_200_300"]) for row in rows])
    y = np.asarray([float(row[key]) for row in rows])
    z = np.asarray([float(row["z"]) for row in rows])
    richness = np.asarray([float(row["richness"]) for row in rows])
    return float(partial_pearson(x, y, z, richness))


def _bootstrap(rows: list[dict[str, str]], satellite_key: str, width_deg: float,
               replicas: int, rng: np.random.Generator, min_blocks: int) -> dict:
    if len(rows) < 5:
        return {"status": "insufficient_pairs", "n_pairs": len(rows)}
    labels = connected_pair_blocks(rows, width_deg)
    blocks = np.unique(labels)
    point = _correlation(rows, satellite_key)
    counts = np.bincount(np.searchsorted(blocks, labels))
    base = {
        "n_pairs": len(rows), "n_spatial_blocks": len(blocks),
        "largest_spatial_block_pairs": int(counts.max()), "partial_r": point,
    }
    if len(blocks) < min_blocks:
        return {**base, "status": "insufficient_spatial_blocks",
                "replicas_requested": replicas, "replicas_valid": 0, "ci95": None}
    indices_by_block = {block: np.flatnonzero(labels == block) for block in blocks}
    values = []
    invalid = 0
    for _ in range(replicas):
        draw = rng.choice(blocks, size=len(blocks), replace=True)
        indices = np.concatenate([indices_by_block[block] for block in draw])
        sampled = [rows[int(index)] for index in indices]
        try:
            values.append(_correlation(sampled, satellite_key))
        except (ValueError, np.linalg.LinAlgError):
            invalid += 1
    return {
        **base,
        "status": "complete_exploratory" if values else "no_valid_replicates",
        "replicas_requested": replicas, "replicas_valid": len(values),
        "replicas_invalid": invalid,
        "ci95": np.quantile(values, [0.025, 0.975]).tolist() if values else None,
    }


def _paired_common(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    no_overlap = {row["cluster_id"]: row for row in rows
                  if row.get("satellite_no_geometric_overlap", "")}
    overlap = {row["cluster_id"]: row for row in rows
               if row.get("satellite_overlap_possible", "")}
    common = []
    for cluster_id in sorted(set(no_overlap) & set(overlap)):
        left, right = no_overlap[cluster_id], overlap[cluster_id]
        if left["filter"] != right["filter"]:
            raise ValueError("paired strata have different filters")
        common.append({
            **left,
            "satellite_no_overlap": left["satellite_no_geometric_overlap"],
            "satellite_overlap_possible": right["satellite_overlap_possible"],
        })
    return common


def _paired_bootstrap(rows: list[dict[str, str]], width_deg: float, replicas: int,
                      rng: np.random.Generator, min_blocks: int) -> dict:
    if len(rows) < 5:
        return {"status": "insufficient_common_pairs", "n_pairs": len(rows)}
    labels = connected_pair_blocks(rows, width_deg)
    blocks = np.unique(labels)
    r_nonoverlap = _correlation(rows, "satellite_no_overlap")
    r_overlap = _correlation(rows, "satellite_overlap_possible")
    point = r_overlap - r_nonoverlap
    base = {"n_pairs": len(rows), "n_spatial_blocks": len(blocks),
            "r_overlap": r_overlap, "r_no_overlap": r_nonoverlap,
            "delta_r_overlap_minus_no_overlap": point}
    if len(blocks) < min_blocks:
        return {**base, "status": "insufficient_spatial_blocks",
                "replicas_requested": replicas, "replicas_valid": 0, "ci95": None}
    indices_by_block = {block: np.flatnonzero(labels == block) for block in blocks}
    values = []
    invalid = 0
    for _ in range(replicas):
        draw = rng.choice(blocks, size=len(blocks), replace=True)
        indices = np.concatenate([indices_by_block[block] for block in draw])
        sampled = [rows[int(index)] for index in indices]
        try:
            values.append(_correlation(sampled, "satellite_overlap_possible")
                          - _correlation(sampled, "satellite_no_overlap"))
        except (ValueError, np.linalg.LinAlgError):
            invalid += 1
    return {
        **base,
        "status": "complete_exploratory" if values else "no_valid_replicates",
        "replicas_requested": replicas, "replicas_valid": len(values),
        "replicas_invalid": invalid,
        "ci95": np.quantile(values, [0.025, 0.975]).tolist() if values else None,
    }


def compare_five_degree_baseline(metrics: dict, baseline: dict) -> dict:
    """Require the 5-degree run to reproduce the archived seeded bootstrap."""
    comparisons = 0

    def same_number(actual, expected, label):
        nonlocal comparisons
        if actual is None or expected is None:
            if actual is not expected:
                raise ValueError(f"5-degree baseline mismatch for {label}")
        elif not math.isclose(float(actual), float(expected), rel_tol=0, abs_tol=1e-12):
            raise ValueError(f"5-degree baseline mismatch for {label}: {actual} != {expected}")
        comparisons += 1

    def same_counts(actual, expected, label):
        nonlocal comparisons
        if int(actual) != int(expected):
            raise ValueError(f"5-degree baseline mismatch for {label}: {actual} != {expected}")
        comparisons += 1

    for band in BANDS:
        for category in CATEGORIES:
            actual = metrics[band][category]
            expected = baseline[band][category]
            same_counts(actual["n_pairs"], expected["n_cluster_filter_pairs"], f"{band}/{category}/n")
            same_counts(actual["n_spatial_blocks"], expected["n_blocks"], f"{band}/{category}/blocks")
            same_number(actual["partial_r"], expected["point_r"], f"{band}/{category}/r")
            same_counts(actual["replicas_valid"], expected["replicas_valid"], f"{band}/{category}/valid")
            same_counts(actual["replicas_invalid"], expected["replicas_invalid"], f"{band}/{category}/invalid")
            if actual["ci95"] is None or expected["ci95"] is None:
                if actual["ci95"] is not expected["ci95"]:
                    raise ValueError(f"5-degree baseline mismatch for {band}/{category}/ci95")
            else:
                for index, (left, right) in enumerate(zip(actual["ci95"], expected["ci95"])):
                    same_number(left, right, f"{band}/{category}/ci95[{index}]")
                if len(actual["ci95"]) != len(expected["ci95"]):
                    raise ValueError(f"5-degree baseline mismatch for {band}/{category}/ci95 length")

        actual = metrics[band]["paired_common_sample_difference"]
        expected = baseline[band]["paired_common_sample_difference"]
        same_counts(actual["n_pairs"], expected["n_cluster_filter_pairs"], f"{band}/paired/n")
        same_counts(actual["n_spatial_blocks"], expected["n_blocks"], f"{band}/paired/blocks")
        for field, baseline_field in (
            ("r_overlap", "r_overlap"),
            ("r_no_overlap", "r_no_overlap"),
            ("delta_r_overlap_minus_no_overlap", "point_delta_r"),
        ):
            same_number(actual[field], expected[baseline_field], f"{band}/paired/{field}")
        same_counts(actual["replicas_valid"], expected["replicas_valid"], f"{band}/paired/valid")
        same_counts(actual["replicas_invalid"], expected["replicas_invalid"], f"{band}/paired/invalid")
        if actual["ci95"] is None or expected["ci95"] is None:
            if actual["ci95"] is not expected["ci95"]:
                raise ValueError(f"5-degree baseline mismatch for {band}/paired/ci95")
        else:
            for index, (left, right) in enumerate(zip(actual["ci95"], expected["ci95"])):
                same_number(left, right, f"{band}/paired/ci95[{index}]")
            if len(actual["ci95"]) != len(expected["ci95"]):
                raise ValueError(f"5-degree baseline mismatch for {band}/paired/ci95 length")
    return {"status": "exact_match_within_1e-12", "numeric_comparisons": comparisons,
            "reference_method": "archived geometric-strata summary, same seed and 5000 replicas"}


def analyze(input_dir: Path, output: Path, widths=DEFAULT_WIDTHS,
            replicas: int = 5000, seed: int = 20261020, min_blocks: int = 6) -> dict:
    input_dir, output = Path(input_dir).resolve(), Path(output).resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite existing output: {output}")
    if replicas < 1000 or min_blocks < 2:
        raise ValueError("replicas must be >=1000 and min_blocks >=2")
    widths = tuple(float(width) for width in widths)
    if not widths or len(set(widths)) != len(widths):
        raise ValueError("block widths must be non-empty and unique")
    summary_path = input_dir / "summary.json"
    table_path = input_dir / "geometric_strata_pairs.tsv"
    source = json.loads(summary_path.read_text(encoding="utf-8"))
    if (source.get("schema") != "hsc-satellite-icl-geometric-strata-v1"
            or source.get("status") != "complete_exploratory"
            or source.get("science_measurement") is not False):
        raise ValueError("unexpected geometric-strata input summary")
    expected = source.get("outputs_sha256", {}).get(table_path.name)
    if not expected or sha256(table_path) != expected:
        raise ValueError("geometric-strata input checksum mismatch")
    raw = read_tsv(table_path)
    if len(raw) != source.get("n_association_pairs"):
        raise ValueError("geometric-strata row count mismatch")
    if len({(row["cluster_id"], row["filter"]) for row in raw}) != len(raw):
        raise ValueError("duplicate cluster/filter row")
    for width in widths:
        if not math.isfinite(width) or not 0 < width <= 180:
            raise ValueError("block widths must lie in (0, 180]")

    output_rows = []
    metrics_by_width = {}
    for width in widths:
        rng = np.random.default_rng(seed)
        metrics = {}
        for band in BANDS:
            metrics[band] = {}
            for category in CATEGORIES:
                field = {
                    "all": "satellite_unweighted",
                    "no_geometric_overlap": "satellite_no_geometric_overlap",
                    "overlap_possible": "satellite_overlap_possible",
                }[category]
                selected = []
                for row in raw:
                    if row["filter"] != band or not row.get(field, ""):
                        continue
                    selected.append({**row, "satellite_delta": row[field]})
                result = _bootstrap(selected, "satellite_delta", width, replicas, rng, min_blocks)
                metrics[band][category] = result
                output_rows.append({
                    "block_width_deg": width, "filter": band, "stratum": category,
                    **result,
                })
            common = [row for row in _paired_common(raw) if row["filter"] == band]
            paired = _paired_bootstrap(common, width, replicas, rng, min_blocks)
            metrics[band]["paired_common_sample_difference"] = paired
            output_rows.append({
                "block_width_deg": width, "filter": band,
                "stratum": "paired_common_sample_difference", **paired,
            })
        metrics_by_width[str(width)] = metrics
    baseline_check = None
    if 5.0 in widths:
        baseline_check = compare_five_degree_baseline(
            metrics_by_width["5.0"], source["metrics_by_filter"])

    output.mkdir(parents=True, exist_ok=False)
    result_path = output / "block_scale_sensitivity.tsv"
    fields = output_fields(output_rows)
    with result_path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t",
                                extrasaction="raise", lineterminator="\n")
        writer.writeheader()
        writer.writerows(output_rows)
    result = {
        "schema": "hsc-icl-satellite-block-scale-sensitivity-v1",
        "status": "complete_exploratory", "science_measurement": False,
        "method": {
            "estimand": "partial Pearson correlation after linear residualization on redshift and richness",
            "bootstrap": "cluster-filter pairs resampled by connected target/control sky cells at each declared RA/Dec cell width",
            "cell_widths_deg": list(widths), "replicas_per_group": replicas,
            "seed_restarted_for_each_width": seed, "min_spatial_blocks": min_blocks,
            "multiplicity_adjustment": "none; this is a block-scale sensitivity audit, not a new significance test",
            "five_degree_stream_order": "matches the archived geometric-strata bootstrap order, including paired stratum draws",
        },
        "input": {"summary": str(summary_path), "summary_sha256": sha256(summary_path),
                  "table": str(table_path), "table_sha256": sha256(table_path),
                  "source_rows": len(raw)},
        "metrics_by_width_deg": metrics_by_width,
        "five_degree_baseline_reproduction": baseline_check,
        "analysis_script_sha256": sha256(Path(__file__).resolve()),
        "results_tsv": str(result_path),
        "results_tsv_sha256": sha256(result_path),
        "limitations": [
            "Changing cell width changes the assumed spatial dependence scale; no width is selected from the observed effect.",
            "RA/Dec bins are not equal-area spherical cells and do not model long-range covariance beyond their connected components.",
            "The selected sample, photometric systematics, mask leakage, PSF wings, coadd noise, and measurement error are unchanged.",
            "This sensitivity analysis does not turn the exploratory association into a confirmed detection.",
        ],
    }
    summary_out = output / "summary.json"
    summary_out.write_text(json.dumps(result, ensure_ascii=False, indent=2,
                                      allow_nan=False) + "\n", encoding="utf-8")
    result["summary_sha256"] = sha256(summary_out)
    print(json.dumps({"status": result["status"], "science_measurement": False,
                      "metrics_by_width_deg": metrics_by_width,
                      "results_tsv_sha256": result["results_tsv_sha256"],
                      "summary_sha256": result["summary_sha256"]},
                     ensure_ascii=False, indent=2))
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--widths-deg", type=float, nargs="+", default=list(DEFAULT_WIDTHS))
    parser.add_argument("--replicas", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=20261020)
    parser.add_argument("--min-blocks", type=int, default=6)
    args = parser.parse_args()
    analyze(args.input, args.output, args.widths_deg, args.replicas, args.seed, args.min_blocks)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
