#!/usr/bin/env python3
"""Measure pair-level coherence between ICL annuli and the 400--500 kpc baseline."""
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
DEFAULT_OUTER = ROOT / "icl_outer_annulus_baseline_2026-10-20_02"
DEFAULT_EXTRACTION = ROOT / "prestack_extraction_2026-10-15_01/extraction_summary.json"
DEFAULT_MANIFEST = ROOT / "sample_manifest.json"
ANNULI = ((50.0, 100.0), (100.0, 200.0), (200.0, 300.0))
FILTERS = ("i", "i2")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def _correlation(x: np.ndarray, y: np.ndarray) -> float:
    if len(x) < 3 or np.std(x) == 0 or np.std(y) == 0:
        return float("nan")
    return float(np.corrcoef(x, y)[0, 1])


def regression_metrics(inner: list[float] | np.ndarray,
                       outer: list[float] | np.ndarray,
                       covariates: list[list[float]] | np.ndarray) -> dict[str, float]:
    """Return raw and covariate-adjusted slope/correlation for inner on outer."""
    y = np.asarray(inner, dtype=float)
    x = np.asarray(outer, dtype=float)
    cov = np.asarray(covariates, dtype=float)
    if y.ndim != 1 or x.shape != y.shape or cov.ndim != 2 or cov.shape[0] != len(y):
        raise ValueError("inner, outer, and covariate arrays have incompatible shapes")
    if not np.isfinite(y).all() or not np.isfinite(x).all() or not np.isfinite(cov).all():
        raise ValueError("regression inputs must be finite")
    design = np.column_stack((np.ones(len(y)), cov))
    if len(y) <= design.shape[1] + 2 or np.linalg.matrix_rank(design) != design.shape[1]:
        raise ValueError("insufficient rows or rank-deficient covariate design")
    beta_y, *_ = np.linalg.lstsq(design, y, rcond=None)
    beta_x, *_ = np.linalg.lstsq(design, x, rcond=None)
    y_resid = y - design @ beta_y
    x_resid = x - design @ beta_x
    x_var = float(np.dot(x_resid, x_resid))
    if x_var <= 0:
        raise ValueError("outer baseline has no residual variation")
    return {
        "raw_pearson_r": _correlation(x, y),
        "raw_slope": float(np.dot(x - x.mean(), y - y.mean()) /
                            np.dot(x - x.mean(), x - x.mean())),
        "partial_pearson_r": _correlation(x_resid, y_resid),
        "partial_slope": float(np.dot(x_resid, y_resid) / x_var),
    }


def _extract_stats(rows: list[dict], indices: np.ndarray) -> dict[str, float]:
    selected = [rows[int(i)] for i in indices]
    return regression_metrics(
        [r["inner"] for r in selected], [r["outer"] for r in selected],
        [[r["z"], r["richness"], r["inner_coverage"], r["outer_coverage"]]
         for r in selected],
    )


def bootstrap_coherence(rows: list[dict], replicas: int, seed: int) -> dict:
    if replicas < 100:
        raise ValueError("at least 100 connected-block replicates are required")
    target = [{"ra_deg": r["ra_deg"], "dec_deg": r["dec_deg"]} for r in rows]
    control = [{"ra_deg": r["control_ra_deg"], "dec_deg": r["control_dec_deg"]}
               for r in rows]
    _, labels, n_target, n_control, n_blocks, largest = pair_block_labels(
        target, control, np.ones(len(rows), dtype=bool))
    unique = np.unique(labels)
    if n_blocks < 3:
        raise ValueError("at least three connected target/control blocks are required")
    members = [np.flatnonzero(labels == label) for label in unique]
    full = _extract_stats(rows, np.arange(len(rows)))
    rng = np.random.default_rng(seed)
    metrics = ("raw_pearson_r", "raw_slope", "partial_pearson_r", "partial_slope")
    draws = {name: [] for name in metrics}
    for _ in range(replicas):
        sampled = rng.integers(0, len(members), size=len(members))
        index = np.concatenate([members[int(i)] for i in sampled])
        try:
            result = _extract_stats(rows, index)
        except ValueError:
            continue
        for name in metrics:
            draws[name].append(result[name])
    intervals = {}
    for name, values in draws.items():
        array = np.asarray(values, dtype=float)
        array = array[np.isfinite(array)]
        if len(array) < 0.8 * replicas:
            raise ValueError(f"too few valid block bootstrap draws for {name}: {len(array)}")
        intervals[name] = {
            "q025": float(np.quantile(array, 0.025)),
            "q975": float(np.quantile(array, 0.975)),
            "valid_replicates": int(len(array)),
        }
    return {
        "n_pairs": len(rows), "n_target_cells": int(n_target),
        "n_control_cells": int(n_control), "n_connected_blocks": int(n_blocks),
        "largest_block_pairs": int(largest), **full,
        "block_bootstrap_95_intervals": intervals,
        "bootstrap_replicates_requested": replicas, "bootstrap_seed": seed,
    }


def analyze(stack_dir: Path, outer_dir: Path, extraction_path: Path,
            manifest_path: Path, output_dir: Path, replicas: int = 5000,
            seed: int = 20261020) -> dict:
    if output_dir.exists():
        raise FileExistsError(f"refusing to overwrite {output_dir}")
    stack_summary_path = stack_dir / "summary.json"
    stack_path = stack_dir / "paired_radial_profiles.tsv"
    outer_summary_path = outer_dir / "summary.json"
    outer_path = outer_dir / "paired_outer_subtracted_profiles.tsv"
    stack_summary = json.loads(stack_summary_path.read_text(encoding="utf-8"))
    outer_summary = json.loads(outer_summary_path.read_text(encoding="utf-8"))
    if stack_summary.get("status") != "complete_exploratory" or stack_summary.get("science_measurement") is not False:
        raise ValueError("source stack is not complete exploratory output")
    if outer_summary.get("status") != "complete_exploratory_sensitivity" or outer_summary.get("science_measurement") is not False:
        raise ValueError("outer-baseline output is not complete exploratory output")
    if outer_summary.get("inputs", {}).get("matched_control_extraction_sha256") != sha256(extraction_path):
        raise ValueError("matched-control extraction checksum differs from the outer-baseline record")
    if stack_summary.get("outputs", {}).get("paired_profiles_sha256") != sha256(stack_path):
        raise ValueError("source paired-profile checksum mismatch")
    if outer_summary.get("output_sha256") != sha256(outer_path):
        raise ValueError("outer-baseline row checksum mismatch")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("science_measurement") is not False or manifest.get("schema") != "hsc-bcg-diffuse-sample-v1":
        raise ValueError("sample manifest is unexpected or science-promoted")
    cluster_by_id = {str(row["id"]): row for row in manifest["clusters"]}
    if len(cluster_by_id) != len(manifest["clusters"]):
        raise ValueError("duplicate target IDs in sample manifest")
    extraction = json.loads(extraction_path.read_text(encoding="utf-8"))
    if extraction.get("status") != "complete" or extraction.get("science_measurement") is not False:
        raise ValueError("matched-control extraction is incomplete or science-promoted")
    controls = {str(row["target_cluster_id"]): row for row in extraction["records"]}
    if len(controls) == 0:
        raise ValueError("matched-control extraction has no records")
    if len(controls) != len(extraction["records"]):
        raise ValueError("duplicate matched-control target IDs")
    stack_rows = {}
    for row in read_tsv(stack_path):
        key = (str(row["target_id"]), row["filter"])
        if key in stack_rows:
            raise ValueError(f"duplicate source profile row {key}")
        stack_rows[key] = row

    records: dict[tuple[str, str], list[dict]] = {}
    baseline_by_key: dict[tuple[str, str], float] = {}
    for row in read_tsv(outer_path):
        band = row["filter"]
        if band not in FILTERS:
            raise ValueError(f"unexpected filter {band}")
        target_id, control_id = str(row["target_id"]), str(row["control_id"])
        key = (target_id, band)
        source = stack_rows.get(key)
        cluster = cluster_by_id.get(target_id)
        control_meta = controls.get(target_id)
        if source is None or cluster is None or control_meta is None:
            raise ValueError(f"missing linked target/control metadata for {key}")
        if str(source["control_id"]) != control_id or str(control_meta["control_id"]) != control_id:
            raise ValueError(f"target/control mapping mismatch for {key}")
        if abs(float(source["z"]) - float(cluster["z"])) > 1e-7:
            raise ValueError(f"target redshift mismatch for {key}")
        if (abs(float(source["ra_deg"]) - float(cluster["ra_deg"])) > 1e-7
                or abs(float(source["dec_deg"]) - float(cluster["dec_deg"])) > 1e-7):
            raise ValueError(f"target coordinates mismatch for {key}")
        lo, hi = float(row["r_min_kpc"]), float(row["r_max_kpc"])
        inner_coverage = float(row["inner_coverage"])
        outer_coverage = float(row["outer_coverage"])
        if inner_coverage < 0.05 or outer_coverage < 0.05:
            raise ValueError(f"outer-baseline table violates its 5% coverage rule for {key}")
        outer_value = float(row["outer_400_500_delta_njy_arcsec2"])
        old_outer = baseline_by_key.setdefault(key, outer_value)
        if abs(old_outer - outer_value) > 1e-10:
            raise ValueError(f"outer baseline varies across inner annuli for {key}")
        records.setdefault((band, f"{lo:g}-{hi:g}"), []).append({
            "target_id": target_id, "control_id": control_id,
            "ra_deg": float(source["ra_deg"]), "dec_deg": float(source["dec_deg"]),
            "control_ra_deg": float(control_meta["ra_deg"]),
            "control_dec_deg": float(control_meta["dec_deg"]),
            "z": float(cluster["z"]), "richness": float(cluster["richness"]),
            "inner": float(row["inner_delta_njy_arcsec2"]),
            "outer": outer_value,
            "inner_coverage": inner_coverage,
            "outer_coverage": outer_coverage,
        })

    metrics = {}
    for band_index, band in enumerate(FILTERS):
        metrics[band] = {}
        for annulus_index, (lo, hi) in enumerate(ANNULI):
            name = f"{lo:g}-{hi:g}"
            rows = records.get((band, name), [])
            if len(rows) < 12:
                raise ValueError(f"insufficient rows for {band}/{name}: {len(rows)}")
            result = bootstrap_coherence(rows, replicas,
                                         seed + band_index * 100 + annulus_index)
            result["r_min_kpc"], result["r_max_kpc"] = lo, hi
            metrics[band][name] = result

    result = {
        "schema": "astrafek-icl-outer-baseline-coherence-v1",
        "status": "complete_exploratory_coherence",
        "science_measurement": False,
        "method": {
            "response": "pair-level raw target-minus-control inner-annulus mean",
            "predictor": "same pair's 400-500 kpc target-minus-control mean",
            "adjusted_covariates": ["cluster redshift", "CAMIRA richness", "inner common-area fraction", "outer common-area fraction"],
            "uncertainty": "95% percentile interval from resampling connected 5-degree target/control spatial blocks",
            "weights": "equal pair weights after the inherited 5% support selection",
        },
        "by_filter": metrics,
        "bootstrap": {"replicates": replicas, "seed": seed},
        "inputs": {
            "stack_summary": str(stack_summary_path.resolve()), "stack_summary_sha256": sha256(stack_summary_path),
            "paired_profiles": str(stack_path.resolve()), "paired_profiles_sha256": sha256(stack_path),
            "outer_summary": str(outer_summary_path.resolve()), "outer_summary_sha256": sha256(outer_summary_path),
            "outer_rows": str(outer_path.resolve()), "outer_rows_sha256": sha256(outer_path),
            "matched_control_extraction": str(extraction_path.resolve()), "matched_control_extraction_sha256": sha256(extraction_path),
            "sample_manifest": str(manifest_path.resolve()), "sample_manifest_sha256": sha256(manifest_path),
        },
        "limitations": [
            "outer-annulus covariance with inner light cannot distinguish residual sky from real extended ICL or cluster-mass variation",
            "the inherited 5% support cut changes the selected pairs across annuli and may select atypical portions of each image",
            "block-bootstrap intervals reflect only the chosen connected 5-degree target/control blocks and omit full pixel/coadd covariance and calibration systematics",
            "six filter-by-annulus comparisons are descriptive; no confirmatory p-values or multiplicity claims are made",
            "no HSC query or download was made; science_measurement=false",
        ],
    }
    result["script_sha256"] = sha256(Path(__file__).resolve())
    output_dir.mkdir(parents=True, exist_ok=False)
    (output_dir / "summary.json").write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    rows_out = []
    for band in FILTERS:
        for annulus, item in metrics[band].items():
            rows_out.append({"filter": band, "annulus_kpc": annulus, **item})
    with (output_dir / "coherence_metrics.jsonl").open("x", encoding="utf-8") as stream:
        for row in rows_out:
            stream.write(json.dumps(row, ensure_ascii=False) + "\n")
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--stack-dir", type=Path, default=DEFAULT_STACK)
    parser.add_argument("--outer-dir", type=Path, default=DEFAULT_OUTER)
    parser.add_argument("--extraction-summary", type=Path, default=DEFAULT_EXTRACTION)
    parser.add_argument("--sample-manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--replicas", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=20261020)
    args = parser.parse_args()
    result = analyze(args.stack_dir.resolve(), args.outer_dir.resolve(),
                     args.extraction_summary.resolve(), args.sample_manifest.resolve(),
                     args.output_dir.resolve(), args.replicas, args.seed)
    print(json.dumps({"status": result["status"], "science_measurement": result["science_measurement"],
                      "by_filter": result["by_filter"], "output": str(args.output_dir.resolve())},
                     ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
