#!/usr/bin/env python3
"""Audit the approved S23B HSC S19A count-only result locally.

The audit verifies all count batches and reports descriptive candidate-count
distribution and sample-selection diagnostics. It never contacts CAS and never
reads or exports source-row, shape, coordinate, or photo-z values.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import re

import numpy as np
from scipy.stats import spearmanr

from submit_s23b_wl_shape_counts import read_count_csv, validate_plan


DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
DEFAULT_PLAN = DATA / "wl_s23b_shape_count_plan_2026-11-07_01"
DEFAULT_RESULTS = DEFAULT_PLAN / "results"
DEFAULT_OUTPUT = DATA / "wl_s23b_shape_count_audit_2026-11-07_01"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def summarize(values: np.ndarray) -> dict:
    return {
        "minimum": int(values.min()),
        "p10": float(np.quantile(values, 0.10)),
        "p25": float(np.quantile(values, 0.25)),
        "median": float(np.median(values)),
        "p75": float(np.quantile(values, 0.75)),
        "p90": float(np.quantile(values, 0.90)),
        "maximum": int(values.max()),
        "zero_count_centers": int(np.count_nonzero(values == 0)),
        "centers_ge_1000": int(np.count_nonzero(values >= 1000)),
        "centers_ge_3000": int(np.count_nonzero(values >= 3000)),
        "centers_ge_5000": int(np.count_nonzero(values >= 5000)),
    }


def run(plan_dir: Path, results_dir: Path, output: Path) -> dict:
    plan_dir, results_dir, output = map(Path.resolve, (plan_dir, results_dir, output))
    if output.exists():
        raise FileExistsError(f"refusing to overwrite {output}")
    plan, plan_hash = validate_plan(plan_dir)
    checkpoint_path = results_dir / "count_status.json"
    summary_path = results_dir / "count_summary.json"
    if not checkpoint_path.is_file() or not summary_path.is_file():
        raise FileNotFoundError("count-only checkpoint or summary is missing")
    state = json.loads(checkpoint_path.read_text(encoding="utf-8"))
    count_summary = json.loads(summary_path.read_text(encoding="utf-8"))
    if (state.get("status") != "complete"
            or state.get("plan_sha256") != plan_hash
            or state.get("summary_sha256") != sha256(summary_path)
            or count_summary.get("status") != "complete"
            or count_summary.get("query_plan_sha256") != plan_hash):
        raise ValueError("count run is incomplete or its plan/summary hashes disagree")
    for key in ("credentials_saved", "source_rows_downloaded", "images_returned",
                "photoz_returned", "science_measurement"):
        if state.get(key) is not False:
            raise ValueError(f"count checkpoint violates the approved scope: {key}")
    for key in ("source_rows_downloaded", "images_returned", "photoz_returned",
                "science_measurement"):
        if count_summary.get(key) is not False:
            raise ValueError(f"count summary violates the approved scope: {key}")

    centers = plan["centers"]
    expected_ids = [str(row["lens_id"]) for row in centers]
    if len(expected_ids) != 427 or len(set(expected_ids)) != 427:
        raise ValueError("the audited plan no longer contains 427 unique centers")
    jobs = state.get("jobs", {})
    if len(jobs) != len(plan["batches"]):
        raise ValueError("checkpoint batch count differs from the approved plan")

    counts: dict[str, int] = {}
    count_hashes: dict[str, str] = {}
    for batch in plan["batches"]:
        stem = Path(batch["name"]).stem
        item = jobs.get(stem, {})
        path = results_dir / f"{stem}.csv"
        if (item.get("status") != "complete" or not path.is_file()
                or sha256(path) != item.get("sha256")):
            raise ValueError(f"completed count batch or checksum is missing: {stem}")
        sql_path = plan_dir / batch["name"]
        batch_ids = re.findall(
            r"'((?:S23B_[0-9]+))'\s+AS\s+LENS_ID",
            sql_path.read_text(encoding="utf-8"), re.IGNORECASE,
        )
        if len(batch_ids) != int(batch["n_lenses"]):
            raise ValueError(f"approved SQL IDs do not match the batch manifest: {stem}")
        rows = read_count_csv(path, batch_ids)
        if item.get("rows") != len(rows) or len(rows) != int(batch["n_lenses"]):
            raise ValueError(f"row count differs from the plan batch: {stem}")
        for row in rows:
            lens_id = str(row["lens_id"])
            if lens_id in counts:
                raise ValueError(f"duplicate center in count batches: {lens_id}")
            counts[lens_id] = int(row["n_shape_sources"])
        count_hashes[path.name] = sha256(path)

    if set(counts) != set(expected_ids) or len(counts) != 427:
        raise ValueError("verified batches do not cover the 427 approved centers exactly")
    values = np.asarray([counts[lens_id] for lens_id in expected_ids], dtype=np.int64)
    if (len(values) != int(count_summary.get("n_count_rows", -1))
            or int(values.sum()) != int(count_summary.get("sum_lens_cone_shape_counts", -1))):
        raise ValueError("recomputed count totals differ from the saved aggregate summary")

    manifest_path = Path(plan["source_target_manifest"]).resolve()
    if sha256(manifest_path) != plan["source_target_manifest_sha256"]:
        raise ValueError("source target manifest differs from approved provenance")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    clusters = {str(row["id"]): row for row in manifest["clusters"]}
    joined = []
    for center in centers:
        lens_id = str(center["lens_id"])
        cluster = clusters.get(lens_id)
        if cluster is None or str(cluster["catalog_id"]) != str(center["catalog_id"]):
            raise ValueError(f"approved center does not match the target manifest: {lens_id}")
        z = float(center["z_lens"])
        if abs(z - float(cluster["z_cl"])) > 5e-7:
            raise ValueError(f"redshift mismatch for {lens_id}")
        joined.append({
            "lens_id": lens_id,
            "catalog_id": str(center["catalog_id"]),
            "z_lens": z,
            "richness": float(cluster["richness"]),
            "n_shape_sources": counts[lens_id],
        })

    count_array = np.asarray([row["n_shape_sources"] for row in joined], dtype=float)
    redshift = np.asarray([row["z_lens"] for row in joined], dtype=float)
    richness = np.asarray([row["richness"] for row in joined], dtype=float)
    rho_z = float(spearmanr(count_array, redshift).statistic)
    rho_richness = float(spearmanr(count_array, richness).statistic)
    if not np.isfinite(rho_z) or not np.isfinite(rho_richness):
        raise ValueError("descriptive rank correlation is non-finite")

    output.mkdir(parents=True, exist_ok=False)
    table_path = output / "per_center_counts.tsv"
    fields = ["lens_id", "catalog_id", "z_lens", "richness", "n_shape_sources"]
    with table_path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t",
                                lineterminator="\n", extrasaction="raise")
        writer.writeheader()
        writer.writerows(joined)

    result = {
        "schema": "astrafek-s23b-wl-shape-count-audit-v1",
        "status": "complete_exploratory_count_support_audit",
        "science_measurement": False,
        "count_only": True,
        "source_rows_downloaded": False,
        "photoz_or_shape_values_returned": False,
        "images_returned": False,
        "n_centers": len(values),
        "n_batches_verified": len(plan["batches"]),
        "sum_lens_cone_shape_counts": int(values.sum()),
        "per_center_count_distribution": summarize(values),
        "spearman_descriptive_only": {
            "count_vs_redshift_rho": rho_z,
            "count_vs_CAMIRA_richness_rho": rho_richness,
            "p_values_reported": False,
        },
        "interpretation": [
            "The count covers entries in the released S19A weak-lensing shape catalog (whose release selection already applies the catalog weak-lensing flag and star mask) inside each 3 proper-Mpc cone, before this project applies any additional source-redshift, tomography, B-mode, or radial cuts.",
            "Summed counts are lens-source candidate pairs; a source may be counted for more than one center and the sum is not a unique-source total.",
            "Redshift/richness correlations are descriptive and mix angular aperture, footprint/depth, and catalog selection; they are not physical trends.",
            "These counts do not estimate final effective source density, lensing signal-to-noise, or a calibrated weak-lensing measurement.",
            "Any source-row export, matched-random count query, or new catalog transfer needs its own approved scope.",
        ],
        "inputs": {
            "plan_sha256": plan_hash,
            "count_summary_sha256": sha256(summary_path),
            "checkpoint_sha256": sha256(checkpoint_path),
            "source_target_manifest_sha256": sha256(manifest_path),
            "count_batch_sha256": count_hashes,
        },
        "output_table": str(table_path),
        "output_table_sha256": sha256(table_path),
        "audit_script_sha256": sha256(Path(__file__).resolve()),
    }
    report_path = output / "assessment.md"
    report_path.write_text(render_assessment(result), encoding="utf-8")
    result["assessment_path"] = str(report_path)
    result["assessment_sha256"] = sha256(report_path)
    summary_output = output / "summary.json"
    summary_output.write_text(json.dumps(result, indent=2, sort_keys=True,
                                         ensure_ascii=False, allow_nan=False) + "\n",
                              encoding="utf-8")
    result["summary_sha256"] = sha256(summary_output)
    return result


def render_assessment(result: dict) -> str:
    distribution = result["per_center_count_distribution"]
    lines = [
        "# S23B S19A weak-lensing count-only support audit",
        "",
        "## Scope",
        "",
        "This local audit verified all 18 batches from the approved count-only plan. It read 427 integer counts and joined them to the already saved target manifest for descriptive redshift and richness summaries. No CAS request was made, and no source rows, coordinates, shape values, photo-z, or images were returned.",
        "",
        "The counts are `object_id` totals from the released S19A weak-lensing shape-catalog subset within 3 proper Mpc. The release catalog has already applied its weak-lensing flag and star mask; this project has not yet applied additional source-redshift, tomographic, B-mode, or radial cuts. They are not a weak-lensing measurement or final effective source counts.",
        "",
        "## Counts",
        "",
        f"- Centers: {result['n_centers']}; verified batches: {result['n_batches_verified']}",
        f"- Lens-source candidate pairs (sum of per-center counts): {result['sum_lens_cone_shape_counts']:,}",
        f"- Per-center count: min {distribution['minimum']:,}, p10 {distribution['p10']:,.1f}, median {distribution['median']:,.1f}, p90 {distribution['p90']:,.1f}, max {distribution['maximum']:,}",
        f"- Zero-count centers: {distribution['zero_count_centers']}; centers >=1,000: {distribution['centers_ge_1000']}; >=3,000: {distribution['centers_ge_3000']}; >=5,000: {distribution['centers_ge_5000']}",
        f"- Descriptive Spearman rho: count versus redshift {result['spearman_descriptive_only']['count_vs_redshift_rho']:.3f}; versus CAMIRA richness {result['spearman_descriptive_only']['count_vs_CAMIRA_richness_rho']:.3f}. No p-values are reported.",
        "",
        "## Interpretation",
        "",
        "The 3 proper-Mpc count totals are an availability screen only. A source can be counted around multiple centers, and the total is a pair count rather than a unique-source count. The HSC release weak-lensing flag and star mask are already built into this catalog table, but no project-specific source-redshift, tomography, B-mode, or radial selections have been applied. Do not infer lensing signal-to-noise or promote `science_measurement`; any source-row export or additional query requires separate approval.",
        "",
        "## Provenance",
        "",
        f"- Approved plan SHA-256: `{result['inputs']['plan_sha256']}`",
        f"- Count summary SHA-256: `{result['inputs']['count_summary_sha256']}`",
        f"- Audit runner SHA-256: `{result['audit_script_sha256']}`",
        f"- Per-center table SHA-256: `{result['output_table_sha256']}`",
        "- HSC release selection and table contract: [S19A Shape Catalog (PDR3)](https://hsc-release.mtk.nao.ac.jp/doc/index.php/s19a-shape-catalog-pdr3/)",
        "- `science_measurement=false`.",
        "",
    ]
    return "\n".join(lines)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan-dir", type=Path, default=DEFAULT_PLAN)
    parser.add_argument("--results-dir", type=Path, default=DEFAULT_RESULTS)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    result = run(args.plan_dir, args.results_dir, args.output)
    print(json.dumps({
        "status": result["status"],
        "science_measurement": False,
        "centers": result["n_centers"],
        "batches_verified": result["n_batches_verified"],
        "candidate_pairs": result["sum_lens_cone_shape_counts"],
        "zero_count_centers": result["per_center_count_distribution"]["zero_count_centers"],
        "output": str(args.output.resolve()),
    }, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
