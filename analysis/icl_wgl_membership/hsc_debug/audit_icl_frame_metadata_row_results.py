#!/usr/bin/env python3
"""Offline integrity and matched-center summary for approved HSC frame metadata."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from prepare_icl_frame_metadata_rows import load_audited_centers
from submit_icl_frame_metadata_rows import EXPECTED_PLAN_SHA256, sha256, validate_plan

ROOT = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
DEFAULT_PLAN = ROOT / "icl_frame_metadata_row_plan_2026-11-05_01"
DEFAULT_RESULTS = ROOT / "icl_frame_metadata_row_results_2026-11-06_01"
DEFAULT_OUTPUT = ROOT / "icl_frame_metadata_row_audit_2026-11-06_01"
CONTROL_MANIFEST = ROOT / "control_download_01/download_manifest.json"
FIELDS = [
    "center_id", "sample_role", "ra_deg", "dec_deg", "coadd_band", "skymap_id",
    "tract", "patch", "patch_s", "n_frame_rows", "n_unique_frames", "n_unique_visits",
    "n_unique_ccds", "sum_exptime_s", "median_sigma_sky", "median_fluxmag0",
    "median_gain1", "median_gain2", "median_gain3", "median_gain4", "n_filter_i", "n_filter_i2",
]


def number(value: str) -> float | None:
    text = str(value).strip()
    if text.lower() in ("", "null", "none", "nan", r"\n", r"\N"):
        return None
    result = float(text)
    if not np.isfinite(result):
        raise ValueError("non-finite frame metadata value")
    return result


def read_csv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream)
        names = [name.strip().removeprefix("#").strip() for name in (reader.fieldnames or [])]
        if len(names) != 21 or len(set(names)) != 21:
            raise ValueError(f"unexpected metadata CSV columns in {path}")
        reader.fieldnames = names
        rows = []
        for raw in reader:
            row = {str(key).strip().removeprefix("#").strip(): value for key, value in raw.items()}
            if None in row or any(value is None for value in row.values()):
                raise ValueError(f"malformed metadata row in {path}")
            rows.append(row)
    return names, rows


def summarize(values: list[float]) -> dict:
    x = np.asarray(values, dtype=float)
    if x.size == 0:
        return {"n": 0, "median": None, "q16": None, "q84": None, "mean": None}
    return {"n": int(x.size), "median": float(np.median(x)),
            "q16": float(np.quantile(x, 0.16)), "q84": float(np.quantile(x, 0.84)),
            "mean": float(np.mean(x))}


def run(plan_dir: Path, results_dir: Path, output_dir: Path) -> dict:
    if output_dir.exists():
        raise FileExistsError(output_dir)
    plan, centers, batches, plan_hash = validate_plan(plan_dir)
    if plan_hash != EXPECTED_PLAN_SHA256:
        raise ValueError("audit plan hash differs from pinned approval")
    status_path = results_dir / "status.json"
    summary_path = results_dir / "summary.json"
    status = json.loads(status_path.read_text(encoding="utf-8"))
    upstream = json.loads(summary_path.read_text(encoding="utf-8"))
    if (status.get("status") != "complete" or status.get("plan_sha256") != plan_hash
            or status.get("query_submitted") is not True or status.get("credentials_saved") is not False
            or upstream.get("status") != "complete" or upstream.get("plan_sha256") != plan_hash
            or upstream.get("n_centers") != 1598 or upstream.get("n_frame_rows") != 67071
            or upstream.get("images_returned") is not False
            or upstream.get("shape_rows_returned") is not False
            or upstream.get("photoz_rows_returned") is not False
            or upstream.get("credentials_saved") is not False
            or upstream.get("science_measurement") is not False):
        raise ValueError("CAS result checkpoint/summary does not meet the approved export scope")

    expected = {str(row["center_id"]): row for row in centers}
    if len(expected) != 1598:
        raise ValueError("audited parent center list changed")
    seen_frames: dict[str, set[str]] = {key: set() for key in expected}
    rows_by_center: dict[str, list[dict[str, str]]] = {key: [] for key in expected}
    batch_hashes = {}
    total_rows = 0
    for batch in batches:
        name = batch["name"]
        target = results_dir / f"{name}.csv"
        job = status.get("jobs", {}).get(name, {})
        if (job.get("status") != "complete" or not target.is_file()
                or sha256(target) != job.get("sha256")):
            raise ValueError(f"missing, incomplete, or changed approved batch: {name}")
        columns, batch_rows = read_csv(target)
        if columns != plan["output_columns"]:
            raise ValueError(f"returned field list differs from approved 21-column contract: {name}")
        for row in batch_rows:
            ident = row["center_id"]
            parent = expected.get(ident)
            if parent is None:
                raise ValueError(f"CAS returned an unrequested center: {ident}")
            if (row["sample_role"] != parent["role"] or row["coadd_band"] != parent["band"]
                    or row["skymap_id"] != str(parent["skymap_id"])
                    or row["tract"] != str(parent["tract"])
                    or row["patch"] != str(parent["patch"])
                    or row["patch_s"] != str(parent["patch_s"])
                    or abs(float(row["ra_deg"]) - float(parent["ra_deg"])) > 1e-8
                    or abs(float(row["dec_deg"]) - float(parent["dec_deg"])) > 1e-8):
                raise ValueError(f"CAS center/patch identity mismatch: {ident}")
            if row["filter"] not in {"i", "i2"}:
                raise ValueError(f"unexpected original frame filter: {row['filter']}")
            frame = row["frame_num"].strip()
            if not frame or frame in seen_frames[ident]:
                raise ValueError(f"duplicate/empty frame membership for {ident}")
            seen_frames[ident].add(frame)
            rows_by_center[ident].append(row)
        if len(batch_rows) != int(batch["expected_rows_from_approved_count_query"]):
            raise ValueError(f"row total differs from the approved count-only parent: {name}")
        if job.get("rows") != len(batch_rows):
            raise ValueError(f"job row count differs from its returned file: {name}")
        total_rows += len(batch_rows)
        batch_hashes[name] = job["sha256"]
    if total_rows != 67071 or any(not rows_by_center[key] for key in expected):
        raise ValueError("returned rows do not reproduce all 1,598 audited centers")
    for ident, parent in expected.items():
        if len(rows_by_center[ident]) != int(parent["n_patch_frame_memberships"]):
            raise ValueError(f"frame membership count differs from the parent audit: {ident}")

    aggregate_rows = []
    lookup = {}
    for ident, parent in expected.items():
        selected = rows_by_center[ident]
        metrics = {
            "center_id": ident, "sample_role": parent["role"], "ra_deg": parent["ra_deg"],
            "dec_deg": parent["dec_deg"], "coadd_band": parent["band"],
            "skymap_id": parent["skymap_id"], "tract": parent["tract"],
            "patch": parent["patch"], "patch_s": parent["patch_s"],
            "n_frame_rows": len(selected),
            "n_unique_frames": len({row["frame_num"] for row in selected}),
            "n_unique_visits": len({row["visit"] for row in selected}),
            "n_unique_ccds": len({row["ccd"] for row in selected}),
            "sum_exptime_s": float(sum(number(row["exptime"]) or 0.0 for row in selected)),
            "median_sigma_sky": float(np.median([number(row["sigma_sky"]) for row in selected
                                                  if number(row["sigma_sky"]) is not None])),
            "median_fluxmag0": float(np.median([number(row["fluxmag0"]) for row in selected
                                                 if number(row["fluxmag0"]) is not None])),
            "n_filter_i": sum(row["filter"] == "i" for row in selected),
            "n_filter_i2": sum(row["filter"] == "i2" for row in selected),
        }
        for amp in range(1, 5):
            values = [number(row[f"gain{amp}"]) for row in selected]
            metrics[f"median_gain{amp}"] = float(np.median([value for value in values if value is not None]))
        aggregate_rows.append(metrics)
        lookup[ident] = metrics

    control_doc = json.loads(CONTROL_MANIFEST.read_text(encoding="utf-8"))
    pairs = []
    for row in control_doc["clusters"]:
        target_id = str(row["target_cluster_id"])
        control = str(row["id"])
        target_key, control_key = f"target:{target_id}", f"control:{control}"
        if target_key not in lookup or control_key not in lookup:
            continue
        a, b = lookup[target_key], lookup[control_key]
        pair = {"target_cluster_id": target_id, "control_id": control,
                "target_center_id": target_key, "control_center_id": control_key}
        for metric in ("n_frame_rows", "n_unique_visits", "n_unique_ccds", "sum_exptime_s",
                       "median_sigma_sky", "median_fluxmag0", "median_gain1", "median_gain2",
                       "median_gain3", "median_gain4", "n_filter_i", "n_filter_i2"):
            pair[f"target_minus_control_{metric}"] = float(a[metric]) - float(b[metric])
        pairs.append(pair)
    if len(pairs) != 798:
        raise ValueError(f"expected 798 audited target/control pairs; found {len(pairs)}")

    output_dir.mkdir(parents=True)
    center_path = output_dir / "center_frame_metadata_aggregates.tsv"
    pair_path = output_dir / "matched_target_control_metadata_differences.tsv"
    for path, fields, rows in ((center_path, FIELDS, aggregate_rows),
                               (pair_path, list(pairs[0]), pairs)):
        with path.open("x", encoding="utf-8", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(sorted(rows, key=lambda item: item.get("center_id", item.get("target_cluster_id"))))

    metrics = ["n_frame_rows", "n_unique_visits", "n_unique_ccds", "sum_exptime_s",
               "median_sigma_sky", "median_fluxmag0", "median_gain1", "median_gain2",
               "median_gain3", "median_gain4", "n_filter_i", "n_filter_i2"]
    differences = {metric: summarize([row[f"target_minus_control_{metric}"] for row in pairs])
                   for metric in metrics}
    result = {
        "schema": "astrafek-icl-frame-metadata-row-audit-v1", "status": "complete_exploratory",
        "science_measurement": False, "network_used_by_auditor": False,
        "approved_query_submitted": True, "plan_sha256": plan_hash,
        "n_centers": len(expected), "n_target_centers": 800, "n_control_centers": 798,
        "n_matched_target_control_pairs": len(pairs), "n_frame_rows": total_rows,
        "source_frame_filter_counts": {name: sum(len([r for r in rows if r["filter"] == name])
                                                   for rows in rows_by_center.values())
                                        for name in ("i", "i2")},
        "metadata_null_counts": upstream["null_counts_by_metadata_column"],
        "matched_target_minus_control_metadata_differences": differences,
        "not_returned": {"images": True, "shape_rows": True, "photoz_rows": True,
                         "object_catalog_rows": True, "credentials": True, "pixel_coadd_weights": True},
        "result_files_sha256": {
            "status.json": sha256(results_dir / "status.json"),
            "summary.json": sha256(summary_path),
            "batch_hashes_sha256": hashlib.sha256("\n".join(
                f"{name}:{digest}" for name, digest in sorted(batch_hashes.items())
            ).encode()).hexdigest(),
            "center_aggregates.tsv": sha256(center_path),
            "matched_differences.tsv": sha256(pair_path),
        },
        "aggregate_tables": {"center_aggregates": str(center_path),
                             "matched_target_control_differences": str(pair_path)},
        "interpretation_limit": "Frame-level metadata summarize contributing visits and CCDs; without pixel-to-frame weights and amplifier maps, they do not establish per-pixel depth, coadd variance, or source-Poisson correction.",
    }
    report = output_dir / "summary.json"
    report.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, ensure_ascii=False, indent=2))
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan-dir", type=Path, default=DEFAULT_PLAN)
    parser.add_argument("--results", type=Path, default=DEFAULT_RESULTS)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    run(args.plan_dir, args.results, args.output)


if __name__ == "__main__":
    main()
