#!/usr/bin/env python3
"""Audit and summarize the exact approved ICL random-point variance export.

This local postprocessor is read-only with respect to HSC. It validates the
approved CAS checkpoints and each returned row, then writes center-level and
matched target/control summaries for the 200--300 proper-kpc annulus.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from fetch_ancillary_catalogues import catalogue_reader
from submit_icl_random_pixel_variance import (
    DEFAULT_PLAN, DEFAULT_RESULTS, EXPECTED_PLAN_SHA256, OUTPUT_COLUMNS,
    load_plan, sha256, validate_rows,
)

ROOT = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
CONTROL_MANIFEST = ROOT / "control_download_01/download_manifest.json"
DEFAULT_OUTPUT = ROOT / "icl_random_pixel_variance_audit_2026-11-06_01"
NUMERIC_FIELDS = (
    "i_pix_variance", "i_sky_mean", "i_sky_std", "i_inputcount_value",
    "i_sdssshape_psf_shape11", "i_sdssshape_psf_shape22",
    "i_sdssshape_psf_shape12",
)
MASK_FIELDS = (
    "i_pixelflags_interpolatedcenter", "i_pixelflags_saturatedcenter",
    "i_pixelflags_crcenter", "i_pixelflags_suspectcenter",
    "i_pixelflags_bright_objectcenter", "i_pixelflags_edge",
)


def finite_number(value: str) -> float | None:
    text = str(value).strip()
    if text.lower() in ("", "null", "none", "nan", r"\n", r"\N"):
        return None
    number = float(text)
    if not math.isfinite(number):
        return None
    return number


def truth(value: str) -> bool:
    text = str(value).strip().lower()
    if text in ("true", "t", "1", "yes"):
        return True
    if text in ("false", "f", "0", "no"):
        return False
    raise ValueError(f"CAS returned invalid boolean mask value: {value!r}")


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def summarize(values: list[float]) -> dict:
    array = np.asarray(values, dtype=float)
    array = array[np.isfinite(array)]
    if array.size == 0:
        return {"n": 0, "median": None, "mean": None, "q16": None, "q84": None}
    return {
        "n": int(array.size), "median": float(np.median(array)),
        "mean": float(np.mean(array)), "q16": float(np.quantile(array, 0.16)),
        "q84": float(np.quantile(array, 0.84)),
    }


def aggregate_center(center_id: str, center: dict,
                     rows: list[dict[str, str]]) -> dict:
    result = {
        "center_id": center_id, "sample_role": center["role"],
        "target_cluster_id": center["target_cluster_id"],
        "n_random_points": len(rows),
        "ra_deg": center["ra_deg"], "dec_deg": center["dec_deg"],
        "cluster_z_for_scale": center["z"],
    }
    for field in NUMERIC_FIELDS:
        values = [value for row in rows
                  if (value := finite_number(row[field])) is not None]
        stats = summarize(values)
        result[f"n_{field}"] = stats["n"]
        for statistic in ("median", "mean", "q16", "q84"):
            result[f"{statistic}_{field}"] = stats[statistic]
    for field in MASK_FIELDS:
        result[f"n_true_{field}"] = sum(truth(row[field]) for row in rows)
        result[f"fraction_true_{field}"] = (
            result[f"n_true_{field}"] / len(rows) if rows else None
        )
    return result


def run(plan_dir: Path, results_dir: Path, output: Path) -> dict:
    plan_dir, results_dir, output = map(Path, (plan_dir, results_dir, output))
    output = output.resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite {output}")
    plan, centers, batches, plan_hash = load_plan(plan_dir)
    if plan_hash != EXPECTED_PLAN_SHA256:
        raise ValueError("approved plan SHA-256 changed")
    status_path = results_dir / "status.json"
    upstream_summary_path = results_dir / "summary.json"
    if not status_path.is_file() or not upstream_summary_path.is_file():
        raise FileNotFoundError("completed approved random-point export is not present")
    status, upstream = read_json(status_path), read_json(upstream_summary_path)
    if (status.get("status") != "complete" or status.get("plan_sha256") != plan_hash
            or status.get("query_submitted") is not True
            or status.get("credentials_saved") is not False
            or upstream.get("status") != "complete_exploratory"
            or upstream.get("plan_sha256") != plan_hash
            or upstream.get("n_requested_centers") != 1598
            or upstream.get("images_returned") is not False
            or upstream.get("shape_rows_returned") is not False
            or upstream.get("photoz_rows_returned") is not False
            or upstream.get("object_catalog_rows_returned") is not False
            or upstream.get("credentials_saved") is not False
            or upstream.get("science_measurement") is not False):
        raise ValueError("CAS checkpoint or result summary violates approved scope")

    control_doc = read_json(CONTROL_MANIFEST)
    controls = {str(row["target_cluster_id"]): row for row in control_doc["clusters"]}
    if len(controls) != 798:
        raise ValueError("expected 798 matched target/control center mappings")
    center_by_id = {str(row["center_id"]): row for row in centers}
    rows_by_center: dict[str, list[dict[str, str]]] = {key: [] for key in center_by_id}
    batch_hashes, total_rows = {}, 0
    unique_ids: dict[str, set[str]] = {key: set() for key in center_by_id}
    for batch in batches:
        name = batch["name"]
        csv_path = results_dir / f"{name}.csv"
        job = status.get("jobs", {}).get(name, {})
        expected = [center_by_id[key] for key in batch["center_ids"]]
        if (job.get("status") != "complete" or not csv_path.is_file()
                or sha256(csv_path) != job.get("sha256")):
            raise ValueError(f"missing, incomplete, or changed approved batch {name}")
        row_count, returned_centers, _ = validate_rows(csv_path, expected)
        if row_count != job.get("rows"):
            raise ValueError(f"returned row count differs from checkpoint for {name}")
        with csv_path.open(encoding="utf-8-sig", newline="") as stream:
            reader = catalogue_reader(stream)
            if reader.fieldnames != OUTPUT_COLUMNS:
                raise ValueError(f"column contract changed in {name}")
            batch_rows = list(reader)
        for row in batch_rows:
            center_id = row["center_id"]
            random_id = row["random_object_id"]
            for field in NUMERIC_FIELDS:
                finite_number(row[field])
            for field in MASK_FIELDS:
                truth(row[field])
            if random_id in unique_ids[center_id]:
                raise ValueError(f"duplicate random ID within center {center_id}")
            unique_ids[center_id].add(random_id)
            rows_by_center[center_id].append(row)
        if len(returned_centers) != len({row["center_id"] for row in batch_rows}):
            raise ValueError(f"center ledger mismatch in {name}")
        total_rows += row_count
        batch_hashes[name] = job["sha256"]

    if (total_rows != upstream.get("n_returned_random_points")
            or len(batch_hashes) != 16):
        raise ValueError("batch totals differ from the completed export summary")
    centers_by_id = {key: aggregate_center(key, center_by_id[key], rows)
                     for key, rows in rows_by_center.items()}
    center_rows = list(centers_by_id.values())
    paired_rows = []
    for target_id, control in controls.items():
        target_key, control_key = f"target:{target_id}", f"control:{control['id']}"
        if target_key not in centers_by_id or control_key not in centers_by_id:
            raise ValueError(f"missing target/control center in approved plan for {target_id}")
        target, random_control = centers_by_id[target_key], centers_by_id[control_key]
        pair = {"target_cluster_id": target_id, "control_id": str(control["id"]),
                "target_n_random_points": target["n_random_points"],
                "control_n_random_points": random_control["n_random_points"]}
        for field in NUMERIC_FIELDS:
            a, b = target[f"median_{field}"], random_control[f"median_{field}"]
            pair[f"target_median_{field}"] = a
            pair[f"control_median_{field}"] = b
            pair[f"target_minus_control_median_{field}"] = (
                float(a) - float(b) if a is not None and b is not None else None
            )
        for field in MASK_FIELDS:
            pair[f"target_fraction_true_{field}"] = target[f"fraction_true_{field}"]
            pair[f"control_fraction_true_{field}"] = random_control[f"fraction_true_{field}"]
        paired_rows.append(pair)
    if len(paired_rows) != 798:
        raise ValueError(f"expected 798 target/control pairs, got {len(paired_rows)}")

    output.mkdir(parents=True, exist_ok=False)
    center_path, pair_path = output / "center_random_point_summary.tsv", output / "matched_target_control_summary.tsv"
    for path, rows in ((center_path, center_rows), (pair_path, paired_rows)):
        with path.open("x", encoding="utf-8", newline="") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t",
                                    lineterminator="\n", extrasaction="raise")
            writer.writeheader()
            writer.writerows(rows)

    n_with_points = sum(row["n_random_points"] > 0 for row in center_rows)
    n_target_with = sum(row["n_random_points"] > 0 for row in center_rows if row["sample_role"] == "target")
    n_control_with = sum(row["n_random_points"] > 0 for row in center_rows if row["sample_role"] == "control")
    valid_pair_counts = {
        field: sum(row[f"target_median_{field}"] is not None
                   and row[f"control_median_{field}"] is not None for row in paired_rows)
        for field in NUMERIC_FIELDS
    }
    summary = {
        "schema": "astrafek-icl-random-point-variance-audit-v1",
        "status": "complete_exploratory_audit", "science_measurement": False,
        "network_used_by_auditor": False, "plan_sha256": plan_hash,
        "n_returned_random_points": total_rows, "n_planned_centers": 1598,
        "n_centers_with_points": n_with_points,
        "n_target_centers_with_points": n_target_with,
        "n_control_centers_with_points": n_control_with,
        "n_matched_target_control_pairs": len(paired_rows),
        "matched_pairs_with_both_endpoint_statistics": valid_pair_counts,
        "random_density_cut": "adjust_density < 0.2",
        "annulus_proper_kpc": [200.0, 300.0],
        "input_sha256": {
            "plan": sha256(plan_dir / "plan.json"),
            "status": sha256(status_path), "upstream_summary": sha256(upstream_summary_path),
            "control_manifest": sha256(CONTROL_MANIFEST),
            "batch_hashes_sha256": hashlib.sha256("\n".join(
                f"{key}:{value}" for key, value in sorted(batch_hashes.items())
            ).encode()).hexdigest(),
        },
        "not_returned": {"images": True, "shape_rows": True, "photoz_rows": True,
                         "object_catalog_rows": True, "credentials": True},
        "limitations": [
            "Random-point pixel variance is a sparse sample, not a pixel variance or coadd weight map.",
            "Patch-level frame metadata and random points do not identify pixel-to-frame or amplifier weights.",
            "The audit does not estimate correlated noise or source-Poisson variance and applies no ICL correction.",
            "A zero-point sample count can leave a center or matched pair without a valid endpoint statistic.",
        ],
        "outputs": {
            center_path.name: sha256(center_path), pair_path.name: sha256(pair_path),
        },
    }
    (output / "summary.json").write_text(
        json.dumps(summary, ensure_ascii=False, indent=2, allow_nan=False) + "\n",
        encoding="utf-8")
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan-dir", type=Path, default=DEFAULT_PLAN)
    parser.add_argument("--results-dir", type=Path, default=DEFAULT_RESULTS)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    summary = run(args.plan_dir, args.results_dir, args.output)
    print(json.dumps({key: summary[key] for key in (
        "status", "science_measurement", "plan_sha256", "n_returned_random_points",
        "n_centers_with_points", "n_target_centers_with_points",
        "n_control_centers_with_points", "matched_pairs_with_both_endpoint_statistics",
    )}, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
