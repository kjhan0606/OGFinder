#!/usr/bin/env python3
"""Verify and compare approved S16A counts at 116 saved ICL control centers.

Local-only QA: the script reads count-only CAS outputs and the hash-pinned
target/control-center table. It reads no source-level rows, shapes, photo-z, or
image pixels and makes no network or CAS request.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.stats import spearmanr

DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
PLAN_DIR = DATA / "wl_s23b_s16a_icl_control_count_plan_2026-11-11_01"
RESULTS = PLAN_DIR / "results"
TARGET_PLAN_RESULTS = DATA / "wl_s23b_s16a_shape_count_plan_2026-11-11_01" / "results"
SELECTION_AUDIT = DATA / "s23b_s16a_shape_coverage_selection_audit_2026-11-11_04"
CONTROL_TABLE = SELECTION_AUDIT / "s16a_ge1000_outside_s19a_icl_control_centers.tsv"
DEFAULT_OUTPUT = DATA / "s23b_s16a_icl_control_count_audit_2026-11-11_01"

EXPECTED_PLAN_SHA256 = "832f2de18f5d512396cc84673e458932635c1f766ea91cdc04c7b47125738dfa"
EXPECTED_CONTROL_TABLE_SHA256 = "81c702166871e75570918f770ffad44702cebe9a29920c80d3d73d2263c21706"
EXPECTED_CONTROL_SUMMARY_SHA256 = "136b86970f15f0d96118420a6ffe9ca9d510ad603b16cdfa28063fb55f4386f2"
EXPECTED_CONTROL_RESULTS_SHA256 = "b72d15093bdc3f4ddc6e03a1f85f9911b8999ed5e463e50c1eb7d6d4d4edbcd7"
EXPECTED_TARGET_RESULTS_SHA256 = "0945329119d1100b783211fc128cd0a329b449f5acf6648313a1484ac55a7b12"
EXPECTED_TARGET_SUMMARY_SHA256 = "0cf90929c2a69eda0e62f3a9350359fd430b37bf869b2ead53246129ac51bb8b"
EXPECTED_AUDIT_SUMMARY_SHA256 = "c08602e1b45492a92389dabd074e5f854050e8e14c2cd55b4683b5ae6fc0b398"
EXPECTED_AUDIT_CONTROL_TABLE_SHA256 = EXPECTED_CONTROL_TABLE_SHA256
EXPECTED_CONTROL_ROWS = 116
EXPECTED_BATCH_SIZE = 25
COUNT_COLUMNS = ["control_id", "n_shape_sources"]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_json(path: Path) -> dict:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def read_tsv(path: Path) -> list[dict[str, str]]:
    with Path(path).open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def read_count_csv(path: Path, expected_ids: list[str]) -> list[dict[str, int | str]]:
    with Path(path).open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.reader(stream)
        header = next(reader)
        header[0] = header[0].lstrip("# ").strip()
        if header != COUNT_COLUMNS:
            raise ValueError(f"unexpected result columns in {path.name}: {header}")
        rows = []
        for values in reader:
            if len(values) != 2:
                raise ValueError(f"malformed result row in {path.name}")
            control_id, raw_count = values
            try:
                count = int(raw_count)
            except ValueError as error:
                raise ValueError(f"non-integer count in {path.name}") from error
            if not control_id or count < 0:
                raise ValueError(f"invalid result row in {path.name}")
            rows.append({"control_id": control_id, "n_shape_sources": count})
    ids = [str(row["control_id"]) for row in rows]
    if len(rows) != len(expected_ids) or len(set(ids)) != len(ids) or set(ids) != set(expected_ids):
        raise ValueError(f"result IDs differ from the approved batch in {path.name}")
    return rows


def summarize(values: list[float | int]) -> dict:
    array = np.asarray(values, dtype=float)
    if not array.size:
        return {"n": 0, "min": None, "q16": None, "median": None, "q84": None, "max": None}
    return {
        "n": int(array.size),
        "min": float(np.min(array)),
        "q16": float(np.quantile(array, 0.16)),
        "median": float(np.median(array)),
        "q84": float(np.quantile(array, 0.84)),
        "max": float(np.max(array)),
    }


def validate_inputs() -> tuple[dict, list[dict], list[dict], dict, dict]:
    plan_path = PLAN_DIR / "shape_count_plan.json"
    if sha256(plan_path) != "78dc14cbdedcd3bc0804fa381f57420b527226720e75bab028df59b4333418f2":
        raise ValueError("approved control-count plan JSON changed")
    plan = read_json(plan_path)
    if (plan.get("query_plan_sha256") != EXPECTED_PLAN_SHA256
            or plan.get("status") != "prepared_not_submitted"
            or plan.get("query_submitted") is not False
            or plan.get("n_controls") != EXPECTED_CONTROL_ROWS
            or plan.get("output_columns") != COUNT_COLUMNS):
        raise ValueError("approved plan scope or provenance changed")
    if (sha256(CONTROL_TABLE) != EXPECTED_CONTROL_TABLE_SHA256
            or plan.get("candidate_control_table_sha256") != EXPECTED_CONTROL_TABLE_SHA256
            or sha256(SELECTION_AUDIT / "summary.json") != EXPECTED_AUDIT_SUMMARY_SHA256):
        raise ValueError("hash-pinned selection-audit inputs changed")
    if (sha256(RESULTS / "count_summary.json") != EXPECTED_CONTROL_SUMMARY_SHA256
            or sha256(RESULTS / "count_results.tsv") != EXPECTED_CONTROL_RESULTS_SHA256):
        raise ValueError("completed control-count outputs changed")
    if (sha256(TARGET_PLAN_RESULTS / "count_results.tsv") != EXPECTED_TARGET_RESULTS_SHA256
            or sha256(TARGET_PLAN_RESULTS / "count_summary.json") != EXPECTED_TARGET_SUMMARY_SHA256):
        raise ValueError("approved target-count results changed")

    audit_summary = read_json(SELECTION_AUDIT / "summary.json")
    if (audit_summary.get("status") != "complete_local_descriptive_audit"
            or audit_summary.get("science_measurement") is not False
            or audit_summary.get("cas_query_submitted") is not False
            or audit_summary.get("source_rows_downloaded") is not False
            or audit_summary.get("outputs_sha256", {}).get(
                "s16a_ge1000_outside_s19a_icl_control_centers.tsv"
            ) != EXPECTED_AUDIT_CONTROL_TABLE_SHA256):
        raise ValueError("selection-audit scope changed")

    run_summary = read_json(RESULTS / "count_summary.json")
    if (run_summary.get("schema") != "astrafek-s23b-s16a-icl-control-count-summary-v1"
            or run_summary.get("status") != "complete_control_count_only"
            or run_summary.get("query_plan_sha256") != EXPECTED_PLAN_SHA256
            or run_summary.get("n_controls") != EXPECTED_CONTROL_ROWS
            or run_summary.get("n_count_rows") != EXPECTED_CONTROL_ROWS
            or run_summary.get("count_results_sha256") != EXPECTED_CONTROL_RESULTS_SHA256
            or run_summary.get("source_rows_downloaded") is not False
            or run_summary.get("images_returned") is not False
            or run_summary.get("photoz_returned") is not False
            or run_summary.get("science_measurement") is not False):
        raise ValueError("completed control-count summary scope changed")

    state = read_json(RESULTS / "count_status.json")
    jobs = state.get("jobs", {})
    if (state.get("schema") != "astrafek-s23b-s16a-icl-control-count-run-v1"
            or state.get("status") != "complete"
            or state.get("plan_sha256") != EXPECTED_PLAN_SHA256
            or state.get("total_count_rows") != EXPECTED_CONTROL_ROWS
            or state.get("credentials_saved") is not False
            or state.get("source_rows_downloaded") is not False
            or state.get("images_returned") is not False
            or state.get("photoz_returned") is not False
            or state.get("science_measurement") is not False
            or len(jobs) != 5
            or any(item.get("status") != "complete" for item in jobs.values())):
        raise ValueError("run checkpoint is incomplete or has an unexpected scope")

    control_rows = read_tsv(CONTROL_TABLE)
    count_rows = read_tsv(RESULTS / "count_results.tsv")
    if (len(control_rows) != EXPECTED_CONTROL_ROWS
            or len(count_rows) != EXPECTED_CONTROL_ROWS
            or list(count_rows[0]) != COUNT_COLUMNS):
        raise ValueError("control or count table has unexpected rows/columns")
    expected_ids = [center["control_id"] for center in plan["centers"]]
    if [row["control_id"] for row in control_rows] != expected_ids:
        raise ValueError("candidate-control table order differs from the approved plan")
    if [row["control_id"] for row in count_rows] != expected_ids:
        raise ValueError("combined CAS result order differs from the approved plan")

    batch_total = 0
    for index, batch in enumerate(plan["batches"]):
        name = Path(batch["name"]).stem
        item = jobs.get(name)
        result_path = RESULTS / f"{name}.csv"
        start = index * EXPECTED_BATCH_SIZE
        expected_batch_ids = expected_ids[start:start + EXPECTED_BATCH_SIZE]
        if (item is None or not result_path.is_file()
                or sha256(result_path) != item.get("sha256")
                or item.get("rows") != len(expected_batch_ids)):
            raise ValueError(f"batch checkpoint or result hash differs for {name}")
        batch_rows = read_count_csv(result_path, expected_batch_ids)
        if sum(int(row["n_shape_sources"]) for row in batch_rows) != item.get("total_sources"):
            raise ValueError(f"batch total differs from checkpoint for {name}")
        batch_total += len(batch_rows)
    if batch_total != EXPECTED_CONTROL_ROWS:
        raise ValueError("batch rows do not sum to 116 controls")

    target_count_map = {}
    for center in plan["centers"]:
        target_count_map[str(center["matched_lens_id"])] = int(
            center["matched_target_s16a_raw_count"]
        )
    target_count_rows = read_tsv(TARGET_PLAN_RESULTS / "count_results.tsv")
    if len(target_count_rows) != 1043 or list(target_count_rows[0]) != ["lens_id", "n_shape_sources"]:
        raise ValueError("approved S16A target counts have unexpected schema")
    observed_target_counts = {row["lens_id"]: int(row["n_shape_sources"]) for row in target_count_rows}
    if any(observed_target_counts.get(lens_id) != count for lens_id, count in target_count_map.items()):
        raise ValueError("candidate target counts disagree with the approved target-count table")
    return plan, control_rows, count_rows, run_summary, state


def analyze(output: Path) -> dict:
    output = output.resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite existing local QA output: {output}")
    plan, controls, counts, run_summary, state = validate_inputs()
    count_by_id = {row["control_id"]: int(row["n_shape_sources"]) for row in counts}
    pairs = []
    for control in controls:
        target_id = str(control["target_id"])
        control_id = str(control["control_id"])
        target_count = int(control["s16a_target_count_3pmpc"])
        control_count = count_by_id[control_id]
        if target_count < 1000:
            raise ValueError(f"candidate target fell below the selection screen: {target_id}")
        pairs.append({
            "target_id": target_id,
            "control_id": control_id,
            "icl_filter": control["icl_filter"],
            "z_cl": float(control["z_cl"]),
            "richness": float(control["richness"]),
            "s16a_target_count_3pmpc": target_count,
            "s16a_control_count_3pmpc": control_count,
            "target_minus_control_count": target_count - control_count,
            "control_to_target_count_ratio": control_count / target_count,
            "delta_ebv_sfd_target_minus_control": float(
                control["delta_ebv_sfd_target_minus_control"]
            ),
        })

    filter_summary = {}
    for band in ("i", "i2"):
        selected = [row for row in pairs if row["icl_filter"] == band]
        targets = [row["s16a_target_count_3pmpc"] for row in selected]
        controls_n = [row["s16a_control_count_3pmpc"] for row in selected]
        ratios = [row["control_to_target_count_ratio"] for row in selected]
        differences = [row["target_minus_control_count"] for row in selected]
        filter_summary[band] = {
            "n_pairs": len(selected),
            "positive_target_centers": sum(value > 0 for value in targets),
            "positive_control_centers": sum(value > 0 for value in controls_n),
            "zero_control_centers": sum(value == 0 for value in controls_n),
            "target_count": summarize(targets),
            "control_count": summarize(controls_n),
            "control_to_target_count_ratio": summarize(ratios),
            "target_minus_control_count": summarize(differences),
        }

    target_counts = [row["s16a_target_count_3pmpc"] for row in pairs]
    control_counts = [row["s16a_control_count_3pmpc"] for row in pairs]
    ratios = [row["control_to_target_count_ratio"] for row in pairs]
    differences = [row["target_minus_control_count"] for row in pairs]
    summary = {
        "schema": "astrafek-s23b-s16a-icl-control-count-audit-v1",
        "status": "complete_local_paired_count_audit",
        "science_measurement": False,
        "network_used": False,
        "cas_query_submitted": True,
        "source_rows_downloaded": False,
        "images_read": False,
        "photoz_or_shape_values_read": False,
        "query_plan_sha256": plan["query_plan_sha256"],
        "cas_release_version": plan["cas_release_version"],
        "n_candidate_targets": plan["n_candidate_targets"],
        "n_target_control_pairs": len(pairs),
        "n_positive_target_centers": sum(value > 0 for value in target_counts),
        "n_positive_control_centers": sum(value > 0 for value in control_counts),
        "n_zero_control_centers": sum(value == 0 for value in control_counts),
        "target_count": summarize(target_counts),
        "control_count": summarize(control_counts),
        "control_to_target_count_ratio": summarize(ratios),
        "target_minus_control_count": summarize(differences),
        "sum_control_cone_shape_counts": sum(control_counts),
        "by_icl_filter": filter_summary,
        "spearman_target_control_count_rho_descriptive_only": float(
            spearmanr(target_counts, control_counts).statistic
        ),
        "p_values_used": False,
        "interpretation_limit": "These are raw pre-cut S16A catalog counts in paired 3 proper-Mpc cones, not usable background-source counts or a weak-lensing estimator; ICL controls do not replace the five-random-centers-per-lens WL design.",
        "credentials_saved": state.get("credentials_saved"),
        "inputs_sha256": {
            "plan_json": sha256(PLAN_DIR / "shape_count_plan.json"),
            "count_plan_sha256": plan["query_plan_sha256"],
            "count_summary": sha256(RESULTS / "count_summary.json"),
            "count_results": sha256(RESULTS / "count_results.tsv"),
            "candidate_control_table": sha256(CONTROL_TABLE),
            "target_count_summary": sha256(TARGET_PLAN_RESULTS / "count_summary.json"),
            "target_count_results": sha256(TARGET_PLAN_RESULTS / "count_results.tsv"),
            "selection_audit_summary": sha256(SELECTION_AUDIT / "summary.json"),
        },
    }
    output.mkdir(parents=True, exist_ok=False)
    table_path = output / "paired_target_control_counts.tsv"
    fields = [
        "target_id", "control_id", "icl_filter", "z_cl", "richness",
        "s16a_target_count_3pmpc", "s16a_control_count_3pmpc",
        "target_minus_control_count", "control_to_target_count_ratio",
        "delta_ebv_sfd_target_minus_control",
    ]
    with table_path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(pairs)
    summary["paired_table_sha256"] = sha256(table_path)
    summary_path = output / "summary.json"
    summary_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    print(json.dumps(analyze(args.output), ensure_ascii=False, indent=2))
    print(f"output={args.output.resolve()}")


if __name__ == "__main__":
    main()
