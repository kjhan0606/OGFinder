#!/usr/bin/env python3
"""Audit overlap between ICL-supported pairs and the existing random-subtracted WL sample.

This is a local, read-only audit. It reads already approved exports and makes no
HSC requests. The resulting eligibility table is exploratory and is not a
weak-lensing measurement.
"""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path


DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
PAIR_FILE = DATA / "icl_satellite_supported_94_association_2026-11-05_01/supported_pairs.tsv"
ICL_METRICS = DATA / "icl_annulus_support_2026-10-18_01/pair_metrics_k2p5_c40_r50_300.tsv"
TARGET_ROWS = DATA / "wl_shape_rows_export_2026-10-18_03"
TARGET_SUMMARY = TARGET_ROWS / "row_summary.json"
TARGET_STATUS = TARGET_ROWS / "row_status.json"
RANDOM_CENTERS = DATA / "wl_random_centers_lens_reuse_2026-11-03_01/random_centers.json"
RANDOM_PLAN = DATA / "wl_random_shape_count_plan_lens_reuse_2026-11-03_01/random_shape_count_plan.json"
RANDOM_COUNTS = DATA / "wl_random_shape_count_results_lens_reuse_2026-11-03_01"
RANDOM_COUNT_STATUS = RANDOM_COUNTS / "count_status.json"
RANDOM_ROW_SUMMARY = DATA / "wl_random_shape_rows_lens_reuse_2026-11-03_01/results/row_summary.json"
WL_STACK = DATA / "wl_random_shape_stack_lens_reuse_2026-11-04_01/summary.json"
OUTPUT = DATA / "icl_wl_sample_overlap_audit_2026-11-05_03"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def verify_hash(path: Path, expected: str, label: str) -> None:
    actual = sha256(path)
    if actual != expected:
        raise ValueError(f"{label} hash mismatch: {path}")


def load_supported_pairs() -> tuple[dict[str, set[str]], dict[str, dict[str, dict[str, str]]]]:
    ids: dict[str, set[str]] = {"i": set(), "i2": set()}
    rows: dict[str, dict[str, dict[str, str]]] = {"i": {}, "i2": {}}
    with PAIR_FILE.open(encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        for row in reader:
            filt, cluster = row["filter"], str(row["cluster_id"])
            ids[filt].add(cluster)
            rows[filt][cluster] = row
    return ids, rows


def load_icl_supported_pairs() -> tuple[dict[str, set[str]], dict[str, str]]:
    ids: dict[str, set[str]] = {"i": set(), "i2": set()}
    with ICL_METRICS.open(encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        for row in reader:
            ids[row["filter"]].add(str(row["id"]))
    if sum(map(len, ids.values())) != 224:
        raise ValueError("expected 224 ICL-supported cluster/filter pairs")
    return ids, {"icl_supported_pairs_sha256": sha256(ICL_METRICS)}


def verify_target_row_export() -> tuple[set[str], dict[str, int], dict[str, str]]:
    summary = read_json(TARGET_SUMMARY)
    status = read_json(TARGET_STATUS)
    if summary.get("status") != "complete" or summary.get("source_rows_downloaded") is not True:
        raise ValueError("target source-row export is not complete")
    jobs = status.get("jobs", {})
    if status.get("status") != "complete" or len(jobs) != 32:
        raise ValueError("target source-row checkpoint is not complete")
    row_counts: dict[str, int] = {}
    provenance: dict[str, str] = {}
    n_rows = 0
    for index in range(32):
        name = f"shape_rows_{index:03d}"
        job = jobs.get(name, {})
        path = TARGET_ROWS / f"{name}.csv"
        if job.get("status") != "complete" or not path.is_file():
            raise ValueError(f"target source-row batch missing: {name}")
        verify_hash(path, job.get("sha256", ""), name)
        provenance[name] = job["sha256"]
        with path.open(encoding="utf-8-sig", newline="") as stream:
            reader = csv.DictReader(stream)
            if reader.fieldnames:
                reader.fieldnames = [name.lstrip("# ").strip() for name in reader.fieldnames]
            if not reader.fieldnames or reader.fieldnames[0] != "lens_id":
                raise ValueError(f"unexpected target source-row schema in {name}")
            for row in reader:
                lens = str(row["lens_id"])
                row_counts[lens] = row_counts.get(lens, 0) + 1
                n_rows += 1
    if n_rows != summary.get("n_pair_rows"):
        raise ValueError("target source-row total differs from its summary")
    provenance["row_summary.json"] = sha256(TARGET_SUMMARY)
    provenance["row_status.json"] = sha256(TARGET_STATUS)
    return set(row_counts), row_counts, provenance


def load_random_support() -> tuple[set[str], dict[str, str], dict[str, int], dict[str, str]]:
    plan = read_json(RANDOM_PLAN)
    status = read_json(RANDOM_COUNT_STATUS)
    centers_doc = read_json(RANDOM_CENTERS)
    row_summary = read_json(RANDOM_ROW_SUMMARY)
    if (plan.get("n_random_centers") != 3100 or plan.get("n_matched_lenses") != 620
            or status.get("status") != "complete" or row_summary.get("status") != "complete"):
        raise ValueError("approved random-center count or row export is incomplete")
    verify_hash(RANDOM_CENTERS, plan["source_centers_sha256"], "random-center manifest")
    plan_hash = sha256(RANDOM_PLAN)
    if (status.get("plan_sha256") != plan_hash
            or row_summary.get("count_plan_sha256") != plan_hash
            or row_summary.get("n_random_centers") != 3100
            or row_summary.get("n_matched_lenses") != 620
            or row_summary.get("n_source_pairs") != 9010324):
        raise ValueError("random export/count checkpoints do not match the approved plans")
    center_to_lens: dict[str, str] = {}
    for center in centers_doc.get("centers", []):
        center_id = str(center["random_center_id"])
        lens_id = str(center["matched_lens_id"])
        if center_id in center_to_lens:
            raise ValueError(f"duplicate random center: {center_id}")
        center_to_lens[center_id] = lens_id
    if len(center_to_lens) != 3100:
        raise ValueError("random-center manifest does not contain exactly 3,100 centers")

    positive_random_lenses: set[str] = set()
    source_counts: dict[str, int] = {}
    count_jobs = status.get("jobs", {})
    if len(count_jobs) != 124:
        raise ValueError("random count checkpoint does not contain 124 batches")
    count_hashes: dict[str, str] = {}
    for index in range(124):
        name = f"random_shape_counts_{index:03d}"
        job = count_jobs.get(name, {})
        path = RANDOM_COUNTS / f"{name}.csv"
        if job.get("status") != "complete" or not path.is_file():
            raise ValueError(f"random count batch missing: {name}")
        verify_hash(path, job.get("sha256", ""), name)
        count_hashes[name] = job["sha256"]
        with path.open(encoding="utf-8-sig", newline="") as stream:
            reader = csv.DictReader(stream)
            if reader.fieldnames:
                reader.fieldnames = [name.lstrip("# ").strip() for name in reader.fieldnames]
            if reader.fieldnames != ["random_center_id", "matched_lens_id", "n_shape_sources"]:
                raise ValueError(f"unexpected random count schema in {name}")
            for row in reader:
                center_id, lens_id = str(row["random_center_id"]), str(row["matched_lens_id"])
                if center_to_lens.get(center_id) != lens_id:
                    raise ValueError(f"random count row mapping mismatch for {center_id}")
                count = int(row["n_shape_sources"])
                if count < 0 or center_id in source_counts:
                    raise ValueError(f"invalid or duplicate random count for {center_id}")
                source_counts[center_id] = count
                if count > 0:
                    positive_random_lenses.add(lens_id)
    if set(source_counts) != set(center_to_lens):
        raise ValueError("random count rows do not cover the exact center manifest")
    if sum(source_counts.values()) != 9010324 or sum(v == 0 for v in source_counts.values()) != 1011:
        raise ValueError("random count totals differ from the approved count summary")
    if not positive_random_lenses <= set(center_to_lens.values()):
        raise ValueError("random support includes an unrequested lens")
    hashes = {
        "random_plan_sha256": sha256(RANDOM_PLAN),
        "random_centers_sha256": sha256(RANDOM_CENTERS),
        "random_count_status_sha256": sha256(RANDOM_COUNT_STATUS),
        "random_row_summary_sha256": sha256(RANDOM_ROW_SUMMARY),
        "random_count_batch_hashes": hashlib.sha256(
            "\n".join(f"{key}:{value}" for key, value in sorted(count_hashes.items())).encode()
        ).hexdigest(),
    }
    return positive_random_lenses, center_to_lens, source_counts, hashes


def main() -> None:
    supported_by_filter, pair_rows = load_supported_pairs()
    icl_supported_by_filter, icl_hashes = load_icl_supported_pairs()
    target_lenses, target_row_counts, target_hashes = verify_target_row_export()
    random_lenses, center_to_lens, random_source_counts, random_hashes = load_random_support()
    stack = read_json(WL_STACK)
    stack_input = stack["input"]
    selected = set(map(str, stack_input["matched_lens_ids"]))
    selected_hash = hashlib.sha256("\n".join(sorted(selected)).encode()).hexdigest()
    if selected_hash != stack_input.get("matched_lens_ids_sha256"):
        raise ValueError("matched-lens IDs disagree with their saved stack hash")
    stack_centers = stack_input.get("random_center_matched_lens_ids", {})
    stack_counts = {key: int(value) for key, value in
                    stack_input.get("random_center_expected_source_counts", {}).items()}
    if stack_centers != center_to_lens or stack_counts != random_source_counts:
        raise ValueError("random mapping/count files disagree with the stored WL stack inputs")
    common = selected & target_lenses & random_lenses
    if len(common) != stack_input["matched_sample"]["n_common_shape_support_lenses"]:
        raise ValueError("reconstructed common WL support disagrees with the stored stack summary")
    common_hash = hashlib.sha256("\n".join(sorted(common)).encode()).hexdigest()
    if common_hash != stack_input["matched_sample"]["common_lens_ids_sha256"]:
        raise ValueError("reconstructed common WL support hash disagrees with the stored stack")

    reports = {}
    icl_reports = {}
    pair_output = []
    for filt, ids in supported_by_filter.items():
        in_selected = ids & selected
        in_target = ids & target_lenses
        in_random = ids & random_lenses
        in_common = ids & common
        reports[filt] = {
            "n_icl_satellite_supported_clusters": len(ids),
            "n_in_wl_matched_620": len(in_selected),
            "n_with_target_shape_rows": len(in_target),
            "n_with_positive_random_shape_counts": len(in_random),
            "n_in_verified_random_subtracted_common_support": len(in_common),
            "n_without_common_wl_support": len(ids - common),
            "cluster_ids_in_common_support": sorted(in_common, key=int),
        }
        icl_ids = icl_supported_by_filter[filt]
        icl_reports[filt] = {
            "n_icl_profile_supported_clusters": len(icl_ids),
            "n_in_wl_matched_620": len(icl_ids & selected),
            "n_with_target_shape_rows": len(icl_ids & target_lenses),
            "n_with_positive_random_shape_counts": len(icl_ids & random_lenses),
            "n_in_verified_random_subtracted_common_support": len(icl_ids & common),
            "n_with_satellite_profile_support": len(ids),
            "n_satellite_supported_in_common_wl_support": len(ids & common),
            "icl_supported_ids_in_common_wl_support": sorted(icl_ids & common, key=int),
        }
        for cluster in sorted(ids, key=int):
            pair = pair_rows[filt][cluster]
            pair_output.append({
                "cluster_id": cluster,
                "filter": filt,
                "wl_matched_sample": cluster in selected,
                "target_shape_rows": cluster in target_lenses,
                "positive_random_shape_count": cluster in random_lenses,
                "random_subtracted_wl_common_support": cluster in common,
                "target_pair_rows": target_row_counts.get(cluster, 0),
                "icl_delta_200_300": pair["icl_delta_200_300"],
                "satellite_delta_50_75": pair["satellite_delta_50_75"],
                "redshift": pair["redshift"],
                "richness": pair["richness"],
            })

    if OUTPUT.exists():
        raise FileExistsError(f"refusing to overwrite audit output: {OUTPUT}")
    OUTPUT.mkdir(parents=True)
    pair_path = OUTPUT / "supported_pair_wl_eligibility.tsv"
    fields = list(pair_output[0])
    with pair_path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(pair_output)
    report = {
        "schema": "astrafek-icl-wl-sample-overlap-audit-v1",
        "status": "complete_exploratory",
        "science_measurement": False,
        "network_used": False,
        "hsc_query_submitted": False,
        "scope": "exact overlap of the existing ICL-satellite supported pairs and approved random-subtracted WL common lens sample",
        "n_icl_profile_supported_pair_rows": sum(map(len, icl_supported_by_filter.values())),
        "n_icl_satellite_pair_rows": sum(map(len, supported_by_filter.values())),
        "wl_sample": {
            "n_matched_lenses": len(selected),
            "n_common_random_subtracted_lenses": len(common),
            "common_lens_ids_sha256": common_hash,
        },
        "icl_profile_support_by_filter": icl_reports,
        "joint_icl_satellite_support_by_filter": reports,
        "inputs": {
            "supported_pairs_sha256": sha256(PAIR_FILE),
            **icl_hashes,
            "target_row_hashes_sha256": hashlib.sha256(
                "\n".join(f"{key}:{value}" for key, value in sorted(target_hashes.items())).encode()
            ).hexdigest(),
            **random_hashes,
            "wl_stack_summary_sha256": sha256(WL_STACK),
            "matched_lens_ids_sha256": selected_hash,
        },
        "eligibility_table": str(pair_path),
        "eligibility_table_sha256": sha256(pair_path),
        "limitations": [
            "This is a sample-overlap audit, not a lensing measurement or an ICL-satellite causal result.",
            "Existing WL stacks remain exploratory; source photo-z PDFs are not propagated.",
            "Filter-specific ICL-satellite rows share clusters and are not independent observations.",
        ],
    }
    report_path = OUTPUT / "summary.json"
    report_path.write_text(json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
