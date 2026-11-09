#!/usr/bin/env python3
"""Validate or submit the exact frame-level HSC calibration metadata plan.

This exports only frame membership and per-CCD-frame calibration metadata for
the existing target/control coadd centers. Submission requires separate user
approval of the pinned plan SHA-256; it never requests images, object rows,
shapes, or photo-z values.
"""
from __future__ import annotations

import argparse
import csv
import fcntl
import getpass
import hashlib
import json
import math
from pathlib import Path
import platform
import sys
import time
import urllib.error

from fetch_ancillary_catalogues import Cas, catalogue_reader, save
from prepare_icl_frame_metadata_rows import (
    COUNT_AUDIT,
    COUNT_PLAN,
    COUNT_RESULTS,
    DEFAULT_OUTPUT as DEFAULT_PLAN,
    OUTPUT_COLUMNS,
    load_audited_centers,
)


DEFAULT_USER = "kjhan0606"
DEFAULT_RESULTS = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_frame_metadata_row_results_2026-11-06_01")
EXPECTED_PLAN_SHA256 = "49aa4b7601dd42243209d7b3b394a4dd3b4cbf00d87d28cbb9467553c29f4ce3"
EXPECTED_SQL_AGGREGATE_SHA256 = "4270f3f0d290ac427b24694f7561f46cdee71cc8c3621a57c3fcf3f034583c76"
MAX_BATCH_BYTES = 25 * 1024 * 1024
STATUS_POLL_INTERVAL_S = 20


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_number(value: str, field: str, *, positive: bool = False,
                 nonnegative: bool = False, nullable: bool = False) -> float | None:
    value = str(value).strip()
    if nullable and value.lower() in ("", "null", "none", "nan", r"\n", r"\N"):
        return None
    try:
        number = float(value)
    except (TypeError, ValueError):
        raise ValueError(f"CAS returned a non-numeric {field}") from None
    if not math.isfinite(number):
        raise ValueError(f"CAS returned non-finite {field}")
    if positive and number <= 0:
        raise ValueError(f"CAS returned non-positive {field}")
    if nonnegative and number < 0:
        raise ValueError(f"CAS returned negative {field}")
    return number


def validate_plan(plan_dir: Path) -> tuple[dict, list[dict], list[dict], str]:
    plan_dir = Path(plan_dir).resolve()
    plan_path = plan_dir / "plan.json"
    plan = json.loads(plan_path.read_text(encoding="utf-8"))
    plan_hash = sha256(plan_path)
    if plan_hash != EXPECTED_PLAN_SHA256:
        raise ValueError("frame-metadata plan differs from the pinned, reviewed plan")
    if (plan.get("schema") != "astrafek-icl-frame-metadata-row-plan-v1"
            or plan.get("status") != "prepared_not_submitted"):
        raise ValueError("unexpected or previously submitted frame-metadata plan")
    if plan.get("tables") != ["pdr3_wide.mosaic", "pdr3_wide.mosaicframe", "pdr3_wide.frame"]:
        raise ValueError("plan table set changed")
    if plan.get("output_columns") != OUTPUT_COLUMNS or len(OUTPUT_COLUMNS) != 21:
        raise ValueError("frame-metadata export schema changed")
    if (plan.get("n_target_centers") != 800 or plan.get("n_control_centers") != 798
            or plan.get("n_centers") != 1598 or plan.get("n_query_batches") != 16
            or plan.get("estimated_rows") != 67071):
        raise ValueError("frame-metadata sample or row estimate changed")
    if (plan.get("execution_requires_explicit_user_approval") is not True
            or plan.get("images_returned") is not False
            or plan.get("shape_rows_returned") is not False
            or plan.get("photoz_rows_returned") is not False
            or plan.get("network_used") is not False
            or plan.get("query_submitted") is not False
            or plan.get("credentials_saved") is not False
            or plan.get("science_measurement") is not False):
        raise ValueError("frame-metadata plan scope/provenance flags changed")
    for path_text, expected in plan.get("input_sha256", {}).items():
        path = Path(path_text).resolve()
        if not path.is_file() or sha256(path) != expected:
            raise ValueError(f"frame-metadata plan input missing or changed: {path}")

    centers, count_audit, _ = load_audited_centers()
    if (count_audit.get("status") != "complete_exploratory"
            or count_audit.get("n_centers") != 1598
            or count_audit.get("centers_without_patch") != 0):
        raise ValueError("parent count audit is not complete")
    batches = plan.get("batches", [])
    if len(batches) != 16:
        raise ValueError("expected exactly 16 pinned SQL batches")
    aggregate = hashlib.sha256("\n".join(row["sha256"] for row in batches).encode()).hexdigest()
    if aggregate != EXPECTED_SQL_AGGREGATE_SHA256 or aggregate != plan.get("sql_aggregate_sha256"):
        raise ValueError("SQL aggregate checksum differs from the reviewed plan")
    total_rows = 0
    for index, batch in enumerate(batches):
        name = f"frame_metadata_rows_{index:03d}"
        path = (plan_dir / f"{name}.sql").resolve()
        if (batch.get("name") != name or Path(batch.get("path", "")).resolve() != path
                or not path.is_file() or path.stat().st_size != batch.get("sql_bytes")
                or sha256(path) != batch.get("sha256")):
            raise ValueError(f"SQL batch changed: {name}")
        sql = path.read_text(encoding="utf-8").lower()
        if (not sql.lstrip().startswith("-- approved-scope candidate")
                or "pdr3_wide.frame" not in sql
                or "gain1" not in sql or "exptime" not in sql
                or any(token in sql for token in ("shape", "photoz", "select *"))):
            raise ValueError(f"SQL batch has unexpected table or column scope: {name}")
        total_rows += int(batch["expected_rows_from_approved_count_query"])
    if total_rows != 67071:
        raise ValueError("batch row totals do not reproduce the approved count query")
    return plan, centers, batches, plan_hash


def read_metadata_rows(path: Path, expected_centers: list[dict]) -> tuple[dict[str, int], dict[str, int]]:
    expected = {str(row["center_id"]): row for row in expected_centers}
    seen: dict[str, set[str]] = {key: set() for key in expected}
    rows_by_center = {key: 0 for key in expected}
    nullable_fields = ("exptime", "gain1", "gain2", "gain3", "gain4",
                       "fluxmag0", "sigma_sky", "mjd")
    null_counts = {name: 0 for name in nullable_fields}
    with Path(path).open(encoding="utf-8-sig", newline="") as stream:
        reader = catalogue_reader(stream)
        if reader.fieldnames != OUTPUT_COLUMNS:
            raise ValueError("CAS returned columns outside the approved 21-column schema")
        for row in reader:
            if None in row or any(value is None for value in row.values()):
                raise ValueError("malformed CAS metadata row")
            center_id = str(row["center_id"])
            parent = expected.get(center_id)
            if parent is None:
                raise ValueError("CAS returned an unrequested coadd center")
            if (row["sample_role"] != parent["role"]
                    or row["coadd_band"] != parent["band"]
                    or row["skymap_id"] != str(parent["skymap_id"])
                    or row["tract"] != str(parent["tract"])
                    or row["patch"] != str(parent["patch"])
                    or row["patch_s"] != str(parent["patch_s"])
                    or abs(float(row["ra_deg"]) - float(parent["ra_deg"])) > 1e-8
                    or abs(float(row["dec_deg"]) - float(parent["dec_deg"])) > 1e-8):
                raise ValueError("CAS center identity or selected mosaic patch differs from audited mapping")
            frame_num = str(row["frame_num"]).strip()
            if not frame_num or frame_num in seen[center_id]:
                raise ValueError("CAS returned an empty or duplicate frame for one center")
            seen[center_id].add(frame_num)
            for field in ("visit", "ccd"):
                value = str(row[field]).strip()
                if not value.isdigit():
                    raise ValueError(f"CAS returned a non-integer {field}")
            # CAS returns the original detector-filter codes here (i/i2),
            # while the coadd selection column uses the HSC-I release label.
            if row["filter"] not in {"i", "i2"}:
                raise ValueError(f"CAS returned an unexpected contributing filter: {row['filter']}")
            for field, constraints in {
                    "exptime": {"positive": True},
                    "gain1": {"positive": True}, "gain2": {"positive": True},
                    "gain3": {"positive": True}, "gain4": {"positive": True},
                    "fluxmag0": {"positive": True},
                    "sigma_sky": {"nonnegative": True}, "mjd": {},
            }.items():
                value = parse_number(row[field], field, nullable=True, **constraints)
                if value is None:
                    null_counts[field] += 1
            rows_by_center[center_id] += 1
    for center_id, row_count in rows_by_center.items():
        expected_count = int(expected[center_id]["n_patch_frame_memberships"])
        if row_count != expected_count:
            raise ValueError(f"frame-row count differs from audited membership count for center {center_id}")
    return rows_by_center, null_counts


def cas_status(client: Cas, job_id: str) -> dict:
    for attempt in range(10):
        try:
            return client.json("status", id=job_id)
        except urllib.error.HTTPError as error:
            if error.code not in (502, 503, 504) or attempt == 9:
                raise
            print(f"CAS status gateway HTTP {error.code}; retrying in 60 seconds ({attempt + 1}/10)", flush=True)
            time.sleep(60)


def run(plan_dir: Path, results_dir: Path, user: str, *, submit: bool,
        approved_plan_sha256: str | None) -> int:
    plan_dir, results_dir = Path(plan_dir).resolve(), Path(results_dir).resolve()
    plan, centers, batches, plan_hash = validate_plan(plan_dir)
    print(json.dumps({
        "plan_status": plan["status"], "plan_sha256": plan_hash,
        "network_used": False, "centers": len(centers), "batches": len(batches),
        "estimated_rows": plan["estimated_rows"],
        "estimated_csv_size": plan["estimated_csv_size"],
        "output_columns": OUTPUT_COLUMNS, "images": False, "shape_rows": False,
        "photoz": False,
    }, ensure_ascii=False), flush=True)
    if not submit:
        if approved_plan_sha256:
            raise ValueError("approval hash is accepted only with --submit")
        return 0
    if str(approved_plan_sha256 or "").lower() != plan_hash:
        raise ValueError("submission requires explicit approval of this exact plan SHA-256")
    if not sys.stdin.isatty():
        raise ValueError("submission requires an interactive hidden credential prompt")
    if not platform.node().lower().startswith("lageunha"):
        raise ValueError("run an approved HSC metadata query from lageunha")
    if results_dir.exists():
        checkpoint = results_dir / "status.json"
        if not checkpoint.is_file():
            raise FileExistsError("refusing existing result directory without a matching checkpoint")
        state = json.loads(checkpoint.read_text(encoding="utf-8"))
        if state.get("plan_sha256") != plan_hash:
            raise ValueError("existing checkpoint belongs to a different metadata plan")
    else:
        results_dir.mkdir(parents=True, exist_ok=False)
        checkpoint = results_dir / "status.json"
        state = {
            "schema": "astrafek-icl-frame-metadata-row-run-v1", "status": "pending",
            "plan_sha256": plan_hash, "jobs": {}, "network_used": False,
            "query_submitted": False, "credentials_saved": False,
            "images_returned": False, "shape_rows_returned": False,
            "photoz_rows_returned": False, "science_measurement": False,
        }
        save(checkpoint, state)
    if any(state.get(key) is not False for key in (
            "credentials_saved", "images_returned", "shape_rows_returned",
            "photoz_rows_returned", "science_measurement")):
        raise ValueError("checkpoint exceeds the approved frame-metadata scope")
    if any(item.get("status") in ("submitting", "error") and "job_id" not in item
           for item in state.get("jobs", {}).values()):
        raise ValueError("submission outcome is ambiguous; inspect CAS before retry")

    with (results_dir / "run.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        password = getpass.getpass("HSC password (hidden; one prompt): ")
        if not password:
            raise ValueError("empty password")
        client = Cas(user, password)
        state.pop("error_type", None)
        state.update(status="running", network_used=True)
        save(checkpoint, state)
        try:
            for index, batch in enumerate(batches):
                name = batch["name"]
                item = state["jobs"].setdefault(name, {"status": "pending"})
                target = results_dir / f"{name}.csv"
                start = index * int(plan["batch_size_centers"])
                expected_centers = centers[start:start + int(plan["batch_size_centers"])]
                if item.get("status") == "complete":
                    if not target.is_file() or sha256(target) != item.get("sha256"):
                        raise ValueError(f"saved metadata rows missing or changed: {name}")
                    read_metadata_rows(target, expected_centers)
                    continue
                if item.get("status") == "error":
                    raise ValueError(f"CAS job failed for {name}; inspect before retry")
                if "job_id" not in item:
                    item["status"] = "submitting"
                    save(checkpoint, state)
                    try:
                        submitted = client.json(
                            "submit", catalog_job=dict(
                                sql=(plan_dir / f"{name}.sql").read_text(encoding="utf-8"),
                                out_format="csv", include_metainfo_to_body=False,
                                release_version="pdr3"),
                            nomail=True, skip_syntax_check=False)
                    except urllib.error.HTTPError:
                        item["status"] = "pending"
                        save(checkpoint, state)
                        raise
                    item.update(status="submitted", job_id=str(submitted["id"]))
                    state["query_submitted"] = True
                    save(checkpoint, state)
                while True:
                    status = cas_status(client, item["job_id"])
                    if status.get("status") == "done":
                        break
                    if status.get("status") == "error":
                        item.update(status="error", query_error=str(status.get("error", "CAS query failed"))[:500])
                        save(checkpoint, state)
                        raise ValueError(f"CAS frame-metadata query failed for {name}")
                    print(f"{name}: {status.get('status', 'unknown')}", flush=True)
                    time.sleep(STATUS_POLL_INTERVAL_S)
                partial = results_dir / f"{name}.csv.part"
                if partial.exists():
                    partial.rename(results_dir / f"{name}.csv.incomplete-{time.time_ns()}")
                with client.post("download", id=item["job_id"]) as response, partial.open("xb") as stream:
                    total_bytes = 0
                    for block in iter(lambda: response.read(1024 * 1024), b""):
                        total_bytes += len(block)
                        if total_bytes > MAX_BATCH_BYTES:
                            raise ValueError(f"CAS response exceeded {MAX_BATCH_BYTES} bytes for {name}")
                        stream.write(block)
                counts, nulls = read_metadata_rows(partial, expected_centers)
                partial.replace(target)
                item.update(status="complete", rows=sum(counts.values()),
                            centers=len(counts), bytes=target.stat().st_size,
                            sha256=sha256(target), metadata_null_counts=nulls)
                save(checkpoint, state)
                print(f"{name}: verified {item['rows']} frame rows", flush=True)

            all_center_counts = {}
            all_metadata_nulls = {name: 0 for name in (
                "exptime", "gain1", "gain2", "gain3", "gain4",
                "fluxmag0", "sigma_sky", "mjd")}
            for index, batch in enumerate(batches):
                target = results_dir / f"{batch['name']}.csv"
                start = index * int(plan["batch_size_centers"])
                expected_centers = centers[start:start + int(plan["batch_size_centers"])]
                counts, nulls = read_metadata_rows(target, expected_centers)
                for center_id, count in counts.items():
                    if center_id in all_center_counts:
                        raise ValueError(f"duplicate center across batches: {center_id}")
                    all_center_counts[center_id] = count
                for name, count in nulls.items():
                    all_metadata_nulls[name] += count
            if set(all_center_counts) != {str(row["center_id"]) for row in centers}:
                raise ValueError("completed frame rows do not cover all 1,598 approved centers")
            summary = {
                "schema": "astrafek-icl-frame-metadata-row-summary-v1",
                "status": "complete", "plan_sha256": plan_hash,
                "n_centers": len(all_center_counts), "n_target_centers": 800,
                "n_control_centers": 798, "n_frame_rows": sum(all_center_counts.values()),
                "null_counts_by_metadata_column": all_metadata_nulls,
                "images_returned": False, "shape_rows_returned": False,
                "photoz_rows_returned": False, "credentials_saved": False,
                "science_measurement": False,
                "scientific_limit": plan["scientific_limit"],
            }
            if summary["n_frame_rows"] != int(plan["estimated_rows"]):
                raise ValueError("returned frame-row total differs from audited parent counts")
            state.update(status="complete", summary=summary)
            save(checkpoint, state)
            (results_dir / "summary.json").write_text(
                json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
        except BaseException as error:
            state.update(status="interrupted" if isinstance(error, KeyboardInterrupt) else "stopped_on_error",
                         error_type=type(error).__name__)
            save(checkpoint, state)
            raise
        finally:
            client.credential.clear()
            password = None
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan-dir", type=Path, default=DEFAULT_PLAN)
    parser.add_argument("--results", type=Path, default=DEFAULT_RESULTS)
    parser.add_argument("--user", default=DEFAULT_USER)
    parser.add_argument("--check", action="store_true", help="validate the pinned plan and SQL without network")
    parser.add_argument("--submit", action="store_true", help="submit the exact approved plan")
    parser.add_argument("--approved-plan-sha256", help="required after explicit approval of this exact plan")
    args = parser.parse_args()
    if args.check == args.submit:
        parser.error("select exactly one of --check or --submit")
    if args.check and args.approved_plan_sha256:
        parser.error("approval hash cannot be used with --check")
    return run(args.plan_dir, args.results, args.user, submit=args.submit,
               approved_plan_sha256=args.approved_plan_sha256)


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except urllib.error.HTTPError as error:
        print(f"HSC CAS HTTP {error.code}; credentials and response body suppressed.", file=sys.stderr)
        raise SystemExit(1)
    except KeyboardInterrupt:
        print("\nStopped; submitted query IDs remain checkpointed.", file=sys.stderr)
        raise SystemExit(130)
