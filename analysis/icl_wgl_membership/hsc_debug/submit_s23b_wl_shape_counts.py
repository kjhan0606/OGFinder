#!/usr/bin/env python3
"""Submit an explicitly approved S23B count-only HSC S19A query plan.

The only returned values are one namespaced lens ID and one integer source
count per S23B member-profile target. Submission requires an explicit,
plan-hash-pinned approval flag and an interactive lageunha terminal.
"""
from __future__ import annotations

import argparse
import csv
import getpass
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import re
import sys
import time
import urllib.error

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from fetch_ancillary_catalogues import Cas, catalogue_reader, save, sha

DEFAULT_PLAN = Path(
    "/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/"
    "wl_s23b_shape_count_plan_2026-11-07_01"
)


def canonical_plan_hash(plan: dict) -> str:
    payload = dict(plan)
    claimed = payload.pop("query_plan_sha256", None)
    digest = hashlib.sha256(json.dumps(
        payload, sort_keys=True, separators=(",", ":"), ensure_ascii=False,
    ).encode("utf-8")).hexdigest()
    if claimed != digest:
        raise ValueError("query-plan hash does not match the saved plan contents")
    return digest


def validate_plan(plan_dir: Path) -> tuple[dict, str]:
    plan_dir = plan_dir.resolve()
    path = plan_dir / "shape_count_plan.json"
    plan = json.loads(path.read_text(encoding="utf-8"))
    if plan.get("schema") != "astrafek-s23b-s19a-wl-shape-count-plan-v1":
        raise ValueError("unexpected S23B shape-count plan schema")
    if plan.get("status") != "prepared_not_submitted":
        raise ValueError("plan is not in prepared_not_submitted state")
    if (plan.get("execution_requires_explicit_user_approval") is not True
            or plan.get("network_used") is not False
            or plan.get("query_submitted") is not False
            or plan.get("source_rows_downloaded") is not False
            or plan.get("images_returned") is not False
            or plan.get("photoz_returned") is not False
            or plan.get("science_measurement") is not False):
        raise ValueError("plan transfer or science-gate metadata is inconsistent")
    plan_hash = canonical_plan_hash(plan)

    for key, hash_key in (
        ("source_target_manifest", "source_target_manifest_sha256"),
        ("source_member_profile_inputs", "source_member_profile_inputs_sha256"),
        ("approved_member_catalog_file", "approved_member_catalog_sha256"),
    ):
        source = Path(plan[key]).resolve()
        if not source.is_file() or sha(source) != plan[hash_key]:
            raise ValueError(f"source provenance changed or is unavailable: {key}")

    centers = plan.get("centers", [])
    lens_ids = [str(row["lens_id"]) for row in centers]
    if (len(centers) != 427 or plan.get("n_lenses") != len(centers)
            or len(set(lens_ids)) != len(lens_ids)
            or any(not re.fullmatch(r"S23B_[0-9]+", value) for value in lens_ids)):
        raise ValueError("expected 427 unique namespaced S23B target centers")
    if (plan.get("output_columns") != ["lens_id", "n_shape_sources"]
            or plan.get("estimated_result_rows") != 427
            or plan.get("query_tables") != [
                "s19a_wide.meas2", "s19a_wide.weaklensing_hsm_regauss", "s19a_wide.meas",
            ]):
        raise ValueError("query output or source-table scope changed")

    expected_batches = math.ceil(len(centers) / int(plan["batch_size"]))
    batches = plan.get("batches", [])
    if len(batches) != expected_batches or plan.get("n_query_batches") != expected_batches:
        raise ValueError("batch count does not cover the approved center list")
    observed_ids = []
    for index, batch in enumerate(batches):
        name = f"shape_counts_{index:03d}.sql"
        query_path = (plan_dir / name).resolve()
        if (Path(batch["name"]).name != name or query_path.parent != plan_dir
                or not query_path.is_file() or sha(query_path) != batch.get("sha256")
                or query_path.stat().st_size != batch.get("sql_bytes")):
            raise ValueError(f"missing or modified approved SQL batch: {name}")
        sql = query_path.read_text(encoding="utf-8")
        upper = sql.upper()
        if (not upper.strip().startswith("SELECT ")
                or upper.count("COUNT(DISTINCT B.OBJECT_ID)") != int(batch["n_lenses"])
                or "SELECT *" in upper
                or any(field in upper for field in (
                    "PHOTOZ", "SOURCE_RA", "SOURCE_DEC", "SHAPE_E1", "PSF", "VARIANCE", "IMAGE",
                ))
                or upper.count("S19A_WIDE.MEAS2") != int(batch["n_lenses"])
                or upper.count("S19A_WIDE.WEAKLENSING_HSM_REGAUSS") != int(batch["n_lenses"])
                or upper.count("S19A_WIDE.MEAS C") != int(batch["n_lenses"])):
            raise ValueError(f"SQL is not the approved minimal count query: {name}")
        batch_ids = re.findall(r"'((?:S23B_[0-9]+))'\s+AS\s+LENS_ID", sql, re.IGNORECASE)
        start = index * int(plan["batch_size"])
        expected = lens_ids[start:start + int(plan["batch_size"])]
        if batch_ids != expected or len(batch_ids) != int(batch["n_lenses"]):
            raise ValueError(f"SQL lens ID list differs from the approved plan: {name}")
        observed_ids.extend(batch_ids)
    if observed_ids != lens_ids:
        raise ValueError("SQL batches do not cover the approved lens list exactly once")
    return plan, plan_hash


def read_count_csv(path: Path, expected_ids: list[str]) -> list[dict]:
    rows = []
    with path.open(encoding="utf-8-sig", newline="") as stream:
        reader = catalogue_reader(stream)
        if reader.fieldnames != ["lens_id", "n_shape_sources"]:
            raise ValueError("CAS count result has unexpected columns")
        for row in reader:
            if None in row or any(value is None for value in row.values()):
                raise ValueError("CAS count result has a malformed row")
            lens_id = str(row["lens_id"])
            count = int(row["n_shape_sources"])
            if count < 0:
                raise ValueError("CAS returned a negative source count")
            rows.append({"lens_id": lens_id, "n_shape_sources": count})
    if ({row["lens_id"] for row in rows} != set(expected_ids)
            or len({row["lens_id"] for row in rows}) != len(expected_ids)):
        raise ValueError("CAS count result IDs differ from the approved batch")
    return rows


def fetch_count(client, name, sql, expected_ids, results, state, checkpoint,
                sleep=time.sleep) -> None:
    item = state["jobs"].setdefault(name, {"status": "pending"})
    target = results / f"{name}.csv"
    if item.get("status") == "complete":
        if not target.is_file() or sha(target) != item.get("sha256"):
            raise ValueError(f"saved count result is missing or changed: {name}")
        read_count_csv(target, expected_ids)
        return
    if item.get("status") == "submitting" and "job_id" not in item:
        raise ValueError(f"CAS submission outcome is unknown for {name}; inspect CAS before retry")
    if item.get("status") == "error":
        raise ValueError(f"CAS query failed for {name}; inspect the checkpoint before retry")

    if "job_id" not in item:
        item["status"] = "submitting"
        save(checkpoint, state)
        try:
            submitted = client.json("submit", catalog_job=dict(
                sql=sql, out_format="csv", include_metainfo_to_body=False,
                release_version="pdr3"), nomail=True, skip_syntax_check=False)
        except urllib.error.HTTPError as error:
            if error.code in (400, 401, 403, 406, 429):
                item["status"] = "pending"
                save(checkpoint, state)
            raise
        item.update(status="submitted", job_id=str(submitted["id"]))
        save(checkpoint, state)

    while True:
        result = client.json("status", id=item["job_id"])
        if result.get("status") == "done":
            break
        if result.get("status") == "error":
            message = str(result.get("error", "CAS query failed"))
            for secret in client.credential.values():
                if secret:
                    message = message.replace(secret, "[redacted]")
            item.update(status="error", query_error=message[:2000])
            save(checkpoint, state)
            raise ValueError(f"{name}: {message[:300]}")
        print(f"{name}: {result.get('status', 'unknown')}", flush=True)
        sleep(60)

    partial = results / f"{name}.csv.part"
    if partial.exists():
        partial.rename(results / f"{name}.csv.incomplete-{time.time_ns()}")
    with client.post("download", id=item["job_id"]) as response, partial.open("xb") as stream:
        for block in iter(lambda: response.read(1024 * 1024), b""):
            stream.write(block)
    rows = read_count_csv(partial, expected_ids)
    partial.replace(target)
    item.update(status="complete", rows=len(rows), bytes=target.stat().st_size,
                sha256=sha(target), total_sources=sum(row["n_shape_sources"] for row in rows),
                zero_source_lenses=sum(row["n_shape_sources"] == 0 for row in rows))
    save(checkpoint, state)
    print(f"{name}: verified {len(rows)} count rows", flush=True)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan-dir", type=Path, default=DEFAULT_PLAN)
    parser.add_argument("--results", type=Path,
                        help="new result directory; defaults to <plan-dir>/results")
    parser.add_argument("--user", default="kjhan0606")
    parser.add_argument("--plan", action="store_true", help="verify the plan; never contact CAS")
    parser.add_argument("--submit", action="store_true", help="submit and retrieve aggregate counts")
    parser.add_argument("--approved-plan-sha256",
                        help="required with --submit; must equal the explicitly approved plan hash")
    args = parser.parse_args()
    if args.plan == args.submit:
        parser.error("select exactly one of --plan or --submit")
    plan_dir = args.plan_dir.resolve()
    plan, plan_hash = validate_plan(plan_dir)
    results = (args.results or plan_dir / "results").resolve()
    print(json.dumps({
        "plan_status": plan["status"], "query_plan_sha256": plan_hash,
        "n_lenses": plan["n_lenses"], "batches": plan["n_query_batches"],
        "result_rows": plan["estimated_result_rows"],
        "result_scope": plan["estimated_result_scope"],
        "network_used": False if args.plan else "pending explicit submission",
        "results": str(results),
    }, ensure_ascii=False), flush=True)
    if args.plan:
        return 0

    if args.approved_plan_sha256 != plan_hash:
        parser.error("--submit requires the exact SHA-256 of the plan explicitly approved by the user")
    if not sys.stdin.isatty():
        parser.error("--submit requires an interactive terminal for one hidden password prompt")
    if not platform.node().lower().startswith("lageunha"):
        parser.error("Run --submit from the established lageunha interactive environment")

    if results.exists():
        checkpoint = results / "count_status.json"
        if not checkpoint.is_file():
            parser.error(f"refusing to reuse an existing result directory without checkpoint: {results}")
        state = json.loads(checkpoint.read_text(encoding="utf-8"))
        if state.get("plan_sha256") != plan_hash:
            parser.error("existing count checkpoint belongs to a different plan")
        if state.get("status") == "complete":
            summary = results / "count_summary.json"
            if not summary.is_file() or sha(summary) != state.get("summary_sha256"):
                parser.error("completed checkpoint summary is missing or changed")
            print(f"Already complete: {summary}")
            return 0
    else:
        results.mkdir(parents=True, exist_ok=False)
        checkpoint = results / "count_status.json"
        state = {
            "schema": "astrafek-s23b-s19a-wl-count-run-v1",
            "status": "pending", "plan_sha256": plan_hash,
            "plan_dir": str(plan_dir), "jobs": {},
            "credentials_saved": False, "source_rows_downloaded": False,
            "images_returned": False, "photoz_returned": False,
            "science_measurement": False,
        }
        save(checkpoint, state)

    import fcntl
    with (results / "run.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        password = getpass.getpass("HSC password (hidden; one prompt): ")
        if not password:
            raise ValueError("Empty password")
        client = Cas(args.user, password)
        state["status"] = "running"
        save(checkpoint, state)
        try:
            for batch in plan["batches"]:
                batch_path = plan_dir / batch["name"]
                ids = re.findall(
                    r"'((?:S23B_[0-9]+))'\s+AS\s+LENS_ID",
                    batch_path.read_text(encoding="utf-8"), re.IGNORECASE,
                )
                state["current_query"] = batch_path.stem
                save(checkpoint, state)
                fetch_count(client, batch_path.stem, batch_path.read_text(encoding="utf-8"),
                            ids, results, state, checkpoint)
            rows = []
            for batch in plan["batches"]:
                name = Path(batch["name"]).stem
                with (results / f"{name}.csv").open(encoding="utf-8-sig", newline="") as stream:
                    rows.extend(dict(row) for row in catalogue_reader(stream))
            expected_ids = [str(center["lens_id"]) for center in plan["centers"]]
            if ({str(row["lens_id"]) for row in rows} != set(expected_ids)
                    or len(set(str(row["lens_id"]) for row in rows)) != len(expected_ids)):
                raise ValueError("combined result does not contain exactly the approved 427 S23B targets")
            row_by_id = {str(row["lens_id"]): row for row in rows}
            rows = [row_by_id[lens_id] for lens_id in expected_ids]
            counts = [int(row["n_shape_sources"]) for row in rows]
            summary_path = results / "count_summary.json"
            summary = {
                "schema": "astrafek-s23b-s19a-wl-count-summary-v1",
                "status": "complete", "n_lenses": len(rows),
                "n_count_rows": len(rows), "sum_lens_cone_shape_counts": sum(counts),
                "zero_source_lenses": sum(value == 0 for value in counts),
                "source_rows_downloaded": False, "images_returned": False,
                "photoz_returned": False, "science_measurement": False,
                "query_plan_sha256": plan_hash,
            }
            with summary_path.open("x", encoding="utf-8") as stream:
                json.dump(summary, stream, ensure_ascii=False, indent=2, allow_nan=False)
                stream.write("\n")
            state.update(status="complete", current_query=None,
                         total_count_rows=len(rows), raw_shape_rows_sum=sum(counts),
                         zero_source_lenses=sum(value == 0 for value in counts),
                         summary_sha256=sha(summary_path), source_rows_downloaded=False,
                         images_returned=False, photoz_returned=False,
                         science_measurement=False)
            save(checkpoint, state)
        except BaseException as error:
            state.update(status="interrupted" if isinstance(error, KeyboardInterrupt) else "stopped_on_error",
                         error_type=type(error).__name__)
            save(checkpoint, state)
            raise
        finally:
            client.credential.clear()
            password = None
    print(f"Complete: {len(rows)} aggregate count rows; no source rows downloaded. Results: {results}")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except urllib.error.HTTPError as error:
        print(f"HSC CAS HTTP {error.code}; no credentials or response body printed.", file=sys.stderr)
        sys.exit(1)
    except KeyboardInterrupt:
        print("\nStopped; submitted jobs and completed batches remain resumable.", file=sys.stderr)
        sys.exit(130)
    except Exception as error:
        message = str(error) if isinstance(error, ValueError) else type(error).__name__
        print("submit_s23b_wl_shape_counts: " + message, file=sys.stderr)
        sys.exit(1)
