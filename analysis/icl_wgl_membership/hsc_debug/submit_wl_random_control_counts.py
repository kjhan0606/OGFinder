#!/usr/bin/env python3
"""Submit approved HSC random-point candidate-count queries.

Only one integer count per approved lens center is retrieved. The script does
not export random coordinates, source rows, shapes, photo-z values, or images.
Credentials are requested once with hidden terminal input and held in memory.
"""
from __future__ import annotations

import argparse
import csv
import getpass
import hashlib
import json
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
from submit_wl_shape_counts import load_plan as load_lens_plan

DEFAULT_PLAN = Path(
    "/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/"
    "wl_random_point_count_plan_2026-10-19_02"
)


def load_random_plan(plan_dir: Path):
    manifest_path = plan_dir / "random_control_count_plan.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("schema") != "hsc-pdr3-random-control-count-plan-v1":
        raise ValueError("Unexpected random-control plan schema")
    if manifest.get("status") != "prepared_not_submitted":
        raise ValueError("Random-control plan is not in prepared state")
    if manifest.get("network_used") is not False or manifest.get("query_submitted") is not False:
        raise ValueError("Plan provenance does not show an unsubmitted plan")
    if manifest.get("output_columns") != ["lens_id", "n_random_candidates"]:
        raise ValueError("Unexpected approved output columns")
    if manifest.get("n_lenses") != 800 or manifest.get("n_query_batches") != 32:
        raise ValueError("Expected 800 lenses split into exactly 32 batches")
    preflight = Path(manifest["source_preflight"])
    if sha(preflight / "preflight_manifest.json") != manifest["source_preflight_sha256"]:
        raise ValueError("Source preflight checksum changed")
    lens_manifest = load_lens_plan(preflight)
    lenses = lens_manifest["lenses"]
    lens_ids = [str(row["lens_id"]) for row in lenses]
    if len(lens_ids) != 800 or len(set(lens_ids)) != 800:
        raise ValueError("Source preflight does not contain 800 unique lenses")
    batches = manifest.get("batches", [])
    if len(batches) != 32:
        raise ValueError("Expected 32 random-count SQL batches")
    for index, item in enumerate(batches):
        path = Path(item["path"])
        if path.parent.resolve() != plan_dir.resolve() or path.name != f"random_counts_{index:03d}.sql":
            raise ValueError("SQL batch path is outside the approved plan directory")
        if sha(path) != item["sha256"]:
            raise ValueError(f"SQL checksum changed: {path.name}")
        sql = path.read_text(encoding="utf-8").strip()
        upper = sql.upper()
        if not upper.startswith("SELECT ") or "COUNT(*) AS N_RANDOM_CANDIDATES" not in upper:
            raise ValueError(f"Batch is not an aggregate count query: {path.name}")
        if "FROM PDR3_WIDE.RANDOM AS R" not in upper or "R.COORD" not in upper:
            raise ValueError(f"Unexpected random table or cone selection: {path.name}")
        if "R.ISPRIMARY IS TRUE" not in upper:
            raise ValueError(f"Batch lacks the documented primary-region filter: {path.name}")
        if "SELECT *" in upper or "R.RA" in upper or "R.DEC" in upper:
            raise ValueError(f"Batch may return random positions: {path.name}")
        start = index * int(manifest["batch_size"])
        expected = lens_ids[start:start + int(manifest["batch_size"])]
        if len(expected) != int(item["n_lenses"]):
            raise ValueError(f"Lens count mismatch in {path.name}")
        for lens_id in expected:
            if f"SELECT '{lens_id}' AS LENS_ID" not in upper:
                raise ValueError(f"Expected lens absent from {path.name}")
    return manifest, lens_ids


def read_counts(path: Path, expected_ids):
    rows = []
    with path.open(encoding="utf-8-sig", newline="") as stream:
        reader = catalogue_reader(stream)
        if reader.fieldnames != ["lens_id", "n_random_candidates"]:
            raise ValueError("CAS returned columns outside the approved two-column schema")
        for row in reader:
            lens_id = str(row["lens_id"])
            count = int(row["n_random_candidates"])
            if count < 0:
                raise ValueError("CAS returned a negative candidate count")
            rows.append({"lens_id": lens_id, "n_random_candidates": count})
    if len(rows) != len(expected_ids) or {row["lens_id"] for row in rows} != set(expected_ids):
        raise ValueError("CAS result did not return exactly one count per requested lens")
    if len({row["lens_id"] for row in rows}) != len(rows):
        raise ValueError("CAS result contains duplicate lens IDs")
    return rows


def safe_error_detail(error, credential):
    try:
        raw = error.read(4096).decode("utf-8", "replace")
    except Exception:
        return "CAS rejected the query; response detail unavailable"
    try:
        payload = json.loads(raw)
        if isinstance(payload, dict):
            detail = next((payload[key] for key in ("error", "message", "detail", "reason")
                           if isinstance(payload.get(key), (str, int, float))), None)
            if detail is None:
                return "CAS rejected the query; response keys: " + ", ".join(sorted(payload)[:12])
            raw = str(detail)
    except json.JSONDecodeError:
        pass
    raw = raw.replace(credential, "[redacted]")
    raw = raw.replace("kjhan0606", "[redacted]")
    if re.search(r"\b(SELECT|FROM|LENS_ID|RANDOM_CANDIDATES)\b", raw, re.I):
        return "CAS rejected the request; SQL-related response detail omitted"
    return re.sub(r"\s+", " ", raw).strip()[:500] or "CAS rejected the query"


def fetch_count(client, name, sql, expected_ids, results, state, checkpoint):
    item = state["jobs"].setdefault(name, {"status": "pending"})
    target = results / f"{name}.csv"
    if item.get("status") == "complete":
        if not target.is_file() or sha(target) != item.get("sha256"):
            raise ValueError(f"Saved count result missing or changed: {name}")
        read_counts(target, expected_ids)
        return
    if item.get("status") in ("submitting", "error"):
        raise ValueError(f"{name} has an ambiguous or failed submission; inspect CAS before retry")
    if "job_id" not in item:
        item["status"] = "submitting"
        save(checkpoint, state)
        try:
            submitted = client.json(
                "submit", catalog_job=dict(sql=sql, out_format="csv",
                    include_metainfo_to_body=False, release_version="pdr3"),
                nomail=True, skip_syntax_check=False,
            )
        except urllib.error.HTTPError as error:
            if error.code in (400, 401, 403, 406, 429):
                item["status"] = "pending"
                save(checkpoint, state)
            if error.code == 406:
                print("CAS 406 detail: " + safe_error_detail(error, client.credential["password"]),
                      file=sys.stderr)
            raise
        item.update(status="submitted", job_id=str(submitted["id"]))
        save(checkpoint, state)
    while True:
        status = client.json("status", id=item["job_id"])
        if status.get("status") == "done":
            break
        if status.get("status") == "error":
            item.update(status="error", query_error=str(status.get("error", "CAS query failed"))[:1000])
            save(checkpoint, state)
            raise ValueError(f"CAS query failed for {name}")
        print(f"{name}: {status.get('status', 'unknown')}", flush=True)
        time.sleep(60)
    partial = results / f"{name}.csv.part"
    if partial.exists():
        partial.rename(results / f"{name}.csv.incomplete-{time.time_ns()}")
    with client.post("download", id=item["job_id"]) as response, partial.open("xb") as stream:
        for block in iter(lambda: response.read(1024 * 1024), b""):
            stream.write(block)
    rows = read_counts(partial, expected_ids)
    partial.replace(target)
    item.update(status="complete", rows=len(rows), bytes=target.stat().st_size,
                sha256=sha(target), total_candidates=sum(r["n_random_candidates"] for r in rows),
                zero_candidate_lenses=sum(r["n_random_candidates"] == 0 for r in rows))
    save(checkpoint, state)
    print(f"{name}: verified {len(rows)} aggregate rows", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", type=Path, default=DEFAULT_PLAN)
    parser.add_argument("--results", type=Path, help="new result directory; defaults to <plan>/cas_counts")
    parser.add_argument("--user", default="kjhan0606")
    parser.add_argument("--check", action="store_true", help="validate plan only; no network")
    parser.add_argument("--submit", action="store_true", help="submit the approved count-only batches")
    args = parser.parse_args()
    if args.check == args.submit:
        parser.error("select exactly one of --check or --submit")
    plan_dir = args.plan.resolve()
    manifest, lens_ids = load_random_plan(plan_dir)
    results = (args.results or plan_dir / "cas_counts").resolve()
    print(json.dumps({"plan_status": manifest["status"], "lenses": len(lens_ids),
        "query_batches": len(manifest["batches"]), "returned_columns": manifest["output_columns"],
        "random_coordinates_returned": False, "shape_rows_returned": False,
        "results": str(results)}, ensure_ascii=False), flush=True)
    if args.check:
        return 0
    if not sys.stdin.isatty():
        parser.error("--submit requires a terminal for one hidden password prompt")
    if not platform.node().lower().startswith("lageunha"):
        parser.error("Run --submit from the established lageunha interactive environment")
    checkpoint = results / "count_status.json"
    plan_sha = sha(plan_dir / "random_control_count_plan.json")
    if results.exists():
        if not checkpoint.is_file():
            parser.error(f"Refusing to reuse a result directory without checkpoint: {results}")
        state = json.loads(checkpoint.read_text(encoding="utf-8"))
        if state.get("plan_sha256") != plan_sha:
            parser.error("Existing count checkpoint belongs to a different plan")
        if any(item.get("status") == "submitting" for item in state.get("jobs", {}).values()):
            parser.error("Submission outcome is ambiguous; inspect CAS job history before retry")
        if any(item.get("status") == "error" for item in state.get("jobs", {}).values()):
            parser.error("A CAS job failed; inspect the checkpoint before retry")
    else:
        results.mkdir(parents=True, exist_ok=False)
        state = {"schema": "hsc-pdr3-random-control-count-run-v1", "status": "pending",
            "plan_sha256": plan_sha, "jobs": {}, "credentials_saved": False,
            "random_coordinates_returned": False, "source_rows_returned": False,
            "science_measurement": False}
    import fcntl
    with (results / "run.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        password = getpass.getpass("HSC password (hidden): ")
        if not password:
            raise ValueError("Empty password")
        client = Cas(args.user, password)
        state["status"] = "running"
        save(checkpoint, state)
        try:
            size = int(manifest["batch_size"])
            for index, item in enumerate(manifest["batches"]):
                start = index * size
                ids = lens_ids[start:start + size]
                sql = Path(item["path"]).read_text(encoding="utf-8")
                state["current_query"] = Path(item["path"]).stem
                save(checkpoint, state)
                fetch_count(client, Path(item["path"]).stem, sql, ids, results, state, checkpoint)
            rows = []
            for index in range(len(manifest["batches"])):
                with (results / f"random_counts_{index:03d}.csv").open(encoding="utf-8-sig", newline="") as stream:
                    rows.extend(dict(row) for row in catalogue_reader(stream))
            if len(rows) != 800 or {str(row["lens_id"]) for row in rows} != set(lens_ids):
                raise ValueError("Combined result does not contain exactly 800 unique lenses")
            counts = [int(row["n_random_candidates"]) for row in rows]
            state.update(status="complete", current_query=None, total_count_rows=800,
                total_random_candidates=sum(counts), zero_candidate_lenses=sum(value == 0 for value in counts))
            save(checkpoint, state)
            summary = {"schema": "hsc-pdr3-random-control-count-summary-v1", "status": "complete",
                "n_lenses": 800, "n_count_rows": len(rows), "sum_candidate_counts": sum(counts),
                "zero_candidate_lenses": sum(value == 0 for value in counts),
                "random_coordinates_returned": False, "shape_rows_returned": False,
                "source_rows_returned": False, "science_measurement": False,
                "query_plan_sha256": state["plan_sha256"]}
            (results / "count_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
        except BaseException as error:
            state.update(status="interrupted" if isinstance(error, KeyboardInterrupt) else "stopped_on_error",
                         error_type=type(error).__name__)
            save(checkpoint, state)
            raise
        finally:
            client.credential.clear()
            password = None
    print(f"Complete: 800 aggregate counts; no random coordinates or source rows. Results: {results}")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except urllib.error.HTTPError as error:
        print(f"HSC CAS HTTP {error.code}; credentials and response body suppressed.", file=sys.stderr)
        sys.exit(1)
    except KeyboardInterrupt:
        print("\nStopped; submitted batches remain checkpointed.", file=sys.stderr)
        sys.exit(130)
