#!/usr/bin/env python3
"""Submit the prepared 800-lens S19A count-only preflight to HSC CAS.

Only one aggregate count row per lens is retrieved. No source shapes, positions,
photo-z values, or image data are returned by these queries. Credentials are
prompted once and kept in memory; job IDs and count results are checkpointed.
"""
import argparse
import csv
import getpass
import hashlib
import json
import os
from pathlib import Path
import platform
import sys
import time
import urllib.error

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from fetch_ancillary_catalogues import Cas, catalogue_reader, save, sha

DEFAULT_PREFLIGHT = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_shape_preflight_2026-10-18_01")


def fingerprint(manifest):
    payload = {
        "source_sha256": manifest["source_sha256"],
        "sql_files": [(item["path"], item["sha256"]) for item in manifest["sql_files"]],
        "release": "pdr3",
    }
    return hashlib.sha256(json.dumps(payload, sort_keys=True).encode()).hexdigest()


def load_plan(preflight):
    manifest_path = preflight / "preflight_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("schema") != "hsc-s19a-shape-count-preflight-v1":
        raise ValueError("Unexpected S19A preflight schema")
    if manifest.get("status") != "prepared_not_submitted" or manifest.get("n_lenses") != 800:
        raise ValueError("Expected the prepared 800-lens count-only query plan")
    if manifest.get("data_transfer_started") is not False:
        raise ValueError("Preflight manifest says data transfer already started")
    if sha(manifest["source_manifest"]) != manifest["source_sha256"]:
        raise ValueError("Lens manifest checksum changed")
    lenses = manifest.get("lenses", [])
    if len(lenses) != 800 or len({str(row["lens_id"]) for row in lenses}) != 800:
        raise ValueError("Lens list must contain 800 unique rows")
    expected_ids = [str(row["lens_id"]) for row in lenses]
    batch_size = int(manifest["batch_size"])
    if len(manifest.get("sql_files", [])) != 32:
        raise ValueError("Expected 32 SQL batches")
    for index, item in enumerate(manifest["sql_files"]):
        path = Path(item["path"])
        if path.parent.resolve() != preflight.resolve() or path.name != f"shape_counts_{index:03d}.sql":
            raise ValueError("SQL batch path is outside the prepared directory")
        if sha(path) != item["sha256"]:
            raise ValueError(f"SQL checksum changed: {path.name}")
        sql = path.read_text(encoding="utf-8").strip()
        upper = sql.upper()
        if not upper.startswith("SELECT ") or "COUNT(DISTINCT B.OBJECT_ID)" not in upper:
            raise ValueError(f"Batch is not the expected count-only query: {path.name}")
        if " B.*" in upper or "I_HSMSHAPEREGAUSS_E1" in upper or "PHOTOZ" in upper:
            raise ValueError(f"Batch would return source data rather than counts: {path.name}")
        start = index * batch_size
        batch_ids = expected_ids[start:start + batch_size]
        if len(batch_ids) != int(item["n_lenses"]):
            raise ValueError(f"Lens count mismatch in {path.name}")
        for lens_id in batch_ids:
            marker = f"SELECT '{lens_id}' AS lens_id".upper()
            if marker not in upper:
                raise ValueError(f"Lens ID is absent from its SQL batch: {path.name}")
    return manifest


def read_count_csv(path, expected_ids):
    rows = []
    with path.open(encoding="utf-8-sig", newline="") as stream:
        reader = catalogue_reader(stream)
        required = {"lens_id", "n_shape_sources"}
        if not required.issubset(reader.fieldnames or []):
            raise ValueError("Unexpected CAS count CSV columns")
        for row in reader:
            lens_id = str(row["lens_id"])
            count = int(row["n_shape_sources"])
            if count < 0:
                raise ValueError("CAS returned a negative source count")
            rows.append({"lens_id": lens_id, "n_shape_sources": count})
    if len(rows) != len(expected_ids) or {row["lens_id"] for row in rows} != set(expected_ids):
        raise ValueError("CAS result did not return exactly one count per requested lens")
    if len({row["lens_id"] for row in rows}) != len(rows):
        raise ValueError("CAS result contains duplicate lens IDs")
    return rows


def fetch_count(client, name, sql, expected_ids, result_dir, state, checkpoint, sleep=time.sleep):
    item = state["jobs"].setdefault(name, {"status": "pending"})
    target = result_dir / f"{name}.csv"
    if item.get("status") == "complete":
        if not target.is_file() or sha(target) != item.get("sha256"):
            raise ValueError(f"Saved count result missing or changed: {name}")
        read_count_csv(target, expected_ids)
        return
    if item.get("status") == "submitting" and "job_id" not in item:
        raise ValueError(f"Submission outcome unknown for {name}; inspect CAS before retry")
    if item.get("status") == "error":
        raise ValueError(f"CAS query failed for {name}; inspect checkpoint before retry")
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

    partial = result_dir / f"{name}.csv.part"
    if partial.exists():
        partial.rename(result_dir / f"{name}.csv.incomplete-{time.time_ns()}")
    with client.post("download", id=item["job_id"]) as response, partial.open("xb") as stream:
        for block in iter(lambda: response.read(1024 * 1024), b""):
            stream.write(block)
    records = read_count_csv(partial, expected_ids)
    partial.replace(target)
    item.update(status="complete", rows=len(records), bytes=target.stat().st_size,
                sha256=sha(target), total_sources=sum(row["n_shape_sources"] for row in records),
                zero_source_lenses=sum(row["n_shape_sources"] == 0 for row in records))
    save(checkpoint, state)
    print(f"{name}: verified {len(records)} count rows", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--preflight", type=Path, default=DEFAULT_PREFLIGHT)
    parser.add_argument("--results", type=Path,
                        help="new result directory; defaults to <preflight>/cas_counts")
    parser.add_argument("--user", default="kjhan0606")
    parser.add_argument("--plan", action="store_true", help="verify plan and print query count; no network")
    parser.add_argument("--submit", action="store_true", help="submit count-only queries and retrieve aggregate counts")
    args = parser.parse_args()
    if args.plan and args.submit:
        parser.error("--plan and --submit are mutually exclusive")
    preflight = args.preflight.resolve()
    manifest = load_plan(preflight)
    results_dir = (args.results or preflight / "cas_counts").resolve()
    plan_hash = fingerprint(manifest)
    print(json.dumps({"preflight_status": manifest["status"], "lenses": 800,
                      "query_batches": len(manifest["sql_files"]),
                      "returned_data": "counts only: one integer per lens; no source rows",
                      "results": str(results_dir)}, ensure_ascii=False), flush=True)
    if args.plan:
        return 0
    if not args.submit:
        parser.error("select --plan or --submit")
    if not sys.stdin.isatty():
        parser.error("--submit requires an interactive terminal for one hidden password prompt")
    if not platform.node().lower().startswith("lageunha"):
        parser.error("Run --submit from the established lageunha interactive environment")
    if results_dir.exists():
        checkpoint = results_dir / "count_status.json"
        if not checkpoint.is_file():
            parser.error(f"Refusing to reuse an existing result directory without checkpoint: {results_dir}")
        state = json.loads(checkpoint.read_text(encoding="utf-8"))
        if state.get("plan_sha256") != plan_hash:
            parser.error("Existing count checkpoint belongs to a different query plan")
    else:
        results_dir.mkdir(parents=True, exist_ok=False)
        checkpoint = results_dir / "count_status.json"
        state = dict(schema="hsc-s19a-count-run-v1", status="pending", plan_sha256=plan_hash,
                     preflight=str(preflight), jobs={}, credentials_saved=False,
                     source_rows_downloaded=False, science_measurement=False)
        save(checkpoint, state)

    import fcntl
    with (results_dir / "run.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        password = getpass.getpass("HSC password (hidden; one prompt): ")
        if not password:
            raise ValueError("Empty password")
        client = Cas(args.user, password)
        state["status"] = "running"
        save(checkpoint, state)
        try:
            batch_size = int(manifest["batch_size"])
            lenses = manifest["lenses"]
            for index, item in enumerate(manifest["sql_files"]):
                start = index * batch_size
                expected_ids = [str(row["lens_id"]) for row in lenses[start:start + batch_size]]
                path = Path(item["path"])
                state["current_query"] = path.stem
                save(checkpoint, state)
                fetch_count(client, path.stem, path.read_text(encoding="utf-8"),
                            expected_ids, results_dir, state, checkpoint)
            all_rows = []
            for index in range(len(manifest["sql_files"])):
                path = results_dir / f"shape_counts_{index:03d}.csv"
                with path.open(encoding="utf-8-sig", newline="") as stream:
                    all_rows.extend(dict(row) for row in catalogue_reader(stream))
            if len(all_rows) != 800 or len({row["lens_id"] for row in all_rows}) != 800:
                raise ValueError("Combined result does not contain exactly 800 unique lenses")
            counts = [int(row["n_shape_sources"]) for row in all_rows]
            state.update(status="complete", current_query=None, total_count_rows=800,
                         raw_shape_rows_sum=sum(counts), zero_source_lenses=sum(value == 0 for value in counts),
                         source_rows_downloaded=False, science_measurement=False)
            save(checkpoint, state)
            summary = results_dir / "count_summary.json"
            summary.write_text(json.dumps({
                "schema": "hsc-s19a-count-summary-v1", "status": "complete",
                "n_lenses": 800, "n_count_rows": len(all_rows),
                "sum_lens_cone_shape_counts": sum(counts),
                "zero_source_lenses": sum(value == 0 for value in counts),
                "source_rows_downloaded": False, "science_measurement": False,
                "query_plan_sha256": plan_hash,
            }, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        except BaseException as error:
            state.update(status="interrupted" if isinstance(error, KeyboardInterrupt) else "stopped_on_error",
                         error_type=type(error).__name__)
            save(checkpoint, state)
            raise
        finally:
            client.credential.clear()
            password = None
    print(f"Complete: 800 aggregate count rows; no source catalog rows downloaded. Results: {results_dir}")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except urllib.error.HTTPError as error:
        print(f"HSC CAS HTTP {error.code}; no credentials or response body printed.", file=sys.stderr)
        sys.exit(1)
    except KeyboardInterrupt:
        print("\nStopped; completed count batches and submitted jobs remain resumable.", file=sys.stderr)
        sys.exit(130)
    except Exception as error:
        message = str(error) if isinstance(error, ValueError) else type(error).__name__
        print("submit_wl_shape_counts: " + message, file=sys.stderr)
        sys.exit(1)
