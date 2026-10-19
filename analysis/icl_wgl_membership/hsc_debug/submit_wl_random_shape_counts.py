#!/usr/bin/env python3
"""Run a separately approved count-only HSC S19A random-lens query plan."""
from __future__ import annotations

import argparse
import csv
import getpass
import json
from pathlib import Path
import platform
import sys
import time
import urllib.error

from prepare_wl_random_shape_count_plan import (
    OUTPUT_COLUMNS, normalize_centers, shape_count_query, sha,
)
from fetch_ancillary_catalogues import Cas, catalogue_reader, save


DEFAULT_USER = "kjhan0606"


def validate_plan(plan_dir):
    plan_dir = Path(plan_dir).resolve()
    path = plan_dir / "random_shape_count_plan.json"
    manifest = json.loads(path.read_text(encoding="utf-8"))
    if manifest.get("schema") != "hsc-s19a-random-shape-count-plan-v1":
        raise ValueError("unexpected random shape-count plan schema")
    if manifest.get("status") != "prepared_not_submitted":
        raise ValueError("random shape-count plan is not in prepared state")
    if manifest.get("output_columns") != OUTPUT_COLUMNS:
        raise ValueError("unexpected random shape-count output columns")
    if any(manifest.get(key) is not False for key in (
            "source_shape_rows_returned", "photoz_rows_returned", "images_returned",
            "network_used", "query_submitted", "science_measurement")):
        raise ValueError("plan requests data outside its count-only scope")
    if manifest.get("execution_requires_explicit_user_approval") is not True:
        raise ValueError("plan does not require explicit approval")
    centers_path = Path(manifest["source_centers"]).resolve()
    if sha(centers_path) != manifest.get("source_centers_sha256"):
        raise ValueError("source random-center manifest changed")
    center_doc = json.loads(centers_path.read_text(encoding="utf-8"))
    centers, h0, omega_m, n_lenses = normalize_centers(
        center_doc, radius_mpc=float(manifest["radius_mpc_proper"]),
        batch_size=int(manifest["batch_size"]),
    )
    if (len(centers) != manifest.get("n_random_centers") or
            n_lenses != manifest.get("n_matched_lenses") or
            h0 != manifest["cosmology"]["h0_km_s_mpc"] or
            omega_m != manifest["cosmology"]["omega_m"]):
        raise ValueError("source-center counts or cosmology disagree with plan")
    batches = manifest.get("batches", [])
    expected_batches = (len(centers) + int(manifest["batch_size"]) - 1) // int(manifest["batch_size"])
    if len(batches) != expected_batches or len(batches) != manifest.get("n_query_batches"):
        raise ValueError("random shape-count batch count mismatch")
    for index, batch in enumerate(batches):
        sql_path = Path(batch["path"]).resolve()
        if sql_path.parent != plan_dir or sql_path.name != f"random_shape_counts_{index:03d}.sql":
            raise ValueError("random shape-count SQL path escapes plan directory")
        start = index * int(manifest["batch_size"])
        chunk = centers[start:start + int(manifest["batch_size"])]
        expected_sql = "\nUNION ALL\n".join(shape_count_query(row) for row in chunk) + "\n"
        if sql_path.read_text(encoding="utf-8") != expected_sql or sha(sql_path) != batch.get("sha256"):
            raise ValueError(f"random shape-count SQL changed: {sql_path.name}")
        if batch.get("n_random_centers") != len(chunk):
            raise ValueError(f"random center count changed in {sql_path.name}")
    return manifest, centers, batches, sha(path)


def read_count_rows(path, expected_centers):
    expected = {row["random_center_id"]: row["matched_lens_id"] for row in expected_centers}
    observed = {}
    with Path(path).open(encoding="utf-8-sig", newline="") as stream:
        reader = catalogue_reader(stream)
        if reader.fieldnames != OUTPUT_COLUMNS:
            raise ValueError("CAS returned columns outside the approved count-only schema")
        for row in reader:
            center_id = str(row["random_center_id"])
            matched_id = str(row["matched_lens_id"])
            if center_id not in expected:
                raise ValueError("CAS returned an unrequested random-center ID")
            if matched_id != expected[center_id]:
                raise ValueError("CAS returned a changed random-center to lens association")
            if center_id in observed:
                raise ValueError("CAS returned duplicate random-center counts")
            try:
                count = int(row["n_shape_sources"])
            except (ValueError, TypeError):
                raise ValueError("CAS returned a non-integer source count") from None
            if count < 0 or str(count) != str(row["n_shape_sources"]).strip():
                raise ValueError("CAS returned an invalid non-negative source count")
            observed[center_id] = {"matched_lens_id": matched_id, "n_shape_sources": count}
    if set(observed) != set(expected):
        raise ValueError("CAS count rows do not cover all requested random centers")
    return observed


def cas_status(client, job_id):
    for attempt in range(10):
        try:
            return client.json("status", id=job_id)
        except urllib.error.HTTPError as error:
            if error.code not in (502, 503, 504) or attempt == 9:
                raise
            print(f"CAS status gateway HTTP {error.code}; retrying in 60 seconds ({attempt + 1}/10)", flush=True)
            time.sleep(60)


def run(plan_dir, results_dir, user, *, submit=False, approved=False):
    plan_dir, results_dir = Path(plan_dir).resolve(), Path(results_dir).resolve()
    manifest, centers, batches, plan_sha = validate_plan(plan_dir)
    print(json.dumps({
        "plan_status": manifest["status"], "network_used": False,
        "random_centers": len(centers), "matched_lenses": manifest["n_matched_lenses"],
        "batches": len(batches), "output_columns": OUTPUT_COLUMNS,
        "source_rows": False, "images": False,
    }, ensure_ascii=False), flush=True)
    if not submit:
        return 0
    if not approved:
        raise ValueError("submission requires --approved-random-shape-count-query after explicit user approval")
    if not sys.stdin.isatty():
        raise ValueError("submission requires one hidden password prompt in an interactive terminal")
    if not platform.node().lower().startswith("lageunha"):
        raise ValueError("run the approved query from the established lageunha interactive environment")
    if results_dir.exists():
        checkpoint = results_dir / "count_status.json"
        if not checkpoint.is_file():
            raise FileExistsError("refusing existing count-results directory without a checkpoint")
        state = json.loads(checkpoint.read_text(encoding="utf-8"))
        if state.get("plan_sha256") != plan_sha:
            raise ValueError("existing checkpoint belongs to another plan")
    else:
        results_dir.mkdir(parents=True, exist_ok=False)
        checkpoint = results_dir / "count_status.json"
        state = {
            "schema": "hsc-s19a-random-shape-count-run-v1", "status": "pending",
            "plan_sha256": plan_sha, "jobs": {}, "credentials_saved": False,
            "source_shape_rows_returned": False, "photoz_rows_returned": False,
            "images_returned": False, "science_measurement": False,
        }
        save(checkpoint, state)
    if state.get("credentials_saved") is not False or any(
            state.get(key) is not False for key in (
                "source_shape_rows_returned", "photoz_rows_returned", "images_returned", "science_measurement")):
        raise ValueError("checkpoint contains state outside the count-only scope")
    if any(row.get("status") in ("submitting", "error") and "job_id" not in row
           for row in state.get("jobs", {}).values()):
        raise ValueError("submission outcome is ambiguous; inspect CAS before retry")

    import fcntl
    with (results_dir / "run.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        password = getpass.getpass("HSC password (hidden; one prompt): ")
        if not password:
            raise ValueError("empty password")
        client = Cas(user, password)
        state["status"] = "running"
        save(checkpoint, state)
        try:
            size = int(manifest["batch_size"])
            for index, batch in enumerate(batches):
                name = f"random_shape_counts_{index:03d}"
                item = state["jobs"].setdefault(name, {"status": "pending"})
                target = results_dir / f"{name}.csv"
                expected = centers[index * size:(index + 1) * size]
                if item.get("status") == "complete":
                    if not target.is_file() or sha(target) != item.get("sha256"):
                        raise ValueError(f"saved count batch missing or changed: {name}")
                    read_count_rows(target, expected)
                    continue
                if item.get("status") == "error":
                    raise ValueError(f"CAS job failed for {name}; inspect before retry")
                if "job_id" not in item:
                    item["status"] = "submitting"
                    save(checkpoint, state)
                    try:
                        submitted = client.json(
                            "submit", catalog_job=dict(sql=Path(batch["path"]).read_text(encoding="utf-8"),
                                out_format="csv", include_metainfo_to_body=False, release_version="pdr3"),
                            nomail=True, skip_syntax_check=False,
                        )
                    except urllib.error.HTTPError:
                        item["status"] = "pending"
                        save(checkpoint, state)
                        raise
                    item.update(status="submitted", job_id=str(submitted["id"]))
                    save(checkpoint, state)
                while True:
                    status = cas_status(client, item["job_id"])
                    if status.get("status") == "done":
                        break
                    if status.get("status") == "error":
                        item.update(status="error", query_error=str(status.get("error", "CAS query failed"))[:500])
                        save(checkpoint, state)
                        raise ValueError(f"CAS count query failed for {name}")
                    print(f"{name}: {status.get('status', 'unknown')}", flush=True)
                    time.sleep(60)
                partial = results_dir / f"{name}.csv.part"
                if partial.exists():
                    partial.rename(results_dir / f"{name}.csv.incomplete-{time.time_ns()}")
                with client.post("download", id=item["job_id"]) as response, partial.open("xb") as stream:
                    for block in iter(lambda: response.read(1024 * 1024), b""):
                        stream.write(block)
                observed = read_count_rows(partial, expected)
                partial.replace(target)
                item.update(status="complete", rows=len(observed), bytes=target.stat().st_size, sha256=sha(target))
                save(checkpoint, state)
                print(f"{name}: verified {len(observed)} counts", flush=True)
            all_counts = {}
            for index in range(len(batches)):
                name = f"random_shape_counts_{index:03d}"
                rows_path = results_dir / f"{name}.csv"
                expected = centers[index * size:(index + 1) * size]
                all_counts.update(read_count_rows(rows_path, expected))
            summary = {
                "schema": "hsc-s19a-random-shape-count-summary-v1", "status": "complete",
                "n_random_centers": len(all_counts), "n_matched_lenses": manifest["n_matched_lenses"],
                "sum_shape_source_counts": sum(row["n_shape_sources"] for row in all_counts.values()),
                "zero_source_centers": sum(row["n_shape_sources"] == 0 for row in all_counts.values()),
                "source_shape_rows_returned": False, "photoz_rows_returned": False,
                "images_returned": False, "science_measurement": False,
                "plan_sha256": plan_sha,
            }
            state.update(status="complete", summary=summary)
            save(checkpoint, state)
            (results_dir / "count_summary.json").write_text(
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


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan-dir", type=Path, required=True)
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--user", default=DEFAULT_USER)
    parser.add_argument("--check", action="store_true", help="validate plan only; no network")
    parser.add_argument("--submit", action="store_true", help="submit the count-only plan")
    parser.add_argument("--approved-random-shape-count-query", action="store_true",
                        help="required after explicit user approval; does not approve source-row export")
    args = parser.parse_args()
    if args.check == args.submit:
        parser.error("select exactly one of --check or --submit")
    if args.check and args.approved_random_shape_count_query:
        parser.error("approval flag is not used with --check")
    return run(args.plan_dir, args.results, args.user, submit=args.submit,
               approved=args.approved_random_shape_count_query)


if __name__ == "__main__":
    try:
        sys.exit(main())
    except urllib.error.HTTPError as error:
        print(f"HSC CAS HTTP {error.code}; credentials and response body suppressed.", file=sys.stderr)
        sys.exit(1)
    except KeyboardInterrupt:
        print("\nStopped; submitted query IDs remain checkpointed.", file=sys.stderr)
        sys.exit(130)
