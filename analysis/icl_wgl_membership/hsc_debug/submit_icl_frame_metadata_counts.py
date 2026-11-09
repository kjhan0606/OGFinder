#!/usr/bin/env python3
"""Validate or submit the separately approved ICL frame-metadata count plan.

The submitted SQL returns one row per existing target/control image center and
only patch identifiers plus contributing frame/visit/CCD counts. It does not
return frame gain/exposure rows, images, shapes, or photo-z values.
"""
from __future__ import annotations

import argparse
import fcntl
import getpass
import hashlib
import json
from pathlib import Path
import platform
import sys
import time
import urllib.error

from fetch_ancillary_catalogues import Cas, catalogue_reader, save
from prepare_icl_frame_metadata_counts import (
    BATCH_SIZE,
    CONTROL_MANIFEST,
    CONTROL_STATUS,
    TARGET_MANIFEST,
    TARGET_STATUS,
    load_group,
)


DEFAULT_USER = "kjhan0606"
OUTPUT_COLUMNS = [
    "center_id", "sample_role", "ra_deg", "dec_deg", "band", "skymap_id",
    "tract", "patch", "patch_s", "n_patch_frame_memberships",
    "n_unique_frames", "n_unique_visits", "n_unique_ccds",
]
COUNT_COLUMNS = OUTPUT_COLUMNS[-4:]
STATUS_POLL_INTERVAL_S = 20


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_expected_centers() -> list[dict]:
    centers = load_group(TARGET_MANIFEST, TARGET_STATUS, "target")
    centers += load_group(CONTROL_MANIFEST, CONTROL_STATUS, "control")
    centers.sort(key=lambda row: row["center_id"])
    if len(centers) != 1598 or len({row["center_id"] for row in centers}) != 1598:
        raise ValueError("expected 800 target and 798 unique matched-control centers")
    return centers


def validate_plan(plan_dir: Path) -> tuple[dict, list[dict], list[dict], str]:
    plan_dir = Path(plan_dir).resolve()
    plan_path = plan_dir / "plan.json"
    plan = json.loads(plan_path.read_text(encoding="utf-8"))
    if (plan.get("schema") != "astrafek-icl-frame-metadata-count-plan-v1"
            or plan.get("status") != "prepared_not_submitted"):
        raise ValueError("unexpected or already-submitted ICL metadata plan")
    if plan.get("tables") != ["pdr3_wide.mosaic", "pdr3_wide.mosaicframe"]:
        raise ValueError("metadata plan references an unexpected table set")
    if plan.get("output_columns") != OUTPUT_COLUMNS:
        raise ValueError("metadata plan output columns differ from the count-only schema")
    if (plan.get("n_target_images") != 800 or plan.get("n_control_images") != 798
            or plan.get("n_centers") != 1598 or plan.get("batch_size_centers") != BATCH_SIZE
            or plan.get("rows_returned_per_complete_query") !=
               "one count-only row per input image center (1,598 expected)"):
        raise ValueError("metadata plan sample or batch contract changed")
    if (plan.get("source_metadata_columns_not_returned_in_this_count_stage") != [
            "frame.exptime", "frame.gain1", "frame.gain2", "frame.gain3",
            "frame.gain4", "frame.fluxmag0", "frame.sigma_sky"]):
        raise ValueError("metadata plan no longer separates count and frame-row exports")
    for key in ("images_returned", "shape_rows_returned", "photoz_rows_returned",
                "network_used", "query_submitted", "credentials_saved", "science_measurement"):
        if plan.get(key) is not False:
            raise ValueError(f"metadata plan scope flag changed: {key}")
    if plan.get("execution_requires_explicit_user_approval") is not True:
        raise ValueError("metadata query plan does not require explicit user approval")

    for path_text, expected_hash in plan.get("input_sha256", {}).items():
        path = Path(path_text).resolve()
        if not path.is_file() or sha256(path) != expected_hash:
            raise ValueError(f"metadata plan input changed or missing: {path}")
    centers = load_expected_centers()
    expected_batches = (len(centers) + BATCH_SIZE - 1) // BATCH_SIZE
    batches = plan.get("batches", [])
    if (len(batches) != expected_batches or len(batches) != plan.get("n_query_batches")
            or len(batches) != 16):
        raise ValueError("metadata plan batch count mismatch")
    total_bytes = 0
    for index, batch in enumerate(batches):
        sql_path = Path(batch.get("path", "")).resolve()
        if sql_path.parent != plan_dir or sql_path.name != f"frame_metadata_counts_{index:03d}.sql":
            raise ValueError("metadata SQL path escapes its approved plan directory")
        start = index * BATCH_SIZE
        expected = centers[start:start + BATCH_SIZE]
        if (batch.get("n_centers") != len(expected)
                or batch.get("first_center_id") != expected[0]["center_id"]
                or batch.get("last_center_id") != expected[-1]["center_id"]
                or not sql_path.is_file()
                or sql_path.stat().st_size != batch.get("sql_bytes")
                or sha256(sql_path) != batch.get("sha256")):
            raise ValueError(f"metadata SQL batch changed: frame_metadata_counts_{index:03d}.sql")
        sql = sql_path.read_text(encoding="utf-8")
        forbidden = ("gain1", "gain2", "gain3", "gain4", "exptime", "fluxmag0", "sigma_sky")
        if any(token in sql.lower() for token in forbidden):
            raise ValueError("count-only SQL contains a frame-level gain/exposure field")
        total_bytes += sql_path.stat().st_size
    if total_bytes != sum(int(batch["sql_bytes"]) for batch in batches):
        raise ValueError("metadata plan SQL byte total is inconsistent")
    return plan, centers, batches, sha256(plan_path)


def parse_count(value: str, column: str) -> int:
    try:
        number = int(value)
    except (TypeError, ValueError):
        raise ValueError(f"CAS returned a non-integer {column}") from None
    if number < 0 or str(number) != str(value).strip():
        raise ValueError(f"CAS returned an invalid non-negative {column}")
    return number


def read_result(path: Path, expected_centers: list[dict]) -> dict[str, dict]:
    expected = {row["center_id"]: row for row in expected_centers}
    observed = {}
    with Path(path).open(encoding="utf-8-sig", newline="") as stream:
        reader = catalogue_reader(stream)
        if reader.fieldnames != OUTPUT_COLUMNS:
            raise ValueError("CAS returned columns outside the approved count-only schema")
        for row in reader:
            center_id = str(row["center_id"])
            source = expected.get(center_id)
            if source is None or center_id in observed:
                raise ValueError("CAS returned an unrequested or duplicate center")
            if (row["sample_role"] != source["role"] or row["band"] != source["band"]
                    or abs(float(row["ra_deg"]) - source["ra_deg"]) > 1e-8
                    or abs(float(row["dec_deg"]) - source["dec_deg"]) > 1e-8):
                raise ValueError("CAS changed a requested center identity or coordinate")
            parsed = dict(row)
            for column in COUNT_COLUMNS:
                parsed[column] = parse_count(row[column], column)
            if (parsed["n_unique_frames"] > parsed["n_patch_frame_memberships"]
                    or parsed["n_unique_visits"] > parsed["n_unique_frames"]
                    or parsed["n_unique_ccds"] > parsed["n_patch_frame_memberships"]):
                raise ValueError("CAS returned inconsistent frame/visit/CCD counts")
            observed[center_id] = parsed
    if set(observed) != set(expected):
        raise ValueError("CAS result does not contain exactly one row per requested center")
    return observed


def cas_status(client: Cas, job_id: str) -> dict:
    for attempt in range(10):
        try:
            return client.json("status", id=job_id)
        except urllib.error.HTTPError as error:
            if error.code not in (502, 503, 504) or attempt == 9:
                raise
            print(f"CAS status gateway HTTP {error.code}; retrying in 60 seconds ({attempt + 1}/10)", flush=True)
            time.sleep(60)


def run(plan_dir: Path, results_dir: Path, user: str, *, submit: bool, approved: bool) -> int:
    plan_dir, results_dir = Path(plan_dir).resolve(), Path(results_dir).resolve()
    plan, centers, batches, plan_sha = validate_plan(plan_dir)
    print(json.dumps({
        "plan_status": plan["status"], "plan_sha256": plan_sha,
        "network_used": False, "centers": len(centers), "batches": len(batches),
        "sql_bytes": sum(int(batch["sql_bytes"]) for batch in batches),
        "output_columns": OUTPUT_COLUMNS, "images": False, "shape_rows": False,
        "photoz_rows": False,
    }, ensure_ascii=False), flush=True)
    if not submit:
        return 0
    if not approved:
        raise ValueError("submission requires --approved-icl-frame-metadata-count-query after explicit user approval")
    if not sys.stdin.isatty():
        raise ValueError("submission requires one hidden password prompt in an interactive terminal")
    if not platform.node().lower().startswith("lageunha"):
        raise ValueError("run the approved query from the established lageunha interactive environment")
    if results_dir.exists():
        checkpoint = results_dir / "status.json"
        if not checkpoint.is_file():
            raise FileExistsError("refusing an existing result directory without a matching checkpoint")
        state = json.loads(checkpoint.read_text(encoding="utf-8"))
        if state.get("plan_sha256") != plan_sha:
            raise ValueError("existing metadata checkpoint belongs to another plan")
    else:
        results_dir.mkdir(parents=True, exist_ok=False)
        checkpoint = results_dir / "status.json"
        state = {
            "schema": "astrafek-icl-frame-metadata-count-run-v1", "status": "pending",
            "plan_sha256": plan_sha, "jobs": {}, "network_used": False,
            "query_submitted": False, "credentials_saved": False,
            "images_returned": False, "shape_rows_returned": False,
            "photoz_rows_returned": False, "science_measurement": False,
        }
        save(checkpoint, state)
    if state.get("credentials_saved") is not False or any(
            state.get(key) is not False for key in (
                "images_returned", "shape_rows_returned", "photoz_rows_returned", "science_measurement")):
        raise ValueError("checkpoint contains state outside the count-only scope")
    if any(item.get("status") in ("submitting", "error") and "job_id" not in item
           for item in state.get("jobs", {}).values()):
        raise ValueError("submission outcome is ambiguous; inspect CAS before retry")

    with (results_dir / "run.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        password = getpass.getpass("HSC password (hidden; one prompt): ")
        if not password:
            raise ValueError("empty password")
        client = Cas(user, password)
        state.update(status="running", network_used=True)
        save(checkpoint, state)
        try:
            for index, batch in enumerate(batches):
                name = batch["name"]
                item = state["jobs"].setdefault(name, {"status": "pending"})
                target = results_dir / f"{name}.csv"
                expected = centers[index * BATCH_SIZE:(index + 1) * BATCH_SIZE]
                if item.get("status") == "complete":
                    if not target.is_file() or sha256(target) != item.get("sha256"):
                        raise ValueError(f"saved metadata batch missing or changed: {name}")
                    read_result(target, expected)
                    continue
                if item.get("status") == "error":
                    raise ValueError(f"CAS job failed for {name}; inspect before retry")
                if "job_id" not in item:
                    item["status"] = "submitting"
                    save(checkpoint, state)
                    try:
                        submitted = client.json(
                            "submit",
                            catalog_job=dict(
                                sql=(plan_dir / f"{name}.sql").read_text(encoding="utf-8"),
                                out_format="csv", include_metainfo_to_body=False,
                                release_version="pdr3",
                            ),
                            nomail=True, skip_syntax_check=False,
                        )
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
                        raise ValueError(f"CAS metadata count query failed for {name}")
                    print(f"{name}: {status.get('status', 'unknown')}", flush=True)
                    time.sleep(STATUS_POLL_INTERVAL_S)
                partial = results_dir / f"{name}.csv.part"
                if partial.exists():
                    partial.rename(results_dir / f"{name}.csv.incomplete-{time.time_ns()}")
                with client.post("download", id=item["job_id"]) as response, partial.open("xb") as stream:
                    for block in iter(lambda: response.read(1024 * 1024), b""):
                        stream.write(block)
                observed = read_result(partial, expected)
                partial.replace(target)
                item.update(status="complete", rows=len(observed), bytes=target.stat().st_size,
                            sha256=sha256(target))
                save(checkpoint, state)
                print(f"{name}: verified {len(observed)} count rows", flush=True)

            all_rows = {}
            for index, batch in enumerate(batches):
                path = results_dir / f"{batch['name']}.csv"
                expected = centers[index * BATCH_SIZE:(index + 1) * BATCH_SIZE]
                all_rows.update(read_result(path, expected))
            if set(all_rows) != {row["center_id"] for row in centers}:
                raise ValueError("completed metadata counts do not cover all 1,598 image centers")
            summary = {
                "schema": "astrafek-icl-frame-metadata-count-summary-v1",
                "status": "complete", "plan_sha256": plan_sha,
                "n_centers": len(all_rows), "n_target_centers": 800,
                "n_control_centers": 798,
                "centers_without_patch": sum(not row["skymap_id"] for row in all_rows.values()),
                "sum_patch_frame_memberships": sum(row["n_patch_frame_memberships"] for row in all_rows.values()),
                "sum_unique_frames_per_center": sum(row["n_unique_frames"] for row in all_rows.values()),
                "sum_unique_visits_per_center": sum(row["n_unique_visits"] for row in all_rows.values()),
                "images_returned": False, "shape_rows_returned": False,
                "photoz_rows_returned": False, "credentials_saved": False,
                "science_measurement": False,
            }
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
    parser.add_argument("--plan-dir", type=Path, required=True)
    parser.add_argument("--results", type=Path, required=True)
    parser.add_argument("--user", default=DEFAULT_USER)
    parser.add_argument("--check", action="store_true", help="validate exact SQL and input hashes; no network")
    parser.add_argument("--submit", action="store_true", help="submit the count-only plan")
    parser.add_argument("--approved-icl-frame-metadata-count-query", action="store_true",
                        help="required after explicit approval of this exact 1,598-center plan")
    args = parser.parse_args()
    if args.check == args.submit:
        parser.error("select exactly one of --check or --submit")
    if args.check and args.approved_icl_frame_metadata_count_query:
        parser.error("approval flag cannot be used with --check")
    return run(args.plan_dir, args.results, args.user, submit=args.submit,
               approved=args.approved_icl_frame_metadata_count_query)


if __name__ == "__main__":
    try:
        sys.exit(main())
    except urllib.error.HTTPError as error:
        print(f"HSC CAS HTTP {error.code}; credentials and response body suppressed.", file=sys.stderr)
        sys.exit(1)
    except KeyboardInterrupt:
        print("\nStopped; submitted query IDs remain checkpointed.", file=sys.stderr)
        sys.exit(130)
