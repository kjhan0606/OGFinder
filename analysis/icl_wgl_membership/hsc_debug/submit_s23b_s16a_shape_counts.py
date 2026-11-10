#!/usr/bin/env python3
"""Submit the exact approved S23B S16A shape-count plan to HSC CAS.

The script returns one cluster ID and one integer source count per center.
It never requests object-level rows, images, photo-z values, or PSF products.
CAS access requires an interactive lageunha terminal and the exact plan hash
passed with ``--approved-plan-sha256``. The password is prompted once and is
never written to disk.
"""
from __future__ import annotations

import argparse
import csv
import getpass
import hashlib
import io
import json
from pathlib import Path
import platform
import sys
import time
import urllib.error

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from fetch_ancillary_catalogues import Cas, catalogue_reader, save, sha
from prepare_s23b_s16a_shape_count_plan import (
    EXPECTED_TARGET_SHA256,
    TARGET_MANIFEST,
    canonical_sha,
    count_query,
)

DEFAULT_PLAN = Path(
    "/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/"
    "wl_s23b_s16a_shape_count_plan_2026-11-11_01"
)
EXPECTED_PLAN_SHA256 = "c426129b150d440c22fca97d1cf2c722414a7f708132a6267563d77ccda997e0"
EXPECTED_TABLES = [
    "s16a_wide.meas2",
    "s16a_wide.weaklensing_hsm_regauss",
    "s16a_wide.meas",
]
EXPECTED_COLUMNS = ["lens_id", "n_shape_sources"]
EXPECTED_CENTERS = 1043
EXPECTED_BATCH_SIZE = 25
EXPECTED_RELEASE = "pdr1"


def validate_plan(plan_dir: Path) -> tuple[dict, str]:
    plan_dir = plan_dir.resolve()
    plan_path = plan_dir / "shape_count_plan.json"
    plan = json.loads(plan_path.read_text(encoding="utf-8"))
    if plan.get("schema") != "astrafek-s23b-s16a-wl-shape-count-plan-v1":
        raise ValueError("unexpected S16A shape-count plan schema")
    if plan.get("status") != "prepared_not_submitted":
        raise ValueError("S16A plan is not in prepared_not_submitted state")
    if (plan.get("execution_requires_explicit_user_approval") is not True
            or plan.get("network_used") is not False
            or plan.get("query_submitted") is not False
            or plan.get("source_rows_downloaded") is not False
            or plan.get("images_returned") is not False
            or plan.get("photoz_returned") is not False
            or plan.get("science_measurement") is not False):
        raise ValueError("S16A plan transfer or science metadata is inconsistent")

    plan_hash = canonical_sha(plan)
    if plan_hash != EXPECTED_PLAN_SHA256 or plan.get("query_plan_sha256") != plan_hash:
        raise ValueError("S16A plan hash differs from the reviewed plan")
    hash_file = plan_dir / "PLAN_SHA256"
    if not hash_file.is_file() or hash_file.read_text(encoding="ascii").split()[0] != plan_hash:
        raise ValueError("PLAN_SHA256 does not match the reviewed plan")

    manifest_path = Path(plan["source_target_manifest"]).resolve()
    if manifest_path != TARGET_MANIFEST.resolve():
        raise ValueError("S16A plan points to an unexpected target manifest")
    if (not manifest_path.is_file()
            or sha(manifest_path) != EXPECTED_TARGET_SHA256
            or plan.get("source_target_manifest_sha256") != EXPECTED_TARGET_SHA256):
        raise ValueError("pinned 1,043-target manifest is missing or changed")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest_rows = manifest.get("clusters", [])
    centers = plan.get("centers", [])
    if (manifest.get("schema") != "hsc-bcg-diffuse-sample-v1"
            or manifest.get("science_measurement") is not False
            or len(manifest_rows) != EXPECTED_CENTERS
            or len(centers) != EXPECTED_CENTERS
            or plan.get("n_centers") != EXPECTED_CENTERS
            or plan.get("estimated_result_rows") != EXPECTED_CENTERS):
        raise ValueError("S16A plan does not cover the pinned 1,043-center sample")

    expected_ids = []
    for row, center in zip(manifest_rows, centers):
        lens_id = str(row.get("id", ""))
        if (row.get("sample_role") != "target"
                or not lens_id.startswith("S23B_")
                or center.get("lens_id") != lens_id
                or float(center.get("ra_deg")) != float(row["ra_deg"])
                or float(center.get("dec_deg")) != float(row["dec_deg"])
                or float(center.get("z_lens")) != float(row["z_cl"])
                or float(center.get("radius_mpc_proper")) != 3.0):
            raise ValueError(f"S16A center differs from the pinned manifest: {lens_id}")
        expected_ids.append(lens_id)
    if len(set(expected_ids)) != EXPECTED_CENTERS:
        raise ValueError("S16A target IDs are not unique")

    if (plan.get("cas_release_version") != EXPECTED_RELEASE
            or plan.get("radius_mpc_proper") != 3.0
            or plan.get("batch_size") != EXPECTED_BATCH_SIZE
            or plan.get("n_query_batches") != 42
            or plan.get("query_tables") != EXPECTED_TABLES
            or plan.get("output_columns") != EXPECTED_COLUMNS
            or plan.get("spatial_function")
            != "coneSearch(c.icoord, ra_deg, dec_deg, radius_arcsec)"):
        raise ValueError("S16A query release, radius, table, or output scope changed")

    batches = plan.get("batches", [])
    if len(batches) != plan["n_query_batches"]:
        raise ValueError("S16A batch metadata is incomplete")
    observed_ids = []
    for index, batch in enumerate(batches):
        name = f"s16a_shape_counts_{index:03d}.sql"
        query_path = (plan_dir / name).resolve()
        start = index * EXPECTED_BATCH_SIZE
        batch_centers = centers[start:start + EXPECTED_BATCH_SIZE]
        batch_ids = [center["lens_id"] for center in batch_centers]
        expected_sql = "\nUNION ALL\n".join(
            count_query(center) for center in batch_centers
        ) + ";\n"
        if (Path(batch.get("name", "")).name != name
                or query_path.parent != plan_dir
                or not query_path.is_file()
                or sha(query_path) != batch.get("sha256")
                or query_path.stat().st_size != batch.get("sql_bytes")
                or query_path.read_text(encoding="utf-8") != expected_sql
                or batch.get("n_centers") != len(batch_ids)
                or batch.get("first_lens_id") != batch_ids[0]
                or batch.get("last_lens_id") != batch_ids[-1]):
            raise ValueError(f"S16A SQL batch differs from the reviewed plan: {name}")
        observed_ids.extend(batch_ids)
    if observed_ids != expected_ids or sum(b["n_centers"] for b in batches) != EXPECTED_CENTERS:
        raise ValueError("S16A SQL batches do not cover every target exactly once")
    return plan, plan_hash


def read_count_csv(path: Path, expected_ids: list[str]) -> list[dict[str, int | str]]:
    rows = []
    with path.open(encoding="utf-8-sig", newline="") as stream:
        reader = catalogue_reader(stream)
        if reader.fieldnames != EXPECTED_COLUMNS:
            raise ValueError("CAS count result has unexpected columns")
        for row in reader:
            if None in row or any(value is None for value in row.values()):
                raise ValueError("CAS returned a malformed count row")
            lens_id = str(row["lens_id"])
            try:
                count = int(row["n_shape_sources"])
            except (TypeError, ValueError) as error:
                raise ValueError("CAS returned a non-integer source count") from error
            if count < 0:
                raise ValueError("CAS returned a negative source count")
            rows.append({"lens_id": lens_id, "n_shape_sources": count})
    observed = [str(row["lens_id"]) for row in rows]
    if (len(rows) != len(expected_ids)
            or len(set(observed)) != len(expected_ids)
            or set(observed) != set(expected_ids)):
        raise ValueError("CAS count result IDs differ from the approved batch")
    return rows


def fetch_count(client, batch: dict, sql: str, expected_ids: list[str], release: str,
                results: Path, state: dict, checkpoint: Path) -> None:
    name = Path(batch["name"]).stem
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
            submitted = client.json(
                "submit",
                catalog_job=dict(
                    sql=sql,
                    out_format="csv",
                    include_metainfo_to_body=False,
                    release_version=release,
                ),
                nomail=True,
                skip_syntax_check=False,
            )
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
        time.sleep(60)

    partial = results / f"{name}.csv.part"
    if partial.exists():
        partial.rename(results / f"{name}.csv.incomplete-{time.time_ns()}")
    with client.post("download", id=item["job_id"]) as response, partial.open("xb") as stream:
        for block in iter(lambda: response.read(1024 * 1024), b""):
            stream.write(block)
    rows = read_count_csv(partial, expected_ids)
    partial.replace(target)
    counts = [int(row["n_shape_sources"]) for row in rows]
    item.update(
        status="complete",
        rows=len(rows),
        bytes=target.stat().st_size,
        sha256=sha(target),
        total_sources=sum(counts),
        zero_source_centers=sum(value == 0 for value in counts),
    )
    save(checkpoint, state)
    print(f"{name}: verified {len(rows)} aggregate count rows", flush=True)


def combined_tsv(rows: list[dict[str, int | str]]) -> bytes:
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(
        stream, fieldnames=EXPECTED_COLUMNS, delimiter="\t", lineterminator="\n"
    )
    writer.writeheader()
    writer.writerows(rows)
    return stream.getvalue().encode("utf-8")


def persist_identical(path: Path, content: bytes) -> None:
    if path.exists():
        if path.read_bytes() != content:
            raise ValueError(f"existing output differs from the approved count result: {path.name}")
        return
    with path.open("xb") as stream:
        stream.write(content)
        stream.flush()


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan-dir", type=Path, default=DEFAULT_PLAN)
    parser.add_argument("--results", type=Path,
                        help="result directory; defaults to <plan-dir>/results")
    parser.add_argument("--user", default="kjhan0606")
    parser.add_argument("--plan", action="store_true",
                        help="validate the plan locally without contacting HSC CAS")
    parser.add_argument("--submit", action="store_true",
                        help="submit the explicitly approved aggregate-count plan")
    parser.add_argument("--approved-plan-sha256",
                        help="required with --submit and must match the approved plan")
    args = parser.parse_args()
    if args.plan == args.submit:
        parser.error("select exactly one of --plan or --submit")

    plan_dir = args.plan_dir.resolve()
    plan, plan_hash = validate_plan(plan_dir)
    results = (args.results or plan_dir / "results").resolve()
    print(json.dumps({
        "plan_status": plan["status"],
        "query_plan_sha256": plan_hash,
        "centers": plan["n_centers"],
        "batches": plan["n_query_batches"],
        "result_rows": plan["estimated_result_rows"],
        "output_columns": plan["output_columns"],
        "cas_release_version": plan["cas_release_version"],
        "network_used": False if args.plan else "pending explicit approval",
        "results": str(results),
    }, ensure_ascii=False), flush=True)
    if args.plan:
        return 0

    if args.approved_plan_sha256 != plan_hash:
        parser.error("--submit requires the exact plan hash approved by the user")
    if not sys.stdin.isatty():
        parser.error("--submit requires an interactive terminal for one hidden password prompt")
    if not platform.node().lower().startswith("lageunha"):
        parser.error("Run --submit from the established lageunha interactive environment")

    if results.exists():
        checkpoint = results / "count_status.json"
        if not checkpoint.is_file():
            parser.error(f"refusing to reuse a result directory without its checkpoint: {results}")
        state = json.loads(checkpoint.read_text(encoding="utf-8"))
        if (state.get("schema") != "astrafek-s23b-s16a-wl-count-run-v1"
                or state.get("plan_sha256") != plan_hash):
            parser.error("existing count checkpoint belongs to a different plan")
        if state.get("status") == "complete":
            summary_path = results / "count_summary.json"
            table_path = results / "count_results.tsv"
            if (not summary_path.is_file()
                    or sha(summary_path) != state.get("summary_sha256")
                    or not table_path.is_file()
                    or sha(table_path) != state.get("count_results_sha256")):
                parser.error("completed count outputs are missing or changed")
            rows = read_count_csv_from_tsv(table_path, plan)
            if len(rows) != EXPECTED_CENTERS:
                parser.error("completed count table has an unexpected number of rows")
            print(f"Already complete: {table_path}")
            return 0
    else:
        results.mkdir(parents=True, exist_ok=False)
        checkpoint = results / "count_status.json"
        state = {
            "schema": "astrafek-s23b-s16a-wl-count-run-v1",
            "status": "pending",
            "plan_sha256": plan_hash,
            "plan_dir": str(plan_dir),
            "cas_release_version": plan["cas_release_version"],
            "jobs": {},
            "credentials_saved": False,
            "source_rows_downloaded": False,
            "images_returned": False,
            "photoz_returned": False,
            "science_measurement": False,
        }
        save(checkpoint, state)

    import fcntl
    with (results / "run.lock").open("a", encoding="utf-8") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        password = getpass.getpass("HSC password (hidden; one prompt): ")
        if not password:
            raise ValueError("Empty password")
        client = Cas(args.user, password)
        state["status"] = "running"
        save(checkpoint, state)
        try:
            for batch_index, batch in enumerate(plan["batches"]):
                start = batch_index * EXPECTED_BATCH_SIZE
                expected_ids = [
                    center["lens_id"]
                    for center in plan["centers"][start:start + EXPECTED_BATCH_SIZE]
                ]
                batch_path = plan_dir / batch["name"]
                state["current_query"] = Path(batch["name"]).stem
                save(checkpoint, state)
                fetch_count(
                    client,
                    batch,
                    batch_path.read_text(encoding="utf-8"),
                    expected_ids,
                    plan["cas_release_version"],
                    results,
                    state,
                    checkpoint,
                )

            expected_ids = [center["lens_id"] for center in plan["centers"]]
            rows_by_id: dict[str, dict[str, int | str]] = {}
            for batch_index, batch in enumerate(plan["batches"]):
                name = Path(batch["name"]).stem
                batch_path = results / f"{name}.csv"
                batch_start = batch_index * EXPECTED_BATCH_SIZE
                batch_ids = [
                    center["lens_id"]
                    for center in plan["centers"][batch_start:batch_start + EXPECTED_BATCH_SIZE]
                ]
                for row in read_count_csv(batch_path, batch_ids):
                    lens_id = str(row["lens_id"])
                    if lens_id in rows_by_id:
                        raise ValueError(f"duplicate count result for {lens_id}")
                    rows_by_id[lens_id] = row
            if set(rows_by_id) != set(expected_ids):
                raise ValueError("combined counts do not match the approved target list")
            rows = [rows_by_id[lens_id] for lens_id in expected_ids]
            table_path = results / "count_results.tsv"
            persist_identical(table_path, combined_tsv(rows))
            counts = [int(row["n_shape_sources"]) for row in rows]
            summary = {
                "schema": "astrafek-s23b-s16a-wl-count-summary-v1",
                "status": "complete",
                "cas_release_version": plan["cas_release_version"],
                "n_centers": len(rows),
                "n_count_rows": len(rows),
                "sum_lens_cone_shape_counts": sum(counts),
                "zero_source_centers": sum(value == 0 for value in counts),
                "positive_source_centers": sum(value > 0 for value in counts),
                "count_results_sha256": sha(table_path),
                "query_plan_sha256": plan_hash,
                "query_submitted": True,
                "source_rows_downloaded": False,
                "images_returned": False,
                "photoz_returned": False,
                "science_measurement": False,
            }
            summary_path = results / "count_summary.json"
            summary_bytes = (json.dumps(summary, ensure_ascii=False, indent=2,
                                        allow_nan=False) + "\n").encode("utf-8")
            persist_identical(summary_path, summary_bytes)
            state.update(
                status="complete",
                current_query=None,
                total_count_rows=len(rows),
                raw_shape_rows_sum=sum(counts),
                zero_source_centers=summary["zero_source_centers"],
                count_results_sha256=sha(table_path),
                summary_sha256=sha(summary_path),
                query_submitted=True,
                source_rows_downloaded=False,
                images_returned=False,
                photoz_returned=False,
                science_measurement=False,
            )
            save(checkpoint, state)
        except BaseException as error:
            state.update(
                status="interrupted" if isinstance(error, KeyboardInterrupt) else "stopped_on_error",
                error_type=type(error).__name__,
            )
            save(checkpoint, state)
            raise
        finally:
            client.credential.clear()
            password = None
    print(f"Complete: {len(rows)} aggregate count rows; no source rows downloaded. Results: {results}")
    return 0


def read_count_csv_from_tsv(path: Path, plan: dict) -> list[dict[str, int | str]]:
    rows = []
    with path.open(encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        if reader.fieldnames != EXPECTED_COLUMNS:
            raise ValueError("saved count table has unexpected columns")
        rows.extend(dict(row) for row in reader)
    ids = [center["lens_id"] for center in plan["centers"]]
    if [str(row["lens_id"]) for row in rows] != ids:
        raise ValueError("saved count table target IDs differ from the approved plan")
    for row in rows:
        try:
            value = int(row["n_shape_sources"])
        except (TypeError, ValueError) as error:
            raise ValueError("saved count table contains a non-integer count") from error
        if value < 0:
            raise ValueError("saved count table contains a negative count")
        row["n_shape_sources"] = value
    return rows


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
        print("submit_s23b_s16a_shape_counts: " + message, file=sys.stderr)
        sys.exit(1)
