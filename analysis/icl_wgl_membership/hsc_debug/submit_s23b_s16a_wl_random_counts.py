#!/usr/bin/env python3
"""Submit only the exact user-approved S16A WL random-center count plan.

The query returns a random-center ID and one integer distinct-object count.
It does not fetch source rows, positions, shapes, photo-z, images, or PSF data.
The credential is requested through one hidden prompt and is never persisted.
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
    canonical_sha,
)
from prepare_s23b_s16a_wl_random_count_plan import (
    EXPECTED_CANDIDATE_AUDIT_SHA256,
    EXPECTED_RANDOM_CENTERS,
    EXPECTED_RANDOM_CENTERS_SHA256,
    EXPECTED_SELECTION_SUMMARY_SHA256,
    EXPECTED_TARGET_SELECTION_SHA256,
    OUTPUT as DEFAULT_PLAN,
    RANDOM_CENTERS,
    SELECTION_SUMMARY,
    TARGET_MANIFEST,
    TARGET_SELECTION,
    BATCH_SIZE,
    OUTPUT_COLUMNS,
    QUERY_TABLES,
    RADIUS_MPC,
    count_query,
)

EXPECTED_PLAN_SHA256 = "dfd2fc75650fc08ce4acfe1eef6d0ed476a7bdf8c84b5d083b607f4009678de0"
EXPECTED_BATCHES = 23
EXPECTED_LENSES = 114
EXPECTED_RELEASE = "pdr1"


def validate_plan(plan_dir: Path) -> tuple[dict, str]:
    plan_dir = plan_dir.resolve()
    plan_path = plan_dir / "shape_count_plan.json"
    plan = json.loads(plan_path.read_text(encoding="utf-8"))
    if (plan.get("schema") != "astrafek-s23b-s16a-wl-random-count-plan-v1"
            or plan.get("status") != "prepared_not_submitted"):
        raise ValueError("unexpected or already submitted random-count plan")
    if (plan.get("execution_requires_explicit_user_approval") is not True
            or plan.get("network_used") is not False
            or plan.get("query_submitted") is not False
            or plan.get("source_rows_downloaded") is not False
            or plan.get("images_returned") is not False
            or plan.get("photoz_returned") is not False
            or plan.get("science_measurement") is not False
            or plan.get("controls_ready_for_science") is not False):
        raise ValueError("random-count plan transfer or science metadata is inconsistent")

    plan_hash = canonical_sha(plan)
    if plan_hash != EXPECTED_PLAN_SHA256 or plan.get("query_plan_sha256") != plan_hash:
        raise ValueError("plan differs from the exact reviewed S16A count-only scope")
    hash_path = plan_dir / "PLAN_SHA256"
    if not hash_path.is_file() or hash_path.read_text(encoding="ascii").split()[0] != plan_hash:
        raise ValueError("PLAN_SHA256 does not match the exact reviewed plan")

    pinned = (
        (Path(plan.get("selection_summary", "")), SELECTION_SUMMARY,
         EXPECTED_SELECTION_SUMMARY_SHA256, "selection summary"),
        (Path(plan.get("random_center_table", "")), RANDOM_CENTERS,
         EXPECTED_RANDOM_CENTERS_SHA256, "random-center table"),
        (Path(plan.get("target_selection_ledger", "")), TARGET_SELECTION,
         EXPECTED_TARGET_SELECTION_SHA256, "target-selection ledger"),
        (Path(plan.get("source_target_manifest", "")), TARGET_MANIFEST,
         EXPECTED_TARGET_SHA256, "target manifest"),
    )
    for saved_path, expected_path, expected_hash, label in pinned:
        if (saved_path.resolve() != expected_path.resolve()
                or not expected_path.is_file()
                or sha(expected_path) != expected_hash):
            raise ValueError(f"pinned {label} is missing or changed")
    if (plan.get("selection_summary_sha256") != EXPECTED_SELECTION_SUMMARY_SHA256
            or plan.get("random_center_table_sha256") != EXPECTED_RANDOM_CENTERS_SHA256
            or plan.get("target_selection_ledger_sha256") != EXPECTED_TARGET_SELECTION_SHA256
            or plan.get("source_target_manifest_sha256") != EXPECTED_TARGET_SHA256
            or plan.get("candidate_audit_csv_sha256") != EXPECTED_CANDIDATE_AUDIT_SHA256):
        raise ValueError("plan provenance hashes differ from the selected random-center sample")

    if (plan.get("cas_release_version") != EXPECTED_RELEASE
            or plan.get("radius_mpc_proper") != RADIUS_MPC
            or plan.get("centers_per_lens") != 5
            or plan.get("n_candidate_lenses") != 116
            or plan.get("n_selected_lenses") != EXPECTED_LENSES
            or plan.get("n_excluded_lenses") != 2
            or plan.get("n_random_centers") != EXPECTED_RANDOM_CENTERS
            or plan.get("n_query_batches") != EXPECTED_BATCHES
            or plan.get("batch_size") != BATCH_SIZE
            or plan.get("query_tables") != QUERY_TABLES
            or plan.get("output_columns") != OUTPUT_COLUMNS
            or plan.get("spatial_function")
            != "coneSearch(c.icoord, ra_deg, dec_deg, radius_arcsec)"
            or plan.get("estimated_result_rows") != EXPECTED_RANDOM_CENTERS):
        raise ValueError("release, center count, batch, or output scope changed")
    if plan.get("estimated_result_scope") != (
            "570 random-center IDs and integer catalog counts only; no source IDs, "
            "coordinates returned, shapes, photo-z values, images, or PSF values"):
        raise ValueError("result-scope safeguard changed")

    centers = plan.get("centers", [])
    if len(centers) != EXPECTED_RANDOM_CENTERS:
        raise ValueError("plan does not contain exactly 570 random centers")
    center_ids = [str(center["random_center_id"]) for center in centers]
    if len(set(center_ids)) != EXPECTED_RANDOM_CENTERS:
        raise ValueError("random-center IDs are not unique")
    for center in centers:
        if (not str(center["random_center_id"]).startswith("RC-")
                or float(center["radius_mpc_proper"]) != RADIUS_MPC
                or float(center["radius_arcsec"]) <= 0
                or not (0 <= float(center["ra_deg"]) < 360)
                or not (-90 <= float(center["dec_deg"]) <= 90)):
            raise ValueError("invalid random-center ID, coordinate, or cone radius")
    per_lens: dict[str, list[dict]] = {}
    for center in centers:
        per_lens.setdefault(str(center["matched_lens_id"]), []).append(center)
    if (len(per_lens) != EXPECTED_LENSES
            or any(len(rows) != 5 for rows in per_lens.values())):
        raise ValueError("plan does not preserve exactly five centers for each selected lens")

    batches = plan.get("batches", [])
    if len(batches) != EXPECTED_BATCHES:
        raise ValueError("random-count SQL batch metadata is incomplete")
    observed_ids = []
    for index, batch in enumerate(batches):
        name = f"s16a_wl_random_counts_{index:03d}.sql"
        query_path = (plan_dir / name).resolve()
        batch_centers = centers[index * BATCH_SIZE:(index + 1) * BATCH_SIZE]
        batch_ids = [str(center["random_center_id"]) for center in batch_centers]
        expected_sql = "\nUNION ALL\n".join(count_query(center) for center in batch_centers) + ";\n"
        if (Path(batch.get("name", "")).name != name
                or query_path.parent != plan_dir
                or not query_path.is_file()
                or sha(query_path) != batch.get("sha256")
                or query_path.stat().st_size != batch.get("sql_bytes")
                or query_path.read_text(encoding="utf-8") != expected_sql
                or batch.get("n_centers") != len(batch_ids)
                or batch.get("first_random_center_id") != batch_ids[0]
                or batch.get("last_random_center_id") != batch_ids[-1]):
            raise ValueError(f"SQL batch differs from the exact count plan: {name}")
        observed_ids.extend(batch_ids)
    if observed_ids != center_ids:
        raise ValueError("SQL batches do not cover every random center exactly once")
    return plan, plan_hash


def read_count_csv(path: Path, expected_ids: list[str]) -> list[dict[str, int | str]]:
    with path.open(encoding="utf-8-sig", newline="") as stream:
        reader = catalogue_reader(stream)
        if reader.fieldnames != OUTPUT_COLUMNS:
            raise ValueError("CAS result columns differ from the approved count-only schema")
        rows = []
        for item in reader:
            if None in item or any(value is None for value in item.values()):
                raise ValueError("CAS returned a malformed count row")
            center_id = str(item["random_center_id"])
            try:
                count = int(item["n_shape_sources"])
            except (TypeError, ValueError) as error:
                raise ValueError("CAS returned a non-integer source count") from error
            if count < 0:
                raise ValueError("CAS returned a negative source count")
            rows.append({"random_center_id": center_id, "n_shape_sources": count})
    ids = [str(row["random_center_id"]) for row in rows]
    if len(ids) != len(expected_ids) or len(set(ids)) != len(ids) or set(ids) != set(expected_ids):
        raise ValueError("CAS result IDs differ from the approved batch")
    return rows


def fetch_count(client, batch: dict, sql: str, expected_ids: list[str], results: Path,
                state: dict, checkpoint: Path) -> None:
    name = Path(batch["name"]).stem
    item = state["jobs"].setdefault(name, {"status": "pending"})
    target = results / f"{name}.csv"
    if item.get("status") == "complete":
        if not target.is_file() or sha(target) != item.get("sha256"):
            raise ValueError(f"saved count result is missing or changed: {name}")
        read_count_csv(target, expected_ids)
        return
    if item.get("status") == "submitting" and "job_id" not in item:
        raise ValueError(f"submission outcome is unknown for {name}; inspect CAS before retry")
    if item.get("status") == "error":
        raise ValueError(f"CAS query failed for {name}; inspect the checkpoint before retry")

    if "job_id" not in item:
        item["status"] = "submitting"
        save(checkpoint, state)
        try:
            submitted = client.json(
                "submit",
                catalog_job=dict(sql=sql, out_format="csv",
                                 include_metainfo_to_body=False,
                                 release_version=EXPECTED_RELEASE),
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
    item.update(status="complete", rows=len(rows), bytes=target.stat().st_size,
                sha256=sha(target), total_sources=sum(counts),
                zero_source_centers=sum(value == 0 for value in counts))
    save(checkpoint, state)
    print(f"{name}: verified {len(rows)} aggregate random-center counts", flush=True)


def combined_tsv(rows: list[dict[str, int | str]]) -> bytes:
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(stream, fieldnames=OUTPUT_COLUMNS, delimiter="\t",
                            lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    return stream.getvalue().encode("utf-8")


def persist_identical(path: Path, content: bytes) -> None:
    if path.exists():
        if path.read_bytes() != content:
            raise ValueError(f"existing output differs from approved results: {path.name}")
        return
    with path.open("xb") as stream:
        stream.write(content)
        stream.flush()


def read_saved_tsv(path: Path, expected_ids: list[str]) -> list[dict[str, int | str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        if reader.fieldnames != OUTPUT_COLUMNS:
            raise ValueError("saved count table has unexpected columns")
        rows = []
        for item in reader:
            rows.append({"random_center_id": str(item["random_center_id"]),
                         "n_shape_sources": int(item["n_shape_sources"])})
    if [row["random_center_id"] for row in rows] != expected_ids:
        raise ValueError("saved result IDs differ from the approved plan")
    return rows


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan-dir", type=Path, default=DEFAULT_PLAN)
    parser.add_argument("--results", type=Path,
                        help="defaults to <plan-dir>/results")
    parser.add_argument("--user", default="kjhan0606")
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--plan", action="store_true",
                      help="validate locally without contacting HSC CAS")
    mode.add_argument("--submit", action="store_true",
                      help="submit only the exact user-approved aggregate-count plan")
    parser.add_argument("--approved-plan-sha256",
                        help="required with --submit")
    args = parser.parse_args()

    plan_dir = args.plan_dir.resolve()
    plan, plan_hash = validate_plan(plan_dir)
    results = (args.results or plan_dir / "results").resolve()
    expected_ids = [str(center["random_center_id"]) for center in plan["centers"]]
    print(json.dumps({
        "plan_status": plan["status"],
        "query_plan_sha256": plan_hash,
        "candidate_lenses": plan["n_candidate_lenses"],
        "selected_lenses": plan["n_selected_lenses"],
        "random_centers": plan["n_random_centers"],
        "batches": plan["n_query_batches"],
        "result_rows": plan["estimated_result_rows"],
        "output_columns": plan["output_columns"],
        "cas_release_version": plan["cas_release_version"],
        "network_used": False if args.plan else "explicitly approved plan only",
        "results": str(results),
    }, ensure_ascii=False), flush=True)
    if args.plan:
        return 0
    if args.approved_plan_sha256 != plan_hash:
        parser.error("--submit requires the exact plan hash approved by the user")
    if not sys.stdin.isatty():
        parser.error("--submit requires a terminal for one hidden password prompt")
    if not platform.node().lower().startswith("lageunha"):
        parser.error("run --submit from the established interactive LagEunha environment")

    checkpoint = results / "count_status.json"
    if results.exists():
        if not checkpoint.is_file():
            parser.error(f"refusing to reuse a result directory without its checkpoint: {results}")
        state = json.loads(checkpoint.read_text(encoding="utf-8"))
        if (state.get("schema") != "astrafek-s23b-s16a-wl-random-count-run-v1"
                or state.get("plan_sha256") != plan_hash):
            parser.error("existing checkpoint belongs to a different count plan")
        if state.get("status") == "complete":
            summary_path, table_path = results / "count_summary.json", results / "count_results.tsv"
            if (not summary_path.is_file() or sha(summary_path) != state.get("summary_sha256")
                    or not table_path.is_file()
                    or sha(table_path) != state.get("count_results_sha256")):
                parser.error("completed count results are missing or changed")
            rows = read_saved_tsv(table_path, expected_ids)
            if len(rows) != EXPECTED_RANDOM_CENTERS:
                parser.error("completed count output has an unexpected row count")
            print(f"Already complete: {table_path}")
            return 0
    else:
        results.mkdir(parents=True, exist_ok=False)
        state = {
            "schema": "astrafek-s23b-s16a-wl-random-count-run-v1",
            "status": "pending",
            "plan_sha256": plan_hash,
            "plan_dir": str(plan_dir),
            "cas_release_version": EXPECTED_RELEASE,
            "jobs": {},
            "credentials_saved": False,
            "source_rows_downloaded": False,
            "images_returned": False,
            "photoz_returned": False,
            "science_measurement": False,
        }
        checkpoint.parent.mkdir(parents=True, exist_ok=True)
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
            for index, batch in enumerate(plan["batches"]):
                batch_centers = plan["centers"][index * BATCH_SIZE:(index + 1) * BATCH_SIZE]
                batch_ids = [str(center["random_center_id"]) for center in batch_centers]
                state["current_query"] = Path(batch["name"]).stem
                save(checkpoint, state)
                fetch_count(client, batch, (plan_dir / batch["name"]).read_text(encoding="utf-8"),
                            batch_ids, results, state, checkpoint)

            rows_by_id = {}
            for index, batch in enumerate(plan["batches"]):
                name = Path(batch["name"]).stem
                batch_ids = [str(center["random_center_id"])
                             for center in plan["centers"][index * BATCH_SIZE:(index + 1) * BATCH_SIZE]]
                for row in read_count_csv(results / f"{name}.csv", batch_ids):
                    center_id = str(row["random_center_id"])
                    if center_id in rows_by_id:
                        raise ValueError(f"duplicate count row for {center_id}")
                    rows_by_id[center_id] = row
            if set(rows_by_id) != set(expected_ids):
                raise ValueError("combined CAS output differs from the exact approved center set")
            rows = [rows_by_id[center_id] for center_id in expected_ids]
            table_path = results / "count_results.tsv"
            persist_identical(table_path, combined_tsv(rows))
            counts = [int(row["n_shape_sources"]) for row in rows]
            summary = {
                "schema": "astrafek-s23b-s16a-wl-random-count-summary-v1",
                "status": "complete_random_center_count_only",
                "query_plan_sha256": plan_hash,
                "cas_release_version": EXPECTED_RELEASE,
                "n_candidate_lenses": plan["n_candidate_lenses"],
                "n_selected_lenses": plan["n_selected_lenses"],
                "n_random_centers": len(rows),
                "n_count_rows": len(rows),
                "sum_random_cone_shape_counts": sum(counts),
                "zero_source_centers": sum(value == 0 for value in counts),
                "positive_source_centers": sum(value > 0 for value in counts),
                "count_results_sha256": sha(table_path),
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
                status="complete", current_query=None,
                total_count_rows=len(rows), sum_random_cone_shape_counts=sum(counts),
                zero_source_centers=summary["zero_source_centers"],
                count_results_sha256=sha(table_path), summary_sha256=sha(summary_path),
                query_submitted=True, source_rows_downloaded=False,
                images_returned=False, photoz_returned=False, science_measurement=False,
            )
            save(checkpoint, state)
        except BaseException as error:
            state.update(status="interrupted" if isinstance(error, KeyboardInterrupt)
                         else "stopped_on_error", error_type=type(error).__name__)
            save(checkpoint, state)
            raise
        finally:
            client.credential.clear()
            password = None

    print(f"Complete: {len(rows)} aggregate counts; no source rows downloaded. Results: {results}")
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
        print("submit_s23b_s16a_wl_random_counts: " + message, file=sys.stderr)
        sys.exit(1)
