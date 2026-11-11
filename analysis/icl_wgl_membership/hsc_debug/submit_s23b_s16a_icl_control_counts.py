#!/usr/bin/env python3
"""Submit the exact approved S16A count-only plan for saved ICL controls.

Only a control ID and one integer distinct-object count are returned for each
of 116 centers. The plan is pinned to its approved SHA-256; passwords are
prompted once, kept in memory, and never written to disk.
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
)
from prepare_s23b_s16a_icl_control_count_plan import (
    AUDIT_SUMMARY,
    CONTROL_TABLE,
    EXPECTED_AUDIT_SUMMARY_SHA256,
    EXPECTED_CONTROL_TABLE_SHA256,
    control_count_query,
)

DEFAULT_PLAN = Path(
    "/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/"
    "wl_s23b_s16a_icl_control_count_plan_2026-11-11_01"
)
EXPECTED_PLAN_SHA256 = "832f2de18f5d512396cc84673e458932635c1f766ea91cdc04c7b47125738dfa"
EXPECTED_TABLES = [
    "s16a_wide.meas2",
    "s16a_wide.weaklensing_hsm_regauss",
    "s16a_wide.meas",
]
EXPECTED_COLUMNS = ["control_id", "n_shape_sources"]
EXPECTED_CONTROLS = 116
EXPECTED_CANDIDATES = 126
EXPECTED_BATCH_SIZE = 25
EXPECTED_RELEASE = "pdr1"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with Path(path).open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def validate_plan(plan_dir: Path) -> tuple[dict, str]:
    plan_dir = plan_dir.resolve()
    plan_path = plan_dir / "shape_count_plan.json"
    plan = json.loads(plan_path.read_text(encoding="utf-8"))
    if (plan.get("schema") != "astrafek-s23b-s16a-icl-control-count-plan-v1"
            or plan.get("status") != "prepared_not_submitted"):
        raise ValueError("unexpected or already submitted ICL control-count plan")
    if (plan.get("execution_requires_explicit_user_approval") is not True
            or plan.get("network_used") is not False
            or plan.get("query_submitted") is not False
            or plan.get("source_rows_downloaded") is not False
            or plan.get("images_returned") is not False
            or plan.get("photoz_returned") is not False
            or plan.get("science_measurement") is not False):
        raise ValueError("control-count transfer or science metadata is inconsistent")

    plan_hash = canonical_sha(plan)
    if plan_hash != EXPECTED_PLAN_SHA256 or plan.get("query_plan_sha256") != plan_hash:
        raise ValueError("control-count plan differs from the exact user-approved plan")
    hash_file = plan_dir / "PLAN_SHA256"
    if not hash_file.is_file() or hash_file.read_text(encoding="ascii").split()[0] != plan_hash:
        raise ValueError("PLAN_SHA256 does not match the approved control-count plan")

    if (plan.get("candidate_control_table") != str(CONTROL_TABLE)
            or plan.get("candidate_control_table_sha256") != EXPECTED_CONTROL_TABLE_SHA256
            or not CONTROL_TABLE.is_file()
            or sha(CONTROL_TABLE) != EXPECTED_CONTROL_TABLE_SHA256
            or plan.get("selection_audit_summary") != str(AUDIT_SUMMARY)
            or plan.get("selection_audit_summary_sha256") != EXPECTED_AUDIT_SUMMARY_SHA256
            or not AUDIT_SUMMARY.is_file()
            or sha(AUDIT_SUMMARY) != EXPECTED_AUDIT_SUMMARY_SHA256):
        raise ValueError("pinned candidate-control inputs are missing or changed")
    if (Path(plan.get("source_target_manifest", "")).resolve() != TARGET_MANIFEST.resolve()
            or not TARGET_MANIFEST.is_file()
            or sha(TARGET_MANIFEST) != EXPECTED_TARGET_SHA256
            or plan.get("source_target_manifest_sha256") != EXPECTED_TARGET_SHA256):
        raise ValueError("pinned S23B target manifest is missing or changed")

    controls = read_tsv(CONTROL_TABLE)
    centers = plan.get("centers", [])
    if (len(controls) != EXPECTED_CONTROLS
            or len(centers) != EXPECTED_CONTROLS
            or plan.get("n_candidate_targets") != EXPECTED_CANDIDATES
            or plan.get("n_controls") != EXPECTED_CONTROLS
            or plan.get("estimated_result_rows") != EXPECTED_CONTROLS):
        raise ValueError("plan does not match the reviewed 126-target/116-control scope")
    expected_ids = []
    for row, center in zip(controls, centers):
        control_id = str(row["control_id"])
        if (center.get("control_id") != control_id
                or center.get("matched_lens_id") != row["target_id"]
                or center.get("icl_filter") != row["icl_filter"]
                or center.get("ra_deg") != float(row["control_ra_deg"])
                or center.get("dec_deg") != float(row["control_dec_deg"])
                or center.get("z_lens") != float(row["z_cl"])
                or center.get("radius_mpc_proper") != 3.0
                or center.get("matched_target_s16a_raw_count") != int(row["s16a_target_count_3pmpc"])
                or center.get("delta_ebv_sfd_target_minus_control")
                != float(row["delta_ebv_sfd_target_minus_control"])):
            raise ValueError(f"control center differs from the approved input table: {control_id}")
        expected_ids.append(control_id)
    if len(set(expected_ids)) != EXPECTED_CONTROLS:
        raise ValueError("approved control IDs are not unique")

    if (plan.get("cas_release_version") != EXPECTED_RELEASE
            or plan.get("radius_mpc_proper") != 3.0
            or plan.get("batch_size") != EXPECTED_BATCH_SIZE
            or plan.get("n_query_batches") != 5
            or plan.get("query_tables") != EXPECTED_TABLES
            or plan.get("output_columns") != EXPECTED_COLUMNS
            or plan.get("spatial_function")
            != "coneSearch(c.icoord, ra_deg, dec_deg, radius_arcsec)"
            or plan.get("estimated_result_bytes_upper_bound") != 4100):
        raise ValueError("release, aperture, batch, or output scope changed")
    if (plan.get("control_interpretation")
            != "This one-to-one ICL-image control screen does not replace the five-random-centers-per-lens weak-lensing control design and is not sufficient for a lensing stack."):
        raise ValueError("control interpretation safeguard changed")

    batches = plan.get("batches", [])
    if len(batches) != 5:
        raise ValueError("control-count SQL batch metadata is incomplete")
    observed_ids = []
    for index, batch in enumerate(batches):
        name = f"s16a_icl_control_counts_{index:03d}.sql"
        query_path = (plan_dir / name).resolve()
        start = index * EXPECTED_BATCH_SIZE
        batch_centers = centers[start:start + EXPECTED_BATCH_SIZE]
        batch_ids = [center["control_id"] for center in batch_centers]
        expected_sql = "\nUNION ALL\n".join(
            control_count_query(center) for center in batch_centers
        ) + ";\n"
        if (Path(batch.get("name", "")).name != name
                or query_path.parent != plan_dir
                or not query_path.is_file()
                or sha(query_path) != batch.get("sha256")
                or query_path.stat().st_size != batch.get("sql_bytes")
                or query_path.read_text(encoding="utf-8") != expected_sql
                or batch.get("n_centers") != len(batch_ids)
                or batch.get("first_control_id") != batch_ids[0]
                or batch.get("last_control_id") != batch_ids[-1]):
            raise ValueError(f"control-count SQL differs from approved plan: {name}")
        observed_ids.extend(batch_ids)
    if observed_ids != expected_ids or sum(batch["n_centers"] for batch in batches) != EXPECTED_CONTROLS:
        raise ValueError("SQL batches do not cover every approved control exactly once")
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
            control_id = str(row["control_id"])
            try:
                count = int(row["n_shape_sources"])
            except (TypeError, ValueError) as error:
                raise ValueError("CAS returned a non-integer count") from error
            if count < 0:
                raise ValueError("CAS returned a negative count")
            rows.append({"control_id": control_id, "n_shape_sources": count})
    observed = [str(row["control_id"]) for row in rows]
    if (len(rows) != len(expected_ids)
            or len(set(observed)) != len(expected_ids)
            or set(observed) != set(expected_ids)):
        raise ValueError("CAS count-result IDs differ from the approved batch")
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
        raise ValueError(f"submission outcome unknown for {name}; inspect CAS before retry")
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
                    release_version="pdr1",
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
    print(f"{name}: verified {len(rows)} integer control counts", flush=True)


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
            raise ValueError(f"existing output differs from approved results: {path.name}")
        return
    with path.open("xb") as stream:
        stream.write(content)
        stream.flush()


def read_saved_tsv(path: Path, expected_ids: list[str]) -> list[dict[str, int | str]]:
    rows = []
    with path.open(encoding="utf-8", newline="") as stream:
        reader = csv.DictReader(stream, delimiter="\t")
        if reader.fieldnames != EXPECTED_COLUMNS:
            raise ValueError("saved count table has unexpected columns")
        for row in reader:
            try:
                count = int(row["n_shape_sources"])
            except (TypeError, ValueError) as error:
                raise ValueError("saved count table has non-integer values") from error
            rows.append({"control_id": str(row["control_id"]), "n_shape_sources": count})
    if [row["control_id"] for row in rows] != expected_ids:
        raise ValueError("saved count IDs differ from the approved plan")
    return rows


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan-dir", type=Path, default=DEFAULT_PLAN)
    parser.add_argument("--results", type=Path, help="defaults to <plan-dir>/results")
    parser.add_argument("--user", default="kjhan0606")
    mode = parser.add_mutually_exclusive_group(required=True)
    mode.add_argument("--plan", action="store_true", help="validate locally without contacting CAS")
    mode.add_argument("--submit", action="store_true", help="submit only the user-approved count plan")
    parser.add_argument("--approved-plan-sha256", help="required with --submit")
    args = parser.parse_args()

    plan_dir = args.plan_dir.resolve()
    plan, plan_hash = validate_plan(plan_dir)
    results = (args.results or plan_dir / "results").resolve()
    expected_ids = [center["control_id"] for center in plan["centers"]]
    print(json.dumps({
        "plan_status": plan["status"],
        "query_plan_sha256": plan_hash,
        "candidate_targets": plan["n_candidate_targets"],
        "control_centers": plan["n_controls"],
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
        parser.error("--submit requires an interactive terminal for one hidden password prompt")
    if not platform.node().lower().startswith("lageunha"):
        parser.error("run --submit from the established interactive LagEunha environment")

    checkpoint = results / "count_status.json"
    if results.exists():
        if not checkpoint.is_file():
            parser.error(f"refusing to reuse result directory without a checkpoint: {results}")
        state = json.loads(checkpoint.read_text(encoding="utf-8"))
        if state.get("schema") != "astrafek-s23b-s16a-icl-control-count-run-v1" or state.get("plan_sha256") != plan_hash:
            parser.error("existing checkpoint belongs to a different plan")
        if state.get("status") == "complete":
            summary_path, table_path = results / "count_summary.json", results / "count_results.tsv"
            if (not summary_path.is_file() or sha(summary_path) != state.get("summary_sha256")
                    or not table_path.is_file() or sha(table_path) != state.get("count_results_sha256")):
                parser.error("completed count outputs are missing or changed")
            rows = read_saved_tsv(table_path, expected_ids)
            if len(rows) != EXPECTED_CONTROLS:
                parser.error("completed result has an unexpected row count")
            print(f"Already complete: {table_path}")
            return 0
    else:
        results.mkdir(parents=True, exist_ok=False)
        state = {
            "schema": "astrafek-s23b-s16a-icl-control-count-run-v1",
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
            for index, batch in enumerate(plan["batches"]):
                start = index * EXPECTED_BATCH_SIZE
                batch_centers = plan["centers"][start:start + EXPECTED_BATCH_SIZE]
                batch_ids = [center["control_id"] for center in batch_centers]
                batch_path = plan_dir / batch["name"]
                state["current_query"] = Path(batch["name"]).stem
                save(checkpoint, state)
                fetch_count(
                    client, batch, batch_path.read_text(encoding="utf-8"),
                    batch_ids, results, state, checkpoint,
                )

            rows_by_id = {}
            for index, batch in enumerate(plan["batches"]):
                name = Path(batch["name"]).stem
                start = index * EXPECTED_BATCH_SIZE
                batch_ids = [
                    center["control_id"]
                    for center in plan["centers"][start:start + EXPECTED_BATCH_SIZE]
                ]
                for row in read_count_csv(results / f"{name}.csv", batch_ids):
                    control_id = str(row["control_id"])
                    if control_id in rows_by_id:
                        raise ValueError(f"duplicate control count result: {control_id}")
                    rows_by_id[control_id] = row
            if set(rows_by_id) != set(expected_ids):
                raise ValueError("combined count results do not match the approved controls")
            rows = [rows_by_id[control_id] for control_id in expected_ids]
            table_path = results / "count_results.tsv"
            persist_identical(table_path, combined_tsv(rows))
            counts = [int(row["n_shape_sources"]) for row in rows]
            summary = {
                "schema": "astrafek-s23b-s16a-icl-control-count-summary-v1",
                "status": "complete_control_count_only",
                "cas_release_version": plan["cas_release_version"],
                "n_candidate_targets": plan["n_candidate_targets"],
                "n_controls": len(rows),
                "n_count_rows": len(rows),
                "sum_control_cone_shape_counts": sum(counts),
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
                sum_control_cone_shape_counts=sum(counts),
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
    print(f"Complete: {len(rows)} aggregate control counts; no source rows downloaded. Results: {results}")
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
        print("submit_s23b_s16a_icl_control_counts: " + message, file=sys.stderr)
        sys.exit(1)
