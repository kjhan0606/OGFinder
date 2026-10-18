#!/usr/bin/env python3
"""Submit an explicitly approved HSC S19A shape-row export plan.

--plan is local and network-free. --submit is guarded by an explicit approval
flag, requires an interactive lageunha terminal, and retrieves only the 22
columns in the prepared plan (no images and no full photo-z PDFs).
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
from prepare_wl_shape_rows import (
    DEFAULT_OUTPUT as DEFAULT_PLAN,
    OUTPUT_COLUMNS,
    load_counts,
)
from submit_wl_shape_counts import load_plan

DEFAULT_PREFLIGHT = Path(
    "/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/"
    "wl_shape_preflight_2026-10-18_01"
)
DEFAULT_RESULTS = Path(
    "/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/"
    "wl_shape_rows_export_2026-10-18_03"
)


def validate_export_plan(plan_dir, preflight):
    plan_dir = plan_dir.resolve()
    preflight = preflight.resolve()
    manifest = json.loads((plan_dir / "row_export_plan.json").read_text(encoding="utf-8"))
    if manifest.get("schema") != "hsc-s19a-shape-row-export-plan-v1":
        raise ValueError("unexpected source-row export plan schema")
    if manifest.get("status") != "prepared_not_submitted":
        raise ValueError("source-row plan is not in prepared_not_submitted state")
    if manifest.get("execution_requires_explicit_user_approval") is not True:
        raise ValueError("source-row plan is missing its approval guard")
    if manifest.get("source_rows_downloaded") is not False or manifest.get("photoz_pdf_rows_included") is not False:
        raise ValueError("plan has an unexpected transfer scope")
    query_plan = load_plan(preflight)
    counts, count_state, summary = load_counts(preflight, query_plan)
    if manifest.get("source_count_plan_sha256") != count_state.get("plan_sha256"):
        raise ValueError("source-row plan is tied to a different count-only run")
    if manifest.get("lens_source_pair_rows_expected") != sum(counts.values()):
        raise ValueError("source-row estimate does not match count-only results")
    batch_size = int(query_plan["batch_size"])
    batches = manifest.get("batches", [])
    if len(batches) != 32 or manifest.get("n_query_batches") != 32:
        raise ValueError("expected exactly 32 source-row SQL batches")
    for index, item in enumerate(batches):
        path = Path(item["path"])
        if path.parent.resolve() != plan_dir or path.name != f"shape_rows_{index:03d}.sql":
            raise ValueError("source-row SQL path escapes the prepared plan directory")
        if sha(path) != item.get("sha256"):
            raise ValueError(f"source-row SQL checksum changed: {path.name}")
        sql = path.read_text(encoding="utf-8").strip()
        upper = sql.upper()
        if not upper.startswith("SELECT DISTINCT") or "SELECT *" in upper or "B.*" in upper:
            raise ValueError(f"unsafe or non-minimal SQL projection: {path.name}")
        if "S19A_WIDE.WEAKLENSING_HSM_REGAUSS" not in upper or "S19A_WIDE.MEAS2" not in upper or "S19A_WIDE.MEAS C" not in upper:
            raise ValueError(f"unexpected source tables: {path.name}")
        if "PHOTOZ" in upper or "PDR3_DNNZ_WIDE" in upper:
            raise ValueError(f"photo-z rows are outside this export plan: {path.name}")
        start = index * batch_size
        lens_batch = query_plan["lenses"][start:start + batch_size]
        expected_rows = sum(counts[str(row["lens_id"])] for row in lens_batch)
        if item.get("n_lenses") != len(lens_batch) or item.get("expected_pair_rows_from_count_only") != expected_rows:
            raise ValueError(f"count-only row estimate mismatch: {path.name}")
        for lens in lens_batch:
            marker = "SELECT DISTINCT\n  '" + str(lens["lens_id"]).replace("'", "''") + "' AS LENS_ID"
            if marker.upper() not in upper:
                raise ValueError(f"lens missing from SQL batch: {path.name}")
    return manifest, counts, summary


def validate_rows(path, counts_for_batch):
    observed = {lens_id: 0 for lens_id in counts_for_batch}
    pair_ids = set()
    with Path(path).open(encoding="utf-8-sig", newline="") as stream:
        reader = catalogue_reader(stream)
        if reader.fieldnames != OUTPUT_COLUMNS:
            raise ValueError(f"unexpected source-row columns: {reader.fieldnames}")
        for row in reader:
            if None in row or not row.get("lens_id") or not row.get("object_id"):
                raise ValueError("source-row CSV has a malformed lens or object ID")
            lens_id = str(row["lens_id"])
            object_id = str(row["object_id"])
            if lens_id not in observed:
                raise ValueError("CAS returned a lens outside the requested batch")
            key = (lens_id, object_id)
            if key in pair_ids:
                raise ValueError("duplicate (lens_id, object_id) pair in CAS output")
            pair_ids.add(key)
            observed[lens_id] += 1
    if observed != counts_for_batch:
        bad = [lens_id for lens_id in observed if observed[lens_id] != counts_for_batch[lens_id]]
        raise ValueError("source-row counts disagree with preflight for " + str(len(bad)) + " lenses")
    return len(pair_ids)


def run(args):
    preflight = args.preflight.resolve()
    plan_dir = args.plan_dir.resolve()
    results_dir = args.results.resolve()
    manifest, lens_counts, count_summary = validate_export_plan(plan_dir, preflight)
    print(json.dumps({
        "plan_status": manifest["status"], "network_used": False,
        "lenses": manifest["n_lenses"], "batches": manifest["n_query_batches"],
        "lens_source_pair_rows_expected": manifest["lens_source_pair_rows_expected"],
        "source_rows_downloaded": False, "photoz_pdf_rows_included": False,
        "results": str(results_dir),
    }, ensure_ascii=False), flush=True)
    if args.plan:
        return 0
    if not args.submit:
        raise ValueError("select --plan or --submit")
    if not args.approved_source_row_export:
        raise ValueError("--submit requires --approved-source-row-export after explicit user approval")
    if not sys.stdin.isatty():
        raise ValueError("--submit requires an interactive terminal")
    if not platform.node().lower().startswith("lageunha"):
        raise ValueError("run row export from the established lageunha interactive environment")

    fingerprint = hashlib.sha256((plan_dir / "row_export_plan.json").read_bytes()).hexdigest()
    if results_dir.exists():
        checkpoint = results_dir / "row_status.json"
        if not checkpoint.is_file():
            raise ValueError(f"refusing existing results directory without checkpoint: {results_dir}")
        state = json.loads(checkpoint.read_text(encoding="utf-8"))
        if state.get("plan_sha256") != fingerprint:
            raise ValueError("existing row-export checkpoint belongs to a different plan")
    else:
        results_dir.mkdir(parents=True, exist_ok=False)
        checkpoint = results_dir / "row_status.json"
        state = {
            "schema": "hsc-s19a-shape-row-run-v1", "status": "pending",
            "plan_sha256": fingerprint, "preflight": str(preflight),
            "jobs": {}, "credentials_saved": False,
            "source_rows_download_started": False,
            "source_rows_downloaded": False, "science_measurement": False,
        }
        save(checkpoint, state)

    import fcntl
    with (results_dir / "run.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        password = getpass.getpass("HSC password (hidden; one prompt): ")
        if not password:
            raise ValueError("empty password")
        client = Cas(args.user, password)
        state["status"] = "running"
        save(checkpoint, state)
        try:
            for index, batch in enumerate(manifest["batches"]):
                name = f"shape_rows_{index:03d}"
                item = state["jobs"].setdefault(name, {"status": "pending"})
                target = results_dir / f"{name}.csv"
                expected_ids = [str(row["lens_id"]) for row in query_plan_lenses(preflight, index)]
                expected = {lens_id: lens_counts[lens_id] for lens_id in expected_ids}
                if item.get("status") == "complete":
                    if not target.is_file() or sha(target) != item.get("sha256"):
                        raise ValueError(f"saved source rows are missing or changed: {name}")
                    validate_rows(target, expected)
                    continue
                if item.get("status") in ("submitting", "error") and "job_id" not in item:
                    raise ValueError(f"submission state is ambiguous for {name}; inspect CAS before retry")
                if item.get("status") == "error":
                    raise ValueError(f"CAS job previously failed for {name}; inspect checkpoint before retry")
                if "job_id" not in item:
                    item["status"] = "submitting"
                    state["source_rows_download_started"] = True
                    save(checkpoint, state)
                    try:
                        submitted = client.json("submit", catalog_job=dict(
                            sql=Path(batch["path"]).read_text(encoding="utf-8"),
                            out_format="csv", include_metainfo_to_body=False,
                            release_version="pdr3"), nomail=True, skip_syntax_check=False)
                    except urllib.error.HTTPError as error:
                        if error.code in (400, 401, 403, 406, 429):
                            item["status"] = "pending"
                            save(checkpoint, state)
                        if error.code == 406:
                            # CAS occasionally returns a useful JSON validation
                            # message with 406. Show only allowlisted fields and
                            # redact credentials before reporting it.
                            try:
                                raw_detail = error.read(4096).decode("utf-8", "replace")
                                payload = json.loads(raw_detail)
                                detail = next((payload[key] for key in ("error", "message", "detail", "reason")
                                               if isinstance(payload, dict) and key in payload), None)
                                if detail is not None:
                                    detail = str(detail)
                                    for secret in (args.user, password):
                                        if secret:
                                            detail = detail.replace(secret, "[redacted]")
                                    print("CAS 406 detail: " + detail[:500], file=sys.stderr)
                                elif isinstance(payload, dict):
                                    print("CAS 406 JSON keys: " + ", ".join(sorted(payload.keys())), file=sys.stderr)
                            except json.JSONDecodeError:
                                detail = raw_detail.strip()
                                for secret in (args.user, password):
                                    if secret:
                                        detail = detail.replace(secret, "[redacted]")
                                if any(token in detail.upper() for token in ("SELECT", "FROM ", "LENS_ID")):
                                    detail = "[response omitted: request-related content]"
                                if detail:
                                    print("CAS 406 detail: " + detail[:300], file=sys.stderr)
                            except Exception:
                                pass
                        raise
                    item.update(status="submitted", job_id=str(submitted["id"]))
                    save(checkpoint, state)
                while True:
                    result = client.json("status", id=item["job_id"])
                    if result.get("status") == "done":
                        break
                    if result.get("status") == "error":
                        item.update(status="error", query_error=str(result.get("error", "CAS query failed"))[:1000])
                        save(checkpoint, state)
                        raise ValueError(f"CAS query failed for {name}")
                    print(f"{name}: {result.get('status', 'unknown')}", flush=True)
                    time.sleep(60)
                partial = results_dir / f"{name}.csv.part"
                if partial.exists():
                    partial.rename(results_dir / f"{name}.csv.incomplete-{time.time_ns()}")
                state["source_rows_download_started"] = True
                save(checkpoint, state)
                with client.post("download", id=item["job_id"]) as response, partial.open("xb") as stream:
                    for block in iter(lambda: response.read(1024 * 1024), b""):
                        if block and not state["source_rows_downloaded"]:
                            state["source_rows_downloaded"] = True
                            save(checkpoint, state)
                        stream.write(block)
                nrows = validate_rows(partial, expected)
                partial.replace(target)
                item.update(status="complete", rows=nrows, bytes=target.stat().st_size, sha256=sha(target))
                state["source_rows_downloaded"] = True
                save(checkpoint, state)
                print(f"{name}: verified {nrows} lens-source rows", flush=True)

            state.update(status="complete", total_pair_rows=manifest["lens_source_pair_rows_expected"],
                         source_rows_downloaded=True, science_measurement=False)
            save(checkpoint, state)
            (results_dir / "row_summary.json").write_text(json.dumps({
                "schema": "hsc-s19a-shape-row-summary-v1", "status": "complete",
                "n_lenses": 800, "n_pair_rows": manifest["lens_source_pair_rows_expected"],
                "columns": OUTPUT_COLUMNS, "source_rows_downloaded": True,
                "photoz_pdf_rows_included": False, "science_measurement": False,
                "count_summary_sum": count_summary["sum_lens_cone_shape_counts"],
            }, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
        except BaseException as error:
            state.update(status="interrupted" if isinstance(error, KeyboardInterrupt) else "stopped_on_error",
                         error_type=type(error).__name__)
            save(checkpoint, state)
            raise
        finally:
            client.credential.clear()
            password = None
    print(f"Complete: {manifest['lens_source_pair_rows_expected']} verified shape pair rows; no photo-z PDFs downloaded.")
    return 0


def query_plan_lenses(preflight, index):
    query_plan = load_plan(preflight)
    batch_size = int(query_plan["batch_size"])
    start = index * batch_size
    return query_plan["lenses"][start:start + batch_size]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--preflight", type=Path, default=DEFAULT_PREFLIGHT)
    parser.add_argument("--plan-dir", type=Path, default=DEFAULT_PLAN)
    parser.add_argument("--results", type=Path, default=DEFAULT_RESULTS)
    parser.add_argument("--user", default="kjhan0606")
    parser.add_argument("--plan", action="store_true", help="validate plan locally; no network")
    parser.add_argument("--submit", action="store_true", help="submit and download the prepared shape rows")
    parser.add_argument("--approved-source-row-export", action="store_true",
                        help="required only after explicit user approval for this exact data export")
    args = parser.parse_args()
    if args.plan and args.submit:
        parser.error("--plan and --submit are mutually exclusive")
    try:
        return run(args)
    except urllib.error.HTTPError as error:
        print(f"CAS HTTP {error.code}; response body and credentials suppressed.", file=sys.stderr)
        return 1
    except KeyboardInterrupt:
        print("Stopped; completed batches remain checkpointed and resumable.", file=sys.stderr)
        return 130
    except Exception as error:
        print("submit_wl_shape_rows: " + (str(error) if isinstance(error, ValueError) else type(error).__name__), file=sys.stderr)
        return 1


if __name__ == "__main__":
    sys.exit(main())
