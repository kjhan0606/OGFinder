#!/usr/bin/env python3
"""Submit the approved HSC PDR3 patch-coverage metadata plan, one batch at a time.

Only the prepared mosaic/skymap/patch_qa query is accepted. The tool returns
patch coverage metadata; it never requests images, object rows, shapes, or photo-z.
"""
from __future__ import annotations

import argparse
import csv
import getpass
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

DEFAULT_PLAN = Path(
    "/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/"
    "icl_expansion_coverage_plan_2026-10-21_01"
)
OUTPUT_COLUMNS = [
    "candidate_id", "probe_id", "skymap_id", "tract", "patch", "band",
    "naxis1", "naxis2", "seeing", "zeropt", "imag_psf_depth", "iseeing",
    "imeas_exists", "idet_bkgd_exists",
]
FORBIDDEN_SQL = (
    "SELECT *", "PHOTOZ", "SHAPE", "HSM", "OBJECT_ID", "SOURCE_ID",
    "PDR3_WIDE.FORCED", "PDR3_WIDE.MERGED", "PDR3_WIDE.RAW",
    "PDR3_WIDE.RANDOM", "PDR3_WIDE.COADD", "FITS", "IMAGE",
)
EXPECTED_S23B_PLAN_SHA256 = "e615ee7c2414109af6d3c156104d71766355d3b4e2a5dba3636a0142995557d7"


def validate_result(path: Path, batch: dict, allowed_probes: set[tuple[str, str]]) -> int:
    rows = 0
    seen = set()
    batch_candidates = set()
    with path.open(encoding="utf-8-sig", newline="") as stream:
        reader = catalogue_reader(stream)
        if reader.fieldnames != OUTPUT_COLUMNS:
            raise ValueError(f"CAS returned columns outside the approved schema: {path.name}")
        for row in reader:
            if None in row or any(value is None for value in row.values()):
                raise ValueError(f"Malformed metadata row in {path.name}")
            candidate_id, probe_id = str(row["candidate_id"]), str(row["probe_id"])
            if (candidate_id, probe_id) not in allowed_probes:
                raise ValueError(f"CAS returned an unrequested probe in {path.name}")
            key = (candidate_id, probe_id, str(row["skymap_id"]))
            if key in seen:
                raise ValueError(f"Duplicate patch metadata row in {path.name}")
            seen.add(key)
            batch_candidates.add(candidate_id)
            if row["band"] != "HSC-I":
                raise ValueError(f"CAS returned a non-HSC-I mosaic in {path.name}")
            if not str(row["tract"]).isdigit() or not str(row["patch"]):
                raise ValueError(f"Invalid tract/patch metadata in {path.name}")
            for key_name in ("naxis1", "naxis2"):
                if int(row[key_name]) <= 0:
                    raise ValueError(f"Invalid mosaic dimensions in {path.name}")
            for key_name in ("seeing", "zeropt"):
                if float(row[key_name]) <= 0:
                    raise ValueError(f"Invalid {key_name} metadata in {path.name}")
            for key_name in ("imag_psf_depth", "iseeing"):
                if row[key_name] not in ("", "NULL", r"\N") and float(row[key_name]) <= 0:
                    raise ValueError(f"Invalid {key_name} QA metadata in {path.name}")
            for key_name in ("imeas_exists", "idet_bkgd_exists"):
                if str(row[key_name]).lower() not in {"true", "false", "t", "f", "1", "0"}:
                    raise ValueError(f"Invalid {key_name} flag in {path.name}")
            rows += 1
    if not batch_candidates.issubset({key[0] for key in allowed_probes}):
        raise ValueError(f"Unexpected candidate IDs in {path.name}")
    if rows and rows > 4 * batch["n_probe_points"]:
        raise ValueError(f"Unexpectedly many patch matches in {path.name}")
    return rows


def validate_plan(plan_dir: Path):
    plan_path = plan_dir / "coverage_query_plan.json"
    manifest = json.loads(plan_path.read_text(encoding="utf-8"))
    schema = manifest.get("schema")
    if schema not in {
        "astrafek-icl-pdr3-coverage-query-plan-v1",
        "astrafek-s23b-pdr3-coverage-query-plan-v1",
    }:
        raise ValueError("Unexpected coverage-plan schema")
    plan_sha = sha(plan_path)
    if schema == "astrafek-icl-pdr3-coverage-query-plan-v1":
        if manifest.get("status") != "prepared_not_submitted" or manifest.get("cas_query_submitted") is not False:
            raise ValueError("Coverage plan is not in reviewed prepared state")
        if any(manifest.get(key) is not False for key in (
            "science_measurement", "images_downloaded", "shapes_downloaded", "photoz_rows_requested"
        )):
            raise ValueError("Coverage plan flags data outside the approved scope")
        expected_batches, expected_candidates, expected_probes = 52, 2596, 23364
        expected_candidate_prefix = None
    else:
        if (plan_sha != EXPECTED_S23B_PLAN_SHA256
                or manifest.get("status") != "prepared_pending_user_approval"
                or manifest.get("user_approval_required_before_submission") is not True
                or manifest.get("cas_query_submitted") is not False):
            raise ValueError("S23B plan differs from the exact scope pending user approval")
        if any(manifest.get(key) is not False for key in (
            "science_measurement", "data_downloaded", "images_downloaded",
            "shape_or_photoz_rows_requested",
        )):
            raise ValueError("S23B plan flags data outside its metadata-only scope")
        expected_batches, expected_candidates, expected_probes = 27, 1333, 11997
        expected_candidate_prefix = "S23B:"
        prior = manifest.get("prior_approved_plan", {})
        if (prior.get("sha256") != "9da220c05f1f24c409a602fcc4ddffc8cbd5bc33f459d335cc8587f1cd235006"
                or prior.get("cas_query_submitted") is not False):
            raise ValueError("S23B plan no longer references the separately prepared S20A plan")
    coverage = manifest.get("coverage_test", {})
    if coverage.get("output_columns") != OUTPUT_COLUMNS:
        raise ValueError("Coverage plan output columns differ from approved schema")
    if coverage.get("images_or_catalog_object_rows") is not False:
        raise ValueError("Coverage plan includes image/object data")
    for raw_path, expected_hash in manifest.get("output_sha256", {}).items():
        path = Path(raw_path)
        if path.resolve().parent != plan_dir.resolve() or not path.is_file() or sha(path) != expected_hash:
            raise ValueError(f"Prepared plan artifact is missing or changed: {path.name}")
    batches = manifest.get("batches", [])
    if (len(batches) != expected_batches
            or manifest.get("n_candidates") != expected_candidates
            or manifest.get("n_probes") != expected_probes):
        raise ValueError("Coverage plan dimensions differ from the approved request")
    centers_path = plan_dir / "candidate_centers.tsv"
    expected_center_key = "id" if expected_candidate_prefix is None else "candidate_id"
    with centers_path.open(encoding="utf-8", newline="") as stream:
        centers = list(csv.DictReader(stream, delimiter="\t"))
    center_ids = [str(row[expected_center_key]) for row in centers]
    if len(center_ids) != expected_candidates or len(set(center_ids)) != expected_candidates:
        raise ValueError("Candidate-center inventory is incomplete or duplicated")
    if expected_candidate_prefix is not None and any(not ident.startswith(expected_candidate_prefix) for ident in center_ids):
        raise ValueError("S23B candidate IDs are not release-prefixed")
    expected_probe_path = plan_dir / "coverage_probes.tsv"
    allowed_probes = set()
    with expected_probe_path.open(encoding="utf-8", newline="") as stream:
        for row in csv.DictReader(stream, delimiter="\t"):
            allowed_probes.add((str(row["candidate_id"]), str(row["probe_id"])))
    if len(allowed_probes) != expected_probes:
        raise ValueError("Coverage probe inventory is incomplete or duplicated")
    if {candidate_id for candidate_id, _ in allowed_probes} != set(center_ids):
        raise ValueError("Probe and candidate-center inventories disagree")
    batch_candidate_total = batch_probe_total = 0
    for index, batch in enumerate(batches):
        sql_path = Path(batch["path"])
        if sql_path.resolve().parent != plan_dir.resolve() or sql_path.name != f"coverage_mosaic_{index:03d}.sql":
            raise ValueError("SQL path is outside the approved plan directory")
        if sha(sql_path) != batch.get("sha256"):
            raise ValueError(f"Approved SQL checksum changed: {sql_path.name}")
        sql = sql_path.read_text(encoding="utf-8").upper()
        if any(token in sql for token in FORBIDDEN_SQL):
            raise ValueError(f"SQL requests fields/data outside metadata scope: {sql_path.name}")
        required = (
            "PDR3_WIDE.MOSAIC", "PUBLIC.SKYMAP", "PDR3_WIDE.PATCH_QA",
            "PATCH_CONTAINS", "M.BAND = 'HSC-I'",
        )
        if any(token not in sql for token in required):
            raise ValueError(f"Unexpected source, geometry, or band in {sql_path.name}")
        if sql.count("SELECT P.CANDIDATE_ID") != 1 or sql.count("FROM PROBES AS P") != 1:
            raise ValueError(f"Unexpected query shape in {sql_path.name}")
        if expected_candidate_prefix is not None and "S23B:" not in sql:
            raise ValueError(f"S23B batch omits release-prefixed IDs: {sql_path.name}")
        if batch.get("expected_rows") is not None:
            raise ValueError("Coverage plan unexpectedly pins a result-row count")
        batch_candidate_total += int(batch.get("n_candidates", 0))
        batch_probe_total += int(batch.get("n_probe_points", 0))
        if expected_candidate_prefix is not None:
            batch_probe_pairs = {
                pair for pair in allowed_probes if f"('{pair[0]}', '{pair[1]}',".upper() in sql
            }
            if len(batch_probe_pairs) != int(batch.get("n_probe_points", -1)):
                raise ValueError(f"S23B SQL/probe inventory mismatch: {sql_path.name}")
    if (len(allowed_probes) != manifest["n_probes"]
            or batch_candidate_total != expected_candidates
            or batch_probe_total != expected_probes):
        raise ValueError("Probe count disagrees with the coverage plan")
    return manifest, batches, allowed_probes, plan_sha


def check_submission_authorization(schema: str, plan_sha: str,
                                   approved_plan_sha256: str | None) -> None:
    """Require a fresh, exact-hash user approval for the separately scoped S23B query."""
    if schema == "astrafek-s23b-pdr3-coverage-query-plan-v1":
        if approved_plan_sha256 != plan_sha:
            raise ValueError("S23B submission requires approval of this exact plan SHA-256")
    elif approved_plan_sha256 is not None and approved_plan_sha256 != plan_sha:
        raise ValueError("approval hash does not match the selected metadata plan")


def status_with_retry(client, job_id):
    for attempt in range(10):
        try:
            return client.json("status", id=job_id)
        except urllib.error.HTTPError as error:
            if error.code not in (502, 503, 504) or attempt == 9:
                raise
            print(f"CAS status gateway HTTP {error.code}; retrying in 60 seconds ({attempt + 1}/10)", flush=True)
            time.sleep(60)


def fetch_batch(client, batch, allowed_probes, results, state, checkpoint):
    name = Path(batch["path"]).stem
    item = state["jobs"].setdefault(name, {"status": "pending"})
    target = results / f"{name}.csv"
    if item.get("status") == "complete":
        if not target.is_file() or sha(target) != item.get("sha256"):
            raise ValueError(f"Saved result missing or changed: {name}")
        validate_result(target, batch, allowed_probes)
        return
    if item.get("status") in ("submitting", "error"):
        raise ValueError(f"{name} has an ambiguous or failed submission; inspect CAS before retry")
    if "job_id" not in item:
        item["status"] = "submitting"
        save(checkpoint, state)
        try:
            submitted = client.json(
                "submit", catalog_job=dict(sql=Path(batch["path"]).read_text(encoding="utf-8"),
                    out_format="csv", include_metainfo_to_body=False, release_version="pdr3"),
                nomail=True, skip_syntax_check=False,
            )
        except urllib.error.HTTPError as error:
            if error.code in (400, 401, 403, 406, 429):
                item["status"] = "pending"
                save(checkpoint, state)
            raise
        item.update(status="submitted", job_id=str(submitted["id"]))
        save(checkpoint, state)
    while True:
        status = status_with_retry(client, item["job_id"])
        item["last_cas_status"] = str(status.get("status", "unknown"))
        item["last_status_poll_utc"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
        item["status_poll_count"] = int(item.get("status_poll_count", 0)) + 1
        save(checkpoint, state)
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
    rows = validate_result(partial, batch, allowed_probes)
    partial.replace(target)
    item.update(status="complete", rows=rows, bytes=target.stat().st_size, sha256=sha(target))
    save(checkpoint, state)
    print(f"{name}: verified {rows} coverage metadata rows", flush=True)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan", type=Path, default=DEFAULT_PLAN)
    parser.add_argument("--check", action="store_true", help="verify the approved offline plan only")
    parser.add_argument("--submit", action="store_true", help="submit the approved metadata batches")
    parser.add_argument("--approved-plan-sha256", help="required exact plan hash for newly scoped S23B submissions")
    args = parser.parse_args(argv)
    if args.check == args.submit:
        parser.error("select exactly one of --check or --submit")
    manifest, batches, allowed_probes, plan_sha = validate_plan(args.plan.resolve())
    if args.submit or args.approved_plan_sha256 is not None:
        try:
            check_submission_authorization(
                manifest.get("schema", ""), plan_sha, args.approved_plan_sha256)
        except ValueError as error:
            parser.error(str(error))
    results = args.plan.resolve() / "cas_coverage_metadata"
    print(json.dumps({"plan_status": manifest["status"], "candidates": manifest["n_candidates"],
        "probes": manifest["n_probes"], "batches": len(batches), "columns": OUTPUT_COLUMNS,
        "images": False, "object_rows": False, "shapes": False, "photoz_rows": False,
        "science_measurement": False, "results": str(results)}, indent=2), flush=True)
    if args.check:
        return 0
    if not sys.stdin.isatty():
        parser.error("--submit requires a terminal for the hidden HSC password prompt")
    if not platform.node().lower().startswith("lageunha"):
        parser.error("Run --submit from the established lageunha interactive environment")
    checkpoint = results / "coverage_status.json"
    if results.exists():
        if not checkpoint.is_file():
            parser.error("Refusing to reuse a result directory without its checkpoint")
        state = json.loads(checkpoint.read_text(encoding="utf-8"))
        if state.get("plan_sha256") != plan_sha:
            parser.error("Existing checkpoint belongs to a different metadata plan")
        if any(item.get("status") in ("submitting", "error") for item in state.get("jobs", {}).values()):
            parser.error("Submission outcome is ambiguous or failed; inspect CAS before retry")
    else:
        results.mkdir(parents=True, exist_ok=False)
        state = {"schema": "astrafek-icl-pdr3-coverage-metadata-run-v1", "status": "pending",
            "plan_sha256": plan_sha, "jobs": {}, "credentials_saved": False,
            "images_returned": False, "object_rows_returned": False, "shape_rows_returned": False,
            "photoz_rows_returned": False, "science_measurement": False}
    import fcntl
    with (results / "run.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        password = getpass.getpass("HSC password (hidden): ")
        if not password:
            raise ValueError("Empty password")
        client = Cas("kjhan0606", password)
        state["status"] = "running"
        save(checkpoint, state)
        try:
            for batch in batches:
                fetch_batch(client, batch, allowed_probes, results, state, checkpoint)
            total_rows = sum(item["rows"] for item in state["jobs"].values())
            state.update(status="complete", total_metadata_rows=total_rows,
                images_returned=False, object_rows_returned=False, shape_rows_returned=False,
                photoz_rows_returned=False, science_measurement=False)
            save(checkpoint, state)
            summary = {"schema": "astrafek-icl-pdr3-coverage-metadata-summary-v1", "status": "complete",
                "n_candidates": manifest["n_candidates"], "n_probes": manifest["n_probes"],
                "n_metadata_rows": total_rows, "output_columns": OUTPUT_COLUMNS,
                "images_returned": False, "object_rows_returned": False,
                "shape_rows_returned": False, "photoz_rows_returned": False,
                "science_measurement": False, "plan_sha256": plan_sha}
            (results / "coverage_summary.json").write_text(
                json.dumps(summary, indent=2) + "\n", encoding="utf-8")
        except BaseException as error:
            state.update(status="interrupted" if isinstance(error, KeyboardInterrupt) else "stopped_on_error",
                error_type=type(error).__name__)
            save(checkpoint, state)
            raise
        finally:
            client.credential.clear()
            password = None
    print(f"Complete: {total_rows} coverage metadata rows; no images or source data. Results: {results}")
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except urllib.error.HTTPError as error:
        print(f"HSC CAS HTTP {error.code}; credentials and response body suppressed.", file=sys.stderr)
        raise SystemExit(1)
