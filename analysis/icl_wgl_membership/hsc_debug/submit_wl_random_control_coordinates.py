#!/usr/bin/env python3
"""Export only HSC random-point coordinates for the approved control pool.

This stage returns lens_id, ra_deg, and dec_deg for the same candidates counted
by submit_wl_random_control_counts.py. It returns no source catalog rows,
shapes, photo-z values, or images. Results and credentials stay outside Git.
"""
from __future__ import annotations

import argparse
import csv
import getpass
import hashlib
import json
import math
from pathlib import Path
import platform
import sys
import time
import urllib.error
from datetime import datetime, timezone

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from fetch_ancillary_catalogues import Cas, catalogue_reader, save, sha
from prepare_wl_random_control_counts import random_count_query
from submit_wl_random_control_counts import (
    DEFAULT_PLAN as DEFAULT_COUNT_PLAN,
    load_random_plan,
    read_counts,
    safe_error_detail,
)

DEFAULT_PLAN = Path(
    "/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/"
    "wl_random_point_coordinate_plan_2026-10-19_01"
)


def coordinate_query(lens, *, inner_radius_arcsec=None,
                     outer_radius_arcsec=1800.0):
    count_sql = random_count_query(
        lens, inner_radius_arcsec=inner_radius_arcsec,
        outer_radius_arcsec=outer_radius_arcsec,
    )
    from_at = count_sql.index("FROM pdr3_wide.random AS r")
    where_at = count_sql.index("\nWHERE ", from_at)
    where_sql = count_sql[where_at:]
    lens_id = str(lens["lens_id"]).replace("'", "''")
    return (
        f"SELECT '{lens_id}' AS lens_id, r.ra AS ra_deg, r.dec AS dec_deg\n"
        "FROM pdr3_wide.random AS r"
        + where_sql
    )


def load_count_inputs(count_plan: Path):
    manifest, lens_ids = load_random_plan(count_plan)
    result_dir = count_plan / "cas_counts"
    checkpoint_path = result_dir / "count_status.json"
    summary_path = result_dir / "count_summary.json"
    checkpoint = json.loads(checkpoint_path.read_text(encoding="utf-8"))
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    plan_sha = sha(count_plan / "random_control_count_plan.json")
    if checkpoint.get("status") != "complete" or checkpoint.get("plan_sha256") != plan_sha:
        raise ValueError("Count-only checkpoint is incomplete or belongs to another plan")
    if summary.get("status") != "complete" or summary.get("query_plan_sha256") != plan_sha:
        raise ValueError("Count-only summary is incomplete or belongs to another plan")
    if summary.get("n_lenses") != 800 or summary.get("n_count_rows") != 800:
        raise ValueError("Count-only summary does not cover exactly 800 lenses")
    counts = {}
    for index, item in enumerate(manifest["batches"]):
        name = f"random_counts_{index:03d}"
        job = checkpoint.get("jobs", {}).get(name, {})
        csv_path = result_dir / f"{name}.csv"
        if job.get("status") != "complete" or sha(csv_path) != job.get("sha256"):
            raise ValueError(f"Count result is missing or changed: {name}")
        start = index * int(manifest["batch_size"])
        ids = lens_ids[start:start + int(manifest["batch_size"])]
        for row in read_counts(csv_path, ids):
            counts[row["lens_id"]] = row["n_random_candidates"]
    if set(counts) != set(lens_ids):
        raise ValueError("Count results do not cover the approved 800-lens input")
    if sum(counts.values()) != summary.get("sum_candidate_counts"):
        raise ValueError("Count results disagree with the aggregate summary")
    lens_manifest = json.loads((Path(manifest["source_preflight"]) / "preflight_manifest.json").read_text(encoding="utf-8"))
    lenses = lens_manifest["lenses"]
    if [str(row["lens_id"]) for row in lenses] != lens_ids:
        raise ValueError("Preflight lens order changed after count query")
    return manifest, lenses, lens_ids, counts, plan_sha


def prepare(count_plan: Path, output: Path):
    if output.exists():
        raise FileExistsError(f"refusing to overwrite existing coordinate plan: {output}")
    count_manifest, lenses, lens_ids, counts, count_plan_sha = load_count_inputs(count_plan)
    batch_size = int(count_manifest["batch_size"])
    output.mkdir(parents=True, exist_ok=False)
    batches = []
    for index, start in enumerate(range(0, len(lenses), batch_size)):
        chunk = lenses[start:start + batch_size]
        candidate_query = count_manifest.get("candidate_query", {})
        sql = "\nUNION ALL\n".join(
            coordinate_query(
                lens,
                inner_radius_arcsec=candidate_query.get("inner_radius_arcsec"),
                outer_radius_arcsec=candidate_query.get("outer_radius_arcsec", 1800.0),
            )
            for lens in chunk
        ) + "\n"
        path = output / f"random_coordinates_{index:03d}.sql"
        path.write_text(sql, encoding="utf-8")
        batches.append({
            "path": str(path),
            "n_lenses": len(chunk),
            "expected_rows": sum(counts[str(lens["lens_id"])] for lens in chunk),
            "first_lens_id": str(chunk[0]["lens_id"]),
            "last_lens_id": str(chunk[-1]["lens_id"]),
            "sha256": sha(path),
        })
    manifest = {
        "schema": "hsc-pdr3-random-control-coordinate-plan-v1",
        "status": "prepared_not_submitted",
        "release": "pdr3",
        "count_plan": str(count_plan.resolve()),
        "count_plan_sha256": count_plan_sha,
        "n_lenses": len(lenses),
        "n_query_batches": len(batches),
        "batch_size": batch_size,
        "candidate_pair_rows_expected": sum(counts.values()),
        "output_columns": ["lens_id", "ra_deg", "dec_deg"],
        "random_points_only": True,
        "source_rows_returned": False,
        "shape_rows_returned": False,
        "photoz_rows_returned": False,
        "images_returned": False,
        "science_measurement": False,
        "post_export_gates": [
            "exclude candidates within twice the 3 Mpc aperture of any of the 800 CAMIRA clusters",
            "check full-aperture HSC footprint and i-band mask coverage",
            "match local HSC-I depth, PSF, sky, and extinction to cluster lenses",
            "deduplicate sky positions and audit overlapping control apertures",
            "measure spatial-block covariance before fixing the final control count",
        ],
        "batches": batches,
    }
    (output / "random_control_coordinate_plan.json").write_text(
        json.dumps(manifest, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    return manifest


def validate_coordinate_plan(plan_dir: Path):
    path = plan_dir / "random_control_coordinate_plan.json"
    manifest = json.loads(path.read_text(encoding="utf-8"))
    if manifest.get("schema") != "hsc-pdr3-random-control-coordinate-plan-v1":
        raise ValueError("Unexpected random-coordinate plan schema")
    if manifest.get("status") != "prepared_not_submitted":
        raise ValueError("Random-coordinate plan is not in prepared state")
    if manifest.get("output_columns") != ["lens_id", "ra_deg", "dec_deg"]:
        raise ValueError("Unexpected coordinate output columns")
    if any(manifest.get(key) is not False for key in (
        "source_rows_returned", "shape_rows_returned", "photoz_rows_returned", "images_returned", "science_measurement"
    )):
        raise ValueError("Plan requests data outside random-point coordinates")
    count_plan = Path(manifest["count_plan"])
    count_manifest, lenses, lens_ids, counts, count_plan_sha = load_count_inputs(count_plan)
    if count_plan_sha != manifest.get("count_plan_sha256"):
        raise ValueError("Count-only input plan changed")
    if manifest.get("candidate_pair_rows_expected") != sum(counts.values()):
        raise ValueError("Expected coordinate-row count changed")
    batches = manifest.get("batches", [])
    if len(batches) != 32 or len(lenses) != 800 or len(lens_ids) != 800:
        raise ValueError("Expected 800 lenses in 32 coordinate batches")
    for index, batch in enumerate(batches):
        sql_path = Path(batch["path"])
        if sql_path.resolve().parent != plan_dir.resolve() or sql_path.name != f"random_coordinates_{index:03d}.sql":
            raise ValueError("SQL path is outside the coordinate plan directory")
        if sha(sql_path) != batch.get("sha256"):
            raise ValueError(f"Coordinate SQL checksum changed: {sql_path.name}")
        sql = sql_path.read_text(encoding="utf-8").upper()
        if "FROM PDR3_WIDE.RANDOM AS R" not in sql or "R.RA AS RA_DEG" not in sql or "R.DEC AS DEC_DEG" not in sql:
            raise ValueError(f"Unexpected coordinate source or selection: {sql_path.name}")
        if "R.ISPRIMARY IS TRUE" not in sql or "R.I_PIXELFLAGS_EDGE IS FALSE" not in sql:
            raise ValueError(f"Required random-point quality filters missing: {sql_path.name}")
        if "PHOTOZ" in sql or "HSM" in sql or "SHAPE" in sql or "SELECT *" in sql:
            raise ValueError(f"Batch requests unrelated catalog fields: {sql_path.name}")
        if "COUNT(" in sql or "CANDIDATE" in sql:
            raise ValueError(f"Batch is not a coordinate-row query: {sql_path.name}")
        start = index * int(manifest["batch_size"])
        ids = lens_ids[start:start + int(manifest["batch_size"])]
        for lens_id in ids:
            if f"SELECT '{lens_id.upper()}' AS LENS_ID" not in sql:
                raise ValueError(f"Expected lens absent from {sql_path.name}")
    return manifest, batches, lenses, lens_ids, counts


def read_coordinates(path: Path, expected_counts):
    expected_ids = list(expected_counts)
    expected_rows = sum(expected_counts.values())
    per_lens = {lens_id: 0 for lens_id in expected_ids}
    rows = 0
    with path.open(encoding="utf-8-sig", newline="") as stream:
        reader = catalogue_reader(stream)
        if reader.fieldnames != ["lens_id", "ra_deg", "dec_deg"]:
            raise ValueError("CAS returned columns outside the approved coordinate schema")
        for row in reader:
            lens_id = str(row["lens_id"])
            if lens_id not in per_lens:
                raise ValueError("CAS returned an unrequested lens ID")
            ra, dec = float(row["ra_deg"]), float(row["dec_deg"])
            if not math.isfinite(ra) or not 0 <= ra < 360:
                raise ValueError("CAS returned invalid RA")
            if not math.isfinite(dec) or not -90 <= dec <= 90:
                raise ValueError("CAS returned invalid Dec")
            per_lens[lens_id] += 1
            rows += 1
    if per_lens != expected_counts:
        raise ValueError("Coordinate rows do not match the approved per-lens candidate counts")
    if rows != expected_rows:
        raise ValueError(f"Expected {expected_rows} lens-candidate rows; got {rows}")
    return rows


def cas_status(client, job_id: str):
    for attempt in range(10):
        try:
            return client.json("status", id=job_id)
        except urllib.error.HTTPError as error:
            if error.code not in (502, 503, 504) or attempt == 9:
                raise
            print(f"CAS status gateway HTTP {error.code}; retrying in 60 seconds "
                  f"({attempt + 1}/10)", flush=True)
            time.sleep(60)


def record_status_poll(state, item, status, timestamp=None):
    """Persist every successful status response so long polls remain auditable."""
    polled_at = timestamp or datetime.now(timezone.utc).isoformat()
    item["last_cas_status"] = str(status.get("status", "unknown"))
    item["last_status_poll_utc"] = polled_at
    item["status_poll_count"] = int(item.get("status_poll_count", 0)) + 1
    state.pop("error_type", None)


def fetch_batch(client, batch, expected_ids, expected_counts, results, state, checkpoint):
    name = Path(batch["path"]).stem
    item = state["jobs"].setdefault(name, {"status": "pending"})
    target = results / f"{name}.csv"
    if item.get("status") == "complete":
        if not target.is_file() or sha(target) != item.get("sha256"):
            raise ValueError(f"Saved coordinate result missing or changed: {name}")
        read_coordinates(target, expected_counts)
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
            if error.code == 406:
                print("CAS 406 detail: " + safe_error_detail(error, client.credential["password"]), file=sys.stderr)
            raise
        item.update(status="submitted", job_id=str(submitted["id"]))
        save(checkpoint, state)
    while True:
        status = cas_status(client, item["job_id"])
        record_status_poll(state, item, status)
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
    rows = read_coordinates(partial, expected_counts)
    partial.replace(target)
    item.update(status="complete", rows=rows, bytes=target.stat().st_size, sha256=sha(target))
    save(checkpoint, state)
    print(f"{name}: verified {rows} coordinate rows", flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--count-plan", type=Path, default=DEFAULT_COUNT_PLAN)
    parser.add_argument("--plan", type=Path, default=DEFAULT_PLAN)
    parser.add_argument("--prepare", action="store_true", help="prepare a coordinate-only plan; no network")
    parser.add_argument("--check", action="store_true", help="validate the prepared plan; no network")
    parser.add_argument("--submit", action="store_true", help="submit the approved coordinate batches")
    args = parser.parse_args()
    if sum((args.prepare, args.check, args.submit)) != 1:
        parser.error("select exactly one of --prepare, --check, or --submit")
    if args.prepare:
        manifest = prepare(args.count_plan.resolve(), args.plan.resolve())
        print(json.dumps({"status": manifest["status"], "lenses": manifest["n_lenses"],
            "expected_rows": manifest["candidate_pair_rows_expected"], "batches": manifest["n_query_batches"],
            "columns": manifest["output_columns"], "plan": str(args.plan.resolve())}, indent=2))
        return 0
    manifest, batches, lenses, lens_ids, counts = validate_coordinate_plan(args.plan.resolve())
    results = args.plan.resolve() / "cas_coordinates"
    print(json.dumps({"plan_status": manifest["status"], "lenses": len(lens_ids),
        "batches": len(batches), "expected_rows": manifest["candidate_pair_rows_expected"],
        "columns": manifest["output_columns"], "source_rows": False, "shapes": False,
        "photoz": False, "images": False, "results": str(results)}, indent=2), flush=True)
    if args.check:
        return 0
    if not sys.stdin.isatty():
        parser.error("--submit requires a terminal for one hidden password prompt")
    if not platform.node().lower().startswith("lageunha"):
        parser.error("Run --submit from the established lageunha interactive environment")
    checkpoint = results / "coordinate_status.json"
    plan_sha = sha(args.plan / "random_control_coordinate_plan.json")
    if results.exists():
        if not checkpoint.is_file():
            parser.error("Refusing to reuse a coordinate result directory without checkpoint")
        state = json.loads(checkpoint.read_text(encoding="utf-8"))
        if state.get("plan_sha256") != plan_sha:
            parser.error("Existing coordinate checkpoint belongs to a different plan")
        if any(item.get("status") in ("submitting", "error") for item in state.get("jobs", {}).values()):
            parser.error("Submission outcome is ambiguous or failed; inspect CAS before retry")
    else:
        results.mkdir(parents=True, exist_ok=False)
        state = {"schema": "hsc-pdr3-random-control-coordinate-run-v1", "status": "pending",
            "plan_sha256": plan_sha, "jobs": {}, "credentials_saved": False,
            "source_rows_returned": False, "shape_rows_returned": False,
            "photoz_rows_returned": False, "images_returned": False, "science_measurement": False}
    import fcntl
    with (results / "run.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        password = getpass.getpass("HSC password (hidden): ")
        if not password:
            raise ValueError("Empty password")
        client = Cas("kjhan0606", password)
        state["status"] = "running"
        state.pop("error_type", None)
        save(checkpoint, state)
        try:
            size = int(manifest["batch_size"])
            for index, batch in enumerate(batches):
                start = index * size
                expected_ids = lens_ids[start:start + size]
                expected_counts = {key: counts[key] for key in expected_ids}
                if sum(expected_counts.values()) != batch["expected_rows"]:
                    raise ValueError(f"Plan count mismatch for {Path(batch['path']).name}")
                fetch_batch(client, batch, expected_ids, expected_counts, results, state, checkpoint)
            total_rows = sum(item["rows"] for item in state["jobs"].values())
            state.update(status="complete", total_coordinate_rows=total_rows,
                no_source_rows=True, no_shape_rows=True, no_photoz_rows=True, no_images=True,
                science_measurement=False)
            save(checkpoint, state)
            summary = {"schema": "hsc-pdr3-random-control-coordinate-summary-v1", "status": "complete",
                "n_lenses": 800, "n_coordinate_rows": total_rows,
                "output_columns": manifest["output_columns"], "source_rows_returned": False,
                "shape_rows_returned": False, "photoz_rows_returned": False,
                "images_returned": False, "science_measurement": False,
                "coordinate_plan_sha256": plan_sha}
            (results / "coordinate_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
        except BaseException as error:
            state.update(status="interrupted" if isinstance(error, KeyboardInterrupt) else "stopped_on_error",
                error_type=type(error).__name__)
            save(checkpoint, state)
            raise
        finally:
            client.credential.clear()
            password = None
    print(f"Complete: {total_rows} coordinate rows; no source data. Results: {results}")
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
