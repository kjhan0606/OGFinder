#!/usr/bin/env python3
"""Validate or submit the pinned S23B HSC random-control candidate plan.

Submission is gated on an explicit approval of the exact plan SHA-256. The
export is limited to random-position identifiers/coordinates and scalar
I-band/patch QA fields. It never requests images, detected-source rows,
shapes, or photo-z values, and it does not save credentials.
"""
from __future__ import annotations

import argparse
import csv
import fcntl
import getpass
import hashlib
import json
import math
from pathlib import Path
import platform
import sys
import time
import urllib.error

from fetch_ancillary_catalogues import Cas, catalogue_reader, save


DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
DEFAULT_PLAN = DATA / "s23b_icl_control_candidate_plan_2026-11-06_01"
DEFAULT_RESULTS = DATA / "s23b_icl_control_candidate_results_2026-11-06_01"
DEFAULT_USER = "kjhan0606"
EXPECTED_PLAN_SHA256 = "143a1bd0a3a3c341576ce420f77ba0089fbfe80fdf58fbfd7afe30ea8e3f6714"
MAX_BATCH_BYTES = 20 * 1024 * 1024
POLL_SECONDS = 20
OUTPUT_COLUMNS = [
    "candidate_id", "random_object_id", "ra", "dec", "skymap_id", "tract", "patch", "patch_s",
    "isprimary", "detect_ispatchinner", "detect_istractinner", "i_pix_variance", "i_sky_mean",
    "i_sky_std", "i_inputcount_value", "i_inputcount_flag", "i_inputcount_flag_noinputs",
    "i_inputcount_flag_badcentroid", "i_pixelflags_offimage", "i_pixelflags_edge",
    "i_pixelflags_bad", "i_pixelflags_interpolatedcenter", "i_pixelflags_saturatedcenter",
    "i_pixelflags_crcenter", "i_pixelflags_suspectcenter", "i_pixelflags_bright_objectcenter",
    "icalexp_exists", "imag_psf_depth", "iseeing",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def parse_float(value: str, field: str, *, nullable: bool = False) -> float | None:
    text = str(value).strip()
    if nullable and text.lower() in ("", "null", "none", "nan", r"\n", r"\N"):
        return None
    try:
        result = float(text)
    except (TypeError, ValueError):
        raise ValueError(f"CAS returned a nonnumeric {field}") from None
    if not math.isfinite(result):
        raise ValueError(f"CAS returned a non-finite {field}")
    return result


def load_plan(plan_dir: Path) -> tuple[dict, list[dict[str, str]], list[dict], str]:
    plan_dir = Path(plan_dir).resolve()
    plan_path = plan_dir / "control_candidate_query_plan.json"
    plan = json.loads(plan_path.read_text(encoding="utf-8"))
    plan_hash = sha256(plan_path)
    if plan_hash != EXPECTED_PLAN_SHA256:
        raise ValueError("query plan differs from the pinned reviewed scope")
    if (plan.get("schema") != "astrafek-s23b-icl-control-candidate-query-plan-v1"
            or plan.get("status") != "prepared_pending_user_approval"
            or plan.get("release") != "pdr3"
            or plan.get("table") != "pdr3_wide.random plus pdr3_wide.patch_qa"
            or plan.get("n_target_candidates") != 1043
            or plan.get("max_random_rows_per_candidate") != 200
            or plan.get("maximum_returned_rows") != 208600
            or plan.get("search_radius_arcsec") != 1800.0
            or plan.get("n_sql_batches") != 105
            or plan.get("estimated_csv_bytes") != 80_000_000
            or plan.get("output_columns") != OUTPUT_COLUMNS):
        raise ValueError("candidate count, radius, row limit, table, or output schema changed")
    if (plan.get("shape_or_photoz_rows_requested") is not False
            or plan.get("images_requested") is not False
            or plan.get("source_detections_requested") is not False
            or plan.get("credentials_saved") is not False
            or plan.get("cas_query_submitted") is not False
            or plan.get("science_measurement") is not False
            or plan.get("user_approval_required_before_submission") is not True):
        raise ValueError("plan provenance or export scope flags changed")
    for text, expected in plan.get("input_sha256", {}).items():
        path = Path(text)
        if not path.is_file() or sha256(path) != expected:
            raise ValueError(f"plan input changed or is missing: {path}")

    candidate_path = Path(plan["candidate_table"])
    if sha256(candidate_path) != plan.get("candidate_table_sha256"):
        raise ValueError("eligible S23B candidate ledger hash changed")
    with candidate_path.open(encoding="utf-8", newline="") as stream:
        candidates = list(csv.DictReader(stream, delimiter="\t"))
    candidate_ids = [row["candidate_id"] for row in candidates]
    if len(candidates) != 1043 or len(set(candidate_ids)) != 1043:
        raise ValueError("eligible candidate ledger is incomplete or duplicated")

    batches = plan.get("batches", [])
    if len(batches) != 105 or sum(int(row["n_candidates"]) for row in batches) != 1043:
        raise ValueError("batch ledger does not cover all 1,043 eligible candidates")
    cursor = 0
    for index, batch in enumerate(batches):
        name = f"random_control_candidates_{index:03d}"
        sql_path = (plan_dir / f"{name}.sql").resolve()
        if (batch.get("name") != name or Path(batch.get("path", "")).resolve() != sql_path
                or not sql_path.is_file() or sql_path.stat().st_size != batch.get("sql_bytes")
                or sha256(sql_path) != batch.get("sha256")):
            raise ValueError(f"SQL batch changed: {name}")
        sql = sql_path.read_text(encoding="utf-8").lower()
        if ("pdr3_wide.random as r" not in sql
                or "pdr3_wide.patch_qa as q" not in sql
                or "conesearch(r.coord," not in sql
                or "limit 200" not in sql
                or "i_pix_variance" not in sql or "imag_psf_depth" not in sql
                or any(token in sql for token in ("photoz", "shape catalog", "forced photometry"))):
            raise ValueError(f"SQL table or science-data scope changed: {name}")
        cursor += int(batch["n_candidates"])
    if cursor != len(candidates):
        raise ValueError("batch candidate count differs from the approved ledger")
    return plan, candidates, batches, plan_hash


def validate_csv(path: Path, expected_ids: set[str]) -> tuple[int, dict[str, int]]:
    counts = {candidate_id: 0 for candidate_id in expected_ids}
    seen: set[tuple[str, str]] = set()
    nullable = {"i_pix_variance", "i_sky_mean", "i_sky_std", "i_inputcount_value",
                "imag_psf_depth", "iseeing"}
    nulls = {name: 0 for name in nullable}
    with Path(path).open(encoding="utf-8-sig", newline="") as stream:
        reader = catalogue_reader(stream)
        if reader.fieldnames != OUTPUT_COLUMNS:
            raise ValueError("CAS returned columns outside the approved 29-column schema")
        for row in reader:
            if None in row or any(value is None for value in row.values()):
                raise ValueError("malformed CAS candidate row")
            candidate_id = row["candidate_id"]
            if candidate_id not in counts:
                raise ValueError("CAS returned an unrequested S23B candidate")
            object_id = str(row["random_object_id"]).strip()
            key = candidate_id, object_id
            if not object_id or key in seen:
                raise ValueError("CAS returned an empty or duplicate random-position ID")
            seen.add(key)
            if counts[candidate_id] >= 200:
                raise ValueError("CAS returned more than 200 positions for one candidate")
            ra, dec = parse_float(row["ra"], "ra"), parse_float(row["dec"], "dec")
            if ra is None or not 0 <= ra < 360 or dec is None or not -90 <= dec <= 90:
                raise ValueError("CAS returned invalid random-position coordinates")
            for field in ("i_pix_variance", "i_sky_mean", "i_sky_std", "i_inputcount_value",
                          "imag_psf_depth", "iseeing"):
                if parse_float(row[field], field, nullable=True) is None:
                    nulls[field] += 1
            if row["i_pixelflags_bright_objectcenter"].lower() not in ("false", "f", "0"):
                raise ValueError("CAS returned a position flagged at a bright-object center")
            if row["i_pixelflags_edge"].lower() not in ("false", "f", "0"):
                raise ValueError("CAS returned an edge position")
            counts[candidate_id] += 1
    return sum(counts.values()), nulls


def cas_status(client: Cas, job_id: str) -> dict:
    for attempt in range(10):
        try:
            return client.json("status", id=job_id)
        except urllib.error.HTTPError as error:
            if error.code not in (502, 503, 504) or attempt == 9:
                raise
            print(f"CAS gateway HTTP {error.code}; retrying in 60 seconds ({attempt + 1}/10)", flush=True)
            time.sleep(60)


def run(plan_dir: Path, results_dir: Path, user: str, *, submit: bool,
        approved_plan_sha256: str | None) -> dict | None:
    plan_dir, results_dir = Path(plan_dir).resolve(), Path(results_dir).resolve()
    plan, candidates, batches, plan_hash = load_plan(plan_dir)
    print(json.dumps({
        "plan_status": plan["status"], "plan_sha256": plan_hash,
        "candidates": len(candidates), "batches": len(batches),
        "maximum_rows": plan["maximum_returned_rows"],
        "estimated_csv_bytes": plan["estimated_csv_bytes"],
        "output_columns": OUTPUT_COLUMNS, "images": False,
        "source_rows": False, "shapes": False, "photoz": False,
        "network_used": False,
    }, ensure_ascii=False), flush=True)
    if not submit:
        if approved_plan_sha256:
            raise ValueError("approval hash is accepted only with --submit")
        return None
    if str(approved_plan_sha256 or "").lower() != plan_hash:
        raise ValueError("submission requires explicit approval of this exact plan SHA-256")
    if not sys.stdin.isatty():
        raise ValueError("submission requires an interactive hidden credential prompt")
    if not platform.node().lower().startswith("lageunha"):
        raise ValueError("run HSC CAS metadata queries from lageunha")

    if results_dir.exists():
        checkpoint = results_dir / "status.json"
        if not checkpoint.is_file():
            raise FileExistsError("refusing an existing results directory without a checkpoint")
        state = json.loads(checkpoint.read_text(encoding="utf-8"))
        if state.get("plan_sha256") != plan_hash:
            raise ValueError("existing checkpoint belongs to a different approved plan")
    else:
        results_dir.mkdir(parents=True, exist_ok=False)
        checkpoint = results_dir / "status.json"
        state = {
            "schema": "astrafek-s23b-icl-control-candidate-run-v1",
            "status": "pending", "plan_sha256": plan_hash, "jobs": {},
            "network_used": False, "query_submitted": False,
            "credentials_saved": False, "images_returned": False,
            "source_detections_returned": False, "shape_rows_returned": False,
            "photoz_rows_returned": False, "science_measurement": False,
        }
        save(checkpoint, state)
    protected = ("credentials_saved", "images_returned", "source_detections_returned",
                 "shape_rows_returned", "photoz_rows_returned", "science_measurement")
    if any(state.get(key) is not False for key in protected):
        raise ValueError("saved checkpoint exceeds the approved export scope")
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
            candidate_ids = [row["candidate_id"] for row in candidates]
            for index, batch in enumerate(batches):
                name = batch["name"]
                item = state["jobs"].setdefault(name, {"status": "pending"})
                target = results_dir / f"{name}.csv"
                ids = set(candidate_ids[index * 10:index * 10 + int(batch["n_candidates"])])
                if item.get("status") == "complete":
                    if not target.is_file() or sha256(target) != item.get("sha256"):
                        raise ValueError(f"saved result missing or changed: {name}")
                    validate_csv(target, ids)
                    continue
                if item.get("status") == "error":
                    raise ValueError(f"CAS job failed for {name}; review checkpoint before retry")
                if "job_id" not in item:
                    item["status"] = "submitting"
                    save(checkpoint, state)
                    try:
                        result = client.json("submit", catalog_job=dict(
                            sql=(plan_dir / f"{name}.sql").read_text(encoding="utf-8"),
                            out_format="csv", include_metainfo_to_body=False,
                            release_version="pdr3"), nomail=True, skip_syntax_check=False)
                    except urllib.error.HTTPError:
                        item["status"] = "pending"
                        save(checkpoint, state)
                        raise
                    item.update(status="submitted", job_id=str(result["id"]))
                    state["query_submitted"] = True
                    save(checkpoint, state)
                while True:
                    response = cas_status(client, item["job_id"])
                    if response.get("status") == "done":
                        break
                    if response.get("status") == "error":
                        item.update(status="error", query_error=str(response.get("error", "CAS query failed"))[:500])
                        save(checkpoint, state)
                        raise ValueError(f"CAS random-candidate query failed for {name}")
                    print(f"{name}: {response.get('status', 'unknown')}", flush=True)
                    time.sleep(POLL_SECONDS)
                partial = target.with_suffix(".csv.part")
                with client.post("download", id=item["job_id"]) as response, partial.open("xb") as stream:
                    total_bytes = 0
                    for block in iter(lambda: response.read(1 << 20), b""):
                        total_bytes += len(block)
                        if total_bytes > MAX_BATCH_BYTES:
                            raise ValueError(f"CAS response exceeded size limit for {name}")
                        stream.write(block)
                rows, nulls = validate_csv(partial, ids)
                partial.replace(target)
                item.update(status="complete", rows=rows, bytes=target.stat().st_size,
                            sha256=sha256(target), nullable_field_null_counts=nulls)
                save(checkpoint, state)
                print(f"{name}: verified {rows} random-position rows", flush=True)

            total = 0
            null_totals = {key: 0 for key in ("i_pix_variance", "i_sky_mean", "i_sky_std",
                                               "i_inputcount_value", "imag_psf_depth", "iseeing")}
            for index, batch in enumerate(batches):
                name = batch["name"]
                target = results_dir / f"{name}.csv"
                ids = set(candidate_ids[index * 10:index * 10 + int(batch["n_candidates"])])
                rows, nulls = validate_csv(target, ids)
                if state["jobs"][name].get("rows") != rows:
                    raise ValueError(f"checkpoint row count differs from CSV: {name}")
                total += rows
                for key, count in nulls.items():
                    null_totals[key] += count
            if total > 208600:
                raise ValueError("returned rows exceed the approved maximum")
            summary = {
                "schema": "astrafek-s23b-icl-control-candidate-summary-v1",
                "status": "complete", "plan_sha256": plan_hash,
                "n_target_candidates": 1043, "n_random_position_rows": total,
                "max_rows_per_target": 200, "nullable_field_null_counts": null_totals,
                "images_returned": False, "source_detections_returned": False,
                "shape_rows_returned": False, "photoz_rows_returned": False,
                "credentials_saved": False, "science_measurement": False,
                "interpretation": "Candidate random positions and scalar QA only; matched controls require downstream quality, contamination, depth, seeing, and spatial-screening audits.",
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
    return summary


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan-dir", type=Path, default=DEFAULT_PLAN)
    parser.add_argument("--results", type=Path, default=DEFAULT_RESULTS)
    parser.add_argument("--user", default=DEFAULT_USER)
    parser.add_argument("--check", action="store_true", help="validate the pinned plan without network")
    parser.add_argument("--submit", action="store_true", help="submit after approval of the exact pinned plan")
    parser.add_argument("--approved-plan-sha256", help="required with --submit")
    args = parser.parse_args()
    if args.check == args.submit:
        parser.error("select exactly one of --check or --submit")
    if args.check and args.approved_plan_sha256:
        parser.error("approval hash cannot be used with --check")
    run(args.plan_dir, args.results, args.user, submit=args.submit,
        approved_plan_sha256=args.approved_plan_sha256)
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except urllib.error.HTTPError as error:
        print(f"HSC CAS HTTP {error.code}; credentials and response body suppressed.", file=sys.stderr)
        raise SystemExit(1)
    except KeyboardInterrupt:
        print("\nStopped; submitted query IDs remain checkpointed.", file=sys.stderr)
        raise SystemExit(130)
