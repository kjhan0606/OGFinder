#!/usr/bin/env python3
"""Validate or submit the pinned HSC random-point variance export plan."""
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


DEFAULT_USER = "kjhan0606"
DEFAULT_PLAN = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_random_pixel_variance_plan_2026-11-06_03")
DEFAULT_RESULTS = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_random_pixel_variance_results_2026-11-06_01")
EXPECTED_PLAN_SHA256 = "eeef98ef2210337f714e36d8a4afcd352fe312a536c90e8e236ff9f843f53e52"
EXPECTED_SQL_AGGREGATE_SHA256 = "5d82df638435ee39029dc51acdcf580ab6e242d21562bc5523ef9c807368f585"
MAX_BATCH_BYTES = 25 * 1024 * 1024
POLL_SECONDS = 20
OUTPUT_COLUMNS = [
    "center_id", "sample_role", "target_cluster_id", "center_ra_deg",
    "center_dec_deg", "cluster_z_for_scale", "random_object_id", "ra", "dec",
    "skymap_id", "tract", "patch", "patch_s", "adjust_density",
    "i_pix_variance", "i_sky_mean", "i_sky_std", "i_inputcount_value",
    "i_pixelflags_interpolatedcenter", "i_pixelflags_saturatedcenter",
    "i_pixelflags_crcenter", "i_pixelflags_suspectcenter",
    "i_pixelflags_bright_objectcenter", "i_pixelflags_edge",
    "i_sdssshape_psf_shape11", "i_sdssshape_psf_shape22",
    "i_sdssshape_psf_shape12",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def num(value: str, name: str, *, nullable: bool = False,
        nonnegative: bool = False) -> float | None:
    value = str(value).strip()
    if nullable and value.lower() in ("", "null", "none", "nan", r"\n", r"\N"):
        return None
    try:
        out = float(value)
    except (TypeError, ValueError):
        raise ValueError(f"CAS returned nonnumeric {name}") from None
    if not math.isfinite(out) or (nonnegative and out < 0):
        raise ValueError(f"CAS returned invalid {name}")
    return out


def load_plan(plan_dir: Path) -> tuple[dict, list[dict], list[dict], str]:
    plan_dir = Path(plan_dir).resolve()
    plan_path = plan_dir / "plan.json"
    plan = json.loads(plan_path.read_text(encoding="utf-8"))
    plan_hash = sha256(plan_path)
    if plan_hash != EXPECTED_PLAN_SHA256:
        raise ValueError("prepared plan differs from the pinned approval scope")
    if (plan.get("schema") != "astrafek-icl-random-pixel-variance-plan-v1"
            or plan.get("status") != "prepared_not_submitted"
            or plan.get("table") != "pdr3_wide.random"
            or plan.get("n_centers") != 1598 or plan.get("n_target_centers") != 800
            or plan.get("n_control_centers") != 798 or plan.get("n_batches") != 16
            or plan.get("annulus_proper_kpc") != [200.0, 300.0]
            or plan.get("output_columns") != OUTPUT_COLUMNS
            or plan.get("estimated_rows") != 46424):
        raise ValueError("plan center count, table, annulus, or column scope changed")
    if (plan.get("execution_requires_explicit_user_approval") is not True
            or plan.get("query_submitted") is not False or plan.get("network_used") is not False
            or plan.get("images_returned") is not False
            or plan.get("shape_rows_returned") is not False
            or plan.get("photoz_rows_returned") is not False
            or plan.get("object_catalog_rows_returned") is not False
            or plan.get("credentials_saved") is not False
            or plan.get("science_measurement") is not False):
        raise ValueError("plan scope/provenance flags changed")
    for filename, digest in plan.get("input_sha256", {}).items():
        path = Path(filename)
        if not path.is_file() or sha256(path) != digest:
            raise ValueError(f"plan input changed: {path}")
    centers = plan.get("centers", [])
    center_ids = {row["center_id"] for row in centers}
    if (len(centers) != 1598 or len(center_ids) != 1598
            or sum(row["role"] == "target" for row in centers) != 800
            or sum(row["role"] == "control" for row in centers) != 798):
        raise ValueError("plan center ledger is incomplete")
    batches = plan.get("batches", [])
    if len(batches) != 16:
        raise ValueError("expected sixteen SQL batches")
    flattened = []
    for batch in batches:
        ids = batch.get("center_ids", [])
        if len(ids) != batch.get("centers"):
            raise ValueError("batch center ledger has an unexpected size")
        flattened.extend(ids)
    if flattened != [row["center_id"] for row in centers] or len(set(flattened)) != len(centers):
        raise ValueError("batch center ledger differs from the plan center table")
    aggregate = hashlib.sha256("\n".join(row["sha256"] for row in batches).encode()).hexdigest()
    if aggregate != EXPECTED_SQL_AGGREGATE_SHA256 or aggregate != plan.get("sql_aggregate_sha256"):
        raise ValueError("SQL aggregate checksum changed")
    for index, batch in enumerate(batches):
        name = f"icl_random_variance_{index:03d}"
        path = (plan_dir / f"{name}.sql").resolve()
        if (batch.get("name") != name or Path(batch.get("path", "")).resolve() != path
                or not path.is_file() or path.stat().st_size != batch.get("sql_bytes")
                or sha256(path) != batch.get("sha256")):
            raise ValueError(f"SQL batch differs from reviewed plan: {name}")
        sql = path.read_text(encoding="utf-8").lower()
        if ("pdr3_wide.random" not in sql or "i_pix_variance" not in sql
                or "adjust_density < 0.20" not in sql or "select *" in sql
                or any(token in sql for token in ("photoz", "shape catalog", "forced"))):
            raise ValueError(f"SQL scope validation failed: {name}")
    return plan, centers, batches, plan_hash


def validate_rows(path: Path, expected_centers: list[dict]) -> tuple[int, set[str], dict[str, int]]:
    expected = {row["center_id"]: row for row in expected_centers}
    seen = set()
    center_ids = set()
    nulls = {"i_pix_variance": 0, "i_sky_std": 0, "i_inputcount_value": 0}
    count = 0
    with Path(path).open(encoding="utf-8-sig", newline="") as stream:
        reader = catalogue_reader(stream)
        if reader.fieldnames != OUTPUT_COLUMNS:
            raise ValueError("CAS returned columns outside the approved schema")
        for row in reader:
            if None in row or any(value is None for value in row.values()):
                raise ValueError("malformed CAS random-point row")
            center = expected.get(row["center_id"])
            if center is None:
                raise ValueError("CAS returned an unrequested center")
            if (row["sample_role"] != center["role"]
                    or row["target_cluster_id"] != center["target_cluster_id"]
                    or abs(float(row["center_ra_deg"]) - center["ra_deg"]) > 1e-8
                    or abs(float(row["center_dec_deg"]) - center["dec_deg"]) > 1e-8
                    or abs(float(row["cluster_z_for_scale"]) - center["z"]) > 1e-8):
                raise ValueError("CAS center identity differs from approved plan")
            object_key = (row["center_id"], row["random_object_id"])
            if not row["random_object_id"] or object_key in seen:
                raise ValueError("duplicate or empty random-point ID at a center")
            seen.add(object_key)
            center_ids.add(row["center_id"])
            adjust = num(row["adjust_density"], "adjust_density", nonnegative=True)
            if adjust is None or adjust >= 0.2:
                raise ValueError("CAS returned point outside approved density thinning")
            ra, dec = num(row["ra"], "ra"), num(row["dec"], "dec")
            if ra is None or dec is None or not -90 <= dec <= 90:
                raise ValueError("CAS returned invalid random-point coordinates")
            dra = math.radians(ra - center["ra_deg"])
            dec0, decr = math.radians(center["dec_deg"]), math.radians(dec)
            cosine = (math.sin(dec0) * math.sin(decr)
                      + math.cos(dec0) * math.cos(decr) * math.cos(dra))
            radius = math.degrees(math.acos(max(-1.0, min(1.0, cosine)))) * 3600.0
            if not center["inner_arcsec"] - .02 <= radius <= center["outer_arcsec"] + .02:
                raise ValueError("CAS returned point outside the approved 200--300 kpc annulus")
            for field in nulls:
                value = num(row[field], field, nullable=True, nonnegative=True)
                if value is None:
                    nulls[field] += 1
            if row["i_pixelflags_bright_objectcenter"].lower() not in ("false", "f", "0"):
                raise ValueError("CAS returned a bright-object-center point")
            if row["i_pixelflags_edge"].lower() not in ("false", "f", "0"):
                raise ValueError("CAS returned an edge point")
            count += 1
    return count, center_ids, nulls


def cas_status(client: Cas, job_id: str) -> dict:
    for attempt in range(10):
        try:
            return client.json("status", id=job_id)
        except urllib.error.HTTPError as error:
            if error.code not in (502, 503, 504) or attempt == 9:
                raise
            print(f"CAS status gateway HTTP {error.code}; retrying in 60 seconds ({attempt + 1}/10)", flush=True)
            time.sleep(60)


def run(plan_dir: Path, results_dir: Path, user: str, *, submit: bool,
        approved_plan_sha256: str | None) -> int:
    plan_dir, results_dir = Path(plan_dir).resolve(), Path(results_dir).resolve()
    plan, centers, batches, plan_hash = load_plan(plan_dir)
    print(json.dumps({"plan_sha256": plan_hash, "centers": len(centers),
        "batches": len(batches), "estimated_rows": plan["estimated_rows"],
        "estimated_csv_size": plan["estimated_csv_size"], "columns": OUTPUT_COLUMNS,
        "images": False, "shapes": False, "photoz": False,
        "network_used": False}, ensure_ascii=False), flush=True)
    if not submit:
        if approved_plan_sha256:
            raise ValueError("approval hash accepted only with --submit")
        return 0
    if str(approved_plan_sha256 or "").lower() != plan_hash:
        raise ValueError("submission requires explicit approval of this exact plan SHA-256")
    if not sys.stdin.isatty():
        raise ValueError("submission requires an interactive hidden credential prompt")
    if not platform.node().lower().startswith("lageunha"):
        raise ValueError("run approved HSC PDR3 query from lageunha")
    if results_dir.exists():
        checkpoint = results_dir / "status.json"
        if not checkpoint.is_file():
            raise FileExistsError("refusing existing result directory without checkpoint")
        state = json.loads(checkpoint.read_text(encoding="utf-8"))
        if state.get("plan_sha256") != plan_hash:
            raise ValueError("existing checkpoint belongs to a different plan")
    else:
        results_dir.mkdir(parents=True, exist_ok=False)
        checkpoint = results_dir / "status.json"
        state = {"schema": "astrafek-icl-random-variance-run-v1", "status": "pending",
            "plan_sha256": plan_hash, "jobs": {}, "network_used": False,
            "query_submitted": False, "credentials_saved": False,
            "images_returned": False, "shape_rows_returned": False,
            "photoz_rows_returned": False, "object_catalog_rows_returned": False,
            "science_measurement": False}
        save(checkpoint, state)
    if any(state.get(k) is not False for k in (
            "credentials_saved", "images_returned", "shape_rows_returned",
            "photoz_rows_returned", "object_catalog_rows_returned", "science_measurement")):
        raise ValueError("checkpoint exceeds approved scope")
    if any(v.get("status") == "submitting" and "job_id" not in v for v in state.get("jobs", {}).values()):
        raise ValueError("submission outcome ambiguous; inspect CAS before retry")

    with (results_dir / "run.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        password = getpass.getpass("HSC password (hidden; one prompt): ")
        if not password:
            raise ValueError("empty password")
        client = Cas(user, password)
        state.update(status="running", network_used=True)
        save(checkpoint, state)
        try:
            center_lookup = {row["center_id"]: row for row in centers}
            for batch in batches:
                name = batch["name"]
                item = state["jobs"].setdefault(name, {"status": "pending"})
                target = results_dir / f"{name}.csv"
                sql_path = plan_dir / f"{name}.sql"
                expected = [center_lookup[key] for key in batch["center_ids"]]
                if item.get("status") == "complete":
                    if not target.is_file() or sha256(target) != item.get("sha256"):
                        raise ValueError(f"saved batch missing or changed: {name}")
                    row_count, _, _ = validate_rows(target, expected)
                    if row_count != item.get("rows"):
                        raise ValueError(f"saved row count changed: {name}")
                    continue
                if item.get("status") == "error":
                    raise ValueError(f"CAS job failed for {name}; inspect before retry")
                if "job_id" not in item:
                    item["status"] = "submitting"
                    save(checkpoint, state)
                    try:
                        reply = client.json("submit", catalog_job=dict(
                            sql=sql_path.read_text(encoding="utf-8"), out_format="csv",
                            include_metainfo_to_body=False, release_version="pdr3"),
                            nomail=True, skip_syntax_check=False)
                    except urllib.error.HTTPError:
                        item["status"] = "pending"
                        save(checkpoint, state)
                        raise
                    item.update(status="submitted", job_id=str(reply["id"]))
                    state["query_submitted"] = True
                    save(checkpoint, state)
                while True:
                    status = cas_status(client, item["job_id"])
                    if status.get("status") == "done":
                        break
                    if status.get("status") == "error":
                        item.update(status="error", query_error=str(status.get("error", "CAS query failed"))[:500])
                        save(checkpoint, state)
                        raise ValueError(f"CAS query failed for {name}")
                    print(f"{name}: {status.get('status', 'unknown')}", flush=True)
                    time.sleep(POLL_SECONDS)
                partial = results_dir / f"{name}.csv.part"
                if partial.exists():
                    partial.rename(results_dir / f"{name}.csv.incomplete-{time.time_ns()}")
                with client.post("download", id=item["job_id"]) as response, partial.open("xb") as stream:
                    total_bytes = 0
                    for block in iter(lambda: response.read(1 << 20), b""):
                        total_bytes += len(block)
                        if total_bytes > MAX_BATCH_BYTES:
                            raise ValueError(f"CAS response exceeded size cap for {name}")
                        stream.write(block)
                expected = [center_lookup[key] for key in batch["center_ids"]]
                count, returned_centers, nulls = validate_rows(partial, expected)
                partial.replace(target)
                item.update(status="complete", rows=count, centers=len(returned_centers),
                            bytes=target.stat().st_size, sha256=sha256(target), null_counts=nulls)
                save(checkpoint, state)
                print(f"{name}: verified {count} random points from {len(returned_centers)} centers", flush=True)

            all_rows = 0
            observed_centers = set()
            hashes = {}
            null_totals = {"i_pix_variance": 0, "i_sky_std": 0, "i_inputcount_value": 0}
            for batch in batches:
                name = batch["name"]
                item = state["jobs"][name]
                target = results_dir / f"{name}.csv"
                if item.get("status") != "complete" or not target.is_file() or sha256(target) != item.get("sha256"):
                    raise ValueError(f"incomplete output batch: {name}")
                expected = [center_lookup[key] for key in batch["center_ids"]]
                count, got, nulls = validate_rows(target, expected)
                if count != item.get("rows"):
                    raise ValueError(f"final row count mismatch: {name}")
                all_rows += count
                observed_centers.update(got)
                for field, value in nulls.items():
                    null_totals[field] += value
                hashes[name] = item["sha256"]
            summary = {"schema": "astrafek-icl-random-variance-summary-v1",
                "status": "complete_exploratory", "plan_sha256": plan_hash,
                "n_requested_centers": len(centers), "n_centers_with_points": len(observed_centers),
                "n_returned_random_points": all_rows, "estimated_rows": plan["estimated_rows"],
                "null_counts": null_totals,
                "batch_hashes_sha256": hashlib.sha256("\n".join(f"{k}:{v}" for k,v in sorted(hashes.items())).encode()).hexdigest(),
                "images_returned": False, "shape_rows_returned": False,
                "photoz_rows_returned": False, "object_catalog_rows_returned": False,
                "credentials_saved": False, "science_measurement": False,
                "limitation": "Random-point variance and patch metadata do not replace coadd weights or provide full correlated-noise/source-Poisson correction."}
            state.update(status="complete", summary=summary)
            save(checkpoint, state)
            (results_dir / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
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
    parser.add_argument("--plan-dir", type=Path, default=DEFAULT_PLAN)
    parser.add_argument("--results", type=Path, default=DEFAULT_RESULTS)
    parser.add_argument("--user", default=DEFAULT_USER)
    action = parser.add_mutually_exclusive_group(required=True)
    action.add_argument("--check", action="store_true")
    action.add_argument("--submit", action="store_true")
    parser.add_argument("--approved-plan-sha256")
    args = parser.parse_args()
    if args.check and args.approved_plan_sha256:
        parser.error("approval hash is only accepted with --submit")
    return run(args.plan_dir, args.results, args.user, submit=args.submit,
               approved_plan_sha256=args.approved_plan_sha256)


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except urllib.error.HTTPError as error:
        print(f"HSC CAS HTTP {error.code}; credentials and response body suppressed.", file=sys.stderr)
        raise SystemExit(1)
    except KeyboardInterrupt:
        print("\nStopped; submitted query IDs remain checkpointed.", file=sys.stderr)
        raise SystemExit(130)
