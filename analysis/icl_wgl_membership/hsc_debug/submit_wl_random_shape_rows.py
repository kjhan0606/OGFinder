#!/usr/bin/env python3
"""Submit an explicitly approved HSC random-center shape-row export plan."""
from __future__ import annotations

import argparse
import csv
import getpass
import json
import math
from pathlib import Path
import platform
import sys
import time
import urllib.error

from fetch_ancillary_catalogues import Cas, catalogue_reader, save
from prepare_wl_random_shape_rows import (
    OUTPUT_COLUMNS, load_verified_counts, shape_rows_query, sha,
)


DEFAULT_USER = "kjhan0606"
MAX_BATCH_BYTES = 1024 * 1024 * 1024
NUMERIC_FIELDS = {
    "random_ra_deg", "random_dec_deg", "lens_z", "source_ra_deg", "source_dec_deg",
    "e1", "e2", "sigma_e", "rms_e", "weight", "m", "c1", "c2", "mag_a10",
    "resolution", "psf_shape11", "psf_shape22", "psf_shape12",
}


def validate_plan(plan_dir):
    plan_dir = Path(plan_dir).resolve()
    manifest_path = plan_dir / "random_shape_row_plan.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("schema") != "hsc-s19a-random-shape-row-export-plan-v1":
        raise ValueError("unexpected random shape-row plan schema")
    if manifest.get("status") != "prepared_not_submitted":
        raise ValueError("random shape-row plan is not in prepared state")
    if manifest.get("output_columns") != OUTPUT_COLUMNS:
        raise ValueError("unexpected random shape-row columns")
    if any(manifest.get(key) is not False for key in (
            "source_shape_rows_returned", "photoz_pdf_rows_included", "images_returned",
            "network_used", "query_submitted", "science_measurement", "credentials_saved")):
        raise ValueError("plan status contains data outside the declared scope")
    if manifest.get("execution_requires_explicit_user_approval") is not True:
        raise ValueError("plan does not require explicit user approval")

    count_plan_path = Path(manifest["source_count_plan"]).resolve()
    count_results_dir = Path(manifest["source_count_results_dir"]).resolve()
    count_manifest, centers, count_batches, counts, count_records, count_plan_sha, summary = \
        load_verified_counts(count_plan_path.parent, count_results_dir)
    if (count_plan_sha != manifest.get("source_count_plan_sha256")
            or sha(count_results_dir / "count_summary.json") != manifest.get("source_count_summary_sha256")
            or len(centers) != manifest.get("n_random_centers")
            or count_manifest.get("n_matched_lenses") != manifest.get("n_matched_lenses")
            or int(summary["sum_shape_source_counts"]) != manifest.get("expected_random_source_pair_rows")
            or count_records != manifest.get("count_batches")):
        raise ValueError("count evidence differs from the immutable row plan")

    batch_size = int(manifest["batch_size"])
    expected_n_batches = (len(centers) + batch_size - 1) // batch_size
    batches = manifest.get("batches", [])
    if len(batches) != expected_n_batches or len(batches) != manifest.get("n_query_batches"):
        raise ValueError("random shape-row batch count mismatch")
    for index, batch in enumerate(batches):
        path = Path(batch["path"]).resolve()
        if path.parent != plan_dir or path.name != f"random_shape_rows_{index:03d}.sql":
            raise ValueError("row SQL path escapes its immutable plan directory")
        chunk = centers[index * batch_size:(index + 1) * batch_size]
        expected_sql = "\nUNION ALL\n".join(shape_rows_query(row) for row in chunk) + "\n"
        expected_rows = sum(counts[row["random_center_id"]]["n_shape_sources"] for row in chunk)
        if (path.read_text(encoding="utf-8") != expected_sql
                or sha(path) != batch.get("sha256")
                or batch.get("n_random_centers") != len(chunk)
                or batch.get("expected_pair_rows_from_count_only") != expected_rows):
            raise ValueError(f"random shape-row SQL or expected count changed: {path.name}")
    return manifest, centers, batches, counts, count_plan_sha, sha(manifest_path)


def read_shape_rows(path, expected_centers, expected_counts):
    expected = {row["random_center_id"]: row for row in expected_centers}
    counts = {center_id: 0 for center_id in expected}
    seen = {center_id: set() for center_id in expected}
    with Path(path).open(encoding="utf-8-sig", newline="") as stream:
        reader = catalogue_reader(stream)
        if reader.fieldnames != OUTPUT_COLUMNS:
            raise ValueError("CAS returned columns outside the approved shape-row schema")
        for row in reader:
            center_id = str(row.get("random_center_id", ""))
            center = expected.get(center_id)
            if center is None:
                raise ValueError("CAS returned an unrequested random center")
            if (str(row["matched_lens_id"]) != str(center["matched_lens_id"])
                    or abs(float(row["random_ra_deg"]) - float(center["ra_deg"])) > 1e-9
                    or abs(float(row["random_dec_deg"]) - float(center["dec_deg"])) > 1e-9
                    or abs(float(row["lens_z"]) - float(center["z_lens"])) > 1e-10):
                raise ValueError("CAS returned a changed random-center/lens association")
            object_id = str(row.get("object_id", ""))
            if not object_id or object_id in seen[center_id]:
                raise ValueError("CAS returned an empty or duplicated center-source pair")
            seen[center_id].add(object_id)
            for field in NUMERIC_FIELDS:
                value = float(row[field])
                if not math.isfinite(value):
                    raise ValueError(f"non-finite value in returned shape field {field}")
            source_ra, source_dec = float(row["source_ra_deg"]), float(row["source_dec_deg"])
            if not 0 <= source_ra < 360 or not -90 <= source_dec <= 90:
                raise ValueError("CAS returned invalid source coordinates")
            mask = int(row["b_mode_mask_raw"])
            zbin = int(row["source_zbin"])
            if str(mask) != str(row["b_mode_mask_raw"]).strip() or mask not in (0, 1):
                raise ValueError("CAS returned a non-binary raw B-mode mask")
            if str(zbin) != str(row["source_zbin"]).strip() or zbin < 0:
                raise ValueError("CAS returned an invalid HSC tomographic bin")
            if float(row["weight"]) <= 0 or not 0 <= float(row["rms_e"]) < 1:
                raise ValueError("CAS returned an invalid source shape weight or rms_e")
            counts[center_id] += 1
    if counts != {key: int(expected_counts[key]["n_shape_sources"]) for key in expected}:
        raise ValueError("shape-row response counts differ from approved per-center preflight")
    return counts


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
    manifest, centers, batches, counts, count_plan_sha, plan_sha = validate_plan(plan_dir)
    print(json.dumps({
        "plan_status": manifest["status"], "network_used": False,
        "random_centers": len(centers), "matched_lenses": manifest["n_matched_lenses"],
        "batches": len(batches), "expected_pair_rows": manifest["expected_random_source_pair_rows"],
        "output_columns": OUTPUT_COLUMNS, "source_shape_rows_requested": True,
        "photoz_pdf_rows": False, "images": False,
    }, ensure_ascii=False), flush=True)
    if not submit:
        return 0
    if not approved:
        raise ValueError("submission requires --approved-random-shape-row-export after explicit user approval")
    if not sys.stdin.isatty():
        raise ValueError("submission requires one hidden password prompt in an interactive terminal")
    if not platform.node().lower().startswith("lageunha"):
        raise ValueError("run approved HSC exports from the established lageunha interactive environment")
    if results_dir.exists():
        checkpoint = results_dir / "row_status.json"
        if not checkpoint.is_file():
            raise FileExistsError("refusing existing row-results directory without a checkpoint")
        state = json.loads(checkpoint.read_text(encoding="utf-8"))
        if state.get("plan_sha256") != plan_sha:
            raise ValueError("existing checkpoint belongs to another row plan")
    else:
        results_dir.mkdir(parents=True, exist_ok=False)
        checkpoint = results_dir / "row_status.json"
        state = {
            "schema": "hsc-s19a-random-shape-row-run-v1", "status": "pending",
            "plan_sha256": plan_sha, "count_plan_sha256": count_plan_sha,
            "jobs": {}, "credentials_saved": False,
            "source_shape_rows_downloaded": False, "photoz_pdf_rows_included": False,
            "images_returned": False, "science_measurement": False,
        }
        save(checkpoint, state)
    if state.get("credentials_saved") is not False or any(
            state.get(key) is not False for key in (
                "source_shape_rows_downloaded", "photoz_pdf_rows_included",
                "images_returned", "science_measurement")):
        raise ValueError("checkpoint contains data outside the approved random shape-row scope")
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
                name = f"random_shape_rows_{index:03d}"
                item = state["jobs"].setdefault(name, {"status": "pending"})
                target = results_dir / f"{name}.csv"
                chunk = centers[index * size:(index + 1) * size]
                expected_counts = {center["random_center_id"]: counts[center["random_center_id"]]
                                   for center in chunk}
                if item.get("status") == "complete":
                    if not target.is_file() or sha(target) != item.get("sha256"):
                        raise ValueError(f"saved shape-row batch missing or changed: {name}")
                    read_shape_rows(target, chunk, expected_counts)
                    continue
                if item.get("status") == "error":
                    raise ValueError(f"CAS query failed for {name}; inspect before retry")
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
                        raise ValueError(f"CAS query failed for {name}")
                    print(f"{name}: {status.get('status', 'unknown')}", flush=True)
                    time.sleep(60)
                partial = results_dir / f"{name}.csv.part"
                if partial.exists():
                    partial.rename(results_dir / f"{name}.csv.incomplete-{time.time_ns()}")
                with client.post("download", id=item["job_id"]) as response, partial.open("xb") as stream:
                    total_bytes = 0
                    for block in iter(lambda: response.read(1024 * 1024), b""):
                        total_bytes += len(block)
                        if total_bytes > MAX_BATCH_BYTES:
                            raise ValueError(f"CAS response exceeds the {MAX_BATCH_BYTES}-byte batch safety limit")
                        stream.write(block)
                read = read_shape_rows(partial, chunk, expected_counts)
                partial.replace(target)
                item.update(status="complete", rows=sum(read.values()),
                            bytes=target.stat().st_size, sha256=sha(target))
                save(checkpoint, state)
                print(f"{name}: verified {item['rows']} random source pairs", flush=True)
            total_rows = 0
            total_bytes = 0
            for index, batch in enumerate(batches):
                name = f"random_shape_rows_{index:03d}"
                item = state["jobs"].get(name, {})
                path = results_dir / f"{name}.csv"
                if item.get("status") != "complete" or not path.is_file() or sha(path) != item.get("sha256"):
                    raise ValueError(f"final random shape-row batch audit failed: {name}")
                total_rows += int(item["rows"])
                total_bytes += int(item["bytes"])
            if total_rows != int(manifest["expected_random_source_pair_rows"]):
                raise ValueError("final random shape-row total differs from count-only plan")
            summary = {
                "schema": "hsc-s19a-random-shape-row-summary-v1", "status": "complete",
                "plan_sha256": plan_sha, "count_plan_sha256": count_plan_sha,
                "n_random_centers": manifest["n_random_centers"],
                "n_matched_lenses": manifest["n_matched_lenses"],
                "n_source_pairs": total_rows, "n_csv_bytes": total_bytes,
                "source_shape_rows_downloaded": True, "photoz_pdf_rows_included": False,
                "images_returned": False, "science_measurement": False,
            }
            state.update(status="complete", summary=summary)
            save(checkpoint, state)
            (results_dir / "row_summary.json").write_text(
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
    parser.add_argument("--check", action="store_true", help="validate immutable plan only; no network")
    parser.add_argument("--submit", action="store_true", help="submit shape-row export")
    parser.add_argument("--approved-random-shape-row-export", action="store_true",
                        help="required after explicit user approval; no image/photo-z PDF export")
    args = parser.parse_args()
    if args.check == args.submit:
        parser.error("select exactly one of --check or --submit")
    if args.check and args.approved_random_shape_row_export:
        parser.error("approval flag is not used with --check")
    return run(args.plan_dir, args.results, args.user, submit=args.submit,
               approved=args.approved_random_shape_row_export)


if __name__ == "__main__":
    try:
        sys.exit(main())
    except urllib.error.HTTPError as error:
        print(f"HSC CAS HTTP {error.code}; credentials and response body suppressed.", file=sys.stderr)
        sys.exit(1)
    except KeyboardInterrupt:
        print("\nStopped; submitted query IDs remain checkpointed.", file=sys.stderr)
        sys.exit(130)
