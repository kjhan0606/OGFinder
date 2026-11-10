#!/usr/bin/env python3
"""Validate and, after exact approval, export the S23B WL source-row plan.

The runner is tied to one immutable post-count plan. It validates both the
80-lens and 400-control SQL against their count-only evidence, then submits
the 20 approved batches sequentially from an interactive lageunha session.
It never requests images or full photo-z PDFs and never promotes a scientific
measurement.
"""
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

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from fetch_ancillary_catalogues import Cas, catalogue_reader, save, sha
from finalize_s23b_wl_common_support_plan import load_lens_counts
from prepare_s23b_wl_common_support_plan import load_verified_counts as load_lens_count_evidence
from prepare_s23b_wl_common_support_plan import COUNT_PLAN as LENS_COUNT_PLAN
from prepare_wl_random_shape_rows import (
    OUTPUT_COLUMNS as RANDOM_COLUMNS,
    load_verified_counts as load_random_count_evidence,
    shape_rows_query as random_shape_rows_query,
)
from prepare_wl_shape_rows import OUTPUT_COLUMNS as LENS_COLUMNS, lens_query
from submit_wl_random_shape_rows import NUMERIC_FIELDS, cas_status, read_shape_rows


DEFAULT_PLAN_DIR = Path(
    "/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/"
    "wl_s23b_wl_common_support_plan_2026-11-10_01"
)
DEFAULT_RESULTS = DEFAULT_PLAN_DIR / "source_rows_post_count_01"
APPROVED_PLAN_SHA256 = "7bfb998149f30a50b5ee7a537bb703ed302b7ea269a3296c35204a836e75aabc"
DEFAULT_USER = "kjhan0606"
MAX_BATCH_BYTES = 256 * 1024 * 1024
LENS_NUMERIC_FIELDS = NUMERIC_FIELDS - {
    "random_ra_deg", "random_dec_deg", "lens_z",
}
COUNT_SUMMARY = "results/count_summary.json"


def read_tsv(path: Path) -> list[dict[str, str]]:
    with Path(path).open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def finite_fields(row: dict[str, str], fields: set[str]) -> None:
    for field in fields:
        if not math.isfinite(float(row[field])):
            raise ValueError(f"non-finite value in returned source field {field}")
    ra, dec = float(row["source_ra_deg"]), float(row["source_dec_deg"])
    if not 0 <= ra < 360 or not -90 <= dec <= 90:
        raise ValueError("CAS returned invalid source coordinates")
    mask_text = str(row["b_mode_mask_raw"]).strip()
    mask = int(mask_text)
    zbin_text = str(row["source_zbin"]).strip()
    zbin = int(zbin_text)
    if str(mask) != mask_text or mask not in (0, 1):
        raise ValueError("CAS returned a non-binary raw B-mode mask")
    if str(zbin) != zbin_text or zbin < 0:
        raise ValueError("CAS returned an invalid HSC tomographic bin")
    if float(row["weight"]) <= 0 or not 0 <= float(row["rms_e"]) < 1:
        raise ValueError("CAS returned an invalid source shape weight or rms_e")


def validate_lens_rows(path: Path, expected_centers: list[dict],
                       expected_counts: dict[str, int]) -> int:
    expected = {str(row["lens_id"]): row for row in expected_centers}
    observed = {lens_id: 0 for lens_id in expected}
    seen: set[tuple[str, str]] = set()
    with Path(path).open(encoding="utf-8-sig", newline="") as stream:
        reader = catalogue_reader(stream)
        if reader.fieldnames != LENS_COLUMNS:
            raise ValueError("CAS returned columns outside the approved 22-column lens schema")
        for row in reader:
            if None in row or any(value is None for value in row.values()):
                raise ValueError("malformed lens source-row CSV record")
            lens_id = str(row.get("lens_id", ""))
            center = expected.get(lens_id)
            object_id = str(row.get("object_id", ""))
            if center is None or not object_id:
                raise ValueError("CAS returned an unrequested or empty lens/source identifier")
            if (abs(float(row["lens_ra_deg"]) - float(center["ra_deg"])) > 1e-9
                    or abs(float(row["lens_dec_deg"]) - float(center["dec_deg"])) > 1e-9
                    or abs(float(row["lens_z"]) - float(center["z_lens"])) > 1e-10):
                raise ValueError("CAS returned changed lens center coordinates or redshift")
            pair = (lens_id, object_id)
            if pair in seen:
                raise ValueError("CAS returned a duplicate (lens_id, object_id) pair")
            seen.add(pair)
            finite_fields(row, LENS_NUMERIC_FIELDS)
            observed[lens_id] += 1
    expected_for_batch = {key: int(expected_counts[key]) for key in expected}
    if observed != expected_for_batch:
        raise ValueError("lens source-row counts differ from the approved count-only results")
    return len(seen)


def validate_plan(plan_dir: Path) -> dict:
    plan_dir = Path(plan_dir).resolve()
    manifest_path = plan_dir / "row_export_plan_post_count.json"
    manifest_bytes = manifest_path.read_bytes()
    plan_sha = __import__("hashlib").sha256(manifest_bytes).hexdigest()
    if plan_sha != APPROVED_PLAN_SHA256:
        raise ValueError("post-count plan SHA-256 differs from the approved immutable plan")
    manifest = json.loads(manifest_bytes)
    if manifest.get("schema") != "astrafek-s23b-wl-common-support-shape-row-plan-v2":
        raise ValueError("unexpected post-count source-row plan schema")
    if manifest.get("status") != "prepared_not_submitted":
        raise ValueError("source-row plan is not in prepared_not_submitted state")
    if manifest.get("execution_requires_explicit_user_approval") is not True:
        raise ValueError("source-row plan is missing its explicit-approval guard")
    for key in ("source_rows_downloaded", "query_submitted", "network_used",
                "credentials_saved", "images_included", "photoz_pdf_rows_included",
                "science_measurement"):
        if manifest.get(key) is not False:
            raise ValueError(f"source-row plan has unexpected scope/status flag: {key}")
    if (manifest.get("query_plan_sha256")
            != __import__("prepare_s23b_wl_common_support_plan").canonical_sha(manifest)):
        raise ValueError("canonical query-plan hash is invalid")
    if manifest.get("lens_output_columns") != LENS_COLUMNS:
        raise ValueError("lens output columns differ from the approved 22-column schema")
    if manifest.get("random_output_columns") != RANDOM_COLUMNS:
        raise ValueError("control output columns differ from the approved 23-column schema")

    lens_plan, lens_summary, lens_counts = load_lens_count_evidence()
    if manifest.get("lens_count_plan_sha256") != "b713eebded489c280fd9ce1b5840bc8d7d3e42d029035f8dbd414b72fb31e984":
        raise ValueError("lens count plan hash differs from its approved source")
    selected_rows = read_tsv(plan_dir / "selected_lenses.tsv")
    selected_ids = [str(row["lens_id"]) for row in selected_rows]
    if len(selected_ids) != 80 or len(set(selected_ids)) != 80:
        raise ValueError("selected lens table no longer contains 80 unique lenses")
    if any(int(row["n_raw_shape_sources"]) != lens_counts.get(str(row["lens_id"]))
           or int(row["n_raw_shape_sources"]) < 1000 for row in selected_rows):
        raise ValueError("selected-lens count table differs from verified count-only evidence")
    if sum(lens_counts[lens_id] for lens_id in selected_ids) != manifest.get(
            "expected_lens_source_pair_rows_from_approved_counts"):
        raise ValueError("lens count-only total differs from the immutable row plan")
    count_centers = {str(row["lens_id"]): row for row in lens_plan["centers"]}
    if any(lens_id not in count_centers for lens_id in selected_ids):
        raise ValueError("selected lens is missing from the approved lens count plan")
    lens_centers = [count_centers[lens_id] for lens_id in selected_ids]
    for selected, center in zip(selected_rows, lens_centers):
        if any(abs(float(selected[key]) - float(center[source])) > tolerance for key, source, tolerance in (
                ("ra_deg", "ra_deg", 1e-10), ("dec_deg", "dec_deg", 1e-10),
                ("z_lens", "z_lens", 1e-12))):
            raise ValueError("selected lens center differs from its approved count query")

    lens_batches = manifest.get("lens_batches", [])
    if (len(lens_batches) != 4 or manifest.get("n_lens_row_sql_batches") != 4
            or manifest.get("n_lenses") != 80):
        raise ValueError("lens batch design changed from 80 centers in four SQL batches")
    lens_batch_size = 25
    for index, batch in enumerate(lens_batches):
        name = f"lens_shape_rows_{index:03d}.sql"
        path = plan_dir / name
        chunk = lens_centers[index * lens_batch_size:(index + 1) * lens_batch_size]
        expected_sql = "\nUNION ALL\n".join(lens_query(center) for center in chunk) + "\n"
        expected_rows = sum(lens_counts[str(row["lens_id"])] for row in chunk)
        if (batch.get("name") != name or not path.is_file()
                or path.read_text(encoding="utf-8") != expected_sql
                or sha(path) != batch.get("sha256")
                or path.stat().st_size != batch.get("sql_bytes")
                or batch.get("n_lenses") != len(chunk)
                or batch.get("expected_rows_from_approved_counts") != expected_rows):
            raise ValueError(f"lens row SQL or count evidence changed: {name}")
        upper = expected_sql.upper()
        if (not upper.startswith("SELECT DISTINCT") or "SELECT *" in upper
                or "PHOTOZ" in upper or "PDR3_DNNZ_WIDE" in upper):
            raise ValueError(f"lens query exceeds the approved minimal source-row scope: {name}")

    count_manifest, random_centers, count_batches, random_counts, count_records, random_plan_sha, random_summary = (
        load_random_count_evidence(plan_dir, plan_dir / "results")
    )
    if random_plan_sha != manifest.get("count_only_preflight_plan_sha256"):
        raise ValueError("control source-row plan differs from the approved count-only plan")
    summary_path = plan_dir / COUNT_SUMMARY
    if sha(summary_path) != manifest.get("count_only_preflight_summary_sha256"):
        raise ValueError("control count summary checksum differs from the approved plan")
    if (len(random_centers) != 400 or len(count_batches) != 16
            or manifest.get("n_random_centers") != 400
            or manifest.get("n_random_row_sql_batches") != 16
            or random_summary.get("status") != "complete"):
        raise ValueError("control count evidence does not cover exactly 400 centers in 16 batches")
    if (sum(int(row["n_shape_sources"]) for row in random_counts.values())
            != manifest.get("expected_random_source_pair_rows")):
        raise ValueError("control count-only total differs from the immutable row plan")
    count_hashes = manifest.get("count_only_preflight_results_sha256", {})
    if {row["name"]: row["sha256"] for row in count_records} != count_hashes:
        raise ValueError("control count batch receipts differ from the approved row plan")
    per_lens = {lens_id: 0 for lens_id in selected_ids}
    for center in random_centers:
        lens_id = str(center["matched_lens_id"])
        if lens_id not in per_lens:
            raise ValueError("control count plan references an unselected lens")
        per_lens[lens_id] += 1
    if any(value != 5 for value in per_lens.values()):
        raise ValueError("control set no longer has exactly five centers per selected lens")

    random_batches = manifest.get("random_batches", [])
    if len(random_batches) != 16:
        raise ValueError("control row SQL batch count changed")
    random_batch_size = int(count_manifest["batch_size"])
    for index, batch in enumerate(random_batches):
        name = f"random_shape_rows_{index:03d}.sql"
        path = plan_dir / name
        chunk = random_centers[index * random_batch_size:(index + 1) * random_batch_size]
        expected_sql = "\nUNION ALL\n".join(random_shape_rows_query(row) for row in chunk) + "\n"
        expected_rows = sum(random_counts[row["random_center_id"]]["n_shape_sources"] for row in chunk)
        if (batch.get("name") != name or not path.is_file()
                or path.read_text(encoding="utf-8") != expected_sql
                or sha(path) != batch.get("sha256")
                or path.stat().st_size != batch.get("sql_bytes")
                or batch.get("n_random_centers") != len(chunk)
                or batch.get("expected_rows_from_count_preflight") != expected_rows
                or batch.get("count_preflight_batch_name") != f"random_shape_counts_{index:03d}"
                or batch.get("count_preflight_batch_sha256") != count_records[index]["sha256"]):
            raise ValueError(f"control row SQL or count evidence changed: {name}")
        upper = expected_sql.upper()
        if (not upper.startswith("SELECT DISTINCT") or "SELECT *" in upper
                or "PHOTOZ" in upper or "PDR3_DNNZ_WIDE" in upper):
            raise ValueError(f"control query exceeds the approved minimal source-row scope: {name}")

    combined = int(manifest["expected_lens_source_pair_rows_from_approved_counts"]
                   + manifest["expected_random_source_pair_rows"])
    if combined != manifest.get("expected_combined_source_pair_rows_from_count_only"):
        raise ValueError("combined pair total differs from the immutable count-only evidence")
    return {
        "manifest": manifest, "plan_sha256": plan_sha, "lens_centers": lens_centers,
        "lens_counts": lens_counts, "random_centers": random_centers,
        "random_counts": random_counts, "lens_batches": lens_batches,
        "random_batches": random_batches, "combined_rows": combined,
        "lens_summary": lens_summary, "random_summary": random_summary,
    }


def response_rows(client: Cas, name: str, item: dict, sql_path: Path,
                  results_dir: Path, validate) -> tuple[int, int, str]:
    target = results_dir / f"{name}.csv"
    partial = results_dir / f"{name}.csv.part"
    if partial.exists():
        partial.rename(results_dir / f"{name}.csv.incomplete-{time.time_ns()}")
    with client.post("download", id=item["job_id"]) as response, partial.open("xb") as stream:
        total_bytes = 0
        for block in iter(lambda: response.read(1024 * 1024), b""):
            if block:
                total_bytes += len(block)
                if total_bytes > MAX_BATCH_BYTES:
                    raise ValueError(f"CAS batch exceeded the {MAX_BATCH_BYTES}-byte transfer safety limit")
                stream.write(block)
    rows = validate(partial)
    partial.replace(target)
    digest = sha(target)
    return rows, target.stat().st_size, digest


def run(args) -> int:
    plan_dir, results_dir = args.plan_dir.resolve(), args.results.resolve()
    plan = validate_plan(plan_dir)
    manifest = plan["manifest"]
    print(json.dumps({
        "plan_status": manifest["status"], "plan_sha256": plan["plan_sha256"],
        "lenses": 80, "random_centers": 400, "batches": 20,
        "expected_lens_pairs": manifest["expected_lens_source_pair_rows_from_approved_counts"],
        "expected_control_pairs": manifest["expected_random_source_pair_rows"],
        "expected_combined_pairs": plan["combined_rows"],
        "columns": {"lens": 22, "control": 23},
        "estimated_uncompressed_csv_bytes": manifest["estimated_uncompressed_csv_bytes"],
        "images": False, "full_photoz_pdfs": False, "science_measurement": False,
        "results": str(results_dir),
    }, ensure_ascii=False), flush=True)
    if args.check:
        return 0
    if not args.approved_source_row_export:
        raise ValueError("submission requires explicit source-row approval")
    if args.approved_plan_sha256 != plan["plan_sha256"]:
        raise ValueError("approval hash does not match the immutable source-row plan")
    if not sys.stdin.isatty():
        raise ValueError("submission requires one hidden password prompt in an interactive terminal")
    if not platform.node().lower().startswith("lageunha"):
        raise ValueError("run the approved HSC source-row export from the established lageunha environment")
    if results_dir.exists():
        checkpoint = results_dir / "row_status.json"
        if not checkpoint.is_file():
            raise FileExistsError("refusing existing results directory without its plan-bound checkpoint")
        state = json.loads(checkpoint.read_text(encoding="utf-8"))
        if state.get("plan_sha256") != plan["plan_sha256"]:
            raise ValueError("existing checkpoint belongs to another source-row plan")
    else:
        results_dir.mkdir(parents=True, exist_ok=False)
        checkpoint = results_dir / "row_status.json"
        state = {
            "schema": "astrafek-s23b-wl-common-source-row-run-v1", "status": "pending",
            "plan_sha256": plan["plan_sha256"], "count_only_plan_sha256": manifest[
                "count_only_preflight_plan_sha256"], "jobs": {}, "credentials_saved": False,
            "source_shape_rows_downloaded": False, "photoz_pdf_rows_included": False,
            "images_returned": False, "science_measurement": False,
        }
        save(checkpoint, state)
    if any(state.get(key) is not False for key in (
            "credentials_saved", "photoz_pdf_rows_included", "images_returned",
            "science_measurement")):
        raise ValueError("saved checkpoint contains data outside the approved source-row scope")

    import fcntl
    with (results_dir / "run.lock").open("a") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        password = getpass.getpass("HSC password (hidden; one prompt): ")
        if not password:
            raise ValueError("empty password")
        client = Cas(args.user, password)
        password = None
        state["status"] = "running"
        state.pop("error_type", None)
        save(checkpoint, state)
        expected_batches = []
        lens_size = 25
        for index, batch in enumerate(plan["lens_batches"]):
            chunk = plan["lens_centers"][index * lens_size:(index + 1) * lens_size]
            ids = [str(row["lens_id"]) for row in chunk]
            expected_batches.append((
                batch["name"].removesuffix(".sql"), batch, chunk,
                {lens_id: int(plan["lens_counts"][lens_id]) for lens_id in ids},
                lambda path, c=chunk, counts={key: int(plan["lens_counts"][key]) for key in ids}:
                    validate_lens_rows(path, c, counts),
            ))
        random_size = int(json.loads((plan_dir / "random_shape_count_plan.json").read_text())[
            "batch_size"])
        for index, batch in enumerate(plan["random_batches"]):
            chunk = plan["random_centers"][index * random_size:(index + 1) * random_size]
            expected = {center["random_center_id"]: plan["random_counts"][center["random_center_id"]]
                        for center in chunk}
            expected_batches.append((
                batch["name"].removesuffix(".sql"), batch, chunk, expected,
                lambda path, c=chunk, counts=expected: sum(read_shape_rows(path, c, counts).values()),
            ))
        total_rows = {"lens": 0, "random": 0}
        total_bytes = {"lens": 0, "random": 0}
        try:
            for name, batch, centers, expected_counts, validate in expected_batches:
                kind = "random" if name.startswith("random_") else "lens"
                item = state["jobs"].setdefault(name, {"status": "pending"})
                target = results_dir / f"{name}.csv"
                if item.get("status") == "complete":
                    if not target.is_file() or sha(target) != item.get("sha256"):
                        raise ValueError(f"saved completed batch is missing or changed: {name}")
                    rows = validate(target)
                    total_rows[kind] += rows
                    total_bytes[kind] += target.stat().st_size
                    continue
                if item.get("status") == "error":
                    raise ValueError(f"CAS marked this batch failed; inspect before retry: {name}")
                if item.get("status") == "submitting" and "job_id" not in item:
                    raise ValueError(f"ambiguous submit outcome; inspect CAS before retry: {name}")
                if "job_id" not in item:
                    item["status"] = "submitting"
                    save(checkpoint, state)
                    sql_path = plan_dir / batch["name"]
                    try:
                        submitted = client.json(
                            "submit", catalog_job=dict(sql=sql_path.read_text(encoding="utf-8"),
                                out_format="csv", include_metainfo_to_body=False,
                                release_version="pdr3"), nomail=True, skip_syntax_check=False,
                        )
                    except urllib.error.HTTPError:
                        item["status"] = "pending"
                        save(checkpoint, state)
                        raise
                    job_id = str(submitted.get("id", ""))
                    if not job_id:
                        item["status"] = "pending"
                        save(checkpoint, state)
                        raise ValueError(f"CAS did not return a job ID for {name}")
                    item.update(status="submitted", job_id=job_id)
                    save(checkpoint, state)
                while True:
                    status = cas_status(client, item["job_id"])
                    if status.get("status") == "done":
                        break
                    if status.get("status") == "error":
                        item.update(status="error", query_error="CAS job error; response body omitted")
                        save(checkpoint, state)
                        raise ValueError(f"CAS query failed for {name}")
                    print(f"{name}: {status.get('status', 'unknown')}", flush=True)
                    time.sleep(60)
                rows, byte_count, digest = response_rows(
                    client, name, item, plan_dir / batch["name"], results_dir, validate,
                )
                item.update(status="complete", rows=rows, bytes=byte_count, sha256=digest)
                state["source_shape_rows_downloaded"] = True
                save(checkpoint, state)
                total_rows[kind] += rows
                total_bytes[kind] += byte_count
                print(f"{name}: verified {rows:,} source pairs ({byte_count:,} bytes)", flush=True)

            expected_totals = {
                "lens": int(manifest["expected_lens_source_pair_rows_from_approved_counts"]),
                "random": int(manifest["expected_random_source_pair_rows"]),
            }
            if total_rows != expected_totals:
                raise ValueError("final source-row totals differ from the approved per-center counts")
            summary = {
                "schema": "astrafek-s23b-wl-common-source-row-summary-v1", "status": "complete",
                "plan_sha256": plan["plan_sha256"],
                "count_only_plan_sha256": manifest["count_only_preflight_plan_sha256"],
                "n_lenses": 80, "n_random_centers": 400,
                "lens_batches": 4, "random_batches": 16,
                "lens_source_pairs": total_rows["lens"],
                "random_source_pairs": total_rows["random"],
                "combined_source_pairs": sum(total_rows.values()),
                "lens_csv_bytes": total_bytes["lens"],
                "random_csv_bytes": total_bytes["random"],
                "combined_csv_bytes": sum(total_bytes.values()),
                "source_shape_rows_downloaded": True, "photoz_pdf_rows_included": False,
                "images_returned": False, "credentials_saved": False,
                "science_measurement": False,
            }
            state.update(status="complete", summary=summary,
                         source_shape_rows_downloaded=True)
            save(checkpoint, state)
            save(results_dir / "row_summary.json", summary)
        except BaseException as error:
            state.update(status="interrupted" if isinstance(error, KeyboardInterrupt) else "stopped_on_error",
                         error_type=type(error).__name__)
            save(checkpoint, state)
            raise
        finally:
            client.credential.clear()
    print(json.dumps(summary, ensure_ascii=False, indent=2), flush=True)
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan-dir", type=Path, default=DEFAULT_PLAN_DIR)
    parser.add_argument("--results", type=Path, default=DEFAULT_RESULTS)
    parser.add_argument("--user", default=DEFAULT_USER)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--check", action="store_true", help="validate all plan inputs without network")
    group.add_argument("--submit", action="store_true", help="submit and retrieve approved source rows")
    parser.add_argument("--approved-source-row-export", action="store_true",
                        help="required after explicit approval of this exact data export")
    parser.add_argument("--approved-plan-sha256", help="exact immutable plan file hash approved by the user")
    args = parser.parse_args()
    if args.check and (args.approved_source_row_export or args.approved_plan_sha256):
        parser.error("approval flags are used only with --submit")
    if args.submit and (not args.approved_source_row_export or not args.approved_plan_sha256):
        parser.error("--submit requires explicit approval flag and exact plan SHA-256")
    return run(args)


if __name__ == "__main__":
    try:
        sys.exit(main())
    except urllib.error.HTTPError as error:
        print(f"HSC CAS HTTP {error.code}; credentials and response body suppressed.", file=sys.stderr)
        sys.exit(1)
    except KeyboardInterrupt:
        print("\nStopped; submitted job IDs and verified batches remain checkpointed.", file=sys.stderr)
        sys.exit(130)
    except Exception as error:
        message = str(error) if isinstance(error, ValueError) else type(error).__name__
        print(f"submit_s23b_wl_common_shape_rows: {message}", file=sys.stderr)
        sys.exit(1)
