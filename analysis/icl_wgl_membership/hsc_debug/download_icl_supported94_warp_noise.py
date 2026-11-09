#!/usr/bin/env python3
"""Download an explicitly approved, hash-pinned 94-pair HSC warp plan.

Default mode is read-only plan inspection. Network access requires --submit and
the exact approved plan SHA-256. One request runs at a time. Per-request and
total byte caps are enforced while streaming; checkpoints support resumption.
Credentials are used in memory only and are never included in output files.
"""
from __future__ import annotations

import argparse
import getpass
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import tempfile
import time
import urllib.error
import urllib.request


HERE = Path(__file__).resolve().parent
PLAN_DIR = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_supported94_warp_noise_plan_2026-11-06_01")
DEFAULT_OUTPUT = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_supported94_warp_noise_download_2026-11-06_01")
APPROVED_PLAN_SHA256 = "bd936a94bcdbff83a607acbfca66616c396db526ada50cc89c9e3e312d37007a"
MAX_REQUEST_BYTES = 160 * 1024 * 1024
MAX_TOTAL_BYTES = 12 * 1024**3
CHUNK_BYTES = 4 * 1024 * 1024


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def json_compatible(value):
    """Convert FITS/Astropy NumPy scalars before checkpoint serialization."""
    if isinstance(value, dict):
        return {str(key): json_compatible(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_compatible(item) for item in value]
    if type(value).__module__.startswith("numpy") and hasattr(value, "item"):
        return value.item()
    return value


def atomic_json(path: Path, value: dict) -> None:
    temp = path.with_name(path.name + ".partial")
    temp.write_text(json.dumps(json_compatible(value), ensure_ascii=False, indent=2, sort_keys=True,
                               allow_nan=False) + "\n", encoding="utf-8")
    temp.replace(path)


def load_pilot_helpers():
    helper_path = HERE / "download_icl_warp_pilot.py"
    spec = importlib.util.spec_from_file_location("astrafek_hsc_warp_pilot_helpers", helper_path)
    if spec is None or spec.loader is None:
        raise RuntimeError("cannot load the audited HSC warp helper")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def load_plan(plan_dir: Path) -> tuple[dict, Path, str]:
    plan_path = Path(plan_dir).resolve() / "plan.json"
    sidecar = plan_path.with_name("plan.sha256")
    if not plan_path.is_file() or not sidecar.is_file():
        raise FileNotFoundError("plan.json or plan.sha256 is missing")
    digest = sha256(plan_path)
    recorded = sidecar.read_text(encoding="ascii").split()[0]
    if digest != recorded or digest != APPROVED_PLAN_SHA256:
        raise ValueError(f"plan does not match the exact reviewable artifact: {digest}")
    plan = json.loads(plan_path.read_text(encoding="utf-8"))
    if (plan.get("schema") != "astrafek-icl-supported94-warp-noise-plan-v1"
            or plan.get("status") != "prepared_not_submitted"
            or plan.get("execution_requires_separate_explicit_user_approval") is not True
            or plan.get("query_submitted") is not False
            or plan.get("network_used") is not False
            or plan.get("science_measurement") is not False
            or plan.get("scope", {}).get("n_pairs") != 94
            or plan.get("scope", {}).get("n_center_filter_requests") != 188
            or plan.get("scope", {}).get("n_warp_visit_products_expected") != 1471
            or plan.get("estimated_payload", {}).get("total_hard_cap_bytes") != MAX_TOTAL_BYTES
            or plan.get("estimated_payload", {}).get("per_request_hard_cap_bytes") != MAX_REQUEST_BYTES):
        raise ValueError("plan scope, guard state, or byte caps differ from the reviewed request")
    for path_text, expected_hash in plan.get("input_sha256", {}).items():
        path = Path(path_text)
        if not path.is_file() or sha256(path) != expected_hash:
            raise ValueError(f"pinned input changed since planning: {path}")
    return plan, plan_path, digest


def inspect_plan(plan: dict, plan_hash: str) -> dict:
    return {
        "status": "prepared_not_submitted",
        "query_submitted": False,
        "plan_sha256": plan_hash,
        "pairs": plan["scope"]["n_pairs"],
        "center_requests": plan["scope"]["n_center_filter_requests"],
        "expected_visit_products": plan["scope"]["n_warp_visit_products_expected"],
        "estimated_gib": plan["estimated_payload"]["estimated_gib"],
        "per_request_cap_mib": plan["estimated_payload"]["per_request_hard_cap_mib"],
        "total_cap_gib": plan["estimated_payload"]["total_hard_cap_gib"],
        "next_action": "network submission requires explicit user approval for this exact plan hash",
    }


def stream_archive(url: str, output_dir: Path, user: str, password: str,
                   opener_factory, *, remaining_bytes: int) -> tuple[Path, int, str]:
    cap = min(MAX_REQUEST_BYTES, remaining_bytes)
    if cap <= 0:
        raise ValueError("total approved byte budget is exhausted")
    request = urllib.request.Request(url, headers={"User-Agent": "Astrafek/ICL-supported94-warp-noise"})
    opener = opener_factory(user, password)
    temp_path = None
    try:
        with opener.open(request, timeout=180) as response:
            status_code = int(response.status)
            content_type = response.headers.get("Content-Type", "")
            content_length = response.headers.get("Content-Length")
            if content_length is not None and int(content_length) > cap:
                raise ValueError(f"response Content-Length {content_length} exceeds remaining approved cap {cap}")
            handle = tempfile.NamedTemporaryFile(
                prefix=".icl-supported94-warp-", suffix=".partial", dir=output_dir, delete=False,
            )
            temp_path = Path(handle.name)
            total = 0
            with handle:
                while True:
                    chunk = response.read(CHUNK_BYTES)
                    if not chunk:
                        break
                    total += len(chunk)
                    if total > cap:
                        raise ValueError(f"stream exceeded remaining approved cap {cap}")
                    handle.write(chunk)
        if status_code != 200:
            raise ValueError(f"HSC endpoint returned HTTP {status_code}")
        return temp_path, total, content_type
    except urllib.error.HTTPError as exc:
        if temp_path and temp_path.exists():
            temp_path.unlink()
        raise RuntimeError(f"HSC endpoint returned HTTP {exc.code}") from None
    except Exception:
        if temp_path and temp_path.exists():
            temp_path.unlink()
        raise


def audit_one(archive_path: Path, request: dict, helper) -> dict:
    audit = helper.inspect_tar(
        archive_path, set(request["expected_visits"]),
        expected_filter=request["filter"], expected_tract=int(request["tract"]),
        expected_rerun=request["rerun"], center_ra=float(request["ra_deg"]),
        center_dec=float(request["dec_deg"]),
    )
    audit["patch_expected_from_metadata"] = int(request["patch"])
    expected = set(request["expected_visits"])
    actual = set(audit["actual_visits"])
    audit["missing_expected_visits"] = sorted(expected - actual)
    audit["additional_cutout_visits"] = sorted(actual - expected)
    audit["all_metadata_expected_visits_returned"] = expected.issubset(actual)
    audit["visit_inventory_note"] = (
        "exact metadata visit set" if actual == expected else
        "all metadata visits returned; additional warp products came with this same spatial cutout"
        if expected.issubset(actual) else
        "one or more metadata-expected visits were absent"
    )
    audit["expected_patch_header_match"] = all(
        member.get("header_patch") is None
        or str(member.get("header_patch")) == str(request["patch"])
        for member in audit["members"]
    )
    audit["all_member_layout_wcs_and_variance_checks_pass"] = bool(
        audit["all_member_layout_wcs_and_variance_checks_pass"]
        and audit["expected_patch_header_match"]
    )
    return audit


def run(plan_dir: Path, output_dir: Path, *, submit: bool,
        confirmed_plan_sha: str | None, user_override: str | None = None) -> dict:
    plan, plan_path, plan_hash = load_plan(plan_dir)
    if not submit:
        payload = inspect_plan(plan, plan_hash)
        print(json.dumps(payload, ensure_ascii=False, indent=2))
        return payload
    if confirmed_plan_sha != plan_hash:
        raise PermissionError("submission requires explicit approval of the exact plan SHA-256")
    output_dir = Path(output_dir).resolve()
    state_path = output_dir / "checkpoint.json"
    summary_path = output_dir / "summary.json"
    if not output_dir.exists():
        output_dir.mkdir(parents=True, exist_ok=False)
    elif not state_path.is_file() or summary_path.exists():
        raise FileExistsError(f"existing output is not a resumable checkpoint: {output_dir}")

    if state_path.exists():
        checkpoint = json.loads(state_path.read_text(encoding="utf-8"))
        if checkpoint.get("plan_sha256") != plan_hash or checkpoint.get("science_measurement") is not False:
            raise ValueError("checkpoint belongs to a different plan or guard state")
    else:
        checkpoint = {
            "schema": "astrafek-icl-supported94-warp-noise-checkpoint-v1",
            "status": "running",
            "plan_sha256": plan_hash,
            "approved_plan_sha256": confirmed_plan_sha,
            "query_submitted": True,
            "network_used": True,
            "credentials_saved": False,
            "science_measurement": False,
            "total_downloaded_bytes": 0,
            "completed_request_count": 0,
            "requests": {},
        }
        atomic_json(state_path, checkpoint)

    downloaded_so_far = int(checkpoint.get("total_downloaded_bytes", 0))
    if downloaded_so_far > MAX_TOTAL_BYTES:
        raise ValueError("checkpoint already exceeds the plan's total byte ceiling")
    helper = load_pilot_helpers()
    auth = helper.load_existing_cutout_helpers()

    # Recover an archive that reached disk just before a previous checkpoint
    # write failed. This avoids an unnecessary duplicate request on resume.
    for index, request_spec in enumerate(plan["requests"]):
        key = f"request_{index:03d}"
        archive_path = output_dir / f"{key}.tar"
        entry = checkpoint["requests"].get(key)
        if entry and entry.get("status") in ("complete", "complete_with_qa_warnings"):
            if not archive_path.is_file() or sha256(archive_path) != entry.get("archive_sha256"):
                raise ValueError(f"completed request archive failed checkpoint hash: {archive_path}")
            continue
        if not archive_path.is_file():
            continue
        nbytes = archive_path.stat().st_size
        if nbytes > MAX_REQUEST_BYTES or downloaded_so_far + nbytes > MAX_TOTAL_BYTES:
            raise ValueError(f"orphan archive exceeds an approved byte cap: {archive_path}")
        audit = audit_one(archive_path, request_spec, helper)
        if not audit["all_metadata_expected_visits_returned"]:
            request_status = "complete_scope_mismatch"
        elif (not audit["all_member_layout_wcs_and_variance_checks_pass"]
              or audit["additional_cutout_visits"]):
            request_status = "complete_with_qa_warnings"
        else:
            request_status = "complete"
        checkpoint["requests"][key] = {
            "request_id": request_spec["request_id"],
            "status": request_status,
            "archive": str(archive_path),
            "archive_bytes": nbytes,
            "archive_sha256": sha256(archive_path),
            "content_type": None,
            "recovered_orphan_archive_without_checkpoint": True,
            "expected_visits": request_spec["expected_visits"],
            "audit": audit,
        }
        downloaded_so_far += nbytes
        checkpoint["total_downloaded_bytes"] = downloaded_so_far
        checkpoint["completed_request_count"] = len(checkpoint["requests"])
        checkpoint["last_request_id"] = request_spec["request_id"]
        atomic_json(state_path, checkpoint)
        if request_status == "complete_scope_mismatch":
            checkpoint["status"] = "stopped_scope_mismatch"
            atomic_json(state_path, checkpoint)
            raise ValueError(f"orphan response did not match the approved visit set: {key}")

    user = user_override or os.environ.get("HSC_USER")
    password = os.environ.get("HSC_PASSWORD")
    has_incomplete = any(
        checkpoint["requests"].get(f"request_{index:03d}", {}).get("status")
        not in ("complete", "complete_with_qa_warnings", "http_404_permanent")
        for index in range(len(plan["requests"]))
    )
    if has_incomplete:
        user = user or input("HSC account name: ").strip()
        password = password or getpass.getpass("HSC password (hidden): ")
        if not user or not password:
            raise ValueError("HSC credentials are required for the approved download")

    for index, request_spec in enumerate(plan["requests"]):
        key = f"request_{index:03d}"
        entry = checkpoint["requests"].get(key)
        archive_path = output_dir / f"{key}.tar"
        if entry and entry.get("status") in ("complete", "complete_with_qa_warnings"):
            if not archive_path.is_file() or sha256(archive_path) != entry.get("archive_sha256"):
                raise ValueError(f"completed request archive failed checkpoint hash: {archive_path}")
            continue
        if entry and entry.get("status") == "http_404_permanent":
            continue
        if entry and entry.get("status") == "complete_scope_mismatch":
            raise ValueError(f"prior response for {key} had a scope mismatch; stopping for review")
        if archive_path.is_file():
            raise ValueError(f"unreconciled archive remains after recovery: {archive_path}")

        remaining = MAX_TOTAL_BYTES - downloaded_so_far
        url = helper.query_url(auth.BASE, request_spec)
        try:
            partial, nbytes, content_type = stream_archive(
                url, output_dir, user, password, auth.opener_for,
                remaining_bytes=remaining,
            )
        except RuntimeError as exc:
            if "HTTP 404" not in str(exc):
                raise
            prior_attempts = int(checkpoint["requests"].get(key, {}).get("attempt_count", 0))
            attempt_count = prior_attempts + 1
            checkpoint["requests"][key] = {
                "request_id": request_spec["request_id"],
                "status": "http_404_permanent" if attempt_count >= 2 else "http_404_retryable",
                "http_status": 404,
                "attempt_count": attempt_count,
                "expected_visits": request_spec["expected_visits"],
                "error": "HSC cutout endpoint returned HTTP 404 for this approved spatial request",
            }
            checkpoint["completed_request_count"] = len(checkpoint["requests"])
            checkpoint["last_request_id"] = request_spec["request_id"]
            atomic_json(state_path, checkpoint)
            if attempt_count == 1:
                time.sleep(1.0)
                try:
                    partial, nbytes, content_type = stream_archive(
                        url, output_dir, user, password, auth.opener_for,
                        remaining_bytes=MAX_TOTAL_BYTES - downloaded_so_far,
                    )
                except RuntimeError as retry_exc:
                    if "HTTP 404" not in str(retry_exc):
                        raise
                    checkpoint["requests"][key]["status"] = "http_404_permanent"
                    checkpoint["requests"][key]["attempt_count"] = 2
                    atomic_json(state_path, checkpoint)
                    continue
            else:
                continue
        partial.replace(archive_path)
        audit = audit_one(archive_path, request_spec, helper)
        if not audit["all_metadata_expected_visits_returned"]:
            request_status = "complete_scope_mismatch"
        elif (not audit["all_member_layout_wcs_and_variance_checks_pass"]
              or audit["additional_cutout_visits"]):
            request_status = "complete_with_qa_warnings"
        else:
            request_status = "complete"
        checkpoint["requests"][key] = {
            "request_id": request_spec["request_id"],
            "status": request_status,
            "archive": str(archive_path),
            "archive_bytes": nbytes,
            "archive_sha256": sha256(archive_path),
            "content_type": content_type,
            "expected_visits": request_spec["expected_visits"],
            "audit": audit,
        }
        downloaded_so_far += nbytes
        checkpoint["total_downloaded_bytes"] = downloaded_so_far
        checkpoint["completed_request_count"] = len(checkpoint["requests"])
        checkpoint["last_request_id"] = request_spec["request_id"]
        checkpoint["last_request_timestamp_utc"] = time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())
        atomic_json(state_path, checkpoint)
        if request_status == "complete_scope_mismatch":
            checkpoint["status"] = "stopped_scope_mismatch"
            atomic_json(state_path, checkpoint)
            raise ValueError(f"HSC response did not match the approved visit set for {key}; checkpoint saved")
        time.sleep(0.25)

    all_rows = list(checkpoint["requests"].values())
    n_archive_complete = sum(row.get("status") in ("complete", "complete_with_qa_warnings")
                             for row in all_rows)
    n_http_404 = sum(row.get("status") == "http_404_permanent" for row in all_rows)
    status = "complete_with_data_gaps" if n_http_404 else "complete"
    if not n_http_404 and any(row["status"] == "complete_with_qa_warnings" for row in all_rows):
        status = "complete_with_qa_warnings"
    summary = {
        "schema": "astrafek-icl-supported94-warp-noise-download-v1",
        "status": status,
        "plan_sha256": plan_hash,
        "query_submitted": True,
        "network_used": True,
        "credentials_saved": False,
        "request_count": len(all_rows),
        "completed_archive_count": n_archive_complete,
        "http_404_request_count": n_http_404,
        "http_404_request_ids": [row["request_id"] for row in all_rows
                                 if row.get("status") == "http_404_permanent"],
        "expected_request_count": len(plan["requests"]),
        "expected_visit_product_count": plan["scope"]["n_warp_visit_products_expected"],
        "downloaded_bytes": downloaded_so_far,
        "total_byte_cap": MAX_TOTAL_BYTES,
        "science_measurement": False,
        "requests": all_rows,
    }
    atomic_json(summary_path, summary)
    checkpoint["status"] = status
    atomic_json(state_path, checkpoint)
    print(json.dumps({"status": status, "requests": len(all_rows),
                      "downloaded_bytes": downloaded_so_far,
                      "science_measurement": False,
                      "summary": str(summary_path)}, ensure_ascii=False, indent=2))
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan-dir", type=Path, default=PLAN_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--submit", action="store_true", help="submit the network requests")
    parser.add_argument("--confirm-plan-sha", help="exact separately approved plan SHA-256")
    parser.add_argument("--user", help="HSC user name; password is prompted securely unless HSC_PASSWORD is already set")
    args = parser.parse_args()
    run(args.plan_dir, args.output_dir, submit=args.submit,
        confirmed_plan_sha=args.confirm_plan_sha, user_override=args.user)


if __name__ == "__main__":
    main()
