#!/usr/bin/env python3
"""Prepare an approval-gated HSC warp transfer plan for the supported ICL set.

This planner reads only pinned local Astrafek summaries, the approved frame
metadata export, and the existing 94-pair ICL/satellite table. It does not
read credentials or access HSC. Its plan covers the target and matched control
positions of every supported pair and is not itself authorization to download.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import defaultdict
from pathlib import Path


ROOT = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
PAIRS = ROOT / "icl_satellite_supported_94_association_2026-11-05_01/supported_pairs.tsv"
FRAME_RESULTS = ROOT / "icl_frame_metadata_row_results_2026-11-06_01"
FRAME_AUDIT = ROOT / "icl_frame_metadata_row_audit_2026-11-06_01/summary.json"
PILOT_SUMMARY = ROOT / "icl_warp_pilot_coadd_variance_2026-11-06_02/summary.json"
PILOT_ARCHIVE = ROOT / "icl_warp_pilot_download_2026-11-06_01/warp_pilot.tar"
DEFAULT_OUTPUT = ROOT / "icl_supported94_warp_noise_plan_2026-11-06_01"

EXPECTED = {
    "pairs": "79f7c4c6b18b29ad837e91da4a4359433b09a564cad62efaf60b4cce655685da",
    "frame_summary": "3925f74be34b19cf16ea0183055ba4d74d91568f66150c569e2d48b6bd2f7b53",
    "frame_status": "ea1e6611390fb05bea79fde928f43c8ac9fd21bcaa67a402b8d3e1abe22fb173",
    "frame_audit": "0f361fcb496fbd7be7c1059f638414f6900535f73a9ac48696030480f400a237",
    "pilot_summary": "ad523613f9efdd5226759d13bef1973cb2add1990614a80b96b27adca0753dae",
    "pilot_archive": "cdc9b642c5c5d0361b58a0bf7c0e733a60cb78caad5ed496ae00cac89b23953c",
}

PIXEL_SCALE_ARCSEC = 0.168
FULL_SIZE_ARCSEC = 120.0
PILOT_ARCHIVE_BYTES = 55_511_040
PILOT_VISITS = 9
MAX_REQUEST_BYTES = 160 * 1024 * 1024
MAX_TOTAL_BYTES = 12 * 1024**3


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def read_metadata(path: Path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        header = stream.readline().lstrip("# ").rstrip("\r\n")
        columns = next(csv.reader([header]))
        yield from csv.DictReader(stream, fieldnames=columns)


def load_inputs() -> tuple[list[dict], dict, dict]:
    required = [PAIRS, FRAME_RESULTS / "summary.json", FRAME_RESULTS / "status.json",
                FRAME_AUDIT, PILOT_SUMMARY, PILOT_ARCHIVE]
    for path in required:
        if not path.is_file():
            raise FileNotFoundError(path)
    hashes = {str(path): sha256(path) for path in required}
    expected = {
        str(PAIRS): EXPECTED["pairs"],
        str(FRAME_RESULTS / "summary.json"): EXPECTED["frame_summary"],
        str(FRAME_RESULTS / "status.json"): EXPECTED["frame_status"],
        str(FRAME_AUDIT): EXPECTED["frame_audit"],
        str(PILOT_SUMMARY): EXPECTED["pilot_summary"],
        str(PILOT_ARCHIVE): EXPECTED["pilot_archive"],
    }
    if hashes != expected:
        raise ValueError(f"pinned planning input changed: {hashes}")

    frame_summary = json.loads((FRAME_RESULTS / "summary.json").read_text(encoding="utf-8"))
    frame_status = json.loads((FRAME_RESULTS / "status.json").read_text(encoding="utf-8"))
    frame_audit = json.loads(FRAME_AUDIT.read_text(encoding="utf-8"))
    pilot = json.loads(PILOT_SUMMARY.read_text(encoding="utf-8"))
    if (frame_summary.get("status") != "complete" or frame_summary.get("n_frame_rows") != 67071
            or frame_summary.get("n_centers") != 1598
            or frame_summary.get("images_returned") is not False
            or frame_summary.get("shape_rows_returned") is not False
            or frame_summary.get("photoz_rows_returned") is not False
            or frame_summary.get("science_measurement") is not False):
        raise ValueError("approved frame metadata export has unexpected scope or guard state")
    if (frame_status.get("status") != "complete" or len(frame_status.get("jobs", {})) != 16
            or frame_audit.get("status") != "complete_exploratory"
            or frame_audit.get("science_measurement") is not False):
        raise ValueError("frame metadata export or independent audit is incomplete")
    if (pilot.get("status") != "complete_exploratory"
            or pilot.get("science_measurement") is not False
            or pilot.get("warp_archive", {}).get("n_visits") != PILOT_VISITS
            or pilot.get("warp_archive", {}).get("bytes") != PILOT_ARCHIVE_BYTES):
        raise ValueError("single-position warp pilot is incomplete or differs from the size reference")

    for batch in range(16):
        job_id = f"frame_metadata_rows_{batch:03d}"
        job = frame_status["jobs"].get(job_id, {})
        path = FRAME_RESULTS / f"{job_id}.csv"
        if job.get("status") != "complete" or not path.is_file() or sha256(path) != job.get("sha256"):
            raise ValueError(f"approved frame metadata batch failed hash/status validation: {path}")

    pairs = list(csv.DictReader(PAIRS.open(encoding="utf-8"), delimiter="\t"))
    if (len(pairs) != 94 or {band: sum(row["filter"] == band for row in pairs)
                             for band in ("i", "i2")} != {"i": 26, "i2": 68}
            or any(row.get("science_measurement", "").lower() != "false" for row in pairs)):
        raise ValueError("the supported 94-pair table has unexpected rows, filter balance, or guard state")
    return pairs, hashes, frame_status


def run(output: Path) -> dict:
    output = Path(output).resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite plan directory: {output}")
    pairs, input_hashes, status = load_inputs()

    request_groups: dict[tuple[str, str], dict] = {}
    pair_records = []
    for pair in pairs:
        band = pair["filter"]
        cluster_id = str(pair["cluster_id"])
        target_id = f"target:{cluster_id}"
        control_id = f"control:{pair['control_id']}"
        pair_records.append({
            "cluster_id": cluster_id,
            "filter_code": band,
            "connected_5deg_block": pair["connected_5deg_block"],
            "target_request_id": f"{band}:{target_id}",
            "control_request_id": f"{band}:{control_id}",
        })
        for role, center_id, ra_key, dec_key in (
            ("target", target_id, "target_ra_deg", "target_dec_deg"),
            ("control", control_id, "control_ra_deg", "control_dec_deg"),
        ):
            key = (band, center_id)
            candidate = {
                "request_id": f"{band}:{center_id}",
                "center_id": center_id,
                "role": role,
                "filter_code": band,
                "ra_deg": float(pair[ra_key]),
                "dec_deg": float(pair[dec_key]),
                "pair_ids": [cluster_id],
            }
            if key in request_groups:
                previous = request_groups[key]
                if (abs(previous["ra_deg"] - candidate["ra_deg"]) > 1e-10
                        or abs(previous["dec_deg"] - candidate["dec_deg"]) > 1e-10
                        or previous["role"] != role):
                    raise ValueError(f"reused center has inconsistent coordinates/role: {key}")
                previous["pair_ids"].append(cluster_id)
            else:
                request_groups[key] = candidate

    selected = set(request_groups)
    metadata_by_center: dict[tuple[str, str], list[dict]] = defaultdict(list)
    metadata_provenance = defaultdict(set)
    for batch in range(16):
        name = f"frame_metadata_rows_{batch:03d}.csv"
        path = FRAME_RESULTS / name
        for row in read_metadata(path):
            key = (row["filter"], row["center_id"])
            if key not in selected:
                continue
            metadata_by_center[key].append(row)
            metadata_provenance[key].add(name)

    if set(metadata_by_center) != selected:
        missing = sorted(selected - set(metadata_by_center))
        raise ValueError(f"selected target/control centers missing from approved metadata: {missing[:10]}")

    for key, request in request_groups.items():
        rows = metadata_by_center[key]
        expected_role = "target" if request["role"] == "target" else "control"
        if any(row["sample_role"] != expected_role or row["coadd_band"] != "HSC-I" for row in rows):
            raise ValueError(f"center role or coadd band does not match approved metadata: {key}")
        for row in rows:
            if (abs(float(row["ra_deg"]) - request["ra_deg"]) > 1e-6
                    or abs(float(row["dec_deg"]) - request["dec_deg"]) > 1e-6):
                raise ValueError(f"pair-table and approved frame-metadata coordinates disagree: {key}")
        tracts = {int(row["tract"]) for row in rows}
        patches = {int(row["patch"]) for row in rows}
        if len(tracts) != 1 or len(patches) != 1:
            raise ValueError(f"a center needs multiple tract/patch requests; refusing ambiguous scope: {key}")
        visits = sorted({int(row["visit"]) for row in rows})
        frames = {str(row["frame_num"]) for row in rows}
        request.update({
            "rerun": "pdr3_wide",
            "type": "warp",
            "filter": "HSC-I",
            "tract": tracts.pop(),
            "patch": patches.pop(),
            "semi_width": "1arcmin",
            "semi_height": "1arcmin",
            "layers": {"image": True, "mask": True, "variance": True},
            "expected_visits": visits,
            "n_expected_visits": len(visits),
            "n_frame_rows": len(rows),
            "n_frames": len(frames),
            "metadata_batches": sorted(metadata_provenance[key]),
        })

    per_visit_bytes = PILOT_ARCHIVE_BYTES / PILOT_VISITS
    estimated_bytes = int(math.ceil(per_visit_bytes * sum(r["n_expected_visits"] for r in request_groups.values())))
    if estimated_bytes >= MAX_TOTAL_BYTES:
        raise ValueError("pilot-scaled estimate exceeds the proposed 12 GiB hard transfer ceiling")
    max_visits = max(row["n_expected_visits"] for row in request_groups.values())
    max_group_estimate = per_visit_bytes * max_visits
    if max_group_estimate >= MAX_REQUEST_BYTES:
        raise ValueError("the maximum per-request estimate exceeds the proposed 160 MiB ceiling")

    visits_by_band = {
        band: sum(row["n_expected_visits"] for row in request_groups.values()
                  if row["filter_code"] == band)
        for band in ("i", "i2")
    }
    center_by_band = {
        band: sum(row["filter_code"] == band for row in request_groups.values())
        for band in ("i", "i2")
    }
    plan = {
        "schema": "astrafek-icl-supported94-warp-noise-plan-v1",
        "status": "prepared_not_submitted",
        "execution_requires_separate_explicit_user_approval": True,
        "query_submitted": False,
        "network_used": False,
        "credentials_read_or_saved": False,
        "science_measurement": False,
        "purpose": "validate delivered coadd variance against visit warp variance across every target/control center in the existing ICL–satellite common-support sample",
        "scope": {
            "supported_pair_table": str(PAIRS),
            "n_pairs": len(pair_records),
            "pairs_by_filter": {"i": 26, "i2": 68},
            "n_center_filter_requests": len(request_groups),
            "center_filter_requests_by_band": center_by_band,
            "n_warp_visit_products_expected": sum(r["n_expected_visits"] for r in request_groups.values()),
            "warp_visit_products_by_band": visits_by_band,
            "unique_visit_counts_per_request": {
                "minimum": min(r["n_expected_visits"] for r in request_groups.values()),
                "median": sorted(r["n_expected_visits"] for r in request_groups.values())[len(request_groups)//2],
                "maximum": max_visits,
            },
            "layers": {"image": True, "mask": True, "variance": True},
            "cutout_size_arcsec": [FULL_SIZE_ARCSEC, FULL_SIZE_ARCSEC],
            "rerun": "pdr3_wide",
            "cutout_filter": "HSC-I",
        },
        "estimated_payload": {
            "pilot_archive_bytes": PILOT_ARCHIVE_BYTES,
            "pilot_visit_count": PILOT_VISITS,
            "measured_bytes_per_returned_visit": per_visit_bytes,
            "estimated_bytes": estimated_bytes,
            "estimated_gib": estimated_bytes / 1024**3,
            "estimate_method": "linear scaling by expected visit count using the measured nine-visit 2x2-arcmin warp pilot; actual FITS/tar payloads vary by position",
            "per_request_hard_cap_bytes": MAX_REQUEST_BYTES,
            "per_request_hard_cap_mib": MAX_REQUEST_BYTES / 1024**2,
            "total_hard_cap_bytes": MAX_TOTAL_BYTES,
            "total_hard_cap_gib": MAX_TOTAL_BYTES / 1024**3,
            "maximum_group_pilot_scaled_estimate_bytes": int(math.ceil(max_group_estimate)),
        },
        "not_requested": ["catalogue rows", "shape rows", "photo-z", "PSF Picker models", "other filters", "coadd images outside these cutouts"],
        "after_receipt_checks": [
            "validate each request archive and exact expected visit set",
            "validate FITS image/mask/variance planes, mask definitions, WCS and pixel scale",
            "compare warp inverse-variance reconstruction with the paired delivered coadd VARIANCE",
            "estimate spatially blocked uncertainty across target/control centers without treating pixels as independent",
            "retain science_measurement=false until source-Poisson, sky, PSF and total covariance gates pass",
        ],
        "pairs": pair_records,
        "requests": sorted(request_groups.values(), key=lambda r: r["request_id"]),
        "input_sha256": input_hashes,
        "pilot_archive_sha256": EXPECTED["pilot_archive"],
    }
    output.mkdir(parents=True, exist_ok=False)
    plan_path = output / "plan.json"
    plan_path.write_text(json.dumps(plan, ensure_ascii=False, indent=2, allow_nan=False) + "\n",
                         encoding="utf-8")
    digest = sha256(plan_path)
    (output / "plan.sha256").write_text(f"{digest}  plan.json\n", encoding="ascii")
    print(json.dumps({
        "status": plan["status"], "query_submitted": False,
        "pairs": len(pair_records), "center_requests": len(request_groups),
        "expected_visit_products": plan["scope"]["n_warp_visit_products_expected"],
        "estimated_gib": plan["estimated_payload"]["estimated_gib"],
        "total_hard_cap_gib": plan["estimated_payload"]["total_hard_cap_gib"],
        "plan": str(plan_path), "sha256": digest,
    }, ensure_ascii=False, indent=2))
    return plan


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    run(args.output)


if __name__ == "__main__":
    main()
