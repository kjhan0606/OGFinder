#!/usr/bin/env python3
"""Prepare, but never submit, a one-center HSC PDR3 warp variance pilot.

The plan is pinned to the previously approved local frame-metadata export.
No credentials are read and no network request is made. A separate explicit
user approval is required before any warp pixels are downloaded.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

ROOT = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
RESULTS = ROOT / "icl_frame_metadata_row_results_2026-11-06_01"
AUDIT = ROOT / "icl_frame_metadata_row_audit_2026-11-06_01"
MANUAL = ROOT.parent / "logs/cutout_global_manual.html"
DEFAULT_OUTPUT = ROOT / "icl_warp_pilot_plan_2026-11-06_02"

EXPECTED_RESULTS_SHA256 = "3925f74be34b19cf16ea0183055ba4d74d91568f66150c569e2d48b6bd2f7b53"
EXPECTED_AUDIT_SHA256 = "0f361fcb496fbd7be7c1059f638414f6900535f73a9ac48696030480f400a237"
EXPECTED_MANUAL_SHA256 = "0022c523a1ea5682904c4fd73ed060a79d7d9f2317536ee55185438f1a70d452"
PILOT_CENTER_ID = "target:6966"
PIXEL_SCALE_ARCSEC = 0.168
FULL_SIZE_ARCSEC = 120.0
BYTES_PER_LAYER_PIXEL = 4
N_LAYERS = 3


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def read_result_rows(path: Path):
    with path.open(encoding="utf-8-sig", newline="") as stream:
        header = stream.readline().lstrip("# ").rstrip("\r\n")
        columns = next(csv.reader([header]))
        yield from csv.DictReader(stream, fieldnames=columns)


def run(output: Path) -> dict:
    output = Path(output).resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite warp pilot plan: {output}")
    results_summary = RESULTS / "summary.json"
    audit_summary = AUDIT / "summary.json"
    if sha256(results_summary) != EXPECTED_RESULTS_SHA256:
        raise ValueError("approved frame-metadata result changed")
    if sha256(audit_summary) != EXPECTED_AUDIT_SHA256:
        raise ValueError("independent frame-metadata audit changed")
    if sha256(MANUAL) != EXPECTED_MANUAL_SHA256:
        raise ValueError("cached PDR3 cutout manual changed")
    cas = json.loads(results_summary.read_text(encoding="utf-8"))
    audit = json.loads(audit_summary.read_text(encoding="utf-8"))
    if (cas.get("status") != "complete" or cas.get("n_frame_rows") != 67071
            or cas.get("images_returned") is not False
            or cas.get("shape_rows_returned") is not False
            or cas.get("photoz_rows_returned") is not False
            or cas.get("science_measurement") is not False):
        raise ValueError("parent export is not the expected metadata-only result")
    if audit.get("status") != "complete_exploratory" or audit.get("science_measurement") is not False:
        raise ValueError("independent metadata audit is incomplete")

    status_path = RESULTS / "status.json"
    status = json.loads(status_path.read_text(encoding="utf-8"))
    if status.get("status") != "complete" or len(status.get("jobs", {})) != 16:
        raise ValueError("the approved 16-batch metadata export is incomplete")
    center_rows = []
    visits = set()
    frames = set()
    ccds = set()
    for batch in range(16):
        job_key = f"frame_metadata_rows_{batch:03d}"
        name = f"{job_key}.csv"
        path = RESULTS / name
        job = status["jobs"].get(job_key, {})
        if job.get("status") != "complete" or not path.is_file() or sha256(path) != job.get("sha256"):
            raise ValueError(f"approved metadata batch is missing or changed: {name}")
        for candidate in read_result_rows(path):
            if candidate.get("center_id") == PILOT_CENTER_ID:
                center_rows.append(candidate)
                frames.add(str(candidate["frame_num"]))
                visits.add(int(candidate["visit"]))
                ccds.add(int(candidate["ccd"]))
    if not center_rows:
        raise ValueError("pilot center is missing from approved metadata")
    row = center_rows[0]
    if (len(center_rows), len(frames), len(visits), len(ccds)) != (43, 43, 9, 40):
        raise ValueError(f"pilot metadata scope changed: frames/visits/CCDs={len(frames)}/{len(visits)}/{len(ccds)}")
    expected_coords = (346.791503, 1.039977)
    if (row["sample_role"] != "target" or row["coadd_band"] != "HSC-I"
            or row["tract"] != "9703" or row["patch"] != "306"
            or abs(float(row["ra_deg"]) - expected_coords[0]) > 1e-8
            or abs(float(row["dec_deg"]) - expected_coords[1]) > 1e-8):
        raise ValueError("pilot center no longer matches the reviewed target and patch")

    pixels = math.ceil(FULL_SIZE_ARCSEC / PIXEL_SCALE_ARCSEC)
    estimated_plane_bytes = pixels * pixels * BYTES_PER_LAYER_PIXEL
    estimate = estimated_plane_bytes * N_LAYERS * len(visits)
    input_hashes = {
        str(results_summary): sha256(results_summary),
        str(status_path): sha256(status_path),
        str(audit_summary): sha256(audit_summary),
        str(MANUAL): sha256(MANUAL),
        **{str(RESULTS / f"frame_metadata_rows_{i:03d}.csv"): status["jobs"][f"frame_metadata_rows_{i:03d}"]["sha256"] for i in range(16)},
    }
    plan = {
        "schema": "astrafek-hsc-pdr3-warp-pilot-plan-v1",
        "status": "prepared_not_submitted",
        "execution_requires_separate_explicit_user_approval": True,
        "query_submitted": False,
        "network_used": False,
        "credentials_read_or_saved": False,
        "science_measurement": False,
        "purpose": "inspect PDR3 single-visit warp image/mask/variance layout and calibrate local variance-to-coadd mapping before considering any larger request",
        "request": {
            "endpoint": "https://hsc-release.mtk.nao.ac.jp/das_cutout/pdr3/cgi-bin/cutout",
            "rerun": "pdr3_wide",
            "type": "warp",
            "filter": "HSC-I",
            "center_id": PILOT_CENTER_ID,
            "ra_deg": expected_coords[0],
            "dec_deg": expected_coords[1],
            "semi_width": "1arcmin",
            "semi_height": "1arcmin",
            "full_tangent_plane_size_arcsec": [FULL_SIZE_ARCSEC, FULL_SIZE_ARCSEC],
            "layers": {"image": True, "mask": True, "variance": True},
            "tract": 9703,
            "patch": 306,
            "coadd_band": "HSC-I",
            "source_filter_codes": ["i2"],
            "frame_count": len(frames),
            "distinct_visits": sorted(visits),
            "distinct_ccds_in_metadata": len(ccds),
        },
        "estimated_payload": {
            "assumption": "full 2x2 arcmin square for each returned visit warp, 0.168 arcsec/pixel, 32-bit image/mask/variance planes; FITS/tar headers and edge overlap are not included",
            "pixels_per_side": pixels,
            "returned_visit_files_expected_from_metadata": len(visits),
            "estimated_uncompressed_plane_bytes": estimate,
            "estimated_uncompressed_plane_mib": estimate / (1024 ** 2),
            "maximum_download_bytes": 70 * 1024 * 1024,
            "maximum_download_mib": 70.0,
            "estimate_is_not_a_server_size_guarantee": True,
        },
        "after_receipt_checks": [
            "archive byte size and member listing",
            "exact visit IDs against the nine metadata visits",
            "FITS HDU 1/2/3 image/mask/variance availability and aligned shapes",
            "tract/patch/filter/WCS headers and pixel-scale coverage",
            "variance positivity, mask occupancy, and image-to-variance alignment",
        ],
        "input_sha256": input_hashes,
    }
    output.mkdir(parents=True, exist_ok=False)
    plan_path = output / "plan.json"
    plan_path.write_text(json.dumps(plan, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    plan_hash = sha256(plan_path)
    # Persist the plan hash in a sidecar to avoid a self-referential hash.
    (output / "plan.sha256").write_text(plan_hash + "  plan.json\n", encoding="ascii")
    print(json.dumps({"status": plan["status"], "query_submitted": False,
                      "visits": len(visits), "metadata_ccds": len(ccds),
                      "estimated_plane_mib": plan["estimated_payload"]["estimated_uncompressed_plane_mib"],
                      "plan": str(plan_path), "sha256": plan_hash}, indent=2))
    return plan


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    run(args.output)


if __name__ == "__main__":
    main()
