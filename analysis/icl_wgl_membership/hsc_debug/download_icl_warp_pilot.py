#!/usr/bin/env python3
"""Download and audit one explicitly approved HSC PDR3 warp pilot.

Default invocation is plan-only. The network path requires both ``--submit``
and the exact reviewed plan SHA. It reuses the existing HSC cutout endpoint and
Basic-Auth helper from ``07_fetch_cutouts.py``. Credentials are read in memory
only and never written to output or logs.
"""
from __future__ import annotations

import argparse
import getpass
import hashlib
import importlib.util
import io
import json
import math
import os
from pathlib import Path
import re
import tarfile
import tempfile
import urllib.error
import urllib.parse
import urllib.request

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales

HERE = Path(__file__).resolve().parent
ROOT = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
DEFAULT_PLAN = ROOT / "icl_warp_pilot_plan_2026-11-06_02"
DEFAULT_OUTPUT = ROOT / "icl_warp_pilot_download_2026-11-06_01"
EXPECTED_PLAN_SHA256 = "a436df245a33bda09cd5599e564b2b9302e5d0d4c006b59fcf6fab1c2c33a65e"
MAX_DOWNLOAD_BYTES = 70 * 1024 * 1024
READ_CHUNK_BYTES = 4 * 1024 * 1024


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def json_compatible(value):
    """Convert NumPy scalar values in FITS QA summaries to JSON primitives."""
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, dict):
        return {str(key): json_compatible(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_compatible(item) for item in value]
    return value


def load_existing_cutout_helpers():
    helper_path = HERE / "07_fetch_cutouts.py"
    spec = importlib.util.spec_from_file_location("astrafek_hsc_cutout_probe", helper_path)
    if spec is None or spec.loader is None:
        raise RuntimeError("cannot load the existing HSC cutout helper")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def load_plan(plan_dir: Path) -> tuple[dict, Path]:
    plan_dir = Path(plan_dir).resolve()
    plan_path = plan_dir / "plan.json"
    hash_path = plan_dir / "plan.sha256"
    if not plan_path.is_file() or not hash_path.is_file():
        raise FileNotFoundError("warp pilot plan.json or plan.sha256 is missing")
    actual_hash = sha256(plan_path)
    recorded_hash = hash_path.read_text(encoding="ascii").split()[0]
    if actual_hash != recorded_hash or actual_hash != EXPECTED_PLAN_SHA256:
        raise ValueError("warp pilot plan differs from the exact reviewed plan")
    plan = json.loads(plan_path.read_text(encoding="utf-8"))
    if (plan.get("schema") != "astrafek-hsc-pdr3-warp-pilot-plan-v1"
            or plan.get("status") != "prepared_not_submitted"
            or plan.get("query_submitted") is not False
            or plan.get("execution_requires_separate_explicit_user_approval") is not True
            or plan.get("network_used") is not False):
        raise ValueError("plan state or approval gate changed")
    if plan.get("estimated_payload", {}).get("maximum_download_bytes") != MAX_DOWNLOAD_BYTES:
        raise ValueError("plan does not match the downloader's explicit 70 MiB ceiling")
    for path_text, expected in plan.get("input_sha256", {}).items():
        path = Path(path_text)
        if not path.is_file() or sha256(path) != expected:
            raise ValueError(f"pinned plan input changed: {path}")
    return plan, plan_path


def query_url(endpoint: str, request: dict) -> str:
    values = {
        "ra": f"{request['ra_deg']:.10f}",
        "dec": f"{request['dec_deg']:.10f}",
        "sw": request["semi_width"],
        "sh": request["semi_height"],
        "type": request["type"],
        "image": "true",
        "mask": "true",
        "variance": "true",
        "filter": request["filter"],
        "rerun": request["rerun"],
        "tract": str(request["tract"]),
    }
    return endpoint.rstrip("?") + "?" + urllib.parse.urlencode(values)


def stream_response(url: str, output_dir: Path, user: str, password: str,
                    opener_factory) -> tuple[Path, int, str]:
    request = urllib.request.Request(url, headers={"User-Agent": "Astrafek/ICL-warp-pilot"})
    opener = opener_factory(user, password)
    temp_path = None
    try:
        with opener.open(request, timeout=180) as response:
            status_code = int(response.status)
            content_type = response.headers.get("Content-Type", "")
            content_length = response.headers.get("Content-Length")
            if content_length is not None and int(content_length) > MAX_DOWNLOAD_BYTES:
                raise ValueError("server Content-Length exceeds the approved 70 MiB pilot ceiling")
            handle = tempfile.NamedTemporaryFile(
                prefix=".warp-pilot-", suffix=".partial", dir=output_dir, delete=False,
            )
            temp_path = Path(handle.name)
            total = 0
            with handle:
                while True:
                    chunk = response.read(READ_CHUNK_BYTES)
                    if not chunk:
                        break
                    total += len(chunk)
                    if total > MAX_DOWNLOAD_BYTES:
                        raise ValueError("streamed response exceeds the approved 70 MiB pilot ceiling")
                    handle.write(chunk)
        if status_code != 200:
            raise ValueError(f"cutout endpoint returned HTTP {status_code}")
        return temp_path, total, content_type
    except urllib.error.HTTPError as exc:
        if temp_path and temp_path.exists():
            temp_path.unlink()
        raise RuntimeError(f"cutout endpoint returned HTTP {exc.code}") from None
    except Exception:
        if temp_path and temp_path.exists():
            temp_path.unlink()
        raise


def inspect_tar(archive_path: Path, expected_visits: set[int], *, expected_filter: str,
                expected_tract: int, expected_rerun: str, center_ra: float,
                center_dec: float, expected_pixel_scale_arcsec: float = 0.168) -> dict:
    visit_pattern = re.compile(r"(?:^|/)warp-(\d+)\.fits$")
    directory_pattern = re.compile(
        r"(?:^|/)warps-([^/]+)-(\d+)-([^/]+)-(\d{6})-(\d{6})/warp-(\d+)\.fits$"
    )
    members = []
    visit_to_member = {}
    with tarfile.open(archive_path, mode="r:*") as archive:
        for member in archive.getmembers():
            if not member.isfile():
                continue
            match = visit_pattern.search(member.name)
            if not match:
                continue
            visit = int(match.group(1))
            directory_match = directory_pattern.search(member.name)
            if not directory_match:
                raise ValueError(f"unrecognized warp archive member path: {member.name}")
            member_filter, member_tract, member_rerun, _date, _clock, member_visit = directory_match.groups()
            if (member_filter != expected_filter or int(member_tract) != expected_tract
                    or member_rerun != expected_rerun or int(member_visit) != visit):
                raise ValueError(f"warp archive member does not match the approved filter/tract/rerun: {member.name}")
            if visit in visit_to_member:
                raise ValueError(f"duplicate warp visit in archive: {visit}")
            file_obj = archive.extractfile(member)
            if file_obj is None:
                raise ValueError(f"cannot read FITS member for visit {visit}")
            data = file_obj.read()
            with fits.open(io.BytesIO(data), memmap=False) as hdul:
                if len(hdul) < 4:
                    raise ValueError(f"visit {visit} FITS lacks fixed image/mask/variance HDUs")
                planes = []
                headers = []
                for hdu_index in (1, 2, 3):
                    plane = hdul[hdu_index].data
                    if plane is None or getattr(plane, "ndim", 0) != 2:
                        raise ValueError(f"visit {visit} HDU {hdu_index} is empty or not 2-D")
                    planes.append(plane)
                    headers.append(hdul[hdu_index].header)
                shapes = [list(plane.shape) for plane in planes]
                if len({tuple(shape) for shape in shapes}) != 1:
                    raise ValueError(f"visit {visit} image/mask/variance plane shapes differ")
                wcs_keys = (
                    "CTYPE1", "CTYPE2", "CRPIX1", "CRPIX2", "CRVAL1", "CRVAL2",
                    "CD1_1", "CD1_2", "CD2_1", "CD2_2", "CDELT1", "CDELT2",
                    "PC1_1", "PC1_2", "PC2_1", "PC2_2", "CROTA2",
                )
                wcs_signatures = [
                    {key: header.get(key) for key in wcs_keys if header.get(key) is not None}
                    for header in headers
                ]
                wcs_aligned = bool(wcs_signatures[0]) and all(
                    signature == wcs_signatures[0] for signature in wcs_signatures[1:]
                )
                try:
                    image_wcs = WCS(headers[0]).celestial
                    if image_wcs.pixel_n_dim != 2 or image_wcs.world_n_dim != 2:
                        raise ValueError("image WCS is not two-dimensional celestial WCS")
                    center_x, center_y = image_wcs.world_to_pixel_values(center_ra, center_dec)
                    pixel_scales = np.asarray(proj_plane_pixel_scales(image_wcs), dtype=float) * 3600.0
                    pixel_scale_arcsec = [float(abs(v)) for v in pixel_scales]
                    center_inside_image = (
                        math.isfinite(float(center_x)) and math.isfinite(float(center_y))
                        and 0 <= center_x < shapes[0][1] and 0 <= center_y < shapes[0][0]
                    )
                    pixel_scale_matches = bool(
                        np.all(np.isfinite(pixel_scales))
                        and np.all(np.abs(pixel_scales / expected_pixel_scale_arcsec - 1.0) <= 0.05)
                    )
                except Exception as exc:
                    center_x = center_y = None
                    pixel_scale_arcsec = []
                    center_inside_image = False
                    pixel_scale_matches = False
                    wcs_error = f"{type(exc).__name__}: {exc}"
                else:
                    wcs_error = None
                variance = planes[2]
                finite_variance = variance[np.isfinite(variance)]
                mask_plane = planes[1]
                variance_check = {
                    "n_pixels": int(variance.size),
                    "n_finite": int(finite_variance.size),
                    "n_positive": int((finite_variance > 0).sum()),
                    "finite_fraction": float(np.isfinite(variance).mean()),
                    "positive_fraction_of_finite": (float((finite_variance > 0).mean())
                                                     if finite_variance.size else 0.0),
                }
                image_check = {
                    "n_pixels": int(planes[0].size),
                    "n_finite": int(np.isfinite(planes[0]).sum()),
                    "finite_fraction": float(np.isfinite(planes[0]).mean()),
                }
                mask_check = {
                    "n_pixels": int(mask_plane.size),
                    "n_nonzero": int(np.count_nonzero(mask_plane)),
                    "nonzero_fraction": float(np.count_nonzero(mask_plane) / mask_plane.size),
                }
                header = hdul[1].header
                header_filter = header.get("FILTER")
                header_tract = header.get("TRACT")
                header_patch = header.get("PATCH")
                header_visit = header.get("VISIT")
                if header_filter is not None and header_filter != expected_filter:
                    raise ValueError(f"visit {visit} FITS FILTER header disagrees with plan")
                if header_tract is not None and int(header_tract) != expected_tract:
                    raise ValueError(f"visit {visit} FITS TRACT header disagrees with plan")
                if header_visit is not None and int(header_visit) != visit:
                    raise ValueError(f"visit {visit} FITS VISIT header disagrees with file name")
                members.append({
                    "member": member.name,
                    "visit": visit,
                    "member_bytes": int(member.size),
                    "hdu_count": len(hdul),
                    "image_mask_variance_shapes": shapes,
                    "image_mask_variance_wcs_aligned": wcs_aligned,
                    "image_finite": image_check,
                    "mask_occupancy": mask_check,
                    "variance_finite_positive": variance_check,
                    "wcs_center_pixel_xy": [None if center_x is None else float(center_x),
                                             None if center_y is None else float(center_y)],
                    "wcs_center_inside_image": center_inside_image,
                    "wcs_pixel_scale_arcsec": pixel_scale_arcsec,
                    "wcs_pixel_scale_matches_0p168_arcsec_within_5pct": pixel_scale_matches,
                    "wcs_error": wcs_error,
                    "header_filter": header_filter,
                    "header_tract": header_tract,
                    "header_patch": header_patch,
                    "header_visit": header_visit,
                    "crval1": header.get("CRVAL1"),
                    "crval2": header.get("CRVAL2"),
                    "cdelt1": header.get("CDELT1"),
                    "cdelt2": header.get("CDELT2"),
                })
            visit_to_member[visit] = member.name
    actual_visits = set(visit_to_member)
    members_ok = bool(members) and all(
        member["image_mask_variance_wcs_aligned"]
        and member["wcs_center_inside_image"]
        and member["wcs_pixel_scale_matches_0p168_arcsec_within_5pct"]
        and member["image_finite"]["finite_fraction"] > 0
        and member["variance_finite_positive"]["n_positive"] > 0
        for member in members
    )
    return {
        "member_count": len(members),
        "expected_visits": sorted(expected_visits),
        "actual_visits": sorted(actual_visits),
        "visit_set_exact_match": actual_visits == expected_visits,
        "all_member_layout_wcs_and_variance_checks_pass": members_ok,
        "members": members,
    }


def run(plan_dir: Path, output_dir: Path, *, submit: bool,
        confirmed_plan_sha: str | None) -> dict:
    plan, plan_path = load_plan(plan_dir)
    request = plan["request"]
    plan_sha = sha256(plan_path)
    if not submit:
        print(json.dumps({"status": plan["status"], "query_submitted": False,
                          "plan_sha256": plan_sha,
                          "expected_visits": request["distinct_visits"],
                          "estimated_payload_mib": plan["estimated_payload"]["estimated_uncompressed_plane_mib"],
                          "next_action": "requires separate explicit user approval"}, indent=2))
        return {"status": "plan_only"}
    if confirmed_plan_sha != plan_sha:
        raise PermissionError("submission requires the exact separately approved plan SHA-256")
    output_dir = Path(output_dir).resolve()
    if output_dir.exists():
        raise FileExistsError(f"refusing to overwrite warp pilot output: {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=False)
    user = os.environ.get("HSC_USER") or input("HSC account name: ").strip()
    password = os.environ.get("HSC_PASSWORD") or getpass.getpass("HSC password (hidden): ")
    if not user or not password:
        raise ValueError("HSC_USER and HSC_PASSWORD are required")

    helpers = load_existing_cutout_helpers()
    url = query_url(helpers.BASE, request)
    archive_partial, nbytes, content_type = stream_response(
        url, output_dir, user, password, helpers.opener_for,
    )
    archive_path = output_dir / "warp_pilot.tar"
    archive_partial.replace(archive_path)
    audit = inspect_tar(
        archive_path, set(request["distinct_visits"]),
        expected_filter=request["filter"], expected_tract=int(request["tract"]),
        expected_rerun=request["rerun"], center_ra=float(request["ra_deg"]),
        center_dec=float(request["dec_deg"]),
    )
    if not audit["visit_set_exact_match"]:
        status_value = "complete_scope_mismatch"
    elif not audit["all_member_layout_wcs_and_variance_checks_pass"]:
        status_value = "complete_with_qa_warnings"
    else:
        status_value = "complete"
    summary = {
        "schema": "astrafek-hsc-warp-pilot-download-v1",
        "status": status_value,
        "plan_sha256": plan_sha,
        "query_submitted": True,
        "network_used": True,
        "credentials_saved": False,
        "downloaded_bytes": nbytes,
        "content_type": content_type,
        "archive": str(archive_path),
        "archive_sha256": sha256(archive_path),
        "archive_audit": audit,
        "science_measurement": False,
    }
    summary_path = output_dir / "summary.json"
    summary_path.write_text(json.dumps(json_compatible(summary), ensure_ascii=False,
                                       indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({"status": summary["status"], "downloaded_bytes": nbytes,
                      "archive_sha256": summary["archive_sha256"],
                      "visit_set_exact_match": audit["visit_set_exact_match"],
                      "member_qa_pass": audit["all_member_layout_wcs_and_variance_checks_pass"],
                      "science_measurement": False}, indent=2))
    return summary


def audit_existing(plan_dir: Path, output_dir: Path) -> dict:
    """Finish FITS/WCS QA for an already downloaded archive without networking."""
    plan, plan_path = load_plan(plan_dir)
    output_dir = Path(output_dir).resolve()
    archive_path = output_dir / "warp_pilot.tar"
    summary_path = output_dir / "summary.json"
    if not output_dir.is_dir() or not archive_path.is_file():
        raise FileNotFoundError("existing warp-pilot output directory or archive is missing")
    if summary_path.exists():
        raise FileExistsError(f"refusing to overwrite existing audit summary: {summary_path}")
    nbytes = archive_path.stat().st_size
    if nbytes > MAX_DOWNLOAD_BYTES:
        raise ValueError("saved archive exceeds the approved 70 MiB pilot ceiling")
    request = plan["request"]
    audit = inspect_tar(
        archive_path, set(request["distinct_visits"]),
        expected_filter=request["filter"], expected_tract=int(request["tract"]),
        expected_rerun=request["rerun"], center_ra=float(request["ra_deg"]),
        center_dec=float(request["dec_deg"]),
    )
    if not audit["visit_set_exact_match"]:
        status_value = "complete_scope_mismatch"
    elif not audit["all_member_layout_wcs_and_variance_checks_pass"]:
        status_value = "complete_with_qa_warnings"
    else:
        status_value = "complete"
    summary = {
        "schema": "astrafek-hsc-warp-pilot-download-v1",
        "status": status_value,
        "plan_sha256": sha256(plan_path),
        "query_submitted": True,
        "network_used": True,
        "credentials_saved": False,
        "downloaded_bytes": nbytes,
        "content_type": None,
        "content_type_note": "not preserved by the initial download attempt; no second request was made",
        "archive": str(archive_path),
        "archive_sha256": sha256(archive_path),
        "archive_audit": audit,
        "audit_execution_mode": "offline_existing_archive",
        "science_measurement": False,
    }
    summary_path.write_text(json.dumps(json_compatible(summary), ensure_ascii=False,
                                       indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({"status": status_value, "downloaded_bytes": nbytes,
                      "archive_sha256": summary["archive_sha256"],
                      "visit_set_exact_match": audit["visit_set_exact_match"],
                      "member_qa_pass": audit["all_member_layout_wcs_and_variance_checks_pass"],
                      "summary": str(summary_path), "science_measurement": False}, indent=2))
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan-dir", type=Path, default=DEFAULT_PLAN)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--submit", action="store_true", help="submit only after separately approved exact plan SHA")
    parser.add_argument("--audit-existing", action="store_true",
                        help="audit a previously downloaded archive without making a network request")
    parser.add_argument("--confirmed-plan-sha", help="exact SHA-256 explicitly approved by the user")
    args = parser.parse_args()
    if args.submit and args.audit_existing:
        parser.error("choose either --submit or --audit-existing")
    if args.audit_existing:
        if args.confirmed_plan_sha is not None:
            parser.error("--confirmed-plan-sha applies only to --submit")
        audit_existing(args.plan_dir, args.output_dir)
    else:
        run(args.plan_dir, args.output_dir, submit=args.submit,
            confirmed_plan_sha=args.confirmed_plan_sha)


if __name__ == "__main__":
    main()
