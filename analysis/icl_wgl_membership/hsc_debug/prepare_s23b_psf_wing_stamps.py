#!/usr/bin/env python3
"""Prepare wide, source-masked HSC stamps for Gaia-supported S23B candidates.

This local CPU workflow uses only existing S23B target or matched-control
images and the separately approved Gaia DR3 match for that population. It
prepares diagnostic PSF-wing stamps; it does not approve a PSF or make a
science measurement.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import platform
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy import units as u
from astropy.wcs import WCS

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[2]))

from download_pilot import BAD_MASK_PLANE_NAMES, mask_plane_bitmask, validate_fits
from diagnose_empirical_psf import objects_from_rows
import prepare_detection_psf as prep
from ogfkit.wlcluster import ARCSEC, Cosmology


DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
PILOT = DATA / "s23b_detection_psf_pilot_2026-11-07_01"
GAIA_AUDIT = DATA / "s23b_psf_gaia_audit_2026-11-07_01"
DEFAULT_OUTPUT = DATA / "s23b_psf_wing_stamps_2026-11-07_01"
CONTROL_PILOT = DATA / "s23b_psf_control_candidates_2026-11-07_01"
CONTROL_GAIA_AUDIT = DATA / "s23b_psf_control_gaia_audit_2026-11-07_01"
CONTROL_GAIA_RUN = DATA / "s23b_psf_control_gaia_match_2026-11-07_01"
CONTROL_GAIA_PLAN = DATA / "s23b_psf_control_gaia_plan_2026-11-07_03"
CONTROL_DEFAULT_OUTPUT = DATA / "s23b_psf_control_wing_stamps_2026-11-07_01"
EXPANSION = DATA / "s23b_psf_expansion_2026-11-07_01"
STAMP_SIDE = 401


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, payload: dict) -> None:
    with path.open("x", encoding="utf-8") as stream:
        json.dump(payload, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")


def build_stamps(output: Path, *, stamp_side: int = STAMP_SIDE,
                 population: str = "target", cohort: str = "pilot") -> dict:
    if not platform.node().lower().startswith("lageunha"):
        raise RuntimeError("Run on lageunha CPU")
    if population not in ("target", "control"):
        raise ValueError("population must be 'target' or 'control'")
    if cohort not in ("pilot", "expansion"):
        raise ValueError("cohort must be 'pilot' or 'expansion'")
    if stamp_side < STAMP_SIDE or stamp_side % 2 != 1:
        raise ValueError("stamp_side must be odd and at least 401 pixels")
    half = stamp_side // 2
    if output.exists():
        raise FileExistsError(f"refusing to overwrite stamp directory: {output}")

    if cohort == "expansion":
        pilot_dir = EXPANSION / f"{population}_preparation"
        pilot_plan_path = None
        pilot_manifest_path = EXPANSION / "plan" / f"{population}_manifest.json"
        selection_path = EXPANSION / "plan/selection.tsv"
        evidence_path = EXPANSION / "gaia_audit_2026-11-07_01/candidate_evidence.tsv"
        gaia_audit_path = EXPANSION / "gaia_audit_2026-11-07_01/summary.json"
        gaia_run_dir = EXPANSION / "gaia_match_2026-11-07_01"
        requests_path = EXPANSION / "gaia_plan/requests.json"
        if population == "target":
            expected_preparation = ("partial", 200, 199)
            expected_candidates, expected_strict = 2069, 225
        else:
            expected_preparation = ("complete", 200, 200)
            expected_candidates, expected_strict = 2230, 223
    else:
        if population == "target":
            pilot_dir = PILOT
            pilot_plan_path = pilot_dir / "plan.json"
            pilot_manifest_path = pilot_dir / "manifest.json"
            selection_path = None
            evidence_path = GAIA_AUDIT / "candidate_evidence.tsv"
            gaia_audit_path = GAIA_AUDIT / "summary.json"
            gaia_run_dir = DATA / "s23b_psf_gaia_match_2026-11-07_01"
            requests_path = DATA / "s23b_psf_gaia_plan_2026-11-07_01/requests.json"
            expected_preparation = ("partial", 50, 49)
            expected_candidates, expected_strict = 492, 48
        else:
            pilot_dir = CONTROL_PILOT
            pilot_plan_path = None
            pilot_manifest_path = pilot_dir / "pilot_manifest.json"
            selection_path = pilot_dir / "selection.tsv"
            evidence_path = CONTROL_GAIA_AUDIT / "control_candidate_evidence.tsv"
            gaia_audit_path = CONTROL_GAIA_AUDIT / "summary.json"
            gaia_run_dir = CONTROL_GAIA_RUN
            requests_path = CONTROL_GAIA_PLAN / "requests.json"
            expected_preparation = ("complete", 49, 49)
            expected_candidates, expected_strict = 458, 38
    preparation_summary_path = (pilot_dir / "summary.json" if cohort == "expansion"
                                else pilot_dir / "preparation/summary.json")
    preparation_root = pilot_dir if cohort == "expansion" else pilot_dir / "preparation"
    gaia_run_path = gaia_run_dir / "status.json"
    gaia_response_path = gaia_run_dir / "gaia_matches.csv"
    gaia_plan_path = requests_path.parent / "plan.json"
    input_paths = [pilot_manifest_path, preparation_summary_path, evidence_path,
                   gaia_audit_path, gaia_run_path, gaia_response_path,
                   gaia_plan_path, requests_path]
    if pilot_plan_path is not None:
        input_paths.append(pilot_plan_path)
    if selection_path is not None:
        input_paths.append(selection_path)
    for path in input_paths:
        if not path.is_file():
            raise FileNotFoundError(path)

    pilot_manifest = json.loads(pilot_manifest_path.read_text(encoding="utf-8"))
    preparation = json.loads(preparation_summary_path.read_text(encoding="utf-8"))
    gaia_audit = json.loads(gaia_audit_path.read_text(encoding="utf-8"))
    gaia_run = json.loads(gaia_run_path.read_text(encoding="utf-8"))
    requests = json.loads(requests_path.read_text(encoding="utf-8"))["candidates"]
    expected_query_candidates = 4299 if cohort == "expansion" else expected_candidates

    if (preparation.get("status") != expected_preparation[0]
            or preparation.get("n_requested") != expected_preparation[1]
            or preparation.get("n_prepared") != expected_preparation[2]
            or preparation.get("science_measurement") is not False):
        raise ValueError(f"unexpected S23B {population} morphology candidate preparation state")
    if (preparation.get("manifest_sha256") != sha256(pilot_manifest_path)
            or gaia_run.get("status") != "complete"
            or gaia_run.get("n_candidates") != expected_query_candidates
            or gaia_audit.get("n_candidate_positions") != expected_query_candidates
            or gaia_audit.get("psf_approved") is not False
            or gaia_audit.get("science_measurement") is not False):
        raise ValueError("Gaia query/audit is incomplete or exceeds its approved scope")
    if (sha256(gaia_plan_path) != gaia_run.get("plan_sha256")
            or sha256(requests_path) != gaia_run.get("request_sha256")
            or sha256(gaia_response_path) != gaia_run.get("response_sha256")):
        raise ValueError("Gaia plan/request/response checksum differs from completed run status")
    for source_path, expected_hash in gaia_audit.get("inputs_sha256", {}).items():
        path = Path(source_path)
        if not path.is_file() or sha256(path) != expected_hash:
            raise ValueError(f"Gaia audit input checksum mismatch: {path}")

    candidate_rows = read_tsv(evidence_path)
    if cohort == "expansion":
        candidate_rows = [row for row in candidate_rows if row.get("population") == population]
        requests_for_population = [row for row in requests if row.get("population") == population]
    else:
        requests_for_population = requests
    strict_rows = [row for row in candidate_rows if row.get("strict_candidate") == "True"]
    if (len(requests_for_population) != expected_candidates or len(candidate_rows) != expected_candidates
            or gaia_audit.get("output_sha256", {}).get(evidence_path.name) != sha256(evidence_path)):
        raise ValueError("Gaia evidence, request, or audit-output checksum differs")
    if len(strict_rows) != expected_strict:
        raise ValueError(f"expected {expected_strict} conservative Gaia candidates, found {len(strict_rows)}")
    if len({row["best_source_id"] for row in strict_rows}) != len(strict_rows):
        raise ValueError("a Gaia source is associated with more than one strict HSC candidate")

    request_by_index = {int(row["local_index"]): row for row in requests}
    frame_by_id = {row["id"]: row for row in pilot_manifest["clusters"]}
    if cohort == "expansion":
        selection_rows = read_tsv(selection_path)
        selection_key = "target_cluster_id" if population == "target" else "control_id"
        selection_hash_key = f"{population}_image_sha256"
        selection_by_id = {
            row[selection_key]: {"id": row[selection_key],
                                 "image_sha256": row[selection_hash_key]}
            for row in selection_rows
        }
        expected_selection_count = 200
    else:
        expected_selection_count = 49
        if pilot_plan_path is not None:
            pilot_plan = json.loads(pilot_plan_path.read_text(encoding="utf-8"))
            selection_by_id = {row["id"]: row for row in pilot_plan["selection"]}
        else:
            selection_by_id = {
                row["control_id"]: {"id": row["control_id"], "image_sha256": row["image_sha256"]}
                for row in read_tsv(selection_path)
            }
    if len(selection_by_id) != expected_selection_count:
        raise ValueError(f"{cohort} {population} selection does not contain {expected_selection_count} unique frames")
    preparation_by_id = {row["id"]: row for row in preparation["records"]}
    if len(request_by_index) != len(requests):
        raise ValueError("approved Gaia candidate indices are duplicated")
    grouped: dict[str, list[dict]] = defaultdict(list)
    for evidence in strict_rows:
        index = int(evidence["local_index"])
        request = request_by_index.get(index)
        if request is None:
            raise ValueError(f"candidate index is absent from approved request: {index}")
        evidence_cluster_id = evidence.get("cluster_id", evidence.get("control_id"))
        if (request["cluster_id"] != evidence_cluster_id
                or int(request["number"]) != int(evidence["number"])):
            raise ValueError(f"Gaia match no longer maps to the same HSC detection: {index}")
        if request["classification"] != "provisional_compact_candidate":
            raise ValueError(f"unexpected input candidate classification at {index}")
        grouped[evidence_cluster_id].append({"evidence": evidence, "request": request})

    input_hashes = {str(path.resolve()): sha256(path) for path in input_paths}
    frame_hashes = {}
    candidate_table_hashes = {}
    output.mkdir(parents=True, exist_ok=False)
    stamps_dir = output / "stamps"
    stamps_dir.mkdir()
    ledger_rows = []
    status_rows = []
    cosmology = Cosmology(H0=pilot_manifest["h0"], Om=pilot_manifest["omega_m"])

    for cluster_id in sorted(grouped, key=lambda value: int(value.rsplit("_", 1)[1])):
        frame = frame_by_id[cluster_id]
        selection = selection_by_id[cluster_id]
        if cluster_id not in preparation_by_id:
            raise ValueError(f"Gaia candidate frame is absent from local preparation: {cluster_id}")
        image_path = Path(frame["image"]).resolve()
        expected_image_sha = selection["image_sha256"]
        actual_image_sha = sha256(image_path)
        if actual_image_sha != expected_image_sha:
            raise ValueError(f"HSC image checksum mismatch for {cluster_id}")
        frame_hashes[str(image_path)] = actual_image_sha

        frame_dir = preparation_root / f"cluster_{cluster_id}"
        detections_path = frame_dir / "detections.tsv"
        masks_path = frame_dir / "masks.fits.gz"
        candidate_path = frame_dir / "psf_candidates.tsv"
        for path in (detections_path, masks_path, candidate_path):
            if not path.is_file():
                raise FileNotFoundError(path)
        candidate_table_hashes[str(detections_path)] = sha256(detections_path)
        candidate_table_hashes[str(masks_path)] = sha256(masks_path)
        candidate_table_hashes[str(candidate_path)] = sha256(candidate_path)

        detections = read_tsv(detections_path)
        candidate_table = read_tsv(candidate_path)
        candidate_by_number = {int(row["NUMBER"]): row for row in candidate_table}
        objects, kron, kron_flags = objects_from_rows(detections, 0, 0)
        member_sky = np.asarray(frame.get("members_radec_deg", []), dtype=float).reshape(-1, 2)

        with fits.open(image_path, memmap=True, do_not_scale_image_data=True, lazy_load_hdus=True) as image_hdus, \
                fits.open(masks_path, memmap=False) as mask_hdus:
            info = validate_fits(image_path, frame)
            filter_name = info["coadd_header_filter"]
            science_hdu = int(info["science_hdu"])
            mask_hdu = int(info["mask_hdu"])
            variance_hdu = int(info["variance_hdu"])
            image_plane = image_hdus[science_hdu]
            variance_plane = image_hdus[variance_hdu]
            bits_plane = image_hdus[mask_hdu]
            shape = image_plane.shape
            if any(plane.shape != shape for plane in (variance_plane, bits_plane)):
                raise ValueError(f"image, variance, and HSC mask shapes differ: {cluster_id}")
            seg = np.asarray(mask_hdus["SEGMENTATION"].data)
            saved_base = np.asarray(mask_hdus["BASE_MASK"].data, dtype=bool)
            if seg.shape != shape or saved_base.shape != shape:
                raise ValueError(f"saved segmentation/base mask shape differs: {cluster_id}")
            image = np.asarray(image_plane.data, dtype=np.float32)
            variance = np.asarray(variance_plane.data, dtype=np.float32)
            bits = np.asarray(bits_plane.data, dtype=np.int32)
            bad_bits = mask_plane_bitmask(info["mask_planes"], BAD_MASK_PLANE_NAMES)
            bright_bit = mask_plane_bitmask(info["mask_planes"], ("MP_BRIGHT_OBJECT",))
            original_base = ((bits & bad_bits) != 0) | ((bits & bright_bit) != 0)
            original_base |= ~np.isfinite(image) | ~np.isfinite(variance) | (variance <= 0)
            if not np.array_equal(original_base, saved_base):
                raise ValueError(f"saved BASE_MASK differs from source HSC quality planes: {cluster_id}")

            header = image_plane.header.copy()
            if image_hdus[0].header.get("INHERIT") is True:
                header.extend(image_hdus[0].header, unique=True)
            wcs = WCS(header).celestial
            if not wcs.has_celestial:
                raise ValueError(f"target image has no celestial WCS: {cluster_id}")
            if len(member_sky):
                member_xy = np.column_stack(wcs.all_world2pix(member_sky[:, 0], member_sky[:, 1], 0))
                inside = (
                    np.isfinite(member_xy).all(axis=1)
                    & (member_xy[:, 0] >= 0) & (member_xy[:, 0] < shape[1])
                    & (member_xy[:, 1] >= 0) & (member_xy[:, 1] < shape[0])
                )
                member_xy = member_xy[inside]
            else:
                member_xy = np.empty((0, 2), dtype=float)

            for pair in grouped[cluster_id]:
                evidence = pair["evidence"]
                request = pair["request"]
                local_index = int(evidence["local_index"])
                number = int(evidence["number"])
                detection = detections[number - 1]
                if int(detection["NUMBER"]) != number or number not in candidate_by_number:
                    raise ValueError(f"candidate detection is missing: {cluster_id}/{number}")
                candidate = candidate_by_number[number]
                if candidate.get("CLASSIFICATION") != "provisional_compact_candidate":
                    raise ValueError(f"unexpected morphology candidate class: {cluster_id}/{number}")
                ra = float(detection["RA"])
                dec = float(detection["DEC"])
                if (abs(ra - float(request["ra_deg"])) > 1e-9
                        or abs(dec - float(request["dec_deg"])) > 1e-9
                        or abs(ra - float(evidence["ra_deg"])) > 1e-9
                        or abs(dec - float(evidence["dec_deg"])) > 1e-9):
                    raise ValueError(f"candidate coordinates disagree across pinned products: {cluster_id}/{number}")

                x = float(detection["X_IMAGE"]) - 1.0
                y = float(detection["Y_IMAGE"]) - 1.0
                ix, iy = int(round(x)), int(round(y))
                x0, y0 = ix - half, iy - half
                if (x0 < 0 or y0 < 0 or x0 + stamp_side > shape[1]
                        or y0 + stamp_side > shape[0]):
                    status_rows.append({
                        "local_index": local_index, "cluster_id": cluster_id,
                        "number": number, "source_id": evidence["best_source_id"],
                        "status": "blocked_stamp_edge", "stamp": "",
                    })
                    continue
                if str(evidence["phot_variable_flag"]).upper() == "VARIABLE":
                    status_rows.append({
                        "local_index": local_index, "cluster_id": cluster_id,
                        "number": number, "source_id": evidence["best_source_id"],
                        "status": "blocked_known_variable", "stamp": "",
                    })
                    continue

                masks, fallback_count = prep.source_masks(
                    objects, seg, kron, kron_flags, saved_base,
                    number - 1, member_xy, float(info["pixel_scale_arcsec"]),
                )
                stamp_slice = np.s_[y0:y0 + stamp_side, x0:x0 + stamp_side]
                science = np.asarray(image[stamp_slice], dtype=np.float32)
                stamp_variance = np.asarray(variance[stamp_slice], dtype=np.float32)
                base_mask = np.asarray(saved_base[stamp_slice], dtype=np.uint8)
                segmentation = np.asarray(seg[stamp_slice], dtype=np.int32)
                mask_25 = np.asarray(masks[0][stamp_slice], dtype=np.uint8)
                mask_35 = np.asarray(masks[1][stamp_slice], dtype=np.uint8)
                if any(array.shape != (stamp_side, stamp_side) for array in (
                    science, stamp_variance, base_mask, segmentation, mask_25, mask_35,
                )):
                    raise ValueError(f"wide PSF stamp has wrong dimensions: {cluster_id}/{number}")
                if np.any(base_mask.astype(bool) & ~mask_25.astype(bool)):
                    raise ValueError("Kron 2.5 mask does not preserve the base mask")
                if np.any(mask_25.astype(bool) & ~mask_35.astype(bool)):
                    raise ValueError("Kron masks are not nested")

                center_x = x - x0
                center_y = y - y0
                pixel_scale = float(info["pixel_scale_arcsec"])
                cluster = frame
                sky_sep = SkyCoord(ra * u.deg, dec * u.deg).separation(
                    SkyCoord(float(cluster["ra_deg"]) * u.deg, float(cluster["dec_deg"]) * u.deg)
                ).arcsec
                radius_kpc = float(sky_sep) * cosmology.angular_diameter(float(cluster["z"])) * 1000 * ARCSEC

                filename = f"star_{local_index:04d}_{evidence['best_source_id']}.fits"
                stamp_path = stamps_dir / filename
                stamp_header = fits.Header()
                stamp_header["SCIVAL"] = False
                stamp_header["PSFOK"] = False
                stamp_header["PSFFIT"] = False
                stamp_header["BGSUB"] = False
                stamp_header["X0"] = x0
                stamp_header["Y0"] = y0
                stamp_header["XCENTER"] = center_x
                stamp_header["YCENTER"] = center_y
                stamp_header["PIXSCALE"] = pixel_scale
                stamp_header["CANDID"] = local_index
                stamp_header["SOURCEID"] = evidence["best_source_id"]
                stamp_header["CLUSTER"] = cluster_id
                stamp_header["POPUL"] = population
                stamp_header["FILTER"] = filter_name
                stamp_header["ZCL"] = float(cluster["z"])
                stamp_header["RCL_KPC"] = radius_kpc
                stamp_header["GAIAEPO"] = 2016.0
                stamp_header["GSEP"] = float(evidence["best_separation_arcsec"])
                stamp_header["RUWE"] = float(evidence["ruwe"])
                if evidence.get("pm_shift_arcsec_at_5yr") not in (None, ""):
                    stamp_header["PM5YR"] = float(evidence["pm_shift_arcsec_at_5yr"])
                stamp_header["COMMENT"] = "Gaia diagnostic candidate; local sky retained; not PSF-approved."
                fits.HDUList([
                    fits.PrimaryHDU(header=stamp_header),
                    fits.ImageHDU(science, name="SCIENCE"),
                    fits.ImageHDU(stamp_variance, name="VARIANCE"),
                    fits.ImageHDU(base_mask, name="BASE_MASK"),
                    fits.ImageHDU(segmentation, name="SEGMENTATION"),
                    fits.ImageHDU(mask_25, name="MASK_K2P5"),
                    fits.ImageHDU(mask_35, name="MASK_K3P5"),
                ]).writeto(stamp_path, overwrite=False)
                stamp_sha = sha256(stamp_path)
                ledger_rows.append({
                    "candidate_id": str(local_index),
                    "cluster_id": cluster_id,
                    "number": number,
                    "source_id": evidence["best_source_id"],
                    "band": filter_name,
                    "cluster_radius_kpc": radius_kpc,
                    "gaia_separation_arcsec": float(evidence["best_separation_arcsec"]),
                    "gaia_star_probability": float(evidence["star_probability"]),
                    "gaia_pm_significance": float(evidence["pm_significance"]),
                    "gaia_ruwe": float(evidence["ruwe"]),
                    "phot_g_mean_mag": float(evidence["phot_g_mean_mag"]),
                    "phot_variable_flag": evidence["phot_variable_flag"],
                    "stamp": str(stamp_path.resolve()),
                    "stamp_sha256": stamp_sha,
                    "image_sha256": actual_image_sha,
                    "mask_fallback_members": fallback_count,
                    "x_center": center_x,
                    "y_center": center_y,
                    "pixel_scale_arcsec": pixel_scale,
                    "status": "complete_diagnostic_only",
                    "science_measurement": False,
                    "psf_approved": False,
                })
                status_rows.append({
                    "local_index": local_index, "cluster_id": cluster_id,
                    "number": number, "source_id": evidence["best_source_id"],
                    "status": "complete_diagnostic_only", "stamp": str(stamp_path.resolve()),
                })

    ledger_fields = list(ledger_rows[0]) if ledger_rows else [
        "candidate_id", "cluster_id", "number", "source_id", "band",
        "cluster_radius_kpc", "stamp", "stamp_sha256",
    ]
    write_tsv(output / "star_ledger.tsv", ledger_rows, ledger_fields)
    write_tsv(output / "star_stamp_status.tsv", status_rows,
              ["local_index", "cluster_id", "number", "source_id", "status", "stamp"])

    subsets = {}
    subset_root = output / "stack_ledgers"
    subset_root.mkdir()
    for band in ("i", "i2"):
        for name, predicate in (
            ("all", lambda row: True),
            ("outer500", lambda row: float(row["cluster_radius_kpc"]) >= 500.0),
        ):
            rows = [row for row in ledger_rows if row["band"] == band and predicate(row)]
            subdir = subset_root / f"{band}_{name}"
            subdir.mkdir()
            write_tsv(subdir / "star_ledger.tsv", rows, ledger_fields)
            subsets[f"{band}_{name}"] = {
                "n_stars": len(rows),
                "n_frames": len({row["cluster_id"] for row in rows}),
                "ledger": str((subdir / "star_ledger.tsv").resolve()),
                "ledger_sha256": sha256(subdir / "star_ledger.tsv"),
            }

    input_hashes.update(frame_hashes)
    input_hashes.update(candidate_table_hashes)
    output_hashes = {
        "star_ledger.tsv": sha256(output / "star_ledger.tsv"),
        "star_stamp_status.tsv": sha256(output / "star_stamp_status.tsv"),
    }
    output_hashes.update({str(Path(row["stamp"]).relative_to(output)): row["stamp_sha256"] for row in ledger_rows})
    output_hashes.update({
        str(Path(value["ledger"]).relative_to(output)): value["ledger_sha256"]
        for value in subsets.values()
    })
    summary = {
        "schema": "astrafek-s23b-psf-wing-stamps-v2",
        "status": "complete_diagnostic_only",
        "cohort": cohort,
        "population": population,
        "n_gaia_candidates_selected": len(strict_rows),
        "n_wide_stamps": len(ledger_rows),
        "n_blocked_candidates": len(strict_rows) - len(ledger_rows),
        "n_sample_frames": len({row["cluster_id"] for row in ledger_rows}),
        "n_target_frames": len({row["cluster_id"] for row in ledger_rows}) if population == "target" else 0,
        "n_control_frames": len({row["cluster_id"] for row in ledger_rows}) if population == "control" else 0,
        "n_filter_i_stamps": sum(row["band"] == "i" for row in ledger_rows),
        "n_filter_i2_stamps": sum(row["band"] == "i2" for row in ledger_rows),
        "n_frames_with_stamps": len({row["cluster_id"] for row in ledger_rows}),
        "subset_ledgers": subsets,
        "stamp_shape": [stamp_side, stamp_side],
        "stamps_are_background_subtracted": False,
        "stamps_are_subpixel_recentered": False,
        "source_masks": ["Kron 2.5", "Kron 3.5"],
        "science_measurement": False,
        "psf_approved": False,
        "external_queries_submitted": False,
        "inputs_sha256": input_hashes,
        "script_sha256": sha256(Path(__file__).resolve()),
        "output_sha256": output_hashes,
        "limitations": [
            f"Only {len(ledger_rows)} of the {len(strict_rows)} Gaia diagnostic candidates yielded usable edge/variability-filtered stamps in {len({row['cluster_id'] for row in ledger_rows})} {population} frames.",
            "S23B image headers in the inspected pilot do not identify the effective coadd epoch; Gaia positions are not propagated to a known HSC epoch.",
            "The mask model removes detected neighbors and defects but cannot guarantee removal of undetected sources.",
            "PSF wings, sky-annulus transfer, and coadd covariance remain uncalibrated.",
            "For controls, RCL_KPC is distance from the random sky center at the paired target redshift; it is not a cluster-centric radius.",
            "Gaia diagnostic thresholds do not establish a calibrated star identity or PSF suitability.",
        ],
    }
    write_json(output / "summary.json", summary)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--population", choices=("target", "control"), default="target")
    parser.add_argument("--cohort", choices=("pilot", "expansion"), default="pilot")
    parser.add_argument("--output", type=Path, default=None)
    parser.add_argument("--stamp-side", type=int, default=STAMP_SIDE,
                        help="odd stamp width in pixels; default 401")
    args = parser.parse_args()
    if args.cohort == "expansion":
        default_output = EXPANSION / f"{args.population}_wing_stamps_2026-11-07_01"
    else:
        default_output = DEFAULT_OUTPUT if args.population == "target" else CONTROL_DEFAULT_OUTPUT
    output = args.output.resolve() if args.output else default_output
    summary = build_stamps(output, stamp_side=args.stamp_side,
                           population=args.population, cohort=args.cohort)
    print(json.dumps({
        key: summary[key] for key in (
            "status", "cohort", "n_gaia_candidates_selected", "n_wide_stamps",
            "n_blocked_candidates", "population", "n_sample_frames", "n_filter_i_stamps",
            "n_filter_i2_stamps", "subset_ledgers", "science_measurement", "psf_approved",
        )
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
