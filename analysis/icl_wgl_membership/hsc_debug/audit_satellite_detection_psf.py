#!/usr/bin/env python3
"""Independently audit the 800-frame detection/segmentation preparation."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import platform

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS

from download_pilot import (BAD_BIT_POSITIONS, BRIGHT_OBJECT_BIT, validate_fits)


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def line_count(path):
    with Path(path).open("rb") as stream:
        return sum(1 for _ in stream)


def image_path_for(frame, manifest_path):
    image_value = frame.get("image")
    if not image_value:
        raise ValueError(f"Missing image path for {frame.get('id')}")
    path = Path(image_value)
    if path.is_absolute():
        return path
    first = manifest_path.parent / "images" / path
    second = manifest_path.parent / path
    return first if first.is_file() else second


def check_wcs(frame, image_path, mask_header):
    info = validate_fits(image_path, frame)
    with fits.open(image_path, memmap=True, checksum=True) as hdul:
        image_header = hdul[info["science_hdu"]].header.copy()
        image_header.extend(hdul[0].header, unique=True)
    image_wcs = WCS(image_header).celestial
    mask_wcs = WCS(mask_header).celestial
    ny, nx = info["shape"]
    pixels = np.array([[0., 0.], [nx - 1., 0.], [0., ny - 1.],
                       [nx - 1., ny - 1.], [(nx - 1.) / 2., (ny - 1.) / 2.]])
    sky = image_wcs.all_pix2world(pixels, 0)
    recovered = mask_wcs.all_world2pix(sky, 0)
    residual = np.hypot(*(recovered - pixels).T)
    if not np.isfinite(residual).all() or float(np.max(residual)) > 1e-5:
        raise AssertionError(f"Mask/image WCS mismatch: {frame['id']} max={np.max(residual)} pixel")
    return info, float(np.max(residual))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--prepared", type=Path, required=True)
    parser.add_argument("--expected", type=int, default=800)
    parser.add_argument("--reference-mask-records", type=Path,
                        default=Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/"
                                    "prestack_cluster_masks_2026-10-16_01/records.tsv"))
    args = parser.parse_args()
    if not platform.node().lower().startswith("lageunha"):
        parser.error("Run this file audit on lageunha CPU")
    manifest_path = args.manifest.resolve()
    prepared = args.prepared.resolve()
    report_path = prepared / "artifact_audit.json"
    if report_path.exists():
        parser.error(f"Refusing to overwrite existing audit: {report_path}")
    manifest = read_json(manifest_path)
    summary_path = prepared / "summary.json"
    summary = read_json(summary_path)
    frames = manifest.get("clusters", [])
    records = summary.get("records", [])
    if (len(frames) != args.expected or summary.get("status") != "complete"
            or summary.get("n_requested") != args.expected
            or summary.get("n_prepared") != args.expected or summary.get("errors")):
        parser.error("Input manifest or preparation summary is not a complete expected sample")
    if summary.get("science_measurement") is not False:
        parser.error("Prepared output unexpectedly activates a science measurement")
    if sha256(manifest_path) != summary.get("manifest_sha256"):
        parser.error("Preparation manifest checksum differs")
    script_path = Path(__file__).with_name("prepare_detection_psf.py")
    if sha256(script_path) != summary.get("preparation_script_sha256"):
        parser.error("Preparation source code checksum differs")
    frame_by_id = {str(row["id"]): row for row in frames}
    record_by_id = {str(row["id"]): row for row in records}
    if len(frame_by_id) != args.expected or len(record_by_id) != args.expected:
        parser.error("Duplicate or missing cluster identifiers")
    if set(frame_by_id) != set(record_by_id):
        parser.error("Manifest and prepared cluster ID sets differ")

    reference = {}
    if args.reference_mask_records.is_file():
        with args.reference_mask_records.open(encoding="utf-8", newline="") as stream:
            reference = {str(row["id"]): row for row in csv.DictReader(stream, delimiter="\t")}
        if set(reference) != set(frame_by_id):
            parser.error("Reference image-hash table does not cover the same 800 clusters")

    mask_hashes = {}
    table_hashes = {}
    records_audit = []
    by_filter = {"i": [], "i2": []}
    for index, ident in enumerate(sorted(frame_by_id, key=lambda x: int(x)), 1):
        frame = dict(frame_by_id[ident])
        record = record_by_id[ident]
        image_path = image_path_for(frame, manifest_path)
        if not image_path.is_file():
            raise FileNotFoundError(image_path)
        frame["image"] = str(image_path.resolve())
        subdir = prepared / ("cluster_" + ident)
        preparation_path = subdir / "preparation.json"
        mask_path = subdir / "masks.fits.gz"
        detections_path = subdir / "detections.tsv"
        members_path = subdir / "members.tsv"
        psf_path = subdir / "psf_candidates.tsv"
        coverage_path = subdir / "mask_radial_coverage.tsv"
        for path in (preparation_path, mask_path, detections_path, members_path, psf_path,
                     coverage_path, subdir / "review.png", subdir / "bcg_review.png"):
            if not path.is_file() or path.stat().st_size == 0:
                raise FileNotFoundError(f"Missing/empty artifact: {path}")

        prep = read_json(preparation_path)
        if (str(prep.get("id")) != ident or prep.get("science_measurement") is not False
                or prep.get("mask_ready_for_science") is not False
                or prep.get("psf_wings_validated") is not False):
            raise AssertionError(f"Unexpected preparation gate for {ident}")
        if prep.get("source_image") != str(image_path.resolve()):
            raise AssertionError(f"Image provenance differs for {ident}")

        detection_count = line_count(detections_path) - 1
        member_count = line_count(members_path) - 1
        psf_count = line_count(psf_path) - 1
        coverage_count = line_count(coverage_path) - 1
        if detection_count != int(record["n_detections"]):
            raise AssertionError(f"Detection table row count mismatch for {ident}")
        if member_count != int(record["n_members"]):
            raise AssertionError(f"Member table row count mismatch for {ident}")
        if psf_count != int(record["n_psf_candidates"]):
            raise AssertionError(f"PSF candidate row count mismatch for {ident}")
        if coverage_count != 60:
            raise AssertionError(f"Radial-coverage table must contain 60 bins for {ident}")

        with fits.open(mask_path, memmap=False, checksum=True) as hdul:
            names = [hdu.name for hdu in hdul]
            expected_names = ["PRIMARY", "SEGMENTATION", "BASE_MASK", "MASK_K2P5",
                              "MASK_K3P5", "MASK_K5P0"]
            if names != expected_names:
                raise AssertionError(f"Unexpected mask HDUs for {ident}: {names}")
            segmentation = np.asarray(hdul["SEGMENTATION"].data)
            base = np.asarray(hdul["BASE_MASK"].data)
            masks = [np.asarray(hdul[name].data) for name in expected_names[3:]]
            headers = [hdul["SEGMENTATION"].header.copy(), hdul["BASE_MASK"].header.copy()]
            header = headers[0]
            if any(arr.shape != segmentation.shape for arr in [base, *masks]):
                raise AssertionError(f"Mask plane shape mismatch for {ident}")
            if segmentation.ndim != 2 or np.any(segmentation < 0) or int(np.max(segmentation)) > detection_count:
                raise AssertionError(f"Invalid segmentation labels for {ident}")
            if not np.isin(base, [0, 1]).all() or any(not np.isin(arr, [0, 1]).all() for arr in masks):
                raise AssertionError(f"Mask plane contains non-binary values for {ident}")
            base_bool = base.astype(bool)
            masks_bool = [arr.astype(bool) for arr in masks]
            if not all(np.all(mask[base_bool]) for mask in masks_bool):
                raise AssertionError(f"Proposed mask cleared a base exclusion for {ident}")
            if not (np.all(~masks_bool[0] | masks_bool[1])
                    and np.all(~masks_bool[1] | masks_bool[2])):
                raise AssertionError(f"Kron masks are not nested for {ident}")
            if header.get("SCIVAL") is not False or header.get("BCGOK") != (record["bcg_status"] == "unique_positional_candidate"):
                raise AssertionError(f"Mask science/BCG flags differ for {ident}")
        info, wcs_residual = check_wcs(frame, image_path, header)
        if tuple(info["shape"]) != tuple(segmentation.shape):
            raise AssertionError(f"Mask/image array dimensions differ for {ident}")
        if record.get("bcg_status") not in ("unique_positional_candidate", "missing_or_ambiguous"):
            raise AssertionError(f"Unknown BCG state for {ident}")
        by_filter[info["coadd_header_filter"]].append({"id": ident, "z": float(frame["z"])})

        mask_hash = sha256(mask_path)
        mask_hashes[str(mask_path)] = mask_hash
        for path in (detections_path, members_path, psf_path, coverage_path, preparation_path):
            table_hashes[str(path)] = sha256(path)
        reference_match = None
        if ident in reference:
            ref = reference[ident]
            if str(Path(ref["image"]).resolve()) != str(image_path.resolve()):
                raise AssertionError(f"Reference image path differs for {ident}")
            if int(ref["image_bytes"]) != image_path.stat().st_size:
                raise AssertionError(f"Reference image size differs for {ident}")
            reference_match = True
        records_audit.append({
            "id": ident, "filter": info["coadd_header_filter"],
            "image_path": str(image_path.resolve()), "image_bytes": image_path.stat().st_size,
            "reference_image_path_size_match": reference_match,
            "mask_file": str(mask_path), "mask_sha256": mask_hash,
            "n_detections": detection_count, "n_members": member_count,
            "n_psf_candidates": psf_count, "bcg_status": record["bcg_status"],
            "segmentation_shape": list(segmentation.shape),
            "segmentation_max_label": int(np.max(segmentation)),
            "mask_wcs_max_residual_pixel": wcs_residual,
        })
        if index % 100 == 0:
            print(f"Audited {index}/{args.expected}", flush=True)

    # Recompute original quality + bright-object masks on a deterministic redshift-stratified sample.
    sampled = []
    for band in ("i", "i2"):
        rows = sorted(by_filter[band], key=lambda row: (row["z"], int(row["id"])))
        indices = np.unique(np.linspace(0, len(rows) - 1, min(8, len(rows))).round().astype(int))
        sampled.extend((band, rows[int(index)]["id"]) for index in indices)
    reference_mask_audits = []
    for band, ident in sampled:
        frame = dict(frame_by_id[ident])
        image_path = image_path_for(frame, manifest_path)
        frame["image"] = str(image_path.resolve())
        info = validate_fits(image_path, frame)
        with fits.open(image_path, memmap=False, checksum=True) as hdul:
            image = np.asarray(hdul[info["science_hdu"]].data, dtype=np.float32)
            variance = np.asarray(hdul[info["variance_hdu"]].data, dtype=np.float32)
            bits = np.asarray(hdul[info["mask_hdu"]].data, dtype=np.int32)
        invalid = ~np.isfinite(image) | ~np.isfinite(variance) | (variance <= 0)
        bad_bits = sum(1 << bit for bit in BAD_BIT_POSITIONS) | BRIGHT_OBJECT_BIT
        expected_base = ((bits & bad_bits) != 0) | invalid
        with fits.open(prepared / ("cluster_" + ident) / "masks.fits.gz", memmap=False) as hdul:
            actual_base = np.asarray(hdul["BASE_MASK"].data, dtype=bool)
        if not np.array_equal(expected_base, actual_base):
            raise AssertionError(f"Base mask differs from original HSC quality bits for {ident}")
        reference_mask_audits.append({"id": ident, "filter": band,
                                      "base_mask_matches_original_bits": True})
        if ident in reference:
            ref = reference[ident]
            if sha256(image_path) != ref["image_sha256"]:
                raise AssertionError(f"Input image checksum differs from prior verified record for {ident}")

    report = {
        "schema": "hsc-satellite-detection-psf-audit-v1",
        "passed": True, "science_measurement": False, "mask_ready_for_science": False,
        "n_frames": len(records_audit), "n_unique_ids": len({r["id"] for r in records_audit}),
        "n_detection_rows": sum(r["n_detections"] for r in records_audit),
        "n_member_positions": sum(r["n_members"] for r in records_audit),
        "n_psf_candidates": sum(r["n_psf_candidates"] for r in records_audit),
        "n_unique_center_candidates": sum(r["bcg_status"] == "unique_positional_candidate" for r in records_audit),
        "n_psf_candidate_frames_ge4": sum(r["n_psf_candidates"] >= 4 for r in records_audit),
        "n_psf_candidate_frames_zero": sum(r["n_psf_candidates"] == 0 for r in records_audit),
        "filter_counts": {band: len(rows) for band, rows in by_filter.items()},
        "input_manifest_sha256": sha256(manifest_path),
        "prepare_script_sha256": sha256(script_path),
        "mask_sha256": mask_hashes,
        "table_sha256": table_hashes,
        "checks": ["800 unique IDs/status", "per-cluster table row counts", "FITS HDUs and array shape",
                   "binary/nested masks and base-exclusion preservation", "segmentation label range",
                   "image-to-mask WCS agreement on all 800", "original quality-bit parity on stratified sample",
                   "cross-check input paths/sizes against prior full-sample mask manifest",
                   "PSF/BCG/science states remain provisional"],
        "base_mask_sample": reference_mask_audits,
        "records": records_audit,
    }
    report_path.write_text(json.dumps(report, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
                           encoding="utf-8")
    print(json.dumps({key: report[key] for key in ("passed", "n_frames", "n_detection_rows",
                                                   "n_unique_center_candidates", "n_psf_candidate_frames_ge4",
                                                   "n_psf_candidate_frames_zero", "filter_counts")},
                     ensure_ascii=False), flush=True)
    print("Base mask parity sample", [row["id"] for row in reference_mask_audits], flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
