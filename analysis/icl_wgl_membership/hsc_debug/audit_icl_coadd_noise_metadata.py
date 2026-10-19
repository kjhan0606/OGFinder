#!/usr/bin/env python3
"""Audit coadd-header support for a source-photon-noise injection model.

The audit follows the 57 HSC target frames in the independently checked
position-dependent-PSF injection run. It reads FITS headers only; it does not
alter or download images and makes no science measurement.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

from astropy.io import fits


DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
DEFAULT_BATCH = DATA / "joint_icl_satellite_psf_redetection_2026-10-19_02"
DEFAULT_AUDIT = DEFAULT_BATCH / "psf_redetection_audit.json"
DEFAULT_MASK_RECORDS = DATA / "prestack_cluster_masks_2026-10-16_01/records.tsv"
DEFAULT_OUTPUT = DATA / "icl_coadd_noise_metadata_2026-10-19_01/summary.json"
KEYWORDS = ("FLUXMAG0", "GAIN", "EXPTIME", "DARKTIME", "BUNIT", "SKYSIGMA", "NCOMBINE")


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def summarize_keyword_presence(frame_headers):
    """Summarize keyword availability across frame -> HDU -> header mappings."""
    if not frame_headers:
        raise ValueError("at least one frame header record is required")
    counts = {key: 0 for key in KEYWORDS}
    examples = {key: set() for key in KEYWORDS}
    for record in frame_headers:
        headers = record.get("headers", [])
        combined = {}
        for header in headers:
            combined.update(header)
        for key in KEYWORDS:
            if key in combined:
                counts[key] += 1
                examples[key].add(str(combined[key]))
    return {
        key: {"frames_with_keyword": counts[key],
              "frames_missing_keyword": len(frame_headers) - counts[key],
              "value_examples": sorted(examples[key])[:5]}
        for key in KEYWORDS
    }


def inspect_frame(path: Path):
    """Read selected metadata from HDU headers without loading pixel arrays."""
    if not path.is_file():
        raise FileNotFoundError(path)
    headers = []
    image_shapes = []
    with fits.open(path, memmap=True, lazy_load_hdus=True) as hdus:
        for index, hdu in enumerate(hdus):
            header = hdu.header
            values = {key: header[key] for key in KEYWORDS if key in header}
            headers.append({"hdu": index, "name": str(hdu.name), "keywords": values})
            naxis = int(header.get("NAXIS", 0))
            if naxis >= 2:
                shape = tuple(int(header[f"NAXIS{axis}"]) for axis in range(naxis, 0, -1))
                image_shapes.append({"hdu": index, "shape": shape})
    merged = {}
    for item in headers:
        merged.update(item["keywords"])
    fluxmag0 = merged.get("FLUXMAG0")
    if fluxmag0 is None or not math.isfinite(float(fluxmag0)) or float(fluxmag0) <= 0:
        raise ValueError(f"missing or invalid FLUXMAG0: {path}")
    if len(image_shapes) < 3:
        raise ValueError(f"expected science, mask and variance image HDUs: {path}")
    return {"path": str(path), "file_bytes": path.stat().st_size,
            "headers": [x["keywords"] for x in headers],
            "hdu_headers": headers, "image_shapes": image_shapes,
            "header_only": True}


def build_audit(batch_dir: Path, audit_path: Path, mask_records_path: Path):
    audit = json.loads(audit_path.read_text(encoding="utf-8"))
    if (audit.get("status") != "verified" or audit.get("science_measurement") is not False
            or audit.get("n_frame_results_checked") != 57 or audit.get("n_errors_reported_by_batch") != 0):
        raise ValueError("PSF injection batch is not the audited 57-frame run")
    with mask_records_path.open(encoding="utf-8", newline="") as stream:
        mask_records = list(csv.DictReader(stream, delimiter="\t"))
    by_id = {str(row["id"]): row for row in mask_records}
    if len(by_id) != len(mask_records):
        raise ValueError("duplicate cluster IDs in target mask records")

    frame_headers = []
    frame_keys = set()
    image_hashes = {}
    for approved in audit["records"]:
        ident, band = str(approved["cluster_id"]), str(approved["filter"])
        key = (ident, band)
        if key in frame_keys:
            raise ValueError(f"duplicate approved frame in PSF audit: {key}")
        frame_keys.add(key)
        result_path = batch_dir / f"cluster_{ident}.json"
        result = json.loads(result_path.read_text(encoding="utf-8"))
        if result.get("status") != "complete" or result.get("science_measurement") is not False:
            raise ValueError(f"frame result is incomplete or promoted: {result_path}")
        if str(result.get("filter")) != band:
            raise ValueError(f"frame filter differs from approved record: {result_path}")
        mask_record = by_id.get(ident)
        if mask_record is None:
            raise ValueError(f"target image record is missing: {ident}")
        image_path = Path(mask_record["image"])
        expected_hash = str(mask_record["image_sha256"])
        if result.get("input_sha256", {}).get(str(image_path)) != expected_hash:
            raise ValueError(f"frame result does not identify the mask-record image hash: {ident}")
        if image_path.name != f"cluster_{ident}.fits":
            raise ValueError(f"unexpected target image filename for {ident}: {image_path.name}")
        frame = inspect_frame(image_path)
        frame["cluster_id"] = ident
        frame["filter"] = band
        frame["source_image_sha256_from_audited_run"] = expected_hash
        frame_headers.append(frame)
        image_hashes[ident] = expected_hash
    if len(frame_keys) != 57:
        raise ValueError(f"expected 57 approved PSF-injection frames; got {len(frame_keys)}")
    keyword_summary = summarize_keyword_presence(frame_headers)
    missing_any = any(keyword_summary[key]["frames_missing_keyword"] for key in ("GAIN", "EXPTIME"))
    return {
        "schema": "astrafek-icl-coadd-noise-metadata-audit-v1",
        "status": "verified_metadata_only",
        "science_measurement": False,
        "n_frames": len(frame_headers),
        "source_psf_injection_audit": str(audit_path.resolve()),
        "source_psf_injection_audit_sha256": sha256(audit_path),
        "target_mask_records": str(mask_records_path.resolve()),
        "target_mask_records_sha256": sha256(mask_records_path),
        "audit_script_sha256": sha256(Path(__file__)),
        "keyword_summary": keyword_summary,
        "poisson_noise_model_from_delivered_coadd_headers_available": not missing_any,
        "variance_plane_present_in_all_frames": all(
            {x["hdu"] for x in frame["image_shapes"]} >= {1, 2, 3}
            and len({tuple(x["shape"]) for x in frame["image_shapes"] if x["hdu"] in {1, 2, 3}}) == 1
            for frame in frame_headers
        ),
        "interpretation": (
            "Delivered coadds contain FLUXMAG0 and a variance image, but lack GAIN and EXPTIME; "
            "the variance plane alone does not specify the source-photon variance for added light. "
            "Do not invent a gain or treat the unchanged variance plane as a photon-noise model."
        ),
        "limitations": [
            "header metadata only; pixel values and variance statistics were not read",
            "noise inputs are representative of the 57 audited target frames, not all 800 targets or controls",
            "this audit does not test noise correlation, source Poisson noise, or upstream HSC sky transfer",
        ],
        "frame_headers": frame_headers,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--batch-dir", type=Path, default=DEFAULT_BATCH)
    parser.add_argument("--psf-audit", type=Path, default=DEFAULT_AUDIT)
    parser.add_argument("--mask-records", type=Path, default=DEFAULT_MASK_RECORDS)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    output = args.output.resolve()
    if output.exists():
        parser.error(f"refusing to overwrite existing output: {output}")
    result = build_audit(args.batch_dir.resolve(), args.psf_audit.resolve(),
                         args.mask_records.resolve())
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("x", encoding="utf-8") as stream:
        json.dump(result, stream, ensure_ascii=False, indent=2, allow_nan=False)
        stream.write("\n")
    print(json.dumps({
        "status": result["status"], "frames": result["n_frames"],
        "gain_present": result["keyword_summary"]["GAIN"]["frames_with_keyword"],
        "exptime_present": result["keyword_summary"]["EXPTIME"]["frames_with_keyword"],
        "variance_present": result["variance_plane_present_in_all_frames"],
        "science_measurement": False, "output": str(output),
    }, ensure_ascii=False))


if __name__ == "__main__":
    main()
