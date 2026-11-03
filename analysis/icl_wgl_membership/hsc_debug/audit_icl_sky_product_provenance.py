#!/usr/bin/env python3
"""Audit the delivered HSC sky-subtraction product for the ICL image sample.

This is a metadata and file-stat audit only. It reads the saved download
manifests/status records, checks the declared product provenance for every
target and matched-control image, and verifies each completed file's existence
and recorded byte size. It does not hash image contents, read image pixels,
contact HSC, or change any science-measurement state.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from collections import Counter
from pathlib import Path


DEFAULT_SAMPLE = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
PRODUCT = "coadd/bg"
SKY = "global-sky"
EVIDENCE_PHRASE = "before aggressive background subtraction"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def validate_collection(name: str, manifest: dict, status: dict, *, expected_count: int):
    frames = manifest.get("clusters")
    completed = status.get("completed")
    if not isinstance(frames, list) or len(frames) != expected_count:
        raise ValueError(f"{name}: expected {expected_count} manifest frames")
    if (status.get("status") != "complete" or status.get("n_expected") != expected_count
            or status.get("downloaded") != expected_count or not isinstance(completed, list)
            or len(completed) != expected_count):
        raise ValueError(f"{name}: download status is not complete")
    if status.get("credentials_saved") is not False:
        raise ValueError(f"{name}: credentials must not be recorded in status")

    frame_by_id = {str(row.get("id", "")): row for row in frames}
    status_by_id = {str(row.get("id", "")): row for row in completed}
    if (len(frame_by_id) != expected_count or "" in frame_by_id
            or len(status_by_id) != expected_count or "" in status_by_id
            or frame_by_id.keys() != status_by_id.keys()):
        raise ValueError(f"{name}: manifest and completed IDs do not match uniquely")

    images = Path(manifest["downloaded_images"])
    products, skies, filters = Counter(), Counter(), Counter()
    total_bytes = 0
    for ident, frame in frame_by_id.items():
        if frame.get("sky") != SKY or frame.get("product") != PRODUCT:
            raise ValueError(f"{name}/{ident}: unexpected sky/product {frame.get('sky')}/{frame.get('product')}")
        if frame.get("product_reviewed") is not True:
            raise ValueError(f"{name}/{ident}: HSC product provenance is not reviewed")
        evidence = str(frame.get("product_evidence", ""))
        if EVIDENCE_PHRASE not in evidence:
            raise ValueError(f"{name}/{ident}: missing pre-background-subtraction provenance")

        item = status_by_id[ident]
        filename = str(item.get("file", ""))
        if not filename or filename != str(frame.get("image", "")):
            raise ValueError(f"{name}/{ident}: image filename mismatch")
        expected_bytes = int(item.get("bytes", -1))
        expected_sha = str(item.get("sha256", ""))
        if expected_bytes <= 0 or len(expected_sha) != 64:
            raise ValueError(f"{name}/{ident}: invalid recorded size or SHA-256")
        image_path = images / filename
        if not image_path.is_file() or image_path.stat().st_size != expected_bytes:
            raise ValueError(f"{name}/{ident}: image is missing or size differs from completed status")
        products[str(frame["product"])] += 1
        skies[str(frame["sky"])] += 1
        filters[str(frame.get("coadd_header_filter", "unknown"))] += 1
        total_bytes += expected_bytes

    return {
        "n_frames": expected_count,
        "products": dict(products),
        "sky_types": dict(skies),
        "coadd_header_filters": dict(filters),
        "recorded_image_bytes": total_bytes,
        "manifest_sha256": manifest.get("_manifest_sha256"),
        "download_status_sha256": status.get("_status_sha256"),
    }


def audit(sample: Path, output: Path):
    sample, output = Path(sample).resolve(), Path(output).resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite {output}")
    collections = {
        "target": (sample / "download_01/download_manifest.json",
                   sample / "download_01/download_status.json",
                   sample / "sample_manifest.json", 800),
        "matched_control": (sample / "control_download_01/download_manifest.json",
                             sample / "control_download_01/download_status.json",
                             sample / "control_patches_2026-10-14_02/control_download_manifest.json", 798),
    }
    results = {}
    input_hashes = {}
    for name, (manifest_path, status_path, source_manifest_path, expected_count) in collections.items():
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        status = json.loads(status_path.read_text(encoding="utf-8"))
        manifest_hash, status_hash = sha256(manifest_path), sha256(status_path)
        source_manifest_hash = sha256(source_manifest_path)
        if status.get("manifest_sha256") != source_manifest_hash:
            raise ValueError(f"{name}: download status references a different source manifest")
        manifest["_manifest_sha256"] = manifest_hash
        status["_status_sha256"] = status_hash
        results[name] = validate_collection(name, manifest, status,
                                            expected_count=expected_count)
        input_hashes[str(manifest_path)] = manifest_hash
        input_hashes[str(status_path)] = status_hash
        input_hashes[str(source_manifest_path)] = source_manifest_hash

    result = {
        "schema": "astrafek-hsc-icl-sky-product-provenance-audit-v1",
        "status": "all_images_global_sky_coadd_bg",
        "sample": str(sample),
        "collections": results,
        "total_images": sum(row["n_frames"] for row in results.values()),
        "total_recorded_image_bytes": sum(row["recorded_image_bytes"] for row in results.values()),
        "input_sha256": input_hashes,
        "product_interpretation": (
            "The recorded coadd/bg product is the HSC PDR3 global-sky patch product, "
            "mapped by the cutout manual to deepCoadd before aggressive background subtraction. "
            "Global-sky processing still subtracts modeled sky and does not prove zero residual transfer."
        ),
        "file_contents_rehashed": False,
        "image_pixels_read": False,
        "hsc_network_used": False,
        "data_requested": False,
        "science_measurement": False,
        "psf_approved": False,
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
                      encoding="utf-8")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", type=Path, default=DEFAULT_SAMPLE)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = audit(args.sample, args.output)
    print(json.dumps({key: result[key] for key in (
        "status", "total_images", "total_recorded_image_bytes", "collections",
        "file_contents_rehashed", "image_pixels_read", "hsc_network_used",
        "science_measurement")}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
