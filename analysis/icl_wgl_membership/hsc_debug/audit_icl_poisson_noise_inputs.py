#!/usr/bin/env python3
"""Audit whether downloaded ICL coadds carry metadata for source shot-noise injection.

This is a read-only header audit. It checks the saved target/control download
records and the IMAGE, MASK, and VARIANCE FITS planes; it does not read image
pixels, contact HSC, or modify the source files. A coadd variance plane alone
does not define the effective gain needed to add a synthetic source's Poisson
variance in the coadd's native units.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

from astropy.io import fits


DEFAULT_DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
DEFAULT_TARGET_STATUS = DEFAULT_DATA / "download_01/download_status.json"
DEFAULT_CONTROL_STATUS = DEFAULT_DATA / "control_download_01/download_status.json"
DEFAULT_OUTPUT = DEFAULT_DATA / "icl_poisson_noise_input_audit_2026-11-04_01/summary.json"

# These names identify explicit scalar detector/calibration metadata, not the
# VARIANCE plane itself. An effective gain map may require upstream visit data.
GAIN_KEYS = ("GAIN", "GAIN1", "GAIN2", "EFFGAIN", "COADDGAIN")
EXPOSURE_KEYS = ("EXPTIME", "EXPOSURE", "EFFEXPT", "SUMEXPT", "NEXP", "NINPUT")
PLANE_TYPES = {"IMAGE": "science_hdu", "MASK": "mask_hdu", "VARIANCE": "variance_hdu"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def audit_group(status_path: Path, image_dir: Path, label: str) -> list[dict]:
    status = json.loads(status_path.read_text(encoding="utf-8"))
    if status.get("status") != "complete":
        raise ValueError(f"{label} download status is not complete")
    records = status.get("completed")
    if not isinstance(records, list) or len(records) != int(status["n_expected"]):
        raise ValueError(f"{label} completed records differ from n_expected")

    output = []
    for record in records:
        review = record.get("review") or {}
        path = image_dir / str(record["file"])
        item = {
            "sample": label,
            "id": str(record.get("id", "")),
            "filter": str(review.get("coadd_header_filter", "")),
            "path": str(path),
            "exists": path.is_file(),
            "expected_bytes": int(record.get("bytes", -1)),
            "actual_bytes": path.stat().st_size if path.is_file() else None,
            "review_hdus": {name: review.get(field) for name, field in PLANE_TYPES.items()},
            "ext_types": {},
            "gain_header_keys": [],
            "exposure_header_keys": [],
            "plane_shapes": {},
            "plane_mapping_matches_review": False,
        }
        if not path.is_file():
            output.append(item)
            continue

        with fits.open(path, memmap=True, lazy_load_hdus=False) as hdus:
            primary = hdus[0].header
            kinds: dict[str, int] = {}
            gain_keys: set[str] = set()
            exposure_keys: set[str] = set()
            for index, hdu in enumerate(hdus):
                header = hdu.header
                kind = str(header.get("EXTTYPE", "")).strip().upper()
                if kind in PLANE_TYPES:
                    if kind in kinds:
                        raise ValueError(f"{path}: duplicate {kind} plane")
                    kinds[kind] = index
                    ndim = int(header.get("NAXIS", 0))
                    if ndim:
                        item["plane_shapes"][kind] = [
                            int(header[f"NAXIS{axis}"]) for axis in range(ndim, 0, -1)
                        ]
                for key in GAIN_KEYS:
                    if key in header:
                        gain_keys.add(key)
                for key in EXPOSURE_KEYS:
                    if key in header:
                        exposure_keys.add(key)
            # The primary header may hold inherited calibration metadata.
            for key in GAIN_KEYS:
                if key in primary:
                    gain_keys.add(key)
            for key in EXPOSURE_KEYS:
                if key in primary:
                    exposure_keys.add(key)
            item["ext_types"] = kinds
            item["gain_header_keys"] = sorted(gain_keys)
            item["exposure_header_keys"] = sorted(exposure_keys)
            item["primary_fluxmag0_present"] = bool(primary.get("FLUXMAG0"))
            item["plane_mapping_matches_review"] = all(
                kinds.get(kind) == review.get(field)
                for kind, field in PLANE_TYPES.items()
            )
        output.append(item)
    return output


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target-status", type=Path, default=DEFAULT_TARGET_STATUS)
    parser.add_argument("--control-status", type=Path, default=DEFAULT_CONTROL_STATUS)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    if args.output.exists():
        raise FileExistsError(f"refusing to overwrite {args.output}")

    records = audit_group(args.target_status, args.target_status.parent / "images", "target")
    records += audit_group(args.control_status, args.control_status.parent / "images", "control")
    groups = {}
    for label in ("target", "control"):
        subset = [row for row in records if row["sample"] == label]
        by_filter = {}
        for filter_name in sorted({row["filter"] for row in subset}):
            rows = [row for row in subset if row["filter"] == filter_name]
            by_filter[filter_name] = {
                "n_images": len(rows),
                "n_files_present_and_size_matched": sum(
                    row["exists"] and row["actual_bytes"] == row["expected_bytes"]
                    for row in rows
                ),
                "n_plane_mappings_matched": sum(row["plane_mapping_matches_review"] for row in rows),
                "n_with_variance_plane": sum("VARIANCE" in row["ext_types"] for row in rows),
                "n_with_explicit_gain_metadata": sum(bool(row["gain_header_keys"]) for row in rows),
                "n_with_exposure_metadata": sum(bool(row["exposure_header_keys"]) for row in rows),
                "gain_keys_seen": sorted({key for row in rows for key in row["gain_header_keys"]}),
                "exposure_keys_seen": sorted({key for row in rows for key in row["exposure_header_keys"]}),
            }
        groups[label] = {"n_images": len(subset), "by_filter": by_filter}

    missing_or_mismatch = [
        {"sample": row["sample"], "id": row["id"], "path": row["path"],
         "exists": row["exists"], "expected_bytes": row["expected_bytes"],
         "actual_bytes": row["actual_bytes"],
         "plane_mapping_matches_review": row["plane_mapping_matches_review"]}
        for row in records
        if (not row["exists"] or row["actual_bytes"] != row["expected_bytes"]
            or not row["plane_mapping_matches_review"])
    ]
    all_have_gain = bool(records) and all(row["gain_header_keys"] for row in records)
    summary = {
        "schema": "astrafek-icl-poisson-noise-input-audit-v1",
        "status": "complete",
        "inputs": {
            str(args.target_status): sha256(args.target_status),
            str(args.control_status): sha256(args.control_status),
        },
        "n_files_audited": len(records),
        "groups": groups,
        "metadata_gate": {
            "all_files_match_saved_byte_size_and_plane_mapping": not missing_or_mismatch,
            "all_files_have_explicit_gain_metadata": all_have_gain,
            "poisson_variance_addition_supported_from_current_coadd_headers": all_have_gain,
            "reason": (
                "Current coadds include a variance plane and FLUXMAG0, but no explicit gain metadata; "
                "the coadd variance does not by itself specify the effective gain for a synthetic source."
                if not all_have_gain else
                "Gain-like header fields exist; their semantics and validity for coadd shot-noise propagation "
                "still require calibration review."
            ),
        },
        "missing_or_mismatch": missing_or_mismatch,
        "pixel_data_read": False,
        "network_used": False,
        "credentials_saved": False,
        "images_modified": False,
        "science_measurement": False,
    }
    if missing_or_mismatch:
        summary["status"] = "audit_incomplete"
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0 if not missing_or_mismatch else 2


if __name__ == "__main__":
    raise SystemExit(main())
