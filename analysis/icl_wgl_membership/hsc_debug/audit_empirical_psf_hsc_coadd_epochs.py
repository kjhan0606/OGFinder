#!/usr/bin/env python3
"""Audit whether existing HSC coadds expose epochs for Gaia position propagation.

Only FITS headers are read. No image arrays, HSC services, or external data are
accessed. A coadd observation epoch is considered available only when a
standard observation date/MJD card is present in one of its HDUs. This is an
audit of Gaia-to-HSC cross-match timing; the empirical PSF stamps themselves
are centered on HSC detection X_IMAGE/Y_IMAGE, not on Gaia coordinates.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import time

from astropy.io import fits


EPOCH_KEYS = (
    "DATE-OBS", "DATE-AVG", "DATE-BEG", "DATE-END",
    "MJD-OBS", "MJD-AVG", "MJD-MID", "MJD-BEG", "MJD-END",
)
DEFAULT_SAMPLE = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
DEFAULT_LEDGER = DEFAULT_SAMPLE / "empirical_psf_gaia_full_sample_2026-10-20_02/star_ledger.tsv"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_frame_ids(ledger: Path) -> list[str]:
    with Path(ledger).open(encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream, delimiter="\t"))
    if not rows or not {"cluster_id", "stamp"}.issubset(rows[0]):
        raise ValueError("Gaia star ledger is empty or missing required columns")
    stamped = [row for row in rows if row["stamp"].strip()]
    if not stamped:
        raise ValueError("Gaia star ledger has no verified stamps")
    frame_ids = sorted({row["cluster_id"].strip() for row in stamped}, key=lambda v: int(v))
    if any(not value.isdigit() for value in frame_ids):
        raise ValueError("non-numeric cluster_id in Gaia star ledger")
    return frame_ids


def count_stamped_rows(ledger: Path) -> int:
    with Path(ledger).open(encoding="utf-8", newline="") as stream:
        return sum(1 for row in csv.DictReader(stream, delimiter="\t")
                   if row.get("stamp", "").strip())


def audit_fits_header(path: Path, cluster_id: str) -> dict:
    """Read all HDU headers, collecting epoch cards and a header-only digest."""
    epoch_values = []
    header_digest = hashlib.sha256()
    n_hdus = 0
    with fits.open(path, mode="readonly", memmap=True, lazy_load_hdus=True) as hdus:
        for index, hdu in enumerate(hdus):
            n_hdus += 1
            header_digest.update(f"HDU={index};EXTNAME={hdu.name}\n".encode("utf-8"))
            header_digest.update(hdu.header.tostring(endcard=True, padding=False).encode("ascii", "replace"))
            header_digest.update(b"\n")
            for key in EPOCH_KEYS:
                if key in hdu.header:
                    epoch_values.append({
                        "hdu": index,
                        "extname": hdu.name,
                        "key": key,
                        "value": str(hdu.header[key]),
                    })
    if n_hdus == 0:
        raise ValueError(f"FITS has no HDUs: {path}")
    return {
        "cluster_id": cluster_id,
        "path": str(path.resolve()),
        "file_bytes": path.stat().st_size,
        "n_hdus": n_hdus,
        "epoch_cards": json.dumps(epoch_values, ensure_ascii=False, separators=(",", ":")),
        "has_standard_epoch": bool(epoch_values),
        "header_sha256": header_digest.hexdigest(),
    }


def run(sample_root: Path, ledger: Path, output: Path) -> dict:
    sample_root, ledger, output = map(Path, (sample_root, ledger, output))
    if output.exists():
        raise FileExistsError(f"refusing to overwrite epoch audit: {output}")
    frame_ids = read_frame_ids(ledger)
    rows = []
    for cluster_id in frame_ids:
        image_path = sample_root / "download_01/images" / f"cluster_{cluster_id}.fits"
        if not image_path.is_file():
            raise FileNotFoundError(f"coadd referenced by Gaia ledger is missing: {image_path}")
        rows.append(audit_fits_header(image_path, cluster_id))

    output.mkdir(parents=True, exist_ok=False)
    table_path = output / "coadd_epoch_headers.tsv"
    fields = list(rows[0])
    with table_path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    n_with_epoch = sum(row["has_standard_epoch"] for row in rows)
    summary = {
        "schema": "astrafek-empirical-psf-coadd-epoch-header-audit-v1",
        "status": "complete_header_only_audit",
        "n_gaia_stamps": count_stamped_rows(ledger),
        "n_unique_target_frames": len(rows),
        "n_frames_with_standard_epoch": n_with_epoch,
        "n_frames_missing_standard_epoch": len(rows) - n_with_epoch,
        "standard_epoch_keys": list(EPOCH_KEYS),
        "can_propagate_gaia_to_actual_coadd_epoch_from_existing_fits": n_with_epoch == len(rows),
        "stamp_center_source": "HSC detection X_IMAGE/Y_IMAGE in native coadd pixels",
        "gaia_epoch_propagation_required_to_recenter_stamp": False,
        "gaia_match_identity_epoch_exactly_verified": n_with_epoch == len(rows),
        "image_pixels_read": False,
        "hsc_requests_made": False,
        "inputs_sha256": {str(ledger.resolve()): sha256(ledger)},
        "header_table_sha256": sha256(table_path),
        "script_sha256": sha256(Path(__file__)),
        "limitations": [
            "The audit checks only standard observation date/MJD cards in the supplied coadd FITS HDUs.",
            "Absence of an epoch card does not imply that HSC exposure metadata are unavailable elsewhere.",
            "Even a coadd DATE/MJD card may not represent a per-object effective epoch for a multi-visit coadd.",
            "The HSC detection centroid centers the extracted stamp, so this missing epoch does not prevent image-space centering.",
            "The missing epoch leaves Gaia cross-match identity at the actual HSC epoch unverified by this audit.",
            "No Gaia motion was applied and no PSF profile was measured.",
        ],
        "science_measurement": False,
        "psf_approved": False,
        "completed_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    }
    summary_path = output / "summary.json"
    summary_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-root", type=Path, default=DEFAULT_SAMPLE)
    parser.add_argument("--ledger", type=Path, default=DEFAULT_LEDGER)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = run(args.sample_root.resolve(), args.ledger.resolve(), args.output.resolve())
    print(json.dumps({key: result[key] for key in (
        "status", "n_unique_target_frames", "n_frames_with_standard_epoch",
        "n_frames_missing_standard_epoch", "can_propagate_gaia_to_actual_coadd_epoch_from_existing_fits",
        "image_pixels_read", "science_measurement",
    )}, ensure_ascii=False, indent=2))
    print(f"summary_sha256={sha256(args.output / 'summary.json')}")


if __name__ == "__main__":
    main()
