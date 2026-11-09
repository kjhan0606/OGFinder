#!/usr/bin/env python3
"""Record pixel-mask evidence for S23B frames rejected by mask preparation.

Local-only audit of named frames in the pinned approved target download. It
does not run detection, estimate surface brightness, or promote a science
measurement.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from astropy.io import fits

from download_pilot import BAD_BIT_POSITIONS
from prepare_s23b_icl_target_image_manifest import sha256


EXPECTED_MANIFEST_SHA256 = "e487582ba40cd3e58da67704242db4d7219e381fb08f4b33fed47ce4ff268b44"
INTRP_BIT = 2
NO_DATA_BIT = 8


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def audit(download_dir: Path, manifest_path: Path, ids: list[str], output: Path) -> dict:
    download_dir, manifest_path, output = map(Path.resolve, (download_dir, manifest_path, output))
    if output.exists():
        raise FileExistsError(f"refusing to overwrite {output}")
    manifest_sha = sha256(manifest_path)
    if manifest_sha != EXPECTED_MANIFEST_SHA256:
        raise ValueError("target manifest is not the pinned approved manifest")
    manifest = read_json(manifest_path)
    status_path = download_dir / "download_status.json"
    downloaded_path = download_dir / "download_manifest.json"
    status, downloaded = map(read_json, (status_path, downloaded_path))
    if (status.get("status") != "complete"
            or status.get("manifest_sha256") != manifest_sha
            or status.get("science_measurement") is not False
            or status.get("credentials_saved") is not False
            or downloaded.get("science_measurement") is not False):
        raise ValueError("download is incomplete, unpinned, or outside the approved scope")

    approved = {str(row["id"]): row for row in manifest.get("clusters", [])}
    completed = {str(row["id"]): row for row in status.get("completed", [])}
    products = {str(row["id"]): row for row in downloaded.get("clusters", [])}
    if len(approved) != 1043 or set(approved) != set(completed) or set(approved) != set(products):
        raise ValueError("download identity sets do not match the approved 1,043-frame sample")
    if len(ids) != len(set(ids)) or not set(ids) <= set(approved):
        raise ValueError("requested frame IDs are duplicated or outside the approved sample")

    rows = []
    bad_bits = sum(1 << bit for bit in BAD_BIT_POSITIONS)
    for ident in ids:
        item, product = completed[ident], products[ident]
        path = download_dir / "images" / str(item["file"])
        if (Path(item["file"]).name != item["file"] or not path.is_file()
                or path.stat().st_size != int(item["bytes"])
                or sha256(path) != item.get("sha256")):
            raise ValueError(f"{ident}: file identity, size, or checksum mismatch")
        review = item.get("review", {})
        with fits.open(path, memmap=False) as hdul:
            image = np.asarray(hdul[int(review["science_hdu"])].data)
            mask = np.asarray(hdul[int(review["mask_hdu"])].data, dtype=np.int32)
            variance = np.asarray(hdul[int(review["variance_hdu"])].data)
        if image.shape != mask.shape or image.shape != variance.shape:
            raise ValueError(f"{ident}: IMAGE/MASK/VARIANCE shapes differ")
        valid = np.isfinite(image) & np.isfinite(variance) & (variance > 0)
        valid &= (mask & bad_bits) == 0
        planes = {
            "MP_INTRP": bool(np.all((mask & (1 << INTRP_BIT)) != 0)),
            "MP_NO_DATA": bool(np.all((mask & (1 << NO_DATA_BIT)) != 0)),
        }
        if valid.any():
            raise ValueError(f"{ident}: requested zero-support exclusion has {int(valid.sum())} valid pixels")
        rows.append({
            "id": ident, "file": path.name, "bytes": path.stat().st_size,
            "sha256": item["sha256"], "shape_y": int(image.shape[0]),
            "shape_x": int(image.shape[1]), "n_pixels": int(image.size),
            "n_finite_image": int(np.isfinite(image).sum()),
            "n_finite_positive_variance": int((np.isfinite(variance) & (variance > 0)).sum()),
            "n_unflagged_finite_positive_variance": int(valid.sum()),
            "entire_frame_has_mask_plane": planes,
            "primary_filter": review.get("coadd_header_filter"),
            "target_ra_deg": approved[ident].get("ra_deg"),
            "target_dec_deg": approved[ident].get("dec_deg"),
            "download_product": product.get("product"),
            "exclusion_reason": "zero unflagged finite positive-variance pixels",
        })

    summary = {
        "schema": "astrafek-s23b-icl-zero-support-frame-audit-v1",
        "status": "complete_zero_pixel_support_audit",
        "n_frames_checked": len(rows), "n_zero_support_frames": len(rows),
        "frames": rows, "science_measurement": False,
        "network_used": False, "credentials_saved": False,
        "method": "finite IMAGE and positive finite VARIANCE, excluding the pinned HSC bad-mask planes",
        "inputs_sha256": {
            str(manifest_path): manifest_sha,
            str(download_dir / "download_status.json"): sha256(status_path),
            str(download_dir / "download_manifest.json"): sha256(downloaded_path),
        },
    }
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"status": summary["status"], "frames": rows}, indent=2))
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--download-dir", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--ids", nargs="+", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    audit(args.download_dir, args.manifest, args.ids, args.output)


if __name__ == "__main__":
    main()
