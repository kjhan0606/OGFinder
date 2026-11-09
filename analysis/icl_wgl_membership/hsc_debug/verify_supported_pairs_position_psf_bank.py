#!/usr/bin/env python3
"""Verify the combined PSF bank for the 224-pair ICL injection sample."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

from position_psf import SupportedPairsPositionPSFBank, sha256


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--injection-plan-dir", type=Path, required=True)
    parser.add_argument("--injection-download-dir", type=Path, required=True)
    parser.add_argument("--test-center-psf", type=Path, required=True)
    parser.add_argument("--target-center-plan-dir", type=Path, required=True)
    parser.add_argument("--target-center-download-dir", type=Path, required=True)
    parser.add_argument("--supported-satellite-plan-dir", type=Path, required=True)
    parser.add_argument("--supported-satellite-download-dir", type=Path, required=True)
    parser.add_argument("--approved-satellite-plan-sha256", required=True,
                        help="exact SHA-256 of the satellite-position plan explicitly approved for download")
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    output = args.output_dir.resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite {output}")

    bank = SupportedPairsPositionPSFBank(
        args.injection_plan_dir.resolve(), args.injection_download_dir.resolve(),
        args.test_center_psf.resolve(), args.target_center_plan_dir.resolve(),
        args.target_center_download_dir.resolve(), args.supported_satellite_plan_dir.resolve(),
        args.supported_satellite_download_dir.resolve(), args.approved_satellite_plan_sha256)
    used_hashes: dict[str, str] = {}
    frames = []
    n_satellite_roles = 0
    for cluster_id, frame_filter in sorted(bank.supported_frames):
        indices = {
            int(target_index)
            for ident, filt, role, target_index in bank._role_coord
            if (ident, filt, role) == (cluster_id, frame_filter, "satellite_injection_center")
            and target_index is not None
        }
        if indices != set(range(len(indices))):
            raise ValueError(f"non-contiguous satellite PSF roles for {cluster_id}/{frame_filter}")
        center, satellites, names = bank.frame(cluster_id, frame_filter, len(indices))
        if len(satellites) != len(indices) or center.ndim != 2:
            raise ValueError(f"incomplete PSF frame bank for {cluster_id}/{frame_filter}")
        role_hashes = bank.file_hashes(names)
        for path, digest in role_hashes.items():
            previous = used_hashes.setdefault(path, digest)
            if previous != digest:
                raise ValueError(f"conflicting hash for PSF file {path}")
        n_satellite_roles += len(indices)
        frames.append({
            "cluster_id": cluster_id,
            "filter": frame_filter,
            "n_satellite_positions": len(indices),
            "center_psf": names[0],
            "psf_names": names,
        })
    if len(frames) != 224 or n_satellite_roles != 2693:
        raise ValueError(f"expected 224 frames and 2,693 satellite roles; got {len(frames)} and {n_satellite_roles}")
    if len(used_hashes) != 2917:
        raise ValueError(f"expected 2,917 unique required PSF coordinates, got {len(used_hashes)}")

    output.mkdir(parents=True, exist_ok=False)
    result = {
        "schema": "astrafek-supported-pairs-position-psf-bank-audit-v1",
        "status": "verified",
        "science_measurement": False,
        "psf_approved": False,
        "request_approval_verified": bank.request_approval_verified,
        "approved_satellite_plan_sha256": bank.approved_satellite_plan_sha256,
        "extended_psf_wings_validated": False,
        "images_read": False,
        "n_supported_frames": len(frames),
        "n_satellite_position_roles": n_satellite_roles,
        "n_unique_psf_files_used": len(used_hashes),
        "input_sha256": {
            str(path): sha256(path) for path in bank.provenance_paths
        },
        "test_center_psf_sha256": sha256(args.test_center_psf),
        "psf_file_sha256": used_hashes,
        "frames": frames,
        "interpretation_limit": "Verified core-PSF role coverage only; not extended-wing, sky-transfer, completeness, or science validation.",
    }
    summary = output / "summary.json"
    summary.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({
        "status": result["status"],
        "n_supported_frames": result["n_supported_frames"],
        "n_satellite_position_roles": result["n_satellite_position_roles"],
        "n_unique_psf_files_used": result["n_unique_psf_files_used"],
        "summary_sha256": sha256(summary),
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
