#!/usr/bin/env python3
"""Verify every approved PSF role resolves through the Astrafek PSF bank."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from position_psf import PositionPSFBank, sha256


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan-dir", type=Path, required=True)
    parser.add_argument("--download-dir", type=Path, required=True)
    parser.add_argument("--test-center-psf", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    output_dir = args.output_dir.resolve()
    if output_dir.exists():
        raise FileExistsError(f"refusing to overwrite output directory: {output_dir}")

    bank = PositionPSFBank(args.plan_dir, args.download_dir,
                           args.test_center_psf)
    mapping_path = args.plan_dir / "position_mapping.json"
    mapping = json.loads(mapping_path.read_text(encoding="utf-8"))
    satellite_indices: dict[tuple[str, str], set[int]] = {}
    n_satellite_roles = 0
    for position in mapping:
        for role in position["roles"]:
            if role["role"] != "satellite_injection_center":
                continue
            key = (str(role["cluster_id"]), str(role["frame_filter"]).lower())
            satellite_indices.setdefault(key, set()).add(int(role["target_index"]))
            n_satellite_roles += 1

    used_hashes: dict[str, str] = {}
    verified_frames = []
    for cluster_id, frame_filter in sorted(bank.supported_frames):
        indices = satellite_indices.get((cluster_id, frame_filter), set())
        if indices and indices != set(range(max(indices) + 1)):
            raise ValueError(f"satellite target indices are not contiguous: {cluster_id}/{frame_filter}")
        n_satellites = len(indices)
        center, satellites, names = bank.frame(cluster_id, frame_filter,
                                               n_satellites)
        if len(center.shape) != 2 or len(satellites) != n_satellites:
            raise AssertionError(f"PSF bank returned an incomplete frame: {cluster_id}/{frame_filter}")
        hashes = bank.file_hashes(names)
        for path, digest in hashes.items():
            previous = used_hashes.setdefault(path, digest)
            if previous != digest:
                raise ValueError(f"inconsistent PSF checksum for {path}")
        verified_frames.append({
            "cluster_id": cluster_id,
            "frame_filter": frame_filter,
            "satellite_psf_count": n_satellites,
            "center_psf": names[0],
        })

    if len(bank.supported_frames) != 57 or n_satellite_roles != 688:
        raise ValueError("role coverage differs from the approved 57-frame/688-position sample")
    if len(used_hashes) != 745:
        raise ValueError(f"expected 745 unique core PSFs in use, got {len(used_hashes)}")

    script_path = Path(__file__).resolve()
    result = {
        "schema": "astrafek-position-psf-bank-validation-v1",
        "status": "verified",
        "plan": str((args.plan_dir / "psf_picker_plan.json").resolve()),
        "plan_sha256": sha256(args.plan_dir / "psf_picker_plan.json"),
        "mapping_sha256": sha256(mapping_path),
        "audit_manifest": str((args.download_dir / "psf_picker_audit.json").resolve()),
        "audit_manifest_sha256": sha256(args.download_dir / "psf_picker_audit.json"),
        "test_center_psf": str(args.test_center_psf.resolve()),
        "test_center_psf_sha256": sha256(args.test_center_psf),
        "n_supported_frames": len(verified_frames),
        "n_satellite_positions": n_satellite_roles,
        "n_unique_psfs_loaded_and_hash_verified": len(used_hashes),
        "science_measurement": False,
        "extended_psf_wings_validated": False,
        "script": str(script_path),
        "script_sha256": sha256(script_path),
        "psf_files": [{"path": path, "sha256": digest}
                       for path, digest in sorted(used_hashes.items())],
        "frames": verified_frames,
    }
    output_dir.mkdir(parents=True, exist_ok=False)
    summary_path = output_dir / "summary.json"
    summary_path.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n",
                            encoding="utf-8")
    summary_hash = hashlib.sha256(summary_path.read_bytes()).hexdigest()
    print(json.dumps({
        "status": result["status"],
        "n_supported_frames": result["n_supported_frames"],
        "n_satellite_positions": result["n_satellite_positions"],
        "n_unique_psfs_loaded_and_hash_verified": result["n_unique_psfs_loaded_and_hash_verified"],
        "science_measurement": result["science_measurement"],
        "extended_psf_wings_validated": result["extended_psf_wings_validated"],
        "summary_sha256": summary_hash,
    }, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
