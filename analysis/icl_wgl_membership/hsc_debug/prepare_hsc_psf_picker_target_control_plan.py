#!/usr/bin/env python3
"""Prepare an offline PSF Picker plan for target/control frame centers.

The planner only reads local manifests and previously received PSF files. It
does not contact HSC, submit a request, or download data.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import statistics
from pathlib import Path


DEFAULT_ROOT = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
DEFAULT_TARGET_MANIFEST = DEFAULT_ROOT / "sample_manifest.json"
DEFAULT_CONTROL_MANIFEST = DEFAULT_ROOT / "control_download_01/download_manifest.json"
DEFAULT_APPROVED_MAPPING = DEFAULT_ROOT / "icl_psf_picker_injection_plan_2026-10-19_02/position_mapping.json"
DEFAULT_PRIOR_AUDIT = DEFAULT_ROOT / "icl_psf_picker_injection_2026-10-19_01/psf_picker_audit.json"
DEFAULT_PRIOR_PSF_DIR = DEFAULT_ROOT / "icl_psf_picker_injection_2026-10-19_01"
DEFAULT_OUTPUT = DEFAULT_ROOT / "icl_psf_picker_target_control_centers_plan_2026-10-21_01"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _load_rows(path: Path, label: str) -> list[dict]:
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    rows = data.get("clusters")
    if not isinstance(rows, list):
        raise ValueError(f"{label} manifest must contain a clusters list")
    ids = [str(row.get("id", row.get("target_cluster_id", ""))) for row in rows]
    if any(not ident for ident in ids) or len(set(ids)) != len(ids):
        raise ValueError(f"{label} rows need nonempty, unique IDs")
    return rows


def _position_key(row: dict) -> tuple[float, float]:
    ra, dec = float(row["ra_deg"]), float(row["dec_deg"])
    if not (math.isfinite(ra) and 0 <= ra < 360 and math.isfinite(dec) and -90 <= dec <= 90):
        raise ValueError("position coordinates must be finite celestial coordinates")
    return round(ra, 9), round(dec, 9)


def collect_positions(targets: list[dict], controls: list[dict], approved: list[dict]):
    positions: dict[tuple[float, float], dict] = {}
    for row in targets:
        key = _position_key(row)
        item = positions.setdefault(key, {"ra_deg": key[0], "dec_deg": key[1], "roles": []})
        item["roles"].append({"role": "target_cluster_center", "cluster_id": str(row["id"])})
    for row in controls:
        key = _position_key(row)
        item = positions.setdefault(key, {"ra_deg": key[0], "dec_deg": key[1], "roles": []})
        item["roles"].append({
            "role": "matched_control_center",
            "control_id": str(row.get("id", row.get("random_object_id", ""))),
            "target_cluster_id": str(row["target_cluster_id"]),
        })
    covered = {_position_key(row) for row in approved}
    all_rows = sorted(positions.values(), key=lambda row: (row["ra_deg"], row["dec_deg"]))
    remaining = [row for row in all_rows if (row["ra_deg"], row["dec_deg"]) not in covered]
    return all_rows, remaining, covered


def prepare(target_manifest: Path, control_manifest: Path, approved_mapping: Path,
            prior_audit: Path, prior_psf_dir: Path, output: Path,
            *, expected_targets: int = 800, expected_controls: int = 798) -> dict:
    paths = [Path(p).resolve() for p in
             (target_manifest, control_manifest, approved_mapping, prior_audit, prior_psf_dir)]
    target_manifest, control_manifest, approved_mapping, prior_audit, prior_psf_dir = paths
    output = Path(output).resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite existing plan directory: {output}")

    targets = _load_rows(target_manifest, "target")
    controls = _load_rows(control_manifest, "control")
    if len(targets) != expected_targets or len(controls) != expected_controls:
        raise ValueError(
            f"unexpected center counts: targets={len(targets)} controls={len(controls)}"
        )
    approved = json.loads(approved_mapping.read_text(encoding="utf-8"))
    if not isinstance(approved, list):
        raise ValueError("approved PSF mapping must be a list")
    audit = json.loads(prior_audit.read_text(encoding="utf-8"))
    sizes = sorted(p.stat().st_size for p in prior_psf_dir.iterdir()
                   if p.is_file() and p.suffix.lower() == ".fits")
    if (len(sizes) != int(audit["n_returned"])
            or sum(sizes) != int(audit["total_bytes"])):
        raise ValueError("received PSF file count/bytes disagree with the saved audit")

    all_rows, remaining, covered = collect_positions(targets, controls, approved)
    output.mkdir(parents=True, exist_ok=False)
    coordinates = output / "psf_picker_coordinates.txt"
    with coordinates.open("x", encoding="ascii", newline="") as stream:
        stream.write("#? ra dec\n")
        for row in remaining:
            stream.write(f"{row['ra_deg']:.9f} {row['dec_deg']:+.9f}\n")
    mapping = output / "position_mapping.json"
    mapping.write_text(json.dumps(remaining, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")

    mean_size = sum(sizes) / len(sizes)
    median_size = statistics.median(sizes)
    plan = {
        "schema": "astrafek-hsc-psf-picker-target-control-center-plan-v1",
        "status": "prepared_not_submitted",
        "purpose": "position-dependent PSF core inventory for existing target/matched-control centers",
        "n_target_centers": len(targets),
        "n_control_centers": len(controls),
        "n_unique_combined_centers": len(all_rows),
        "n_already_covered_by_approved_psfs": len(
            {(row["ra_deg"], row["dec_deg"]) for row in all_rows} & covered
        ),
        "n_requests_remaining": len(remaining),
        "rerun": "pdr3_wide",
        "type": "coadd",
        "filter": "HSC-I",
        "centered": True,
        "coordinate_list": str(coordinates),
        "coordinate_list_sha256": sha256(coordinates),
        "position_mapping": str(mapping),
        "position_mapping_sha256": sha256(mapping),
        "target_manifest": str(target_manifest),
        "target_manifest_sha256": sha256(target_manifest),
        "control_manifest": str(control_manifest),
        "control_manifest_sha256": sha256(control_manifest),
        "previous_approved_mapping": str(approved_mapping),
        "previous_approved_mapping_sha256": sha256(approved_mapping),
        "previous_psf_audit": str(prior_audit),
        "previous_psf_audit_sha256": sha256(prior_audit),
        "previous_psf_bytes_min": sizes[0],
        "previous_psf_bytes_median": median_size,
        "previous_psf_bytes_mean": mean_size,
        "previous_psf_bytes_max": sizes[-1],
        "estimated_total_bytes_median": len(remaining) * median_size,
        "estimated_total_bytes_mean": round(len(remaining) * mean_size),
        "estimated_total_bytes_max": len(remaining) * sizes[-1],
        "images_requested": False,
        "shape_rows_requested": False,
        "photoz_requested": False,
        "credentials_saved": False,
        "network_used": False,
        "request_submitted": False,
        "explicit_request_approval_required": True,
        "science_measurement": False,
    }
    plan_path = output / "psf_picker_plan.json"
    plan_path.write_text(json.dumps(plan, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return plan


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--target-manifest", type=Path, default=DEFAULT_TARGET_MANIFEST)
    parser.add_argument("--control-manifest", type=Path, default=DEFAULT_CONTROL_MANIFEST)
    parser.add_argument("--approved-mapping", type=Path, default=DEFAULT_APPROVED_MAPPING)
    parser.add_argument("--prior-audit", type=Path, default=DEFAULT_PRIOR_AUDIT)
    parser.add_argument("--prior-psf-dir", type=Path, default=DEFAULT_PRIOR_PSF_DIR)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--expected-targets", type=int, default=800)
    parser.add_argument("--expected-controls", type=int, default=798)
    args = parser.parse_args()
    plan = prepare(args.target_manifest, args.control_manifest, args.approved_mapping,
                   args.prior_audit, args.prior_psf_dir, args.output,
                   expected_targets=args.expected_targets,
                   expected_controls=args.expected_controls)
    print(json.dumps({key: plan[key] for key in (
        "status", "n_target_centers", "n_control_centers", "n_unique_combined_centers",
        "n_already_covered_by_approved_psfs", "n_requests_remaining", "filter",
        "estimated_total_bytes_mean", "network_used", "request_submitted",
    )}, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
