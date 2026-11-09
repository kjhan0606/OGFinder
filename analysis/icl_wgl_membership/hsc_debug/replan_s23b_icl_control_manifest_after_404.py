#!/usr/bin/env python3
"""Create a revised, offline control manifest after one approved cutout returns 404.

The script keeps the approved 1,041-target scope and all unchanged control
positions, excludes the failed random-object ID, and selects the best remaining
metadata-matched, CAMIRA-clean position for that same target. It requires new
user approval for the revised manifest hash and performs no network requests.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

from download_pilot import estimate_bytes
from prepare_s23b_icl_target_image_manifest import sha256, write_json


DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
DEFAULT_PLAN = DATA / "s23b_icl_control_image_plan_2026-11-07_01"
DEFAULT_AUDIT = DATA / "icl_s23b_control_candidate_audit_2026-11-07_01"
DEFAULT_OUTPUT = DATA / "s23b_icl_control_image_plan_2026-11-07_02"
EXPECTED_ORIGINAL_MANIFEST_SHA256 = (
    "4e182a3160055f1dfec4171035bceb1f4a8b574d10f0b12660f67a39dd037c35"
)


def angular_separation_deg(ra1: float, dec1: float, ra2: float, dec2: float) -> float:
    ra1, dec1, ra2, dec2 = map(math.radians, (ra1, dec1, ra2, dec2))
    cosine = (math.sin(dec1) * math.sin(dec2)
              + math.cos(dec1) * math.cos(dec2) * math.cos(ra1 - ra2))
    return math.degrees(math.acos(max(-1.0, min(1.0, cosine))))


def footprint_cap_deg(frame: dict) -> float:
    return math.sqrt(2.0) * float(frame["requested_half_width_arcsec"]) / 3600.0 + 1.0 / 3600.0


def read_candidates(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream))


def prepare(plan_dir: Path, audit_dir: Path, failed_control_id: str,
            failed_random_object_id: str, output_dir: Path) -> dict:
    plan_dir, audit_dir = plan_dir.resolve(), audit_dir.resolve()
    output_dir = output_dir.resolve()
    if output_dir.exists():
        raise FileExistsError(f"refusing to overwrite {output_dir}")
    manifest_path = plan_dir / "control_manifest.json"
    plan_summary_path = plan_dir / "summary.json"
    candidate_path = audit_dir / "candidate_audit.csv"
    audit_summary_path = audit_dir / "summary.json"
    for path in (manifest_path, plan_summary_path, candidate_path, audit_summary_path):
        if not path.is_file():
            raise FileNotFoundError(path)

    original_sha = sha256(manifest_path)
    if original_sha != EXPECTED_ORIGINAL_MANIFEST_SHA256:
        raise ValueError("original approved control manifest SHA-256 changed")
    old_summary = json.loads(plan_summary_path.read_text(encoding="utf-8"))
    audit_summary = json.loads(audit_summary_path.read_text(encoding="utf-8"))
    if (old_summary.get("status") != "prepared_pending_specific_control_image_approval"
            or old_summary.get("n_controls") != 1041
            or audit_summary.get("status") != "complete_provisional_candidate_audit"
            or audit_summary.get("science_measurement") is not False
            or sha256(candidate_path) != audit_summary.get("outputs_sha256", {}).get("candidate_audit.csv")):
        raise ValueError("control-plan or candidate-audit provenance failed")

    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    frames = manifest.get("clusters") or []
    if len(frames) != 1041 or any(frame.get("image") for frame in frames):
        raise ValueError("expected 1,041 undownloaded control frames")
    failed = [frame for frame in frames if str(frame.get("id")) == failed_control_id]
    if len(failed) != 1 or str(failed[0].get("random_object_id")) != failed_random_object_id:
        raise ValueError("failed control ID/random-object pair does not match the approved manifest")
    failed_frame = failed[0]
    target_id = str(failed_frame["target_cluster_id"])
    used_random_ids = {str(frame["random_object_id"]) for frame in frames}
    used_random_ids.discard(failed_random_object_id)
    existing = [frame for frame in frames if frame is not failed_frame]

    eligible = []
    for row in read_candidates(candidate_path):
        if row.get("target_cluster_id", "").replace(":", "_") != target_id:
            continue
        if (row.get("metadata_match") != "True"
                or int(row.get("cluster_overlap_1mpc", -1)) != 0
                or str(row.get("random_object_id")) in used_random_ids
                or str(row.get("random_object_id")) == failed_random_object_id):
            continue
        candidate_cap = float(row["footprint_cap_deg"])
        if any(angular_separation_deg(
                    float(row["ra_deg"]), float(row["dec_deg"]),
                    float(frame["ra_deg"]), float(frame["dec_deg"]))
               <= candidate_cap + footprint_cap_deg(frame) for frame in existing):
            continue
        eligible.append(row)
    eligible.sort(key=lambda row: (float(row["match_score"]),
                                   float(row["offset_from_target_arcmin"]),
                                   str(row["random_object_id"])))
    if not eligible:
        raise ValueError(f"no globally non-overlapping replacement candidate for {target_id}")
    replacement = eligible[0]

    # Change only the failed position's metadata; preserve target ID and output filename.
    updated = dict(failed_frame)
    updates = {
        "random_object_id": str(replacement["random_object_id"]),
        "ra_deg": float(replacement["ra_deg"]),
        "dec_deg": float(replacement["dec_deg"]),
        "candidate_depth_mag": float(replacement["candidate_depth_mag"]),
        "delta_depth_mag": float(replacement["delta_depth_mag"]),
        "candidate_seeing_arcsec": float(replacement["candidate_seeing_arcsec"]),
        "seeing_ratio": float(replacement["seeing_ratio"]),
        "candidate_ebv_sfd": float(replacement["candidate_ebv_sfd"]),
        "delta_ebv_sfd": float(replacement["delta_ebv_sfd"]),
        "match_score": float(replacement["match_score"]),
        "offset_from_target_arcmin": float(replacement["offset_from_target_arcmin"]),
        "skymap_id": replacement["skymap_id"],
        "tract": replacement["tract"],
        "patch": replacement["patch"],
    }
    updated.update(updates)
    frames = [updated if frame is failed_frame else frame for frame in frames]
    manifest["clusters"] = frames
    manifest["selection_limits"] = (
        str(manifest.get("selection_limits", ""))
        + "; control for " + target_id + " replaced after HTTP 404, using the next-ranked "
        + "globally non-overlapping 1-Mpc CAMIRA-clean metadata match"
    ).strip("; ")

    estimate, estimate_rows = estimate_bytes(frames)
    if estimate != int(old_summary["maximum_estimated_bytes"]):
        raise ValueError("replacement unexpectedly changed the approved cutout-size estimate")
    output_dir.mkdir(parents=True, exist_ok=False)
    revised_manifest = output_dir / "control_manifest.json"
    write_json(revised_manifest, manifest)
    estimate_path = output_dir / "per_frame_estimate.tsv"
    with estimate_path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(estimate_rows[0]), delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(estimate_rows)
    summary = {
        "schema": "astrafek-s23b-icl-control-image-plan-v1",
        "status": "prepared_pending_reapproval_after_http_404",
        "user_approval_required_before_resume": True,
        "n_controls": len(frames),
        "n_targets_with_control": len({str(frame["target_cluster_id"]) for frame in frames}),
        "failed_control_id": failed_control_id,
        "failed_random_object_id": failed_random_object_id,
        "failed_http_status": 404,
        "replacement_target_cluster_id": target_id,
        "replacement_random_object_id": str(replacement["random_object_id"]),
        "replacement_ra_deg": float(replacement["ra_deg"]),
        "replacement_dec_deg": float(replacement["dec_deg"]),
        "replacement_match_score": float(replacement["match_score"]),
        "replacement_tract": replacement["tract"],
        "replacement_patch": replacement["patch"],
        "original_approved_manifest_sha256": original_sha,
        "maximum_estimated_bytes": estimate,
        "maximum_estimated_gib": estimate / 1024**3,
        "minimum_free_space_with_30pct_margin_bytes": math.ceil(estimate * 1.3),
        "filter_family": "HSC-I", "product": "coadd/bg",
        "planes": ["IMAGE", "MASK", "VARIANCE"],
        "half_width_kpc_proper": 1000,
        "shape_rows": False, "photoz_rows": False,
        "science_measurement": False, "images_downloaded": False,
        "inputs_sha256": {str(path): sha256(path) for path in (
            manifest_path, plan_summary_path, candidate_path, audit_summary_path)},
        "outputs_sha256": {"control_manifest.json": sha256(revised_manifest),
                           "per_frame_estimate.tsv": sha256(estimate_path)},
        "limitations": [
            "all 1,040 unchanged controls retain their original positions",
            "replacement passed metadata and global-footprint non-overlap screens, but its image is unverified",
            "the revised manifest requires explicit user approval before resume",
        ],
    }
    write_json(output_dir / "summary.json", summary)
    print(json.dumps({key: summary[key] for key in (
        "status", "n_controls", "replacement_target_cluster_id",
        "failed_random_object_id", "replacement_random_object_id",
        "replacement_ra_deg", "replacement_dec_deg", "replacement_match_score",
        "original_approved_manifest_sha256", "maximum_estimated_bytes",
        "minimum_free_space_with_30pct_margin_bytes", "outputs_sha256")}, indent=2))
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan-dir", type=Path, default=DEFAULT_PLAN)
    parser.add_argument("--audit-dir", type=Path, default=DEFAULT_AUDIT)
    parser.add_argument("--failed-control-id", default="CTRL_S23B_456")
    parser.add_argument("--failed-random-object-id", default="42740933254447653")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    prepare(args.plan_dir, args.audit_dir, args.failed_control_id,
            args.failed_random_object_id, args.output_dir)


if __name__ == "__main__":
    main()
