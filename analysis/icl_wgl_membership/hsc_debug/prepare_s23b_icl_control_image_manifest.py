#!/usr/bin/env python3
"""Prepare a checksum-pinned control-cutout manifest from the local audit.

This offline stage consumes only the completed random-position metadata audit
and the already approved target manifest. It does not retrieve images or make
any network request. Its output supports a concrete control-image request for
user approval after the target/control match yield and storage estimate are
known.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

from download_pilot import MAX_FILE_BYTES, MAX_SAMPLE_FRAMES, estimate_bytes
from fetch_ancillary_catalogues import catalogue_reader
from prepare_s23b_icl_target_image_manifest import sha256, write_json


DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
DEFAULT_AUDIT = DATA / "icl_s23b_control_candidate_audit_2026-11-07_01"
DEFAULT_TARGET_PLAN = DATA / "s23b_icl_target_image_plan_2026-11-07_01"
DEFAULT_OUTPUT = DATA / "s23b_icl_control_image_plan_2026-11-07_01"
USER_APPROVAL_REQUIRED = True


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as stream:
        return list(catalogue_reader(stream))


def prepare(audit_dir: Path, target_plan_dir: Path, output_dir: Path) -> dict:
    audit_dir, target_plan_dir = Path(audit_dir).resolve(), Path(target_plan_dir).resolve()
    output_dir = Path(output_dir).resolve()
    if output_dir.exists():
        raise FileExistsError(f"refusing to overwrite {output_dir}")
    audit_summary_path = audit_dir / "summary.json"
    candidate_audit_path = audit_dir / "candidate_audit.csv"
    controls_path = audit_dir / "provisional_controls.csv"
    target_manifest_path = target_plan_dir / "target_manifest.json"
    target_summary_path = target_plan_dir / "summary.json"
    for path in (audit_summary_path, candidate_audit_path, controls_path,
                 target_manifest_path, target_summary_path):
        if not path.is_file():
            raise FileNotFoundError(path)
    audit = json.loads(audit_summary_path.read_text(encoding="utf-8"))
    target_summary = json.loads(target_summary_path.read_text(encoding="utf-8"))
    if (audit.get("schema") != "astrafek-s23b-icl-random-control-candidate-audit-v1"
            or audit.get("status") != "complete_provisional_candidate_audit"
            or audit.get("science_measurement") is not False
            or audit.get("images_downloaded") != 0
            or target_summary.get("schema") != "astrafek-s23b-icl-target-image-plan-v1"
            or target_summary.get("status") != "prepared_after_user_image_approval"
            or target_summary.get("science_measurement") is not False):
        raise ValueError("candidate audit or target-image plan provenance is not complete")
    for path, key in ((candidate_audit_path, "candidate_audit.csv"),
                      (controls_path, "provisional_controls.csv")):
        if sha256(path) != audit.get("outputs_sha256", {}).get(key):
            raise ValueError(f"candidate audit output checksum mismatch: {path}")
    if sha256(target_manifest_path) != target_summary.get("outputs_sha256", {}).get("target_manifest.json"):
        raise ValueError("approved target manifest changed")

    target_doc = json.loads(target_manifest_path.read_text(encoding="utf-8"))
    targets = {str(row["id"]): row for row in target_doc.get("clusters", [])}
    if len(targets) != 1043:
        raise ValueError("expected all 1,043 user-approved target centers")
    matched_rows = [row for row in read_csv(controls_path)
                    if float(row["cluster_exclusion_mpc"]) == 1.0]
    target_ids, random_ids = set(), set()
    frames = []
    for row in matched_rows:
        target_id = row["target_cluster_id"].replace(":", "_")
        target = targets.get(target_id)
        random_id = str(row["random_object_id"])
        if target is None or target_id in target_ids or random_id in random_ids:
            raise ValueError("unknown target or duplicate matched control in 1-Mpc selection")
        target_ids.add(target_id); random_ids.add(random_id)
        control_id = f"CTRL_{target_id}"
        frame = {
            "id": control_id, "target_cluster_id": target_id,
            "random_object_id": random_id, "sample_role": "matched_control",
            "frame_role": "matched_control", "center_definition": "HSC PDR3 random-table sky position",
            "ra_deg": float(row["ra_deg"]), "dec_deg": float(row["dec_deg"]),
            "z": float(target["z"]), "z_cl": float(target["z_cl"]),
            "z_spec_bcg": None, "richness": float(target["richness"]),
            "requested_half_width_arcsec": float(target["requested_half_width_arcsec"]),
            "requested_half_width_kpc_proper": 1000.0,
            "target_depth_mag": float(row["target_depth_mag"]),
            "candidate_depth_mag": float(row["candidate_depth_mag"]),
            "delta_depth_mag": float(row["delta_depth_mag"]),
            "target_seeing_arcsec": float(row["target_seeing_arcsec"]),
            "candidate_seeing_arcsec": float(row["candidate_seeing_arcsec"]),
            "seeing_ratio": float(row["seeing_ratio"]),
            "target_ebv_sfd": float(row["target_ebv_sfd"]),
            "candidate_ebv_sfd": float(row["candidate_ebv_sfd"]),
            "delta_ebv_sfd": float(row["delta_ebv_sfd"]),
            "match_score": float(row["match_score"]),
            "offset_from_target_arcmin": float(row["offset_from_target_arcmin"]),
            "cluster_exclusion_mpc": 1.0,
            "skymap_id": row["skymap_id"], "tract": row["tract"], "patch": row["patch"],
            "product": "coadd/bg", "band": "HSC-I", "bcg_policy": "none",
            "members_radec_deg": [], "randoms": [], "image": None,
            "science_measurement": False,
        }
        frames.append(frame)
    if not frames:
        raise ValueError("no 1-Mpc exclusion matched controls were selected")
    if len(frames) > MAX_SAMPLE_FRAMES:
        raise ValueError("control count exceeds downloader capacity")
    estimate, estimate_rows = estimate_bytes(frames)
    if any(row["estimated_bytes"] > MAX_FILE_BYTES for row in estimate_rows):
        raise ValueError("one or more per-frame estimates exceed 512 MiB")

    manifest = {key: value for key, value in target_doc.items() if key != "clusters"}
    manifest.update({
        "product": "matched-control diffuse background",
        "center_definition": "HSC PDR3 random-table sky position matched to one S23B target",
        "selection": "one provisional HSC PDR3 random position per S23B CAMIRA target, matched in I-band depth/seeing/SFD E(B-V), excluding CAMIRA centers within 1 proper Mpc",
        "selection_limits": "provisional metadata/geometry match; full image WCS, masks, sky, variance, and PSF remain unverified",
        "controls_ready_for_science": False,
        "controls_downloaded": 0,
        "science_measurement": False,
        "clusters": frames,
    })
    output_dir.mkdir(parents=True, exist_ok=False)
    manifest_path = output_dir / "control_manifest.json"
    write_json(manifest_path, manifest)
    estimate_path = output_dir / "per_frame_estimate.tsv"
    with estimate_path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(estimate_rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(estimate_rows)
    summary = {
        "schema": "astrafek-s23b-icl-control-image-plan-v1",
        "status": "prepared_pending_specific_control_image_approval",
        "user_approval_required_before_download": USER_APPROVAL_REQUIRED,
        "n_controls": len(frames), "n_targets_with_control": len(target_ids),
        "maximum_estimated_bytes": estimate, "maximum_estimated_gib": estimate / 1024**3,
        "minimum_free_space_with_30pct_margin_bytes": math.ceil(estimate * 1.3),
        "filter_family": "HSC-I", "product": "coadd/bg",
        "planes": ["IMAGE", "MASK", "VARIANCE"],
        "half_width_kpc_proper": 1000,
        "shape_rows": False, "photoz_rows": False,
        "science_measurement": False, "images_downloaded": False,
        "inputs_sha256": {str(path): sha256(path) for path in (
            audit_summary_path, candidate_audit_path, controls_path,
            target_manifest_path, target_summary_path)},
        "outputs_sha256": {"control_manifest.json": sha256(manifest_path),
                           "per_frame_estimate.tsv": sha256(estimate_path)},
        "limitations": [
            "selection uses metadata proxies and CAMIRA exclusion, not full cutout validation",
            "CAMIRA cluster catalogs are incomplete tracers of projected contamination",
            "matched controls are provisional until the image-level QA and stack pipeline are run",
        ],
    }
    write_json(output_dir / "summary.json", summary)
    print(json.dumps({key: summary[key] for key in (
        "status", "n_controls", "n_targets_with_control", "maximum_estimated_bytes",
        "maximum_estimated_gib", "minimum_free_space_with_30pct_margin_bytes",
        "outputs_sha256", "science_measurement")}, indent=2))
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--audit-dir", type=Path, default=DEFAULT_AUDIT)
    parser.add_argument("--target-plan-dir", type=Path, default=DEFAULT_TARGET_PLAN)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    prepare(args.audit_dir, args.target_plan_dir, args.output_dir)


if __name__ == "__main__":
    main()
