#!/usr/bin/env python3
"""Prepare a hash-pinned matched-control subset for the S23B PSF pilot.

This is a local-only manifest builder. It pairs the 49 target frames that
completed the S23B morphology pilot with their already downloaded HSC PDR3
matched-control images. It does not query Gaia/CAS or download data.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path


DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
TARGET_PILOT = DATA / "s23b_detection_psf_pilot_2026-11-07_01"
CONTROL_PLAN = DATA / "s23b_icl_control_image_plan_2026-11-07_02"
CONTROL_AUDIT = DATA / "s23b_icl_control_image_audit_2026-11-07_01"
CONTROL_DOWNLOAD = (
    DATA / "s23b_icl_control_image_plan_2026-11-07_01/control_download_01/images"
)
DEFAULT_OUTPUT = DATA / "s23b_psf_control_candidates_2026-11-07_01"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, payload: dict) -> None:
    with path.open("x", encoding="utf-8") as stream:
        json.dump(payload, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write("\n")


def prepare(output: Path) -> dict:
    if output.exists():
        raise FileExistsError(f"refusing to overwrite control pilot directory: {output}")

    target_plan_path = TARGET_PILOT / "plan.json"
    target_manifest_path = TARGET_PILOT / "manifest.json"
    target_summary_path = TARGET_PILOT / "preparation/summary.json"
    control_manifest_path = CONTROL_PLAN / "control_manifest.json"
    control_plan_summary_path = CONTROL_PLAN / "summary.json"
    control_audit_summary_path = CONTROL_AUDIT / "summary.json"
    control_frames_path = CONTROL_AUDIT / "verified_control_frames.tsv"
    for path in (
        target_plan_path, target_manifest_path, target_summary_path,
        control_manifest_path, control_plan_summary_path,
        control_audit_summary_path, control_frames_path,
    ):
        if not path.is_file():
            raise FileNotFoundError(path)

    target_plan = json.loads(target_plan_path.read_text(encoding="utf-8"))
    target_manifest = json.loads(target_manifest_path.read_text(encoding="utf-8"))
    target_summary = json.loads(target_summary_path.read_text(encoding="utf-8"))
    control_manifest = json.loads(control_manifest_path.read_text(encoding="utf-8"))
    control_plan_summary = json.loads(control_plan_summary_path.read_text(encoding="utf-8"))
    control_audit_summary = json.loads(control_audit_summary_path.read_text(encoding="utf-8"))

    if target_summary.get("n_prepared") != 49 or target_summary.get("science_measurement") is not False:
        raise ValueError("target pilot must contain the 49 prepared, diagnostic-only frames")
    if control_plan_summary.get("approved_manifest_sha256") != sha256(control_manifest_path):
        raise ValueError("control manifest SHA-256 differs from the approved plan summary")
    if control_audit_summary.get("manifest_sha256") != sha256(control_manifest_path):
        raise ValueError("control audit did not pin this control manifest")
    if control_audit_summary.get("n_verified_frames") != 1041:
        raise ValueError("expected 1,041 integrity-audited matched-control images")
    if sha256(control_frames_path) != control_audit_summary["output_sha256"]["verified_control_frames.tsv"]:
        raise ValueError("verified-control frame table hash mismatch")

    target_selection = {row["id"]: row for row in target_plan["selection"]}
    target_frames = {row["id"]: row for row in target_manifest["clusters"]}
    target_prepared = [row["id"] for row in target_summary["records"]]
    if len(target_prepared) != 49 or len(set(target_prepared)) != 49:
        raise ValueError("target pilot prepared-frame list is incomplete or duplicated")
    if set(target_prepared) - set(target_selection) or set(target_prepared) - set(target_frames):
        raise ValueError("target preparation references an unselected frame")

    controls_by_target: dict[str, list[dict]] = {}
    for row in control_manifest["clusters"]:
        controls_by_target.setdefault(row["target_cluster_id"], []).append(row)
    verified_by_id = {row["id"]: row for row in read_tsv(control_frames_path)}

    selected_controls = []
    selection_rows = []
    image_hashes = {}
    for target_id in target_prepared:
        matches = controls_by_target.get(target_id, [])
        if len(matches) != 1:
            raise ValueError(f"expected one matched control for {target_id}, found {len(matches)}")
        control = dict(matches[0])
        verified = verified_by_id.get(control["id"])
        if verified is None:
            raise ValueError(f"control frame is absent from the integrity audit: {control['id']}")
        target = target_selection[target_id]
        target_filter = target["primary_filter"]
        if verified["filter"] != target_filter:
            raise ValueError(
                f"filter mismatch for target/control pair {target_id}: "
                f"{target_filter} != {verified['filter']}"
            )
        if control.get("members_radec_deg"):
            raise ValueError(f"matched control unexpectedly carries member coordinates: {control['id']}")
        image_path = CONTROL_DOWNLOAD / verified["file"]
        if not image_path.is_file():
            raise FileNotFoundError(image_path)
        if image_path.stat().st_size != int(verified["bytes"]):
            raise ValueError(f"control image byte count changed: {image_path}")
        image_sha = sha256(image_path)
        if image_sha != verified["sha256"]:
            raise ValueError(f"control image checksum differs from the integrity audit: {image_path}")

        control["image"] = str(image_path.resolve())
        control["preparation_cohort"] = "matched_control_psf_pilot"
        control["paired_target_cluster_id"] = target_id
        control["primary_filter"] = target_filter
        selected_controls.append(control)
        image_hashes[str(image_path.resolve())] = image_sha
        selection_rows.append({
            "target_cluster_id": target_id,
            "target_primary_filter": target_filter,
            "control_id": control["id"],
            "control_primary_filter": verified["filter"],
            "control_image": str(image_path.resolve()),
            "image_bytes": verified["bytes"],
            "image_sha256": image_sha,
            "match_score": control["match_score"],
            "target_seeing_arcsec": control["target_seeing_arcsec"],
            "control_seeing_arcsec": control["candidate_seeing_arcsec"],
            "center_definition": control["center_definition"],
        })

    if len({row["id"] for row in selected_controls}) != len(selected_controls):
        raise ValueError("one matched control was reused for multiple target frames")
    if len(selected_controls) != 49:
        raise ValueError("control pilot does not contain exactly 49 frames")

    output.mkdir(parents=True, exist_ok=False)
    subset_manifest = dict(control_manifest)
    subset_manifest["clusters"] = selected_controls
    subset_manifest["science_measurement"] = False
    subset_manifest["controls_ready_for_science"] = False
    manifest_path = output / "pilot_manifest.json"
    write_json(manifest_path, subset_manifest)
    selection_path = output / "selection.tsv"
    write_tsv(selection_path, selection_rows, list(selection_rows[0]))

    input_paths = (
        target_plan_path, target_manifest_path, target_summary_path,
        control_manifest_path, control_plan_summary_path,
        control_audit_summary_path, control_frames_path,
    )
    summary = {
        "schema": "astrafek-s23b-matched-control-psf-pilot-plan-v1",
        "status": "prepared_local_only",
        "n_target_frames": len(target_prepared),
        "n_control_frames": len(selected_controls),
        "filter_counts": dict(Counter(row["target_primary_filter"] for row in selection_rows)),
        "n_image_hashes_verified": len(image_hashes),
        "input_sha256": {str(path.resolve()): sha256(path) for path in input_paths},
        "verified_image_sha256": image_hashes,
        "pilot_manifest_sha256": sha256(manifest_path),
        "selection_sha256": sha256(selection_path),
        "science_measurement": False,
        "images_downloaded": False,
        "external_queries_submitted": False,
        "notes": [
            "Pairs are the exact matched controls for the 49 target frames successfully prepared in the target pilot.",
            "All control images already existed and matched the prior 1,041-frame byte/hash/WCS audit.",
            "The control center is a random sky position, not a BCG or member-galaxy position.",
            "Any central-position diagnostic from the legacy candidate runner is not a BCG identification.",
            "This manifest only prepares local morphology candidates; it submits no Gaia query.",
        ],
    }
    write_json(output / "summary.json", summary)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    summary = prepare(args.output.resolve())
    print(json.dumps({key: summary[key] for key in (
        "status", "n_target_frames", "n_control_frames", "filter_counts",
        "n_image_hashes_verified", "pilot_manifest_sha256",
        "selection_sha256", "science_measurement", "external_queries_submitted",
    )}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
