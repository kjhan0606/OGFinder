#!/usr/bin/env python3
"""Prepare a stratified, offline S23B target/control PSF-candidate expansion.

The script selects already-downloaded target/control image pairs, attaches the
approved CAMIRA member positions to target frames for compact-source screening,
and writes hash-pinned manifests for the existing local CPU preparation tool.
It does not submit Gaia/CAS queries or download data.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import random
from collections import Counter, defaultdict
from pathlib import Path


DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
TARGET_PLAN = DATA / "s23b_icl_target_image_plan_2026-11-07_01"
TARGET_AUDIT = DATA / "s23b_icl_target_image_audit_2026-11-07_01"
TARGET_IMAGES = TARGET_PLAN / "target_download_01/images"
CONTROL_PLAN = DATA / "s23b_icl_control_image_plan_2026-11-07_02"
CONTROL_AUDIT = DATA / "s23b_icl_control_image_audit_2026-11-07_01"
CONTROL_IMAGES = DATA / "s23b_icl_control_image_plan_2026-11-07_01/control_download_01/images"
TARGET_PILOT = DATA / "s23b_detection_psf_pilot_2026-11-07_01"
MEMBER_IMPORT = DATA / "s23b_camira_member_catalog_2026-11-07_01/import_2026-11-07_02"
DEFAULT_OUTPUT = DATA / "s23b_psf_expansion_2026-11-07_01/plan"
PER_FILTER_PAIRS = 100
N_SEEING_STRATA = 5
SEED = 20261107


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


def load_json(path: Path) -> dict:
    with path.open(encoding="utf-8") as stream:
        return json.load(stream)


def strata_sample(rows: list[dict], *, filter_name: str, n: int,
                  seeing_key: str, seed: int) -> list[tuple[dict, int]]:
    ordered = sorted(rows, key=lambda row: (float(row[seeing_key]), row["target_cluster_id"]))
    if len(ordered) < n or n % N_SEEING_STRATA:
        raise ValueError(f"cannot draw {n} pairs across {N_SEEING_STRATA} seeing strata for {filter_name}")
    rng = random.Random(seed + (0 if filter_name == "i" else 1))
    per_stratum = n // N_SEEING_STRATA
    selected: list[tuple[dict, int]] = []
    for index in range(N_SEEING_STRATA):
        lo = math.floor(index * len(ordered) / N_SEEING_STRATA)
        hi = math.floor((index + 1) * len(ordered) / N_SEEING_STRATA)
        group = ordered[lo:hi]
        selected.extend((row, index + 1) for row in rng.sample(group, per_stratum))
    return selected


def build(output: Path) -> dict:
    output = output.resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite PSF expansion plan: {output}")

    paths = {
        "target_manifest": TARGET_PLAN / "target_manifest.json",
        "target_plan_summary": TARGET_PLAN / "summary.json",
        "target_audit_summary": TARGET_AUDIT / "summary.json",
        "target_verified_frames": TARGET_AUDIT / "verified_target_frames.tsv",
        "control_manifest": CONTROL_PLAN / "control_manifest.json",
        "control_plan_summary": CONTROL_PLAN / "summary.json",
        "control_audit_summary": CONTROL_AUDIT / "summary.json",
        "control_verified_frames": CONTROL_AUDIT / "verified_control_frames.tsv",
        "target_pilot_manifest": TARGET_PILOT / "manifest.json",
        "target_pilot_preparation": TARGET_PILOT / "preparation/summary.json",
        "member_positions": MEMBER_IMPORT / "target_member_positions.tsv",
        "member_import_summary": MEMBER_IMPORT / "summary.json",
    }
    for path in paths.values():
        if not path.is_file():
            raise FileNotFoundError(path)

    target_manifest = load_json(paths["target_manifest"])
    target_plan_summary = load_json(paths["target_plan_summary"])
    target_audit = load_json(paths["target_audit_summary"])
    control_manifest = load_json(paths["control_manifest"])
    control_plan_summary = load_json(paths["control_plan_summary"])
    control_audit = load_json(paths["control_audit_summary"])
    pilot_manifest = load_json(paths["target_pilot_manifest"])
    pilot_preparation = load_json(paths["target_pilot_preparation"])
    member_summary = load_json(paths["member_import_summary"])

    target_rows = read_tsv(paths["target_verified_frames"])
    control_rows = read_tsv(paths["control_verified_frames"])
    member_rows = read_tsv(paths["member_positions"])

    if (target_audit.get("status") != "complete_header_wcs_file_integrity_audit"
            or target_audit.get("n_verified_frames") != 1043
            or target_audit.get("output_sha256", {}).get("verified_target_frames.tsv")
            != sha256(paths["target_verified_frames"])):
        raise ValueError("S23B target image audit does not match its verified frame table")
    if (control_audit.get("status") != "complete_header_wcs_file_integrity_audit"
            or control_audit.get("n_verified_frames") != 1041
            or control_audit.get("output_sha256", {}).get("verified_control_frames.tsv")
            != sha256(paths["control_verified_frames"])):
        raise ValueError("S23B control image audit does not match its verified frame table")
    target_manifest_sha = sha256(paths["target_manifest"])
    if member_summary.get("status") != "complete_catalog_join_and_image_footprint_audit":
        raise ValueError("approved CAMIRA member import is incomplete")
    if member_summary.get("catalog_join", {}).get("target_manifest_sha256") != target_manifest_sha:
        raise ValueError("CAMIRA member positions were joined against a different target manifest")
    if (member_summary.get("output_sha256", {}).get("target_member_positions.tsv")
            != sha256(paths["member_positions"])):
        raise ValueError("CAMIRA member-position table checksum differs from its import audit")
    if (control_plan_summary.get("approved_manifest_sha256") != sha256(paths["control_manifest"])
            or control_audit.get("manifest_sha256") != sha256(paths["control_manifest"])):
        raise ValueError("control images are not pinned to the approved control manifest")
    if (pilot_preparation.get("n_prepared") != 49
            or pilot_preparation.get("science_measurement") is not False
            or pilot_preparation.get("manifest_sha256") != sha256(paths["target_pilot_manifest"])):
        raise ValueError("existing target PSF pilot does not match its preparation manifest")

    target_verified = {row["id"]: row for row in target_rows}
    control_verified = {row["id"]: row for row in control_rows}
    target_by_id = {row["id"]: row for row in target_manifest["clusters"]}
    controls_by_target: dict[str, list[dict]] = defaultdict(list)
    for row in control_manifest["clusters"]:
        controls_by_target[row["target_cluster_id"]].append(row)
    members_by_target: dict[str, list[list[float]]] = defaultdict(list)
    member_counts: Counter[str] = Counter()
    for row in member_rows:
        target_id = row["target_cluster_id"]
        if row["inside_image_pixel_centers"] == "True":
            members_by_target[target_id].append([float(row["member_ra_deg"]), float(row["member_dec_deg"])])
            member_counts[target_id] += 1

    prepared_ids = {row["id"] for row in pilot_preparation["records"]}
    if len(prepared_ids) != 49:
        raise ValueError("expected 49 distinct already-prepared target frames")

    pairs_by_filter: dict[str, list[dict]] = defaultdict(list)
    pair_ids = set()
    filter_mismatch_ids: list[dict[str, str]] = []
    for target_id, target in target_by_id.items():
        control_matches = controls_by_target.get(target_id, [])
        if len(control_matches) != 1:
            continue
        target_audit_row = target_verified.get(target_id)
        if target_audit_row is None:
            raise ValueError(f"target image audit lacks {target_id}")
        control = control_matches[0]
        control_audit_row = control_verified.get(control["id"])
        if control_audit_row is None:
            raise ValueError(f"control image audit lacks {control['id']}")
        filt = target_audit_row["primary_filter"].lower()
        control_filter = control_audit_row["filter"].lower()
        if filt not in ("i", "i2") or control_filter not in ("i", "i2"):
            raise ValueError(f"unrecognized HSC-I family filter for {target_id}")
        if control_filter != filt:
            filter_mismatch_ids.append({"target_cluster_id": target_id,
                                        "target_filter": filt,
                                        "control_filter": control_filter})
            continue
        if target_id in prepared_ids:
            continue
        seeing = float(control["target_seeing_arcsec"])
        pair = {
            "target_cluster_id": target_id,
            "target_filter": filt,
            "target_seeing_arcsec": seeing,
            "control_id": control["id"],
            "control_seeing_arcsec": float(control["candidate_seeing_arcsec"]),
            "control_match_score": float(control["match_score"]),
            "target": target,
            "control": control,
            "target_audit": target_audit_row,
            "control_audit": control_audit_row,
        }
        if target_id in pair_ids or control["id"] in pair_ids:
            raise ValueError("duplicate target/control frame in pair inventory")
        pair_ids.add(target_id)
        pair_ids.add(control["id"])
        pairs_by_filter[filt].append(pair)

    selected: list[tuple[dict, int]] = []
    for filt in ("i", "i2"):
        selected.extend(strata_sample(pairs_by_filter[filt], filter_name=filt,
                                      n=PER_FILTER_PAIRS,
                                      seeing_key="target_seeing_arcsec", seed=SEED))
    if len({pair["target_cluster_id"] for pair, _ in selected}) != 2 * PER_FILTER_PAIRS:
        raise ValueError("selected pair list contains duplicates")

    selected.sort(key=lambda item: (item[0]["target_filter"], item[0]["target_cluster_id"]))
    target_frames = []
    control_frames = []
    selection_rows = []
    for pair, stratum in selected:
        target = dict(pair["target"])
        target_audit_row = pair["target_audit"]
        target_path = TARGET_IMAGES / target_audit_row["file"]
        if not target_path.is_file() or target_path.stat().st_size != int(target_audit_row["bytes"]):
            raise ValueError(f"audited target image is missing or has changed size: {target_path}")
        target.update(
            image=str(target_path.resolve()),
            filter="HSC-I",
            coadd_header_filter=target_audit_row["primary_filter"].lower(),
            input_image_sha256=target_audit_row["sha256"],
            input_image_bytes=int(target_audit_row["bytes"]),
            members_radec_deg=members_by_target.get(target["id"], []),
            member_catalog_inside_frame_count=member_counts.get(target["id"], 0),
            preparation_cohort="s23b_matched_psf_expansion_target",
            paired_control_id=pair["control"]["id"],
            center_is_bcg=False,
            science_measurement=False,
        )

        control = dict(pair["control"])
        control_audit_row = pair["control_audit"]
        control_path = CONTROL_IMAGES / control_audit_row["file"]
        if not control_path.is_file() or control_path.stat().st_size != int(control_audit_row["bytes"]):
            raise ValueError(f"audited control image is missing or has changed size: {control_path}")
        control.update(
            image=str(control_path.resolve()),
            filter="HSC-I",
            coadd_header_filter=control_audit_row["filter"].lower(),
            input_image_sha256=control_audit_row["sha256"],
            input_image_bytes=int(control_audit_row["bytes"]),
            members_radec_deg=[],
            preparation_cohort="s23b_matched_psf_expansion_control",
            paired_target_cluster_id=pair["target_cluster_id"],
            center_is_bcg=False,
            science_measurement=False,
        )
        target_frames.append(target)
        control_frames.append(control)
        selection_rows.append({
            "target_cluster_id": pair["target_cluster_id"],
            "control_id": pair["control_id"],
            "filter": pair["target_filter"],
            "seeing_stratum": stratum,
            "stratification_seed": SEED,
            "target_seeing_arcsec": pair["target_seeing_arcsec"],
            "control_seeing_arcsec": pair["control_seeing_arcsec"],
            "control_match_score": pair["control_match_score"],
            "target_image": str(target_path.resolve()),
            "target_image_sha256": target_audit_row["sha256"],
            "control_image": str(control_path.resolve()),
            "control_image_sha256": control_audit_row["sha256"],
            "target_member_positions_in_image": member_counts.get(pair["target_cluster_id"], 0),
            "z_cl": pair["target"].get("z_cl"),
            "richness": pair["target"].get("richness"),
        })

    output.mkdir(parents=True, exist_ok=False)
    target_doc = dict(target_manifest)
    control_doc = dict(control_manifest)
    for doc, frames, population in (
        (target_doc, target_frames, "target"),
        (control_doc, control_frames, "control"),
    ):
        doc["clusters"] = frames
        doc["science_measurement"] = False
        doc["psf_approved"] = False
        doc["preparation_cohort"] = f"s23b_matched_psf_expansion_{population}"
        doc["selection_seed"] = SEED
        doc["selection_per_filter"] = PER_FILTER_PAIRS
        doc["center_is_bcg"] = False
        doc["member_catalog_import"] = str(paths["member_positions"].resolve()) if population == "target" else None
        doc["member_catalog_sha256"] = sha256(paths["member_positions"]) if population == "target" else None
        doc["notes"] = [
            "Selected from already received and audited HSC PDR3 coadd/bg images.",
            "Target/control pairs are matched by the previously screened depth/seeing/extinction procedure.",
            "Selection is stratified by target seeing in five equal-count strata per i-family filter.",
            "CAMIRA cluster centers and random control centers are not identified BCG positions.",
            "Target member coordinates are used only to screen nearby compact-source candidates.",
            "No Gaia/CAS request, image download, PSF approval, or science measurement is performed.",
        ]
    target_manifest_out = output / "target_manifest.json"
    control_manifest_out = output / "control_manifest.json"
    write_json(target_manifest_out, target_doc)
    write_json(control_manifest_out, control_doc)
    selection_out = output / "selection.tsv"
    write_tsv(selection_out, selection_rows, list(selection_rows[0]))

    source_hashes = {str(path.resolve()): sha256(path) for path in paths.values()}
    summary = {
        "schema": "astrafek-s23b-matched-psf-expansion-plan-v1",
        "status": "prepared_local_cpu_analysis_pending",
        "selection_seed": SEED,
        "stratification": "five equal-count target-seeing strata within each actual HSC-I family filter",
        "n_pairs_selected": len(selected),
        "n_frames_selected": len(target_frames) + len(control_frames),
        "n_pairs_by_filter": dict(Counter(pair["target_filter"] for pair, _ in selected)),
        "n_candidate_pairs_available_after_pilot_exclusion": {
            filt: len(pairs_by_filter[filt]) for filt in ("i", "i2")
        },
        "n_existing_target_pilot_pairs_excluded": len(prepared_ids),
        "n_filter_mismatched_pairs_excluded": len(filter_mismatch_ids),
        "filter_mismatch_examples": filter_mismatch_ids[:20],
        "n_member_positions_attached_to_targets": sum(len(row["members_radec_deg"]) for row in target_frames),
        "n_member_positions_inside_targets": sum(row["member_catalog_inside_frame_count"] for row in target_frames),
        "input_sha256": source_hashes,
        "target_manifest": str(target_manifest_out),
        "target_manifest_sha256": sha256(target_manifest_out),
        "control_manifest": str(control_manifest_out),
        "control_manifest_sha256": sha256(control_manifest_out),
        "selection_table": str(selection_out),
        "selection_table_sha256": sha256(selection_out),
        "target_image_sha256": {row["image"]: row["input_image_sha256"] for row in target_frames},
        "control_image_sha256": {row["image"]: row["input_image_sha256"] for row in control_frames},
        "external_queries_submitted": False,
        "images_downloaded": False,
        "science_measurement": False,
        "psf_approved": False,
        "limitations": [
            "selection excludes the 49 previously prepared target/control pilot pairs",
            "the 200 new pairs are a stratified diagnostic expansion, not the full S23B sample",
            "Gaia matches are required before calling compact candidates stars",
            "no coadd-epoch propagation, PSF wing calibration, sky transfer, or covariance correction is supplied",
        ],
    }
    summary_path = output / "summary.json"
    write_json(summary_path, summary)
    summary["summary_sha256"] = sha256(summary_path)
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    result = build(args.output)
    print(json.dumps({key: result[key] for key in (
        "status", "n_pairs_selected", "n_frames_selected", "n_pairs_by_filter",
        "n_candidate_pairs_available_after_pilot_exclusion",
        "n_member_positions_attached_to_targets", "external_queries_submitted",
        "science_measurement", "psf_approved", "target_manifest_sha256",
        "control_manifest_sha256", "summary_sha256",
    )}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
