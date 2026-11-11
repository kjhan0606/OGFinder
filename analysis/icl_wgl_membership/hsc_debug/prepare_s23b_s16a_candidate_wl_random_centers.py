#!/usr/bin/env python3
"""Select provisional five-per-lens WL random centers for S16A candidates.

This local-only selector reuses the established S23B random-candidate audit
and deterministic global-uniqueness matcher. It prepares coordinates for a
later count-only plan; it does not contact HSC CAS or certify controls as
science-ready.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path

from astropy.cosmology import FlatLambdaCDM

from prepare_s23b_wl_common_support_plan import select_unique_controls

DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
CONTROL_AUDIT = DATA / "icl_s23b_control_candidate_audit_2026-11-07_01"
CANDIDATE_AUDIT = CONTROL_AUDIT / "candidate_audit.csv"
CONTROL_AUDIT_SUMMARY = CONTROL_AUDIT / "summary.json"
SELECTION_AUDIT = DATA / "s23b_s16a_shape_coverage_selection_audit_2026-11-11_04"
SELECTION_SUMMARY = SELECTION_AUDIT / "summary.json"
ICL_CONTROL_TABLE = SELECTION_AUDIT / "s16a_ge1000_outside_s19a_icl_control_centers.tsv"
TARGET_MANIFEST = DATA / "s23b_icl_target_image_plan_2026-11-07_01" / "target_manifest.json"
DEFAULT_OUTPUT = DATA / "s23b_s16a_expansion_wl_random_centers_2026-11-11_01"

EXPECTED_CANDIDATE_AUDIT_PLAN_SHA256 = "143a1bd0a3a3c341576ce420f77ba0089fbfe80fdf58fbfd7afe30ea8e3f6714"
EXPECTED_CANDIDATE_AUDIT_SHA256 = "ac9cb72ca769b19ace23e2b90de3f6dd5676c9f3813a525e8f8319f6fae6cb14"
EXPECTED_SELECTION_SUMMARY_SHA256 = "c08602e1b45492a92389dabd074e5f854050e8e14c2cd55b4683b5ae6fc0b398"
EXPECTED_ICL_CONTROL_TABLE_SHA256 = "81c702166871e75570918f770ffad44702cebe9a29920c80d3d73d2263c21706"
EXPECTED_TARGET_MANIFEST_SHA256 = "e487582ba40cd3e58da67704242db4d7219e381fb08f4b33fed47ce4ff268b44"
MIN_TARGET_CONTROL_SEPARATION_MPC = 6.0
CONTROLS_PER_LENS = 5
H0 = 70.0
OMEGA_M = 0.3


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_json(path: Path) -> dict:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def read_tsv(path: Path) -> list[dict[str, str]]:
    with Path(path).open(encoding="utf-8-sig", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def build(output: Path) -> dict:
    output = output.resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite existing random-center output: {output}")
    if (sha256(CANDIDATE_AUDIT) != EXPECTED_CANDIDATE_AUDIT_SHA256
            or sha256(TARGET_MANIFEST) != EXPECTED_TARGET_MANIFEST_SHA256
            or sha256(SELECTION_SUMMARY) != EXPECTED_SELECTION_SUMMARY_SHA256
            or sha256(ICL_CONTROL_TABLE) != EXPECTED_ICL_CONTROL_TABLE_SHA256):
        raise ValueError("one or more pinned input files changed")

    candidate_summary = read_json(CONTROL_AUDIT_SUMMARY)
    if (candidate_summary.get("schema") != "astrafek-s23b-icl-random-control-candidate-audit-v1"
            or candidate_summary.get("status") != "complete_provisional_candidate_audit"
            or candidate_summary.get("plan_sha256") != EXPECTED_CANDIDATE_AUDIT_PLAN_SHA256
            or candidate_summary.get("science_measurement") is not False
            or candidate_summary.get("controls_ready_for_science") != 0
            or candidate_summary.get("outputs_sha256", {}).get("candidate_audit.csv")
            != EXPECTED_CANDIDATE_AUDIT_SHA256):
        raise ValueError("saved random-candidate audit is incomplete or changed")
    selection_summary = read_json(SELECTION_SUMMARY)
    selected_scope = selection_summary.get("expansion_candidate_icl_control_pairs", {})
    if (selection_summary.get("status") != "complete_local_descriptive_audit"
            or selection_summary.get("science_measurement") is not False
            or selection_summary.get("cas_query_submitted") is not False
            or selected_scope.get("n_candidate_targets") != 126
            or selected_scope.get("n_targets_with_saved_icl_control_pair") != 116):
        raise ValueError("S16A selection audit does not match the reviewed 116-pair sample")

    targets_doc = read_json(TARGET_MANIFEST)
    if (targets_doc.get("schema") != "hsc-bcg-diffuse-sample-v1"
            or targets_doc.get("science_measurement") is not False):
        raise ValueError("unexpected S23B target manifest")
    targets = {str(row["id"]): row for row in targets_doc.get("clusters", [])}
    if len(targets) != 1043:
        raise ValueError("pinned target manifest does not contain 1,043 clusters")
    with ICL_CONTROL_TABLE.open(encoding="utf-8-sig", newline="") as stream:
        selected_rows = list(csv.DictReader(stream, delimiter="\t"))
    if len(selected_rows) != 116:
        raise ValueError("ICL control table does not contain 116 paired targets")

    target_rows = {}
    for row in selected_rows:
        lens_id = str(row["target_id"])
        if lens_id in target_rows or lens_id not in targets:
            raise ValueError(f"unexpected or repeated S16A candidate target: {lens_id}")
        target = targets[lens_id]
        if (str(row["control_s16a_count_status"]) != "not_queried"
                or int(row["s16a_target_count_3pmpc"]) < 1000
                or float(row["z_cl"]) != float(target["z_cl"])
                or str(target.get("sample_role")) != "target"):
            raise ValueError(f"candidate target metadata or raw-count screen changed: {lens_id}")
        catalog_id = str(target.get("catalog_id", ""))
        if not catalog_id:
            raise ValueError(f"manifest lacks catalog ID for {lens_id}")
        target_rows[lens_id] = {
            "lens_id": lens_id,
            "catalog_id": catalog_id,
            "z_lens": float(target["z_cl"]),
            "ra_deg": float(target["ra_deg"]),
            "dec_deg": float(target["dec_deg"]),
            "richness": float(target["richness"]),
            "icl_filter": str(row["icl_filter"]),
            "target_s16a_raw_count": int(row["s16a_target_count_3pmpc"]),
        }

    cosmology = FlatLambdaCDM(H0=H0, Om0=OMEGA_M, Tcmb0=0)
    groups: dict[str, list[dict]] = {lens_id: [] for lens_id in target_rows}
    filter_stage_counts = {
        "candidate_rows": 0,
        "center_flags_pass": 0,
        "metadata_match": 0,
        "cluster_overlap_2mpc_zero": 0,
        "separation_gt_6_mpc": 0,
    }
    with CANDIDATE_AUDIT.open(encoding="utf-8", newline="") as stream:
        for row in csv.DictReader(stream):
            lens_id = str(row["target_cluster_id"]).replace(":", "_")
            if lens_id not in target_rows:
                continue
            filter_stage_counts["candidate_rows"] += 1
            target = target_rows[lens_id]
            if str(row.get("target_catalog_id")) != target["catalog_id"]:
                raise ValueError(f"candidate catalog ID does not match target: {lens_id}")
            if not math.isclose(float(row["target_z"]), target["z_lens"], rel_tol=0, abs_tol=1e-9):
                raise ValueError(f"candidate redshift does not match target: {lens_id}")
            if row.get("center_flags_pass") != "True":
                continue
            filter_stage_counts["center_flags_pass"] += 1
            if row.get("metadata_match") != "True":
                continue
            filter_stage_counts["metadata_match"] += 1
            if int(row["cluster_overlap_2mpc"]) != 0:
                continue
            filter_stage_counts["cluster_overlap_2mpc_zero"] += 1
            z = target["z_lens"]
            da = float(cosmology.angular_diameter_distance(z).value)
            offset_mpc = float(row["offset_from_target_arcmin"]) * da * math.pi / 10800.0
            if not offset_mpc > MIN_TARGET_CONTROL_SEPARATION_MPC:
                continue
            filter_stage_counts["separation_gt_6_mpc"] += 1
            groups[lens_id].append({
                "random_object_id": str(row["random_object_id"]),
                "ra_deg": float(row["ra_deg"]),
                "dec_deg": float(row["dec_deg"]),
                "target_z": z,
                "offset_from_target_mpc": offset_mpc,
                "offset_from_target_arcmin": float(row["offset_from_target_arcmin"]),
                "match_score": float(row["match_score"]),
                "tract": str(row["tract"]),
                "patch": str(row["patch"]),
                "skymap_id": str(row["skymap_id"]),
                "delta_depth_mag": float(row["delta_depth_mag"]),
                "seeing_ratio": float(row["seeing_ratio"]),
                "delta_ebv_sfd": float(row["delta_ebv_sfd"]),
                "cluster_overlap_2mpc": int(row["cluster_overlap_2mpc"]),
            })

    for lens_id, options in groups.items():
        unique = {}
        for option in options:
            position = (round(option["ra_deg"], 12), round(option["dec_deg"], 12))
            if position not in unique or (option["match_score"], option["random_object_id"]) < (
                    unique[position]["match_score"], unique[position]["random_object_id"]):
                unique[position] = option
        groups[lens_id] = sorted(unique.values(), key=lambda item: (
            item["match_score"], item["random_object_id"]))

    eligible = {lens_id for lens_id, options in groups.items()
                if len(options) >= CONTROLS_PER_LENS}
    assignment = select_unique_controls(groups, CONTROLS_PER_LENS)
    assigned_lenses = {slot[0] for slot in assignment}
    if assigned_lenses != eligible:
        raise ValueError("global unique matching did not preserve every eligible target")

    rows = []
    random_ids = set()
    for lens_id in sorted(target_rows):
        if lens_id not in assigned_lenses:
            continue
        target = target_rows[lens_id]
        for slot_index in range(CONTROLS_PER_LENS):
            candidate = assignment[(lens_id, slot_index)]
            random_object_id = candidate["random_object_id"]
            if random_object_id in random_ids:
                raise ValueError(f"random-table object ID reused: {random_object_id}")
            random_ids.add(random_object_id)
            center_id = "RC-" + hashlib.sha256(
                f"{lens_id}\0{random_object_id}\0{candidate['ra_deg']:.12f}\0{candidate['dec_deg']:.12f}".encode()
            ).hexdigest()[:20]
            rows.append({
                "random_center_id": center_id,
                "matched_lens_id": lens_id,
                "matched_catalog_id": target["catalog_id"],
                "icl_filter": target["icl_filter"],
                "ra_deg": candidate["ra_deg"],
                "dec_deg": candidate["dec_deg"],
                "z_lens": target["z_lens"],
                "target_s16a_raw_count": target["target_s16a_raw_count"],
                "random_object_id": random_object_id,
                "control_slot": slot_index + 1,
                "offset_from_target_mpc": candidate["offset_from_target_mpc"],
                "offset_from_target_arcmin": candidate["offset_from_target_arcmin"],
                "match_score": candidate["match_score"],
                "tract": candidate["tract"],
                "patch": candidate["patch"],
                "skymap_id": candidate["skymap_id"],
                "delta_depth_mag": candidate["delta_depth_mag"],
                "seeing_ratio": candidate["seeing_ratio"],
                "delta_ebv_sfd": candidate["delta_ebv_sfd"],
                "cluster_overlap_2mpc": candidate["cluster_overlap_2mpc"],
            })
    selected_lenses = {row["matched_lens_id"] for row in rows}
    if len(rows) != len(selected_lenses) * CONTROLS_PER_LENS:
        raise ValueError("selected random-control rows do not have exactly five centers per target")

    output.mkdir(parents=True, exist_ok=False)
    fields = list(rows[0]) if rows else [
        "random_center_id", "matched_lens_id", "matched_catalog_id", "icl_filter",
        "ra_deg", "dec_deg", "z_lens", "target_s16a_raw_count", "random_object_id",
        "control_slot", "offset_from_target_mpc", "offset_from_target_arcmin", "match_score",
        "tract", "patch", "skymap_id", "delta_depth_mag", "seeing_ratio", "delta_ebv_sfd",
        "cluster_overlap_2mpc",
    ]
    centers_path = output / "random_centers.tsv"
    with centers_path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    selection_path = output / "target_selection.tsv"
    selection_fields = [
        "matched_lens_id", "matched_catalog_id", "icl_filter", "z_lens",
        "target_s16a_raw_count", "n_eligible_random_candidates",
        "n_selected_random_centers", "selection_status", "selection_reason",
    ]
    with selection_path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=selection_fields, delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        for lens_id in sorted(target_rows):
            target = target_rows[lens_id]
            n_eligible = len(groups[lens_id])
            n_selected = CONTROLS_PER_LENS if lens_id in selected_lenses else 0
            if n_selected:
                status = "selected"
                reason = "five eligible globally unique random centers"
            else:
                status = "excluded"
                reason = f"only {n_eligible} eligible candidates; need {CONTROLS_PER_LENS}"
            writer.writerow({
                "matched_lens_id": lens_id,
                "matched_catalog_id": target["catalog_id"],
                "icl_filter": target["icl_filter"],
                "z_lens": target["z_lens"],
                "target_s16a_raw_count": target["target_s16a_raw_count"],
                "n_eligible_random_candidates": n_eligible,
                "n_selected_random_centers": n_selected,
                "selection_status": status,
                "selection_reason": reason,
            })
    by_filter = {}
    for band in ("i", "i2"):
        band_candidates = {lens_id for lens_id, target in target_rows.items()
                           if target["icl_filter"] == band}
        ids = {lens_id for lens_id in selected_lenses
               if target_rows[lens_id]["icl_filter"] == band}
        by_filter[band] = {
            "n_candidate_lenses": len(band_candidates),
            "n_lenses_with_five_eligible_candidates": sum(
                len(groups[lens_id]) >= CONTROLS_PER_LENS for lens_id in band_candidates),
            "n_lenses_selected": len(ids),
            "n_lenses_excluded": len(band_candidates) - len(ids),
            "n_random_centers": sum(row["icl_filter"] == band for row in rows),
        }
    summary = {
        "schema": "astrafek-s23b-s16a-expansion-wl-random-centers-v2",
        "status": "complete_local_provisional_selection",
        "science_measurement": False,
        "network_used": False,
        "cas_query_submitted": False,
        "source_rows_downloaded": False,
        "images_read": False,
        "n_candidate_targets": len(target_rows),
        "n_targets_with_at_least_five_eligible_candidates": len(eligible),
        "n_targets_with_five_globally_unique_random_centers": len(selected_lenses),
        "n_random_centers": len(rows),
        "controls_per_lens": CONTROLS_PER_LENS,
        "minimum_target_control_separation_mpc": MIN_TARGET_CONTROL_SEPARATION_MPC,
        "candidate_filter_stages": filter_stage_counts,
        "eligible_candidates_per_target": {
            "min": min((len(groups[lens_id]) for lens_id in target_rows), default=0),
            "median": float(__import__("numpy").median([len(groups[lens_id]) for lens_id in target_rows])) if target_rows else None,
            "max": max((len(groups[lens_id]) for lens_id in target_rows), default=0),
        },
        "by_icl_filter": by_filter,
        "candidate_audit_plan_sha256": EXPECTED_CANDIDATE_AUDIT_PLAN_SHA256,
        "candidate_audit_csv_sha256": sha256(CANDIDATE_AUDIT),
        "selection_summary_sha256": sha256(SELECTION_SUMMARY),
        "icl_control_table_sha256": sha256(ICL_CONTROL_TABLE),
        "target_manifest_sha256": sha256(TARGET_MANIFEST),
        "random_center_table_sha256": sha256(centers_path),
        "target_selection_table_sha256": sha256(selection_path),
        "excluded_target_ids": sorted(set(target_rows) - selected_lenses),
        "controls_ready_for_science": False,
        "control_quality_limit": "The upstream provisional control audit has controls_ready_for_science=0. These centers pass its local metadata/flag, CAMIRA 2-Mpc, >6 proper-Mpc separation, and global-uniqueness rules only; they are not image-QA-certified WL controls.",
        "existing_saved_wl_random_centers_covering_candidates": 4,
        "later_gates": [
            "prepare an exact S16A count-only plan for these centers and obtain separate approval before CAS submission",
            "verify local source support and sample-wide common cuts after the count-only result",
            "perform image/PSF/mask/control quality gates before any weak-lensing shape-row export",
            "obtain separate explicit approval before source-row transfer",
        ],
    }
    summary_path = output / "summary.json"
    summary_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    return summary


def main() -> None:
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    print(json.dumps(build(args.output), ensure_ascii=False, indent=2))
    print(f"output={args.output.resolve()}")


if __name__ == "__main__":
    main()
