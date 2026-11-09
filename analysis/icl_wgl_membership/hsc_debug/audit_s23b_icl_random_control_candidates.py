#!/usr/bin/env python3
"""Audit returned S23B random positions and prepare provisional ICL controls.

This is a local-only postprocessor for the explicitly approved 29-column HSC
random-position export. It verifies every batch, estimates target/control
metadata compatibility, removes positions near known CAMIRA clusters, and
selects a reproducible one-control-per-target non-overlapping queue. Its output
is a retrieval candidate list, not science-ready imagery or a measurement.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import platform

import numpy as np
from astropy.coordinates import SkyCoord
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
from scipy.spatial import cKDTree

from fetch_ancillary_catalogues import catalogue_reader
from prepare_ancillary import query_sfd
from submit_s23b_icl_control_candidates import (
    DEFAULT_PLAN, DEFAULT_RESULTS, EXPECTED_PLAN_SHA256, OUTPUT_COLUMNS,
    load_plan, sha256, validate_csv,
)


DATA = Path("/gpfs/kjhan/Astrafek/data")
DEFAULT_S20A = DATA / "catalog/camira_s20a_wide_sm.tbl"
DEFAULT_S23B = DATA / "catalogue_s23b_camira_2026-10-20_01/camira_s23b_wide_sm.tbl"
DEFAULT_ANCILLARY = DATA / "sample_800_2026-10-13/ancillary_2026-10-13_01/preparation_status.json"
DEFAULT_OUTPUT = DATA / "sample_800_2026-10-13/icl_s23b_control_candidate_audit_2026-11-07_01"
EXPECTED_S20A_SHA256 = "d929985d7f6ee38c4435fc12099e9df1e5dbfd9ce9603cea319d5704129d22aa"
EXPECTED_S23B_SHA256 = "5af057028feab7976cd3ed389a7246cbe73add8bbb96891101f0f692b1de4d71"
RADII_MPC = (0.5, 1.0, 2.0)
DEPTH_TOLERANCE_MAG = 0.30
SEEING_RATIO_MAX = 1.20
EBV_TOLERANCE = 0.02
BAD_FLAGS = (
    "i_inputcount_flag", "i_inputcount_flag_noinputs", "i_inputcount_flag_badcentroid",
    "i_pixelflags_offimage", "i_pixelflags_edge", "i_pixelflags_bad",
    "i_pixelflags_interpolatedcenter", "i_pixelflags_saturatedcenter",
    "i_pixelflags_crcenter", "i_pixelflags_suspectcenter",
    "i_pixelflags_bright_objectcenter",
)


def unit_vectors(ra: np.ndarray, dec: np.ndarray) -> np.ndarray:
    ra, dec = np.deg2rad(ra), np.deg2rad(dec)
    return np.column_stack((np.cos(dec) * np.cos(ra),
                            np.cos(dec) * np.sin(ra), np.sin(dec)))


def separation_deg(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    return np.rad2deg(2 * np.arcsin(np.clip(np.linalg.norm(a - b, axis=-1) / 2, 0, 1)))


def number(value: object) -> float | None:
    try:
        result = float(value)
    except (TypeError, ValueError):
        return None
    return result if math.isfinite(result) else None


def as_bool(value: str) -> bool | None:
    value = str(value).strip().lower()
    if value in {"true", "t", "1"}:
        return True
    if value in {"false", "f", "0"}:
        return False
    return None


def center_flags_pass(row: dict[str, str]) -> bool:
    return (
        all(as_bool(row.get(key, "")) is True for key in
            ("isprimary", "detect_ispatchinner", "detect_istractinner"))
        and all(as_bool(row.get(key, "")) is False for key in BAD_FLAGS)
        and (number(row.get("i_inputcount_value")) or 0) > 0
        and (number(row.get("i_pix_variance")) or 0) > 0
    )


def read_camira(path: Path) -> list[dict]:
    rows = []
    with path.open(encoding="utf-8") as stream:
        for line_no, line in enumerate(stream, 1):
            if not line.strip() or line.lstrip().startswith("#"):
                continue
            fields = line.split()
            if len(fields) != 7:
                raise ValueError(f"{path}:{line_no}: expected 7 CAMIRA columns")
            ident, _name = fields[:2]
            ra, dec, z, richness, _z_bcg = map(float, fields[2:])
            if (not all(math.isfinite(v) for v in (ra, dec, z, richness))
                    or not (0 <= ra < 360 and -90 <= dec <= 90 and z > 0 and richness > 0)):
                raise ValueError(f"{path}:{line_no}: invalid CAMIRA coordinate or property")
            rows.append({"id": ident, "ra_deg": ra, "dec_deg": dec, "z": z,
                         "richness": richness})
    if len({row["id"] for row in rows}) != len(rows):
        raise ValueError(f"duplicate CAMIRA IDs in {path}")
    return rows


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def read_batch(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as stream:
        reader = catalogue_reader(stream)
        if reader.fieldnames != OUTPUT_COLUMNS:
            raise ValueError(f"columns differ from the approved 29-column schema: {path.name}")
        rows = list(reader)
    return rows


def candidate_match_rows(rows: list[dict], radius_mpc: float) -> list[dict]:
    groups: dict[str, list[dict]] = {}
    for row in rows:
        if row["metadata_match"] and row[f"cluster_overlap_{radius_mpc:g}mpc"] == 0:
            groups.setdefault(row["target_cluster_id"], []).append(row)
    chosen, used_ids, points, caps = [], set(), [], []
    for cluster_id in sorted(groups, key=lambda key: (len(groups[key]), key)):
        choices = sorted(groups[cluster_id], key=lambda r: (
            r["match_score"], r["offset_from_target_arcmin"], r["random_object_id"]))
        for row in choices:
            ident = row["random_object_id"]
            if ident in used_ids:
                continue
            point = unit_vectors(np.array([row["ra_deg"]]), np.array([row["dec_deg"]]))[0]
            cap = row["footprint_cap_deg"]
            if points and np.any(separation_deg(point, np.asarray(points)) <= cap + np.asarray(caps)):
                continue
            chosen.append(row)
            used_ids.add(ident)
            points.append(point)
            caps.append(cap)
            break
    return chosen


def run(plan_dir: Path, results_dir: Path, s20a_path: Path, s23b_path: Path,
        ancillary_path: Path, output_dir: Path) -> dict:
    if not platform.node().lower().startswith("lageunha"):
        raise RuntimeError("Run this local metadata audit on LagEunha CPU")
    if output_dir.exists():
        raise FileExistsError(f"refusing to overwrite {output_dir}")
    plan, targets, batches, plan_hash = load_plan(plan_dir)
    status_path = results_dir / "status.json"
    status = json.loads(status_path.read_text(encoding="utf-8"))
    if (status.get("schema") != "astrafek-s23b-icl-control-candidate-run-v1"
            or status.get("status") != "complete"
            or status.get("plan_sha256") != EXPECTED_PLAN_SHA256
            or status.get("network_used") is not True
            or status.get("query_submitted") is not True
            or len(status.get("jobs", {})) != len(batches)
            or any(job.get("status") != "complete" for job in status["jobs"].values())):
        raise ValueError("approved S23B metadata export is incomplete or has an unexpected scope")
    if any(status.get(key) is not False for key in (
            "credentials_saved", "images_returned", "source_detections_returned",
            "shape_rows_returned", "photoz_rows_returned", "science_measurement")):
        raise ValueError("metadata export provenance exceeds the approved scope")
    for path, expected in ((s20a_path, EXPECTED_S20A_SHA256), (s23b_path, EXPECTED_S23B_SHA256)):
        if not path.is_file() or sha256(path) != expected:
            raise ValueError(f"CAMIRA input hash mismatch: {path}")

    candidate_by_id = {row["candidate_id"]: row for row in targets}
    if len(candidate_by_id) != 1043:
        raise ValueError("expected 1,043 screened S23B targets")
    raw_rows: list[dict[str, str]] = []
    plan_json = Path(plan_dir) / "control_candidate_query_plan.json"
    candidate_table = Path(plan["candidate_table"])
    inputs = {
        str(status_path): sha256(status_path),
        str(plan_json): sha256(plan_json),
        str(candidate_table): sha256(candidate_table),
    }
    inputs.update({str(Path(batch["path"])): batch["sha256"] for batch in batches})
    for index, batch in enumerate(batches):
        name = batch["name"]
        result_path = results_dir / f"{name}.csv"
        job = status["jobs"].get(name)
        if (not job or not result_path.is_file() or sha256(result_path) != job.get("sha256")
                or job.get("status") != "complete"):
            raise ValueError(f"missing or changed returned batch {name}")
        expected_ids = {row["candidate_id"] for row in targets[
            index * 10:index * 10 + int(batch["n_candidates"])]}
        n_rows, _nulls = validate_csv(result_path, expected_ids)
        if n_rows != int(job.get("rows", -1)):
            raise ValueError(f"checkpoint row count differs for {name}")
        raw_rows.extend(read_batch(result_path))
        inputs[str(result_path)] = job["sha256"]
    if len(raw_rows) != sum(int(job["rows"]) for job in status["jobs"].values()):
        raise ValueError("independent aggregate row count differs from checkpoint")

    catalogues = read_camira(s20a_path) + read_camira(s23b_path)
    cluster_vectors = unit_vectors(np.array([r["ra_deg"] for r in catalogues]),
                                   np.array([r["dec_deg"] for r in catalogues]))
    cluster_tree = cKDTree(cluster_vectors)
    cosmo = FlatLambdaCDM(H0=70.0, Om0=0.3)
    cluster_deg_per_mpc = np.array([
        math.degrees(1 / float(cosmo.angular_diameter_distance(row["z"]).value))
        for row in catalogues
    ])
    map_status = json.loads(ancillary_path.read_text(encoding="utf-8"))
    maps = {}
    for pole, item in map_status.get("maps", {}).items():
        path = Path(item["file"])
        if not path.is_file() or sha256(path) != item["sha256"]:
            raise ValueError(f"SFD map missing or changed: {path}")
        maps[pole] = path
        inputs[str(path)] = item["sha256"]
    if set(maps) != {"ngp", "sgp"}:
        raise ValueError("both verified SFD maps are required")
    inputs[str(s20a_path)] = EXPECTED_S20A_SHA256
    inputs[str(s23b_path)] = EXPECTED_S23B_SHA256
    inputs[str(ancillary_path)] = sha256(ancillary_path)
    target_ebv = query_sfd(SkyCoord(
        [float(r["ra_deg"]) for r in targets] * u.deg,
        [float(r["dec_deg"]) for r in targets] * u.deg), maps)
    target_ebv = {r["candidate_id"]: float(v) for r, v in zip(targets, target_ebv)}
    point_ebv = query_sfd(SkyCoord(
        [float(r["ra"]) for r in raw_rows] * u.deg,
        [float(r["dec"]) for r in raw_rows] * u.deg), maps)

    unique: dict[str, tuple[float, float, str]] = {}
    records, seen_pair = [], set()
    point_vectors = unit_vectors(np.array([float(r["ra"]) for r in raw_rows]),
                                 np.array([float(r["dec"]) for r in raw_rows]))
    for index, (row, ebv) in enumerate(zip(raw_rows, point_ebv)):
        target_id, object_id = row["candidate_id"], str(row["random_object_id"])
        target = candidate_by_id.get(target_id)
        if target is None:
            raise ValueError(f"unrequested target ID {target_id}")
        pair = (target_id, object_id)
        if pair in seen_pair:
            raise ValueError(f"duplicate random object within target {target_id}")
        seen_pair.add(pair)
        ra, dec = float(row["ra"]), float(row["dec"])
        if object_id in unique and unique[object_id] != (ra, dec, row["skymap_id"]):
            raise ValueError(f"random object {object_id} has inconsistent sky position")
        unique[object_id] = (ra, dec, row["skymap_id"])
        target_vec = unit_vectors(np.array([float(target["ra_deg"])]),
                                  np.array([float(target["dec_deg"])]))[0]
        offset_arcmin = float(separation_deg(point_vectors[index], target_vec)) * 60
        if offset_arcmin > 30.00001:
            raise ValueError("returned random point lies outside the approved 30 arcmin search radius")
        depth, seeing = number(row.get("imag_psf_depth")), number(row.get("iseeing"))
        target_depth, target_seeing = number(target.get("median_imag_psf_depth")), number(target.get("median_seeing_arcsec"))
        qa_pass = (center_flags_pass(row) and as_bool(row.get("icalexp_exists", "")) is True
                   and depth is not None and depth > 0 and seeing is not None and seeing > 0)
        dd = (depth - target_depth) if depth is not None and target_depth is not None else None
        ratio = (seeing / target_seeing) if seeing is not None and target_seeing and target_seeing > 0 else None
        de = float(ebv) - target_ebv[target_id] if math.isfinite(float(ebv)) else None
        metadata_match = bool(qa_pass and dd is not None and ratio is not None and de is not None
                             and abs(dd) <= DEPTH_TOLERANCE_MAG
                             and 1 / SEEING_RATIO_MAX <= ratio <= SEEING_RATIO_MAX
                             and abs(de) <= EBV_TOLERANCE)
        score = ((dd / DEPTH_TOLERANCE_MAG) ** 2 +
                 (math.log(ratio) / math.log(SEEING_RATIO_MAX)) ** 2 +
                 (de / EBV_TOLERANCE) ** 2) if metadata_match else None
        z = float(target["z_cl"])
        half_width_arcsec = 1000.0 / float(cosmo.kpc_proper_per_arcmin(z).value) * 60.0
        footprint_cap_deg = math.sqrt(2) * half_width_arcsec / 3600 + 1 / 3600
        search_deg = footprint_cap_deg + max(RADII_MPC) * float(np.max(cluster_deg_per_mpc))
        nearby = cluster_tree.query_ball_point(
            point_vectors[index], 2 * math.sin(math.radians(search_deg) / 2))
        overlap = {}
        for radius in RADII_MPC:
            n = 0
            for ci in nearby:
                cluster = catalogues[ci]
                one_mpc_deg = cluster_deg_per_mpc[ci]
                sep = float(separation_deg(point_vectors[index], cluster_vectors[ci]))
                if sep <= footprint_cap_deg + radius * one_mpc_deg:
                    n += 1
            overlap[f"cluster_overlap_{radius:g}mpc"] = n
        records.append({
            "target_cluster_id": target_id, "target_catalog_id": target["catalog_id"],
            "random_object_id": object_id, "ra_deg": ra, "dec_deg": dec,
            "skymap_id": row["skymap_id"], "tract": row["tract"], "patch": row["patch"],
            "offset_from_target_arcmin": offset_arcmin, "target_z": z,
            "target_richness": float(target["richness"]), "target_depth_mag": target_depth,
            "candidate_depth_mag": depth, "delta_depth_mag": dd,
            "target_seeing_arcsec": target_seeing, "candidate_seeing_arcsec": seeing,
            "seeing_ratio": ratio, "target_ebv_sfd": target_ebv[target_id],
            "candidate_ebv_sfd": float(ebv), "delta_ebv_sfd": de,
            "center_flags_pass": qa_pass, "metadata_match": metadata_match,
            "match_score": score, "footprint_cap_deg": footprint_cap_deg,
            **overlap,
        })
    if len(unique) > len(raw_rows):
        raise AssertionError("unique position count cannot exceed row count")

    sensitivity, plans = [], []
    for radius in RADII_MPC:
        candidates = [r for r in records if r["metadata_match"]
                      and r[f"cluster_overlap_{radius:g}mpc"] == 0]
        chosen = candidate_match_rows(records, radius)
        sensitivity.append({
            "cluster_exclusion_mpc": radius,
            "candidate_positions": len(candidates),
            "targets_with_candidate_positions": len({r["target_cluster_id"] for r in candidates}),
            "assigned_nonoverlap_controls": len(chosen),
        })
        for row in chosen:
            plans.append({**row, "cluster_exclusion_mpc": radius,
                          "status": "provisional_control_pending_image_QA",
                          "science_ready": False})

    output_dir.mkdir(parents=True, exist_ok=False)
    fields = list(records[0]) if records else ["target_cluster_id", "random_object_id"]
    with (output_dir / "candidate_audit.csv").open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
        writer.writeheader(); writer.writerows(records)
    plan_fields = list(dict.fromkeys(key for row in plans for key in row)) or fields
    with (output_dir / "provisional_controls.csv").open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=plan_fields, lineterminator="\n")
        writer.writeheader(); writer.writerows(plans)
    summary = {
        "schema": "astrafek-s23b-icl-random-control-candidate-audit-v1",
        "status": "complete_provisional_candidate_audit",
        "plan_sha256": plan_hash,
        "n_targets": len(targets), "n_batches": len(batches), "n_returned_rows": len(raw_rows),
        "n_unique_random_objects": len(unique),
        "n_center_qa_pass": sum(bool(r["center_flags_pass"]) for r in records),
        "n_metadata_match": sum(bool(r["metadata_match"]) for r in records),
        "match_thresholds": {"depth_abs_mag": DEPTH_TOLERANCE_MAG,
                             "seeing_ratio_max": SEEING_RATIO_MAX,
                             "sfd_ebv_abs": EBV_TOLERANCE},
        "sensitivity": sensitivity,
        "science_measurement": False, "controls_ready_for_science": 0,
        "images_downloaded": 0, "network_used_by_audit": False,
        "inputs_sha256": inputs,
        "outputs_sha256": {},
        "limitations": [
            "matching tolerances and 0.5/1/2 proper-Mpc exclusion radii are declared screening choices, not fitted results",
            "HSC random-table center flags and patch-level QA are metadata; full cutout WCS, masks, backgrounds, variance and PSF remain unverified",
            "the 1-Mpc square-cutout exclusion is conservatively represented by a circular cap around each image center",
            "CAMIRA catalogs are incomplete tracers of projected cluster contamination; the S20A and S23B entries are unioned for exclusion",
            "one-per-target greedy assignment is provisional and may favor dense, well-covered regions",
            "no cutouts were requested or downloaded; each selected target/control image requires a separate exact-scope approval",
        ],
    }
    for name in ("candidate_audit.csv", "provisional_controls.csv"):
        summary["outputs_sha256"][name] = sha256(output_dir / name)
    (output_dir / "summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({key: summary[key] for key in (
        "status", "n_targets", "n_returned_rows", "n_unique_random_objects",
        "n_center_qa_pass", "n_metadata_match", "sensitivity", "science_measurement")}, indent=2))
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan-dir", type=Path, default=DEFAULT_PLAN)
    parser.add_argument("--results-dir", type=Path, default=DEFAULT_RESULTS)
    parser.add_argument("--s20a-catalog", type=Path, default=DEFAULT_S20A)
    parser.add_argument("--s23b-catalog", type=Path, default=DEFAULT_S23B)
    parser.add_argument("--ancillary-status", type=Path, default=DEFAULT_ANCILLARY)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    run(args.plan_dir, args.results_dir, args.s20a_catalog, args.s23b_catalog,
        args.ancillary_status, args.output_dir)


if __name__ == "__main__":
    main()
