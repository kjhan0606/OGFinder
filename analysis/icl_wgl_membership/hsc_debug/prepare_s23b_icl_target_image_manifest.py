#!/usr/bin/env python3
"""Prepare the approved S23B HSC-I cluster-center cutout manifest locally.

No network access occurs here. The manifest contains the 1,043 S23B CAMIRA
centers that passed the nine-probe HSC-I coverage screen. Each requested
coadd/bg cutout includes image, mask, and variance planes over a 1-Mpc proper
half-width. Matched random-control cutouts are handled by a separate plan.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import re

import numpy as np
from astropy.cosmology import FlatLambdaCDM
from scipy.spatial import cKDTree

from download_pilot import MAX_FILE_BYTES, MAX_SAMPLE_FRAMES, estimate_bytes
from submit_s23b_icl_control_candidates import (
    DEFAULT_PLAN, EXPECTED_PLAN_SHA256, load_plan, sha256,
)


DATA = Path("/gpfs/kjhan/Astrafek/data")
DEFAULT_PARENT_MANIFEST = DATA / "sample_800_2026-10-13/sample_manifest.json"
DEFAULT_S23B_CATALOG = DATA / "catalogue_s23b_camira_2026-10-20_01/camira_s23b_wide_sm.tbl"
DEFAULT_OUTPUT = DATA / "sample_800_2026-10-13/s23b_icl_target_image_plan_2026-11-07_01"
EXPECTED_PARENT_SHA256 = "95388814062e95f6e8429e2543a5432433515c9d9607680d2e873488d2758401"
EXPECTED_S23B_SHA256 = "5af057028feab7976cd3ed389a7246cbe73add8bbb96891101f0f692b1de4d71"
EXPECTED_PLAN_SHA256 = "143a1bd0a3a3c341576ce420f77ba0089fbfe80fdf58fbfd7afe30ea8e3f6714"
HALF_WIDTH_KPC = 1000.0
MIN_SEPARATION_ARCSEC = 60.0


def unit_vectors(ra: np.ndarray, dec: np.ndarray) -> np.ndarray:
    ra, dec = np.deg2rad(ra), np.deg2rad(dec)
    return np.column_stack((np.cos(dec) * np.cos(ra),
                            np.cos(dec) * np.sin(ra), np.sin(dec)))


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def write_json(path: Path, value: dict) -> None:
    with path.open("x", encoding="utf-8") as stream:
        json.dump(value, stream, ensure_ascii=False, indent=2, allow_nan=False)
        stream.write("\n")


def prepare(plan_dir: Path, parent_manifest: Path, s23b_catalog: Path,
            output_dir: Path) -> dict:
    output_dir = output_dir.resolve()
    if output_dir.exists():
        raise FileExistsError(f"refusing to overwrite {output_dir}")
    plan, candidates, _batches, plan_hash = load_plan(plan_dir)
    if plan_hash != EXPECTED_PLAN_SHA256 or len(candidates) != 1043:
        raise ValueError("eligible S23B candidate plan differs from the approved expansion sample")
    for path, expected in ((parent_manifest, EXPECTED_PARENT_SHA256),
                           (s23b_catalog, EXPECTED_S23B_SHA256)):
        if not path.is_file() or sha256(path) != expected:
            raise ValueError(f"pinned input checksum mismatch: {path}")
    parent = json.loads(parent_manifest.read_text(encoding="utf-8"))
    if (parent.get("schema") != "hsc-bcg-diffuse-sample-v1"
            or parent.get("science_measurement") is not False
            or len(parent.get("clusters", [])) != 800):
        raise ValueError("unexpected original Astrafek ICL sample manifest")
    candidate_path = Path(plan["candidate_table"]).resolve()
    if sha256(candidate_path) != plan.get("candidate_table_sha256"):
        raise ValueError("eligible S23B candidate table hash differs from the approved plan")

    # Reject coordinates already represented in the 800-cluster HSC sample.
    old = parent["clusters"]
    old_vectors = unit_vectors(np.array([float(row["ra_deg"]) for row in old]),
                               np.array([float(row["dec_deg"]) for row in old]))
    old_tree = cKDTree(old_vectors)
    new_vectors = unit_vectors(np.array([float(row["ra_deg"]) for row in candidates]),
                               np.array([float(row["dec_deg"]) for row in candidates]))
    chord, nearest = old_tree.query(new_vectors, k=1)
    nearest_arcsec = np.rad2deg(2 * np.arcsin(np.clip(chord / 2, 0, 1))) * 3600
    if np.any(nearest_arcsec <= MIN_SEPARATION_ARCSEC):
        bad = int(np.argmin(nearest_arcsec))
        raise ValueError(f"S23B target is within 60 arcsec of existing sample: {candidates[bad]['candidate_id']}")

    cosmo = FlatLambdaCDM(H0=float(parent["h0"]), Om0=float(parent["omega_m"]))
    frames = []
    for row in candidates:
        catalog_id = str(row["catalog_id"])
        target_id = f"S23B_{catalog_id}"
        if not re.fullmatch(r"[A-Za-z0-9_-]+", target_id):
            raise ValueError(f"unsafe S23B target identifier: {target_id}")
        z = float(row["z_cl"])
        ra, dec, richness = float(row["ra_deg"]), float(row["dec_deg"]), float(row["richness"])
        half_width = HALF_WIDTH_KPC / float(cosmo.kpc_proper_per_arcmin(z).value) * 60.0
        if not (0 <= ra < 360 and -90 <= dec <= 90 and z > 0 and richness > 0
                and 30 <= half_width <= 600):
            raise ValueError(f"invalid center, redshift, richness, or cutout size for {target_id}")
        frames.append({
            "id": target_id, "catalog_id": catalog_id, "name": row["name"],
            "sample_role": "target", "center_definition": (
                "S23B CAMIRA catalog center; correspondence to algorithm-identified BCG "
                "coordinates is not independently verified"
            ),
            "ra_deg": ra, "dec_deg": dec, "z": z, "z_cl": z,
            "z_spec_bcg": None, "richness": richness,
            "requested_half_width_arcsec": half_width,
            "requested_half_width_kpc_proper": HALF_WIDTH_KPC,
            "redshift_bin": "S23B CAMIRA; continuous z retained",
            "richness_bin": "S23B CAMIRA; continuous richness retained",
            "image": None, "science_hdu": None, "sky": None,
            "product": "coadd/bg", "band": "HSC-I", "bcg_policy": "not_applicable_cluster_center",
            "members_radec_deg": [], "randoms": [], "reused_from_previous": False,
        })
    ids = [frame["id"] for frame in frames]
    if len(frames) != 1043 or len(set(ids)) != 1043 or len(frames) > MAX_SAMPLE_FRAMES:
        raise ValueError("S23B target manifest is incomplete or exceeds downloader capacity")

    estimate, estimate_rows = estimate_bytes(frames)
    if any(row["estimated_bytes"] > MAX_FILE_BYTES for row in estimate_rows):
        raise ValueError("one or more per-frame estimates exceed the downloader's 512-MiB limit")
    output_dir.mkdir(parents=True, exist_ok=False)
    manifest = {key: value for key, value in parent.items() if key != "clusters"}
    manifest.update({
        "selection": "S23B CAMIRA centers with center and eight 1-Mpc HSC-I footprint/QA probes all covered; no S20A sample match within 60 arcsec",
        "selection_limits": (
            "S23B CAMIRA candidate sample; catalog completeness, the RA/Dec-to-BCG "
            "coordinate mapping, and member selection remain unverified"
        ),
        "geometry_redshift": "z_cl",
        "catalog_sha256": EXPECTED_S23B_SHA256,
        "member_catalog_sha256": None,
        "parent_selection": str(candidate_path),
        "parent_selection_sha256": sha256(candidate_path),
        "sample_origin": "S23B CAMIRA Wide candidate expansion",
        "center_definition": (
            "CAMIRA catalog center; correspondence to the algorithm-identified BCG "
            "coordinates is not independently verified"
        ),
        "requested_outer_radius_kpc": HALF_WIDTH_KPC,
        "product": "cluster-center diffuse light",
        "science_measurement": False,
        "clusters": frames,
    })
    manifest_path = output_dir / "target_manifest.json"
    write_json(manifest_path, manifest)
    estimate_path = output_dir / "per_frame_estimate.tsv"
    with estimate_path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(estimate_rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(estimate_rows)

    summary = {
        "schema": "astrafek-s23b-icl-target-image-plan-v1",
        "status": "prepared_after_user_image_approval",
        "user_approval_received": "2026-11-07: image receipt approved for the active S23B ICL expansion",
        "approved_scope": {
            "targets": len(frames), "release": "HSC PDR3 pdr3_wide",
            "filter_family": "HSC-I", "product": "coadd/bg",
            "planes": ["IMAGE", "MASK", "VARIANCE"],
            "half_width_kpc_proper": HALF_WIDTH_KPC,
            "shape_rows": False, "photoz_rows": False,
        },
        "maximum_estimated_bytes": estimate,
        "maximum_estimated_gib": estimate / 1024**3,
        "minimum_free_space_with_30pct_margin_bytes": math.ceil(estimate * 1.3),
        "n_frames_over_512_mib_estimate": sum(row["estimated_bytes"] > MAX_FILE_BYTES for row in estimate_rows),
        "n_half_width_over_600_arcsec": sum(frame["requested_half_width_arcsec"] > 600 for frame in frames),
        "minimum_half_width_arcsec": min(frame["requested_half_width_arcsec"] for frame in frames),
        "median_half_width_arcsec": float(np.median([frame["requested_half_width_arcsec"] for frame in frames])),
        "maximum_half_width_arcsec": max(frame["requested_half_width_arcsec"] for frame in frames),
        "existing_sample_nearest_separation_min_arcsec": float(np.min(nearest_arcsec)),
        "credentials_saved": False, "science_measurement": False,
        "images_downloaded": False, "controls_included": False,
        "network_used": False,
        "inputs_sha256": {
            str(parent_manifest.resolve()): sha256(parent_manifest),
            str(s23b_catalog.resolve()): EXPECTED_S23B_SHA256,
            str(candidate_path): sha256(candidate_path),
            str(Path(plan_dir).resolve() / "control_candidate_query_plan.json"):
                sha256(Path(plan_dir).resolve() / "control_candidate_query_plan.json"),
        },
        "outputs_sha256": {
            "target_manifest.json": sha256(manifest_path),
            "per_frame_estimate.tsv": sha256(estimate_path),
        },
        "limitations": [
            "target centers use the CAMIRA catalog coordinates; their exact mapping to the algorithm-identified BCG coordinates has not been independently verified",
            "member catalog positions are not included; satellite matching remains a later pipeline stage",
            "metadata footprint QA does not establish full pixel-level WCS, mask, background, variance, or PSF quality",
            "no matched-control images are included in this target-only request",
            "coadd/bg IMAGE/MASK/VARIANCE is not raw visit imaging and its sky transfer remains under audit",
        ],
    }
    summary_path = output_dir / "summary.json"
    write_json(summary_path, summary)
    print(json.dumps({key: summary[key] for key in (
        "status", "approved_scope", "maximum_estimated_bytes", "maximum_estimated_gib",
        "minimum_free_space_with_30pct_margin_bytes", "minimum_half_width_arcsec",
        "median_half_width_arcsec", "maximum_half_width_arcsec", "outputs_sha256",
        "science_measurement")}, indent=2))
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--plan-dir", type=Path, default=DEFAULT_PLAN)
    parser.add_argument("--parent-manifest", type=Path, default=DEFAULT_PARENT_MANIFEST)
    parser.add_argument("--s23b-catalog", type=Path, default=DEFAULT_S23B_CATALOG)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    prepare(args.plan_dir, args.parent_manifest, args.s23b_catalog, args.output_dir)


if __name__ == "__main__":
    main()
