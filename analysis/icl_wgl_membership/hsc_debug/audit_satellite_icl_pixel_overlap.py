#!/usr/bin/env python3
"""Measure actual shared valid pixels in ICL and satellite annuli.

This is a support audit on already downloaded HSC images. It applies the saved
K2.5 ICL exclusion mask and the same SEP/Kron satellite mask used by the
satellite profile estimator. It does not query HSC or measure flux covariance.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from astropy.cosmology import FlatLambdaCDM
from astropy.io import fits
from astropy import units as u
from scipy.spatial import cKDTree

import measure_satellite_outer_profiles as satellite

ROOT = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
DEFAULT_MASKS = ROOT / "prestack_cluster_masks_2026-10-16_01"
DEFAULT_PREP = ROOT / "satellite_detection_psf_2026-10-18_01"
DEFAULT_PROFILES = ROOT / "satellite_outer_profile_perbin_full_2026-10-18_01"
DEFAULT_ICL = ROOT / "prestack_flux_profiles_2026-10-18_02/pair_metrics_k2p5_c40_r50_300.tsv"
ICL_INNER, ICL_OUTER = 200.0, 300.0
SAT_INNER, SAT_OUTER = 50.0, 75.0
ICL_COVERAGE_FLOOR = 0.40
SAT_COVERAGE_FLOOR = 0.60


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with Path(path).open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def selected_pairs(icl_pairs_path: Path, satellite_profiles_path: Path) -> dict[str, list[dict[str, str]]]:
    icl_rows = read_tsv(icl_pairs_path)
    sat_rows = read_tsv(satellite_profiles_path)
    if not icl_rows or not {"id", "filter", "coverage_200_300"}.issubset(icl_rows[0]):
        raise ValueError("ICL table lacks required pair-support fields")
    required = {"cluster_id", "filter", "member_index", "detection_number",
                "radial_bin_kpc", "target_coverage", "control_coverage", "n_controls"}
    if not sat_rows or not required.issubset(sat_rows[0]):
        raise ValueError("satellite table lacks required per-annulus fields")
    icl_support = {(row["id"], row["filter"])
                   for row in icl_rows
                   if float(row["coverage_200_300"]) >= ICL_COVERAGE_FLOOR}
    selected: dict[str, list[dict[str, str]]] = {}
    for row in sat_rows:
        if (row["cohort"] != "full_sample" or row["radial_bin_kpc"] != "50-75"
                or float(row["target_coverage"]) < SAT_COVERAGE_FLOOR
                or float(row["control_coverage"]) < SAT_COVERAGE_FLOOR
                or int(row["n_controls"]) < 3
                or (row["cluster_id"], row["filter"]) not in icl_support):
            continue
        selected.setdefault(row["cluster_id"], []).append(row)
    return selected


def shared_pixel_counts(satellite_annulus: np.ndarray, icl_annulus: np.ndarray,
                        satellite_valid: np.ndarray, icl_valid: np.ndarray) -> dict[str, int]:
    """Count geometrical and estimator-valid support in aligned image crops."""
    arrays = (satellite_annulus, icl_annulus, satellite_valid, icl_valid)
    if any(array.shape != satellite_annulus.shape for array in arrays):
        raise ValueError("overlap masks must have identical shapes")
    geometric = satellite_annulus & icl_annulus
    sat_support = satellite_annulus & satellite_valid
    icl_support = icl_annulus & icl_valid
    common = geometric & satellite_valid & icl_valid
    return {
        "n_satellite_annulus_geometric_pixels": int(satellite_annulus.sum()),
        "n_satellite_annulus_valid_pixels": int(sat_support.sum()),
        "n_icl_annulus_valid_pixels_local": int(icl_support.sum()),
        "n_geometric_intersection_pixels": int(geometric.sum()),
        "n_shared_valid_measurement_pixels": int(common.sum()),
    }


def process_cluster(ident: str, target_rows: list[dict[str, str]], prep_dir: Path,
                    manifest: dict, mask_summary: dict, mask_records: list[dict[str, str]],
                    icl_pairs_by_key: dict[tuple[str, str], dict[str, str]]) -> list[dict]:
    mask_matches = [row for row in mask_records if row["id"] == ident]
    if len(mask_matches) != 1:
        raise ValueError(f"cluster {ident}: expected one saved ICL mask record")
    mask_record = mask_matches[0]
    filt = mask_record["primary_filter"]
    target_rows = [row for row in target_rows if row["filter"] == filt]
    if not target_rows:
        return []
    frame_matches = [row for row in manifest["clusters"] if str(row["id"]) == ident]
    if len(frame_matches) != 1:
        raise ValueError(f"cluster {ident}: manifest identity is not unique")
    frame = dict(frame_matches[0])
    image_path = Path(mask_record["image"])
    frame["image"] = str(image_path)
    icl_mask_path = Path(mask_record["mask_file"])
    if sha256(image_path) != mask_record["image_sha256"]:
        raise ValueError(f"cluster {ident}: image SHA-256 mismatch")
    if (sha256(icl_mask_path) != mask_record["mask_sha256"]
            or mask_summary.get("mask_output_sha256", {}).get(str(icl_mask_path))
            != mask_record["mask_sha256"]):
        raise ValueError(f"cluster {ident}: ICL mask SHA-256 mismatch")
    prep = prep_dir / f"cluster_{ident}"
    detection_path, member_path = prep / "detections.tsv", prep / "members.tsv"
    sat_mask_path, prep_summary_path = prep / "masks.fits.gz", prep / "preparation.json"
    preparation = satellite.read_json(prep_summary_path)
    if preparation.get("science_measurement") is not False:
        raise ValueError(f"cluster {ident}: unexpected preparation science gate")
    info = satellite.validate_fits(image_path, frame)
    if info["coadd_header_filter"] != filt:
        raise ValueError(f"cluster {ident}: image/mask filter mismatch")
    detections_list = satellite.parse_detections(satellite.read_tsv(detection_path))
    member_rows = satellite.read_tsv(member_path)
    matched, _ = satellite.match_members(member_rows, detections_list, info["pixel_scale_arcsec"])
    with fits.open(image_path, memmap=False, checksum=True) as hdul:
        image = np.asarray(hdul[info["science_hdu"]].data, dtype=np.float32)
        variance = np.asarray(hdul[info["variance_hdu"]].data, dtype=np.float32)
    with fits.open(icl_mask_path, memmap=False, checksum=True) as hdul:
        icl_excluded = np.asarray(hdul["MASK_K2P5"].data, dtype=bool)
    with fits.open(sat_mask_path, memmap=False, checksum=True) as hdul:
        segmentation = np.asarray(hdul["SEGMENTATION"].data, dtype=np.int32)
        base = np.asarray(hdul["BASE_MASK"].data, dtype=bool)
    if not (image.shape == variance.shape == icl_excluded.shape == segmentation.shape == base.shape):
        raise ValueError(f"cluster {ident}: image and mask shapes differ")
    finite = np.isfinite(image) & np.isfinite(variance) & (variance > 0)
    if np.any(~finite & ~base):
        raise ValueError(f"cluster {ident}: BASE_MASK misses invalid image/variance pixels")

    pixel_scale = float(info["pixel_scale_arcsec"])
    center_xy = np.asarray(info["center_pixel_xy"], dtype=float)
    kpc_per_pixel = float(FlatLambdaCDM(
        H0=float(manifest.get("h0", 70.)), Om0=float(manifest.get("omega_m", .3)),
        Tcmb0=0*u.K).kpc_proper_per_arcmin(float(frame["z"])).value * pixel_scale / 60.)
    detections = satellite.detection_arrays(detections_list, pixel_scale)
    matched_indices = {item["member_index"] for item in matched}
    in_frame_members = [row for row in member_rows if row["IN_FRAME"] == "True"]
    member_xy = np.asarray([[float(row["X_IMAGE"])-1., float(row["Y_IMAGE"])-1.]
                            for row in in_frame_members], dtype=float).reshape((-1, 2))
    fallback_xy = np.asarray([[float(row["X_IMAGE"])-1., float(row["Y_IMAGE"])-1.]
                              for row in in_frame_members
                              if int(row["INDEX"]) not in matched_indices], dtype=float).reshape((-1, 2))
    fallback_tree = cKDTree(fallback_xy) if len(fallback_xy) else None
    yy, xx = np.ogrid[:image.shape[0], :image.shape[1]]
    cluster_radius = np.hypot((xx-center_xy[0])*kpc_per_pixel,
                              (yy-center_xy[1])*kpc_per_pixel)
    icl_valid = ~icl_excluded & finite
    n_icl_valid = int(np.count_nonzero((cluster_radius >= ICL_INNER)
                                       & (cluster_radius < ICL_OUTER) & icl_valid))

    target_by_number = {det["NUMBER"]: det for det in detections_list}
    member_by_index = {item["member_index"]: item for item in matched}
    output = []
    for row in target_rows:
        number = int(row["detection_number"])
        member_index = int(row["member_index"])
        if member_index not in member_by_index or number not in target_by_number:
            raise ValueError(f"cluster {ident}: profile target is absent from saved match tables")
        item = member_by_index[member_index]
        if int(detections_list[item["detection_index"]]["NUMBER"]) != number:
            raise ValueError(f"cluster {ident}: member-to-detection identity mismatch")
        target = target_by_number[number]
        cx, cy = float(target["x"]), float(target["y"])
        half_size = int(math.ceil(SAT_OUTER / kpc_per_pixel)) + 2
        x0 = max(0, int(math.floor(cx-half_size)))
        x1 = min(image.shape[1], int(math.ceil(cx+half_size))+1)
        y0 = max(0, int(math.floor(cy-half_size)))
        y1 = min(image.shape[0], int(math.ceil(cy+half_size))+1)
        sat_mask = satellite.make_stamp_mask(
            base, segmentation, detections, fallback_tree, (cx, cy), x0, y0,
            number, pixel_scale, half_size, satellite.KRON_FACTOR)
        sy, sx = np.ogrid[y0:y1, x0:x1]
        sat_radius = np.hypot((sx-cx)*kpc_per_pixel, (sy-cy)*kpc_per_pixel)
        annulus = (sat_radius >= SAT_INNER) & (sat_radius < SAT_OUTER)
        cluster_slice = cluster_radius[y0:y1, x0:x1]
        icl_annulus = (cluster_slice >= ICL_INNER) & (cluster_slice < ICL_OUTER)
        finite_slice = finite[y0:y1, x0:x1]
        sat_valid = ~sat_mask & finite_slice
        icl_valid_local = icl_valid[y0:y1, x0:x1]
        geometric = annulus & icl_annulus
        counts = shared_pixel_counts(annulus, icl_annulus, sat_valid, icl_valid_local)
        n_sat_geom = counts["n_satellite_annulus_geometric_pixels"]
        n_sat_valid = counts["n_satellite_annulus_valid_pixels"]
        n_geometric = counts["n_geometric_intersection_pixels"]
        n_common = counts["n_shared_valid_measurement_pixels"]
        pair = icl_pairs_by_key[(ident, filt)]
        output.append({
            "cluster_id": ident, "filter": filt, "member_index": member_index,
            "detection_number": number, "cluster_radius_kpc": float(row["cluster_radius_kpc"]),
            **counts,
            "n_icl_annulus_valid_pixels_global": n_icl_valid,
            "geometric_overlap_fraction_of_satellite_annulus": n_geometric/max(1, n_sat_geom),
            "shared_valid_fraction_of_satellite_valid": n_common/max(1, n_sat_valid),
            "shared_valid_fraction_of_geometric_overlap": n_common/max(1, n_geometric),
            "target_coverage_recorded": float(row["target_coverage"]),
            "icl_coverage_recorded": float(pair["coverage_200_300"]),
        })
    return output


def run(data_root: Path, masks_dir: Path, prep_dir: Path, profile_dir: Path,
        icl_pairs_path: Path, output_dir: Path, cluster_ids: list[str] | None = None) -> dict:
    if output_dir.exists():
        raise FileExistsError(f"refusing to overwrite {output_dir}")
    manifest_path = data_root / "sample_manifest.json"
    mask_summary_path, mask_records_path = masks_dir / "summary.json", masks_dir / "records.tsv"
    profile_summary_path = profile_dir / "summary.json"
    satellite_profiles_path = profile_dir / "satellite_profiles.tsv"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    mask_summary = json.loads(mask_summary_path.read_text(encoding="utf-8"))
    profile_summary = json.loads(profile_summary_path.read_text(encoding="utf-8"))
    if (mask_summary.get("science_measurement") is not False
            or profile_summary.get("science_measurement") is not False
            or profile_summary.get("n_frames_completed") != 800):
        raise ValueError("expected complete exploratory 800-frame mask/profile inputs")
    grouped = selected_pairs(icl_pairs_path, satellite_profiles_path)
    icl_pairs_by_key = {(row["id"], row["filter"]): row for row in read_tsv(icl_pairs_path)}
    if cluster_ids:
        missing = sorted(set(cluster_ids)-set(grouped))
        if missing:
            raise ValueError(f"requested clusters are outside the matched support: {missing}")
        grouped = {ident: grouped[ident] for ident in cluster_ids}
    mask_records = read_tsv(mask_records_path)
    records = []
    for index, (ident, rows) in enumerate(sorted(grouped.items()), 1):
        print(f"[{index}/{len(grouped)}] pixel overlap cluster {ident} ({len(rows)} satellites)",
              flush=True)
        records.extend(process_cluster(
            ident, rows, prep_dir, manifest, mask_summary, mask_records,
            icl_pairs_by_key))
    output_dir.mkdir(parents=True, exist_ok=False)
    table_path = output_dir / "pixel_overlap.tsv"
    if records:
        with table_path.open("x", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(stream, fieldnames=list(records[0]), delimiter="\t")
            writer.writeheader(); writer.writerows(records)
    by_filter = {}
    for filt in sorted({row["filter"] for row in records}):
        subset = [row for row in records if row["filter"] == filt]
        fractions = [row["shared_valid_fraction_of_satellite_valid"] for row in subset]
        by_filter[filt] = {
            "n_satellite_profiles": len(subset),
            "n_clusters": len({row["cluster_id"] for row in subset}),
            "n_with_geometric_overlap": sum(row["n_geometric_intersection_pixels"] > 0 for row in subset),
            "n_with_shared_valid_pixels": sum(row["n_shared_valid_measurement_pixels"] > 0 for row in subset),
            "total_shared_valid_pixels_across_profiles": sum(row["n_shared_valid_measurement_pixels"] for row in subset),
            "median_shared_fraction_of_satellite_valid": float(np.median(fractions)),
            "q90_shared_fraction_of_satellite_valid": float(np.quantile(fractions, .9)),
        }
    summary = {
        "schema": "hsc-satellite-icl-pixel-overlap-v1",
        "status": "complete_exploratory",
        "science_measurement": False,
        "flux_covariance_available": False,
        "selection": {"icl_annulus_kpc": [ICL_INNER, ICL_OUTER],
                      "icl_coverage_floor": ICL_COVERAGE_FLOOR,
                      "satellite_annulus_kpc": [SAT_INNER, SAT_OUTER],
                      "satellite_target_and_control_coverage_floor": SAT_COVERAGE_FLOOR,
                      "n_clusters": len({row["cluster_id"] for row in records}),
                      "n_satellite_profiles": len(records)},
        "method": "Counts the pixel intersection of ICL annulus pixels retained by MASK_K2P5 and satellite annulus pixels retained by the saved BASE_MASK, SEGMENTATION, SEP/Kron neighbor masks and fallback-member masks; finite positive-variance pixels only.",
        "by_filter": by_filter,
        "input_sha256": {
            str(manifest_path): sha256(manifest_path),
            str(mask_summary_path): sha256(mask_summary_path),
            str(mask_records_path): sha256(mask_records_path),
            str(profile_summary_path): sha256(profile_summary_path),
            str(satellite_profiles_path): sha256(satellite_profiles_path),
            str(icl_pairs_path): sha256(icl_pairs_path),
            str(Path(__file__).resolve()): sha256(Path(__file__).resolve()),
        },
        "output_sha256": {str(table_path): sha256(table_path)} if records else {},
        "limitations": [
            "Shared pixel support is not a flux covariance estimate and does not model PSF wings or upstream sky subtraction.",
            "The ICL and satellite estimators use different masks and controls; this audit quantifies shared valid pixels only in each satellite target annulus, not all source/control covariance.",
            "The coverage cuts and masks are exploratory and do not certify a science sample.",
        ],
    }
    (output_dir / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True, allow_nan=False)+"\n", encoding="utf-8")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=ROOT)
    parser.add_argument("--icl-masks", type=Path, default=DEFAULT_MASKS)
    parser.add_argument("--satellite-prep", type=Path, default=DEFAULT_PREP)
    parser.add_argument("--satellite-profiles", type=Path, default=DEFAULT_PROFILES)
    parser.add_argument("--icl-pairs", type=Path, default=DEFAULT_ICL)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--cluster-id", action="append", default=[])
    args = parser.parse_args()
    summary = run(args.data.resolve(), args.icl_masks.resolve(), args.satellite_prep.resolve(),
                  args.satellite_profiles.resolve(), args.icl_pairs.resolve(),
                  args.output_dir.resolve(), args.cluster_id or None)
    print(json.dumps({"status": summary["status"], "science_measurement": False,
                      "selection": summary["selection"], "by_filter": summary["by_filter"],
                      "output": str(args.output_dir.resolve())}, indent=2))


if __name__ == "__main__":
    main()
