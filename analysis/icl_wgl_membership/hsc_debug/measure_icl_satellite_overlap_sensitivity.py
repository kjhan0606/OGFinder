#!/usr/bin/env python3
"""Mask satellite 50--75 kpc target pixels from the ICL ring as a sensitivity.

For the fixed 74-pair association sample, reconstruct the actual valid satellite
target-aperture pixels that overlap the 200--300 kpc ICL annulus. Recompute the
raw target-minus-control ICL brightness after excluding the union of those
pixels, then compare redshift/richness-adjusted correlations using the same
filter-specific connected spatial blocks. This is a de-overlap sensitivity,
not a replacement science estimator.
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
from audit_satellite_icl_pixel_overlap import selected_pairs
from audit_satellite_icl_shared_pixel_covariance import annulus_pixel_indices


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError("refusing to write an empty sensitivity table")
    with path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t", extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def partial_correlation(rows, indices=None):
    data = [rows[index] for index in (range(len(rows)) if indices is None else indices)]
    if len(data) < 8:
        return None
    x = np.asarray([row["icl_delta_njy_arcsec2"] for row in data], dtype=float)
    y = np.asarray([row["satellite_delta_njy_arcsec2"] for row in data], dtype=float)
    z = np.asarray([row["z"] for row in data], dtype=float)
    richness = np.asarray([row["richness"] for row in data], dtype=float)
    if np.std(z) <= 0 or np.std(richness) <= 0:
        return None
    z = (z - z.mean()) / z.std(ddof=1)
    richness = (richness - richness.mean()) / richness.std(ddof=1)
    design = np.column_stack((np.ones(len(data)), z, richness))
    if np.linalg.matrix_rank(design) < 3:
        return None
    rx = x - design @ np.linalg.lstsq(design, x, rcond=None)[0]
    ry = y - design @ np.linalg.lstsq(design, y, rcond=None)[0]
    denom = float(np.linalg.norm(rx) * np.linalg.norm(ry))
    return float(rx @ ry / denom) if denom > 0 else None


def block_bootstrap(rows, replicas, seed):
    labels = np.asarray([row["connected_5deg_block"] for row in rows], dtype=str)
    unique_blocks = np.unique(labels)
    indices_by_block = {block: np.flatnonzero(labels == block) for block in unique_blocks}
    rng = np.random.default_rng(seed)
    original, masked, difference = [], [], []
    for _ in range(replicas):
        sampled = rng.choice(unique_blocks, size=len(unique_blocks), replace=True)
        indices = np.concatenate([indices_by_block[block] for block in sampled])
        x = partial_correlation(rows, indices)
        trimmed_rows = [{**row, "icl_delta_njy_arcsec2": row["icl_delta_without_satellite_shells_njy_arcsec2"]}
                        for row in rows]
        y = partial_correlation(trimmed_rows, indices)
        if x is None or y is None:
            continue
        original.append(x)
        masked.append(y)
        difference.append(y - x)

    def interval(values):
        if not values:
            return None
        return np.quantile(np.asarray(values, dtype=float), [0.025, 0.5, 0.975]).tolist()

    return {
        "n_blocks": int(len(unique_blocks)),
        "replicas_requested": int(replicas),
        "replicas_valid": len(original),
        "original_partial_r_2p5_50_97p5": interval(original),
        "satellite_shells_masked_partial_r_2p5_50_97p5": interval(masked),
        "paired_change_masked_minus_original_2p5_50_97p5": interval(difference),
    }


def run(args):
    manifest_path = args.data / "sample_manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    mask_summary_path = args.masks / "summary.json"
    mask_records_path = args.masks / "records.tsv"
    profile_summary_path = args.profiles / "summary.json"
    profile_table_path = args.profiles / "satellite_profiles.tsv"
    icl_support_path = args.icl_support / "pair_metrics_k2p5_c40_r50_300.tsv"
    icl_annulus_table = args.icl_support / "annulus_pair_profiles.tsv"
    assoc_summary_path = args.association_dir / "summary.json"
    association_path = args.association_dir / "pixel_extinction_satellite_association.tsv"
    reblocked_path = args.reblocked_association / "association_rows.tsv"
    for path in (manifest_path, mask_summary_path, mask_records_path, profile_summary_path,
                 profile_table_path, icl_support_path, icl_annulus_table,
                 assoc_summary_path, association_path, reblocked_path):
        if not path.is_file():
            raise FileNotFoundError(path)
    mask_summary = json.loads(mask_summary_path.read_text(encoding="utf-8"))
    profile_summary = json.loads(profile_summary_path.read_text(encoding="utf-8"))
    icl_summary_path = args.icl_support / "summary.json"
    icl_summary = json.loads(icl_summary_path.read_text(encoding="utf-8"))
    assoc_summary = json.loads(assoc_summary_path.read_text(encoding="utf-8"))
    if (mask_summary.get("science_measurement") is not False
            or profile_summary.get("science_measurement") is not False
            or profile_summary.get("n_frames_completed") != 800
            or icl_summary.get("schema") != "hsc-icl-per-annulus-support-v1"
            or icl_summary.get("status") != "complete"
            or assoc_summary.get("science_measurement") is not False):
        raise ValueError("one or more exploratory source products are incomplete or have unexpected gates")
    for file_path, summary, hash_field in (
        (profile_table_path, profile_summary, "output_sha256"),
        (icl_support_path, icl_summary, "output_sha256"),
        (association_path, assoc_summary, "output_sha256"),
    ):
        if summary.get(hash_field, {}).get(str(file_path)) != sha256(file_path):
            raise ValueError(f"input table checksum mismatch: {file_path}")

    profiles_by_pair = selected_pairs(icl_support_path, profile_table_path)
    association = read_tsv(association_path)
    reblocked = read_tsv(reblocked_path)
    association_by_key = {(str(row["cluster_id"]), str(row["filter"])): row for row in association}
    block_by_key = {(str(row["target_id"]), str(row["filter"])): str(row["connected_5deg_block"])
                    for row in reblocked}
    if len(association_by_key) != 74 or len(block_by_key) != 74:
        raise ValueError("expected the fixed 74-pair ICL--satellite association sample")
    if set(association_by_key) != set(block_by_key):
        raise ValueError("association values and corrected spatial block labels differ")
    selected_pair_keys = {
        (ident, str(row["filter"]))
        for ident, rows in profiles_by_pair.items() for row in rows
        if (ident, str(row["filter"])) in association_by_key
    }
    if selected_pair_keys != set(association_by_key):
        raise ValueError("the fixed association sample is not fully covered by satellite profile support")

    support_rows = read_tsv(icl_support_path)
    annulus_rows = read_tsv(icl_annulus_table)
    support_by_key = {(str(row["id"]), str(row["filter"])): row for row in support_rows}
    annulus_by_key = {(str(row["id"]), str(row["filter"])): row for row in annulus_rows
                      if row["radial_bin_kpc"] == "200-300"}
    frame_by_id = {str(row["id"]): row for row in manifest["clusters"]}
    mask_records = read_tsv(mask_records_path)
    mask_by_key = {(str(row["id"]), str(row["primary_filter"])): row for row in mask_records}
    all_satellite_rows = read_tsv(profile_table_path)
    selected_by_key = {}
    for ident, selected in profiles_by_pair.items():
        for row in selected:
            key = (ident, str(row["filter"]))
            if key in association_by_key:
                selected_by_key.setdefault(key, []).append(row)
    if sum(map(len, selected_by_key.values())) != 206:
        raise ValueError("expected 206 selected satellite profiles in the fixed association sample")

    output_rows = []
    for index, key in enumerate(sorted(association_by_key), 1):
        ident, filt = key
        print(f"[{index}/74] exclude satellite-shell pixels from ICL ring {ident} {filt}", flush=True)
        frame = frame_by_id[ident]
        mask_record = mask_by_key.get(key)
        pair_icl = support_by_key.get(key)
        annulus_measurement = annulus_by_key.get(key)
        if mask_record is None or pair_icl is None or annulus_measurement is None:
            raise ValueError(f"missing frame, mask, or ICL annulus support for {key}")
        image_path, icl_mask_path = Path(mask_record["image"]), Path(mask_record["mask_file"])
        if sha256(image_path) != mask_record["image_sha256"]:
            raise ValueError(f"target image hash changed for {key}")
        if sha256(icl_mask_path) != mask_record["mask_sha256"]:
            raise ValueError(f"ICL mask hash changed for {key}")
        prep = args.satellite_prep / f"cluster_{ident}"
        detections_list = satellite.parse_detections(satellite.read_tsv(prep / "detections.tsv"))
        member_rows = satellite.read_tsv(prep / "members.tsv")
        preparation = satellite.read_json(prep / "preparation.json")
        if preparation.get("science_measurement") is not False:
            raise ValueError(f"satellite preparation has an unexpected science gate for {key}")
        info = satellite.validate_fits(image_path, frame)
        if info["coadd_header_filter"] != filt:
            raise ValueError(f"target image filter mismatch for {key}")
        matched, _ = satellite.match_members(member_rows, detections_list, info["pixel_scale_arcsec"])
        matched_by_index = {item["member_index"]: item for item in matched}
        member_by_index = {int(row["INDEX"]): row for row in member_rows if row["IN_FRAME"] == "True"}
        matched_indices = set(matched_by_index)
        member_xy = np.asarray([[float(row["X_IMAGE"])-1.0, float(row["Y_IMAGE"])-1.0]
                                for row in member_by_index.values()], dtype=float).reshape((-1, 2))
        fallback_xy = np.asarray([[float(row["X_IMAGE"])-1.0, float(row["Y_IMAGE"])-1.0]
                                  for member_index, row in member_by_index.items()
                                  if member_index not in matched_indices], dtype=float).reshape((-1, 2))
        fallback_tree = cKDTree(fallback_xy) if len(fallback_xy) else None
        pixel_scale = float(info["pixel_scale_arcsec"])
        center_xy = np.asarray(info["center_pixel_xy"], dtype=float)
        kpc_per_pixel = float(FlatLambdaCDM(
            H0=float(manifest.get("h0", 70.0)), Om0=float(manifest.get("omega_m", 0.3)),
            Tcmb0=0*u.K).kpc_proper_per_arcmin(float(frame["z"])).value * pixel_scale / 60.0)

        with fits.open(image_path, memmap=True, checksum=True) as hdul:
            image = np.asarray(hdul[info["science_hdu"]].data)
            variance = np.asarray(hdul[info["variance_hdu"]].data)
        with fits.open(icl_mask_path, memmap=True, checksum=True) as hdul:
            icl_excluded = np.asarray(hdul["MASK_K2P5"].data, dtype=bool)
        with fits.open(prep / "masks.fits.gz", memmap=True, checksum=True) as hdul:
            segmentation = np.asarray(hdul["SEGMENTATION"].data, dtype=np.int32)
            base = np.asarray(hdul["BASE_MASK"].data, dtype=bool)
        if not (image.shape == variance.shape == icl_excluded.shape == segmentation.shape == base.shape):
            raise ValueError(f"image/mask shape mismatch for {key}")
        finite = np.isfinite(image) & np.isfinite(variance) & (variance > 0)
        if np.any(~finite & ~base):
            raise ValueError(f"base mask misses invalid pixels for {key}")
        yy, xx = np.ogrid[:image.shape[0], :image.shape[1]]
        cluster_radius = np.hypot((xx-center_xy[0])*kpc_per_pixel,
                                  (yy-center_xy[1])*kpc_per_pixel)
        geometric = (cluster_radius >= 200.0) & (cluster_radius < 300.0)
        icl_valid = ~icl_excluded & finite
        full_ids = np.flatnonzero(geometric & icl_valid)
        union = np.zeros(image.shape, dtype=bool)
        detections = satellite.detection_arrays(detections_list, pixel_scale)
        detections_by_number = {int(row["NUMBER"]): row for row in detections_list}
        selected_profiles = selected_by_key[key]
        for profile in selected_profiles:
            member_index = int(profile["member_index"])
            number = int(profile["detection_number"])
            match = matched_by_index.get(member_index)
            if match is None or int(detections_list[match["detection_index"]]["NUMBER"]) != number:
                raise ValueError(f"member/detection mapping mismatch for {key}/{member_index}")
            target = detections_by_number.get(number)
            if target is None:
                raise ValueError(f"profile detection missing for {key}/{member_index}")
            target_xy = np.asarray([float(target["x"]), float(target["y"])])
            pixel_ids = annulus_pixel_indices(
                target_xy, number, image.shape, base, segmentation, detections,
                fallback_tree, pixel_scale, kpc_per_pixel, finite,
            )
            if len(pixel_ids) != int(profile["target_n_pixels"]):
                raise ValueError(f"satellite target support changed for {key}/{member_index}")
            in_icl = geometric.ravel()[pixel_ids] & icl_valid.ravel()[pixel_ids]
            union.ravel()[pixel_ids[in_icl]] = True

        overlap_ids = np.flatnonzero(union)
        retained_ids = full_ids[~union.ravel()[full_ids]]
        if len(retained_ids) == 0:
            raise ValueError(f"satellite-shell exclusion removes the entire ICL annulus for {key}")
        conversion = float(satellite.nJy_factor(info))
        full_target = float(np.mean(image.ravel()[full_ids], dtype=np.float64) * conversion)
        trimmed_target = float(np.mean(image.ravel()[retained_ids], dtype=np.float64) * conversion)
        control_sb = float(annulus_measurement["control_sb_njy_arcsec2"])
        if not math.isclose(full_target, float(annulus_measurement["target_sb_njy_arcsec2"]),
                            rel_tol=0, abs_tol=1e-5):
            raise ValueError(f"reconstructed original ICL target mean differs for {key}")
        original_delta = full_target - control_sb
        trimmed_delta = trimmed_target - control_sb
        recorded_delta = float(association_by_key[key]["icl_delta_raw_njy_arcsec2"])
        if not math.isclose(original_delta, recorded_delta, rel_tol=0, abs_tol=1e-5):
            raise ValueError(f"reconstructed original ICL pair does not match association for {key}")
        original_coverage = len(full_ids) / int(np.count_nonzero(geometric))
        trimmed_coverage = len(retained_ids) / int(np.count_nonzero(geometric))
        association_row = association_by_key[key]
        output_rows.append({
            "cluster_id": ident, "filter": filt,
            "icl_delta_njy_arcsec2": recorded_delta,
            "icl_delta_raw_njy_arcsec2": recorded_delta,
            "icl_delta_without_satellite_shells_njy_arcsec2": trimmed_delta,
            "delta_icl_njy_arcsec2": trimmed_delta - recorded_delta,
            "satellite_delta_njy_arcsec2": float(association_row["satellite_delta_raw_njy_arcsec2"]),
            "satellite_delta_raw_njy_arcsec2": float(association_row["satellite_delta_raw_njy_arcsec2"]),
            "z": float(association_row["z"]), "richness": float(association_row["richness"]),
            "connected_5deg_block": block_by_key[key],
            "icl_annulus_geometric_pixels": int(np.count_nonzero(geometric)),
            "icl_pixels_before": int(len(full_ids)),
            "icl_pixels_removed_by_satellite_shell_union": int(len(overlap_ids)),
            "icl_pixels_after": int(len(retained_ids)),
            "icl_coverage_before": float(original_coverage),
            "icl_coverage_after": float(trimmed_coverage),
            "n_satellite_profiles": int(len(selected_profiles)),
        })

    correlations = {}
    for band_index, band in enumerate(("i", "i2")):
        rows = [row for row in output_rows if row["filter"] == band]
        original = partial_correlation(rows)
        trimmed_rows = [{**row, "icl_delta_njy_arcsec2": row["icl_delta_without_satellite_shells_njy_arcsec2"]}
                        for row in rows]
        trimmed = partial_correlation(trimmed_rows)
        if original is None or trimmed is None:
            raise ValueError(f"partial correlation failed for {band}")
        bootstrap = block_bootstrap(rows, args.replicas, args.seed + band_index)
        correlations[band] = {
            "n_pairs": len(rows),
            "n_pairs_remaining_at_icl_coverage_floor_0p4": sum(row["icl_coverage_after"] >= 0.4 for row in rows),
            "original_partial_r_z_richness_adjusted": original,
            "satellite_shells_masked_partial_r_z_richness_adjusted": trimmed,
            "point_change_masked_minus_original": trimmed - original,
            "block_bootstrap": bootstrap,
            "median_removed_icl_support_fraction": float(np.median([
                row["icl_pixels_removed_by_satellite_shell_union"] / row["icl_pixels_before"]
                for row in rows
            ])),
            "median_delta_icl_njy_arcsec2": float(np.median([row["delta_icl_njy_arcsec2"] for row in rows])),
        }
        retained = [row for row in rows if row["icl_coverage_after"] >= 0.4]
        if len(retained) >= 8:
            retained_trimmed = [{**row, "icl_delta_njy_arcsec2": row["icl_delta_without_satellite_shells_njy_arcsec2"]}
                                for row in retained]
            correlations[band]["same_coverage_floor_subset"] = {
                "n_pairs": len(retained),
                "original_partial_r": partial_correlation(retained),
                "satellite_shells_masked_partial_r": partial_correlation(retained_trimmed),
                "block_bootstrap": block_bootstrap(retained, args.replicas, args.seed + 100 + band_index),
            }
        else:
            correlations[band]["same_coverage_floor_subset"] = {
                "n_pairs": len(retained), "status": "too_few_pairs_for_partial_correlation",
            }

    return {
        "schema": "astrafek-icl-satellite-shell-overlap-sensitivity-v1",
        "status": "exploratory",
        "science_measurement": False,
        "method": "recompute ICL target-minus-control 200-300 kpc raw surface brightness after excluding the union of valid satellite target 50-75 kpc annulus pixels from the target ICL ring; keep the matched-control ICL value and satellite profile fixed",
        "n_association_pairs": len(output_rows),
        "n_satellite_profiles": sum(row["n_satellite_profiles"] for row in output_rows),
        "profiles_by_filter": {
            band: sum(row["n_satellite_profiles"] for row in output_rows if row["filter"] == band)
            for band in ("i", "i2")
        },
        "correlations": correlations,
        "replicas": args.replicas,
        "seed": args.seed,
        "limitations": [
            "this is a fixed-sample aperture-exclusion sensitivity, not a PSF/deblending correction",
            "only supported satellite target annuli are excluded; satellite light outside 50-75 kpc and mask leakage remain",
            "the ICL annulus support changes locally after masking; matched-control and all satellite measurements remain fixed",
            "selected sample and 5-degree spatial blocks are exploratory and do not establish independent HSC fields",
            "large-scale sky covariance, extended PSF wings, extinction/bandpass, completeness, and calibration errors remain unpropagated",
            "no HSC query or download was performed",
        ],
        "input_sha256": {
            str(path.resolve()): sha256(path) for path in (
                manifest_path, mask_summary_path, mask_records_path, profile_summary_path,
                profile_table_path, icl_summary_path, icl_support_path, icl_annulus_table,
                assoc_summary_path, association_path, reblocked_path, Path(__file__).resolve(),
            )
        },
        "output_rows": output_rows,
    }


def main():
    root = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=root)
    parser.add_argument("--masks", type=Path, default=root / "prestack_cluster_masks_2026-10-16_01")
    parser.add_argument("--satellite-prep", type=Path, default=root / "satellite_detection_psf_2026-10-18_01")
    parser.add_argument("--profiles", type=Path, default=root / "satellite_outer_profile_perbin_full_2026-10-18_01")
    parser.add_argument("--icl-support", type=Path, default=root / "icl_annulus_support_2026-10-18_01")
    parser.add_argument("--association-dir", type=Path, default=root / "icl_pixel_extinction_2026-10-19_03")
    parser.add_argument("--reblocked-association", type=Path,
                        default=root / "icl_outer_sky_association_reblocked_2026-11-03_01")
    parser.add_argument("--replicas", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=20261105)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    if args.replicas < 1000:
        parser.error("at least 1000 spatial block resamples are required")
    if args.output_dir.exists():
        parser.error(f"refusing to overwrite {args.output_dir}")
    result = run(args)
    args.output_dir.mkdir(parents=True, exist_ok=False)
    write_tsv(args.output_dir / "pair_sensitivity.tsv", result.pop("output_rows"))
    summary_path = args.output_dir / "summary.json"
    result["pair_table_sha256"] = sha256(args.output_dir / "pair_sensitivity.tsv")
    summary_path.write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n",
                            encoding="utf-8")
    print(json.dumps({
        "status": result["status"],
        "science_measurement": result["science_measurement"],
        "n_pairs": result["n_association_pairs"],
        "correlations": result["correlations"],
        "output_dir": str(args.output_dir.resolve()),
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
