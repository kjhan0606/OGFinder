#!/usr/bin/env python3
"""Separate image-edge, base-quality-mask, and source-mask losses at 50–75 kpc.

The audit uses the 40%-floor satellite candidate list, existing HSC image
headers, and saved preparation masks. It reads no science-image pixels and
does not contact HSC or change the measurement.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from astropy.cosmology import FlatLambdaCDM
import astropy.units as u
from astropy.io import fits
from scipy.spatial import cKDTree

import sys
HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
import measure_satellite_outer_profiles as core
from download_pilot import validate_fits


DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
MANIFEST = DATA / "download_01/download_manifest.json"
PREP = DATA / "satellite_detection_psf_2026-10-18_01"
PREP_SUMMARY = PREP / "summary.json"
PREP_AUDIT = PREP / "artifact_audit.json"
ICL_SUPPORT = DATA / "icl_annulus_support_2026-10-18_01/pair_metrics_k2p5_c40_r50_300.tsv"
SATELLITE_40 = DATA / "satellite_outer_profile_perbin_cov040_2026-11-05_01"
OUTER_BIN = "50-75"
COVERAGE_FLOOR = 0.40
MIN_CONTROLS = 3


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty table: {path}")
    with path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t",
                                extrasaction="raise", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def coverage_parts(base: np.ndarray, segmentation: np.ndarray, detections: dict,
                   fallback_tree: cKDTree | None, center_xy: tuple[float, float],
                   target_number: int, pixel_scale: float, kpc_per_pixel: float) -> dict:
    cx, cy = map(float, center_xy)
    lo, hi = map(float, core.R_EDGES_KPC[-2:])
    half = int(math.ceil(hi / kpc_per_pixel)) + 2
    x0, x1 = max(0, int(math.floor(cx-half))), min(base.shape[1], int(math.ceil(cx+half))+1)
    y0, y1 = max(0, int(math.floor(cy-half))), min(base.shape[0], int(math.ceil(cy+half))+1)
    if x0 >= x1 or y0 >= y1:
        return {"footprint": 0., "after_base_mask": 0., "after_source_masks": 0.}
    # Match profile_at's float32 pixel-center arithmetic exactly. Using float64
    # here can move boundary pixels between radial bins and disagree with the
    # saved three-decimal coverage even when the mask is identical.
    yy, xx = np.indices((y1-y0, x1-x0), dtype=np.float32)
    radius = np.hypot(xx+x0-cx, yy+y0-cy) * kpc_per_pixel
    annulus = (radius >= lo) & (radius < hi)
    expected = math.pi * (hi**2-lo**2) / kpc_per_pixel**2
    geometric_n = int(annulus.sum())
    mask = np.array(base[y0:y1, x0:x1], dtype=bool, copy=True)
    base_n = int((annulus & ~mask).sum())

    # Apply each source-mask layer in the same order as make_stamp_mask so the
    # marginal area loss can identify which prepared mask drives low coverage.
    labels = segmentation[y0:y1, x0:x1]
    mask |= (labels > 0) & (labels != target_number)
    after_segmentation_n = int((annulus & ~mask).sum())

    xs, ys = detections["x"], detections["y"]
    near = ((xs + detections["extent"] >= x0) & (xs - detections["extent"] < x1)
            & (ys + detections["extent"] >= y0) & (ys - detections["extent"] < y1)
            & (detections["number"] != target_number))
    selected = np.flatnonzero(near)
    valid_selected = selected[detections["shape_valid"][selected]]
    invalid_selected = selected[~detections["shape_valid"][selected]]
    if len(valid_selected):
        aa = np.maximum(detections["a"][valid_selected], .5)
        bb = np.maximum(detections["b"][valid_selected], .5)
        core.sep.mask_ellipse(mask, xs[valid_selected] - x0, ys[valid_selected] - y0,
                              aa, bb, detections["theta"][valid_selected],
                              r=detections["kron_scale"][valid_selected])
    fallback_radius = max(2. / pixel_scale, 1.5)
    for index in invalid_selected:
        core.add_circle(mask, float(xs[index] - x0), float(ys[index] - y0),
                        fallback_radius)
    after_detections_n = int((annulus & ~mask).sum())

    if fallback_tree is not None:
        for index in fallback_tree.query_ball_point(center_xy, half + fallback_radius):
            point = fallback_tree.data[index]
            core.add_circle(mask, float(point[0] - x0), float(point[1] - y0),
                            fallback_radius)
    full_n = int((annulus & ~mask).sum())
    return {
        "footprint": min(1., geometric_n/expected),
        "after_base_mask": min(1., base_n/expected),
        "after_segmentation_mask": min(1., after_segmentation_n/expected),
        "after_detection_neighbor_masks": min(1., after_detections_n/expected),
        "after_source_masks": min(1., full_n/expected),
        "footprint_pixels": geometric_n,
        "base_valid_pixels": base_n,
        "segmentation_valid_pixels": after_segmentation_n,
        "detection_masks_valid_pixels": after_detections_n,
        "source_mask_valid_pixels": full_n,
    }


def run(output: Path) -> dict:
    for path in (MANIFEST, PREP_SUMMARY, PREP_AUDIT, ICL_SUPPORT,
                 SATELLITE_40/"summary.json", SATELLITE_40/"candidate_diagnostics.tsv",
                 SATELLITE_40/"satellite_profiles.tsv"):
        if not path.is_file():
            raise FileNotFoundError(path)
    if output.exists():
        raise FileExistsError(f"refusing to overwrite {output}")
    manifest = json.loads(MANIFEST.read_text(encoding="utf-8"))
    prep_summary = json.loads(PREP_SUMMARY.read_text(encoding="utf-8"))
    prep_audit = json.loads(PREP_AUDIT.read_text(encoding="utf-8"))
    low_summary = json.loads((SATELLITE_40/"summary.json").read_text(encoding="utf-8"))
    if (prep_summary.get("status") != "complete" or prep_summary.get("n_prepared") != 800
            or prep_audit.get("passed") is not True or prep_audit.get("science_measurement") is not False
            or low_summary.get("status") != "complete" or low_summary.get("n_frames_completed") != 800
            or low_summary.get("method", {}).get("coverage_floor") != COVERAGE_FLOOR
            or low_summary.get("method", {}).get("minimum_controls_per_annulus") != MIN_CONTROLS):
        raise ValueError("expected verified 800-frame preparation and completed 40% diagnostic")

    frames = core.resolve_image_paths(manifest["clusters"], MANIFEST)
    frame_by_key = {(str(frame["id"]), frame["coadd_header_filter"]): frame for frame in frames}
    if len(frame_by_key) != 800:
        raise ValueError("expected 800 unique cluster/filter image frames")
    prep_records = {str(row["id"]): row for row in prep_summary["records"]}
    original_records = {(str(row["cluster_id"]), row["filter"]): row
                       for row in low_summary["frame_records"]}
    if len(original_records) != 800:
        raise ValueError("40% profile frame records do not cover the 800-image sample")

    icl_keys = {(row["id"], row["filter"]) for row in read_tsv(ICL_SUPPORT)}
    if len(icl_keys) != 224:
        raise ValueError("expected 224 ICL 200–300 kpc supported cluster/filter pairs")
    candidate_rows = [row for row in read_tsv(SATELLITE_40/"candidate_diagnostics.tsv")
                      if (row["cluster_id"], row["filter"]) in icl_keys]
    profiles = read_tsv(SATELLITE_40/"satellite_profiles.tsv")
    outer_profile_keys = {(row["cluster_id"], row["filter"], row["member_index"])
                          for row in profiles if row["radial_bin_kpc"] == OUTER_BIN}
    profile_by_frame = defaultdict(set)
    for ident, band, member_index in outer_profile_keys:
        profile_by_frame[(ident, band)].add(member_index)

    cosmology = FlatLambdaCDM(H0=float(manifest.get("h0", 70.)),
                              Om0=float(manifest.get("omega_m", .3)), Tcmb0=0*u.K)
    candidate_output, frame_output, errors = [], [], []
    rows_by_frame = defaultdict(list)
    for row in candidate_rows:
        rows_by_frame[(row["cluster_id"], row["filter"])].append(row)

    for ordinal, key in enumerate(sorted(icl_keys), 1):
        frame = frame_by_key.get(key)
        source_record = original_records.get(key)
        if frame is None or source_record is None:
            raise ValueError(f"missing target image/provenance for ICL support {key}")
        ident, band = key
        prep_path = PREP / f"cluster_{ident}"
        paths = {"detections": prep_path/"detections.tsv", "members": prep_path/"members.tsv",
                 "masks": prep_path/"masks.fits.gz", "preparation": prep_path/"preparation.json"}
        expected_hashes = {
            "detections": source_record["detection_table_sha256"],
            "members": source_record["member_table_sha256"],
            "masks": source_record["mask_sha256"],
            "image": source_record["image_sha256"],
        }
        if any(not path.is_file() for path in paths.values()):
            raise FileNotFoundError(f"missing preparation artifact for {ident}")
        if any(sha256(paths[name]) != expected_hashes[name]
               for name in ("detections", "members", "masks")):
            raise ValueError(f"preparation artifact hash mismatch for {ident}")
        if sha256(Path(frame["image"])) != expected_hashes["image"]:
            raise ValueError(f"source image hash mismatch for {ident}")

        info = validate_fits(frame["image"], frame)
        det_rows = core.read_tsv(paths["detections"])
        member_rows = core.read_tsv(paths["members"])
        detections = core.parse_detections(det_rows)
        matched, _ = core.match_members(member_rows, detections, info["pixel_scale_arcsec"])
        det_arrays = core.detection_arrays(detections, info["pixel_scale_arcsec"])
        with fits.open(paths["masks"], memmap=False, checksum=True) as hdul:
            segmentation = np.asarray(hdul["SEGMENTATION"].data, dtype=np.int32)
            base = np.asarray(hdul["BASE_MASK"].data, dtype=bool)
        if (segmentation.shape != info["shape"] or base.shape != info["shape"]):
            raise ValueError(f"saved mask/image shape mismatch for {ident}")
        unmatched_indices = {item["member_index"] for item in matched}
        in_frame_members = [r for r in member_rows if r["IN_FRAME"] == "True"]
        fallback_xy = np.asarray([
            [float(r["X_IMAGE"])-1., float(r["Y_IMAGE"])-1.]
            for r in in_frame_members if int(r["INDEX"]) not in unmatched_indices
        ], dtype=float)
        fallback_tree = cKDTree(fallback_xy) if len(fallback_xy) else None
        det_by_number = {int(det["NUMBER"]): det for det in detections}
        matched_by_index = {int(item["member_index"]): item for item in matched}
        kpc_per_pixel = float(cosmology.kpc_proper_per_arcmin(float(frame["z"])).value
                              * info["pixel_scale_arcsec"] / 60.)
        outer_candidates = [r for r in rows_by_frame[key]
                            if len(r["target_coverage_by_bin"].split(",")) == len(core.R_EDGES_KPC)-1]
        n_coverage_pass = n_controls_pass = n_outer_profiles = 0
        for row in outer_candidates:
            member_index = int(row["member_index"])
            member_match = matched_by_index.get(member_index)
            det = det_by_number.get(int(row["detection_number"]))
            if member_match is None or det is None or member_match["detection_number"] != det["NUMBER"]:
                raise ValueError(f"candidate/member/detection mismatch for {key}/{member_index}")
            parts = coverage_parts(base, segmentation, det_arrays,
                                   fallback_tree, (det["x"], det["y"]), det["NUMBER"],
                                   info["pixel_scale_arcsec"], kpc_per_pixel)
            target_cov = float(row["target_coverage_by_bin"].split(",")[-1])
            if not math.isclose(parts["after_source_masks"], target_cov, rel_tol=0, abs_tol=0.00051):
                raise ValueError(f"reconstructed 50–75 kpc target coverage differs for {key}/{member_index}")
            valid_controls = ([int(v) for v in row["valid_controls_by_bin"].split(",") if v]
                              if row["valid_controls_by_bin"] else [])
            controls_n = valid_controls[-1] if valid_controls else 0
            profile_retained = (ident, band, str(member_index)) in outer_profile_keys
            exact_target_cov = parts["after_source_masks"]
            n_coverage_pass += exact_target_cov >= COVERAGE_FLOOR
            n_controls_pass += exact_target_cov >= COVERAGE_FLOOR and controls_n >= MIN_CONTROLS
            n_outer_profiles += profile_retained
            candidate_output.append({
                "cluster_id": ident, "filter": band, "member_index": member_index,
                "detection_number": int(row["detection_number"]),
                "cluster_radius_kpc": float(row["cluster_radius_kpc"]),
                "target_coverage_saved_40pct": target_cov,
                "footprint_coverage_50_75": parts["footprint"],
                "after_base_quality_mask_coverage_50_75": parts["after_base_mask"],
                "after_segmentation_mask_coverage_50_75": parts["after_segmentation_mask"],
                "after_detection_neighbor_masks_coverage_50_75": parts["after_detection_neighbor_masks"],
                "after_source_neighbor_masks_coverage_50_75": parts["after_source_masks"],
                "coverage_loss_to_image_edge": 1.-parts["footprint"],
                "coverage_loss_to_base_quality_mask": parts["footprint"]-parts["after_base_mask"],
                "coverage_loss_to_segmentation": parts["after_base_mask"]-parts["after_segmentation_mask"],
                "coverage_loss_to_detection_neighbor_masks": parts["after_segmentation_mask"]-parts["after_detection_neighbor_masks"],
                "coverage_loss_to_fallback_member_circles": parts["after_detection_neighbor_masks"]-parts["after_source_masks"],
                "coverage_loss_to_all_source_neighbor_masks": parts["after_base_mask"]-parts["after_source_masks"],
                "n_valid_azimuthal_controls": controls_n,
                "outer_annulus_profile_retained": profile_retained,
                "science_measurement": False,
            })
        frame_output.append({
            "cluster_id": ident, "filter": band, "n_candidate_members": len(outer_candidates),
            "n_target_coverage_pass_40pct": n_coverage_pass,
            "n_target_and_control_pass": n_controls_pass,
            "n_outer_annulus_profiles_retained": n_outer_profiles,
            "frame_support_50_75": bool(profile_by_frame.get(key)),
            "science_measurement": False,
        })
        if ordinal % 50 == 0:
            print(f"audited {ordinal}/224 ICL-supported frames", flush=True)

    if len(frame_output) != 224 or len(candidate_output) == 0:
        raise ValueError("coverage-loss audit did not cover the expected ICL support")
    output.mkdir(parents=True, exist_ok=False)
    candidates_path = output/"candidate_coverage_components.tsv"
    frames_path = output/"frame_support_components.tsv"
    write_tsv(candidates_path, candidate_output)
    write_tsv(frames_path, frame_output)
    components = {name: sum(max(0., float(row[name])) for row in candidate_output)
                  for name in ("coverage_loss_to_image_edge", "coverage_loss_to_base_quality_mask",
                               "coverage_loss_to_segmentation", "coverage_loss_to_detection_neighbor_masks",
                               "coverage_loss_to_fallback_member_circles",
                               "coverage_loss_to_all_source_neighbor_masks")}
    count_causes = Counter()
    for row in candidate_output:
        # Candidate diagnostics serialize coverage to three decimals. Use the
        # unrounded pixel reconstruction for threshold classification so a
        # true value such as 0.39998 is not promoted by formatting to 0.400.
        if row["footprint_coverage_50_75"] < COVERAGE_FLOOR:
            count_causes["footprint_below_40pct"] += 1
        elif row["after_base_quality_mask_coverage_50_75"] < COVERAGE_FLOOR:
            count_causes["base_quality_mask_below_40pct"] += 1
        elif row["after_segmentation_mask_coverage_50_75"] < COVERAGE_FLOOR:
            count_causes["segmentation_mask_below_40pct"] += 1
        elif row["after_detection_neighbor_masks_coverage_50_75"] < COVERAGE_FLOOR:
            count_causes["detection_neighbor_masks_below_40pct"] += 1
        elif row["after_source_neighbor_masks_coverage_50_75"] < COVERAGE_FLOOR:
            count_causes["source_neighbor_masks_reduce_below_40pct"] += 1
        else:
            count_causes["target_pass"] += 1

    inputs = (MANIFEST, PREP_SUMMARY, PREP_AUDIT, ICL_SUPPORT,
              SATELLITE_40/"summary.json", SATELLITE_40/"candidate_diagnostics.tsv",
              SATELLITE_40/"satellite_profiles.tsv")
    summary = {
        "schema": "astrafek-icl-satellite-outer-coverage-loss-v3",
        "status": "complete_offline_diagnostic", "science_measurement": False,
        "hsc_network_used": False, "science_image_pixels_read": False,
        "prepared_mask_pixels_read": True, "n_icl_supported_frames": len(frame_output),
        "n_member_candidates_audited": len(candidate_output),
        "coverage_floor": COVERAGE_FLOOR, "minimum_azimuthal_controls": MIN_CONTROLS,
        "candidate_outcome_counts": dict(count_causes),
        "summed_coverage_deficit_diagnostics_not_additive": components,
        "interpretation": "candidate-level projected annulus area decomposition; stage-specific floor outcomes use unrounded pixel coverage, while saved candidate diagnostics are rounded to three decimals; losses are sequential and do not define independent uncertainties",
        "limitations": [
            "The area decomposition follows saved BASE_MASK, SEP segmentation, Kron neighbor masks, and fallback member circles; it is not a flux-bias or completeness correction.",
            "The 50–75 kpc annulus is audited only for candidates in ICL 200–300 kpc supported cluster/filter frames.",
            "No unmasked/relaxed-floor measurements are promoted; all results are diagnostic and science_measurement remains false.",
        ],
        "input_sha256": {str(path): sha256(path) for path in inputs},
        "script_sha256": sha256(Path(__file__).resolve()),
        "output_sha256": {str(candidates_path): sha256(candidates_path),
                          str(frames_path): sha256(frames_path)},
    }
    summary_path = output/"summary.json"
    summary_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True)+"\n",
                            encoding="utf-8")
    print(json.dumps({"status": summary["status"], "candidate_count": len(candidate_output),
                      "outcome_counts": dict(count_causes), "output": str(output)},
                     ensure_ascii=False, indent=2))
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    run(args.output.resolve())


if __name__ == "__main__":
    main()
