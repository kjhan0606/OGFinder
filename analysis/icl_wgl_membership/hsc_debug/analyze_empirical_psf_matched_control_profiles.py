#!/usr/bin/env python3
"""Compare Gaia-star wings with same-offset matched-control sky profiles.

This diagnostic reads only existing global-sky target/control cutouts and saved
star stamps. It estimates local residual profiles; target-minus-control is not
interpreted as a PSF-only signal because cluster ICL curvature can remain.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import platform
from collections import defaultdict
from pathlib import Path

import numpy as np
from astropy.io import fits

from fetch_ancillary_catalogues import sha
from stack_empirical_psf_star_profiles import RADIAL_EDGES, SKY_ANNULI, measure_stamp
import prepare_detection_psf as prep

BASE = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
SOURCE_PLANES = ("MP_DETECTED", "MP_DETECTED_NEGATIVE", "MP_NOT_DEBLENDED")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def control_profile(image, variance, mask, x: float, y: float, core_flux: float,
                    mask_status: str = "ok") -> list[dict]:
    """Measure target-core-normalized local residuals at a masked control point."""
    image = np.asarray(image, dtype=float)
    variance = np.asarray(variance, dtype=float)
    mask = np.asarray(mask, dtype=bool)
    if image.shape != (401, 401) or variance.shape != image.shape or mask.shape != image.shape:
        raise ValueError("control stamp planes must all be 401x401")
    if not math.isfinite(core_flux) or core_flux <= 0:
        raise ValueError("target star core flux must be finite and positive")
    yy, xx = np.indices(image.shape, dtype=float)
    radius = np.hypot(xx - x, yy - y)
    angle = np.arctan2(yy - y, xx - x)
    sectors = np.clip(np.floor((angle + np.pi) * (4 / np.pi)).astype(int), 0, 7)
    good = ~mask & np.isfinite(image) & np.isfinite(variance) & (variance > 0)
    core_region = radius < 12
    core_fraction = float(np.count_nonzero(core_region & good) / np.count_nonzero(core_region))
    output = []
    edges = np.asarray(RADIAL_EDGES, dtype=float)
    for sky_lo, sky_hi in SKY_ANNULI:
        sky_region = (radius >= sky_lo) & (radius < sky_hi)
        sector_means, sector_counts, sector_fractions = [], [], []
        for sector_id in range(8):
            region = sky_region & (sectors == sector_id)
            valid = region & good
            count = int(valid.sum())
            sector_counts.append(count)
            sector_fractions.append(count / max(1, int(region.sum())))
            if count:
                sector_means.append(float(np.mean(image[valid])))
        sky_good = sky_region & good
        sky_fraction = float(sky_good.sum() / max(1, sky_region.sum()))
        sky_ok = (sky_good.sum() >= 300 and sky_fraction >= 0.30
                  and min(sector_counts) >= 20 and min(sector_fractions) >= 0.25)
        background = float(np.mean(sector_means)) if sky_ok else None
        for lo, hi in zip(edges[:-1], edges[1:]):
            ring = (radius >= lo) & (radius < hi)
            valid = ring & good
            area = int(ring.sum())
            count = int(valid.sum())
            covered = sum(
                np.count_nonzero(valid & (sectors == sector_id))
                >= max(1, 0.25 * np.count_nonzero(ring & (sectors == sector_id)))
                for sector_id in range(8)
            )
            status = "ok"
            if mask_status != "ok":
                status = mask_status
            elif lo < 12 and core_fraction < 0.99:
                status = "blocked_control_center_source_or_mask"
            elif not sky_ok:
                status = "blocked_background_coverage"
            elif count / max(1, area) < 0.50 or covered != 8:
                status = "blocked_radial_coverage"
            value = None
            raw = float(np.mean(image[valid])) if count else None
            if status == "ok":
                value = (raw - background) / core_flux
            output.append({
                "sky_r_min_pix": sky_lo, "sky_r_max_pix": sky_hi,
                "r_min_pix": float(lo), "r_max_pix": float(hi),
                "n_valid": count, "area_fraction": count / max(1, area),
                "n_covered_sectors": covered, "sky_area_fraction": sky_fraction,
                "background": background, "raw_mean": raw,
                "core_open_fraction": core_fraction,
                "profile_per_target_core_per_pixel": value, "status": status,
            })
    return output


def detect_control_masks(image, variance, base_mask, x: float, y: float, pixel_scale: float):
    """Build target-like SEP masks; return explicit per-stamp support failures."""
    try:
        objects, segmentation, kron, kron_flags, _, _, _ = prep.detect(
            image, variance, base_mask, np.asarray([x, y]), np.empty((0, 2)), pixel_scale,
        )
    except ValueError as error:
        if "Too few valid pixels" not in str(error):
            raise
        blocked = np.ones((401, 401), dtype=bool)
        return [blocked.copy(), blocked.copy(), blocked.copy()], "blocked_control_too_few_quality_pixels"
    masks, _ = prep.source_masks(
        objects, segmentation, kron, kron_flags, base_mask, -1,
        np.empty((0, 2)), pixel_scale,
    )
    return masks, "ok"


def _quality_mask_bitmask(control_manifest_row: dict) -> int:
    bits = 0
    for value in control_manifest_row["bad_mask"]["exclude_bits"]:
        bits |= int(value)
    for value in control_manifest_row["bright_star_mask"]["exclude_bits"]:
        bits |= int(value)
    return bits


def _write_tsv(path: Path, rows: list[dict]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    with path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t", extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def frame_summary(rows: list[dict]) -> list[dict]:
    grouped = defaultdict(lambda: defaultdict(list))
    for row in rows:
        if row["pair_status"] != "complete_profile_pair":
            continue
        key = (row["filter"], float(row["kron_factor"]), int(row["sky_r_min_pix"]),
               int(row["sky_r_max_pix"]), float(row["r_min_pix"]), float(row["r_max_pix"]))
        grouped[key][str(row["target_cluster_id"])].append(row)
    summary = []
    for (band, kron, sky_lo, sky_hi, rlo, rhi), frames in sorted(grouped.items()):
        triples = [tuple(
            float(np.mean([float(row[name]) for row in values]))
            for name in ("target_profile", "control_profile", "target_minus_control")
        ) for values in frames.values()]
        values = np.asarray(list(triples), dtype=float)
        summary.append({
            "filter": band, "kron_factor": kron, "sky_r_min_pix": sky_lo, "sky_r_max_pix": sky_hi,
            "r_min_pix": rlo, "r_max_pix": rhi, "n_frames": len(values),
            "target_frame_equal_mean": float(np.mean(values[:, 0])),
            "control_frame_equal_mean": float(np.mean(values[:, 1])),
            "paired_difference_frame_equal_mean": float(np.mean(values[:, 2])),
            "paired_difference_frame_q16": float(np.quantile(values[:, 2], 0.16)),
            "paired_difference_frame_median": float(np.median(values[:, 2])),
            "paired_difference_frame_q84": float(np.quantile(values[:, 2], 0.84)),
            "interval_type": "descriptive_distribution_across_frames_not_confidence_interval",
            "science_measurement": False, "psf_approved": False,
        })
    return summary


def run(sample_root: Path, support_dir: Path, ledger_path: Path, output_dir: Path) -> dict:
    if output_dir.exists():
        raise FileExistsError(f"refusing to overwrite {output_dir}")
    support_summary_path = support_dir / "summary.json"
    support_rows_path = support_dir / "matched_control_support.tsv"
    target_manifest_path = sample_root / "download_01/download_manifest.json"
    control_manifest_path = sample_root / "control_download_01/download_manifest.json"
    target_status_path = sample_root / "download_01/download_status.json"
    control_status_path = sample_root / "control_download_01/download_status.json"
    inputs = (support_summary_path, support_rows_path, target_manifest_path, control_manifest_path,
              target_status_path, control_status_path, ledger_path)
    for path in inputs:
        if not path.is_file():
            raise FileNotFoundError(path)
    support_summary = json.loads(support_summary_path.read_text(encoding="utf-8"))
    if support_summary.get("status") != "complete_header_wcs_preflight":
        raise ValueError("matched-control support preflight is not complete")
    support = [row for row in read_tsv(support_rows_path) if row["pair_status"] == "same_filter_offset_available"]
    target_doc = json.loads(target_manifest_path.read_text(encoding="utf-8"))
    control_doc = json.loads(control_manifest_path.read_text(encoding="utf-8"))
    target_status = json.loads(target_status_path.read_text(encoding="utf-8"))
    control_status = json.loads(control_status_path.read_text(encoding="utf-8"))
    target_by_id = {str(row["id"]): row for row in target_doc["clusters"]}
    control_by_id = {str(row["id"]): row for row in control_doc["clusters"]}
    control_by_target = {str(row["target_cluster_id"]): row for row in control_doc["clusters"]}
    control_status_by_id = {str(row["id"]): row for row in control_status["completed"]}
    ledger_by_id = {row["candidate_id"]: row for row in read_tsv(ledger_path) if row.get("stamp")}
    if len(ledger_by_id) != len([row for row in read_tsv(ledger_path) if row.get("stamp")]):
        raise ValueError("duplicate candidate ID in Gaia stamp ledger")

    grouped = defaultdict(list)
    for row in support:
        if row["candidate_id"] not in ledger_by_id:
            raise ValueError(f"support row has no Gaia stamp: {row['candidate_id']}")
        grouped[str(row["target_cluster_id"])].append(row)

    pair_rows = []
    for target_id, candidates in grouped.items():
        target = target_by_id[target_id]
        control = control_by_target[target_id]
        control_status_row = control_status_by_id[str(control["id"])]
        control_path = sample_root / "control_download_01/images" / control_status_row["file"]
        science_hdu = int(control_status_row["review"]["science_hdu"])
        variance_hdu = int(control_status_row["review"]["variance_hdu"])
        mask_hdu = int(control_status_row["review"]["mask_hdu"])
        expected_shape = tuple(int(v) for v in control_status_row["review"]["shape"])
        quality_mask_bits = _quality_mask_bitmask(control)
        pixel_scale = float(control_status_row["review"]["pixel_scale_arcsec"])
        with fits.open(control_path, memmap=True, lazy_load_hdus=True) as hdus:
            for item in candidates:
                stamp_row = ledger_by_id[item["candidate_id"]]
                stamp_path = Path(stamp_row["stamp"])
                if sha(stamp_path) != stamp_row["stamp_sha256"]:
                    raise ValueError(f"star stamp checksum mismatch: {stamp_path}")
                target_profiles = measure_stamp(stamp_path, target_id, item["candidate_id"])
                x, y = float(item["control_x"]), float(item["control_y"])
                x0, y0 = math.floor(x) - 200, math.floor(y) - 200
                x1, y1 = x0 + 401, y0 + 401
                if x0 < 0 or y0 < 0 or x1 > expected_shape[1] or y1 > expected_shape[0]:
                    raise ValueError(f"control stamp edge changed after support audit: {item['candidate_id']}")
                sl = np.s_[y0:y1, x0:x1]
                image = np.asarray(hdus[science_hdu].section[sl], dtype=float)
                variance = np.asarray(hdus[variance_hdu].section[sl], dtype=float)
                raw_mask = np.asarray(hdus[mask_hdu].section[sl])
                center_x, center_y = x - x0, y - y0
                base_mask = (raw_mask & quality_mask_bits) != 0
                base_mask |= ~np.isfinite(image) | ~np.isfinite(variance) | (variance <= 0)
                control_masks, mask_status = detect_control_masks(
                    image, variance, base_mask, center_x, center_y, pixel_scale,
                )
                target_by_method = {}
                for row in target_profiles:
                    key = (float(row["kron_factor"]), int(row["sky_r_min_pix"]), int(row["sky_r_max_pix"]))
                    target_by_method.setdefault(key, {})[(float(row["r_min_pix"]), float(row["r_max_pix"]))] = row
                control_by_method = {}
                for (kron, sky_lo, sky_hi), target_rings in target_by_method.items():
                    if (sky_lo, sky_hi) not in SKY_ANNULI:
                        raise ValueError(f"unexpected sky annulus: {(sky_lo, sky_hi)}")
                    representative = next(iter(target_rings.values()))
                    core_flux = representative.get("core_flux")
                    mask_index = {2.5: 0, 3.5: 1}.get(kron)
                    if mask_index is None:
                        raise ValueError(f"unexpected target Kron mask: {kron}")
                    control_mask = control_masks[mask_index]
                    if core_flux is not None and math.isfinite(float(core_flux)) and float(core_flux) > 0:
                        profile_rows = control_profile(image, variance, control_mask, center_x, center_y,
                                                       float(core_flux), mask_status=mask_status)
                        control_by_method[(kron, sky_lo, sky_hi)] = {
                            (float(r["r_min_pix"]), float(r["r_max_pix"])): r for r in profile_rows
                            if (int(r["sky_r_min_pix"]), int(r["sky_r_max_pix"])) == (sky_lo, sky_hi)
                        }
                    else:
                        control_by_method[(kron, sky_lo, sky_hi)] = {
                            (float(r["r_min_pix"]), float(r["r_max_pix"])): {
                                "status": "blocked_target_core_flux", "profile_per_target_core_per_pixel": None,
                                "area_fraction": 0.0, "sky_area_fraction": 0.0, "core_open_fraction": 0.0,
                            }
                            for r in target_rings.values()
                        }
                for (kron, sky_lo, sky_hi), target_rings in target_by_method.items():
                    for ring_key, target_ring in target_rings.items():
                        control_ring = control_by_method[(kron, sky_lo, sky_hi)][ring_key]
                        target_value = target_ring.get("profile_per_core_per_pixel")
                        control_value = control_ring.get("profile_per_target_core_per_pixel")
                        pair_status = "complete_profile_pair"
                        if target_ring.get("status") != "ok":
                            pair_status = "blocked_target_" + str(target_ring.get("status"))
                        elif control_ring.get("status") != "ok":
                            pair_status = str(control_ring.get("status"))
                        if pair_status == "complete_profile_pair" and (target_value is None or control_value is None):
                            raise ValueError("valid paired profile has a null value")
                        pair_rows.append({
                            "candidate_id": item["candidate_id"], "target_cluster_id": target_id,
                            "control_id": control["id"], "filter": item["target_filter"],
                            "kron_factor": kron, "sky_r_min_pix": sky_lo, "sky_r_max_pix": sky_hi,
                            "r_min_pix": ring_key[0], "r_max_pix": ring_key[1],
                            "target_profile": target_value, "control_profile": control_value,
                            "target_minus_control": (float(target_value) - float(control_value)
                                                     if target_value is not None and control_value is not None else None),
                            "control_area_fraction": control_ring["area_fraction"],
                            "control_sky_area_fraction": control_ring["sky_area_fraction"],
                            "control_core_open_fraction": control_ring["core_open_fraction"],
                            "pair_status": pair_status,
                            "science_measurement": False, "psf_approved": False,
                        })
        print(f"target {target_id}: {len(candidates)} paired star stamps", flush=True)

    stack = frame_summary(pair_rows)
    output_dir.mkdir(parents=True, exist_ok=False)
    pair_path = output_dir / "paired_profiles.tsv"
    stack_path = output_dir / "frame_weighted_profiles.tsv"
    _write_tsv(pair_path, pair_rows)
    _write_tsv(stack_path, stack)
    status_counts = defaultdict(int)
    for row in pair_rows:
        status_counts[row["pair_status"]] += 1
    summary = {
        "schema": "astrafek-empirical-psf-matched-control-profiles-v1",
        "status": "complete_exploratory_profile_diagnostic",
        "n_same_filter_star_stamps": len(support), "n_target_frames": len(grouped),
        "n_profile_pair_rows": len(pair_rows), "n_frame_weighted_profile_bins": len(stack),
        "pair_status_counts": dict(status_counts),
        "control_source_mask": "Astrafek SEP detections with matching Kron 2.5/3.5 masks; HSC bad/bright pixels remain excluded",
        "control_center_gate": "center openness is required only for bins inside 12 pixels; outer bins are gated by their own masked area and 8-sector coverage",
        "bootstrap": False,
        "profile_intervals": "q16/q84 are descriptive distributions across frames, not confidence intervals",
        "limitations": [
            "The control position has no matched star; control profiles measure local background residuals/placebos.",
            "Target-minus-control can retain target ICL curvature, source selection differences, and local sky mismatch; it is not a PSF-only estimate.",
            "Control SEP masks are redetected on each matched-control stamp; the segmentation can differ from the original target detection catalog.",
            "No coadd covariance, photon-noise injection, Gaia epoch correction, or PSF approval is included.",
        ],
        "input_sha256": {str(path.resolve()): sha(path) for path in inputs},
        "script_sha256": sha(Path(__file__).resolve()),
        "paired_profiles_sha256": sha(pair_path), "frame_weighted_profiles_sha256": sha(stack_path),
        "science_measurement": False, "psf_approved": False,
    }
    summary_path = output_dir / "summary.json"
    summary_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-root", type=Path, default=BASE)
    parser.add_argument("--support", type=Path, default=BASE / "empirical_psf_matched_control_support_2026-10-20_01")
    parser.add_argument("--ledger", type=Path, default=BASE / "empirical_psf_gaia_full_sample_2026-10-20_02/star_ledger.tsv")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if not platform.node().lower().startswith("lageunha"):
        parser.error("Run this existing-image CPU profile diagnostic on lageunha")
    result = run(args.sample_root.resolve(), args.support.resolve(), args.ledger.resolve(), args.output.resolve())
    print(json.dumps({key: result[key] for key in (
        "status", "n_same_filter_star_stamps", "n_target_frames", "n_profile_pair_rows",
        "pair_status_counts", "science_measurement", "psf_approved",
    )}, ensure_ascii=False, indent=2))
    print(f"summary_sha256={sha(args.output / 'summary.json')}")


if __name__ == "__main__":
    main()
