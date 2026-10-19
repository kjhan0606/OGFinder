#!/usr/bin/env python3
"""Build exploratory HSC-I 2-D ICL target-minus-matched-control stacks.

This offline Astrafek analysis tool uses the existing 800 target coadds, their
reviewed K=2.5 masks, and 798 downloaded matched-control coadds. It performs no
network access and never modifies source data. Maps are resampled to a common
proper-kpc grid, paired on the intersection of valid target/control pixels,
and averaged with equal pair weight. The product is exploratory, not a final
science measurement.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path
import sys
import time

import numpy as np
from astropy.io import fits

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))
if str(HERE.parents[2]) not in sys.path:
    sys.path.insert(0, str(HERE.parents[2]))

from download_pilot import validate_fits
from icl.wgl import warp_metric
from ogfkit.wlcluster import Cosmology
from validate_control_images import control_masks, image_path as control_image_path, load_controls
from validate_control_images import nJy_factor


DEFAULT_SAMPLE = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
DEFAULT_WORK = DEFAULT_SAMPLE / "icl_2d_matched_stack_2026-10-19_01"
FILTERS = ("i", "i2")
PROFILE_EDGES_KPC = np.arange(0.0, 525.0, 25.0)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_json(path: Path) -> dict:
    with Path(path).open(encoding="utf-8") as stream:
        return json.load(stream)


def read_tsv(path: Path) -> list[dict[str, str]]:
    with Path(path).open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def json_safe(value):
    if isinstance(value, dict):
        return {str(k): json_safe(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [json_safe(v) for v in value]
    if isinstance(value, np.ndarray):
        return json_safe(value.tolist())
    if isinstance(value, (np.bool_, bool)):
        return bool(value)
    if isinstance(value, (np.integer, int)):
        return int(value)
    if isinstance(value, (np.floating, float)):
        return float(value) if np.isfinite(value) else None
    return value


def write_json(path: Path, value: dict) -> None:
    path.write_text(json.dumps(json_safe(value), ensure_ascii=False, indent=2,
                               sort_keys=True, allow_nan=False) + "\n", encoding="utf-8")


class OnlineMap:
    """Pixelwise Welford mean/SEM for a changing set of valid input maps."""

    def __init__(self, shape: tuple[int, int]):
        self.count = np.zeros(shape, dtype=np.uint16)
        self.mean = np.zeros(shape, dtype=np.float64)
        self.m2 = np.zeros(shape, dtype=np.float64)

    def add(self, values: np.ndarray, valid: np.ndarray) -> None:
        values = np.asarray(values, dtype=np.float64)
        valid = np.asarray(valid, dtype=bool) & np.isfinite(values)
        old_count = self.count[valid].astype(np.float64)
        new_count = old_count + 1.0
        delta = values[valid] - self.mean[valid]
        self.mean[valid] += delta / new_count
        delta2 = values[valid] - self.mean[valid]
        self.m2[valid] += delta * delta2
        if np.any(self.count[valid] == np.iinfo(np.uint16).max):
            raise OverflowError("pair count exceeds uint16 capacity")
        self.count[valid] += 1

    def result(self) -> tuple[np.ndarray, np.ndarray]:
        mean = np.full(self.mean.shape, np.nan, dtype=np.float32)
        sem = np.full(self.mean.shape, np.nan, dtype=np.float32)
        occupied = self.count > 0
        mean[occupied] = self.mean[occupied].astype(np.float32)
        many = self.count > 1
        sem[many] = np.sqrt(self.m2[many] / (self.count[many] - 1) / self.count[many]).astype(np.float32)
        return mean, sem


def pair_common_support(target: np.ndarray, control: np.ndarray,
                        target_weight: np.ndarray, control_weight: np.ndarray):
    """Return finite pair maps on the strict intersection of both warped supports."""
    target = np.asarray(target, dtype=np.float64)
    control = np.asarray(control, dtype=np.float64)
    valid = (np.asarray(target_weight) > 0) & (np.asarray(control_weight) > 0)
    valid &= np.isfinite(target) & np.isfinite(control)
    delta = np.full(target.shape, np.nan, dtype=np.float64)
    delta[valid] = target[valid] - control[valid]
    return target, control, delta, valid


def pair_radial_profile(target: np.ndarray, control: np.ndarray, delta: np.ndarray,
                        valid: np.ndarray, axis_kpc: np.ndarray,
                        edges_kpc: np.ndarray = PROFILE_EDGES_KPC) -> list[dict]:
    yy, xx = np.meshgrid(axis_kpc, axis_kpc, indexing="ij")
    radius = np.hypot(xx, yy)
    rows = []
    for lo, hi in zip(edges_kpc[:-1], edges_kpc[1:]):
        annulus = (radius >= lo) & (radius < hi)
        geometric = int(annulus.sum())
        selected = annulus & valid
        n_valid = int(selected.sum())
        rows.append({
            "r_min_kpc": float(lo), "r_max_kpc": float(hi),
            "n_valid_pixels": n_valid,
            "geometric_pixels": geometric,
            "coverage_fraction": n_valid / geometric if geometric else 0.0,
            "target_mean_njy_arcsec2": float(np.mean(target[selected])) if n_valid else None,
            "control_mean_njy_arcsec2": float(np.mean(control[selected])) if n_valid else None,
            "delta_mean_njy_arcsec2": float(np.mean(delta[selected])) if n_valid else None,
        })
    return rows


def summarize_profiles(pair_rows: list[dict], filters: tuple[str, ...] = FILTERS) -> list[dict]:
    """Equal-weight pair-level radial profiles and their between-pair SEM."""
    output = []
    for band in filters:
        selected = [row for row in pair_rows if row["filter"] == band]
        for bin_index, (lo, hi) in enumerate(zip(PROFILE_EDGES_KPC[:-1], PROFILE_EDGES_KPC[1:])):
            base = {"filter": band, "r_min_kpc": float(lo), "r_max_kpc": float(hi)}
            for key, stem in (("target_mean_njy_arcsec2", "target"),
                              ("control_mean_njy_arcsec2", "control"),
                              ("delta_mean_njy_arcsec2", "delta")):
                vals = np.asarray([row["profile"][bin_index][key] for row in selected
                                   if row["profile"][bin_index][key] is not None], dtype=float)
                base[f"n_pairs_{stem}"] = int(vals.size)
                base[f"mean_{stem}_njy_arcsec2"] = float(vals.mean()) if vals.size else None
                base[f"pair_sd_{stem}_njy_arcsec2"] = float(vals.std(ddof=1)) if vals.size > 1 else None
                base[f"pair_sem_{stem}_njy_arcsec2"] = float(vals.std(ddof=1) / math.sqrt(vals.size)) if vals.size > 1 else None
            cov = np.asarray([row["profile"][bin_index]["coverage_fraction"] for row in selected], dtype=float)
            base["median_pair_common_coverage"] = float(np.median(cov)) if cov.size else None
            output.append(base)
    return output


def pixel_flux_factor(info: dict) -> float:
    """Native input flux/pixel to nJy/pixel; warp_metric divides by pixel area."""
    return float(nJy_factor(info) * info["pixel_scale_arcsec"] ** 2)


def target_and_control_map(image: np.ndarray, mask: np.ndarray, info: dict,
                           redshift: float, cosmology: Cosmology,
                           kpc_max: float, n_side: int) -> dict:
    frame = warp_metric(image * pixel_flux_factor(info), info["center_pixel_xy"],
                        info["pixel_scale_arcsec"], redshift, cosmo=cosmology,
                        mask=mask, kpc_max=kpc_max, n_side=n_side)
    return frame


def process_pair(target: dict, target_mask_record: dict, target_mask_summary: dict,
                 control: dict, control_download_record: dict, control_base: Path,
                 control_status_sha: str, cosmology: Cosmology,
                 kpc_max: float, n_side: int) -> tuple[dict, dict]:
    ident = str(target["id"])
    control_id = str(control["id"])
    target_path = Path(target_mask_record["image"])
    target_mask_path = Path(target_mask_record["mask_file"])
    control_path = control_image_path(control_base, control)
    if not target_path.is_file() or not target_mask_path.is_file() or not control_path.is_file():
        raise FileNotFoundError(f"{ident}: target, target mask, or control image missing")
    for path, expected, label in (
        (target_path, target_mask_record["image_sha256"], "target image"),
        (target_mask_path, target_mask_record["mask_sha256"], "target mask"),
        (control_path, control_download_record["sha256"], "control image"),
    ):
        actual = sha256(path)
        if actual != expected:
            raise ValueError(f"{ident}: {label} SHA-256 mismatch")

    target_info = validate_fits(target_path, target)
    control_info = validate_fits(control_path, control)
    band = target_info["coadd_header_filter"]
    if band not in FILTERS or control_info["coadd_header_filter"] != band:
        raise ValueError(f"{ident}: target/control HSC-I family filter mismatch")
    if target_mask_record["primary_filter"] != band or control.get("coadd_header_filter") != band:
        raise ValueError(f"{ident}: stored target/control filter evidence mismatch")
    if (str(control.get("target_cluster_id")) != ident
            or control.get("frame_role") != "matched_control"
            or control.get("bcg_policy") != "none"
            or control.get("science_measurement") is not False
            or control.get("members_radec_deg") != []
            or control.get("randoms") != []
            or float(control["z"]) != float(target["z"])):
        raise ValueError(f"{ident}: control role, pairing, or assigned redshift mismatch")

    with fits.open(target_path, memmap=False, checksum=True) as hdul:
        target_image = np.asarray(hdul[target_info["science_hdu"]].data, dtype=np.float32)
    with fits.open(target_mask_path, memmap=False, checksum=True) as hdul:
        target_mask = np.asarray(hdul["MASK_K2P5"].data, dtype=bool)
    if target_mask_summary.get("mask_output_sha256", {}).get(str(target_mask_path)) != target_mask_record["mask_sha256"]:
        raise ValueError(f"{ident}: target mask is absent from verified mask summary")
    if target_image.shape != target_mask.shape:
        raise ValueError(f"{ident}: target image/mask shapes differ")

    with fits.open(control_path, memmap=False, checksum=True) as hdul:
        control_image = np.asarray(hdul[control_info["science_hdu"]].data, dtype=np.float32)
        control_variance = np.asarray(hdul[control_info["variance_hdu"]].data, dtype=np.float32)
        control_bits = np.asarray(hdul[control_info["mask_hdu"]].data, dtype=np.int32)
    if control_image.shape != control_variance.shape or control_bits.shape != control_image.shape:
        raise ValueError(f"{ident}: control image/variance/mask shapes differ")
    _, control_mask_variants, _, _, _ = control_masks(
        control_image, control_variance, control_bits,
        control_info["center_pixel_xy"], control_info["pixel_scale_arcsec"],
    )
    control_mask = np.asarray(control_mask_variants[0], dtype=bool)

    target_warp = target_and_control_map(target_image, target_mask, target_info,
                                         float(target["z"]), cosmology, kpc_max, n_side)
    control_warp = target_and_control_map(control_image, control_mask, control_info,
                                          float(target["z"]), cosmology, kpc_max, n_side)
    tmap, cmap, dmap, valid = pair_common_support(
        target_warp["sb"], control_warp["sb"], target_warp["weight"], control_warp["weight"])
    if not np.any(valid):
        raise ValueError(f"{ident}: target/control warped maps have no common pixels")
    axis = np.asarray(target_warp["kpc"], dtype=float)
    pair = {
        "id": ident, "control_id": control_id, "filter": band,
        "z": float(target["z"]), "ra_deg": float(target["ra_deg"]),
        "dec_deg": float(target["dec_deg"]), "n_common_pixels": int(valid.sum()),
        "common_grid_fraction": float(valid.mean()),
        "target_image": str(target_path), "target_image_sha256": target_mask_record["image_sha256"],
        "target_mask": str(target_mask_path), "target_mask_sha256": target_mask_record["mask_sha256"],
        "control_image": str(control_path), "control_image_sha256": control_download_record["sha256"],
        "control_download_status_sha256": control_status_sha,
        "profile": pair_radial_profile(tmap, cmap, dmap, valid, axis),
    }
    maps = {"target": tmap, "control": cmap, "delta": dmap, "valid": valid, "axis_kpc": axis}
    return pair, maps


def write_pair_table(path: Path, rows: list[dict]) -> None:
    fields = ["target_id", "control_id", "filter", "z", "ra_deg", "dec_deg",
              "n_common_pixels", "common_grid_fraction", "profile_json"]
    with path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", extrasaction="raise")
        writer.writeheader()
        for row in rows:
            writer.writerow({"target_id": row["id"], "control_id": row["control_id"],
                             "filter": row["filter"], "z": row["z"], "ra_deg": row["ra_deg"],
                             "dec_deg": row["dec_deg"], "n_common_pixels": row["n_common_pixels"],
                             "common_grid_fraction": row["common_grid_fraction"],
                             "profile_json": json.dumps(row["profile"], separators=(",", ":"),
                                                        allow_nan=False)})


def save_band_fits(path: Path, band: str, kpc: np.ndarray, maps: dict,
                   accumulators: dict[str, OnlineMap], n_pairs: int,
                   cosmology: dict, kpc_max: float) -> dict:
    target, target_sem = accumulators["target"].result()
    control, control_sem = accumulators["control"].result()
    delta, delta_sem = accumulators["delta"].result()
    n_pair = accumulators["delta"].count
    y, x = np.meshgrid(kpc, kpc, indexing="ij")
    radius = np.hypot(x, y).astype(np.float32)
    header = fits.Header()
    header["FILTER"] = band
    header["NPAIR"] = int(n_pairs)
    header["KPCMAX"] = float(kpc_max)
    header["GRIDSTEP"] = float(kpc[1] - kpc[0])
    header["BUNIT"] = "nJy/arcsec2"
    header["SCIVAL"] = False
    header["H0"] = float(cosmology["H0"])
    header["OMEGAM"] = float(cosmology["Om0"])
    header["MASK"] = "Kron-2.5 target/control; common support"
    hdus = [fits.PrimaryHDU(header=header)]
    for name, data, unit in (
        ("TARGET", target, "nJy/arcsec2"), ("CONTROL", control, "nJy/arcsec2"),
        ("DELTA", delta, "nJy/arcsec2"), ("NPAIR", n_pair.astype(np.uint16), "pair"),
        ("TARGET_SEM", target_sem, "nJy/arcsec2"), ("CONTROL_SEM", control_sem, "nJy/arcsec2"),
        ("DELTA_SEM", delta_sem, "nJy/arcsec2"), ("GRID_KPC", radius, "kpc"),
    ):
        h = fits.Header(); h["BUNIT"] = unit; h["FILTER"] = band; h["SCIVAL"] = False
        hdus.append(fits.ImageHDU(data=data, header=h, name=name))
    fits.HDUList(hdus).writeto(path, overwrite=False, checksum=True)
    finite = np.isfinite(delta)
    return {"file": str(path), "sha256": sha256(path), "n_pairs": n_pairs,
            "shape": list(delta.shape), "covered_pixels": int(finite.sum()),
            "n_pair_min": int(n_pair[n_pair > 0].min()) if np.any(n_pair) else 0,
            "n_pair_max": int(n_pair.max()), "delta_finite_range_njy_arcsec2":
                [float(np.min(delta[finite])), float(np.max(delta[finite]))] if finite.any() else None}


def run(sample: Path, work: Path, kpc_max: float, n_side: int) -> dict:
    if work.exists():
        raise FileExistsError(f"refusing to overwrite existing work directory: {work}")
    work.mkdir(parents=True, exist_ok=False)
    started = time.monotonic()
    sample_manifest_path = sample / "sample_manifest.json"
    target_mask_dir = sample / "prestack_cluster_masks_2026-10-16_01"
    target_mask_summary_path = target_mask_dir / "summary.json"
    target_mask_records_path = target_mask_dir / "records.tsv"
    sample_doc = read_json(sample_manifest_path)
    mask_summary = read_json(target_mask_summary_path)
    targets = sample_doc.get("clusters", [])
    mask_rows = read_tsv(target_mask_records_path)
    target_by_id = {str(frame["id"]): frame for frame in targets}
    mask_by_id = {str(row["id"]): row for row in mask_rows}
    if len(targets) != 800 or len(target_by_id) != 800 or len(mask_by_id) != 800:
        raise ValueError("expected 800 unique targets and 800 target-mask records")
    if (mask_summary.get("status") != "complete" or mask_summary.get("n_success") != 800
            or mask_summary.get("science_measurement") is not False
            or mask_summary.get("mask_ready_for_science") is not False
            or sha256(target_mask_records_path) != mask_summary.get("output_sha256", {}).get(str(target_mask_records_path))):
        raise ValueError("target mask manifest is incomplete or provenance-invalid")

    controls_doc, control_downloads, control_base, control_source_hashes = load_controls(sample)
    control_status_path = sample / "control_download_01/download_status.json"
    control_manifest_path = sample / "control_download_01/download_manifest.json"
    controls_by_target = {}
    for frame in controls_doc["clusters"]:
        target_id = str(frame["target_cluster_id"])
        if target_id in controls_by_target:
            raise ValueError(f"duplicate matched control for target {target_id}")
        controls_by_target[target_id] = frame
    if len(controls_by_target) != 798 or len(controls_doc["clusters"]) != 798:
        raise ValueError("expected 798 unique matched-control fields")
    unknown = set(controls_by_target) - set(target_by_id)
    if unknown:
        raise ValueError(f"control manifest has unknown target ids: {sorted(unknown)[:5]}")
    controls_download_by_id = {str(row["id"]): row for row in control_downloads.values()}
    if len(controls_download_by_id) != 798:
        raise ValueError("control download status does not contain exactly 798 products")
    control_status_sha = control_source_hashes[str(control_status_path)]
    control_manifest_sha = control_source_hashes[str(control_manifest_path)]
    if (sample_doc.get("science_measurement") is not False
            or controls_doc.get("science_measurement") is not False
            or controls_doc.get("controls_ready_for_science") is not False):
        raise ValueError("input research gates are not explicitly false")

    cosmology_values = {"H0": float(sample_doc["h0"]), "Om0": float(sample_doc["omega_m"])}
    if (cosmology_values["H0"] != float(controls_doc["h0"])
            or cosmology_values["Om0"] != float(controls_doc["omega_m"])):
        raise ValueError("target and matched-control cosmologies differ")
    cosmology = Cosmology(H0=cosmology_values["H0"], Om=cosmology_values["Om0"])

    accumulators = {}
    pair_rows = []
    errors = []
    counts = {band: 0 for band in FILTERS}
    exclusions = []
    for index, target in enumerate(targets, 1):
        ident = str(target["id"])
        control = controls_by_target.get(ident)
        if control is None:
            continue
        control_id = str(control["id"])
        download_row = controls_download_by_id.get(control_id)
        if download_row is None:
            raise ValueError(f"{ident}: matched control is absent from completed download status")
        target_band = mask_by_id[ident]["primary_filter"].strip().lower()
        control_band = str(control.get("coadd_header_filter", "")).strip().lower()
        if control_band != str(download_row["review"].get("coadd_header_filter", "")).strip().lower():
            raise ValueError(f"{ident}: control manifest and download-status filter evidence disagree")
        if target_band not in FILTERS or control_band not in FILTERS:
            raise ValueError(f"{ident}: target/control filter metadata is outside HSC-I family")
        if target_band != control_band:
            exclusions.append({"target_id": ident, "control_id": control_id,
                               "reason": "target_control_primary_filter_mismatch",
                               "target_filter": target_band, "control_filter": control_band})
            continue
        try:
            pair, maps = process_pair(
                target, mask_by_id[ident], mask_summary, control, download_row,
                control_base, control_status_sha, cosmology, kpc_max, n_side,
            )
            band = pair["filter"]
            if band not in accumulators:
                shape = maps["delta"].shape
                accumulators[band] = {key: OnlineMap(shape) for key in ("target", "control", "delta")}
            for key in ("target", "control", "delta"):
                accumulators[band][key].add(maps[key], maps["valid"])
            pair_rows.append(pair)
            counts[band] += 1
        except Exception as exc:
            errors.append({"target_id": ident, "control_id": control_id,
                           "error_type": type(exc).__name__, "error": str(exc)})
            raise
        if index % 20 == 0 or index == len(targets):
            progress = {"status": "running", "targets_seen": index,
                        "pairs_complete": len(pair_rows), "by_filter": counts,
                        "elapsed_seconds": time.monotonic() - started}
            write_json(work / "progress.json", progress)
            print(f"Stacked {index}/{len(targets)} targets; pairs={len(pair_rows)}; bands={counts}", flush=True)

    expected_eligible = len(controls_doc["clusters"]) - len(exclusions)
    if errors or len(pair_rows) != expected_eligible or sum(counts.values()) != expected_eligible:
        raise ValueError(f"incomplete target-control stack: {len(pair_rows)}/{expected_eligible} eligible pairs, errors={len(errors)}")
    if any(counts[band] == 0 for band in FILTERS):
        raise ValueError(f"one or more HSC-I-family filters have no frames: {counts}")

    outputs = {}
    axis = np.linspace(-kpc_max, kpc_max, n_side, dtype=float)
    for band in FILTERS:
        band_pairs = [row for row in pair_rows if row["filter"] == band]
        band_acc = accumulators[band]
        outputs[band] = save_band_fits(
            work / f"icl_target_control_stack_{band}.fits", band, axis,
            {"n_pairs": len(band_pairs)}, band_acc, len(band_pairs),
            cosmology_values, kpc_max,
        )
    pair_table_path = work / "paired_radial_profiles.tsv"
    write_pair_table(pair_table_path, pair_rows)
    radial_rows = summarize_profiles(pair_rows)
    radial_path = work / "stacked_radial_profiles.tsv"
    with radial_path.open("x", encoding="utf-8", newline="") as stream:
        fields = list(radial_rows[0])
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t", extrasaction="raise")
        writer.writeheader(); writer.writerows(radial_rows)

    repo_root = HERE.parents[2]
    script_paths = [Path(__file__).resolve(), repo_root / "icl/wgl.py",
                    HERE / "download_pilot.py", HERE / "validate_control_images.py",
                    HERE / "prepare_detection_psf.py"]
    result = {
        "schema": "astrafek-hsc-icl-target-control-2d-stack-v1",
        "status": "complete_exploratory",
        "sample": str(sample.resolve()), "n_targets_available": len(targets),
        "n_matched_controls_available": len(controls_doc["clusters"]),
        "n_pairs_excluded_filter_mismatch": len(exclusions),
        "filter_mismatch_exclusions": exclusions,
        "n_pairs_stacked": len(pair_rows), "n_by_filter": counts,
        "grid": {"coordinate": "proper_kpc", "kpc_max": kpc_max,
                 "n_side": n_side, "pixel_step_kpc": float(axis[1] - axis[0])},
        "stack_definition": {
            "calibration": "HSC FLUXMAG0 to nJy/pixel, then pixel-area division in icl.wgl.warp_metric",
            "mask": "target MASK_K2P5; matched-control SEP/Kron 2.5; native bad/bright-object pixels retained",
            "support": "strict pairwise intersection of fully valid linearly interpolated target and control pixels",
            "weighting": "equal pair weight independently at each grid pixel; no inverse-variance weighting",
            "maps": "target, matched-control, and target-minus-control means; pixelwise pair count and between-pair SEM",
            "radial_profiles": "each matched pair contributes one annular mean per radial bin; between-pair SD/SEM",
        },
        "outputs": {"maps": outputs, "paired_profiles": str(pair_table_path),
                    "paired_profiles_sha256": sha256(pair_table_path),
                    "stacked_profiles": str(radial_path), "stacked_profiles_sha256": sha256(radial_path)},
        "input_sha256": {
            str(sample_manifest_path): sha256(sample_manifest_path),
            str(target_mask_summary_path): sha256(target_mask_summary_path),
            str(target_mask_records_path): sha256(target_mask_records_path),
            str(control_manifest_path): control_manifest_sha,
            str(control_status_path): control_status_sha,
            **{row["target_image"]: row["target_image_sha256"] for row in pair_rows},
            **{row["target_mask"]: row["target_mask_sha256"] for row in pair_rows},
            **{row["control_image"]: row["control_image_sha256"] for row in pair_rows},
        },
        "script_sha256": {str(path): sha256(path) for path in script_paths},
        "software": {"numpy": np.__version__, "python": sys.version.split()[0]},
        "science_measurement": False,
        "limitations": [
            "Exploratory HSC-I-family stack; i and i2 are reported separately.",
            "No extinction correction, residual sky-transfer correction, or PSF-core homogenization is applied.",
            "PSF Picker cores do not constrain extended wings; the stack does not establish a diffuse ICL detection.",
            "Between-pair SEM is not the full uncertainty: it omits correlated coadd noise, sky, PSF, calibration, and astrophysical covariance.",
            "The target/control coadds are assumed independent for paired subtraction; cross-image covariance was not measured.",
            "Pair count varies with position; inspect NPAIR and radial-profile n_pairs before interpreting outer radii.",
            "This stack is not the weak-lensing background-source stack and has no bearing on lensing calibration.",
        ],
        "elapsed_seconds": time.monotonic() - started,
    }
    write_json(work / "summary.json", result)
    (work / "progress.json").unlink(missing_ok=True)
    print(json.dumps({k: result[k] for k in ("status", "n_pairs_stacked", "n_by_filter",
                                              "grid", "science_measurement", "elapsed_seconds")},
                     indent=2), flush=True)
    return result


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", type=Path, default=DEFAULT_SAMPLE)
    parser.add_argument("--work", type=Path, default=DEFAULT_WORK)
    parser.add_argument("--kpc-max", type=float, default=500.0)
    parser.add_argument("--n-side", type=int, default=201)
    args = parser.parse_args(argv)
    if not os.environ.get("SLURM_JOB_ID"):
        parser.error("run the full 800-pair image stack inside a Slurm allocation")
    if not np.isfinite(args.kpc_max) or args.kpc_max <= 0 or args.n_side < 3:
        parser.error("grid size must be positive and n-side at least 3")
    run(args.sample.resolve(), args.work.resolve(), args.kpc_max, args.n_side)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
