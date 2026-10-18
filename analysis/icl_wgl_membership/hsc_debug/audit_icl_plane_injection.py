#!/usr/bin/env python3
"""Measure local-plane transfer of diffuse light in paired HSC target/control images.

This offline audit injects analytic diffuse profiles into the already downloaded
coadd/bg science arrays, before Astrafek's 600--900 kpc tile-plane estimator.
The existing K2.5 cluster masks and re-derived control source masks are held
fixed. The same 74 target/control pairs, filter assignments, annulus and
connected spatial blocks are used throughout. It does not test the HSC pipeline
sky subtraction that occurred before these coadd/bg products were delivered.
"""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import sys
import time

import numpy as np
from astropy.cosmology import FlatLambdaCDM
from astropy.io import fits

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import download_pilot
import prepare_detection_psf
import prestack_core
import validate_control_images
from evaluate_prestack_pair_support import pair_block_labels
from validate_control_images import nJy_factor


DEFAULT_SAMPLE = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
DEFAULT_ASSOCIATION = DEFAULT_SAMPLE / "icl_pixel_extinction_2026-10-19_03/pixel_extinction_satellite_association.tsv"
SHAPES = ((100.0, 1.0), (100.0, 0.65), (250.0, 1.0), (250.0, 0.65))
AMPLITUDE_NJY_ARCSEC2 = 1.0
BIN_INDEX = 3  # 200--300 proper kpc, matching the existing ICL–satellite table.


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def read_tsv(path):
    with Path(path).open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def write_json(path, value):
    path = Path(path)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, ensure_ascii=False, indent=2,
                                    sort_keys=True, allow_nan=False) + "\n",
                         encoding="utf-8")
    temporary.replace(path)


def write_tsv(path, rows):
    if not rows:
        raise ValueError("cannot write an empty injection table")
    path = Path(path)
    temporary = Path(str(path) + ".tmp")
    with temporary.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t",
                                extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)
    temporary.replace(path)


def block_bootstrap_mean(values, labels, replicas, seed):
    values = np.asarray(values, dtype=float)
    labels = np.asarray(labels, dtype=str)
    if len(values) != len(labels) or not np.isfinite(values).all():
        raise ValueError("invalid bootstrap vector or labels")
    blocks = np.unique(labels)
    if len(blocks) < 6:
        raise ValueError("at least six connected spatial blocks are required")
    per_block = {block: np.flatnonzero(labels == block) for block in blocks}
    rng = np.random.default_rng(seed)
    draws = np.empty(replicas, dtype=float)
    for index in range(replicas):
        chosen = rng.choice(blocks, size=len(blocks), replace=True)
        rows = np.concatenate([per_block[block] for block in chosen])
        draws[index] = float(np.mean(values[rows]))
    return {
        "mean": float(np.mean(values)),
        "median": float(np.median(values)),
        "ci95": [float(np.percentile(draws, 2.5)), float(np.percentile(draws, 97.5))],
        "replicas": int(replicas), "n_blocks": int(len(blocks)),
    }


def select_shape_rows(rows, band, scale_kpc, axis_ratio, expected_pairs):
    subset = [row for row in rows if row["filter"] == band
              and float(row["scale_kpc"]) == float(scale_kpc)
              and float(row["axis_ratio"]) == float(axis_ratio)]
    pair_ids = [(row["cluster_id"], row["control_id"]) for row in subset]
    if len(subset) != expected_pairs or len(set(pair_ids)) != expected_pairs:
        raise ValueError("shape group differs from the fixed paired sample")
    return subset


def load_frame_arrays(path, frame, expected_sha):
    if not Path(path).is_file() or sha256(path) != expected_sha:
        raise ValueError(f"source image checksum mismatch: {path}")
    info = download_pilot.validate_fits(path, frame)
    with fits.open(path, memmap=False, checksum=True) as hdul:
        image = np.asarray(hdul[info["science_hdu"]].data, dtype=np.float32)
        variance = np.asarray(hdul[info["variance_hdu"]].data, dtype=np.float32)
        bits = np.asarray(hdul[info["mask_hdu"]].data, dtype=np.int32)
    if not (image.shape == variance.shape == bits.shape):
        raise ValueError("science, variance and mask planes have different shapes")
    conversion = float(nJy_factor(info))
    if not np.isfinite(conversion) or conversion <= 0:
        raise ValueError("invalid FLUXMAG0 surface-brightness conversion")
    return info, image, variance, bits, conversion


def target_mask(mask_record, mask_summary):
    path = Path(mask_record["mask_file"])
    expected = mask_record["mask_sha256"]
    if (mask_summary.get("mask_output_sha256", {}).get(str(path)) != expected
            or not path.is_file() or sha256(path) != expected):
        raise ValueError(f"target mask provenance mismatch: {path}")
    with fits.open(path, memmap=False, checksum=True) as hdul:
        mask = np.asarray(hdul["MASK_K2P5"].data, dtype=bool)
    return mask


def measure_transfer(image, variance, mask, info, conversion, z,
                     scale_kpc, axis_ratio, cosmology,
                     amplitude=AMPLITUDE_NJY_ARCSEC2):
    center = np.asarray(info["center_pixel_xy"], dtype=float)
    kpc_per_pixel = float(cosmology.kpc_proper_per_arcmin(z).value
                          * info["pixel_scale_arcsec"] / 60.0)
    geometry = prestack_core.geometry(image.shape, center, kpc_per_pixel)
    bins = geometry[3]
    valid_annulus = (~mask & np.isfinite(image) & np.isfinite(variance)
                     & (variance > 0) & (bins == BIN_INDEX))
    if valid_annulus.sum() < 10:
        raise ValueError("insufficient unmasked pixels in ICL 200--300 kpc annulus")
    baseline = prestack_core.measure(image, variance, mask, geometry, center, kpc_per_pixel)
    if baseline["n_tiles"] < 10 or not np.isfinite(baseline["corrected"][BIN_INDEX]):
        raise ValueError("baseline has no valid local-plane estimate")
    model = prestack_core.diffuse_model(geometry[0], geometry[1], scale_kpc,
                                        axis_ratio, amplitude).astype(np.float32)
    true_mean = float(np.mean(model[valid_annulus]))
    if not np.isfinite(true_mean) or true_mean <= 0:
        raise ValueError("invalid visible injected annulus mean")
    injected_image = image + model / conversion
    injected = prestack_core.measure(injected_image, variance, mask, geometry, center,
                                     kpc_per_pixel)
    if injected["n_tiles"] < 10 or not np.isfinite(injected["corrected"][BIN_INDEX]):
        raise ValueError("injected frame has no valid local-plane estimate")
    raw_delta = float((injected["raw"][BIN_INDEX] - baseline["raw"][BIN_INDEX]) * conversion)
    corrected_delta = float((injected["corrected"][BIN_INDEX]
                             - baseline["corrected"][BIN_INDEX]) * conversion)
    return {
        "annulus_pixels": int(valid_annulus.sum()),
        "truth_visible_mean_njy_arcsec2": true_mean,
        "raw_delta_njy_arcsec2": raw_delta,
        "plane_corrected_delta_njy_arcsec2": corrected_delta,
        "raw_recovery_fraction": raw_delta / true_mean,
        "plane_corrected_recovery_fraction": corrected_delta / true_mean,
        "baseline_plane_njy_arcsec2_per_unit": (baseline["plane"] * conversion).tolist(),
        "injected_plane_njy_arcsec2_per_unit": (injected["plane"] * conversion).tolist(),
        "n_plane_tiles_baseline": int(baseline["n_tiles"]),
        "n_plane_tiles_injected": int(injected["n_tiles"]),
        "kpc_per_pixel": kpc_per_pixel,
    }


def run(args):
    if not os.environ.get("SLURM_JOB_ID"):
        raise ValueError("image-level injection requires a Slurm allocation on syntax")
    sample, input_path, output = args.sample.resolve(), args.input.resolve(), args.work.resolve()
    if output.exists():
        raise ValueError(f"refusing to overwrite existing output: {output}")
    if args.replicas < 1000:
        raise ValueError("at least 1000 block-bootstrap replicas are required")

    manifests = {
        "sample": sample / "sample_manifest.json",
        "mask_summary": sample / "prestack_cluster_masks_2026-10-16_01/summary.json",
        "mask_records": sample / "prestack_cluster_masks_2026-10-16_01/records.tsv",
        "control_manifest": sample / "control_download_01/download_manifest.json",
        "control_status": sample / "control_download_01/download_status.json",
        "control_extraction": sample / "prestack_extraction_2026-10-15_01/extraction_summary.json",
    }
    for path in manifests.values():
        if not path.is_file():
            raise ValueError(f"required provenance file missing: {path}")
    sample_doc = read_json(manifests["sample"])
    mask_summary = read_json(manifests["mask_summary"])
    control_doc = read_json(manifests["control_manifest"])
    control_status = read_json(manifests["control_status"])
    control_extraction = read_json(manifests["control_extraction"])
    if len(sample_doc.get("clusters", [])) != 800 or mask_summary.get("n_success") != 800:
        raise ValueError("expected the complete exploratory 800-target mask set")
    if mask_summary.get("status") != "complete" or mask_summary.get("mask_ready_for_science") is not False:
        raise ValueError("unexpected target mask science status")
    if (mask_summary.get("output_sha256", {}).get(str(manifests["mask_records"]))
            != sha256(manifests["mask_records"])):
        raise ValueError("target mask records checksum mismatch")
    if control_status.get("status") != "complete" or control_status.get("downloaded") != 798:
        raise ValueError("expected the existing 798 downloaded matched controls")
    if control_extraction.get("status") != "complete" or control_extraction.get("n_measured") != 798:
        raise ValueError("expected the existing 798 control profile extraction")

    association = read_tsv(input_path)
    if len(association) != 74:
        raise ValueError(f"expected the fixed 74-pair sample, got {len(association)}")
    pair_keys = [(str(r["cluster_id"]), r["filter"]) for r in association]
    if len(set(pair_keys)) != 74:
        raise ValueError("duplicate cluster/filter pair in ICL association input")
    target_manifest = {str(r["id"]): r for r in sample_doc["clusters"]}
    mask_records = {str(r["id"]): r for r in read_tsv(manifests["mask_records"])}
    control_frames = {str(r["id"]): r for r in control_doc["clusters"]}
    controls_by_target = {str(r["target_cluster_id"]): r for r in control_extraction["records"]}
    if len(mask_records) != 800 or len(controls_by_target) != 798:
        raise ValueError("target mask/control link records are incomplete or duplicated")
    if control_doc.get("frame_role") != "matched_control" or control_doc.get("bcg_policy") != "none":
        raise ValueError("control manifest is not marked as matched random fields")
    control_downloaded = {str(r["id"]): r for r in control_status.get("completed", [])}
    if len(control_downloaded) != 798:
        raise ValueError("control download completion IDs are incomplete or duplicated")

    cosmo = FlatLambdaCDM(H0=float(sample_doc["h0"]), Om0=float(sample_doc["omega_m"]), Tcmb0=0)
    # Validate block IDs once from the fixed 74 matched target/control coordinates.
    target_coords = [{"ra_deg": float(r["target_ra_deg"]), "dec_deg": float(r["target_dec_deg"])}
                     for r in association]
    control_coords = [{"ra_deg": float(r["control_ra_deg"]), "dec_deg": float(r["control_dec_deg"])}
                      for r in association]
    _, labels, _, _, n_blocks, largest_block = pair_block_labels(
        target_coords, control_coords, np.ones(len(association), dtype=bool))
    output.mkdir(parents=True, exist_ok=False)
    progress_path = output / "progress.json"
    source_hashes = {str(path): sha256(path) for path in (*manifests.values(), input_path)}
    rows = []
    started = time.monotonic()
    write_json(progress_path, {"status": "running", "n_pairs_done": 0, "n_pairs_expected": 74,
                               "slurm_job_id": os.environ["SLURM_JOB_ID"]})
    try:
        for pair_index, pair in enumerate(association):
            ident = str(pair["cluster_id"])
            band = pair["filter"]
            if ident not in target_manifest or ident not in mask_records or ident not in controls_by_target:
                raise ValueError(f"missing target/mask/control for pair {ident}")
            frame = dict(target_manifest[ident])
            mask_record = mask_records[ident]
            control_record = controls_by_target[ident]
            control_id = str(control_record["control_id"])
            control_frame = control_frames.get(control_id)
            if control_frame is None:
                raise ValueError(f"control {control_id} missing from download manifest")
            if (control_id not in control_downloaded
                    or control_downloaded[control_id].get("sha256") != control_record["source_sha256"]):
                raise ValueError(f"control download/extraction provenance differs for {control_id}")
            if (mask_record["primary_filter"] != band or control_record["header_filter"] != band
                    or control_frame["coadd_header_filter"] != band):
                raise ValueError(f"filter mismatch in target/control pair {ident}")
            frame["image"] = mask_record["image"]
            target_path = Path(mask_record["image"])
            target_expected = mask_record["image_sha256"]
            target_mask_array = target_mask(mask_record, mask_summary)
            target_info, target_image, target_var, _, target_conversion = load_frame_arrays(
                target_path, frame, target_expected)
            if target_info["coadd_header_filter"] != band:
                raise ValueError(f"target image filter mismatch for {ident}")
            source_hashes[str(target_path)] = target_expected

            control_path = Path(control_record["source_image"])
            control_expected = control_record["source_sha256"]
            control_info, control_image, control_var, control_bits, control_conversion = load_frame_arrays(
                control_path, control_frame, control_expected)
            if control_info["coadd_header_filter"] != band:
                raise ValueError(f"control image filter mismatch for {control_id}")
            _, control_masks, _, _, _ = validate_control_images.control_masks(
                control_image, control_var, control_bits, control_info["center_pixel_xy"],
                control_info["pixel_scale_arcsec"])
            control_mask = np.asarray(control_masks[prepare_detection_psf.FACTORS.index(2.5)],
                                      dtype=bool)
            source_hashes[str(control_path)] = control_expected

            target_z = float(frame["z"])
            control_z = float(control_frame["z"])
            if not np.isclose(target_z, control_z, rtol=0, atol=1e-12):
                raise ValueError(f"control angular scale differs from paired target {ident}")
            for scale_kpc, axis_ratio in SHAPES:
                shape_name = f"exp{int(scale_kpc)}_q{axis_ratio:g}"
                target_result = measure_transfer(
                    target_image, target_var, target_mask_array, target_info, target_conversion,
                    target_z, scale_kpc, axis_ratio, cosmo, args.amplitude)
                control_result = measure_transfer(
                    control_image, control_var, control_mask, control_info, control_conversion,
                    control_z, scale_kpc, axis_ratio, cosmo, args.amplitude)
                rows.append({
                    "pair_index": pair_index, "cluster_id": ident, "filter": band,
                    "control_id": control_id, "block_id": str(labels[pair_index]),
                    "scale_kpc": scale_kpc, "axis_ratio": axis_ratio,
                    "amplitude_njy_arcsec2_at_100kpc": args.amplitude,
                    **{f"target_{key}": value for key, value in target_result.items()},
                    **{f"control_{key}": value for key, value in control_result.items()},
                })
            if (pair_index + 1) % 5 == 0 or pair_index + 1 == len(association):
                write_json(progress_path, {"status": "running", "n_pairs_done": pair_index + 1,
                                           "n_pairs_expected": 74,
                                           "slurm_job_id": os.environ["SLURM_JOB_ID"],
                                           "elapsed_seconds": time.monotonic() - started})
                print(f"Injected and measured {pair_index + 1}/74 matched fields", flush=True)

        # Paired spatial block bootstrap for the difference in recovery between
        # real cluster and random control fields, for every input shape/filter.
        grouped = {}
        for band in ("i2", "i"):
            band_rows = [r for r in rows if r["filter"] == band]
            expected_pairs = sum(r["filter"] == band for r in association)
            for scale_kpc, axis_ratio in SHAPES:
                shape_name = f"exp{int(scale_kpc)}_q{axis_ratio:g}"
                subset = select_shape_rows(rows, band, scale_kpc, axis_ratio, expected_pairs)
                subset_labels = np.asarray([r["block_id"] for r in subset], dtype=str)
                rng_seed = args.seed + (1000 if band == "i2" else 2000) + int(scale_kpc) + int(axis_ratio * 100)
                metrics = {}
                for estimator in ("raw", "plane_corrected"):
                    target_values = [r[f"target_{estimator}_recovery_fraction"] for r in subset]
                    control_values = [r[f"control_{estimator}_recovery_fraction"] for r in subset]
                    difference = np.asarray(target_values) - np.asarray(control_values)
                    metrics[f"target_{estimator}"] = block_bootstrap_mean(
                        target_values, subset_labels, args.replicas, rng_seed)
                    metrics[f"control_{estimator}"] = block_bootstrap_mean(
                        control_values, subset_labels, args.replicas, rng_seed + 1)
                    metrics[f"target_minus_control_{estimator}"] = block_bootstrap_mean(
                        difference, subset_labels, args.replicas, rng_seed + 2)
                grouped[f"{band}_{shape_name}"] = {"n_pairs": len(subset), "results": metrics}

        rows_path = output / "transfer_rows.tsv"
        write_tsv(rows_path, rows)
        summary = {
            "schema": "hsc-icl-plane-transfer-paired-injection-v1",
            "status": "complete_exploratory", "science_measurement": False,
            "observational_validation_passed": False,
            "n_pairs": len(association), "filter_counts": {band: sum(r["filter"] == band for r in association)
                                                               for band in ("i", "i2")},
            "n_connected_target_control_blocks": int(n_blocks),
            "largest_block_pairs": int(largest_block),
            "method": {
                "injection": "elliptical exponential diffuse profile, fixed angle 37 degrees, observed amplitude at elliptical radius 100 kpc",
                "shapes": [{"scale_kpc": scale, "axis_ratio": q} for scale, q in SHAPES],
                "amplitude_njy_arcsec2": args.amplitude,
                "measurement": "same K2.5 valid-pixel mask and Astrafek 600--900 kpc 64-pixel tile-plane estimator; measure 200--300 kpc annulus before and after injection",
                "control": "the corresponding matched random coadd/bg frame; control source masks are re-derived on the unmodified frame and held fixed for injection",
                "uncertainty": "paired target/control differences resampled over identical connected joint 5-degree sky blocks",
                "bootstrap_replicas": args.replicas, "seed": args.seed,
                "noise_and_psf": "no source Poisson noise or PSF convolution; original diagonal variance held fixed",
                "image_product_limit": "tests Astrafek's downstream local plane estimator on delivered global-sky coadd/bg data; does not test upstream HSC pipeline sky subtraction",
            },
            "response_by_filter_shape": grouped,
            "limitations": [
                "Only the existing 74-pair ICL/satellite support sample is used; inference does not represent the full 800-cluster parent sample.",
                "K2.5 source masks are fixed after injection; source detection, segmentation, and mask growth response are not measured.",
                "Analytic exponential profiles do not span the full observed ICL morphology distribution.",
                "This validates transfer through one downstream plane estimator only; upstream sky subtraction, extended PSF wings, extinction, and correlated coadd noise remain unvalidated.",
                "The result is an exploratory transfer calibration and must not be used to correct the ICL signal or covariance without broader validation.",
            ],
            "input_sha256": source_hashes,
            "code_sha256": {str(Path(__file__).resolve()): sha256(__file__),
                            str(Path(prestack_core.__file__).resolve()): sha256(prestack_core.__file__),
                            str(Path(download_pilot.__file__).resolve()): sha256(download_pilot.__file__),
                            str(Path(prepare_detection_psf.__file__).resolve()): sha256(prepare_detection_psf.__file__),
                            str(Path(validate_control_images.__file__).resolve()): sha256(validate_control_images.__file__),
                            str(Path(validate_control_images.quality.__file__).resolve()): sha256(validate_control_images.quality.__file__)},
            "output_sha256": {str(rows_path): sha256(rows_path)},
            "elapsed_seconds": time.monotonic() - started,
        }
        write_json(output / "summary.json", summary)
        write_json(progress_path, {"status": "complete", "n_pairs_done": 74,
                                   "n_pairs_expected": 74, "slurm_job_id": os.environ["SLURM_JOB_ID"],
                                   "elapsed_seconds": summary["elapsed_seconds"]})
    except BaseException as error:
        write_json(progress_path, {"status": "stopped_on_error", "n_pairs_done": len(rows) // len(SHAPES),
                                   "n_pairs_expected": 74,
                                   "slurm_job_id": os.environ["SLURM_JOB_ID"],
                                   "error_type": type(error).__name__,
                                   "elapsed_seconds": time.monotonic() - started})
        raise
    print(json.dumps({"status": summary["status"], "science_measurement": False,
                      "output": str(output), "response_by_filter_shape": grouped},
                     ensure_ascii=False))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", type=Path, default=DEFAULT_SAMPLE)
    parser.add_argument("--input", type=Path, default=DEFAULT_ASSOCIATION)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--replicas", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=20261019)
    parser.add_argument("--amplitude", type=float, default=AMPLITUDE_NJY_ARCSEC2)
    args = parser.parse_args()
    if not np.isfinite(args.amplitude) or args.amplitude <= 0:
        parser.error("injection amplitude must be positive and finite")
    try:
        run(args)
    except (OSError, ValueError) as error:
        parser.error(str(error))


if __name__ == "__main__":
    main()
