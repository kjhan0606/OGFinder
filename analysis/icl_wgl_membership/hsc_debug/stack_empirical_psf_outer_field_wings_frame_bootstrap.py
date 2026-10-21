#!/usr/bin/env python3
"""Measure matched Gaia-star wings with frame-level resampling for multi-star frames."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np
from astropy.cosmology import FlatLambdaCDM

import stack_empirical_psf_outer_field_wings as core


DEFAULT_SAMPLE = core.BASE
DEFAULT_OUTER_RADIUS_KPC = 600.0


def select_candidates(sample: Path, *, outer_radius_kpc: float) -> tuple[list[dict], dict]:
    """Same physical/quality gates as the 800-kpc pilot, allowing multiple stars per frame."""
    sample = Path(sample).resolve()
    target_dir, control_dir = sample / "download_01", sample / "control_download_01"
    target_manifest_path = target_dir / "download_manifest.json"
    target_status_path = target_dir / "download_status.json"
    control_manifest_path = control_dir / "download_manifest.json"
    control_status_path = control_dir / "download_status.json"
    ledger_path = sample / "empirical_psf_gaia_full_sample_2026-10-20_02/star_ledger.tsv"
    quality_path = sample / "empirical_psf_star_stamp_quality_audit_2026-10-20_02/star_stamp_quality.tsv"
    support_path = sample / "empirical_psf_matched_control_support_2026-10-20_01/matched_control_support.tsv"
    source_paths = (target_manifest_path, target_status_path, control_manifest_path,
                    control_status_path, ledger_path, quality_path, support_path)
    if not all(path.is_file() for path in source_paths):
        raise FileNotFoundError("one or more existing HSC/PSF input ledgers are missing")
    target_manifest = json.loads(target_manifest_path.read_text(encoding="utf-8"))
    target_status = json.loads(target_status_path.read_text(encoding="utf-8"))
    control_manifest = json.loads(control_manifest_path.read_text(encoding="utf-8"))
    control_status = json.loads(control_status_path.read_text(encoding="utf-8"))
    if target_manifest.get("science_measurement") is not False or control_manifest.get("science_measurement") is not False:
        raise ValueError("input manifests have unexpected science-measurement status")
    cosmology = FlatLambdaCDM(H0=float(target_manifest["h0"]), Om0=float(target_manifest["omega_m"]), Tcmb0=0)
    targets = {str(row["id"]): row for row in target_manifest["clusters"]}
    target_done = {str(row["id"]): row for row in target_status["completed"]}
    controls = {str(row["id"]): row for row in control_manifest["clusters"]}
    control_done = {str(row["id"]): row for row in control_status["completed"]}
    ledger = core.read_tsv(ledger_path)
    quality = {str(row["candidate_id"]): row for row in core.read_tsv(quality_path)}
    support = {str(row["candidate_id"]): row for row in core.read_tsv(support_path)}
    if len(targets) != 800 or len(target_done) != 800 or len(controls) != 798 or len(control_done) != 798:
        raise ValueError("the pinned 800-target/798-control data set is incomplete")
    stamped = [row for row in ledger if row.get("stamp")]
    if len(ledger) != 480 or len(stamped) != 333 or len(quality) != 333 or len(support) != 333:
        raise ValueError("unexpected Gaia star-ledger, quality, or control-support cardinality")
    counts = {key: 0 for key in ("complete_stamps", "not_flagged_variable", "outer_field",
        "target_300kpc_footprint", "same_filter_control", "control_300kpc_footprint", "k2p5_core_99pct", "selected")}
    selected = []
    for star in stamped:
        counts["complete_stamps"] += 1
        if core.sha(Path(star["stamp"])) != star["stamp_sha256"]:
            raise ValueError(f"star-stamp hash mismatch: {star['candidate_id']}")
        if str(star.get("gaia_variability", "")).upper() == "VARIABLE":
            continue
        counts["not_flagged_variable"] += 1
        cid = str(star["cluster_id"])
        target, target_record = targets[cid], target_done[cid]
        target_review = target_record["review"]
        z = float(target["z"])
        pixel_scale = float(target_review["pixel_scale_arcsec"])
        kpp = core.kpc_per_pixel(z, pixel_scale, cosmology)
        x, y = float(star["native_x"]), float(star["native_y"])
        cx, cy = map(float, target_review["center_pixel_xy"])
        cluster_radius_kpc = float(np.hypot(x - cx, y - cy) * kpp)
        if cluster_radius_kpc <= outer_radius_kpc:
            continue
        counts["outer_field"] += 1
        target_shape = tuple(int(v) for v in target_review["shape"])
        if core.edge_radius(target_shape, x, y) * kpp < 300.0:
            continue
        counts["target_300kpc_footprint"] += 1
        pair = support[str(star["candidate_id"])]
        if (pair["pair_status"] != "same_filter_offset_available"
                or pair["target_filter"] != pair["control_filter"]):
            continue
        counts["same_filter_control"] += 1
        control_id = str(pair["control_id"])
        control, control_record = controls[control_id], control_done[control_id]
        if control.get("target_cluster_id") != cid:
            raise ValueError(f"matched control belongs to a different target for {star['candidate_id']}")
        control_review = control_record["review"]
        if control_review["coadd_header_filter"] != target_review["coadd_header_filter"]:
            raise ValueError(f"target/control filter mismatch for {star['candidate_id']}")
        if abs(float(control_review["pixel_scale_arcsec"]) - pixel_scale) > 1e-5:
            raise ValueError(f"target/control pixel scales differ for {star['candidate_id']}")
        control_x, control_y = float(pair["control_x"]), float(pair["control_y"])
        control_shape = tuple(int(v) for v in control_review["shape"])
        if core.edge_radius(control_shape, control_x, control_y) * kpp < 300.0:
            continue
        counts["control_300kpc_footprint"] += 1
        q = quality[str(star["candidate_id"])]
        if (int(q["detection_flags"]) != 0 or int(q["hsc_bright_object_pixels_core_0_12"]) != 0
                or float(q["mask_k2p5_core_0_12_open_fraction"]) < 0.99):
            continue
        counts["k2p5_core_99pct"] += 1
        selected.append({
            "candidate_id": str(star["candidate_id"]), "cluster_id": cid,
            "detection_number": int(star["number"]), "source_id": str(star["source_id"]),
            "gaia_variability": str(star["gaia_variability"]), "epoch_review_required": True,
            "target_image": str((target_dir / "images" / target_record["file"]).resolve()),
            "control_id": control_id,
            "control_image": str((control_dir / "images" / control_record["file"]).resolve()),
            "target_x": x, "target_y": y, "control_x": control_x, "control_y": control_y,
            "z_lens": z, "cluster_centric_radius_kpc": cluster_radius_kpc,
            "pixel_scale_arcsec": pixel_scale, "kpc_per_pixel": kpp,
            "target_filter": str(target_review["coadd_header_filter"]),
            "k2p5_core_open_fraction": float(q["mask_k2p5_core_0_12_open_fraction"]),
            "target_shape": target_shape, "control_shape": control_shape,
            "target_status_row": target_record, "control_status_row": control_record,
        })
    counts["selected"] = len(selected)
    if not selected:
        raise ValueError("no matched target/control Gaia stars passed the geometry gates")
    provenance = {
        "sample_root": str(sample),
        "cosmology": {"H0_km_s_Mpc": float(target_manifest["h0"]),
                      "Omega_m": float(target_manifest["omega_m"]), "Tcmb0_K": 0.0},
        "outer_radius_kpc": float(outer_radius_kpc), "counts": counts,
        "n_unique_target_frames": len({row["cluster_id"] for row in selected}),
        "input_sha256": {str(path.resolve()): core.sha(path) for path in source_paths},
    }
    return selected, provenance


def stack_frame_pairs(rows: list[dict], *, n_bootstrap: int, seed: int) -> list[dict]:
    """Average multiple stars within a frame, then bootstrap independent frames."""
    per_star: dict[tuple, dict[str, dict[str, float]]] = {}
    for row in rows:
        if row["status"] != "ok" or row["profile_per_target_core_per_pixel"] is None:
            continue
        key = (str(row["filter"]), int(row["sky_r_min_pix"]), int(row["sky_r_max_pix"]),
               float(row["r_min_kpc"]), float(row["r_max_kpc"]), str(row["cluster_id"]))
        per_star.setdefault(key, {}).setdefault(str(row["candidate_id"]), {})[row["role"]] = float(row["profile_per_target_core_per_pixel"])
    frame_values: dict[tuple, dict[str, list[float]]] = {}
    for (band, sky_lo, sky_hi, r_min, r_max, frame), candidates in per_star.items():
        paired = [value["target"] - value["control"] for value in candidates.values()
                  if "target" in value and "control" in value]
        if paired:
            frame_values.setdefault((band, sky_lo, sky_hi, r_min, r_max), {})[frame] = paired
    all_keys = sorted({(str(row["filter"]), int(row["sky_r_min_pix"]), int(row["sky_r_max_pix"]),
                        float(row["r_min_kpc"]), float(row["r_max_kpc"])) for row in rows})
    rng = np.random.default_rng(seed)
    output = []
    for key in all_keys:
        grouped = frame_values.get(key, {})
        frames = sorted(grouped)
        frame_means = np.asarray([np.mean(grouped[frame]) for frame in frames], dtype=float)
        n_frames, n_pairs = len(frames), sum(len(grouped[frame]) for frame in frames)
        if n_frames and n_bootstrap:
            draws = rng.integers(0, n_frames, size=(n_bootstrap, n_frames))
            boot = np.mean(frame_means[draws], axis=1)
            q16, q84 = map(float, np.quantile(boot, (0.16, 0.84)))
        else:
            q16 = q84 = None
        output.append({
            "filter": key[0], "sky_r_min_pix": key[1], "sky_r_max_pix": key[2],
            "r_min_kpc": key[3], "r_max_kpc": key[4], "r_mid_kpc": (key[3] + key[4]) / 2,
            "n_complete_star_pairs": n_pairs, "n_independent_target_frames": n_frames,
            "n_frames_with_multiple_complete_stars": sum(len(grouped[frame]) > 1 for frame in frames),
            "frame_equal_paired_difference_mean": float(np.mean(frame_means)) if n_frames else None,
            "frame_equal_paired_difference_median": float(np.median(frame_means)) if n_frames else None,
            "bootstrap_q16": q16, "bootstrap_q84": q84,
            "bootstrap_unit": "target/control HSC coadd frame; stars averaged within frame",
            "science_measurement": False, "psf_approved": False,
        })
    return output


def run(sample: Path, output_dir: Path, *, outer_radius_kpc=DEFAULT_OUTER_RADIUS_KPC,
        n_bootstrap=10000, seed=20261024):
    sample, output_dir = Path(sample).resolve(), Path(output_dir).resolve()
    if output_dir.exists():
        raise FileExistsError(f"refusing to overwrite output directory: {output_dir}")
    if n_bootstrap < 1000:
        raise ValueError("at least 1000 frame-bootstrap replicates are required")
    candidates, provenance = select_candidates(sample, outer_radius_kpc=outer_radius_kpc)
    ledger_summary_path = sample / "empirical_psf_gaia_full_sample_2026-10-20_02/summary.json"
    ledger_summary = json.loads(ledger_summary_path.read_text(encoding="utf-8"))
    detection_dir = Path(ledger_summary["detection_dir"])
    rows = []
    hashes: dict[str, str] = dict(provenance["input_sha256"])
    hashes[str(ledger_summary_path.resolve())] = core.sha(ledger_summary_path)
    quality_summary_path = sample / "empirical_psf_star_stamp_quality_audit_2026-10-20_02/summary.json"
    hashes[str(quality_summary_path.resolve())] = core.sha(quality_summary_path)
    # Full image SHA-256 is expensive; compute it once per frame rather than once per star.
    verified_image_hashes = {}
    original_sha = core.sha
    def cached_sha(path):
        key = str(Path(path).resolve())
        if key not in verified_image_hashes:
            verified_image_hashes[key] = original_sha(path)
        return verified_image_hashes[key]
    core.sha = cached_sha
    try:
        for candidate in candidates:
            measured = core.measure_candidate(candidate, sample, detection_dir)
            rows.extend(measured)
            hashes[candidate["target_image"]] = candidate["target_status_row"]["sha256"]
            hashes[candidate["control_image"]] = candidate["control_status_row"]["sha256"]
            frame_dir = detection_dir / f"cluster_{candidate['cluster_id']}"
            for path in (frame_dir / "detections.tsv", frame_dir / "masks.fits.gz"):
                hashes[str(path.resolve())] = original_sha(path)
            print(f"measured candidate={candidate['candidate_id']} frame={candidate['cluster_id']}", flush=True)
    finally:
        core.sha = original_sha
    stack = stack_frame_pairs(rows, n_bootstrap=n_bootstrap, seed=seed)
    output_dir.mkdir(parents=True, exist_ok=False)
    selected = [{key: value for key, value in row.items() if key not in ("target_status_row", "control_status_row")}
                for row in candidates]
    core.write_tsv(output_dir / "selected_candidates.tsv", selected)
    core.write_tsv(output_dir / "matched_radial_profiles.tsv", rows)
    core.write_tsv(output_dir / "frame_bootstrap_paired_stack.tsv", stack)
    summary = {
        "schema": "astrafek-empirical-psf-outer-field-frame-bootstrap-v1",
        "status": "complete_exploratory_diagnostic",
        "counts": provenance["counts"], "n_unique_target_frames": provenance["n_unique_target_frames"],
        "n_profile_rows": len(rows), "n_bootstrap_draws": n_bootstrap, "seed": seed,
        "outer_radius_kpc": float(outer_radius_kpc),
        "target_filter_star_counts": {band: sum(row["target_filter"] == band for row in candidates)
                                      for band in sorted({row["target_filter"] for row in candidates})},
        "radial_edges_kpc": core.RADIAL_EDGES_KPC.tolist(),
        "sky_annuli_pix": [list(item) for item in core.SKY_ANNULI_PIX],
        "input_sha256": hashes,
        "output_sha256": {name: core.sha(output_dir / name) for name in (
            "selected_candidates.tsv", "matched_radial_profiles.tsv", "frame_bootstrap_paired_stack.tsv")},
        "bootstrap_unit": "target/control HSC coadd frame; multiple stars averaged within frame",
        "science_measurement": False, "psf_approved": False,
        "limitations": [
            f"The {outer_radius_kpc:g} kpc outer-star threshold keeps the measured 0-100 kpc satellite-scale stellar profile outside the 200-300 kpc ICL annulus; the 200-300 kpc stellar wing may overlap the cluster ICL when the threshold is below 600 kpc.",
            "Gaia epoch matching remains unresolved; UNKNOWN/NOT_AVAILABLE variability flags are not evidence of non-variability.",
            "Source masks, local-sky choices, target/control residual differences, and coadd covariance remain systematic limitations.",
            "This is a PSF-wing diagnostic only; no ICL correction or science measurement is produced.",
        ],
    }
    summary_path = output_dir / "summary.json"
    with summary_path.open("x", encoding="utf-8") as stream:
        json.dump(summary, stream, ensure_ascii=False, indent=2, allow_nan=False)
        stream.write("\n")
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", type=Path, default=DEFAULT_SAMPLE)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--outer-radius-kpc", type=float, default=DEFAULT_OUTER_RADIUS_KPC)
    parser.add_argument("--bootstrap-draws", type=int, default=10000)
    parser.add_argument("--seed", type=int, default=20261024)
    args = parser.parse_args()
    result = run(args.sample, args.output_dir, outer_radius_kpc=args.outer_radius_kpc,
                 n_bootstrap=args.bootstrap_draws, seed=args.seed)
    print(json.dumps({"status": result["status"], "counts": result["counts"],
        "n_unique_target_frames": result["n_unique_target_frames"],
        "output_sha256": result["output_sha256"], "science_measurement": False,
        "psf_approved": False}, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
