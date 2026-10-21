#!/usr/bin/env python3
"""Sweep outer-field Gaia-star support thresholds without reading science pixels."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from astropy.cosmology import FlatLambdaCDM


BASE = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
DEFAULT_THRESHOLDS_KPC = (500.0, 600.0, 700.0, 800.0, 900.0, 1000.0)
ARCSEC_PER_RADIAN = 206264.80624709636


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with Path(path).open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def kpc_per_pixel(z: float, pixel_scale_arcsec: float, cosmology) -> float:
    if not np.isfinite([z, pixel_scale_arcsec]).all() or z <= 0 or pixel_scale_arcsec <= 0:
        raise ValueError("redshift and pixel scale must be finite and positive")
    return float(cosmology.angular_diameter_distance(z).to_value("kpc")
                 * pixel_scale_arcsec / ARCSEC_PER_RADIAN)


def edge_radius(shape: tuple[int, int], x: float, y: float) -> float:
    height, width = shape
    if not np.isfinite([x, y]).all() or not (0 <= x < width and 0 <= y < height):
        return -1.0
    return float(min(x, y, width - 1 - x, height - 1 - y))


def summarize_thresholds(candidates: list[dict], thresholds_kpc=DEFAULT_THRESHOLDS_KPC) -> list[dict]:
    thresholds = tuple(float(value) for value in thresholds_kpc)
    if not thresholds or any(not np.isfinite(value) or value < 300 for value in thresholds):
        raise ValueError("outer-star thresholds must be finite and at least 300 kpc")
    if len(set(thresholds)) != len(thresholds):
        raise ValueError("outer-star thresholds must be unique")
    rows = []
    for threshold in thresholds:
        selected = [row for row in candidates if row["eligible"] and row["cluster_centric_radius_kpc"] > threshold]
        distances = np.asarray([row["cluster_centric_radius_kpc"] for row in selected], dtype=float)
        frame_counts: dict[str, int] = {}
        for row in selected:
            frame_counts[row["cluster_id"]] = frame_counts.get(row["cluster_id"], 0) + 1
        rows.append({
            "minimum_cluster_centric_radius_kpc": threshold,
            "n_star_control_pairs": len(selected),
            "n_independent_target_frames": len(frame_counts),
            "n_frames_with_multiple_stars": sum(count > 1 for count in frame_counts.values()),
            "max_stars_in_one_frame": max(frame_counts.values(), default=0),
            "star_radius_q10_kpc": float(np.quantile(distances, 0.10)) if len(distances) else None,
            "star_radius_median_kpc": float(np.median(distances)) if len(distances) else None,
            "star_radius_q90_kpc": float(np.quantile(distances, 0.90)) if len(distances) else None,
            "minimum_clearance_inside_300kpc_profile_kpc": float(np.min(distances) - 300.0) if len(distances) else None,
            "same_frame_candidates_are_not_independent": True,
        })
    return rows


def load_candidates(sample: Path) -> tuple[list[dict], dict]:
    target_dir = sample / "download_01"
    control_dir = sample / "control_download_01"
    ledger_path = sample / "empirical_psf_gaia_full_sample_2026-10-20_02/star_ledger.tsv"
    quality_path = sample / "empirical_psf_star_stamp_quality_audit_2026-10-20_02/star_stamp_quality.tsv"
    support_path = sample / "empirical_psf_matched_control_support_2026-10-20_01/matched_control_support.tsv"
    paths = [target_dir / "download_manifest.json", target_dir / "download_status.json",
             control_dir / "download_manifest.json", control_dir / "download_status.json",
             ledger_path, quality_path, support_path]
    if any(not path.is_file() for path in paths):
        raise FileNotFoundError("existing HSC image or PSF-support audit inputs are missing")
    target_manifest, target_status, control_manifest, control_status = [
        json.loads(path.read_text(encoding="utf-8")) for path in paths[:4]
    ]
    if target_manifest.get("band") != "HSC-I" or control_manifest.get("band") != "HSC-I":
        raise ValueError("expected HSC-I target and matched-control manifests")
    targets = {str(row["id"]): row for row in target_manifest["clusters"]}
    target_done = {str(row["id"]): row for row in target_status["completed"]}
    controls = {str(row["id"]): row for row in control_manifest["clusters"]}
    control_done = {str(row["id"]): row for row in control_status["completed"]}
    if len(targets) != 800 or len(target_done) != 800 or len(controls) != 798 or len(control_done) != 798:
        raise ValueError("expected the verified 800-target/798-control HSC-I sample")
    ledger, quality_rows, support_rows = [read_tsv(path) for path in (ledger_path, quality_path, support_path)]
    stamped = [row for row in ledger if row.get("stamp", "").strip()]
    quality = {str(row["candidate_id"]): row for row in quality_rows}
    support = {str(row["candidate_id"]): row for row in support_rows}
    if len(ledger) != 480 or len(stamped) != 333 or len(quality) != 333 or len(support) != 333:
        raise ValueError("unexpected coverage in saved Gaia star, quality, or control-support tables")
    cosmology = FlatLambdaCDM(H0=float(target_manifest["h0"]), Om0=float(target_manifest["omega_m"]), Tcmb0=0)
    candidates = []
    stage_counts = {key: 0 for key in ("stamped", "not_flagged_variable", "target_300kpc_geometry",
        "same_filter_control", "control_300kpc_geometry", "central_mask_quality", "fully_eligible")}
    for star in stamped:
        stage_counts["stamped"] += 1
        if sha256(Path(star["stamp"])) != star["stamp_sha256"]:
            raise ValueError(f"star stamp checksum changed: {star['candidate_id']}")
        if str(star.get("gaia_variability", "")).upper() == "VARIABLE":
            continue
        stage_counts["not_flagged_variable"] += 1
        cluster_id = str(star["cluster_id"])
        target = targets[cluster_id]
        target_record = target_done[cluster_id]
        target_review = target_record["review"]
        scale = kpc_per_pixel(float(target["z"]), float(target_review["pixel_scale_arcsec"]), cosmology)
        x, y = float(star["native_x"]), float(star["native_y"])
        center_x, center_y = map(float, target_review["center_pixel_xy"])
        cluster_radius = math.hypot(x - center_x, y - center_y) * scale
        if edge_radius(tuple(target_review["shape"]), x, y) * scale < 300.0:
            continue
        stage_counts["target_300kpc_geometry"] += 1
        pair = support[str(star["candidate_id"])]
        if pair["pair_status"] != "same_filter_offset_available" or pair["target_filter"] != pair["control_filter"]:
            continue
        stage_counts["same_filter_control"] += 1
        control = controls[str(pair["control_id"])]
        control_record = control_done[str(pair["control_id"])]
        if str(control.get("target_cluster_id")) != cluster_id:
            raise ValueError(f"matched control maps to wrong target: {star['candidate_id']}")
        if control_record["review"]["coadd_header_filter"] != target_review["coadd_header_filter"]:
            raise ValueError(f"matched control filter changed: {star['candidate_id']}")
        if edge_radius(tuple(control_record["review"]["shape"]), float(pair["control_x"]),
                       float(pair["control_y"])) * scale < 300.0:
            continue
        stage_counts["control_300kpc_geometry"] += 1
        q = quality[str(star["candidate_id"])]
        eligible = (int(q["detection_flags"]) == 0
                    and int(q["hsc_bright_object_pixels_core_0_12"]) == 0
                    and float(q["mask_k2p5_core_0_12_open_fraction"]) >= 0.99)
        if eligible:
            stage_counts["central_mask_quality"] += 1
        candidates.append({
            "candidate_id": str(star["candidate_id"]), "cluster_id": cluster_id,
            "cluster_centric_radius_kpc": cluster_radius, "eligible": eligible,
        })
        stage_counts["fully_eligible"] += int(eligible)
    hashes = {str(path.resolve()): sha256(path) for path in paths}
    hashes[str(ledger_path.resolve())] = sha256(ledger_path)
    return candidates, {"stage_counts": stage_counts, "input_sha256": hashes,
        "n_stamped_stars": len(stamped), "n_unique_stamped_frames": len({row["cluster_id"] for row in stamped}),
        "cosmology": {"H0_km_s_Mpc": float(target_manifest["h0"]), "Omega_m": float(target_manifest["omega_m"]), "Tcmb0_K": 0.0}}


def run(sample: Path, output_dir: Path, thresholds_kpc=DEFAULT_THRESHOLDS_KPC) -> dict:
    sample, output_dir = Path(sample).resolve(), Path(output_dir).resolve()
    if output_dir.exists():
        raise FileExistsError(f"refusing to overwrite output directory: {output_dir}")
    candidates, provenance = load_candidates(sample)
    sweep = summarize_thresholds(candidates, thresholds_kpc)
    result = {
        "schema": "astrafek-empirical-psf-outer-star-support-sweep-v1",
        "status": "complete_pixel_free_support_audit",
        "science_measurement": False,
        "psf_approved": False,
        "pixel_data_read": False,
        "hsc_query_submitted": False,
        "data_downloaded": False,
        **provenance,
        "threshold_sweep": sweep,
        "interpretation": "Counts treat matched target/control star pairs as candidate measurements, but independent resampling must group all stars from one target coadd frame.",
        "limitations": [
            "This audit measures selection support only; it does not read pixels or estimate PSF wings.",
            "The 500 kpc threshold places the inner edge of a 300 kpc star profile at the 200 kpc ICL annulus boundary; it does not prove zero cluster-correlated background.",
            "Gaia epoch matching, outer source masks, local sky, and matched-control pixel residuals remain unresolved.",
            "science_measurement remains false and psf_approved remains false.",
        ],
    }
    output_dir.mkdir(parents=True, exist_ok=False)
    result["analysis_script_sha256"] = sha256(Path(__file__).resolve())
    out_path = output_dir / "support_sweep.json"
    with out_path.open("x", encoding="utf-8") as stream:
        json.dump(result, stream, ensure_ascii=False, indent=2, allow_nan=False)
        stream.write("\n")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample", type=Path, default=BASE)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--thresholds-kpc", type=float, nargs="+", default=DEFAULT_THRESHOLDS_KPC)
    args = parser.parse_args()
    result = run(args.sample, args.output_dir, args.thresholds_kpc)
    print(json.dumps({"status": result["status"], "stage_counts": result["stage_counts"],
        "threshold_sweep": result["threshold_sweep"], "output": str(args.output_dir.resolve()),
        "science_measurement": False, "psf_approved": False}, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
