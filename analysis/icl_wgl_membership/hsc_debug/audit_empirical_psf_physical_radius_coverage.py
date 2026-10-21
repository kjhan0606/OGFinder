#!/usr/bin/env python3
"""Compare empirical-star radial support to the cluster ICL measurement scales.

Uses existing manifests, hash-recorded star ledgers and audit summaries only.
It reads no science pixels and submits no HSC requests.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import statistics

import numpy as np
from astropy.cosmology import FlatLambdaCDM


PROFILE_MAX_RADIUS_PIX = 96.0
PROFILE_WING_INNER_RADIUS_PIX = 48.0
STAMP_SIDE_PIX = 401
ICL_INNER_KPC = 200.0
ICL_OUTER_KPC = 300.0
ARCSEC_PER_RADIAN = 206264.80624709636


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def kpc_per_pixel(z: float, pixel_scale_arcsec: float, cosmology: FlatLambdaCDM) -> float:
    z, scale = float(z), float(pixel_scale_arcsec)
    if not math.isfinite(z) or z <= 0 or not math.isfinite(scale) or scale <= 0:
        raise ValueError("redshift and pixel scale must be finite and positive")
    return float(cosmology.angular_diameter_distance(z).to_value("kpc") * scale / ARCSEC_PER_RADIAN)


def inscribed_stamp_radius_pixels(center_x: float, center_y: float, side: int = STAMP_SIDE_PIX) -> float:
    if side < 2 or not np.isfinite([center_x, center_y]).all():
        raise ValueError("invalid stamp geometry")
    upper = side - 1
    if not (0 <= center_x <= upper and 0 <= center_y <= upper):
        raise ValueError("stamp center lies outside the stamp")
    return float(min(center_x, center_y, upper - center_x, upper - center_y))


def _quantiles(values: list[float]) -> dict[str, float]:
    array = np.asarray(values, dtype=float)
    return {key: float(np.quantile(array, q)) for key, q in (
        ("q10", 0.10), ("q50", 0.50), ("q90", 0.90),
    )}


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with Path(path).open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def run(sample_root: Path, ledger_path: Path, profile_summary_path: Path,
        quality_summary_path: Path, output_dir: Path) -> dict:
    sample_root, ledger_path, output_dir = map(Path, (sample_root, ledger_path, output_dir))
    profile_summary_path, quality_summary_path = map(Path, (profile_summary_path, quality_summary_path))
    if output_dir.exists():
        raise FileExistsError(f"refusing to overwrite PSF radius audit: {output_dir}")
    manifest_path = sample_root / "download_01/download_manifest.json"
    status_path = sample_root / "download_01/download_status.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    status = json.loads(status_path.read_text(encoding="utf-8"))
    profile_summary = json.loads(profile_summary_path.read_text(encoding="utf-8"))
    quality_summary = json.loads(quality_summary_path.read_text(encoding="utf-8"))
    ledger = _read_tsv(ledger_path)
    stamped = [row for row in ledger if row.get("stamp", "").strip()]
    if len(stamped) != 333 or len({row["candidate_id"] for row in stamped}) != 333:
        raise ValueError("expected 333 unique complete Gaia-supported stamps")
    cluster_ids = {row["cluster_id"] for row in stamped}
    if len(cluster_ids) != 227:
        raise ValueError("expected 227 HSC frames in the star ledger")
    if (profile_summary.get("n_stamps") != len(stamped)
            or profile_summary.get("n_unique_frames") != len(cluster_ids)
            or profile_summary.get("science_measurement") is not False
            or profile_summary.get("psf_approved") is not False):
        raise ValueError("profile stack summary does not match the audited ledger or status")
    if quality_summary.get("n_stamps") != len(stamped) or quality_summary.get("psf_approved") is not False:
        raise ValueError("quality audit summary does not match the audited ledger")
    profile_inputs = profile_summary.get("input_sha256", {})
    quality_inputs = quality_summary.get("source_hashes", {})
    ledger_abs = str(ledger_path.resolve())
    ledger_digest = sha256(ledger_path)
    if profile_inputs.get(ledger_abs) != ledger_digest or quality_inputs.get(ledger_abs) != ledger_digest:
        raise ValueError("profile/quality provenance does not match the star ledger")
    for row in stamped:
        stamp = str(Path(row["stamp"]).resolve())
        expected = row["stamp_sha256"]
        if profile_inputs.get(stamp) != expected or quality_inputs.get(stamp) != expected:
            raise ValueError(f"profile/quality provenance does not match star ledger: {stamp}")

    h0, omega_m = float(manifest["h0"]), float(manifest["omega_m"])
    cosmology = FlatLambdaCDM(H0=h0, Om0=omega_m, Tcmb0=0)
    frames = {str(row["id"]): row for row in manifest["clusters"]}
    reviews = {str(row["id"]): row["review"] for row in status["completed"]}
    if len(frames) != 800 or set(frames) != set(reviews):
        raise ValueError("target image manifests do not cover the same 800 cluster IDs")
    profile_max = float(PROFILE_MAX_RADIUS_PIX)
    rows = []
    for star in stamped:
        cluster_id = str(star["cluster_id"])
        frame, review = frames.get(cluster_id), reviews.get(cluster_id)
        if frame is None or review is None:
            raise ValueError(f"star frame is absent from the verified 800-image sample: {cluster_id}")
        if frame.get("band") != "HSC-I" or review.get("filter") != "HSC-I":
            raise ValueError(f"PSF star frame is outside the HSC-I sample: {cluster_id}")
        z = float(frame["z"])
        scale = float(review["pixel_scale_arcsec"])
        kpp = kpc_per_pixel(z, scale, cosmology)
        x = float(star["native_x"]) - int(star["x0"])
        y = float(star["native_y"]) - int(star["y0"])
        stamp_radius_pix = inscribed_stamp_radius_pixels(x, y)
        rows.append({
            "candidate_id": star["candidate_id"], "cluster_id": cluster_id,
            "z_lens": z, "pixel_scale_arcsec": scale, "kpc_per_pixel": kpp,
            "profile_wing_inner_radius_pix": PROFILE_WING_INNER_RADIUS_PIX,
            "profile_max_radius_pix": profile_max,
            "profile_max_radius_kpc": profile_max * kpp,
            "profile_48_96_kpc_min": PROFILE_WING_INNER_RADIUS_PIX * kpp,
            "profile_48_96_kpc_max": profile_max * kpp,
            "stamp_inscribed_radius_pix": stamp_radius_pix,
            "stamp_inscribed_radius_kpc": stamp_radius_pix * kpp,
            "icl_200kpc_radius_pix": ICL_INNER_KPC / kpp,
            "icl_300kpc_radius_pix": ICL_OUTER_KPC / kpp,
            "profile_reaches_50kpc": profile_max * kpp >= 50.0,
            "profile_reaches_75kpc": profile_max * kpp >= 75.0,
            "profile_reaches_200kpc": profile_max * kpp >= ICL_INNER_KPC,
            "stamp_reaches_200kpc": stamp_radius_pix * kpp >= ICL_INNER_KPC,
            "stamp_reaches_300kpc": stamp_radius_pix * kpp >= ICL_OUTER_KPC,
        })

    frame_rows = {}
    for row in rows:
        frame_rows.setdefault(row["cluster_id"], []).append(row)
    profile_radii_frame = [max(r["profile_max_radius_kpc"] for r in v) for v in frame_rows.values()]
    stamp_radii_frame = [max(r["stamp_inscribed_radius_kpc"] for r in v) for v in frame_rows.values()]
    output_dir.mkdir(parents=True, exist_ok=False)
    table_path = output_dir / "psf_physical_radius_coverage.tsv"
    with table_path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)

    input_hashes = {
        str(path.resolve()): sha256(path) for path in (
            manifest_path, status_path, ledger_path, profile_summary_path, quality_summary_path,
        )
    }
    summary = {
        "schema": "astrafek-empirical-psf-physical-radius-coverage-v1",
        "status": "complete_physical_radius_coverage_audit",
        "n_stamps": len(rows), "n_frames": len(frame_rows),
        "cosmology": {"H0_km_s_Mpc": h0, "Omega_m": omega_m, "Tcmb0_K": 0.0},
        "profile_max_radius_pix": profile_max,
        "profile_max_radius_kpc_by_frame": _quantiles(profile_radii_frame),
        "profile_48_96_kpc_by_stamp": {
            "inner_radius": _quantiles([row["profile_48_96_kpc_min"] for row in rows]),
            "outer_radius": _quantiles([row["profile_48_96_kpc_max"] for row in rows]),
        },
        "stamp_inscribed_radius_kpc_by_frame": _quantiles(stamp_radii_frame),
        "profile_reaches_50kpc_stamps": sum(row["profile_reaches_50kpc"] for row in rows),
        "profile_reaches_75kpc_stamps": sum(row["profile_reaches_75kpc"] for row in rows),
        "profile_reaches_200kpc_stamps": sum(row["profile_reaches_200kpc"] for row in rows),
        "profile_reaches_200kpc_frames": sum(value >= ICL_INNER_KPC for value in profile_radii_frame),
        "stamp_reaches_200kpc_stamps": sum(row["stamp_reaches_200kpc"] for row in rows),
        "stamp_reaches_200kpc_frames": sum(value >= ICL_INNER_KPC for value in stamp_radii_frame),
        "stamp_reaches_300kpc_stamps": sum(row["stamp_reaches_300kpc"] for row in rows),
        "stamp_reaches_300kpc_frames": sum(value >= ICL_OUTER_KPC for value in stamp_radii_frame),
        "icl_annulus_kpc": [ICL_INNER_KPC, ICL_OUTER_KPC],
        "science_measurement": False, "psf_approved": False,
        "input_sha256": input_hashes,
        "table_sha256": sha256(table_path),
        "audit_script_sha256": sha256(Path(__file__)),
        "science_image_pixels_read": False, "hsc_requests_made": False,
        "limitations": [
            "This is angular-to-proper-radius coverage only, not an empirical PSF-wing measurement.",
            "The measured profile stack ends at 96 pixels although each source stamp is 401x401 pixels.",
            "Stamp radius is the inscribed circular radius, not corner reach; it still does not test usable masked-flux support.",
            "A radial range reaching a science aperture does not certify PSF matching, sky transfer, or unbiased recovery.",
        ],
    }
    summary_path = output_dir / "summary.json"
    summary_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-root", type=Path, required=True)
    parser.add_argument("--ledger", type=Path, required=True)
    parser.add_argument("--profile-summary", type=Path, required=True)
    parser.add_argument("--quality-summary", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = run(args.sample_root, args.ledger, args.profile_summary, args.quality_summary, args.output)
    print(json.dumps({key: result[key] for key in (
        "status", "n_stamps", "n_frames", "profile_max_radius_kpc_by_frame",
        "stamp_inscribed_radius_kpc_by_frame", "profile_reaches_200kpc_frames",
        "stamp_reaches_200kpc_frames", "stamp_reaches_300kpc_frames",
        "science_measurement", "psf_approved",
    )}, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
