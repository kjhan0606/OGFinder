#!/usr/bin/env python3
"""Audit the radial support and normalization of approved HSC extended PSFs.

This is a local-file audit only. It never downloads data and never approves a
PSF correction. FITS input is expected to be the authors' two-HDU 4001-square
PDR3 reference model; the image is read from extension 1.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import re
import tempfile

import numpy as np
from astropy.io import fits
from astropy.cosmology import FlatLambdaCDM


DEFAULT_PIXEL_SCALE_ARCSEC = 0.168
DEFAULT_RADIAL_EDGES_PIX = np.concatenate(([0.0], np.geomspace(1.0, 2000.0, 49)))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def analyze_image(data: np.ndarray, pixel_scale_arcsec: float = DEFAULT_PIXEL_SCALE_ARCSEC,
                  radial_edges_pix: np.ndarray = DEFAULT_RADIAL_EDGES_PIX):
    """Return annular profile rows and geometry/normalization diagnostics."""
    image = np.asarray(data, dtype=np.float64)
    edges = np.asarray(radial_edges_pix, dtype=np.float64)
    if image.ndim != 2 or image.shape[0] != image.shape[1] or image.shape[0] % 2 != 1:
        raise ValueError("PSF image must be square with odd side length")
    if not np.all(np.isfinite(image)):
        raise ValueError("PSF image contains non-finite pixels")
    if not np.isfinite(pixel_scale_arcsec) or pixel_scale_arcsec <= 0:
        raise ValueError("pixel scale must be finite and positive")
    if (edges.ndim != 1 or len(edges) < 2 or not np.all(np.isfinite(edges))
            or np.any(np.diff(edges) <= 0) or edges[0] < 0):
        raise ValueError("radial edges must be finite, nonnegative, and increasing")

    total = float(np.sum(image, dtype=np.float64))
    if not np.isfinite(total) or total <= 0:
        raise ValueError("PSF total flux must be finite and positive")
    side = image.shape[0]
    center = (side - 1) / 2.0
    yy, xx = np.indices(image.shape, dtype=np.float32)
    radius = np.hypot(xx - center, yy - center)
    peak_y, peak_x = np.unravel_index(np.argmax(image), image.shape)
    peak_offset = float(np.hypot(peak_x - center, peak_y - center))
    annulus_index = np.searchsorted(edges, radius, side="right") - 1
    in_profile = (annulus_index >= 0) & (annulus_index < len(edges) - 1)
    flat_index = annulus_index[in_profile]
    flat_flux = image[in_profile]
    sums = np.bincount(flat_index, weights=flat_flux, minlength=len(edges) - 1)
    counts = np.bincount(flat_index, minlength=len(edges) - 1)
    profile = []
    cumulative = np.cumsum(sums, dtype=np.float64) / total
    for index, (lo, hi) in enumerate(zip(edges[:-1], edges[1:])):
        values = image[in_profile & (annulus_index == index)]
        profile.append({
            "r_inner_pix": float(lo), "r_outer_pix": float(hi),
            "r_mid_pix": float(np.sqrt(lo * hi)),
            "r_inner_arcsec": float(lo * pixel_scale_arcsec),
            "r_outer_arcsec": float(hi * pixel_scale_arcsec),
            "n_pixels": int(counts[index]),
            "mean_flux_per_pixel": float(sums[index] / counts[index]) if counts[index] else None,
            "median_flux_per_pixel": float(np.median(values)) if len(values) else None,
            "annular_flux_fraction": float(sums[index] / total),
            "encircled_flux_fraction_at_outer_edge": float(cumulative[index]),
        })

    return profile, {
        "shape": [int(side), int(side)],
        "pixel_scale_arcsec": float(pixel_scale_arcsec),
        "total_flux_raw": total,
        "minimum_pixel": float(np.min(image)),
        "maximum_pixel": float(np.max(image)),
        "peak_xy_pixel": [int(peak_x), int(peak_y)],
        "geometric_center_xy_pixel": [float(center), float(center)],
        "peak_offset_from_center_pixel": peak_offset,
        "profile_outer_radius_pix": float(edges[-1]),
        "profile_outer_radius_arcmin": float(edges[-1] * pixel_scale_arcsec / 60.0),
        "profile_covered_flux_fraction": float(np.sum(sums) / total),
        "science_measurement": False,
        "psf_approved": False,
    }


def parse_input(value: str):
    if "=" not in value:
        raise argparse.ArgumentTypeError("use LABEL=/absolute/path/to/model.fits")
    label, raw_path = value.split("=", 1)
    if not re.fullmatch(r"[A-Za-z0-9_-]+", label):
        raise argparse.ArgumentTypeError("label may contain letters, digits, _ and - only")
    return label, Path(raw_path).expanduser().resolve()


def physical_aperture_diagnostics(data: np.ndarray, redshifts=(0.30, 0.55),
                                 pixel_scale_arcsec=DEFAULT_PIXEL_SCALE_ARCSEC,
                                 h0=70.0, omega_m=0.3):
    image = np.asarray(data, dtype=np.float64)
    side = image.shape[0]
    center = (side - 1) / 2.0
    yy, xx = np.indices(image.shape, dtype=np.float32)
    radius_pix = np.hypot(xx - center, yy - center)
    total = float(np.sum(image, dtype=np.float64))
    cosmology = FlatLambdaCDM(H0=h0, Om0=omega_m)
    results = []
    for redshift in redshifts:
        kpc_per_arcsec = float(cosmology.kpc_proper_per_arcmin(redshift).value / 60.0)
        kpc_per_pixel = kpc_per_arcsec * pixel_scale_arcsec
        radii = {}
        for radius_kpc in (50, 75, 100, 200, 250, 300):
            radius_pixels = radius_kpc / kpc_per_pixel
            inside = radius_pix <= radius_pixels
            radii[str(radius_kpc)] = {
                "radius_arcsec": float(radius_kpc / kpc_per_arcsec),
                "radius_pixels": float(radius_pixels),
                "encircled_flux_fraction": float(np.sum(image[inside], dtype=np.float64) / total),
            }
        annuli = {}
        sectors = np.floor((np.arctan2(yy - center, xx - center) + np.pi)
                           * (8.0 / (2.0 * np.pi))).astype(np.int8)
        sectors = np.clip(sectors, 0, 7)
        for name, (inner_kpc, outer_kpc) in {
            "satellite_50_100_kpc": (50, 100),
            "icl_200_300_kpc": (200, 300),
        }.items():
            lo, hi = inner_kpc / kpc_per_pixel, outer_kpc / kpc_per_pixel
            annulus = (radius_pix > lo) & (radius_pix <= hi)
            flux = float(np.sum(image[annulus], dtype=np.float64))
            sector_flux = np.array([
                np.sum(image[annulus & (sectors == index)], dtype=np.float64)
                for index in range(8)
            ])
            sector_shares = sector_flux / flux if flux != 0 else np.full(8, np.nan)
            annuli[name] = {
                "inner_radius_kpc": inner_kpc, "outer_radius_kpc": outer_kpc,
                "annular_flux_fraction_total": flux / total,
                "eight_sector_flux_shares": [float(v) for v in sector_shares],
                "sector_share_rms_about_equal": float(np.sqrt(np.mean((sector_shares - 0.125) ** 2))),
            }
        results.append({"redshift": float(redshift),
                        "proper_kpc_per_arcsec": kpc_per_arcsec,
                        "proper_kpc_per_pixel": kpc_per_pixel,
                        "apertures": radii, "annuli": annuli})
    return results


def audit_fits(label: str, path: Path, pixel_scale_arcsec: float,
               redshifts=(0.30, 0.55), h0=70.0, omega_m=0.3):
    path = Path(path).resolve()
    if not path.is_file():
        raise FileNotFoundError(path)
    with fits.open(path, memmap=True, checksum=True) as hdus:
        if len(hdus) < 2 or hdus[1].data is None:
            raise ValueError(f"{path}: expected image data in FITS extension 1")
        data = np.asarray(hdus[1].data)
        profile, metrics = analyze_image(data, pixel_scale_arcsec)
        physical = physical_aperture_diagnostics(
            data, redshifts, pixel_scale_arcsec, h0=h0, omega_m=omega_m)
        primary_header = hdus[0].header.tostring(sep="\n", endcard=False, padding=False)
        image_header = hdus[1].header.tostring(sep="\n", endcard=False, padding=False)
        hdu_shapes = [list(hdu.data.shape) if hdu.data is not None else None for hdu in hdus]
    return {
        "label": label,
        "path": str(path),
        "bytes": path.stat().st_size,
        "sha256": sha256(path),
        "hdu_shapes": hdu_shapes,
        "primary_header_cards": primary_header,
        "image_header_cards": image_header,
        "physical_aperture_diagnostics": physical,
        "cosmology": {"model": "FlatLambdaCDM", "H0_km_s_Mpc": h0, "Omega_m": omega_m},
        **metrics,
    }, profile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fits", action="append", type=parse_input, required=True,
                        metavar="LABEL=PATH", help="repeat for each locally available model")
    parser.add_argument("--pixel-scale-arcsec", type=float, default=DEFAULT_PIXEL_SCALE_ARCSEC)
    parser.add_argument("--redshift", type=float, action="append", default=None,
                        help="cluster redshift for physical aperture diagnostics (repeatable; default: 0.30, 0.55)")
    parser.add_argument("--h0", type=float, default=70.0)
    parser.add_argument("--omega-m", type=float, default=0.3)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    labels = [label for label, _ in args.fits]
    if len(set(labels)) != len(labels):
        parser.error("input labels must be unique")
    redshifts = tuple(args.redshift if args.redshift is not None else (0.30, 0.55))
    if any(not np.isfinite(z) or z <= 0 for z in redshifts):
        parser.error("redshifts must be finite and positive")
    output = args.output.resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite {output}")
    output.parent.mkdir(parents=True, exist_ok=True)
    artifacts = []
    with tempfile.TemporaryDirectory(prefix=output.name + ".tmp-", dir=output.parent) as tmp:
        stage = Path(tmp)
        for label, path in args.fits:
            metadata, profile = audit_fits(label, path, args.pixel_scale_arcsec,
                                           redshifts, args.h0, args.omega_m)
            profile_path = stage / f"{label}_radial_profile.csv"
            with profile_path.open("x", newline="", encoding="utf-8") as stream:
                writer = csv.DictWriter(stream, fieldnames=list(profile[0]))
                writer.writeheader()
                writer.writerows(profile)
            metadata["profile_csv"] = str(output / profile_path.name)
            metadata["profile_csv_sha256"] = sha256(profile_path)
            artifacts.append(metadata)
        summary = {
            "schema": "astrafek-hsc-public-extended-psf-reference-audit-v1",
            "status": "complete_reference_characterization_only",
            "inputs": artifacts,
            "psf_approved": False,
            "science_measurement": False,
            "limitations": [
                "A survey-level reference PSF is not a local per-frame PSF calibration.",
                "This audit does not correct sky transfer, masks, or source scattering in science images.",
                "No ICL measurement or PSF approval is produced.",
            ],
        }
        summary_path = stage / "summary.json"
        summary_path.write_text(json.dumps(summary, indent=2, allow_nan=False) + "\n", encoding="utf-8")
        if output.exists():
            raise FileExistsError(f"refusing to overwrite {output}")
        os.rename(stage, output)
    summary_path = output / "summary.json"
    print(json.dumps({"status": summary["status"], "inputs": len(artifacts),
                      "output": str(output), "summary_sha256": sha256(summary_path)}, indent=2))


if __name__ == "__main__":
    main()
