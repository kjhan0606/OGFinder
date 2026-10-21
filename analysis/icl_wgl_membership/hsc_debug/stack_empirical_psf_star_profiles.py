#!/usr/bin/env python3
"""Frame-weighted stack of existing Gaia-supported HSC star diagnostics.

This is a diagnostic of the available stellar wings, not a certified PSF model
or a science measurement. It reads only stamps listed in one immutable ledger.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
from astropy.io import fits

from fetch_ancillary_catalogues import sha

RADIAL_EDGES = np.asarray((0, 2, 4, 6, 8, 12, 16, 24, 32, 48, 64, 96), dtype=float)
MASKS = ((2.5, "MASK_K2P5"), (3.5, "MASK_K3P5"))
SKY_ANNULI = ((120, 150), (160, 190))


def _table(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def measure_stamp(stamp: Path, frame_id: str, candidate_id: str) -> list[dict]:
    """Return one row per mask, sky annulus, and radial bin."""
    with fits.open(stamp, memmap=False) as hdus:
        header = hdus[0].header
        if header.get("SCIVAL") is not False or header.get("PSFOK") is not False:
            raise ValueError(f"stamp has unexpected science/PSF state: {stamp}")
        image = np.asarray(hdus["SCIENCE"].data, dtype=float)
        base = np.asarray(hdus["BASE_MASK"].data, dtype=bool)
        center_x, center_y = float(header["XCENTER"]), float(header["YCENTER"])
        pixel_scale = float(header["PIXSCALE"])
        if image.shape != (401, 401) or base.shape != image.shape:
            raise ValueError(f"expected 401x401 stamp and base mask: {stamp}")
        variance = np.asarray(hdus["VARIANCE"].data, dtype=float)
        if variance.shape != image.shape:
            raise ValueError(f"variance shape mismatch: {stamp}")
        yy, xx = np.indices(image.shape, dtype=float)
        radius = np.hypot(xx - center_x, yy - center_y)
        sector = np.floor((np.arctan2(yy - center_y, xx - center_x) + np.pi) * (4 / np.pi)).astype(int)
        sector = np.clip(sector, 0, 7)

        output = []
        for factor, mask_name in MASKS:
            if mask_name not in hdus:
                raise ValueError(f"missing {mask_name} in {stamp}")
            source_mask = np.asarray(hdus[mask_name].data, dtype=bool)
            good = (~base & ~source_mask & np.isfinite(image) & np.isfinite(variance) & (variance > 0))
            core = (radius < 12) & good
            core_area = np.count_nonzero(radius < 12)
            core_fraction = float(np.count_nonzero(core) / core_area)
            for sky_lo, sky_hi in SKY_ANNULI:
                sky_region = (radius >= sky_lo) & (radius < sky_hi)
                sky_good = sky_region & good
                sector_means = []
                sector_counts = []
                sector_fractions = []
                for s in range(8):
                    region = sky_region & (sector == s)
                    valid = region & good
                    count = int(valid.sum())
                    sector_counts.append(count)
                    sector_fractions.append(count / max(1, int(region.sum())))
                    if count:
                        sector_means.append(float(np.mean(image[valid])))
                sky_ok = (sky_good.sum() >= 300 and sky_good.sum() / max(1, sky_region.sum()) >= 0.3
                          and min(sector_counts) >= 20 and min(sector_fractions) >= 0.25)
                status = "ok"
                if not sky_ok:
                    status = "blocked_background_coverage"
                    background = core_flux = None
                else:
                    # Same equal-azimuth weighting used in the original diagnostic.
                    background = float(np.mean(sector_means))
                    core_flux = float(np.sum(image[core] - background))
                    core_noise = float(np.sqrt(np.sum(variance[core])))
                    if core_fraction < 0.99:
                        status = "blocked_core_mask"
                    elif not np.isfinite(core_flux) or core_flux <= 0 or core_noise <= 0 or core_flux / core_noise < 30:
                        status = "blocked_core_flux"
                for lo, hi in zip(RADIAL_EDGES[:-1], RADIAL_EDGES[1:]):
                    ring = (radius >= lo) & (radius < hi)
                    valid = ring & good
                    area = int(ring.sum())
                    count = int(valid.sum())
                    covered_sectors = 0
                    for s in range(8):
                        sector_area = int(np.count_nonzero(ring & (sector == s)))
                        sector_valid = int(np.count_nonzero(valid & (sector == s)))
                        if sector_valid >= max(1, 0.25 * sector_area):
                            covered_sectors += 1
                    radial_ok = count / max(1, area) >= 0.5 and covered_sectors == 8
                    value = None
                    if status == "ok" and radial_ok and count:
                        value = float((np.mean(image[valid]) - background) / core_flux)
                    row_status = status if status != "ok" else ("ok" if radial_ok else "blocked_radial_coverage")
                    output.append({
                        "candidate_id": str(candidate_id), "frame_id": str(frame_id),
                        "kron_factor": factor, "mask_name": mask_name,
                        "sky_r_min_pix": sky_lo, "sky_r_max_pix": sky_hi,
                        "r_min_pix": float(lo), "r_max_pix": float(hi),
                        "r_mid_arcsec": float((lo + hi) * 0.5 * pixel_scale),
                        "n_valid": count, "area_fraction": count / max(1, area),
                        "n_covered_sectors": covered_sectors, "core_flux": core_flux,
                        "background": background, "profile_per_core_per_pixel": value,
                        "status": row_status, "science_measurement": False, "psf_approved": False,
                    })
        return output


def frame_weighted_stack(rows: list[dict], *, n_bootstrap: int = 10000, seed: int = 20261020) -> list[dict]:
    """Aggregate stars within frame first, then give each frame equal weight."""
    keys = sorted({(float(r["kron_factor"]), int(r["sky_r_min_pix"]), int(r["sky_r_max_pix"]),
                    float(r["r_min_pix"]), float(r["r_max_pix"])) for r in rows})
    rng = np.random.default_rng(seed)
    stacked = []
    for factor, sky_lo, sky_hi, rlo, rhi in keys:
        frame_values: dict[str, list[tuple[float, float]]] = defaultdict(list)
        for row in rows:
            if (row["status"] == "ok" and float(row["kron_factor"]) == factor
                    and int(row["sky_r_min_pix"]) == sky_lo and int(row["sky_r_max_pix"]) == sky_hi
                    and float(row["r_min_pix"]) == rlo and float(row["r_max_pix"]) == rhi):
                frame_values[str(row["frame_id"])].append((float(row["profile_per_core_per_pixel"]), float(row["r_mid_arcsec"])))
        frames = sorted(frame_values)
        per_frame = np.asarray([np.mean([v for v, _ in frame_values[f]]) for f in frames], dtype=float)
        per_frame_radius = np.asarray([np.mean([x for _, x in frame_values[f]]) for f in frames], dtype=float)
        if len(per_frame):
            mean = float(np.mean(per_frame))
            median = float(np.median(per_frame))
            if len(per_frame) >= 2 and n_bootstrap > 0:
                draws = rng.integers(0, len(per_frame), size=(n_bootstrap, len(per_frame)))
                boot_mean = np.mean(per_frame[draws], axis=1)
                low, high = (float(v) for v in np.quantile(boot_mean, (0.16, 0.84)))
            else:
                low = high = None
            mean_radius = float(np.mean(per_frame_radius))
        else:
            mean = median = low = high = mean_radius = None
        stacked.append({
            "kron_factor": factor, "sky_r_min_pix": sky_lo, "sky_r_max_pix": sky_hi,
            "r_min_pix": rlo, "r_max_pix": rhi, "r_mid_arcsec": mean_radius,
            "n_stars": sum(len(v) for v in frame_values.values()), "n_frames": len(frames),
            "frame_equal_mean": mean, "frame_equal_median": median,
            "bootstrap_mean_q16": low, "bootstrap_mean_q84": high,
            "science_measurement": False, "psf_approved": False,
        })
    return stacked


def _write_tsv(path: Path, rows: list[dict]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    with path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def _plot(stack: list[dict], destination: Path) -> None:
    from matplotlib import pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharex=True)
    for ax, (factor, sky) in zip(axes.flat, ((2.5, (120, 150)), (2.5, (160, 190)),
                                             (3.5, (120, 150)), (3.5, (160, 190)))):
        selected = [r for r in stack if r["kron_factor"] == factor and r["sky_r_min_pix"] == sky[0]
                    and r["sky_r_max_pix"] == sky[1] and r["n_frames"] >= 2]
        selected.sort(key=lambda r: r["r_mid_arcsec"])
        if selected:
            x = np.asarray([r["r_mid_arcsec"] for r in selected])
            y = np.asarray([r["frame_equal_mean"] for r in selected])
            lo = np.asarray([r["bootstrap_mean_q16"] for r in selected])
            hi = np.asarray([r["bootstrap_mean_q84"] for r in selected])
            ax.plot(x, y, marker="o", ms=3, lw=1)
            if np.all(np.isfinite(lo)) and np.all(np.isfinite(hi)):
                ax.fill_between(x, lo, hi, alpha=0.2)
        ax.axhline(0, color="0.4", lw=0.8)
        ax.set_yscale("symlog", linthresh=1e-7)
        ax.set_title(f"Kron {factor:g}; sky {sky[0]}–{sky[1]} px")
        ax.set_ylabel("mean / core flux per pixel")
        ax.grid(alpha=0.18)
    for ax in axes[-1]:
        ax.set_xlabel("Radius [arcsec]")
    fig.suptitle("Frame-weighted empirical-star wing diagnostic (not PSF-approved)")
    fig.tight_layout()
    fig.savefig(destination, dpi=160)
    plt.close(fig)


def run(input_dir: Path, output_dir: Path, *, n_bootstrap: int = 10000, seed: int = 20261020) -> dict:
    input_dir, output_dir = input_dir.resolve(), output_dir.resolve()
    ledger_path = input_dir / "star_ledger.tsv"
    if not ledger_path.is_file():
        raise FileNotFoundError(ledger_path)
    if output_dir.exists():
        raise FileExistsError(f"refusing to overwrite {output_dir}")
    ledger = _table(ledger_path)
    if not ledger or not {"candidate_id", "cluster_id", "stamp", "stamp_sha256"} <= set(ledger[0]):
        raise ValueError("unexpected star ledger schema")
    stamped = [row for row in ledger if row.get("stamp")]
    if not stamped:
        raise ValueError("the star ledger contains no complete stamps")
    input_hashes = {str(ledger_path): sha(ledger_path)}
    profiles = []
    for row in stamped:
        stamp = Path(row["stamp"]).resolve()
        if not stamp.is_file() or sha(stamp) != row["stamp_sha256"]:
            raise ValueError(f"stamp checksum mismatch: {stamp}")
        input_hashes[str(stamp)] = row["stamp_sha256"]
        profiles.extend(measure_stamp(stamp, row["cluster_id"], row["candidate_id"]))
    stack = frame_weighted_stack(profiles, n_bootstrap=n_bootstrap, seed=seed)
    output_dir.mkdir(parents=True, exist_ok=False)
    _write_tsv(output_dir / "star_radial_profiles.tsv", profiles)
    _write_tsv(output_dir / "frame_weighted_stack.tsv", stack)
    _plot(stack, output_dir / "frame_weighted_psf_wings.png")
    summary = {
        "schema": "astrafek-empirical-psf-frame-weighted-stack-v1", "status": "complete_diagnostic_only",
        "n_ledger_rows": len(ledger), "n_stamps": len(stamped),
        "n_unique_frames": len({row["cluster_id"] for row in stamped}),
        "n_profile_rows": len(profiles), "n_stacked_radial_bins": len(stack),
        "bootstrap_unit": "HSC frame (cluster_id); stars in a frame are averaged before resampling",
        "bootstrap_draws": int(n_bootstrap), "random_seed": int(seed),
        "script_sha256": sha(Path(__file__)),
        "input_sha256": input_hashes,
        "science_measurement": False, "psf_approved": False,
        "limitations": [
            "Gaia proper-motion propagation to the effective HSC coadd epoch is unresolved.",
            f"Only {len(stamped)} complete stamps from {len({row['cluster_id'] for row in stamped})} frames in this ledger are used; they are not representative of all 800 cluster frames.",
            "The two measured sky annuli may contain PSF wings, ICL, undetected sources, or residual sky structure.",
            "The radial support/sector rules are geometry gates; bootstrap intervals describe frame resampling, not total uncertainty.",
            "No effective coadd PSF, extended-wing correction, or science surface-brightness correction is approved.",
        ],
        "stack_tsv_sha256": sha(output_dir / "frame_weighted_stack.tsv"),
        "profile_tsv_sha256": sha(output_dir / "star_radial_profiles.tsv"),
        "plot_sha256": sha(output_dir / "frame_weighted_psf_wings.png"),
    }
    (output_dir / "summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True, help="existing empirical_psf_... ledger directory")
    parser.add_argument("--output", type=Path, required=True, help="new output directory; must not exist")
    parser.add_argument("--bootstrap-draws", type=int, default=10000)
    parser.add_argument("--seed", type=int, default=20261020)
    args = parser.parse_args()
    if args.bootstrap_draws < 100:
        parser.error("--bootstrap-draws must be at least 100")
    print(json.dumps(run(args.input, args.output, n_bootstrap=args.bootstrap_draws, seed=args.seed),
                     ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
