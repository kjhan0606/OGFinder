#!/usr/bin/env python3
"""Compare the published HSC extended-PSF core with approved local PSF cores.

The public model is cropped and renormalized only to compare core *shape*.
The cropped fraction is recorded separately; this comparison is not a hybrid
kernel, wing correction, or science measurement.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from astropy.io import fits

from summarize_hsc_psf_core_variation import measure_core


EXPECTED_PUBLIC = {
    "HSC_PSF_i.fits": {
        "label": "symmetric",
        "bytes": 9809280,
        "sha256": "8897af89767e5fbb8d5f9c70965126bc89409f6c65ecd8d117c7995f34be2d85",
    },
    "HSC_PSF_i_asymmetric.fits": {
        "label": "asymmetric",
        "bytes": 9809280,
        "sha256": "b2ab473c6af4da11c1df74850bdd0a8cc0c1da4d109c118d1f4033b4235027fd",
    },
}
EXPECTED_LOCAL_SUMMARY_SHA256 = "2a1a0f63f7cbabadaf6f1a99f27d53ea9b87d8df2232305c73a89e8ce75807ad"
PIXEL_SCALE_ARCSEC = 0.168


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def measure_public_core(data: np.ndarray, side: int = 43) -> tuple[dict, float]:
    """Measure shape metrics on a centered odd crop, separately recording flux."""
    image = np.asarray(data, dtype=np.float64)
    if image.ndim != 2 or image.shape[0] != image.shape[1] or image.shape[0] % 2 != 1:
        raise ValueError("public PSF must be an odd square image")
    if not np.isfinite(image).all() or not np.isfinite(image.sum()) or image.sum() <= 0:
        raise ValueError("public PSF must have finite positive total flux")
    if side <= 0 or side % 2 != 1 or side > image.shape[0]:
        raise ValueError("core crop side must be positive, odd, and no larger than the image")
    center = image.shape[0] // 2
    half = side // 2
    crop = image[center - half:center + half + 1, center - half:center + half + 1].copy()
    crop_flux = float(crop.sum())
    total_flux = float(image.sum())
    if not np.isfinite(crop_flux) or crop_flux <= 0:
        raise ValueError("public PSF core crop has nonpositive flux")
    metrics = measure_core(crop / crop_flux, PIXEL_SCALE_ARCSEC)
    return metrics, crop_flux / total_flux


def distribution(records: list[dict], key: str) -> dict:
    values = np.asarray([float(row[key]) for row in records], dtype=float)
    if not len(values) or not np.isfinite(values).all():
        raise ValueError(f"no finite values for {key}")
    return {"n": int(len(values)), "median": float(np.median(values)),
            "q16": float(np.quantile(values, 0.16)), "q84": float(np.quantile(values, 0.84)),
            "q05": float(np.quantile(values, 0.05)), "q95": float(np.quantile(values, 0.95))}


def compare(public_files: dict[str, Path], local_summary_path: Path) -> dict:
    local_summary_path = Path(local_summary_path).resolve()
    if sha256(local_summary_path) != EXPECTED_LOCAL_SUMMARY_SHA256:
        raise ValueError("approved local PSF-core summary checksum mismatch")
    local_summary = json.loads(local_summary_path.read_text(encoding="utf-8"))
    if (local_summary.get("schema") != "astrafek-hsc-position-psf-core-variation-v1"
            or local_summary.get("status") != "complete_exploratory_psf_core_audit"
            or local_summary.get("n_total_models_including_test") != 745
            or local_summary.get("science_measurement") is not False):
        raise ValueError("local position-dependent PSF-core audit is incomplete or promoted")
    local_records = [row for row in local_summary["records"]
                     if row.get("source") == "approved_744_position_batch"]
    if len(local_records) != 744:
        raise ValueError(f"expected the approved 744 local models, found {len(local_records)}")
    local_distributions = {
        key: distribution(local_records, key)
        for key in ("moment_fwhm_major_arcsec", "moment_fwhm_minor_arcsec",
                    "encircled_energy_r50_arcsec", "encircled_energy_r80_arcsec",
                    "encircled_energy_r90_arcsec")
    }

    public_metrics = {}
    for name, expected in EXPECTED_PUBLIC.items():
        path = Path(public_files[expected["label"]]).resolve()
        if not path.is_file() or path.stat().st_size != expected["bytes"]:
            raise ValueError(f"public PSF file missing or has wrong size: {path}")
        file_hash = sha256(path)
        if file_hash != expected["sha256"]:
            raise ValueError(f"public PSF checksum mismatch: {path}")
        with fits.open(path, memmap=True, checksum=True) as hdus:
            if len(hdus) < 2 or hdus[1].data is None or hdus[1].data.shape != (4001, 4001):
                raise ValueError(f"unexpected public PSF FITS shape: {path}")
            image = np.asarray(hdus[1].data, dtype=np.float64)
        if not np.isfinite(image).all():
            raise ValueError(f"non-finite public PSF pixels: {path}")
        core, crop_fraction = measure_public_core(image)
        public_metrics[expected["label"]] = {
            "path": str(path), "bytes": path.stat().st_size, "sha256": file_hash,
            "full_image_shape": list(image.shape), "core_crop_side": 43,
            "signed_flux_fraction_inside_core_crop": crop_fraction,
            "core_shape_metrics": core,
            "core_to_local_median_ratios": {
                key: core[key] / local_distributions[key]["median"]
                for key in local_distributions
            },
        }
    return {
        "schema": "astrafek-hsc-public-extended-psf-local-core-comparison-v1",
        "status": "complete_core_shape_comparison_only",
        "public_psf_files": public_metrics,
        "local_picker_core_summary": {
            "path": str(local_summary_path), "sha256": EXPECTED_LOCAL_SUMMARY_SHA256,
            "n_position_models": len(local_records), "distributions": local_distributions,
        },
        "comparison_definition": (
            "The central 43x43-pixel public model is renormalized to unit signed sum for core shape metrics; "
            "the pre-renormalization crop fraction is retained. Local PSF Picker models are not modified."
        ),
        "psf_approved": False,
        "science_measurement": False,
        "limitations": [
            "Core-shape mismatch does not test outer-wing transfer or target-image contamination.",
            "The survey-level public model is not a local per-frame PSF.",
            "No hybrid PSF, correction, or ICL measurement is produced.",
        ],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--symmetric", type=Path, required=True)
    parser.add_argument("--asymmetric", type=Path, required=True)
    parser.add_argument("--local-core-summary", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    output = args.output.resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite {output}")
    result = compare({"symmetric": args.symmetric, "asymmetric": args.asymmetric},
                     args.local_core_summary)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({"status": result["status"], "output": str(output),
                      "sha256": sha256(output)}, indent=2))


if __name__ == "__main__":
    main()
