#!/usr/bin/env python3
"""Independently validate the approved public HSC extended-PSF audit."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.cosmology import FlatLambdaCDM


EXPECTED = {
    "HSC_PSF_i.fits": {
        "bytes": 9809280,
        "md5": "5bf3317bcc41a83c81c25acab2faab09",
        "sha256": "8897af89767e5fbb8d5f9c70965126bc89409f6c65ecd8d117c7995f34be2d85",
        "label": "symmetric",
    },
    "HSC_PSF_i_asymmetric.fits": {
        "bytes": 9809280,
        "md5": "665e7b459cda203d6610133c3626194d",
        "sha256": "b2ab473c6af4da11c1df74850bdd0a8cc0c1da4d109c118d1f4033b4235027fd",
        "label": "asymmetric",
    },
}


def digest(path: Path, algorithm="sha256") -> str:
    h = hashlib.new(algorithm)
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def _close(observed, expected, *, rtol=2e-9, atol=1e-12):
    if not np.isclose(float(observed), float(expected), rtol=rtol, atol=atol):
        raise ValueError(f"independent value mismatch: {observed} != {expected}")


def validate(data_dir: Path, audit_dir: Path) -> dict:
    data_dir, audit_dir = Path(data_dir).resolve(), Path(audit_dir).resolve()
    manifest_path = data_dir / "manifest.json"
    summary_path = audit_dir / "summary.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    if (manifest.get("doi") != "10.5281/zenodo.14915578"
            or manifest.get("record_version") != "v4.0.0"
            or manifest.get("science_measurement") is not False):
        raise ValueError("acquisition manifest provenance/status mismatch")
    if (summary.get("schema") != "astrafek-hsc-public-extended-psf-reference-audit-v1"
            or summary.get("status") != "complete_reference_characterization_only"
            or summary.get("psf_approved") is not False
            or summary.get("science_measurement") is not False):
        raise ValueError("audit schema/status mismatch")
    listed = {row["path"]: row for row in manifest.get("files", [])}
    if set(listed) != set(EXPECTED):
        raise ValueError("manifest does not contain exactly the two approved FITS files")
    artifact_rows = {row["label"]: row for row in summary.get("inputs", [])}
    if set(artifact_rows) != {"symmetric", "asymmetric"}:
        raise ValueError("audit summary does not contain both PSF variants")

    cosmo = FlatLambdaCDM(H0=70.0, Om0=0.3)
    z_values = (0.30, 0.55)
    physical_radii_kpc = (50, 75, 100, 200, 250, 300)
    checks = {}
    for name, expected in EXPECTED.items():
        path = data_dir / name
        meta = listed[name]
        if (not path.is_file() or path.stat().st_size != expected["bytes"]
                or digest(path, "md5") != expected["md5"]
                or digest(path) != expected["sha256"]):
            raise ValueError(f"source file size/checksum mismatch: {name}")
        for key, value in (("bytes", expected["bytes"]),
                           ("zenodo_md5", expected["md5"]),
                           ("sha256", expected["sha256"])):
            if meta.get(key) != value:
                raise ValueError(f"manifest {key} mismatch: {name}")
        record = artifact_rows[expected["label"]]
        if record.get("sha256") != expected["sha256"] or record.get("bytes") != expected["bytes"]:
            raise ValueError(f"audit source provenance mismatch: {name}")
        profile_path = audit_dir / f"{expected['label']}_radial_profile.csv"
        if digest(profile_path) != record.get("profile_csv_sha256"):
            raise ValueError(f"radial profile checksum mismatch: {profile_path.name}")
        with profile_path.open(newline="", encoding="utf-8") as stream:
            profile = list(csv.DictReader(stream))
        if len(profile) != 49 or float(profile[-1]["r_outer_pix"]) != 2000.0:
            raise ValueError(f"radial profile coverage mismatch: {profile_path.name}")
        _close(sum(float(row["annular_flux_fraction"]) for row in profile),
               record["profile_covered_flux_fraction"])

        with fits.open(path, memmap=True, checksum=True) as hdus:
            if len(hdus) < 2 or hdus[1].data is None or hdus[1].data.shape != (4001, 4001):
                raise ValueError(f"unexpected FITS image layout: {name}")
            image = np.asarray(hdus[1].data, dtype=np.float64)
        if not np.all(np.isfinite(image)):
            raise ValueError(f"non-finite PSF pixel: {name}")
        total = float(np.sum(image, dtype=np.float64))
        _close(total, record["total_flux_raw"], rtol=1e-11)
        if np.unravel_index(np.argmax(image), image.shape) != (2000, 2000):
            raise ValueError(f"PSF peak is not centered: {name}")
        y, x = np.indices(image.shape, dtype=np.float32)
        radius_pix = np.hypot(x - 2000.0, y - 2000.0)
        per_z = []
        for redshift in z_values:
            kpc_per_arcsec = float(cosmo.kpc_proper_per_arcmin(redshift).value / 60.0)
            kpc_per_pix = kpc_per_arcsec * 0.168
            apertures = {}
            for radius_kpc in physical_radii_kpc:
                radius_pixels = radius_kpc / kpc_per_pix
                fraction = float(np.sum(image[radius_pix <= radius_pixels], dtype=np.float64) / total)
                saved = record["physical_aperture_diagnostics"][z_values.index(redshift)]
                _close(fraction, saved["apertures"][str(radius_kpc)]["encircled_flux_fraction"],
                       rtol=1e-10)
                apertures[str(radius_kpc)] = fraction
            lo, hi = 200.0 / kpc_per_pix, 300.0 / kpc_per_pix
            annular_fraction = float(np.sum(image[(radius_pix > lo) & (radius_pix <= hi)],
                                            dtype=np.float64) / total)
            saved_annulus = record["physical_aperture_diagnostics"][z_values.index(redshift)]["annuli"]["icl_200_300_kpc"]
            _close(annular_fraction, saved_annulus["annular_flux_fraction_total"], rtol=1e-10)
            per_z.append({"redshift": redshift, "encircled_fraction_by_kpc": apertures,
                          "icl_200_300_kpc_fraction": annular_fraction})
        checks[expected["label"]] = {
            "bytes": path.stat().st_size, "md5": expected["md5"],
            "sha256": expected["sha256"], "hdu1_shape": [4001, 4001],
            "raw_total_flux": total, "peak_pixel_yx": [2000, 2000],
            "independent_physical_recalculation": per_z,
        }
    return {
        "schema": "astrafek-hsc-public-extended-psf-reference-validation-v1",
        "status": "validated_reference_characterization_only",
        "science_measurement": False,
        "psf_approved": False,
        "source_manifest_sha256": digest(manifest_path),
        "source_summary_sha256": digest(summary_path),
        "variants": checks,
        "limitations": [
            "Validation does not establish local-frame PSF transfer or sky correction.",
            "Point-source annular leakage is not the ICL contamination fraction.",
            "No science measurement or PSF approval is produced.",
        ],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-dir", type=Path, required=True)
    parser.add_argument("--audit-dir", type=Path, required=True)
    args = parser.parse_args()
    output = args.audit_dir.resolve() / "validation.json"
    if output.exists():
        raise FileExistsError(f"refusing to overwrite {output}")
    result = validate(args.data_dir, args.audit_dir)
    output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({"status": result["status"], "output": str(output),
                      "sha256": digest(output)}, indent=2))


if __name__ == "__main__":
    main()
