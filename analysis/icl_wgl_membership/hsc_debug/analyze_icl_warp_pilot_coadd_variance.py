#!/usr/bin/env python3
"""Compare one approved HSC warp pilot's inverse-variance map with its coadd.

This is a local consistency diagnostic for one 2-arcmin cutout. It does not
infer source-Poisson noise, coadd weights for the 800-cluster sample, or a
science-measurement uncertainty. No network access is performed.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import tarfile
from pathlib import Path

import numpy as np
from astropy import units as u
from astropy.cosmology import FlatLambdaCDM
from astropy.io import fits
from astropy.wcs import WCS


DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
PLAN_DIR = DATA / "icl_warp_pilot_plan_2026-11-06_02"
DOWNLOAD_DIR = DATA / "icl_warp_pilot_download_2026-11-06_01"
COADD_PATH = DATA / "download_01/images/cluster_6966.fits"
SAMPLE_MANIFEST = DATA / "sample_manifest.json"
PLAN_SHA256 = "a436df245a33bda09cd5599e564b2b9302e5d0d4c006b59fcf6fab1c2c33a65e"
ARCHIVE_SHA256 = "cdc9b642c5c5d0361b58a0bf7c0e733a60cb78caad5ed496ae00cac89b23953c"
DOWNLOAD_SUMMARY_SHA256 = "858b26f411cda8131662fceb4b6eb62e913ff53eb19f26831f11439a67c1865c"
COADD_SHA256 = "66aa1e34fcc07f0335364bbc63fbbee192bdc33256ed109ad9023a77e10c83b9"
SAMPLE_MANIFEST_SHA256 = "95388814062e95f6e8429e2543a5432433515c9d9607680d2e873488d2758401"
DEFAULT_OUTPUT = DATA / "icl_warp_pilot_coadd_variance_2026-11-06_02"
NO_DATA_NAME = "MP_NO_DATA"
MASK_EXCLUDE_NAMES = (
    "MP_BAD", "MP_SAT", "MP_INTRP", "MP_EDGE", "MP_CR", "MP_SUSPECT",
    "MP_NO_DATA", "MP_BRIGHT_OBJECT", "MP_CROSSTALK", "MP_UNMASKEDNAN",
    "MP_CLIPPED", "MP_REJECTED", "MP_SENSOR_EDGE", "MP_INEXACT_PSF",
)
RADIUS_BINS_KPC = ((0, 50), (50, 100), (100, 200), (200, 300), (300, 400))


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def mask_bits(header) -> dict[str, int]:
    result = {}
    for key, value in header.items():
        if key.startswith("MP_"):
            result[key] = int(value)
    return result


def bitmask(names: tuple[str, ...], mapping: dict[str, int]) -> int:
    missing = sorted(set(names) - set(mapping))
    if missing:
        raise ValueError(f"mask plane is missing required bits: {missing}")
    return sum(1 << mapping[name] for name in names)


def quantile_summary(values: np.ndarray) -> dict:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if not values.size:
        return {"n_pixels": 0, "q16_median_q84": None, "median_relative_offset_percent": None,
                "fraction_within_1pct_of_unity": None, "fraction_within_2pct_of_unity": None}
    q16, median, q84 = np.quantile(values, [0.16, 0.5, 0.84])
    return {
        "n_pixels": int(values.size),
        "q16_median_q84": [float(q16), float(median), float(q84)],
        "median_relative_offset_percent": float(100 * (median - 1)),
        "fraction_within_1pct_of_unity": float(np.mean(np.abs(values - 1) <= 0.01)),
        "fraction_within_2pct_of_unity": float(np.mean(np.abs(values - 1) <= 0.02)),
    }


def load_pinned_inputs() -> tuple[dict, dict, dict, dict]:
    plan_path = PLAN_DIR / "plan.json"
    download_summary_path = DOWNLOAD_DIR / "summary.json"
    archive_path = DOWNLOAD_DIR / "warp_pilot.tar"
    for path in (plan_path, download_summary_path, archive_path, COADD_PATH, SAMPLE_MANIFEST):
        if not path.is_file():
            raise FileNotFoundError(path)
    hashes = {
        "plan_json": sha256(plan_path),
        "download_summary": sha256(download_summary_path),
        "warp_archive": sha256(archive_path),
        "coadd_fits": sha256(COADD_PATH),
        "sample_manifest": sha256(SAMPLE_MANIFEST),
    }
    expected = {
        "plan_json": PLAN_SHA256,
        "download_summary": DOWNLOAD_SUMMARY_SHA256,
        "warp_archive": ARCHIVE_SHA256,
        "coadd_fits": COADD_SHA256,
        "sample_manifest": SAMPLE_MANIFEST_SHA256,
    }
    if hashes != expected:
        raise ValueError(f"one or more pinned inputs changed: {hashes}")
    plan = json.loads(plan_path.read_text(encoding="utf-8"))
    download = json.loads(download_summary_path.read_text(encoding="utf-8"))
    sample = json.loads(SAMPLE_MANIFEST.read_text(encoding="utf-8"))
    if (plan.get("status") != "prepared_not_submitted"
            or plan.get("query_submitted") is not False
            or plan.get("science_measurement") is not False
            or plan.get("request", {}).get("center_id") != "target:6966"
            or plan.get("request", {}).get("type") != "warp"
            or plan.get("request", {}).get("layers") != {"image": True, "mask": True, "variance": True}
            or download.get("plan_sha256") != PLAN_SHA256
            or download.get("query_submitted") is not True
            or download.get("science_measurement") is not False
            or download.get("archive_sha256") != ARCHIVE_SHA256
            or not download.get("archive_audit", {}).get("visit_set_exact_match")):
        raise ValueError("warp-plan or downloaded-audit provenance/guard state is inconsistent")
    cluster = next((row for row in sample.get("clusters", []) if str(row.get("id")) == "6966"), None)
    if cluster is None or cluster.get("z_spec_bcg") is None:
        raise ValueError("target 6966 is not represented in the pinned parent sample manifest")
    return plan, download, sample, {"hashes": hashes, "cluster": cluster}


def run(output_dir: Path) -> dict:
    output_dir = Path(output_dir).resolve()
    if output_dir.exists():
        raise FileExistsError(f"refusing to overwrite {output_dir}")
    plan, download, sample, provenance = load_pinned_inputs()
    request = plan["request"]
    expected_visits = {int(value) for value in request["distinct_visits"]}
    archive_path = DOWNLOAD_DIR / "warp_pilot.tar"
    archive_records = []
    with fits.open(COADD_PATH, memmap=True) as coadd:
        if len(coadd) < 4 or str(coadd[0].header.get("FILTER", "")).lower() != "i2":
            raise ValueError("pinned target coadd is not the expected HSC-I2 image/mask/variance FITS")
        image = np.asarray(coadd[1].data)
        coadd_mask = np.asarray(coadd[2].data)
        coadd_variance = np.asarray(coadd[3].data)
        if image.ndim != 2 or image.shape != coadd_mask.shape or image.shape != coadd_variance.shape:
            raise ValueError("coadd image/mask/variance shapes differ")
        coadd_wcs = WCS(coadd[1].header).celestial
        center_x_raw, center_y_raw = coadd_wcs.world_to_pixel_values(
            float(request["ra_deg"]), float(request["dec_deg"]))
        center_x, center_y = float(np.asarray(center_x_raw)), float(np.asarray(center_y_raw))
        coadd_bits = mask_bits(coadd[2].header)
        if NO_DATA_NAME not in coadd_bits:
            raise ValueError("coadd mask does not define MP_NO_DATA")

        with tarfile.open(archive_path, mode="r:*") as archive:
            members = [member for member in archive.getmembers() if member.isfile()
                       and member.name.endswith(".fits")]
            if len(members) != len(expected_visits):
                raise ValueError(f"expected {len(expected_visits)} warp FITS members; got {len(members)}")
            visits_seen = set()
            inverse_variance_layers, valid_layers, clean_layers = [], [], []
            map_offsets = []
            mask_bit_maps = []
            strict_mask_names_used = set(MASK_EXCLUDE_NAMES)
            warp_footprint = np.zeros(image.shape, dtype=bool)
            for member in members:
                try:
                    visit = int(Path(member.name).stem.removeprefix("warp-"))
                except ValueError as exc:
                    raise ValueError(f"unrecognized visit file name: {member.name}") from exc
                if visit not in expected_visits or visit in visits_seen:
                    raise ValueError(f"unexpected or duplicate visit in approved archive: {visit}")
                visits_seen.add(visit)
                stream = archive.extractfile(member)
                if stream is None:
                    raise ValueError(f"cannot read warp member {member.name}")
                with fits.open(stream, memmap=False) as warp:
                    if len(warp) < 4:
                        raise ValueError(f"visit {visit} has no image/mask/variance planes")
                    w_image = np.asarray(warp[1].data)
                    w_mask = np.asarray(warp[2].data)
                    w_var = np.asarray(warp[3].data)
                    if w_image.ndim != 2 or w_image.shape != w_mask.shape or w_image.shape != w_var.shape:
                        raise ValueError(f"visit {visit} plane dimensions do not agree")
                    warp_wcs = WCS(warp[1].header).celestial
                    sample_x = np.asarray([0, w_image.shape[1] - 1, 0, w_image.shape[1] - 1,
                                           (w_image.shape[1] - 1) / 2.0])
                    sample_y = np.asarray([0, 0, w_image.shape[0] - 1, w_image.shape[0] - 1,
                                           (w_image.shape[0] - 1) / 2.0])
                    world = warp_wcs.pixel_to_world_values(sample_x, sample_y)
                    coadd_x, coadd_y = coadd_wcs.world_to_pixel_values(*world)
                    offsets = np.column_stack([coadd_x - sample_x, coadd_y - sample_y])
                    integer_offset = np.rint(np.median(offsets, axis=0)).astype(int)
                    residual = float(np.max(np.abs(offsets - integer_offset)))
                    if residual > 1e-3:
                        raise ValueError(f"visit {visit} warp is not on the coadd pixel grid ({residual} pixel)")
                    map_offsets.append(integer_offset.tolist())
                    y0, x0 = map(int, integer_offset[::-1])
                    y1, x1 = y0 + w_image.shape[0], x0 + w_image.shape[1]
                    if y0 < 0 or x0 < 0 or y1 > image.shape[0] or x1 > image.shape[1]:
                        raise ValueError(f"visit {visit} footprint lies outside the downloaded coadd")
                    w_mask_bits = mask_bits(warp[2].header)
                    if w_mask_bits.get(NO_DATA_NAME) != coadd_bits.get(NO_DATA_NAME):
                        raise ValueError(f"visit {visit} MP_NO_DATA bit does not match the coadd")
                    no_data_bit = 1 << w_mask_bits[NO_DATA_NAME]
                    visit_exclusions = tuple(name for name in MASK_EXCLUDE_NAMES if name in w_mask_bits)
                    if not visit_exclusions:
                        raise ValueError(f"visit {visit} has no supported strict-mask planes")
                    excluded = bitmask(visit_exclusions, w_mask_bits)
                    strict_mask_names_used.intersection_update(visit_exclusions)
                    valid = (np.isfinite(w_image) & np.isfinite(w_var) & (w_var > 0)
                             & ((w_mask & no_data_bit) == 0))
                    clean = valid & ((w_mask & excluded) == 0)
                    full_inverse_variance = np.zeros(image.shape, dtype=np.float64)
                    full_valid = np.zeros(image.shape, dtype=bool)
                    full_clean = np.zeros(image.shape, dtype=bool)
                    inv_local = np.zeros(w_var.shape, dtype=np.float64)
                    inv_local[valid] = 1.0 / w_var[valid]
                    full_inverse_variance[y0:y1, x0:x1] = inv_local
                    full_valid[y0:y1, x0:x1] = valid
                    full_clean[y0:y1, x0:x1] = clean
                    warp_footprint[y0:y1, x0:x1] = True
                    inverse_variance_layers.append(full_inverse_variance)
                    valid_layers.append(full_valid)
                    clean_layers.append(full_clean)
                    archive_records.append({
                        "visit": visit,
                        "member": member.name,
                        "bytes": int(member.size),
                        "valid_image_variance_fraction": float(valid.mean()),
                        "clean_pixel_fraction": float(clean.mean()),
                        "median_positive_variance": float(np.median(w_var[valid])) if valid.any() else None,
                        "median_nonzero_mask_fraction": float(np.mean(w_mask != 0)),
                        "strict_mask_exclusions_used": list(visit_exclusions),
                        "strict_mask_planes_absent": sorted(set(MASK_EXCLUDE_NAMES) - set(w_mask_bits)),
                        "coadd_integer_pixel_offset_xy": integer_offset.tolist(),
                        "wcs_grid_max_residual_pixel": residual,
                    })
            if visits_seen != expected_visits:
                raise ValueError("warp visit set differs from the approved visit inventory")

    n_valid = np.sum(valid_layers, axis=0, dtype=np.int16)
    n_clean = np.sum(clean_layers, axis=0, dtype=np.int16)
    inverse_variance_sum = np.sum(inverse_variance_layers, axis=0, dtype=np.float64)
    predicted_variance = np.divide(
        1.0, inverse_variance_sum,
        out=np.full(image.shape, np.nan, dtype=np.float64),
        where=inverse_variance_sum > 0,
    )
    ratio = np.divide(
        coadd_variance, predicted_variance,
        out=np.full(image.shape, np.nan, dtype=np.float64),
        where=(np.isfinite(coadd_variance) & (coadd_variance > 0)
               & np.isfinite(predicted_variance) & (predicted_variance > 0)),
    )
    coadd_no_data = 1 << coadd_bits[NO_DATA_NAME]
    coadd_exclusions = tuple(name for name in sorted(strict_mask_names_used) if name in coadd_bits)
    coadd_excluded = bitmask(coadd_exclusions, coadd_bits)
    coadd_available = (np.isfinite(image) & np.isfinite(coadd_variance)
                       & (coadd_variance > 0) & ((coadd_mask & coadd_no_data) == 0))
    clean_coadd = coadd_available & ((coadd_mask & coadd_excluded) == 0)
    if not warp_footprint.any():
        raise ValueError("approved warp archive has an empty footprint")

    center_sky_x_raw, center_sky_y_raw = coadd_wcs.world_to_pixel_values(
        float(request["ra_deg"]), float(request["dec_deg"]))
    center_sky_x, center_sky_y = float(np.asarray(center_sky_x_raw)), float(np.asarray(center_sky_y_raw))
    kpc_per_pixel = float(
        FlatLambdaCDM(H0=float(sample["h0"]), Om0=float(sample["omega_m"]), Tcmb0=0 * u.K)
        .angular_diameter_distance(float(provenance["cluster"]["z_spec_bcg"]))
        .to_value(u.kpc) * 0.168 / (180.0 * 3600.0 / math.pi)
    )
    grid_y, grid_x = np.indices(image.shape)
    radius_kpc = np.hypot((grid_x + 0.5 - center_sky_x) * kpc_per_pixel,
                          (grid_y + 0.5 - center_sky_y) * kpc_per_pixel)
    radial_rows = []
    for inner, outer in RADIUS_BINS_KPC:
        radial = (radius_kpc >= inner) & (radius_kpc < outer)
        primary = radial & coadd_available & (n_valid >= 3) & np.isfinite(ratio) & (ratio > 0)
        strict = radial & clean_coadd & (n_clean >= 3) & np.isfinite(ratio) & (ratio > 0)
        radial_rows.append({
            "inner_kpc": inner,
            "outer_kpc": outer,
            "geometric_pixels": int(radial.sum()),
            "primary_pixels": int(primary.sum()),
            "primary_coverage_fraction": float(primary.sum() / radial.sum()) if radial.any() else 0.0,
            "primary_median_n_valid_visits": float(np.median(n_valid[primary])) if primary.any() else None,
            "primary_coadd_to_inverse_variance_ratio": quantile_summary(ratio[primary]),
            "strict_pixels": int(strict.sum()),
            "strict_coverage_fraction": float(strict.sum() / radial.sum()) if radial.any() else 0.0,
            "strict_median_n_valid_visits": float(np.median(n_valid[strict])) if strict.any() else None,
            "strict_median_n_clean_visits": float(np.median(n_clean[strict])) if strict.any() else None,
            "strict_coadd_to_inverse_variance_ratio": quantile_summary(ratio[strict]),
        })

    coadd_center_x_raw, coadd_center_y_raw = coadd_wcs.world_to_pixel_values(
        float(request["ra_deg"]), float(request["dec_deg"]))
    coadd_center_x = float(np.asarray(coadd_center_x_raw))
    coadd_center_y = float(np.asarray(coadd_center_y_raw))
    ix, iy = int(round(coadd_center_x)), int(round(coadd_center_y))
    summary = {
        "schema": "astrafek-icl-warp-pilot-coadd-variance-v1",
        "status": "complete_exploratory",
        "science_measurement": False,
        "hsc_network_used_by_analysis": False,
        "image_pixels_read": True,
        "target": {
            "center_id": request["center_id"],
            "filter": request["filter"],
            "redshift": float(provenance["cluster"]["z_spec_bcg"]),
            "coadd_center_pixel_xy": [float(coadd_center_x), float(coadd_center_y)],
            "coadd_center_nearest_pixel": [ix, iy],
            "kpc_per_pixel": kpc_per_pixel,
            "radius_200_300_kpc_arcsec": [float(200 * 0.168 / kpc_per_pixel),
                                           float(300 * 0.168 / kpc_per_pixel)],
        },
        "warp_archive": {
            "bytes": int((DOWNLOAD_DIR / "warp_pilot.tar").stat().st_size),
            "expected_visits": sorted(expected_visits),
            "n_visits": len(archive_records),
            "exact_visit_set": True,
            "per_visit": sorted(archive_records, key=lambda row: row["visit"]),
            "median_n_valid_visits_per_pixel_in_warp_footprint": float(np.median(n_valid[n_valid > 0])),
            "median_n_clean_visits_per_pixel_in_warp_footprint": float(np.median(n_clean[n_clean > 0])),
            "footprint_pixels": int(warp_footprint.sum()),
            "coadd_available_fraction_in_warp_footprint": float(np.mean(coadd_available[warp_footprint])),
            "at_least_6_valid_visit_fraction_in_warp_footprint": float(np.mean(n_valid[warp_footprint] >= 6)),
        },
        "coadd_variance_comparison": {
            "predictor": "1 / sum(1 / warp_variance) over finite, positive-variance warp pixels not marked MP_NO_DATA",
            "ratio": "delivered coadd VARIANCE / inverse-variance-combined warp variance",
            "minimum_valid_visits": 3,
            "strict_mask_exclusions": sorted(strict_mask_names_used),
            "coadd_strict_mask_exclusions_used": list(coadd_exclusions),
            "radial_bins": radial_rows,
        },
        "inputs_sha256": provenance["hashes"],
        "limitations": [
            "This is one HSC-I2 target position, not the 800-cluster sample or a field-level validation.",
            "Two requested visit warp files contain no finite image/variance pixels; one more has partial spatial coverage.",
            "The inverse-variance reconstruction uses delivered warp VARIANCE planes and assumes unit relative weights; the close agreement is a local coadd-product consistency check, not an independent noise calibration.",
            "The comparison does not establish source-Poisson variance, effective gain, correlated resampling noise, sky-subtraction covariance, or an ICL estimator uncertainty.",
            "Mask-clean radial coverage is incomplete and the pilot footprint is limited to 2 arcmin square.",
        ],
    }
    output_dir.mkdir(parents=True, exist_ok=False)
    summary_path = output_dir / "summary.json"
    summary_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2,
                                       allow_nan=False) + "\n", encoding="utf-8")
    radial_path = output_dir / "radial_comparison.tsv"
    with radial_path.open("x", encoding="utf-8", newline="") as stream:
        fieldnames = [
            "inner_kpc", "outer_kpc", "geometric_pixels", "primary_pixels",
            "primary_coverage_fraction", "primary_median_n_valid_visits",
            "primary_q16_ratio", "primary_median_ratio", "primary_q84_ratio",
            "primary_fraction_within_1pct", "strict_pixels", "strict_coverage_fraction",
            "strict_median_n_clean_visits", "strict_q16_ratio", "strict_median_ratio",
            "strict_q84_ratio", "strict_fraction_within_1pct",
        ]
        writer = csv.DictWriter(stream, fieldnames=fieldnames, delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        for row in radial_rows:
            primary_q = row["primary_coadd_to_inverse_variance_ratio"]["q16_median_q84"] or [None] * 3
            strict_q = row["strict_coadd_to_inverse_variance_ratio"]["q16_median_q84"] or [None] * 3
            writer.writerow({
                "inner_kpc": row["inner_kpc"], "outer_kpc": row["outer_kpc"],
                "geometric_pixels": row["geometric_pixels"], "primary_pixels": row["primary_pixels"],
                "primary_coverage_fraction": row["primary_coverage_fraction"],
                "primary_median_n_valid_visits": row["primary_median_n_valid_visits"],
                "primary_q16_ratio": primary_q[0], "primary_median_ratio": primary_q[1],
                "primary_q84_ratio": primary_q[2],
                "primary_fraction_within_1pct": row["primary_coadd_to_inverse_variance_ratio"]["fraction_within_1pct_of_unity"],
                "strict_pixels": row["strict_pixels"], "strict_coverage_fraction": row["strict_coverage_fraction"],
                "strict_median_n_clean_visits": row["strict_median_n_clean_visits"],
                "strict_q16_ratio": strict_q[0], "strict_median_ratio": strict_q[1],
                "strict_q84_ratio": strict_q[2],
                "strict_fraction_within_1pct": row["strict_coadd_to_inverse_variance_ratio"]["fraction_within_1pct_of_unity"],
            })
    print(json.dumps({"status": summary["status"], "science_measurement": False,
                      "summary": str(summary_path), "n_visits": len(archive_records),
                      "radial_comparison": str(radial_path)}, indent=2))
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    run(args.output_dir)


if __name__ == "__main__":
    main()
