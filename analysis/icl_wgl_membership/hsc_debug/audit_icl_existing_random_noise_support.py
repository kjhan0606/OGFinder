#!/usr/bin/env python3
"""Audit whether the saved HSC random points sample ICL radial apertures.

This is an offline support diagnostic. It uses only the random-point rows
already saved for the 800 target centers; it does not query HSC or read image
pixels. HSC PDR3 documents ``i_pix_variance`` as pixel variance at each random
point, but sparse radial support cannot be treated as a pixel-level coadd
variance map.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from astropy.cosmology import FlatLambdaCDM


ROOT = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
CATALOG_DIR = ROOT / "catalogues_2026-10-13_01"
ANCILLARY_DIR = ROOT / "ancillary_2026-10-13_01"
REQUESTS = ANCILLARY_DIR / "hsc_catalogue_requests.json"
STATUS = CATALOG_DIR / "catalogue_status.json"
TARGET_MANIFEST = ROOT / "sample_manifest.json"
DEFAULT_OUTPUT = ROOT / "icl_existing_random_noise_support_2026-11-06_01"
EXPECTED_CATALOG_PLAN_SHA256 = "eda38489262eed7467746255673e9abf3fd2dadddfbdf4235ff5acd2bc46f504"
COSMOLOGY = FlatLambdaCDM(H0=70.0, Om0=0.3, Tcmb0=2.725)
ANNULI_KPC = ((50.0, 75.0), (75.0, 150.0), (150.0, 200.0),
              (200.0, 300.0), (300.0, 500.0))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_csv(path: Path) -> tuple[list[str], list[dict[str, str]]]:
    with path.open(encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream)
        names = [name.strip().removeprefix("#").strip() for name in (reader.fieldnames or [])]
        reader.fieldnames = names
        if len(names) != len(set(names)):
            raise ValueError(f"duplicate columns in {path}")
        rows = list(reader)
    return names, rows


def angular_separation_arcsec(ra: np.ndarray, dec: np.ndarray,
                              ra0: float, dec0: float) -> np.ndarray:
    rar, decr = np.deg2rad(ra), np.deg2rad(dec)
    ra0r, dec0r = math.radians(ra0), math.radians(dec0)
    cosine = (np.sin(decr) * math.sin(dec0r)
              + np.cos(decr) * math.cos(dec0r) * np.cos(rar - ra0r))
    return np.rad2deg(np.arccos(np.clip(cosine, -1.0, 1.0))) * 3600.0


def summarize(values: list[float]) -> dict:
    array = np.asarray(values, dtype=float)
    if array.size == 0:
        return {"n": 0, "median": None, "q16": None, "q84": None}
    return {"n": int(array.size), "median": float(np.median(array)),
            "q16": float(np.quantile(array, .16)),
            "q84": float(np.quantile(array, .84))}


def run(output: Path) -> dict:
    output = Path(output).resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite {output}")
    status = json.loads(STATUS.read_text(encoding="utf-8"))
    if status.get("status") != "complete" or status.get("plan_sha256") != EXPECTED_CATALOG_PLAN_SHA256:
        raise ValueError("saved ancillary catalog checkpoint is incomplete or changed")
    requests = json.loads(REQUESTS.read_text(encoding="utf-8"))
    manifest = json.loads(TARGET_MANIFEST.read_text(encoding="utf-8"))
    if len(requests) != 800 or len(manifest.get("clusters", [])) != 800:
        raise ValueError("expected random-point requests and target manifest for 800 centers")
    centers = {str(row["id"]): row for row in manifest["clusters"]}
    request_ids = {str(row["cluster_id"]) for row in requests}
    if request_ids != set(centers):
        raise ValueError("saved random-point request IDs do not match the 800-center sample")

    rows_by_center: dict[str, list[dict[str, str]]] = {key: [] for key in centers}
    file_hashes = {}
    file_row_total = 0
    for batch_index in range(80):
        path = CATALOG_DIR / f"random_{batch_index:03d}.csv"
        if not path.is_file():
            raise FileNotFoundError(path)
        names, rows = read_csv(path)
        required = {"target_cluster_id", "object_id", "ra", "dec", "adjust_density",
                    "i_pix_variance", "i_sky_mean", "i_sky_std", "i_inputcount_value",
                    "i_pixelflags_bright_objectcenter", "i_pixelflags_edge"}
        if not required.issubset(names):
            raise ValueError(f"saved random-point columns missing from {path.name}")
        if len(rows) != 2000:
            raise ValueError(f"expected 2,000 saved random points per ten-request batch in {path.name}; got {len(rows)}")
        for row in rows:
            ident = str(row["target_cluster_id"])
            if ident not in rows_by_center:
                raise ValueError(f"unrequested center in {path.name}: {ident}")
            rows_by_center[ident].append(row)
        file_row_total += len(rows)
        file_hashes[path.name] = sha256(path)
    if len(file_hashes) != 80 or file_row_total != 160000:
        raise ValueError("saved random-point batch inventory changed")

    detailed = []
    annulus_values: dict[str, dict[str, list[float]]] = {
        f"{lo:g}-{hi:g}": {"pix_variance": [], "sky_std": [], "inputcount": []}
        for lo, hi in ANNULI_KPC
    }
    summary_by_annulus = {}
    center_annulus_counts = {f"{lo:g}-{hi:g}": [] for lo, hi in ANNULI_KPC}
    center_summary: dict[str, dict] = {}
    for ident, source in centers.items():
        selected = rows_by_center[ident]
        if len(selected) != 200:
            raise ValueError(f"center {ident} has {len(selected)} rows; expected 200")
        ids = [row["object_id"] for row in selected]
        if len(ids) != len(set(ids)):
            raise ValueError(f"duplicate random point for center {ident}")
        ra = np.asarray([float(row["ra"]) for row in selected])
        dec = np.asarray([float(row["dec"]) for row in selected])
        separation = angular_separation_arcsec(ra, dec, float(source["ra_deg"]),
                                                float(source["dec_deg"]))
        if np.any(separation > 1800.01):
            raise ValueError(f"random point outside requested 1800-arcsec cone for center {ident}")
        z = float(source["z"])
        proper_kpc_per_arcsec = COSMOLOGY.kpc_proper_per_arcmin(z).to_value("kpc/arcsec")
        radius_kpc = separation * proper_kpc_per_arcsec
        center_row = {"cluster_id": ident, "z": z, "n_random_points": len(selected),
                      "proper_kpc_per_arcsec": proper_kpc_per_arcsec}
        for lo, hi in ANNULI_KPC:
            label = f"{lo:g}-{hi:g}"
            mask = (radius_kpc >= lo) & (radius_kpc < hi)
            in_ring = [row for row, keep in zip(selected, mask, strict=True) if keep]
            count = len(in_ring)
            center_annulus_counts[label].append(count)
            center_row[f"n_{label}_kpc"] = count
            for field, out_field in (("i_pix_variance", "pix_variance"),
                                     ("i_sky_std", "sky_std"),
                                     ("i_inputcount_value", "inputcount")):
                finite = []
                for row in in_ring:
                    try:
                        value = float(row[field])
                    except (TypeError, ValueError):
                        continue
                    if math.isfinite(value):
                        finite.append(value)
                annulus_values[label][out_field].extend(finite)
        detailed.append(center_row)
        center_summary[ident] = center_row

    for label, counts in center_annulus_counts.items():
        arr = np.asarray(counts, dtype=int)
        pix = annulus_values[label]["pix_variance"]
        sky = annulus_values[label]["sky_std"]
        inputs = annulus_values[label]["inputcount"]
        summary_by_annulus[label] = {
            "n_centers": len(arr), "total_random_points": int(arr.sum()),
            "points_per_center": summarize(arr.astype(float).tolist()),
            "centers_with_any_point": int(np.count_nonzero(arr)),
            "centers_with_at_least_10_points": int(np.count_nonzero(arr >= 10)),
            "centers_with_at_least_30_points": int(np.count_nonzero(arr >= 30)),
            "pix_variance_at_sampled_points": summarize(pix),
            "sky_std_at_sampled_points": summarize(sky),
            "inputcount_at_sampled_points": summarize(inputs),
        }

    output.mkdir(parents=True, exist_ok=False)
    table_path = output / "random_point_radial_support_by_cluster.tsv"
    with table_path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(detailed[0]), delimiter="\t",
                                lineterminator="\n")
        writer.writeheader()
        writer.writerows(sorted(detailed, key=lambda row: int(row["cluster_id"])))
    result = {
        "schema": "astrafek-icl-existing-random-noise-support-v1",
        "status": "complete_offline_support_audit",
        "science_measurement": False,
        "hsc_network_used": False,
        "image_pixels_read": False,
        "n_target_centers": 800,
        "n_saved_random_points": file_row_total,
        "random_points_per_center": 200,
        "original_query_radius_arcsec": 1800.0,
        "annuli_proper_kpc": [[lo, hi] for lo, hi in ANNULI_KPC],
        "cosmology": {"H0_km_s_Mpc": 70.0, "Om0": 0.3, "geometry_redshift": "cluster z from sample_manifest"},
        "radial_support": summary_by_annulus,
        "matched_control_random_points_available": False,
        "interpretation": "Existing 200-point, 30-arcmin random samples can characterize only the sampled target locations. Sparse annular counts and no corresponding matched-control random samples cannot provide paired target/control pixel-variance correction for the ICL stack.",
        "input_sha256": {str(path): sha256(path) for path in
                         (REQUESTS, STATUS, TARGET_MANIFEST)},
        "random_file_sha256": file_hashes,
        "center_table": str(table_path),
        "center_table_sha256": sha256(table_path),
        "not_returned_or_read": ["new HSC data", "coadd images", "image pixels", "shape rows", "photo-z", "object catalog rows"],
    }
    summary_path = output / "summary.json"
    summary_path.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    result["summary_sha256"] = sha256(summary_path)
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    print(json.dumps(run(args.output), ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
