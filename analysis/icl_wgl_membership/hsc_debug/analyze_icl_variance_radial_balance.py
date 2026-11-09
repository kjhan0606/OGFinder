#!/usr/bin/env python3
"""Compare delivered coadd-variance profiles in supported target/control pairs.

This is an offline diagnostic on already downloaded HSC coadds. It measures
the delivered VARIANCE plane by proper radius, compares matched target and
control fields, and tests whether their variance ratio covaries with the saved
ICL contrast. It does not infer source-Poisson variance or promote a science
measurement.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import sys

import numpy as np
from astropy.cosmology import FlatLambdaCDM
from astropy.io import fits

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from download_pilot import file_sha256, validate_fits
from evaluate_prestack_pair_support import pair_block_labels
from audit_icl_target_control_outer_sky import partial_pearson


DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
SUPPORT = DATA / "icl_satellite_supported_94_association_2026-11-05_01/supported_pairs.tsv"
TARGET_MANIFEST = DATA / "download_01/download_manifest.json"
TARGET_STATUS = DATA / "download_01/download_status.json"
CONTROL_MANIFEST = DATA / "control_download_01/download_manifest.json"
CONTROL_STATUS = DATA / "control_download_01/download_status.json"
OUTPUT = DATA / "icl_variance_radial_balance_supported94_2026-11-06_01"
OUTPUT_ALL = DATA / "icl_variance_radial_balance_full783_2026-11-06_01"
EXPECTED_SUPPORT_SHA256 = "79f7c4c6b18b29ad837e91da4a4359433b09a564cad62efaf60b4cce655685da"
RADIAL_EDGES_KPC = (0, 50, 75, 100, 150, 200, 300, 400, 500, 600, 900)
N_BOOT = 5000
SEED = 20261106
BAD_AND_STAR_BIT_VALUES = [1, 2, 4, 8, 16, 128, 256, 512, 1024, 4096,
                           8192, 16384, 32768, 65536]


def sha256(path: Path) -> str:
    return file_sha256(path)


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f"refusing to write empty table: {path}")
    with path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t", extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def variance_profile(path: Path, frame: dict, download: dict,
                     redshift: float, cosmology: FlatLambdaCDM) -> list[dict]:
    if not path.is_file() or path.stat().st_size != int(download["bytes"]):
        raise ValueError(f"missing image or saved byte-size mismatch: {path}")
    digest = sha256(path)
    if digest != download["sha256"]:
        raise ValueError(f"image SHA-256 differs from download record: {path}")
    info = validate_fits(path, frame)
    expected_filter = str(download.get("review", {}).get("coadd_header_filter", "")).lower()
    if info["coadd_header_filter"] != expected_filter:
        raise ValueError(f"filter mismatch for {path}: {info['coadd_header_filter']} vs {expected_filter}")
    with fits.open(path, memmap=False, checksum=True) as hdus:
        variance = np.asarray(hdus[info["variance_hdu"]].data, dtype=np.float64)
        bits = np.asarray(hdus[info["mask_hdu"]].data, dtype=np.uint64)
    if variance.shape != bits.shape:
        raise ValueError(f"variance and mask arrays differ in shape: {path}")
    if info["input_unit"] == "instrumental_flux/pixel":
        conversion = 3.631e12 / float(info["fluxmag0"]) / float(info["pixel_scale_arcsec"]) ** 2
    elif info["input_unit"] == "nJy/pixel":
        conversion = 1.0 / float(info["pixel_scale_arcsec"]) ** 2
    else:
        raise ValueError(f"unsupported science unit for variance conversion: {info['input_unit']}")
    variance_njy = variance * conversion ** 2
    exclude_mask = 0
    for bit in BAD_AND_STAR_BIT_VALUES:
        exclude_mask |= bit
    valid = np.isfinite(variance_njy) & (variance_njy > 0) & ((bits & exclude_mask) == 0)
    kpc_per_pixel = float(cosmology.kpc_proper_per_arcmin(redshift).value) \
        * float(info["pixel_scale_arcsec"]) / 60.0
    yy, xx = np.ogrid[:variance.shape[0], :variance.shape[1]]
    center_x, center_y = map(float, info["center_pixel_xy"])
    radius = np.hypot(xx - center_x, yy - center_y) * kpc_per_pixel
    result = []
    for lo, hi in zip(RADIAL_EDGES_KPC[:-1], RADIAL_EDGES_KPC[1:]):
        annulus = (radius >= lo) & (radius < hi)
        selected = variance_njy[annulus & valid]
        enough_pixels = selected.size >= 50
        result.append({
            "r_min_kpc": lo,
            "r_max_kpc": hi,
            "median_variance_njy2_arcsec4": float(np.median(selected)) if enough_pixels else None,
            "q16_variance_njy2_arcsec4": float(np.quantile(selected, 0.16)) if enough_pixels else None,
            "q84_variance_njy2_arcsec4": float(np.quantile(selected, 0.84)) if enough_pixels else None,
            "median_meets_50_pixel_floor": enough_pixels,
            "n_geometric_pixels": int(annulus.sum()),
            "n_valid_pixels": int(selected.size),
            "valid_fraction": float(selected.size / annulus.sum()),
            "pixel_scale_arcsec": float(info["pixel_scale_arcsec"]),
            "conversion_njy_arcsec2_per_native_unit": float(conversion),
            "file_sha256": digest,
        })
    return result


def block_bootstrap(rows: list[dict], replicas: int, seed: int) -> tuple[dict, np.ndarray]:
    target_coords = [{"ra_deg": float(r["target_ra_deg"]), "dec_deg": float(r["target_dec_deg"])} for r in rows]
    control_coords = [{"ra_deg": float(r["control_ra_deg"]), "dec_deg": float(r["control_dec_deg"])} for r in rows]
    _, labels, n_target_cells, n_control_cells, n_blocks, largest = pair_block_labels(
        target_coords, control_coords, np.ones(len(rows), dtype=bool)
    )
    blocks = np.unique(labels)
    if n_blocks < 6:
        raise ValueError(f"only {n_blocks} connected spatial blocks")
    members = {label: np.flatnonzero(labels == label) for label in blocks}
    log_ratio = np.asarray([float(row["log_target_control_variance_ratio"]) for row in rows])
    point = float(np.median(log_ratio))
    rng = np.random.default_rng(seed)
    boot_median = np.empty(replicas)
    for rep in range(replicas):
        draw = rng.choice(blocks, size=len(blocks), replace=True)
        idx = np.concatenate([members[label] for label in draw])
        boot_median[rep] = np.median(log_ratio[idx])

    icl_rows = [row for row in rows if row.get("icl_delta_200_300") not in (None, "")]
    corr_point = None
    corr_interval = [None, None]
    corr_valid_fraction = 0.0
    corr_blocks = 0
    if len(icl_rows) >= 6:
        corr_t = [{"ra_deg": float(r["target_ra_deg"]), "dec_deg": float(r["target_dec_deg"])} for r in icl_rows]
        corr_c = [{"ra_deg": float(r["control_ra_deg"]), "dec_deg": float(r["control_dec_deg"])} for r in icl_rows]
        _, corr_labels, _, _, corr_blocks, _ = pair_block_labels(
            corr_t, corr_c, np.ones(len(icl_rows), dtype=bool))
        corr_groups = np.unique(corr_labels)
        if corr_blocks >= 6:
            corr_members = {label: np.flatnonzero(corr_labels == label) for label in corr_groups}
            x = np.asarray([float(row["icl_delta_200_300"]) for row in icl_rows])
            y = np.asarray([float(row["log_target_control_variance_ratio"]) for row in icl_rows])
            z = np.asarray([float(row["redshift"]) for row in icl_rows])
            logr = np.asarray([float(row["log_richness"]) for row in icl_rows])
            corr_point = partial_pearson(x, y, [z, logr])
            corr_rng = np.random.default_rng(seed + 500000)
            boot_corr = np.empty(replicas)
            for rep in range(replicas):
                draw = corr_rng.choice(corr_groups, size=len(corr_groups), replace=True)
                idx = np.concatenate([corr_members[label] for label in draw])
                boot_corr[rep] = partial_pearson(x[idx], y[idx], [z[idx], logr[idx]])
            finite_corr = boot_corr[np.isfinite(boot_corr)]
            if finite_corr.size:
                corr_interval = [float(v) for v in np.quantile(finite_corr, [0.025, 0.975])]
                corr_valid_fraction = float(finite_corr.size / replicas)
    summary = {
        "n_pairs": len(rows),
        "n_pairs_with_icl_support": len(icl_rows),
        "n_connected_5deg_blocks": int(n_blocks),
        "largest_block_pairs": int(largest),
        "n_target_sky_cells": int(n_target_cells),
        "n_control_sky_cells": int(n_control_cells),
        "median_log_target_control_variance_ratio": point,
        "target_control_variance_ratio_median": float(np.exp(point)),
        "median_log_ratio_block_bootstrap_95": [float(v) for v in np.quantile(boot_median, [0.025, 0.975])],
        "icl_partial_correlation_controlling_z_logrichness": None if corr_point is None else float(corr_point),
        "icl_partial_correlation_block_bootstrap_95": corr_interval,
        "correlation_connected_5deg_blocks": int(corr_blocks),
        "correlation_bootstrap_valid_fraction": corr_valid_fraction,
    }
    return summary, boot_median


def analyze(output: Path, replicas: int = N_BOOT, scope: str = "supported94") -> dict:
    output = output.resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite {output}")
    if replicas < 1000:
        raise ValueError("use at least 1000 connected-block bootstrap replicates")
    if sha256(SUPPORT) != EXPECTED_SUPPORT_SHA256:
        raise ValueError("pinned ICL–satellite support table changed")
    target_doc, control_doc = read_json(TARGET_MANIFEST), read_json(CONTROL_MANIFEST)
    target_state, control_state = read_json(TARGET_STATUS), read_json(CONTROL_STATUS)
    if any(state.get("status") != "complete" for state in (target_state, control_state)):
        raise ValueError("target/control image download states are not complete")
    targets = {str(row["id"]): row for row in target_doc["clusters"]}
    controls = {str(row["id"]): row for row in control_doc["clusters"]}
    target_downloads = {str(row["id"]): row for row in target_state["completed"]}
    control_downloads = {str(row["id"]): row for row in control_state["completed"]}
    controls_by_target = {str(row["target_cluster_id"]): row for row in control_doc["clusters"]}
    target_base = Path(target_doc["downloaded_images"])
    control_base = Path(control_doc["downloaded_images"])
    if not target_base.is_absolute():
        target_base = TARGET_MANIFEST.parent / target_base
    if not control_base.is_absolute():
        control_base = CONTROL_MANIFEST.parent / control_base
    support_rows = read_tsv(SUPPORT)
    if len(support_rows) != 94 or any(row["science_measurement"].lower() != "false" for row in support_rows):
        raise ValueError("expected the existing 94-row exploratory support table")
    support_by_id = {str(row["cluster_id"]): row for row in support_rows}
    if scope == "supported94":
        selected_pairs = support_rows
    elif scope == "all783":
        selected_pairs = []
        for ident, target in targets.items():
            control = controls_by_target.get(ident)
            if control is None:
                continue
            control_id = str(control["id"])
            t_download = target_downloads.get(ident)
            c_download = control_downloads.get(control_id)
            if t_download is None or c_download is None:
                continue
            target_filter = str(t_download["review"]["coadd_header_filter"]).lower()
            control_filter = str(c_download["review"]["coadd_header_filter"]).lower()
            if target_filter not in ("i", "i2") or target_filter != control_filter:
                continue
            selected_pairs.append({
                "cluster_id": ident,
                "control_id": control_id,
                "filter": target_filter,
                "icl_delta_200_300": support_by_id.get(ident, {}).get("icl_delta_200_300", ""),
                "satellite_delta_50_75": support_by_id.get(ident, {}).get("satellite_delta_50_75", ""),
                "redshift": str(target["z"]),
                "richness": str(target["richness"]),
                "science_measurement": "false",
            })
        if len(selected_pairs) != 783:
            raise ValueError(f"expected 783 same-filter target/control pairs, got {len(selected_pairs)}")
    else:
        raise ValueError(f"unknown sample scope: {scope}")
    cosmology = FlatLambdaCDM(H0=float(target_doc["h0"]), Om0=float(target_doc["omega_m"]), Tcmb0=0)
    profile_rows = []
    per_pair: dict[tuple[str, str], dict] = {}
    for row in selected_pairs:
        ident = str(row["cluster_id"])
        control_id = str(row["control_id"])
        band = str(row["filter"])
        target = targets.get(ident)
        control = controls.get(control_id)
        if target is None or control is None or controls_by_target.get(ident, {}).get("id") != control_id:
            raise ValueError(f"matched target/control identity mismatch for {ident}")
        if str(target_downloads[ident]["review"]["coadd_header_filter"]).lower() != band \
                or str(control_downloads[control_id]["review"]["coadd_header_filter"]).lower() != band:
            raise ValueError(f"target/control filter mismatch for {ident}")
        if not math.isclose(float(target["z"]), float(control["z"]), abs_tol=1e-12):
            raise ValueError(f"control angular-scale redshift differs from target for {ident}")
        target_path = target_base / target_downloads[ident]["file"]
        control_path = control_base / control_downloads[control_id]["file"]
        target_prof = variance_profile(target_path, target, target_downloads[ident], float(target["z"]), cosmology)
        control_prof = variance_profile(control_path, control, control_downloads[control_id], float(target["z"]), cosmology)
        per_pair[(ident, band)] = {
            "support": row, "target": target, "control": control,
            "target_profile": target_prof, "control_profile": control_prof,
            "target_sha256": target_downloads[ident]["sha256"],
            "control_sha256": control_downloads[control_id]["sha256"],
        }
        for role, profile in (("target", target_prof), ("matched_control", control_prof)):
            for value in profile:
                profile_rows.append({
                    "cluster_id": ident, "control_id": control_id, "filter": band,
                    "has_icl_satellite_support": ident in support_by_id,
                    "sample_role": role, "icl_delta_200_300": row["icl_delta_200_300"],
                    "satellite_delta_50_75": row["satellite_delta_50_75"],
                    "redshift": row["redshift"], "richness": row["richness"],
                    **value,
                })

    comparison_rows = []
    summary_by_filter = {}
    for band in ("i", "i2"):
        band_pairs = [(key, value) for key, value in per_pair.items() if key[1] == band]
        expected_counts = {"supported94": {"i": 26, "i2": 68},
                           "all783": {"i": 302, "i2": 481}}[scope]
        if len(band_pairs) != expected_counts[band]:
            raise ValueError(f"unexpected {band} supported-pair count: {len(band_pairs)}")
        band_results = []
        for radial_index, (lo, hi) in enumerate(zip(RADIAL_EDGES_KPC[:-1], RADIAL_EDGES_KPC[1:])):
            rows = []
            for key, values in band_pairs:
                target_metric = values["target_profile"][radial_index]
                control_metric = values["control_profile"][radial_index]
                target_var = target_metric["median_variance_njy2_arcsec4"]
                control_var = control_metric["median_variance_njy2_arcsec4"]
                if target_var is None or control_var is None:
                    continue
                item = {
                    "cluster_id": key[0], "control_id": values["support"]["control_id"],
                    "filter": band, "r_min_kpc": lo, "r_max_kpc": hi,
                    "has_icl_satellite_support": key[0] in support_by_id,
                    "icl_delta_200_300": values["support"]["icl_delta_200_300"],
                    "satellite_delta_50_75": values["support"]["satellite_delta_50_75"],
                    "redshift": values["support"]["redshift"],
                    "log_richness": math.log(float(values["support"]["richness"])),
                    "target_ra_deg": values["target"]["ra_deg"],
                    "target_dec_deg": values["target"]["dec_deg"],
                    "control_ra_deg": values["control"]["ra_deg"],
                    "control_dec_deg": values["control"]["dec_deg"],
                    "target_variance": target_var, "control_variance": control_var,
                    "log_target_control_variance_ratio": math.log(target_var / control_var),
                    "target_valid_fraction": target_metric["valid_fraction"],
                    "control_valid_fraction": control_metric["valid_fraction"],
                    "target_n_valid": target_metric["n_valid_pixels"],
                    "control_n_valid": control_metric["n_valid_pixels"],
                }
                rows.append(item)
                comparison_rows.append(item)
            target_coords = [{"ra_deg": float(r["target_ra_deg"]), "dec_deg": float(r["target_dec_deg"])}
                             for r in rows]
            control_coords = [{"ra_deg": float(r["control_ra_deg"]), "dec_deg": float(r["control_dec_deg"])}
                              for r in rows]
            _, labels, _, _, n_blocks, _ = pair_block_labels(
                target_coords, control_coords, np.ones(len(rows), dtype=bool))
            for item, label in zip(rows, labels):
                item["pair_block_5deg"] = str(label)
            stats, _ = block_bootstrap(rows, replicas, SEED + radial_index + (0 if band == "i" else 1000))
            stats.update({"filter": band, "r_min_kpc": lo, "r_max_kpc": hi,
                          "median_target_variance": float(np.median([r["target_variance"] for r in rows])),
                          "median_control_variance": float(np.median([r["control_variance"] for r in rows])),
                          "median_target_valid_fraction": float(np.median([r["target_valid_fraction"] for r in rows])),
                          "median_control_valid_fraction": float(np.median([r["control_valid_fraction"] for r in rows]))})
            band_results.append(stats)
        summary_by_filter[band] = {
            "n_pairs": len(band_pairs),
            "n_pairs_with_icl_satellite_support": sum(1 for _, values in band_pairs
                                                       if str(values["support"]["icl_delta_200_300"]) != ""),
            "radial_results": band_results,
        }

    output.mkdir(parents=True, exist_ok=False)
    write_tsv(output / "per_pair_radial_variance.tsv", profile_rows)
    write_tsv(output / "target_control_radial_comparison.tsv", comparison_rows)
    image_hash_rows = []
    for (ident, band), values in sorted(per_pair.items(), key=lambda x: (x[0][1], int(x[0][0]))):
        image_hash_rows.append({"cluster_id": ident, "control_id": values["support"]["control_id"],
                                "filter": band, "target_image_sha256": values["target_sha256"],
                                "control_image_sha256": values["control_sha256"]})
    write_tsv(output / "image_hashes.tsv", image_hash_rows)
    paths = [SUPPORT, TARGET_MANIFEST, TARGET_STATUS, CONTROL_MANIFEST, CONTROL_STATUS]
    result = {
        "schema": "astrafek-icl-radial-variance-balance-v2",
        "status": "complete_diagnostic",
        "science_measurement": False,
        "network_used": False,
        "new_hsc_data_requested": False,
        "science_image_values_used": False,
        "variance_and_mask_pixels_read": True,
        "sample_scope": scope,
        "purpose": "check whether delivered coadd-variance radial structure differs between matched targets and controls and, where available, tracks saved ICL contrast",
        "n_pairs": len(selected_pairs),
        "n_pairs_with_icl_satellite_support": len(support_rows),
        "radial_edges_proper_kpc": list(RADIAL_EDGES_KPC),
        "variance_definition": "median HSC VARIANCE plane in nJy^2 arcsec^-4 after excluding non-finite/non-positive pixels and reviewed bad/bright-star mask bits; no object-detection masks applied",
        "bootstrap": {"replicates": replicas, "method": "resample connected 5-degree target/control sky blocks within filter",
                      "seed": SEED},
        "filters": summary_by_filter,
        "limitations": [
            "delivered coadd variance omits the unknown effective-gain term for source shot noise in synthetic injections",
            "median variance maps are not the covariance of the masked radial ICL estimator",
            "mask-valid fractions can differ between target and control and are reported alongside the variance medians",
            "ICL correlation is restricted to the 94-pair support; connected equal-area cells are approximate, not validated HSC fields",
            "no correction to the measured ICL profile or its science gate is applied",
        ],
        "inputs_sha256": {str(path): sha256(path) for path in paths},
        "output_sha256": {name: sha256(output / name) for name in
                          ("per_pair_radial_variance.tsv", "target_control_radial_comparison.tsv", "image_hashes.tsv")},
    }
    (output / "summary.json").write_text(json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n")
    print(json.dumps({"status": result["status"], "science_measurement": False,
                      "sample_scope": scope, "n_pairs": result["n_pairs"],
                      "summary": str(output / "summary.json")}, indent=2))
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scope", choices=("supported94", "all783"), default="supported94")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--bootstrap", type=int, default=N_BOOT)
    args = parser.parse_args()
    output = args.output or (OUTPUT if args.scope == "supported94" else OUTPUT_ALL)
    analyze(output, args.bootstrap, args.scope)


if __name__ == "__main__":
    main()
