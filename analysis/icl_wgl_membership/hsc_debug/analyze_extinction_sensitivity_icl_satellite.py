#!/usr/bin/env python3
"""Apply center-extinction corrections to saved ICL/satellite profiles.

This is an offline sensitivity analysis. It does not read image pixels or
request data from HSC CAS. It compares per-band SFD center corrections with a
PDR3 forced-catalogue a_i correction for target frames, and keeps the science
gate closed because extinction is not corrected pixel-by-pixel and other sky,
PSF, masking, and injection/recovery gates remain open.
"""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import sys

import numpy as np

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

from evaluate_prestack_pair_support import pair_block_labels
from fetch_ancillary_catalogues import catalogue_reader


DEFAULT_DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_json(path):
    return json.loads(Path(path).read_text(encoding="utf-8"))


def read_tsv(path):
    with Path(path).open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def write_tsv(path, rows, fields):
    with Path(path).open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter="\t",
                                extrasaction="raise")
        writer.writeheader()
        writer.writerows(rows)


def angular_separation_arcsec(ra1, dec1, ra2, dec2):
    r1, d1, r2, d2 = np.deg2rad([ra1, dec1, ra2, dec2])
    sd = math.sin((d2 - d1) / 2)
    sr = math.sin((r2 - r1) / 2)
    hav = sd * sd + math.cos(d1) * math.cos(d2) * sr * sr
    return 2 * math.asin(math.sqrt(min(1.0, max(0.0, hav)))) * 206264.806247


def read_forced_extinction(catalogue_dir, status):
    """Read only checksum-verified PDR3 forced-source extinction results."""
    if status.get("status") != "complete":
        raise ValueError("PDR3 catalogue query checkpoint is not complete")
    output = {}
    hashes = {}
    for name, record in status.get("jobs", {}).items():
        if not name.startswith("extinction_"):
            continue
        if record.get("status") != "complete":
            raise ValueError(f"PDR3 forced query incomplete: {name}")
        path = catalogue_dir / (name + ".csv")
        if not path.is_file() or sha256(path) != record.get("sha256"):
            raise ValueError(f"PDR3 forced query checksum mismatch: {path}")
        hashes[str(path)] = record["sha256"]
        with path.open(encoding="utf-8-sig", newline="") as stream:
            reader = catalogue_reader(stream)
            expected = {"target_cluster_id", "object_id", "ra", "dec", "a_i"}
            if not expected.issubset(reader.fieldnames or []):
                raise ValueError(f"Unexpected PDR3 columns in {path.name}")
            for row in reader:
                cluster_id = str(row["target_cluster_id"])
                values = tuple(float(row[key]) for key in ("ra", "dec", "a_i"))
                if not all(math.isfinite(value) for value in values) or values[2] < 0:
                    raise ValueError(f"Invalid PDR3 forced extinction row in {path.name}")
                output.setdefault(cluster_id, []).append(values)
    if not output:
        raise ValueError("No verified PDR3 forced-source extinction rows found")
    return output, hashes


def choose_nearest_forced(cluster_rows, pdr3_rows):
    chosen = {}
    audit = {}
    for row in cluster_rows:
        cluster_id = str(row["id"])
        candidates = pdr3_rows.get(cluster_id, [])
        near = sorted((angular_separation_arcsec(
            float(row["ra_deg"]), float(row["dec_deg"]), ra, dec), ai)
            for ra, dec, ai in candidates)
        if not near or near[0][0] > 3.0:
            raise ValueError(f"No PDR3 forced object within 3 arcsec of target {cluster_id}")
        if len(near) > 1 and near[1][0] - near[0][0] < 1e-5:
            raise ValueError(f"Ambiguous nearest PDR3 forced object for target {cluster_id}")
        chosen[cluster_id] = near[0][1]
        audit[cluster_id] = {"separation_arcsec": near[0][0], "a_i": near[0][1],
                             "n_candidates_within_3arcsec": sum(item[0] <= 3 for item in near)}
    return chosen, audit


def partial_pearson(x, y, z, richness):
    arrays = [np.asarray(v, dtype=float) for v in (x, y, z, richness)]
    if len(arrays[0]) < 5 or any(not np.all(np.isfinite(v)) for v in arrays):
        return np.nan
    design = np.column_stack((np.ones(len(arrays[0])), arrays[2], arrays[3]))
    if np.linalg.matrix_rank(design) < 3:
        return np.nan
    rx = arrays[0] - design @ np.linalg.lstsq(design, arrays[0], rcond=None)[0]
    ry = arrays[1] - design @ np.linalg.lstsq(design, arrays[1], rcond=None)[0]
    denom = np.linalg.norm(rx) * np.linalg.norm(ry)
    return float(rx @ ry / denom) if denom > 0 else np.nan


def percentile(values, q):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    return float(np.percentile(values, q)) if len(values) else None


def block_bootstrap_mean(values, labels, replicas, seed):
    labels = np.asarray(labels, dtype=str)
    blocks = np.unique(labels)
    if len(blocks) < 6:
        raise ValueError("At least six connected spatial blocks are required")
    per_block = {block: np.flatnonzero(labels == block) for block in blocks}
    rng = np.random.default_rng(seed)
    result = np.empty(replicas, dtype=float)
    for i in range(replicas):
        draw = rng.choice(blocks, size=len(blocks), replace=True)
        index = np.concatenate([per_block[item] for item in draw])
        result[i] = float(np.mean(np.asarray(values)[index]))
    return {"ci95_low": percentile(result, 2.5), "ci95_high": percentile(result, 97.5),
            "replicas": replicas, "n_blocks": int(len(blocks))}


def block_bootstrap_correlations(rows, labels, replicas, seed):
    labels = np.asarray(labels, dtype=str)
    blocks = np.unique(labels)
    if len(blocks) < 6:
        raise ValueError("At least six connected spatial blocks are required")
    per_block = {block: np.flatnonzero(labels == block) for block in blocks}
    xraw = np.array([row["icl_delta_raw_njy_arcsec2"] for row in rows], dtype=float)
    x_sfd = np.array([row["icl_delta_sfd_corrected_njy_arcsec2"] for row in rows], dtype=float)
    x_pdr3 = np.array([row["icl_delta_pdr3_target_corrected_njy_arcsec2"] for row in rows], dtype=float)
    yraw = np.array([row["satellite_delta_raw_njy_arcsec2"] for row in rows], dtype=float)
    y_sfd = np.array([row["satellite_delta_sfd_corrected_njy_arcsec2"] for row in rows], dtype=float)
    y_pdr3 = np.array([row["satellite_delta_pdr3_target_corrected_njy_arcsec2"] for row in rows], dtype=float)
    z = np.array([row["z"] for row in rows], dtype=float)
    richness = np.array([row["richness"] for row in rows], dtype=float)
    point = {
        "raw": partial_pearson(xraw, yraw, z, richness),
        "sfd": partial_pearson(x_sfd, y_sfd, z, richness),
        "pdr3_target": partial_pearson(x_pdr3, y_pdr3, z, richness),
    }
    rng = np.random.default_rng(seed)
    draws = {key: np.full(replicas, np.nan) for key in point}
    for replica in range(replicas):
        selected = rng.choice(blocks, size=len(blocks), replace=True)
        indices = np.concatenate([per_block[item] for item in selected])
        draws["raw"][replica] = partial_pearson(xraw[indices], yraw[indices], z[indices], richness[indices])
        draws["sfd"][replica] = partial_pearson(x_sfd[indices], y_sfd[indices], z[indices], richness[indices])
        draws["pdr3_target"][replica] = partial_pearson(x_pdr3[indices], y_pdr3[indices], z[indices], richness[indices])
    output = {}
    for key in point:
        output[key] = {"partial_pearson_r": point[key],
                       "block_bootstrap_ci95_low": percentile(draws[key], 2.5),
                       "block_bootstrap_ci95_high": percentile(draws[key], 97.5),
                       "valid_bootstrap_replicas": int(np.isfinite(draws[key]).sum())}
    return output, int(len(blocks))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data", type=Path, default=DEFAULT_DATA)
    parser.add_argument("--work", type=Path, required=True)
    parser.add_argument("--replicas", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=20261019)
    args = parser.parse_args()
    data, work = args.data.resolve(), args.work.resolve()
    if work.exists():
        parser.error(f"Refusing to overwrite existing output: {work}")
    if args.replicas < 1000:
        parser.error("Use at least 1000 block-bootstrap replicas")

    sample_path = data / "sample_manifest.json"
    candidates_path = data / "ancillary_2026-10-13_01/extinction_candidates.tsv"
    extinction_validation_path = data / "ancillary_2026-10-13_01/extinction_validation.json"
    catalogue_dir = data / "catalogues_2026-10-13_01"
    catalogue_status_path = catalogue_dir / "catalogue_status.json"
    control_path = data / "control_download_01/download_manifest.json"
    annulus_dir = data / "icl_annulus_sensitivity_2026-10-18_02"
    annulus_summary_path = annulus_dir / "summary.json"
    annulus_path = annulus_dir / "selected_pair_profiles.tsv"
    association_dir = data / "icl_satellite_cut_sensitivity_2026-10-18_01"
    association_summary_path = association_dir / "summary.json"
    association_path = association_dir / "matched_cluster_values.tsv"
    required = (sample_path, candidates_path, extinction_validation_path,
                catalogue_status_path, control_path, annulus_summary_path,
                annulus_path, association_summary_path, association_path)
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        parser.error("Missing inputs: " + ", ".join(missing))

    sample = read_json(sample_path)
    sample_rows = sample.get("clusters", [])
    control_manifest = read_json(control_path)
    controls = control_manifest.get("clusters", [])
    ext_validation = read_json(extinction_validation_path)
    catalogue_status = read_json(catalogue_status_path)
    ann_summary = read_json(annulus_summary_path)
    assoc_summary = read_json(association_summary_path)
    if len(sample_rows) != 800 or len(controls) != 798:
        parser.error("Unexpected target/control manifest size")
    if ext_validation.get("status") != "pass" or ext_validation.get("science_measurement") is not False:
        parser.error("Expected passed exploratory SFD-map/WCS validation")
    if ann_summary.get("status") != "complete" or ann_summary.get("science_measurement") is not False:
        parser.error("Unexpected annulus sensitivity status/science gate")
    if assoc_summary.get("status") != "complete" or assoc_summary.get("science_measurement") is not False:
        parser.error("Unexpected satellite association status/science gate")
    if ann_summary.get("output_sha256", {}).get(str(annulus_path)) != sha256(annulus_path):
        parser.error("Annulus profile table hash mismatch")
    if assoc_summary.get("output_sha256", {}).get(str(association_path)) != sha256(association_path):
        parser.error("Satellite association table hash mismatch")

    candidate_rows = read_tsv(candidates_path)
    if len(candidate_rows) != 800:
        parser.error("Expected one validated SFD extinction row per target")
    candidates = {str(row["cluster_id"]): row for row in candidate_rows}
    if len(candidates) != 800:
        parser.error("Duplicate target IDs in SFD extinction table")
    coeffs = {}
    for band in ("i", "i2"):
        values = [float(row["coefficient"]) for row in candidate_rows if row["header_filter"] == band]
        if not values or max(values) - min(values) > 1e-10:
            parser.error(f"No unique extinction coefficient for {band}")
        coeffs[band] = float(values[0])

    forced_rows, forced_hashes = read_forced_extinction(catalogue_dir, catalogue_status)
    pdr3_ai, pdr3_audit = choose_nearest_forced(sample_rows, forced_rows)
    control_by_id = {str(row["id"]): row for row in controls}
    if len(control_by_id) != 798:
        parser.error("Duplicate control IDs in control manifest")

    annulus_all = read_tsv(annulus_path)
    annulus = [row for row in annulus_all
               if float(row["kron_factor"]) == 2.5
               and float(row["coverage_floor"]) == 0.40
               and row["radial_bin_kpc"] == "200-300"]
    assoc_all = read_tsv(association_path)
    assoc = [row for row in assoc_all
             if float(row["kron_factor"]) == 2.5
             and float(row["coverage_floor"]) == 0.40]
    annulus_by_key = {(str(row["id"]), row["filter"]): row for row in annulus}
    if len(annulus_by_key) != len(annulus) or len(annulus) != 224:
        parser.error("Unexpected or duplicate fixed-cut 200-300 kpc annulus sample")
    assoc_by_key = {(str(row["cluster_id"]), row["filter"]): row for row in assoc}
    if len(assoc_by_key) != len(assoc) or len(assoc) != 74:
        parser.error("Unexpected or duplicate fixed-cut satellite-association sample")
    if not set(assoc_by_key).issubset(annulus_by_key):
        parser.error("Satellite association rows are not a subset of the selected ICL annulus")

    paired = []
    annulus_corrected = []
    for key, row in annulus_by_key.items():
        cluster_id, band = key
        if band not in coeffs:
            raise ValueError(f"Unexpected target filter label: {band}")
        candidate = candidates[cluster_id]
        if candidate["header_filter"] != band:
            raise ValueError(f"Target filter metadata mismatch for {cluster_id}: {band}")
        control_id = str(row["control_id"])
        control = control_by_id.get(control_id)
        if control is None or str(control.get("target_cluster_id")) != cluster_id:
            raise ValueError(f"Matched-control provenance mismatch for {cluster_id}")
        ebv_target = float(candidate["ebv_sfd"])
        ebv_control = float(control["ebv_sfd"])
        a_target_sfd = ebv_target * coeffs[band]
        a_control_sfd = ebv_control * coeffs[band]
        a_target_pdr3 = float(pdr3_ai[cluster_id])
        f_target_sfd = 10 ** (0.4 * a_target_sfd)
        f_control_sfd = 10 ** (0.4 * a_control_sfd)
        f_target_pdr3 = 10 ** (0.4 * a_target_pdr3)
        t, c = float(row["target_flux"]), float(row["control_flux"])
        delta_raw = float(row["delta_flux"])
        if not math.isclose(delta_raw, t - c, rel_tol=0, abs_tol=1e-8):
            raise ValueError(f"Raw ICL target-control subtraction failed for {key}")
        delta_sfd = t * f_target_sfd - c * f_control_sfd
        delta_pdr3 = t * f_target_pdr3 - c * f_control_sfd
        common = {
            "cluster_id": cluster_id, "filter": band, "control_id": control_id,
            "radial_bin_kpc": "200-300", "target_flux_raw_njy_arcsec2": t,
            "control_flux_raw_njy_arcsec2": c, "icl_delta_raw_njy_arcsec2": delta_raw,
            "target_ebv_sfd": ebv_target, "control_ebv_sfd": ebv_control,
            "extinction_coefficient": coeffs[band], "target_a_sfd_mag": a_target_sfd,
            "target_a_pdr3_forced_i_mag": a_target_pdr3,
            "control_a_sfd_mag": a_control_sfd,
            "pdr3_nearest_separation_arcsec": pdr3_audit[cluster_id]["separation_arcsec"],
            "icl_delta_sfd_corrected_njy_arcsec2": delta_sfd,
            "icl_delta_pdr3_target_corrected_njy_arcsec2": delta_pdr3,
            "target_coverage": row["target_coverage"],
            "control_coverage": row["control_coverage"],
            "common_coverage": row["common_coverage"],
        }
        annulus_corrected.append(common)
        if key not in assoc_by_key:
            continue
        sr = assoc_by_key[key]
        if not math.isclose(float(sr["icl_delta_200_300"]), delta_raw, rel_tol=0, abs_tol=1e-6):
            raise ValueError(f"Satellite table ICL value does not reproduce annulus input for {key}")
        satellite_raw = float(sr["satellite_delta_50_75"])
        if not all(math.isfinite(x) for x in (satellite_raw, float(sr["z"]), float(sr["richness"]))):
            raise ValueError(f"Nonfinite association data for {key}")
        pair_record = dict(common)
        pair_record.update({
            "satellite_delta_raw_njy_arcsec2": satellite_raw,
            "satellite_delta_sfd_corrected_njy_arcsec2": satellite_raw * f_target_sfd,
            "satellite_delta_pdr3_target_corrected_njy_arcsec2": satellite_raw * f_target_pdr3,
            "z": float(sr["z"]), "richness": float(sr["richness"]),
            "target_ra_deg": float(sr["target_ra_deg"]), "target_dec_deg": float(sr["target_dec_deg"]),
            "control_ra_deg": float(sr["control_ra_deg"]), "control_dec_deg": float(sr["control_dec_deg"]),
        })
        paired.append(pair_record)

    work.mkdir(parents=True, exist_ok=False)
    ann_fields = list(annulus_corrected[0])
    pair_fields = list(paired[0])
    ann_out = work / "annulus_extinction_corrected.tsv"
    pair_out = work / "icl_satellite_extinction_corrected.tsv"
    write_tsv(ann_out, annulus_corrected, ann_fields)
    write_tsv(pair_out, paired, pair_fields)

    annulus_metrics = {}
    association_metrics = {}
    for index, band in enumerate(("i", "i2")):
        selected_ann = [row for row in annulus_corrected if row["filter"] == band]
        targets = [{"ra_deg": float(row["target_ra_deg"]), "dec_deg": float(row["target_dec_deg"])}
                   for row in annulus_by_key.values() if row["filter"] == band]
        controls_coords = [{"ra_deg": float(row["control_ra_deg"]), "dec_deg": float(row["control_dec_deg"])}
                           for row in annulus_by_key.values() if row["filter"] == band]
        _, ann_labels, _, _, n_ann_blocks, largest_ann_block = pair_block_labels(
            targets, controls_coords, np.ones(len(targets), dtype=bool))
        raw_values = np.array([row["icl_delta_raw_njy_arcsec2"] for row in selected_ann])
        sfd_values = np.array([row["icl_delta_sfd_corrected_njy_arcsec2"] for row in selected_ann])
        pdr3_values = np.array([row["icl_delta_pdr3_target_corrected_njy_arcsec2"] for row in selected_ann])
        annulus_metrics[band] = {
            "n_pairs": len(selected_ann), "n_joint_5deg_blocks": n_ann_blocks,
            "largest_joint_block": largest_ann_block,
            "mean_raw_delta_njy_arcsec2": float(np.mean(raw_values)),
            "mean_sfd_corrected_delta_njy_arcsec2": float(np.mean(sfd_values)),
            "mean_sfd_correction_shift_njy_arcsec2": float(np.mean(sfd_values - raw_values)),
            "mean_pdr3_target_corrected_delta_njy_arcsec2": float(np.mean(pdr3_values)),
            "mean_pdr3_target_correction_shift_njy_arcsec2": float(np.mean(pdr3_values - raw_values)),
            "mean_raw_delta_spatial_block_bootstrap_ci95": block_bootstrap_mean(
                raw_values, ann_labels, args.replicas, args.seed + 10 * index),
            "mean_sfd_corrected_delta_spatial_block_bootstrap_ci95": block_bootstrap_mean(
                sfd_values, ann_labels, args.replicas, args.seed + 10 * index + 1),
            "mean_pdr3_target_corrected_delta_spatial_block_bootstrap_ci95": block_bootstrap_mean(
                pdr3_values, ann_labels, args.replicas, args.seed + 10 * index + 2),
        }

        selected_assoc = [row for row in paired if row["filter"] == band]
        target_coords = [{"ra_deg": row["target_ra_deg"], "dec_deg": row["target_dec_deg"]}
                         for row in selected_assoc]
        control_assoc_coords = [{"ra_deg": row["control_ra_deg"], "dec_deg": row["control_dec_deg"]}
                                for row in selected_assoc]
        _, labels, _, _, n_blocks, largest = pair_block_labels(
            target_coords, control_assoc_coords, np.ones(len(target_coords), dtype=bool))
        metrics, _ = block_bootstrap_correlations(
            selected_assoc, labels, args.replicas, args.seed + 100 + index)
        association_metrics[band] = {"n_pairs": len(selected_assoc),
                                     "n_joint_5deg_blocks": n_blocks,
                                     "largest_joint_block": largest,
                                     "correlations": metrics}

    # Summarize how closely the nearest forced PDR3 catalogue a_i agrees with
    # the per-target SFD estimate; this is a cross-check, not the correction map.
    comparison = []
    for row in sample_rows:
        cluster_id = str(row["id"])
        candidate = candidates[cluster_id]
        band = candidate["header_filter"]
        a_sfd = float(candidate["ebv_sfd"]) * coeffs[band]
        a_cat = pdr3_ai[cluster_id]
        comparison.append((band, a_cat - a_sfd))
    pdr3_crosscheck = {}
    for band in ("i", "i2"):
        delta = np.array([value for label, value in comparison if label == band], dtype=float)
        pdr3_crosscheck[band] = {
            "n_targets": int(len(delta)),
            "median_signed_pdr3_minus_sfd_a_mag": float(np.median(delta)),
            "median_absolute_difference_mag": float(np.median(np.abs(delta))),
            "p95_absolute_difference_mag": float(np.percentile(np.abs(delta), 95)),
            "max_absolute_difference_mag": float(np.max(np.abs(delta))),
        }

    source_hashes = {str(path): sha256(path) for path in required}
    source_hashes.update(forced_hashes)
    summary = {
        "schema": "hsc-icl-satellite-center-extinction-sensitivity-v1",
        "status": "complete_exploratory", "science_measurement": False,
        "method": {
            "primary": "apply per-band SFD center extinction A=coefficient*E(B-V) separately to target and matched-control annulus fluxes",
            "sensitivity": "replace target-center SFD A with nearest HSC PDR3 forced-source a_i; retain SFD for matched control",
            "satellite_profile": "multiply target-image satellite delta by the same target-center factor",
            "association": "partial Pearson residualized on z and richness; paired spatial block bootstrap over connected joint 5-degree target/control blocks",
            "bootstrap_replicas": args.replicas, "seed": args.seed,
            "sfd_coefficients_by_filter": coeffs,
        },
        "pdr3_nearest_forced_source_a_i_crosscheck": pdr3_crosscheck,
        "annulus_200_300_kpc": annulus_metrics,
        "satellite_association_icl200_300_vs_satellite50_75": association_metrics,
        "limitations": [
            "Extinction factors are evaluated at target/control field centers, not per pixel or at each profile location; the correction does not account for a spatially varying dust map across the annulus.",
            "PDR3 forced-source a_i is an i-band catalogue correction and is a sensitivity check for the target image; it is not a pixel-level correction and does not independently calibrate every profile pixel.",
            "The input annulus and satellite samples were selected before extinction correction; no selection correction or multiplicity-adjusted confirmatory inference is performed.",
            "Block bootstrap reflects the chosen connected 5-degree blocks only; it does not include sky subtraction, PSF wings, mask leakage, injection/recovery, or map-systematic uncertainties.",
            "Existing ICL profiles remain exploratory; no final ICL detection or physical satellite–ICL association is claimed.",
        ],
        "source_sha256": source_hashes,
        "program_sha256": sha256(__file__),
        "output_sha256": {str(ann_out): sha256(ann_out), str(pair_out): sha256(pair_out)},
    }
    summary_path = work / "summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False,
                                       allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({"status": summary["status"], "science_measurement": False,
                      "annulus_counts": {band: annulus_metrics[band]["n_pairs"] for band in annulus_metrics},
                      "association_counts": {band: association_metrics[band]["n_pairs"] for band in association_metrics},
                      "output": str(work)}, ensure_ascii=False))


if __name__ == "__main__":
    main()
