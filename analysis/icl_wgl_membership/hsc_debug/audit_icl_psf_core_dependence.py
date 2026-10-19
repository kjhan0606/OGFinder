#!/usr/bin/env python3
"""Audit whether measured ICL excess covaries with local model-PSF core size.

This is a diagnostic on the fixed 57-frame joint-injection support sample. It
does not infer PSF-wing corrections or promote the ICL signal to a detection.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

import numpy as np
from astropy.io import fits

from evaluate_prestack_pair_support import pair_block_labels
from summarize_hsc_psf_core_variation import measure_core, sha256

ROOT = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
DEFAULT_ICL = ROOT / "icl_annulus_support_2026-10-18_01"
DEFAULT_EXTRACTION = ROOT / "prestack_extraction_2026-10-15_01/extraction_summary.json"
DEFAULT_JOINT = ROOT / "joint_icl_satellite_psf_redetection_2026-10-19_02"
DEFAULT_PSF = ROOT / "psf_core_variation_2026-10-20_01/summary.json"
EXPECTED_PSF_SUMMARY_SHA256 = "2a1a0f63f7cbabadaf6f1a99f27d53ea9b87d8df2232305c73a89e8ce75807ad"
MODEL_COVARIATES = {
    "unadjusted": (),
    "redshift_richness_common_coverage": ("redshift", "richness", "common_coverage"),
    "redshift_richness_target_control_coverage": (
        "redshift", "richness", "target_coverage", "control_coverage"),
}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with Path(path).open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def partial_correlation(rows: list[dict], covariates: tuple[str, ...]) -> dict:
    if len(rows) < 5:
        raise ValueError("partial correlation requires at least five rows")
    x = np.asarray([float(row["icl_delta"]) for row in rows])
    y = np.asarray([float(row["psf_fwhm_major_arcsec"]) for row in rows])
    other = [np.asarray([float(row[field]) for row in rows]) for field in covariates]
    if not all(np.isfinite(values).all() for values in (x, y, *other)):
        raise ValueError("correlation inputs must be finite")
    design = np.column_stack((np.ones(len(rows)), *other))
    rank = int(np.linalg.matrix_rank(design))
    if rank != design.shape[1]:
        raise ValueError("nuisance design is rank deficient")
    residual_dof = len(rows) - rank
    if residual_dof < 3:
        raise ValueError("fewer than three residual degrees of freedom remain")
    x_resid = x - design @ np.linalg.lstsq(design, x, rcond=None)[0]
    y_resid = y - design @ np.linalg.lstsq(design, y, rcond=None)[0]
    denominator = float(np.linalg.norm(x_resid) * np.linalg.norm(y_resid))
    if not math.isfinite(denominator) or denominator <= 0:
        raise ValueError("partial correlation is undefined for a constant residual")
    return {"partial_r": float(x_resid @ y_resid / denominator),
            "n_rows": len(rows), "design_rank": rank,
            "residual_degrees_of_freedom": residual_dof}


def connected_block_labels(rows: list[dict]) -> np.ndarray:
    targets = [{"ra_deg": float(row["target_ra_deg"]),
                "dec_deg": float(row["target_dec_deg"])} for row in rows]
    controls = [{"ra_deg": float(row["control_ra_deg"]),
                 "dec_deg": float(row["control_dec_deg"])} for row in rows]
    selected, labels, *_ = pair_block_labels(targets, controls,
                                               np.ones(len(rows), dtype=bool))
    if len(selected) != len(rows) or len(labels) != len(rows):
        raise ValueError("connected-block labeling changed the sample")
    return labels


def block_bootstrap(rows: list[dict], labels: np.ndarray, covariates: tuple[str, ...],
                    replicas: int, seed: int) -> dict:
    labels = np.asarray(labels, dtype=str)
    if len(labels) != len(rows):
        raise ValueError("one connected spatial-block label is required per row")
    blocks = np.unique(labels)
    if len(blocks) < 6:
        return {"status": "insufficient_spatial_blocks", "n_blocks": len(blocks),
                "replicas_requested": replicas, "replicas_valid": 0,
                "percentiles_2_5_50_97_5": None}
    indices_by_block = {block: np.flatnonzero(labels == block) for block in blocks}
    rng = np.random.default_rng(seed)
    draws = []
    invalid = 0
    for _ in range(replicas):
        chosen = rng.choice(blocks, size=len(blocks), replace=True)
        indices = np.concatenate([indices_by_block[block] for block in chosen])
        sample = [rows[int(index)] for index in indices]
        try:
            draws.append(partial_correlation(sample, covariates)["partial_r"])
        except ValueError:
            invalid += 1
    if not draws:
        return {"status": "no_valid_replicas", "n_blocks": len(blocks),
                "replicas_requested": replicas, "replicas_valid": 0,
                "replicas_invalid": invalid, "percentiles_2_5_50_97_5": None}
    return {"status": "complete_descriptive_bootstrap", "n_blocks": len(blocks),
            "largest_block_rows": max(len(indices) for indices in indices_by_block.values()),
            "replicas_requested": replicas, "replicas_valid": len(draws),
            "replicas_invalid": invalid, "seed": seed,
            "percentiles_2_5_50_97_5": np.percentile(draws, [2.5, 50, 97.5]).tolist()}


def leave_one_block_out(rows: list[dict], labels: np.ndarray,
                        covariates: tuple[str, ...]) -> dict:
    values = []
    blocks = np.unique(labels)
    for block in blocks:
        sample = [row for row, label in zip(rows, labels) if label != block]
        try:
            result = partial_correlation(sample, covariates)
        except ValueError:
            continue
        values.append({"omitted_block": str(block), "n_rows": len(sample),
                       "partial_r": result["partial_r"]})
    if not values:
        return {"n_blocks": int(len(blocks)), "n_valid_omissions": 0,
                "partial_r_range": None, "records": []}
    estimates = [row["partial_r"] for row in values]
    return {"n_blocks": int(len(blocks)), "n_valid_omissions": len(values),
            "partial_r_range": [float(min(estimates)), float(max(estimates))],
            "records": values}


def load_verified_inputs(icl_dir: Path, extraction_path: Path,
                         joint_dir: Path, psf_summary_path: Path) -> tuple[list[dict], dict]:
    icl_dir, extraction_path = Path(icl_dir).resolve(), Path(extraction_path).resolve()
    joint_dir, psf_summary_path = Path(joint_dir).resolve(), Path(psf_summary_path).resolve()
    icl_summary_path = icl_dir / "summary.json"
    annulus_path = icl_dir / "annulus_pair_profiles.tsv"
    joint_summary_path = joint_dir / "summary.json"
    joint_audit_path = joint_dir / "psf_redetection_audit.json"
    for path in (icl_summary_path, annulus_path, extraction_path, joint_summary_path,
                 joint_audit_path, psf_summary_path):
        if not path.is_file():
            raise FileNotFoundError(path)

    icl_summary = json.loads(icl_summary_path.read_text(encoding="utf-8"))
    if (icl_summary.get("schema") != "hsc-icl-per-annulus-support-v1"
            or icl_summary.get("status") != "complete"
            or icl_summary.get("science_gates", {}).get("science_measurement") is not False
            or icl_summary.get("output_sha256", {}).get(str(annulus_path)) != sha256(annulus_path)):
        raise ValueError("ICL annulus table is incomplete, changed, or science-promoted")
    if icl_summary.get("input_sha256", {}).get(str(extraction_path)) != sha256(extraction_path):
        raise ValueError("matched-control extraction summary hash differs from the ICL run record")

    extraction = json.loads(extraction_path.read_text(encoding="utf-8"))
    if (extraction.get("status") != "complete"
            or str(extraction.get("science_measurement", "")).lower() != "false"):
        raise ValueError("matched-control coordinate extraction is incomplete or science-promoted")
    controls = {}
    for row in extraction.get("records", []):
        ident = str(row["control_id"])
        if ident in controls or row.get("frame_role") != "matched_control":
            raise ValueError(f"duplicate or non-control extraction record: {ident}")
        controls[ident] = row

    joint_summary = json.loads(joint_summary_path.read_text(encoding="utf-8"))
    joint_audit = json.loads(joint_audit_path.read_text(encoding="utf-8"))
    if (joint_summary.get("schema") != "hsc-joint-icl-satellite-psf-redetection-batch-v1"
            or joint_summary.get("status") != "complete"
            or joint_summary.get("n_completed") != 57
            or joint_summary.get("n_errors") != 0
            or joint_summary.get("science_measurement") is not False
            or joint_audit.get("status") != "verified"
            or joint_audit.get("n_frame_results_checked") != 57
            or joint_audit.get("batch_summary_sha256") != sha256(joint_summary_path)):
        raise ValueError("joint PSF-redetection run is not fully hash-audited")

    psf_digest = sha256(psf_summary_path)
    if psf_digest != EXPECTED_PSF_SUMMARY_SHA256:
        raise ValueError("PSF core summary differs from its audited checksum")
    psf_summary = json.loads(psf_summary_path.read_text(encoding="utf-8"))
    if (psf_summary.get("schema") != "astrafek-hsc-position-psf-core-variation-v1"
            or psf_summary.get("status") != "complete_exploratory_psf_core_audit"
            or psf_summary.get("n_total_models_including_test") != 745
            or psf_summary.get("science_measurement") is not False):
        raise ValueError("PSF core summary is incomplete or science-promoted")

    centers = {}
    psf_download_dir = Path(psf_summary["inputs"]["batch_audit"]).resolve().parent
    test_psf = Path(psf_summary["inputs"]["test_psf"]).resolve()
    for psf_row in psf_summary["records"]:
        for role in psf_row["roles"]:
            if role["role"] != "cluster_center":
                continue
            key = (str(role["cluster_id"]), str(role["frame_filter"]))
            if key in centers:
                raise ValueError(f"duplicate center PSF for frame: {key}")
            path = (test_psf if psf_row["source"] == "approved_single_position_test"
                    else psf_download_dir / psf_row["file"])
            digest = sha256(path)
            if digest != psf_row["sha256"]:
                raise ValueError(f"center PSF file checksum changed: {path.name}")
            with fits.open(path, memmap=False, checksum=True) as hdus:
                if len(hdus) != 1 or hdus[0].data is None:
                    raise ValueError(f"unexpected center PSF FITS structure: {path}")
                metric = measure_core(hdus[0].data)
            if not math.isclose(metric["moment_fwhm_major_arcsec"],
                                float(psf_row["moment_fwhm_major_arcsec"]),
                                rel_tol=0, abs_tol=1e-12):
                raise ValueError(f"center PSF metric does not reproduce: {path.name}")
            centers[key] = {"path": str(path), "sha256": digest,
                            "moment_fwhm_major_arcsec": metric["moment_fwhm_major_arcsec"]}
    if len(centers) != 57:
        raise ValueError(f"expected 57 center PSFs, found {len(centers)}")

    annulus_rows = {}
    for row in read_tsv(annulus_path):
        if row["radial_bin_kpc"] != "200-300":
            continue
        key = (str(row["id"]), row["filter"])
        if key in annulus_rows:
            raise ValueError(f"duplicate ICL target-control annulus row: {key}")
        annulus_rows[key] = row

    output, seen = [], set()
    for record in joint_summary["records"]:
        key = (str(record["cluster_id"]), record["filter"])
        if key in seen or key not in annulus_rows or key not in centers:
            raise ValueError(f"duplicate or unmatched joint-injection frame: {key}")
        seen.add(key)
        annulus = annulus_rows[key]
        redetected_coverage = float(record["icl_coverage_200_300"])
        if not math.isfinite(redetected_coverage) or not 0 <= redetected_coverage <= 1:
            raise ValueError(f"invalid redetected ICL coverage for {key}")
        frame_path = Path(record["path"])
        if not frame_path.is_file() or sha256(frame_path) != record["sha256"]:
            raise ValueError(f"joint-injection frame result hash mismatch: {key}")
        control = controls.get(str(annulus["control_id"]))
        if (control is None or control.get("header_filter") != record["filter"]
                or str(control.get("target_cluster_id")) != key[0]):
            raise ValueError(f"matched-control extraction is missing or filter-mismatched: {key}")
        control_ra, control_dec = float(control["ra_deg"]), float(control["dec_deg"])
        if not (math.isfinite(control_ra) and 0 <= control_ra < 360
                and math.isfinite(control_dec) and -90 <= control_dec <= 90):
            raise ValueError(f"invalid matched-control coordinate for {key}")
        values = {
            "cluster_id": key[0], "filter": key[1],
            "icl_delta": float(annulus["delta_sb_njy_arcsec2"]),
            "psf_fwhm_major_arcsec": centers[key]["moment_fwhm_major_arcsec"],
            "redshift": float(annulus["z"]), "richness": float(annulus["richness"]),
            "target_coverage": float(annulus["target_coverage"]),
            "control_coverage": float(annulus["control_coverage"]),
            "common_coverage": float(annulus["common_coverage"]),
            "joint_redetected_icl_coverage": redetected_coverage,
            "redetected_minus_saved_coverage": redetected_coverage - float(annulus["common_coverage"]),
            "target_ra_deg": float(annulus["ra_deg"]),
            "target_dec_deg": float(annulus["dec_deg"]),
            "control_ra_deg": control_ra,
            "control_dec_deg": control_dec,
            "frame_result_sha256": record["sha256"],
            "center_psf_sha256": centers[key]["sha256"],
        }
        if not all(math.isfinite(float(values[name])) for name in (
                "icl_delta", "psf_fwhm_major_arcsec", "redshift", "richness",
                "target_coverage", "control_coverage", "common_coverage")):
            raise ValueError(f"non-finite joined frame data: {key}")
        output.append(values)
    if len(output) != 57 or set(row["filter"] for row in output) != {"i", "i2"}:
        raise ValueError(f"fixed PSF frame sample changed: {len(output)}")
    if len({row["cluster_id"] for row in output}) != 57:
        raise ValueError("fixed PSF frame sample contains repeated clusters across filters")
    provenance = {
        "icl_summary_sha256": sha256(icl_summary_path),
        "icl_annulus_table_sha256": sha256(annulus_path),
        "extraction_summary_sha256": sha256(extraction_path),
        "joint_summary_sha256": sha256(joint_summary_path),
        "joint_audit_sha256": sha256(joint_audit_path),
        "psf_core_summary_sha256": psf_digest,
        "n_joined_frames": len(output),
    }
    return output, provenance


def analyze(icl_dir: Path, extraction_path: Path, joint_dir: Path,
            psf_summary_path: Path, output_path: Path, replicas: int = 10000,
            seed: int = 20261020) -> dict:
    if replicas < 1000:
        raise ValueError("at least 1000 spatial-block bootstrap replicas are required")
    output_path = Path(output_path).resolve()
    if output_path.exists():
        raise FileExistsError(f"refusing to overwrite existing output: {output_path}")
    rows, provenance = load_verified_inputs(icl_dir, extraction_path, joint_dir, psf_summary_path)
    by_filter = {}
    for band in ("i", "i2"):
        selected = [row for row in rows if row["filter"] == band]
        expected_n = 13 if band == "i" else 44
        if len(selected) != expected_n:
            raise ValueError(f"fixed {band} PSF support changed: {len(selected)} != {expected_n}")
        labels = connected_block_labels(selected)
        model_results = {}
        for name, covariates in MODEL_COVARIATES.items():
            point = partial_correlation(selected, covariates)
            loo = leave_one_block_out(selected, labels, covariates)
            boot = block_bootstrap(selected, labels, covariates, replicas,
                                   seed + (1 if band == "i" else 0) + list(MODEL_COVARIATES).index(name) * 101)
            model_results[name] = {**point, "covariates": list(covariates),
                                   "leave_one_block_out": loo,
                                   "spatial_block_bootstrap": boot}
        by_filter[band] = {"n_frames": len(selected),
                           "n_connected_target_control_blocks": int(len(np.unique(labels))),
                           "psf_fwhm_major_distribution_arcsec": {
                               "median": float(np.median([row["psf_fwhm_major_arcsec"] for row in selected])),
                               "q16": float(np.quantile([row["psf_fwhm_major_arcsec"] for row in selected], 0.16)),
                               "q84": float(np.quantile([row["psf_fwhm_major_arcsec"] for row in selected], 0.84)),
                           },
                           "models": model_results}
    result = {
        "schema": "astrafek-icl-psf-core-dependence-v1",
        "status": "complete_exploratory_sensitivity",
        "science_measurement": False,
        "psf_wings_validated": False,
        "sample": "fixed 57-frame PSF-redetection support joined one-to-one to hash-verified target-minus-control ICL annuli and center PSFs",
        "method": {
            "association": "partial Pearson correlation of ICL 200-300 kpc target-minus-control brightness and cluster-center model-PSF moment FWHM",
            "adjustments": "redshift, richness, and ICL target/control area support",
            "spatial_resampling": "paired connected 5-degree target/control sky-cell bootstrap; percentile ranges are descriptive, not calibrated confidence intervals",
            "interpretation": "correlation with compact model-PSF core only; no causal attribution or correction",
        },
        "by_filter": by_filter,
        "records": rows,
        "provenance": provenance,
        "limitations": [
            "57 frames are a support-selected subset (13 i and 44 i2), not the full 800-cluster sample",
            "model PSF stamps cover only about 7 arcsec and do not measure extended wings on ICL scales",
            "spatial-block intervals are unstable with only 10 i blocks and 24 i2 blocks",
            "regression adjustment cannot remove selection or unmeasured sky, mask, noise, and field effects",
            "the joint-run redetected mask can yield a different area support than the saved ICL annulus mask; both are retained separately",
            "this association audit does not replace matched-random lens subtraction or observational transfer validation",
            "no new HSC request or download was made",
        ],
    }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("x", encoding="utf-8") as stream:
        json.dump(result, stream, indent=2, ensure_ascii=False, allow_nan=False)
        stream.write("\n")
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--icl", type=Path, default=DEFAULT_ICL)
    parser.add_argument("--extraction-summary", type=Path, default=DEFAULT_EXTRACTION)
    parser.add_argument("--joint", type=Path, default=DEFAULT_JOINT)
    parser.add_argument("--psf-summary", type=Path, default=DEFAULT_PSF)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--replicas", type=int, default=10000)
    parser.add_argument("--seed", type=int, default=20261020)
    args = parser.parse_args()
    result = analyze(args.icl, args.extraction_summary, args.joint,
                     args.psf_summary, args.output, args.replicas, args.seed)
    print(json.dumps({"status": result["status"], "science_measurement": result["science_measurement"],
                      "by_filter": {band: {"n_frames": v["n_frames"],
                                           "n_blocks": v["n_connected_target_control_blocks"],
                                           "models": {name: {"partial_r": model["partial_r"],
                                                             "ci": model["spatial_block_bootstrap"]["percentiles_2_5_50_97_5"],
                                                             "loo_range": model["leave_one_block_out"]["partial_r_range"]}
                                                       for name, model in v["models"].items()}}
                                    for band, v in result["by_filter"].items()}}, indent=2))


if __name__ == "__main__":
    main()
