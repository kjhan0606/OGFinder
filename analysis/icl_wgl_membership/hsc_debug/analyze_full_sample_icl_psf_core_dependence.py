#!/usr/bin/env python3
"""Test whether paired PSF-core mismatch tracks the existing ICL excess.

The sample is restricted to the existing common-area-qualified 200--300 kpc
target/control profile rows. The PSF mismatch is a descriptive moment-width
difference, not a wing correction or an image-level homogenization.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import platform

import numpy as np

from audit_icl_psf_core_dependence import connected_block_labels


ROOT = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
DEFAULT_PROFILES = ROOT / "prestack_flux_profiles_2026-10-18_02"
DEFAULT_PSF_DIAGNOSTIC = ROOT / "target_control_psf_core_diagnostic_2026-10-21_01"
MODEL_COVARIATES = {
    "unadjusted": (),
    "redshift_richness_common_coverage": ("redshift", "richness", "common_coverage"),
    "plus_center_mask_flag": ("redshift", "richness", "common_coverage", "center_masked"),
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with Path(path).open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def partial_correlation(rows: list[dict], covariates: tuple[str, ...]) -> dict:
    if len(rows) < 5:
        raise ValueError("partial correlation requires at least five rows")
    x = np.asarray([float(row["icl_delta_njy_arcsec2"]) for row in rows])
    y = np.asarray([float(row["psf_core_delta_fwhm_major_arcsec"]) for row in rows])
    other = [np.asarray([float(row[name]) for row in rows]) for name in covariates]
    if not all(np.isfinite(value).all() for value in (x, y, *other)):
        raise ValueError("correlation inputs must be finite")
    design = np.column_stack((np.ones(len(rows)), *other))
    rank = int(np.linalg.matrix_rank(design))
    if rank != design.shape[1] or len(rows) - rank < 3:
        raise ValueError("nuisance design is rank deficient or leaves too few degrees of freedom")
    x_resid = x - design @ np.linalg.lstsq(design, x, rcond=None)[0]
    y_resid = y - design @ np.linalg.lstsq(design, y, rcond=None)[0]
    denom = float(np.linalg.norm(x_resid) * np.linalg.norm(y_resid))
    if not np.isfinite(denom) or denom <= 0:
        raise ValueError("partial correlation is undefined for a constant residual")
    return {"partial_r": float(x_resid @ y_resid / denom),
            "n_rows": len(rows), "design_rank": rank,
            "residual_degrees_of_freedom": len(rows) - rank}


def block_bootstrap(rows: list[dict], labels: np.ndarray,
                    covariates: tuple[str, ...], replicas: int, seed: int) -> dict:
    labels = np.asarray(labels, dtype=str)
    if len(labels) != len(rows):
        raise ValueError("one connected block label is required per row")
    blocks = np.unique(labels)
    if len(blocks) < 6:
        return {"status": "insufficient_spatial_blocks", "n_blocks": len(blocks),
                "replicas_requested": replicas, "replicas_valid": 0,
                "replicas_invalid": 0, "percentiles_2_5_50_97_5": None}
    by_block = {block: np.flatnonzero(labels == block) for block in blocks}
    rng = np.random.default_rng(seed)
    draws, invalid = [], 0
    for _ in range(replicas):
        selected = rng.choice(blocks, size=len(blocks), replace=True)
        indices = np.concatenate([by_block[key] for key in selected])
        try:
            draws.append(partial_correlation([rows[int(i)] for i in indices], covariates)["partial_r"])
        except ValueError:
            invalid += 1
    if not draws:
        return {"status": "no_valid_replicas", "n_blocks": len(blocks),
                "replicas_requested": replicas, "replicas_valid": 0,
                "replicas_invalid": invalid, "percentiles_2_5_50_97_5": None}
    return {"status": "complete_descriptive_block_bootstrap", "n_blocks": len(blocks),
            "largest_block_rows": max(map(len, by_block.values())),
            "replicas_requested": replicas, "replicas_valid": len(draws),
            "replicas_invalid": invalid, "seed": seed,
            "percentiles_2_5_50_97_5": np.percentile(draws, [2.5, 50, 97.5]).tolist()}


def leave_one_block_out(rows: list[dict], labels: np.ndarray,
                        covariates: tuple[str, ...]) -> dict:
    blocks = np.unique(labels)
    values = []
    for block in blocks:
        selected = [row for row, label in zip(rows, labels) if label != block]
        try:
            values.append(partial_correlation(selected, covariates)["partial_r"])
        except ValueError:
            continue
    return {"n_blocks": int(len(blocks)), "n_valid_omissions": len(values),
            "partial_r_range": ([float(min(values)), float(max(values))] if values else None)}


def join_rows(profile_rows: list[dict[str, str]], pair_rows: list[dict[str, str]],
              psf_pairs: list[dict[str, str]], targets: dict[str, dict],
              controls: dict[str, dict]) -> list[dict]:
    psf_by_pair = {}
    for row in psf_pairs:
        key = (str(row["target_cluster_id"]), str(row["control_id"]))
        if key in psf_by_pair:
            raise ValueError(f"duplicate PSF pair: {key}")
        psf_by_pair[key] = row
    profile_by_key = {}
    for row in profile_rows:
        key = (str(row["id"]), str(row["filter"]))
        if key in profile_by_key:
            raise ValueError(f"duplicate target profile: {key}")
        profile_by_key[key] = row
    seen = set()
    output = []
    for row in pair_rows:
        ident, control_id, band = str(row["id"]), str(row["control_id"]), str(row["filter"])
        key = (ident, control_id, band)
        if key in seen:
            raise ValueError(f"duplicate ICL target/control/filter row: {key}")
        seen.add(key)
        psf = psf_by_pair.get((ident, control_id))
        profile = profile_by_key.get((ident, band))
        control = controls.get(control_id)
        target = targets.get(ident)
        if psf is None or profile is None or control is None or target is None:
            raise ValueError(f"unmatched ICL/PSF/manifest row: {key}")
        if str(control.get("target_cluster_id")) != ident:
            raise ValueError(f"control-to-target manifest link differs for {key}")
        if str(target.get("id")) != ident:
            raise ValueError(f"target manifest ID differs for {key}")
        common = float(row["coverage_200_300"])
        values = {
            "cluster_id": ident, "control_id": control_id, "filter": band,
            "icl_delta_njy_arcsec2": float(row["delta_sb_200_300"]),
            "psf_core_delta_fwhm_major_arcsec": (
                float(psf["target_minus_control_moment_fwhm_major_arcsec"])),
            "psf_target_fwhm_major_arcsec": float(psf["target_moment_fwhm_major_arcsec"]),
            "psf_control_fwhm_major_arcsec": float(psf["control_moment_fwhm_major_arcsec"]),
            "redshift": float(row["z"]), "richness": float(row["richness"]),
            "common_coverage": common,
            "center_masked": int(str(row["center_masked_k2p5"]).lower() == "true"),
            "target_ra_deg": float(target["ra_deg"]),
            "target_dec_deg": float(target["dec_deg"]),
            "control_ra_deg": float(control["ra_deg"]),
            "control_dec_deg": float(control["dec_deg"]),
            "coverage_floor_met": common >= 0.40,
        }
        if not all(np.isfinite(float(values[name])) for name in (
                "icl_delta_njy_arcsec2", "psf_core_delta_fwhm_major_arcsec",
                "redshift", "richness", "common_coverage", "target_ra_deg",
                "target_dec_deg", "control_ra_deg", "control_dec_deg")):
            raise ValueError(f"non-finite joined row: {key}")
        if values["filter"] not in {"i", "i2"}:
            raise ValueError(f"unexpected diagnostic filter label: {band}")
        output.append(values)
    if len(output) != len(pair_rows):
        raise ValueError("not all ICL profile rows were joined")
    return output


def analyze(profile_dir: Path, psf_dir: Path, output_dir: Path,
            replicas: int = 5000, seed: int = 20261021) -> dict:
    profile_dir, psf_dir = Path(profile_dir).resolve(), Path(psf_dir).resolve()
    output_dir = Path(output_dir).resolve()
    if output_dir.exists():
        raise FileExistsError(f"refusing to overwrite {output_dir}")
    if replicas < 1000:
        raise ValueError("at least 1000 block-bootstrap replicas are required")
    profile_summary_path = profile_dir / "summary.json"
    profile_path = profile_dir / "pair_metrics_k2p5_c40_r50_300.tsv"
    target_profile_path = profile_dir / "cluster_profiles.tsv"
    psf_summary_path = psf_dir / "summary.json"
    psf_pairs_path = psf_dir / "paired_target_control_psf_core_metrics.tsv"
    for path in (profile_summary_path, profile_path, target_profile_path,
                 psf_summary_path, psf_pairs_path):
        if not path.is_file():
            raise FileNotFoundError(path)
    profile_summary = json.loads(profile_summary_path.read_text(encoding="utf-8"))
    if (profile_summary.get("schema") != "hsc-icl-full-sample-profiles-v1"
            or profile_summary.get("status") != "complete"
            or profile_summary.get("science_gates", {}).get("science_measurement") is not False
            or profile_summary.get("output_sha256", {}).get(str(profile_path)) != sha256(profile_path)
            or profile_summary.get("output_sha256", {}).get(str(target_profile_path)) != sha256(target_profile_path)):
        raise ValueError("full-sample profile inputs are incomplete, changed, or promoted")
    psf_summary = json.loads(psf_summary_path.read_text(encoding="utf-8"))
    if (psf_summary.get("schema") != "astrafek-hsc-target-control-psf-core-diagnostic-v1"
            or psf_summary.get("status") != "complete_descriptive_core_mismatch_diagnostic"
            or psf_summary.get("science_measurement") is not False
            or psf_summary.get("psf_approved") is not False
            or psf_summary.get("support", {}).get("paired_target_control_centers") != 798
            or psf_summary.get("tables", {}).get("paired_target_control_psf_core_metrics", {}).get("sha256") != sha256(psf_pairs_path)):
        raise ValueError("full-sample PSF-core diagnostic is incomplete or changed")
    target_plan_path = Path(psf_summary["inputs"]["target_control_plan"])
    target_plan = json.loads(target_plan_path.read_text(encoding="utf-8"))
    target_manifest = Path(target_plan["target_manifest"])
    control_manifest = Path(target_plan["control_manifest"])
    if (sha256(target_plan_path) != psf_summary["inputs"]["target_control_plan_sha256"]
            or sha256(target_manifest) != target_plan["target_manifest_sha256"]
            or sha256(control_manifest) != target_plan["control_manifest_sha256"]):
        raise ValueError("sample manifest provenance differs from the PSF plan")
    targets = {str(row["id"]): row for row in json.loads(target_manifest.read_text())["clusters"]}
    controls = {str(row["id"]): row for row in json.loads(control_manifest.read_text())["clusters"]}
    rows = join_rows(read_tsv(target_profile_path), read_tsv(profile_path),
                     read_tsv(psf_pairs_path), targets, controls)
    if len(rows) != 137:
        raise ValueError(f"fixed full-sample profile support changed: {len(rows)}")

    by_filter = {}
    for band, expected_n in (("i", 44), ("i2", 93)):
        selected = [row for row in rows if row["filter"] == band]
        if len(selected) != expected_n:
            raise ValueError(f"{band} support changed: {len(selected)} != {expected_n}")
        labels = connected_block_labels(selected)
        models = {}
        for index, (name, covariates) in enumerate(MODEL_COVARIATES.items()):
            point = partial_correlation(selected, covariates)
            models[name] = {
                **point, "covariates": list(covariates),
                "spatial_block_bootstrap": block_bootstrap(
                    selected, labels, covariates, replicas,
                    seed + (0 if band == "i2" else 10_000) + index * 101),
                "leave_one_block_out": leave_one_block_out(selected, labels, covariates),
            }
        by_filter[band] = {
            "n_pairs": len(selected),
            "n_connected_target_control_spatial_blocks": int(len(np.unique(labels))),
            "icl_delta_distribution_njy_arcsec2": {
                "median": float(np.median([r["icl_delta_njy_arcsec2"] for r in selected])),
                "q16": float(np.quantile([r["icl_delta_njy_arcsec2"] for r in selected], .16)),
                "q84": float(np.quantile([r["icl_delta_njy_arcsec2"] for r in selected], .84)),
            },
            "psf_core_mismatch_distribution_arcsec": {
                "median": float(np.median([r["psf_core_delta_fwhm_major_arcsec"] for r in selected])),
                "q16": float(np.quantile([r["psf_core_delta_fwhm_major_arcsec"] for r in selected], .16)),
                "q84": float(np.quantile([r["psf_core_delta_fwhm_major_arcsec"] for r in selected], .84)),
                "median_absolute_difference": float(np.median([abs(r["psf_core_delta_fwhm_major_arcsec"]) for r in selected])),
            },
            "models": models,
        }

    output_dir.mkdir(parents=True, exist_ok=False)
    joined_path = output_dir / "joined_icl_psf_core_rows.tsv"
    with joined_path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    script_path = Path(__file__).resolve()
    result = {
        "schema": "astrafek-full-sample-icl-psf-core-dependence-v1",
        "status": "complete_exploratory_psf_mismatch_sensitivity",
        "science_measurement": False,
        "psf_approved": False,
        "sample": "existing 200-300 kpc K2.5 matched target-control profile rows passing the 40% common-coverage floor; joined to the 800/798 target/control model-PSF center bank",
        "method": {
            "association": "partial Pearson correlation between target-minus-control ICL surface brightness and target-minus-control moment-based major-axis PSF-core FWHM",
            "covariates": "redshift, richness, common annular coverage; an additional sensitivity model includes the central K2.5 mask flag",
            "spatial_resampling": "connected 5-degree target/control sky-cell block bootstrap; percentile intervals are descriptive, not calibrated confidence intervals",
            "filter_labels": "i/i2 follow the existing science-frame strata; both use the HSC-I model coadd PSF product",
        },
        "support": {"n_joined_pairs": len(rows),
                    "filter_counts": {band: sum(r["filter"] == band for r in rows) for band in ("i", "i2")},
                    "center_masked_k2p5": sum(r["center_masked"] for r in rows),
                    "coverage_floor_minimum": min(r["common_coverage"] for r in rows),
                    "coverage_floor_maximum": max(r["common_coverage"] for r in rows)},
        "by_filter": by_filter,
        "provenance": {
            "profile_summary_sha256": sha256(profile_summary_path),
            "profile_pair_table_sha256": sha256(profile_path),
            "target_profile_table_sha256": sha256(target_profile_path),
            "psf_core_summary_sha256": sha256(psf_summary_path),
            "psf_core_pair_table_sha256": sha256(psf_pairs_path),
            "target_manifest_sha256": sha256(target_manifest),
            "control_manifest_sha256": sha256(control_manifest),
        },
        "joined_table": {"path": str(joined_path), "sha256": sha256(joined_path), "rows": len(rows)},
        "limitations": [
            "137 profile rows are a common-coverage-selected subset of the 800-cluster parent (44 i, 93 i2).",
            "Core moment width is a seeing/PSF-core proxy, not an extended-wing measurement.",
            "The HSC-I model PSF cores do not measure the local sky-subtraction transfer at 200-300 kpc.",
            "Block-bootstrap ranges are descriptive sensitivity intervals; no p-values or detection claim are produced.",
            "The diagnostic does not subtract a PSF correction from ICL data or approve science measurement.",
        ],
        "script": str(script_path), "script_sha256": sha256(script_path),
        "python": platform.python_version(), "numpy": np.__version__,
        "bootstrap_replicas": replicas, "bootstrap_seed": seed,
    }
    summary_path = output_dir / "summary.json"
    summary_path.write_text(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False) + "\n",
                            encoding="utf-8")
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile-dir", type=Path, default=DEFAULT_PROFILES)
    parser.add_argument("--psf-diagnostic-dir", type=Path, default=DEFAULT_PSF_DIAGNOSTIC)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--replicas", type=int, default=5000)
    parser.add_argument("--seed", type=int, default=20261021)
    args = parser.parse_args()
    result = analyze(args.profile_dir, args.psf_diagnostic_dir,
                     args.output_dir, args.replicas, args.seed)
    print(json.dumps({"status": result["status"], "support": result["support"],
                      "by_filter": result["by_filter"],
                      "science_measurement": result["science_measurement"],
                      "psf_approved": result["psf_approved"],
                      "summary_sha256": sha256(Path(args.output_dir) / "summary.json")}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
