#!/usr/bin/env python3
"""Independently reproduce the full-sample ICL–PSF-core sensitivity result."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np

from audit_icl_psf_core_dependence import connected_block_labels


COVARIATES = {
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


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with Path(path).open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def _independent_partial_r(rows: list[dict], covariates: tuple[str, ...]) -> float:
    x = np.array([float(row["icl_delta_njy_arcsec2"]) for row in rows])
    y = np.array([float(row["psf_core_delta_fwhm_major_arcsec"]) for row in rows])
    z = np.column_stack((np.ones(len(rows)), *[
        np.array([float(row[name]) for row in rows]) for name in covariates
    ]))
    # QR residualization is independent of the analysis script's least-squares implementation.
    q, _ = np.linalg.qr(z, mode="reduced")
    xr, yr = x - q @ (q.T @ x), y - q @ (q.T @ y)
    return float(np.dot(xr, yr) / (np.linalg.norm(xr) * np.linalg.norm(yr)))


def _recompute_bootstrap(rows: list[dict], labels: np.ndarray,
                         covariates: tuple[str, ...], replicas: int,
                         seed: int) -> tuple[list[float], list[float]]:
    blocks = np.unique(np.asarray(labels, dtype=str))
    indices = {block: np.flatnonzero(labels == block) for block in blocks}
    rng = np.random.default_rng(seed)
    estimates = []
    for _ in range(replicas):
        draw = rng.choice(blocks, size=len(blocks), replace=True)
        selected = np.concatenate([indices[block] for block in draw])
        estimates.append(_independent_partial_r([rows[int(i)] for i in selected], covariates))
    leave = []
    for block in blocks:
        selected = [row for row, label in zip(rows, labels) if label != block]
        leave.append(_independent_partial_r(selected, covariates))
    return np.percentile(estimates, [2.5, 50, 97.5]).tolist(), [min(leave), max(leave)]


def validate(profile_dir: Path, psf_dir: Path, result_dir: Path) -> dict:
    profile_dir, psf_dir, result_dir = map(Path, (profile_dir, psf_dir, result_dir))
    summary_path = result_dir / "summary.json"
    joined_path = result_dir / "joined_icl_psf_core_rows.tsv"
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    rows = _read_tsv(joined_path)
    if (summary.get("status") != "complete_exploratory_psf_mismatch_sensitivity"
            or summary.get("science_measurement") is not False
            or summary.get("psf_approved") is not False):
        raise ValueError("result status or science gates changed")
    if sha256(joined_path) != summary.get("joined_table", {}).get("sha256"):
        raise ValueError("joined table checksum mismatch")
    if len(rows) != 137 or len({(r["cluster_id"], r["control_id"], r["filter"]) for r in rows}) != 137:
        raise ValueError("joined row count or pair uniqueness changed")

    profile_table = profile_dir / "pair_metrics_k2p5_c40_r50_300.tsv"
    target_table = profile_dir / "cluster_profiles.tsv"
    profile_summary = profile_dir / "summary.json"
    psf_summary_path = psf_dir / "summary.json"
    psf_pair_table = psf_dir / "paired_target_control_psf_core_metrics.tsv"
    for path, key in ((profile_summary, "profile_summary_sha256"),
                      (profile_table, "profile_pair_table_sha256"),
                      (target_table, "target_profile_table_sha256"),
                      (psf_summary_path, "psf_core_summary_sha256"),
                      (psf_pair_table, "psf_core_pair_table_sha256")):
        if sha256(path) != summary["provenance"][key]:
            raise ValueError(f"input hash mismatch: {path}")
    psf_summary = json.loads(psf_summary_path.read_text(encoding="utf-8"))
    plan_path = Path(psf_summary["inputs"]["target_control_plan"])
    plan = json.loads(plan_path.read_text(encoding="utf-8"))
    target_manifest, control_manifest = Path(plan["target_manifest"]), Path(plan["control_manifest"])
    if (sha256(target_manifest) != summary["provenance"]["target_manifest_sha256"]
            or sha256(control_manifest) != summary["provenance"]["control_manifest_sha256"]):
        raise ValueError("sample manifest hash mismatch")

    icl = {(r["id"], r["control_id"], r["filter"]): r for r in _read_tsv(profile_table)}
    psf = {(r["target_cluster_id"], r["control_id"]): r for r in _read_tsv(psf_pair_table)}
    target_profiles = {(r["id"], r["filter"]): r for r in _read_tsv(target_table)}
    targets = {str(r["id"]): r for r in json.loads(target_manifest.read_text())["clusters"]}
    controls = {str(r["id"]): r for r in json.loads(control_manifest.read_text())["clusters"]}
    for row in rows:
        ident, cid, band = row["cluster_id"], row["control_id"], row["filter"]
        ikey = (ident, cid, band)
        pkey = (ident, cid)
        if ikey not in icl or pkey not in psf or (ident, band) not in target_profiles:
            raise ValueError(f"joined row does not resolve to source records: {ikey}")
        if ident not in targets or cid not in controls:
            raise ValueError(f"joined row does not resolve to source manifests: {ikey}")
        if str(controls[cid]["target_cluster_id"]) != ident:
            raise ValueError(f"control target link changed: {ikey}")
        source, psf_row = icl[ikey], psf[pkey]
        checks = (
            (float(row["icl_delta_njy_arcsec2"]), float(source["delta_sb_200_300"])),
            (float(row["common_coverage"]), float(source["coverage_200_300"])),
            (float(row["psf_core_delta_fwhm_major_arcsec"]),
             float(psf_row["target_minus_control_moment_fwhm_major_arcsec"])),
            (float(row["redshift"]), float(source["z"])),
            (float(row["richness"]), float(source["richness"])),
        )
        if any(not np.isclose(a, b, rtol=0, atol=1e-12) for a, b in checks):
            raise ValueError(f"joined values do not reproduce from source records: {ikey}")

    replicas = int(summary["bootstrap_replicas"])
    verified = {}
    for band, expected_n in (("i", 44), ("i2", 93)):
        selected = [row for row in rows if row["filter"] == band]
        if len(selected) != expected_n:
            raise ValueError(f"filter support changed for {band}: {len(selected)}")
        labels = connected_block_labels(selected)
        band_summary = summary["by_filter"][band]
        if int(band_summary["n_connected_target_control_spatial_blocks"]) != len(np.unique(labels)):
            raise ValueError(f"connected block count differs for {band}")
        models = {}
        for index, (name, covariates) in enumerate(COVARIATES.items()):
            observed = band_summary["models"][name]
            point = _independent_partial_r(selected, covariates)
            if not np.isclose(point, observed["partial_r"], rtol=0, atol=1e-12):
                raise ValueError(f"partial correlation does not reproduce: {band}/{name}")
            seed = int(observed["spatial_block_bootstrap"]["seed"])
            boot_ci, loo_range = _recompute_bootstrap(
                selected, labels, covariates, replicas, seed
            )
            saved_ci = observed["spatial_block_bootstrap"]["percentiles_2_5_50_97_5"]
            saved_loo = observed["leave_one_block_out"]["partial_r_range"]
            if not np.allclose(boot_ci, saved_ci, rtol=0, atol=1e-12):
                raise ValueError(f"block bootstrap does not reproduce: {band}/{name}")
            if not np.allclose(loo_range, saved_loo, rtol=0, atol=1e-12):
                raise ValueError(f"leave-one-block-out does not reproduce: {band}/{name}")
            models[name] = {"partial_r": point, "bootstrap_percentiles": boot_ci,
                            "leave_one_block_out_range": loo_range,
                            "n_replicas": replicas}
        verified[band] = {"n": len(selected), "n_blocks": len(np.unique(labels)),
                          "models": models}
    return {"status": "independently_verified", "rows": len(rows),
            "filters": verified, "science_measurement": False,
            "summary_sha256": sha256(summary_path),
            "joined_table_sha256": sha256(joined_path)}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--profile-dir", type=Path, required=True)
    parser.add_argument("--psf-diagnostic-dir", type=Path, required=True)
    parser.add_argument("--result-dir", type=Path, required=True)
    args = parser.parse_args()
    result = validate(args.profile_dir, args.psf_diagnostic_dir, args.result_dir)
    print(json.dumps(result, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
