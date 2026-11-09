#!/usr/bin/env python3
"""Audit satellite-light leakage into the measured ICL--satellite association.

Join the observed 94-pair association sample to the independently audited
position-dependent-PSF injection/redetection responses. Injecting cluster light
into target images while holding matched control images fixed measures the
incremental response of the target-minus-control ICL contrast. On the exact
overlap, compare transfer from injected satellite envelopes into the ICL
annulus with transfer into the measured satellite profile. A paired
block-bootstrap stress test asks how the association changes under a
median-response leakage subtraction. The stress test is not a calibrated
correction because one median coefficient does not propagate the broad,
pair-dependent transfer response or the support selection.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np

DATA = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
INJECTION_DIR = DATA / "joint_icl_satellite_psf_redetection_icl_supported224_2026-11-06_01"
INJECTION_SUMMARY = INJECTION_DIR / "summary.json"
INJECTION_AUDIT = INJECTION_DIR / "psf_redetection_audit.json"
ASSOCIATION = DATA / "icl_satellite_pixel_noise_bias_supported_94_2026-11-06_03/association_rows.tsv"
SUPPORT = DATA / "icl_satellite_ipw_variance_sensitivity_2026-11-06_01/support_variance_join.tsv"
DEFAULT_OUTPUT = DATA / "icl_satellite_injection_cross_talk_2026-11-06_02"
EXPECTED_INJECTION_SUMMARY_SHA = "be3c2e8cb22185c06768d9022ae3e6192072f180fa54d55e21b85920e5572750"
EXPECTED_INJECTION_AUDIT_SHA = "b91783ee4e1582c24676d9a9bd85aacec2cebafa46d76980e909fcbfecbc0bd1"
BOOTSTRAP_REPLICAS = 10000
BOOTSTRAP_SEED = 20261106


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_json(path: Path) -> dict:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def read_tsv(path: Path) -> list[dict[str, str]]:
    with Path(path).open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f"refusing empty output table: {path}")
    with Path(path).open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]), delimiter="\t",
                                extrasaction="raise", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def partial_r(x: np.ndarray, y: np.ndarray, z: np.ndarray,
              log_richness: np.ndarray) -> float:
    if len(x) < 6 or not all(np.isfinite(v).all() for v in (x, y, z, log_richness)):
        return float("nan")
    design = np.column_stack((np.ones(len(x)), z, log_richness))
    if np.linalg.matrix_rank(design) != 3:
        return float("nan")
    rx = x - design @ np.linalg.lstsq(design, x, rcond=None)[0]
    ry = y - design @ np.linalg.lstsq(design, y, rcond=None)[0]
    denom = float(np.linalg.norm(rx) * np.linalg.norm(ry))
    return float(rx @ ry / denom) if denom > 0 else float("nan")


def interval(values: list[float]) -> list[float] | None:
    array = np.asarray(values, dtype=float)
    array = array[np.isfinite(array)]
    return np.quantile(array, [0.025, 0.975]).tolist() if len(array) else None


def robust_summary(values) -> dict:
    array = np.asarray(values, dtype=float)
    array = array[np.isfinite(array)]
    if len(array) == 0:
        return {"n": 0, "median": None, "q16": None, "q84": None}
    return {"n": int(len(array)), "median": float(np.median(array)),
            "q16": float(np.quantile(array, .16)),
            "q84": float(np.quantile(array, .84))}


def summarize_filter(rows: list[dict], full_rows: list[dict], band: str,
                     replicas: int, seed: int) -> dict:
    selected = [row for row in rows if row["filter"] == band]
    full = [row for row in full_rows if row["filter"] == band]
    if len(selected) < 6:
        raise ValueError(f"too few common injection/association pairs for {band}: {len(selected)}")
    blocks = sorted({row["connected_5deg_block"] for row in selected})
    if len(blocks) < 6:
        raise ValueError(f"too few connected blocks for {band}: {len(blocks)}")
    indices = {block: np.asarray([i for i, row in enumerate(selected)
                                  if row["connected_5deg_block"] == block], dtype=int)
               for block in blocks}
    rng = np.random.default_rng(seed)
    base_draws, leak_draws, delta_draws, k_draws = [], [], [], []
    for _ in range(replicas):
        chosen_blocks = rng.choice(blocks, size=len(blocks), replace=True)
        idx = np.concatenate([indices[str(block)] for block in chosen_blocks])
        sample = [selected[int(i)] for i in idx]
        x = np.asarray([r["icl_delta_njy_arcsec2"] for r in sample])
        y = np.asarray([r["satellite_delta_njy_arcsec2"] for r in sample])
        z = np.asarray([r["redshift"] for r in sample])
        rich = np.asarray([r["log_richness"] for r in sample])
        base = partial_r(x, y, z, rich)
        b = np.asarray([r["icl_from_satellite_response"] for r in sample])
        d = np.asarray([r["satellite_from_satellite_response"] for r in sample])
        denominator = float(np.median(d))
        numerator = float(np.median(b))
        if not math.isfinite(denominator) or abs(denominator) < 1e-8:
            continue
        k = numerator / denominator
        stress = partial_r(x-k*y, y, z, rich)
        if np.isfinite(base) and np.isfinite(stress):
            base_draws.append(float(base)); leak_draws.append(float(stress))
            delta_draws.append(float(stress-base)); k_draws.append(float(k))

    sel_arrays = {
        "icl_delta": np.asarray([row["icl_delta_njy_arcsec2"] for row in selected], dtype=float),
        "satellite_delta": np.asarray([row["satellite_delta_njy_arcsec2"] for row in selected], dtype=float),
        "redshift": np.asarray([row["redshift"] for row in selected], dtype=float),
        "log_richness": np.asarray([row["log_richness"] for row in selected], dtype=float),
        **{key: np.asarray([row[key] for row in selected], dtype=float)
           for key in ("icl_from_diffuse_response", "icl_from_satellite_response",
                       "satellite_from_diffuse_response", "satellite_from_satellite_response",
                       "matrix_determinant", "matrix_condition", "leakage_ratio")},
    }
    b_med = float(np.median(sel_arrays["icl_from_satellite_response"]))
    d_med = float(np.median(sel_arrays["satellite_from_satellite_response"]))
    k_ratio_of_medians = b_med/d_med if abs(d_med) > 1e-8 else None
    x, y = sel_arrays["icl_delta"], sel_arrays["satellite_delta"]
    z, rich = sel_arrays["redshift"], sel_arrays["log_richness"]
    base_r = partial_r(x, y, z, rich)
    stress_r = (partial_r(x-float(k_ratio_of_medians)*y, y, z, rich)
                if k_ratio_of_medians is not None else float("nan"))
    full_r = partial_r(
        np.asarray([float(r["icl_delta_raw_njy_arcsec2"]) for r in full]),
        np.asarray([float(r["satellite_delta_raw_njy_arcsec2"]) for r in full]),
        np.asarray([float(r["z"]) for r in full]),
        np.log(np.asarray([float(r["richness"]) for r in full])),
    )
    valid_fraction = len(base_draws)/replicas
    if valid_fraction < 0.8:
        raise ValueError(f"too few valid paired bootstrap draws for {band}: {len(base_draws)}/{replicas}")
    return {
        "n_association_pairs_full_filter": len(full),
        "n_common_injection_pairs": len(selected),
        "n_connected_5deg_blocks": len(blocks),
        "full_94_sample_partial_r_z_logrichness": full_r,
        "common_injection_subset_partial_r_z_logrichness": base_r,
        "median_response_coefficients": {
            "icl_from_diffuse_injection": robust_summary(sel_arrays["icl_from_diffuse_response"]),
            "icl_from_satellite_injection": robust_summary(sel_arrays["icl_from_satellite_response"]),
            "satellite_from_diffuse_injection": robust_summary(sel_arrays["satellite_from_diffuse_response"]),
            "satellite_from_satellite_injection": robust_summary(sel_arrays["satellite_from_satellite_response"]),
        },
        "satellite_to_icl_leakage_ratio_b_over_d": {
            "median_of_pairwise_ratios_q16_q84": robust_summary(sel_arrays["leakage_ratio"]),
            "ratio_of_median_responses": k_ratio_of_medians,
        },
        "two_component_transfer_matrix": {
            "matrix_order": "[[ICL_annulus<-diffuse, ICL_annulus<-satellite], [satellite_profile<-diffuse, satellite_profile<-satellite]]",
            "median_determinant": float(np.median(sel_arrays["matrix_determinant"])),
            "condition_number_q16_median_q84": [float(np.quantile(sel_arrays["matrix_condition"], q)) for q in (.16, .5, .84)],
        },
        "leakage_stress_test_not_a_correction": {
            "definition": "ICL target-minus-control contrast minus ratio-of-median target-side satellite-to-ICL response times measured satellite value; matched controls remain fixed because the injected cluster signal is absent from random control fields",
            "stress_test_partial_r_with_global_ratio_of_median_response": stress_r,
            "paired_block_bootstrap_base_r_ci95": interval(base_draws),
            "paired_block_bootstrap_stress_r_ci95": interval(leak_draws),
            "paired_block_bootstrap_stress_minus_base_ci95": interval(delta_draws),
            "bootstrap_ratio_of_median_responses_ci95": interval(k_draws),
            "valid_replicas": len(base_draws), "replicas_requested": replicas,
            "seed": seed,
        },
        "interpretation": "Target-only injections measure the incremental response of the target-minus-control ICL contrast when matched controls remain fixed. The median-ratio stress test is not a calibrated correction because it assumes one leakage coefficient despite broad pair-dependent responses and selected injection support.",
    }


def run(output: Path, replicas: int = BOOTSTRAP_REPLICAS,
        seed: int = BOOTSTRAP_SEED) -> dict:
    output = Path(output).resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite {output}")
    if sha256(INJECTION_SUMMARY) != EXPECTED_INJECTION_SUMMARY_SHA:
        raise ValueError("injection batch summary hash changed")
    if sha256(INJECTION_AUDIT) != EXPECTED_INJECTION_AUDIT_SHA:
        raise ValueError("independent injection audit hash changed")
    summary, audit = read_json(INJECTION_SUMMARY), read_json(INJECTION_AUDIT)
    if (summary.get("status") != "complete"
            or summary.get("science_measurement") is not False
            or summary.get("observational_validation_passed") is not False
            or audit.get("n_frame_results_checked") != 137
            or audit.get("n_errors_reported_by_batch") != 0):
        raise ValueError("PSF/redetection batch is not independently verified")
    injection_records = {}
    for record in summary.get("records", []):
        key = (str(record["cluster_id"]), str(record["filter"]))
        path = Path(record["path"])
        if not path.is_file() or sha256(path) != record.get("sha256"):
            raise ValueError(f"per-cluster injection record changed: {key}")
        detail = read_json(path)
        if (detail.get("cluster_id") != key[0] or detail.get("filter") != key[1]
                or detail.get("science_measurement") is not False):
            raise ValueError(f"per-cluster injection scope mismatch: {key}")
        injection_records[key] = {**record, "detail": detail}
    if len(injection_records) != 137:
        raise ValueError("expected 137 unique audited injection records")

    associations = read_tsv(ASSOCIATION)
    support_rows = read_tsv(SUPPORT)
    if len(associations) != 94 or len(support_rows) != 224:
        raise ValueError("association or parent support sample size differs from pinned audit")
    support_by_key = {(str(r["cluster_id"]), str(r["filter"])): r for r in support_rows}
    association_by_key = {(str(r["cluster_id"]), str(r["filter"])): r for r in associations}
    if len(association_by_key) != 94 or len(support_by_key) != 224:
        raise ValueError("duplicate cluster/filter keys in association/support tables")
    rows, unmatched = [], []
    for key, assoc in association_by_key.items():
        parent = support_by_key.get(key)
        if (parent is None or parent.get("selected_50_75") != "True"
                or abs(float(assoc["z"])-float(parent["redshift"])) > 1e-8
                or abs(float(assoc["richness"])-float(parent["richness"])) > 1e-8):
            raise ValueError(f"association key is not the exact audited satellite support row: {key}")
        injection = injection_records.get(key)
        if injection is None:
            unmatched.append({"cluster_id": key[0], "filter": key[1], "reason": "no approved PSF injection record"})
            continue
        detail = injection["detail"]
        gaps = detail.get("satellite_support_gaps", [])
        if gaps:
            unmatched.append({"cluster_id": key[0], "filter": key[1], "reason": "satellite annulus support gap"})
            continue
        t = detail["transfer_response"]
        names = ("icl_response_to_diffuse_unit_injection",
                 "icl_response_to_all_satellite_envelopes",
                 "satellite_response_to_smooth_diffuse",
                 "satellite_response_to_satellite_unit_injection")
        values = [t.get(name) for name in names]
        if any(value is None or not math.isfinite(float(value)) for value in values):
            unmatched.append({"cluster_id": key[0], "filter": key[1], "reason": "non-finite transfer response"})
            continue
        a, b, c, d = map(float, values)
        ratio = b/d if abs(d) > 1e-10 else float("nan")
        matrix = np.asarray([[a, b], [c, d]], dtype=float)
        rows.append({
            "cluster_id": key[0], "filter": key[1],
            "connected_5deg_block": str(parent["connected_5deg_block"]),
            "icl_delta_njy_arcsec2": float(assoc["icl_delta_raw_njy_arcsec2"]),
            "satellite_delta_njy_arcsec2": float(assoc["satellite_delta_raw_njy_arcsec2"]),
            "redshift": float(assoc["z"]), "log_richness": math.log(float(assoc["richness"])),
            "icl_from_diffuse_response": a, "icl_from_satellite_response": b,
            "satellite_from_diffuse_response": c, "satellite_from_satellite_response": d,
            "matrix_determinant": float(np.linalg.det(matrix)),
            "matrix_condition": float(np.linalg.cond(matrix)),
            "leakage_ratio": ratio,
            "injection_record_sha256": injection["sha256"],
            "target_only_injection": True,
        })
    if not rows:
        raise ValueError("no supported association pairs overlap valid injection responses")
    input_hashes = {str(path.resolve()): sha256(path)
                    for path in (INJECTION_SUMMARY, INJECTION_AUDIT, ASSOCIATION, SUPPORT)}
    results = {}
    for index, band in enumerate(("i", "i2")):
        results[band] = summarize_filter(rows, list(association_by_key.values()), band,
                                         replicas, seed+index*1000)
    output.mkdir(parents=True, exist_ok=False)
    row_path = output/"association_injection_cross_talk.tsv"
    unmatched_path = output/"unmatched_association_pairs.tsv"
    write_tsv(row_path, sorted(rows, key=lambda r: (r["filter"], int(r["cluster_id"]))))
    write_tsv(unmatched_path, unmatched or [{"cluster_id": "", "filter": "", "reason": "none"}])
    result = {
        "schema": "astrafek-icl-satellite-injection-cross-talk-v2",
        "status": "complete_exploratory_diagnostic", "science_measurement": False,
        "hsc_network_used": False, "data_requested": False,
        "n_observed_association_pairs": len(associations),
        "n_pairs_with_valid_injection_transfer": len(rows),
        "n_unmatched_association_pairs": len(unmatched),
        "unmatched_reasons": {reason: sum(item["reason"] == reason for item in unmatched)
                               for reason in sorted({item["reason"] for item in unmatched})},
        "n_frames_in_independent_injection_audit": audit["n_frame_results_checked"],
        "by_filter": results,
        "method": "join observed 94-pair ICL/satellite target-minus-control measurements to the independently audited target-only PSF-convolved redetection responses; calculate the target-side leakage ratio b/d and paired 5-degree-block bootstrap under a median-response leakage stress test",
        "limitations": [
            "Target-only injections are appropriate for the incremental cluster-signal response of the target-minus-control estimator when matched random controls remain fixed; the current response support covers only 58 of 94 observed association pairs.",
            "Only 58/94 observed satellite-supported pairs have valid injection responses (14 i, 44 i2); this overlap is selected by PSF availability and image coverage.",
            "Satellite response distributions are broad and cross zero, making b/d an unstable de-mixing coefficient.",
            "The leakage stress test is not a correction or a population-calibrated transfer function; it omits injected Poisson noise, coadd weights, extended PSF wings, and HSC sky transfer.",
            "The observed association still has unverified missingness and spatial covariance; no science gate is promoted.",
        ],
        "input_sha256": input_hashes,
        "output_sha256": {"association_injection_cross_talk.tsv": sha256(row_path),
                           "unmatched_association_pairs.tsv": sha256(unmatched_path)},
        "analysis_script_sha256": sha256(Path(__file__).resolve()),
    }
    summary_path = output/"summary.json"
    summary_path.write_text(json.dumps(result, ensure_ascii=False, indent=2,
                                       sort_keys=True, allow_nan=False)+"\n", encoding="utf-8")
    (output/"output_hashes.json").write_text(json.dumps({"summary_sha256": sha256(summary_path),
        "association_injection_cross_talk_sha256": sha256(row_path),
        "unmatched_association_pairs_sha256": sha256(unmatched_path)}, indent=2)+"\n", encoding="utf-8")
    print(json.dumps({"status": result["status"], "science_measurement": False,
                      "pairs": result["n_pairs_with_valid_injection_transfer"],
                      "by_filter": {band: {"n": value["n_common_injection_pairs"],
                          "r": value["common_injection_subset_partial_r_z_logrichness"],
                          "k_median_ratio": value["satellite_to_icl_leakage_ratio_b_over_d"]["ratio_of_median_responses"],
                          "stress_r": value["leakage_stress_test_not_a_correction"]["stress_test_partial_r_with_global_ratio_of_median_response"],
                          "delta_ci95": value["leakage_stress_test_not_a_correction"]["paired_block_bootstrap_stress_minus_base_ci95"]}
                          for band, value in results.items()},
                      "output": str(output)}, ensure_ascii=False, indent=2), flush=True)
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    parser.add_argument("--replicas", type=int, default=BOOTSTRAP_REPLICAS)
    parser.add_argument("--seed", type=int, default=BOOTSTRAP_SEED)
    args = parser.parse_args()
    if args.replicas < 1000:
        parser.error("at least 1000 block-bootstrap replicates are required")
    run(args.output, replicas=args.replicas, seed=args.seed)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
