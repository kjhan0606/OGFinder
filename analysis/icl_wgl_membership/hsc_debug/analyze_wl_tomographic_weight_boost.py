#!/usr/bin/env python3
"""Decompose HSC S19A source weight boosts and test a low-z-bin exclusion.

Uses only the approved local source-row exports and saved random-center
catalog. This is a diagnostic of redshift-dependent source-density/weight
selection; it is not a member-contamination correction or a science result.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import sys

import numpy as np
from scipy.stats import norm

HERE = Path(__file__).resolve().parent
if str(HERE) not in sys.path:
    sys.path.insert(0, str(HERE))

import analyze_wl_source_density_boost_blocks as base


DATA = base.DATA
DEFAULT_OUTPUT = DATA / "wl_tomographic_weight_boost_2026-11-05_04"
SOURCE_BINS = (1, 2, 3, 4)
METRICS_OF_INTEREST = ("sigma_crit_weight_boost", "estimator_weight_boost")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def subset_source_bin(pairs: dict, source_bin: int) -> dict:
    bins = np.asarray(pairs["source_zbin"], dtype=int)
    keep = bins == source_bin
    if not np.any(keep):
        raise ValueError(f"no source rows in HSC Y3 bin {source_bin}")
    return {key: values[keep] for key, values in pairs.items()}


def expand_tensor(dataset: dict, global_blocks: np.ndarray) -> np.ndarray:
    tensor = dataset["tensor"]
    local_blocks = np.asarray(dataset["blocks"], dtype=np.int64)
    local_to_global = np.searchsorted(global_blocks, local_blocks)
    expanded = np.zeros((tensor.shape[0], len(global_blocks), len(global_blocks), tensor.shape[-1]),
                        dtype=np.float64)
    for i, gi in enumerate(local_to_global):
        expanded[:, gi, local_to_global, :] = tensor[:, i, :, :]
    return expanded


def holm_adjust(p_values: list[float]) -> list[float]:
    p = np.asarray(p_values, dtype=float)
    order = np.argsort(p)
    adjusted = np.empty_like(p)
    running = 0.0
    for rank, index in enumerate(order):
        running = max(running, min(1.0, (len(p) - rank) * p[index]))
        adjusted[index] = running
    return adjusted.tolist()


def run(output: Path) -> dict:
    output = Path(output).resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite {output}")

    stack_25 = json.loads(base.STACK_2P5.read_text(encoding="utf-8"))
    stack_5 = json.loads(base.STACK_5.read_text(encoding="utf-8"))
    for saved in (stack_25, stack_5):
        if saved.get("science_measurement") is not False:
            raise ValueError("input stack must retain science_measurement=false")
        source = saved.get("input", {})
        if (source.get("n_z_sacc_sha256") != sha256(base.NZ_SACC)
                or source.get("b_mode_masked_value") != base.MASKED_BMODE_VALUE
                or source.get("cosmology") != {"H0_km_s_Mpc": 70.0, "Omega_m": 0.3}):
            raise ValueError("saved stack calibration differs from the pinned audit")

    random_pairs, random_prov = base.load_random_pairs(
        base.RANDOM_ROWS, base.RANDOM_PLAN, base.MASKED_BMODE_VALUE)
    lens_pairs, lens_prov = base.load_cluster_pairs(
        base.CLUSTER_ROWS, base.CLUSTER_PREFLIGHT, base.CLUSTER_PLAN,
        base.MASKED_BMODE_VALUE)
    matched_lens_ids = set(random_prov["matched_lens_ids"])
    lens_rows_ids = set(lens_pairs["lens_id"].astype(str))
    random_rows_lenses = {
        random_prov["random_center_matched_lens_ids"][center]
        for center, count in random_prov["random_center_expected_source_counts"].items()
        if int(count) > 0
    }
    common = matched_lens_ids & lens_rows_ids & random_rows_lenses
    sample_summary = stack_25["input"]["matched_sample"]
    if len(common) != int(sample_summary["n_common_shape_support_lenses"]):
        raise ValueError("source-row common lens sample differs from the saved 435-lens stack")
    common_hash = hashlib.sha256("\n".join(sorted(common)).encode("utf-8")).hexdigest()
    if common_hash != sample_summary.get("common_lens_ids_sha256"):
        raise ValueError("common lens IDs differ from the approved saved stack")

    center_doc = json.loads(base.RANDOM_CENTERS.read_text(encoding="utf-8"))
    centers = base.source_centers_for_sample(
        center_doc, random_prov["random_center_matched_lens_ids"], common)
    manifest = json.loads(base.SAMPLE_MANIFEST.read_text(encoding="utf-8"))
    lens_by_id = {str(row["id"]): row for row in manifest["clusters"]}
    lens_positions = [lens_by_id[ident] for ident in sorted(common)]
    center_to_lens = random_prov["random_center_matched_lens_ids"]
    z_nz, nz_by_bin = base.load_hsc_y3_nz_sacc(base.NZ_SACC)
    from ogfkit.wlcluster import Cosmology
    cosmology = Cosmology(H0=70.0, Om=0.3)

    rows = []
    summaries = {}
    for width in (5.0, 2.5):
        saved = stack_5 if width == 5.0 else stack_25
        raw_datasets = {}
        block_parts = [
            base.equal_area_block_ids(
                [row["ra_deg"] for row in lens_positions],
                [row["dec_deg"] for row in lens_positions], cell_width_deg=width),
            base.equal_area_block_ids(
                [row["ra_deg"] for row in centers],
                [row["dec_deg"] for row in centers], cell_width_deg=width),
        ]
        for source_bin in SOURCE_BINS:
            lens = base.reduce_pairs(
                subset_source_bin(lens_pairs, source_bin), common, None,
                z_nz, nz_by_bin, cosmology, width)
            random = base.reduce_pairs(
                subset_source_bin(random_pairs, source_bin), common, center_to_lens,
                z_nz, nz_by_bin, cosmology, width)
            raw_datasets[source_bin] = (lens, random)
            block_parts.extend([lens["blocks"], random["blocks"]])
        global_blocks = np.unique(np.concatenate(block_parts)).astype(np.int64)
        lens_center_blocks = np.searchsorted(global_blocks, block_parts[0])
        random_center_blocks = np.searchsorted(global_blocks, block_parts[1])

        per_bin_draws = []
        per_bin_values = []
        expanded_lens_tensors = []
        expanded_random_tensors = []
        lens_counts_by_bin = []
        random_counts_by_bin = []
        for source_bin in SOURCE_BINS:
            lens, random = raw_datasets[source_bin]
            lens_jk = {
                "tensor": expand_tensor(lens, global_blocks),
                "cluster_center_blocks": lens_center_blocks,
                "n_lens_centers": len(common),
                "block_ids": global_blocks,
            }
            random_jk = {
                "tensor": expand_tensor(random, global_blocks),
                "random_center_blocks": random_center_blocks,
                "n_random_centers": len(centers),
                "block_ids": global_blocks,
            }
            expanded_lens_tensors.append(lens_jk["tensor"])
            expanded_random_tensors.append(random_jk["tensor"])
            metrics, joint, draws = base.jackknife_summary(
                lens_jk, random_jk, len(global_blocks))
            per_bin_draws.append((joint, draws))
            per_bin_values.append(metrics)
            lens_counts_by_bin.append(lens["selected_pair_counts"])
            random_counts_by_bin.append(random["selected_pair_counts"])
            expected_lens = sum(int(row["n_lens_pairs"]) for row in
                                saved["paired_spatial_jackknife"]["profile"])
            rows_for_bin = int(lens["n_selected_rows"] + random["n_selected_rows"])
            rows.append({
                "block_width_deg": width,
                "source_zbin": source_bin,
                "selected_cluster_rows": int(lens["n_selected_rows"]),
                "selected_random_rows": int(random["n_selected_rows"]),
                "selected_rows_combined": rows_for_bin,
                "lens_pair_counts_by_radius": lens["selected_pair_counts"].tolist(),
                "random_pair_counts_by_radius": random["selected_pair_counts"].tolist(),
                "full_sample_lens_pair_count_across_all_source_bins": expected_lens,
                "count_boost": metrics["count_boost"],
                "sigma_crit_weight_boost": metrics["sigma_crit_weight_boost"],
                "estimator_weight_boost": metrics["estimator_weight_boost"],
            })

        # The saved four-bin stack excludes zbin 0. Require the four separate
        # source-bin row counts to reproduce its selected pair totals exactly.
        expected_rows = saved["paired_spatial_jackknife"]["profile"]
        lens_bin_sum = np.sum(np.stack(lens_counts_by_bin), axis=0)
        random_bin_sum = np.sum(np.stack(random_counts_by_bin), axis=0)
        if lens_bin_sum.tolist() != [int(row["n_lens_pairs"]) for row in expected_rows]:
            raise ValueError(f"four source bins do not reproduce saved lens pair counts at {width:g} degrees")
        if random_bin_sum.tolist() != [int(row["n_random_pairs"]) for row in expected_rows]:
            raise ValueError(f"four source bins do not reproduce saved random pair counts at {width:g} degrees")

        # Paired cross-bin covariance for the exact estimator-weight boost.
        selected_replicas = []
        parameter_order = []
        for source_bin, (joint, draws) in zip(SOURCE_BINS, per_bin_draws):
            order = joint["parameter_order"]
            keep = [i for i, item in enumerate(order)
                    if item["metric"] in METRICS_OF_INTEREST]
            selected_replicas.append(draws[:, keep])
            parameter_order.extend({"source_zbin": source_bin, **order[i]} for i in keep)
        combined = np.concatenate(selected_replicas, axis=1)
        covariance = base.jackknife_covariance(combined)
        sigma = np.sqrt(np.diag(covariance))
        denom = np.outer(sigma, sigma)
        correlation = np.divide(covariance, denom, out=np.zeros_like(covariance), where=denom > 0)
        contrasts = []
        contrast_vectors = []
        p_values = []
        for radial_bin in range(len(base.RADIAL_EDGES) - 1):
            for comparison_bin in (2, 3, 4):
                z1_point = per_bin_values[0]["estimator_weight_boost"][radial_bin]["point"]
                other_point = per_bin_values[comparison_bin - 1]["estimator_weight_boost"][radial_bin]["point"]
                z1_replicas = selected_replicas[0]
                other_replicas = selected_replicas[comparison_bin - 1]
                # Selected columns are SigmaCrit weight bins 0..6, then exact
                # estimator-weight bins 0..6 for each source-redshift bin.
                idx_z1 = 7 + radial_bin
                idx_other = (comparison_bin - 1) * 14 + 7 + radial_bin
                vector = np.zeros(len(parameter_order), dtype=float)
                vector[idx_z1] = 1.0
                vector[idx_other] = -1.0
                contrast_vectors.append(vector)
                variance = float(vector @ covariance @ vector)
                if variance <= 0 or not np.isfinite(variance):
                    raise ValueError("source-bin estimator-weight contrast has invalid jackknife variance")
                sigma_contrast = float(np.sqrt(variance))
                difference = float(z1_point - other_point)
                z_score = difference / sigma_contrast
                p_value = float(2.0 * norm.sf(abs(z_score)))
                p_values.append(p_value)
                contrasts.append({
                    "radial_bin": radial_bin,
                    "radius_mpc": [float(base.RADIAL_EDGES[radial_bin]),
                                   float(base.RADIAL_EDGES[radial_bin + 1])],
                    "contrast": "source_zbin_1_minus_source_zbin_%d" % comparison_bin,
                    "difference": difference,
                    "paired_jackknife_sigma": sigma_contrast,
                    "nominal_z_score": z_score,
                    "nominal_two_sided_normal_p": p_value,
                    "nominal_95pct_interval": [difference - 1.96 * sigma_contrast,
                                               difference + 1.96 * sigma_contrast],
                    "holm_family_p": None,
                })
        adjusted = holm_adjust(p_values)
        for item, value in zip(contrasts, adjusted):
            item["holm_family_p"] = value
        contrast_matrix = np.asarray(contrast_vectors)
        contrast_covariance = contrast_matrix @ covariance @ contrast_matrix.T

        def stack_profile(tensor: np.ndarray, center_blocks: np.ndarray,
                          n_centers: int) -> dict:
            return {"tensor": tensor, "cluster_center_blocks": center_blocks,
                    "random_center_blocks": center_blocks,
                    "n_lens_centers": n_centers, "n_random_centers": n_centers,
                    "block_ids": global_blocks}

        all_lens_tensor = np.sum(np.stack(expanded_lens_tensors), axis=0)
        all_random_tensor = np.sum(np.stack(expanded_random_tensors), axis=0)
        bins_234_lens_tensor = np.sum(np.stack(expanded_lens_tensors[1:]), axis=0)
        bins_234_random_tensor = np.sum(np.stack(expanded_random_tensors[1:]), axis=0)
        lens_profile_jk = stack_profile(all_lens_tensor, lens_center_blocks, len(common))
        random_profile_jk = stack_profile(all_random_tensor, random_center_blocks, len(centers))
        lens_234_jk = stack_profile(bins_234_lens_tensor, lens_center_blocks, len(common))
        random_234_jk = stack_profile(bins_234_random_tensor, random_center_blocks, len(centers))
        lens_full, lens_full_replicas = base.delta_sigma_leave_one_block(lens_profile_jk)
        random_full, random_full_replicas = base.delta_sigma_leave_one_block(random_profile_jk)
        full_profile = lens_full - random_full
        full_replicas = lens_full_replicas - random_full_replicas
        lens_234, lens_234_replicas = base.delta_sigma_leave_one_block(lens_234_jk)
        random_234, random_234_replicas = base.delta_sigma_leave_one_block(random_234_jk)
        profile_234 = lens_234 - random_234
        replicas_234 = lens_234_replicas - random_234_replicas
        profile_difference = profile_234 - full_profile
        difference_replicas = replicas_234 - full_replicas
        full_covariance = base.jackknife_covariance(full_replicas.reshape(len(global_blocks), -1))
        covariance_234 = base.jackknife_covariance(replicas_234.reshape(len(global_blocks), -1))
        difference_covariance = base.jackknife_covariance(
            difference_replicas.reshape(len(global_blocks), -1))
        saved_profile = np.asarray([
            [row["delta_sigma_random_subtracted_msun_mpc2"],
             row["cross_delta_sigma_random_subtracted_msun_mpc2"]]
            for row in saved["paired_spatial_jackknife"]["profile"]
        ], dtype=float)
        if not np.allclose(full_profile, saved_profile, rtol=5e-9, atol=1e-8):
            raise ValueError(f"tomographic source-bin sum does not reproduce saved profile at {width:g} degrees")

        from scipy.stats import chi2
        def null_summary(profile: np.ndarray, cov: np.ndarray, component: int) -> dict:
            indices = np.arange(component, 2 * len(profile), 2)
            selected_cov = cov[np.ix_(indices, indices)]
            values = profile[:, component]
            rank = int(np.linalg.matrix_rank(selected_cov))
            statistic = float(values @ np.linalg.pinv(selected_cov, hermitian=True) @ values)
            return {"chi2": statistic, "dof_covariance_rank": rank,
                    "nominal_chi2_p": float(chi2.sf(statistic, rank)),
                    "covariance_condition_number": float(np.linalg.cond(selected_cov))}

        bins_234_sensitivity = {
            "omitted_source_zbin": 1,
            "retained_source_zbins": [2, 3, 4],
            "profile_msun_mpc2": profile_234.tolist(),
            "profile_covariance_msun2_mpc4": covariance_234.tolist(),
            "difference_from_all_bins_msun_mpc2": profile_difference.tolist(),
            "difference_covariance_msun2_mpc4": difference_covariance.tolist(),
            "all_bins_profile_reproduces_saved_stack": True,
            "null_checks": {
                "all_bins_tangential": null_summary(full_profile, full_covariance, 0),
                "all_bins_cross": null_summary(full_profile, full_covariance, 1),
                "bins_2_to_4_tangential": null_summary(profile_234, covariance_234, 0),
                "bins_2_to_4_cross": null_summary(profile_234, covariance_234, 1),
                "bins_2_to_4_minus_all_bins_tangential_change": null_summary(
                    profile_difference, difference_covariance, 0),
                "bins_2_to_4_minus_all_bins_cross_change": null_summary(
                    profile_difference, difference_covariance, 1),
            },
            "difference_delete_one_block_sigma": np.sqrt(
                np.diag(difference_covariance)).reshape(len(profile_234), 2).tolist(),
            "retained_cluster_pairs_by_radius": np.sum(np.stack(lens_counts_by_bin[1:]), axis=0).tolist(),
            "retained_random_pairs_by_radius": np.sum(np.stack(random_counts_by_bin[1:]), axis=0).tolist(),
            "omitted_zbin_1_cluster_pairs_by_radius": lens_counts_by_bin[0].tolist(),
            "omitted_zbin_1_random_pairs_by_radius": random_counts_by_bin[0].tolist(),
            "interpretation_limit": "source-z-bin exclusion is a sample-sensitivity test; it does not correct residual photo-z dilution or validate random-aperture matching",
        }
        summaries[str(width)] = {
            "n_lenses": len(common), "n_random_centers": len(centers),
            "n_endpoint_blocks": int(len(global_blocks)),
            "block_ids": global_blocks.tolist(),
            "saved_all_bin_pair_counts_referenced": [
                int(row["n_lens_pairs"]) for row in saved["paired_spatial_jackknife"]["profile"]
            ],
            "joint_cross_tomographic_covariance": {
                "method": "paired delete-one equal-area sky-cell jackknife; same lens/random support and spatial deletion across all four source bins",
                "parameter_order": parameter_order,
                "covariance": covariance.tolist(),
                "correlation": correlation.tolist(),
            },
            "source_zbin_1_vs_bins_2_to_4_estimator_weight_contrasts": {
                "family": "three comparisons across seven radial bins; Holm correction across 21 tests",
                "covariance": contrast_covariance.tolist(),
                "results": contrasts,
                "interpretation_limit": "nominal normal p-values use the jackknife covariance and do not propagate n(z) calibration uncertainty or validate full-aperture random matching",
            },
            "four_source_bin_pair_counts_reproduce_saved_stack": True,
            "source_zbin_1_exclusion_sensitivity": bins_234_sensitivity,
        }

    output.mkdir(parents=True, exist_ok=False)
    import csv
    table_path = output / "tomographic_weight_boost.tsv"
    columns = [
        "block_width_deg", "source_zbin", "selected_cluster_rows", "selected_random_rows",
        "selected_rows_combined", "lens_pair_counts_by_radius", "random_pair_counts_by_radius",
        "full_sample_lens_pair_count_across_all_source_bins", "count_boost",
        "sigma_crit_weight_boost", "estimator_weight_boost",
    ]
    with table_path.open("x", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=columns, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: json.dumps(row[key], separators=(",", ":"))
                             if isinstance(row[key], (dict, list)) else row[key]
                             for key in columns})

    input_paths = [
        base.NZ_SACC, base.RANDOM_CENTERS, base.SAMPLE_MANIFEST,
        base.RANDOM_PLAN / "random_shape_row_plan.json",
        base.RANDOM_ROWS / "row_summary.json", base.RANDOM_ROWS / "row_status.json",
        base.CLUSTER_PLAN / "row_export_plan.json",
        base.CLUSTER_ROWS / "row_summary.json", base.CLUSTER_ROWS / "row_status.json",
        base.STACK_5, base.STACK_2P5,
    ]
    result = {
        "schema": "astrafek-wl-tomographic-weight-boost-v1",
        "status": "complete_diagnostic",
        "science_measurement": False,
        "boost_correction_applied": False,
        "data_requested": False,
        "hsc_network_used": False,
        "images_read": False,
        "matched_sample": {
            "n_common_lenses": len(common), "n_random_centers": len(centers),
            "random_centers_per_lens": 5,
            "common_lens_ids_sha256": common_hash,
        },
        "source_redshift_bins": list(SOURCE_BINS),
        "source_bin_zero_excluded": True,
        "estimator": "per-bin source count, SigmaCritInv-squared weight, and exact S19A response-adjusted estimator weight per planned lens/random center",
        "block_width_results": summaries,
        "profile_rows": rows,
        "profile_table": str(table_path),
        "profile_table_sha256": sha256(table_path),
        "input_sha256": {str(path.resolve()): sha256(path) for path in input_paths},
        "limitations": [
            "A redshift-bin weight excess does not identify cluster-member contamination; magnification, depth, masks, PSF, sky response, and source-selection differences can also change these ratios.",
            "Random-center full-aperture depth, mask, PSF, and sky matching are not established.",
            "HSC Y3 n(z) uncertainty and source-redshift calibration uncertainty are not propagated.",
            "Do not use these ratios as a multiplicative correction to DeltaSigma; the saved science gate remains false.",
            "Equal-area RA/sin(Dec) blocks are not verified contiguous HSC survey regions.",
        ],
    }
    summary_path = output / "summary.json"
    summary_path.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    result["summary_path"] = str(summary_path)
    result["summary_sha256"] = sha256(summary_path)
    return result


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    result = run(args.output)
    summary = {
        "status": result["status"],
        "science_measurement": result["science_measurement"],
        "matched_sample": result["matched_sample"],
        "estimator_weight_boost_by_source_bin": {
            width: {
                str(row["source_zbin"]): [
                    round(float(value["point"]), 5) for value in row["estimator_weight_boost"]
                ] for row in result["profile_rows"] if float(row["block_width_deg"]) == float(width)
            } for width in ("5.0", "2.5")
        },
        "output": result["summary_path"],
        "summary_sha256": result["summary_sha256"],
    }
    print(json.dumps(summary, ensure_ascii=False, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
