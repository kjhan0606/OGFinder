#!/usr/bin/env python3
"""Build a hash-provenance ICL and weak-lensing evidence report.

The builder reads only already saved Astrafek summaries and compact audit
tables. It makes no HSC requests, does not read image pixels, and preserves
the exploratory science status of all component analyses.
"""
from __future__ import annotations

import argparse
import csv
import json
import hashlib
import ast
from pathlib import Path

DATA_ROOT = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
QUERY_AUDIT = DATA_ROOT / "icl_expansion_coverage_plan_2026-11-06_01/cas_coverage_metadata/submission_attempt_audit.json"
S23B_QUERY_AUDIT = DATA_ROOT / "icl_s23b_coverage_plan_2026-11-06_04/cas_coverage_metadata/submission_attempt_audit.json"
SOURCES = {
    "icl_support94": "icl_satellite_supported_94_association_2026-11-05_01/summary.json",
    "icl_pixelnoise94": "icl_satellite_pixel_noise_bias_supported_94_2026-11-06_04/summary.json",
    "icl_ipw_pixelnoise": "icl_satellite_ipw_pixelnoise_joint_2026-11-06_01/summary.json",
    "icl_membercount": "icl_satellite_membercount_adjustment_2026-11-06_02/summary.json",
    "icl_positionmix": "icl_satellite_clustercentric_position_adjustment_2026-11-06_01/summary.json",
    "icl_aperture_covariance": "icl_aperture_covariance_200_300_2026-11-06_01/summary.json",
    "icl_injection_crosstalk": "icl_satellite_injection_cross_talk_2026-11-06_02/summary.json",
    "icl_transfer_grid": "icl_satellite_amplitude_grid_transfer_matrix_2026-11-06_01/summary.json",
    "icl_shell_exclusion": "icl_hybrid_shell_exclusion_audit_2026-11-06_01/summary.json",
    "frame_metadata": "icl_frame_metadata_variance_association_2026-11-06_02/summary.json",
    "frame_metadata_satellite": "icl_frame_metadata_satellite_relation_2026-11-06_03/summary.json",
    "icl_nonlinear_adjustment": "icl_satellite_nonlinear_adjustment_2026-11-06_01/summary.json",
    "icl_spatial_cv": "icl_satellite_spatial_cv_2026-11-06_01/summary.json",
    "icl_fixed_profile": "icl_satellite_fixed_profile_2026-11-06_02/summary.json",
    "icl_fixed_profile_split": "icl_satellite_fixed_profile_split_2026-11-06_02/summary.json",
    "wl_stack": "wl_random_shape_stack_lens_reuse_2026-11-04_01/summary.json",
    "wl_icl_relation": "wl_continuous_icl_relation_i2_2026-11-06_02/summary.json",
    "wl_region_influence": "wl_random_tract_region_influence_2026-11-06_01/summary.json",
    "wl_boost": "wl_source_density_boost_spatial_2026-11-05_06/summary.json",
    "wl_random_placebo": "wl_random_only_placebo_2026-11-06_01/summary.json",
    "icl_warp_variance": "icl_warp_pilot_coadd_variance_2026-11-06_02/summary.json",
    "icl_warp_variance_supported94": "icl_supported94_warp_noise_analysis_affine_03_2026-11-06/summary.json",
    "icl_warp_lanczos5": "icl_supported94_warp_noise_lanczos5_2026-11-06_02/summary.json",
    "icl_warp_mask_assoc94": "icl_supported94_warp_mask_association_2026-11-06_02/summary.json",
    "icl_annular_covariance": "icl_annular_coadd_covariance_2026-11-06_01/summary.json",
    "icl_annular_full783_covariance_64px": "icl_annular_coadd_covariance_full783_2026-11-06_02/summary.json",
    "icl_annular_full783_support_64px": "icl_annular_covariance_full783_support_audit_2026-11-06_02/summary.json",
    "icl_annular_full783_covariance_32px": "icl_annular_coadd_covariance_full783_32px_2026-11-06_01/summary.json",
    "icl_annular_full783_support_32px": "icl_annular_covariance_full783_32px_support_audit_2026-11-06_02/summary.json",
    "icl_random_variance_export": "icl_random_pixel_variance_results_2026-11-06_01/summary.json",
    "icl_random_variance_audit": "icl_random_pixel_variance_audit_2026-11-06_01/summary.json",
    "icl_random_variance_satellite_relation": "icl_random_point_variance_satellite_relation_2026-11-06_01/summary.json",
    "icl_wl_triple_max_t": "icl_satellite_wl_triple_overlap_maxt_2026-11-06_01/summary.json",
    "icl_wl_overlap_selection": "icl_wl_overlap_selection_i2_2026-11-06_03/summary.json",
    "icl_forced_coverage": "icl_satellite_forced_coverage_both_sides_2026-11-06_03/summary.json",
    "icl_forced_coverage_blocks": "icl_forced_coverage_selection_blocks_2026-11-06_02/summary.json",
    "icl_forced_coverage_adjusted": "icl_forced_coverage_selection_adjusted_2026-11-06_01/summary.json",
    "icl_forced_position_radial": "icl_forced_position_radial_relations_2026-11-06_01/summary.json",
    "icl_radial_outer_contrast": "icl_radial_outer_baseline_contrast_2026-11-06_01/summary.json",
    "icl_satellite_outer_profile_relation": "icl_satellite_outer_profile_relation_2026-11-06_01/summary.json",
    "s23b_coverage_audit": "icl_s23b_coverage_audit_2026-11-06_02/summary.json",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def load_sources(data_root: Path) -> tuple[dict, dict]:
    payloads, provenance = {}, {}
    for name, relative in SOURCES.items():
        path = data_root / relative
        if not path.is_file():
            raise FileNotFoundError(f"required Astrafek summary missing: {path}")
        payload = json.loads(path.read_text(encoding="utf-8"))
        if payload.get("science_measurement") is not False:
            raise ValueError(f"{name} has unexpected science gate: {payload.get('science_measurement')}")
        payloads[name] = payload
        provenance[str(path)] = sha256(path)
    shell = payloads["icl_shell_exclusion"]
    if (shell.get("status") != "complete"
            or shell.get("n_association_psf_supported_frames") != 48
            or shell.get("observational_validation_passed") is not False
            or shell.get("flux_covariance_available") is not False
            or any(row.get("psf_approved") is not False
                   for row in shell.get("filters", []))
            or {row.get("filter") for row in shell.get("filters", [])} != {"i", "i2"}):
        raise ValueError("hybrid-PSF shell-exclusion audit is incomplete or has unexpected guard flags")
    warp = payloads["icl_warp_variance"]
    if (warp.get("status") != "complete_exploratory"
            or warp.get("science_measurement") is not False
            or warp.get("hsc_network_used_by_analysis") is not False
            or warp.get("warp_archive", {}).get("exact_visit_set") is not True
            or warp.get("warp_archive", {}).get("n_visits") != 9
            or warp.get("inputs_sha256", {}).get("warp_archive")
            != "cdc9b642c5c5d0361b58a0bf7c0e733a60cb78caad5ed496ae00cac89b23953c"):
        raise ValueError("approved warp variance pilot is incomplete or has unexpected provenance/guard state")
    warp94 = payloads["icl_warp_variance_supported94"]
    if (warp94.get("schema") != "astrafek-icl-supported94-warp-noise-analysis-v1"
            or warp94.get("status") != "complete_with_data_gaps"
            or warp94.get("science_measurement") is not False
            or warp94.get("hsc_network_used") is not False
            or warp94.get("parent_plan_sha256")
            != "bd936a94bcdbff83a607acbfca66616c396db526ada50cc89c9e3e312d37007a"
            or warp94.get("parent_download_summary_sha256")
            != "0039321462f8158d1c006a1df64bb3e96ae988744ed780c845f74a60c79b95b4"
            or warp94.get("planned_pair_count") != 94
            or warp94.get("planned_center_count") != 188
            or warp94.get("analyzed_center_count") != 184
            or warp94.get("failed_center_count") != 4
            or warp94.get("radial_200_300_center_count") != 184):
        raise ValueError("supported-94 warp variance audit has unexpected scope or guard state")
    lanczos5 = payloads["icl_warp_lanczos5"]
    lanczos_dir = data_root / "icl_supported94_warp_noise_lanczos5_2026-11-06_02"
    lanczos_table = lanczos_dir / "per_center_lanczos5.tsv"
    if (lanczos5.get("schema") != "astrafek-icl-supported94-warp-noise-lanczos5-v1"
            or lanczos5.get("status") != "complete_with_data_gaps"
            or lanczos5.get("science_measurement") is not False
            or lanczos5.get("hsc_network_used") is not False
            or lanczos5.get("planned_pairs") != 94
            or lanczos5.get("planned_centers") != 188
            or lanczos5.get("analyzed_centers") != 184
            or lanczos5.get("failed_count") != 4
            or lanczos5.get("output_sha256", {}).get("per_center_lanczos5.tsv") != sha256(lanczos_table)):
        raise ValueError("supported-94 Lanczos-5 audit has unexpected scope, status, or provenance")
    lanczos_rows = read_tsv(lanczos_table)
    if len(lanczos_rows) != 184:
        raise ValueError("Lanczos-5 per-center table does not contain the 184 returned centers")
    lanczos5["validated_per_center_rows"] = lanczos_rows
    provenance[str(lanczos_table)] = sha256(lanczos_table)
    nearest_table = data_root / "icl_supported94_warp_noise_analysis_affine_03_2026-11-06/per_center_radial.tsv"
    nearest_rows = [row for row in read_tsv(nearest_table)
                    if row.get("inner_kpc") == "200" and row.get("outer_kpc") == "300"]
    if len(nearest_rows) != 184:
        raise ValueError("nearest-neighbor baseline table does not contain the expected annulus rows")
    lanczos5["validated_nearest_neighbor_rows"] = nearest_rows
    provenance[str(nearest_table)] = sha256(nearest_table)
    mask_assoc = payloads["icl_warp_mask_assoc94"]
    lanczos5 = payloads["icl_warp_lanczos5"]
    if (mask_assoc.get("schema") != "astrafek-icl-supported94-warp-mask-association-v1"
            or mask_assoc.get("status") != "complete_with_data_gaps"
            or mask_assoc.get("science_measurement") is not False
            or mask_assoc.get("hsc_network_used") is not False
            or mask_assoc.get("approved_plan_sha256")
            != "bd936a94bcdbff83a607acbfca66616c396db526ada50cc89c9e3e312d37007a"
            or mask_assoc.get("download_summary_sha256")
            != "0039321462f8158d1c006a1df64bb3e96ae988744ed780c845f74a60c79b95b4"
            or mask_assoc.get("planned_pairs") != 94
            or mask_assoc.get("planned_centers") != 188
            or mask_assoc.get("warp_complete_pair_count_before_mask_coverage_floor") != 90
            or mask_assoc.get("strict_pair_count_after_original_40pct_floor") != 81
            or mask_assoc.get("per_filter_coverage", {}).get("i", {}).get("n_strict_pairs_after_original_40pct_floor") != 23
            or mask_assoc.get("per_filter_coverage", {}).get("i2", {}).get("n_strict_pairs_after_original_40pct_floor") != 58):
        raise ValueError("supported-94 strict warp-mask association audit has unexpected scope or guard state")
    for band in ("i", "i2"):
        ci = mask_assoc["association_by_filter"][band]["paired_delta_ci95"]
        if not (ci[0] <= 0 <= ci[1]):
            raise ValueError(f"strict-mask paired interval unexpectedly excludes zero for {band}")
    annular = payloads["icl_annular_covariance"]
    annular_dir = data_root / "icl_annular_coadd_covariance_2026-11-06_01"
    annular_table = annular_dir / "paired_annular_covariance_comparison.tsv"
    if (annular.get("status") != "complete_exploratory"
            or annular.get("science_measurement") is not False
            or annular.get("hsc_network_used") is not False
            or annular.get("input_pair_count") != 94
            or annular.get("successful_center_annulus_count") != 264
            or annular.get("failed_inputs")
            or annular.get("outputs_sha256", {}).get(str(annular_table)) != sha256(annular_table)):
        raise ValueError("annular covariance audit has unexpected scope, status, or provenance")
    for band in ("i", "i2"):
        for metric in ("zero_lag", "sum8", "coverage"):
            for region in ("icl_200_300", "outer_350_500", "target_control_radial_difference_in_differences"):
                ci = annular["paired_block_bootstrap"][band][metric][region]["ci95"]
                if not ci[0] <= 0 <= ci[1]:
                    raise ValueError(f"unexpected annular target/control interval for {band}/{metric}/{region}")

    covariance_summaries = {
        64: payloads["icl_annular_full783_covariance_64px"],
        32: payloads["icl_annular_full783_covariance_32px"],
    }
    support_summaries = {
        64: payloads["icl_annular_full783_support_64px"],
        32: payloads["icl_annular_full783_support_32px"],
    }
    for patch_size in (64, 32):
        covariance = covariance_summaries[patch_size]
        support = support_summaries[patch_size]
        directory = data_root / (
            "icl_annular_coadd_covariance_full783_2026-11-06_02" if patch_size == 64
            else "icl_annular_coadd_covariance_full783_32px_2026-11-06_01")
        support_directory = data_root / (
            "icl_annular_covariance_full783_support_audit_2026-11-06_01" if patch_size == 64
            else "icl_annular_covariance_full783_32px_support_audit_2026-11-06_02")
        if (covariance.get("schema") != "astrafek-icl-full783-annular-coadd-covariance-v1"
                or covariance.get("status") != "complete_exploratory"
                or covariance.get("science_measurement") is not False
                or covariance.get("hsc_network_used") is not False
                or covariance.get("n_same_filter_pairs") != 783
                or covariance.get("n_image_records") != 1566
                or covariance.get("failed_inputs")
                or covariance.get("method", {}).get("patch_size_pixels") != patch_size):
            raise ValueError(f"full-783 {patch_size}px covariance audit has unexpected scope or status")
        for name in ("per_center_icl_annulus_covariance.tsv", "paired_icl_annulus_covariance.tsv"):
            table = directory / name
            if covariance.get("outputs_sha256", {}).get(str(table)) != sha256(table):
                raise ValueError(f"full-783 {patch_size}px covariance table hash mismatch: {table}")
            provenance[str(table)] = sha256(table)
        if (support.get("schema") != "astrafek-icl-full783-annular-covariance-support-audit-v1"
                or support.get("status") != "complete_exploratory"
                or support.get("science_measurement") is not False
                or support.get("hsc_network_used") is not False
                or support.get("image_pixels_read") is not False
                or support.get("method", {}).get("patch_size_pixels") != patch_size
                or support.get("input", {}).get("summary_sha256") != sha256(directory / "summary.json")):
            raise ValueError(f"full-783 {patch_size}px support audit has unexpected provenance or scope")
        for name, digest in support.get("outputs_sha256", {}).items():
            table = Path(name)
            if not table.is_file() or sha256(table) != digest:
                raise ValueError(f"full-783 {patch_size}px support table hash mismatch: {table}")
            provenance[str(table)] = digest
    for patch_size in (64, 32):
        for band, expected_pairs in (("i", 302), ("i2", 481)):
            row = next(item for item in support_summaries[patch_size]["support_by_filter"]
                       if item["filter"] == band)
            if row.get("n_same_filter_pairs") != expected_pairs:
                raise ValueError(f"full-783 {patch_size}px {band}-band pair count changed")

    variance_export = payloads["icl_random_variance_export"]
    variance_audit = payloads["icl_random_variance_audit"]
    variance_results_dir = data_root / "icl_random_pixel_variance_results_2026-11-06_01"
    variance_audit_dir = data_root / "icl_random_pixel_variance_audit_2026-11-06_01"
    variance_status_path = variance_results_dir / "status.json"
    variance_center_table = variance_audit_dir / "center_random_point_summary.tsv"
    variance_pair_table = variance_audit_dir / "matched_target_control_summary.tsv"
    variance_status = json.loads(variance_status_path.read_text(encoding="utf-8"))
    if (variance_export.get("schema") != "astrafek-icl-random-variance-summary-v1"
            or variance_export.get("status") != "complete_exploratory"
            or variance_export.get("plan_sha256") != "eeef98ef2210337f714e36d8a4afcd352fe312a536c90e8e236ff9f843f53e52"
            or variance_export.get("n_requested_centers") != 1598
            or variance_export.get("n_returned_random_points") != 40048
            or variance_export.get("images_returned") is not False
            or variance_export.get("shape_rows_returned") is not False
            or variance_export.get("photoz_rows_returned") is not False
            or variance_export.get("object_catalog_rows_returned") is not False
            or variance_export.get("credentials_saved") is not False
            or variance_status.get("status") != "complete"
            or variance_status.get("query_submitted") is not True
            or variance_status.get("plan_sha256") != variance_export.get("plan_sha256")
            or len(variance_status.get("jobs", {})) != 16
            or any(job.get("status") != "complete" for job in variance_status.get("jobs", {}).values())):
        raise ValueError("approved random-point variance export has unexpected scope or incomplete batches")
    if (variance_audit.get("schema") != "astrafek-icl-random-point-variance-audit-v1"
            or variance_audit.get("status") != "complete_exploratory_audit"
            or variance_audit.get("science_measurement") is not False
            or variance_audit.get("plan_sha256") != variance_export.get("plan_sha256")
            or variance_audit.get("n_returned_random_points") != variance_export.get("n_returned_random_points")
            or variance_audit.get("n_planned_centers") != 1598
            or variance_audit.get("n_centers_with_points") != 1587
            or variance_audit.get("n_target_centers_with_points") != 791
            or variance_audit.get("n_control_centers_with_points") != 796
            or variance_audit.get("n_matched_target_control_pairs") != 798
            or variance_audit.get("matched_pairs_with_both_endpoint_statistics", {}).get("i_pix_variance") != 787
            or variance_audit.get("input_sha256", {}).get("status") != sha256(variance_status_path)
            or variance_audit.get("outputs", {}).get(variance_center_table.name) != sha256(variance_center_table)
            or variance_audit.get("outputs", {}).get(variance_pair_table.name) != sha256(variance_pair_table)):
        raise ValueError("random-point variance audit has unexpected support, status, or provenance")
    variance_audit["validated_center_rows"] = read_tsv(variance_center_table)
    variance_audit["validated_pair_rows"] = read_tsv(variance_pair_table)
    point_counts = sorted(int(row["n_random_points"]) for row in variance_audit["validated_center_rows"])
    variance_audit["random_points_per_center_q16_median_q84"] = [
        point_counts[int(q * (len(point_counts) - 1))] for q in (0.16, 0.50, 0.84)
    ]
    provenance[str(variance_status_path)] = sha256(variance_status_path)
    provenance[str(variance_center_table)] = sha256(variance_center_table)
    provenance[str(variance_pair_table)] = sha256(variance_pair_table)
    for name, job in sorted(variance_status["jobs"].items()):
        csv_path = variance_results_dir / f"{name}.csv"
        if not csv_path.is_file() or sha256(csv_path) != job.get("sha256"):
            raise ValueError(f"random-point variance batch hash mismatch: {csv_path}")
        provenance[str(csv_path)] = job["sha256"]

    variance_relation = payloads["icl_random_variance_satellite_relation"]
    variance_relation_dir = data_root / "icl_random_point_variance_satellite_relation_2026-11-06_01"
    variance_relation_table = variance_relation_dir / "joined_supported_pairs.tsv"
    if (variance_relation.get("schema") != "astrafek-icl-random-point-variance-satellite-relation-v1"
            or variance_relation.get("status") != "complete_exploratory_sensitivity"
            or variance_relation.get("science_measurement") is not False
            or variance_relation.get("hsc_network_used") is not False
            or variance_relation.get("image_pixels_read") is not False
            or variance_relation.get("data_requested") is not False
            or variance_relation.get("plan_sha256") != variance_export.get("plan_sha256")
            or variance_relation.get("n_pairs") != 94
            or variance_relation.get("filter_counts") != {"i": 26, "i2": 68}
            or variance_relation.get("paired_block_resampling", {}).get("replicas_valid_joint") != 5000
            or len(variance_relation.get("random_point_features", [])) != 8
            or variance_relation.get("joined_pairs_sha256") != sha256(variance_relation_table)):
        raise ValueError("random-point ICL/satellite sensitivity has unexpected scope or provenance")
    for path_text, digest in variance_relation.get("inputs_sha256", {}).items():
        if path_text.startswith("/"):
            path = Path(path_text)
            if not path.is_file() or sha256(path) != digest:
                raise ValueError(f"random-point ICL/satellite input hash mismatch: {path}")
    for band in ("i", "i2"):
        for metric in variance_relation["random_point_features"]:
            interval = variance_relation["results"][band][metric]["paired_change_maxT_ci95"]
            if interval[0] > 0 or interval[1] < 0:
                raise ValueError(f"unexpected nonzero random-point adjustment interval: {band}/{metric}")
    provenance[str(variance_relation_table)] = sha256(variance_relation_table)

    triple = payloads["icl_wl_triple_max_t"]
    triple_dir = data_root / "icl_satellite_wl_triple_overlap_maxt_2026-11-06_01"
    triple_table = triple_dir / "max_t_metrics.tsv"
    triple_rows = read_tsv(triple_table)
    if (triple.get("schema") != "astrafek-icl-satellite-wl-triple-overlap-max-t-v1"
            or triple.get("status") != "complete_exploratory"
            or triple.get("science_measurement") is not False
            or triple.get("network_used") is not False
            or triple.get("family_size") != 29
            or triple.get("n_triple_overlap_lenses") != 27
            or triple.get("replicates_valid") != 10000
            or len(triple_rows) != 29
            or triple.get("output_sha256", {}).get(triple_table.name) != sha256(triple_table)):
        raise ValueError("triple-overlap max-T audit has unexpected scope, status, or provenance")
    triple["validated_metrics"] = triple_rows
    provenance[str(annular_table)] = sha256(annular_table)
    provenance[str(triple_table)] = sha256(triple_table)

    selection = payloads["icl_wl_overlap_selection"]
    selection_dir = data_root / "icl_wl_overlap_selection_i2_2026-11-06_03"
    selection_table = selection_dir / "selection_balance.tsv"
    block_table = selection_dir / "block_support.tsv"
    selection_rows = read_tsv(selection_table)
    block_rows = read_tsv(block_table)
    if (selection.get("schema") != "astrafek-icl-wl-overlap-selection-audit-v1"
            or selection.get("status") != "complete_exploratory"
            or selection.get("science_measurement") is not False
            or selection.get("network_used") is not False
            or selection.get("parent_i2_satellite_supported_n") != 68
            or selection.get("common_wl_supported_n") != 27
            or selection.get("spatial_blocks") != 33
            or selection.get("mixed_selection_blocks") != 9
            or len(selection_rows) != 4 or len(block_rows) != 33
            or selection.get("output_sha256", {}).get(selection_table.name) != sha256(selection_table)
            or selection.get("output_sha256", {}).get(block_table.name) != sha256(block_table)):
        raise ValueError("WL overlap-selection audit has unexpected scope, status, or provenance")
    mixed_blocks = [r for r in block_rows
                    if int(r["n_common_wl"]) > 0 and int(r["n_without_common_wl"]) > 0]
    selection["mixed_block_selected_n"] = sum(int(r["n_common_wl"]) for r in mixed_blocks)
    selection["mixed_block_excluded_n"] = sum(int(r["n_without_common_wl"]) for r in mixed_blocks)
    selection["familywise_maxT_p_min"] = min(float(r["familywise_maxT_p_4"]) for r in selection_rows)
    if (selection["mixed_block_selected_n"] != 10
            or selection["mixed_block_excluded_n"] != 16):
        raise ValueError("WL overlap-selection spatial positivity counts changed")
    selection["validated_metrics"] = selection_rows
    provenance[str(selection_table)] = sha256(selection_table)
    provenance[str(block_table)] = sha256(block_table)

    forced_coverage = payloads["icl_forced_coverage"]
    forced_dir = data_root / "icl_satellite_forced_coverage_both_sides_2026-11-06_03"
    forced_table = forced_dir / "paired_cluster_comparison.tsv"
    forced_rows = read_tsv(forced_table)
    if (forced_coverage.get("schema") != "astrafek-icl-satellite-forced-coverage-paired-comparison-v1"
            or forced_coverage.get("status") != "complete_exploratory"
            or forced_coverage.get("science_measurement") is not False
            or forced_coverage.get("replicas_requested") != 10000
            or forced_coverage.get("scope") != "existing HSC cutouts only; target and rotated-control annular coverage floors both 20 percent"
            or len(forced_rows) != 224
            or forced_coverage.get("output_sha256", {}).get(forced_table.name) != sha256(forced_table)):
        raise ValueError("20/20 forced-profile coverage comparison has unexpected scope or provenance")
    for band, expected in (("i", (71, 50, 21)), ("i2", (153, 110, 43))):
        support_row = forced_coverage["selection_support"][band]
        association = forced_coverage["association_by_filter"][band]
        observed = (support_row.get("icl_supported_pairs"),
                    support_row.get("20_20_selected_pairs"),
                    support_row.get("excluded_pairs"))
        if observed != expected or association.get("replicas_valid", 0) < 9000:
            raise ValueError(f"20/20 forced-profile support changed for {band}: {observed}")
        for key in ("paired_change_sample_expansion_full_minus_common",
                    "paired_change_added_forced_measurements_common_minus_detected",
                    "paired_change_total_full_minus_detected"):
            interval = association[key]["ci95"]
            if not interval[0] <= 0 <= interval[1]:
                raise ValueError(f"20/20 paired interval unexpectedly excludes zero: {band}/{key}")
    for path_text, digest in forced_coverage.get("inputs_sha256", {}).items():
        path = Path(path_text)
        if not path.is_file() or sha256(path) != digest:
            raise ValueError(f"20/20 forced-profile input hash mismatch: {path}")
        provenance[str(path)] = digest
    provenance[str(forced_table)] = sha256(forced_table)

    block_selection = payloads["icl_forced_coverage_blocks"]
    block_selection_dir = data_root / "icl_forced_coverage_selection_blocks_2026-11-06_02"
    block_support_table = block_selection_dir / "spatial_block_support.tsv"
    mixed_contrast_table = block_selection_dir / "mixed_block_icl_contrasts.tsv"
    if (block_selection.get("schema") != "astrafek-icl-forced-coverage-spatial-selection-audit-v1"
            or block_selection.get("status") != "complete_exploratory"
            or block_selection.get("science_measurement") is not False
            or block_selection.get("network_used") is not False
            or block_selection.get("image_pixels_read") is not False
            or block_selection.get("inputs_sha256", {}).get(str(forced_table)) != sha256(forced_table)
            or block_selection.get("output_sha256", {}).get(block_support_table.name) != sha256(block_support_table)
            or block_selection.get("output_sha256", {}).get(mixed_contrast_table.name) != sha256(mixed_contrast_table)):
        raise ValueError("20/20 spatial selection audit has unexpected scope or provenance")
    expected_block_selection = {
        "i": (18, 9, 40, 19),
        "i2": (40, 21, 77, 40),
    }
    for band, expected in expected_block_selection.items():
        row = block_selection["by_filter"][band]
        observed = (row.get("connected_paired_sky_blocks"), row.get("mixed_selection_blocks"),
                    row.get("selected_pairs_in_mixed_blocks"), row.get("excluded_pairs_in_mixed_blocks"))
        if observed != expected or row.get("block_bootstrap_replicates") != 10000:
            raise ValueError(f"20/20 spatial overlap counts changed for {band}: {observed}")
    for path_text, digest in block_selection.get("inputs_sha256", {}).items():
        path = Path(path_text)
        if not path.is_file() or sha256(path) != digest:
            raise ValueError(f"20/20 spatial selection input hash mismatch: {path}")
        provenance[str(path)] = digest
    for table in (block_support_table, mixed_contrast_table):
        provenance[str(table)] = sha256(table)

    adjusted_selection = payloads["icl_forced_coverage_adjusted"]
    adjusted_dir = data_root / "icl_forced_coverage_selection_adjusted_2026-11-06_01"
    adjusted_table = adjusted_dir / "mixed_block_model_support.tsv"
    if (adjusted_selection.get("schema") != "astrafek-icl-forced-coverage-adjusted-selection-v1"
            or adjusted_selection.get("status") != "complete_exploratory"
            or adjusted_selection.get("science_measurement") is not False
            or adjusted_selection.get("network_used") is not False
            or adjusted_selection.get("image_pixels_read") is not False
            or adjusted_selection.get("output_sha256", {}).get(adjusted_table.name) != sha256(adjusted_table)):
        raise ValueError("adjusted 20/20 selection analysis has unexpected scope or provenance")
    for band, expected in (("i", (9, 59, 40, 19)), ("i2", (21, 117, 77, 40))):
        row = adjusted_selection["by_filter"][band]
        observed = (row.get("n_mixed_blocks"), row.get("n_pairs_in_mixed_blocks"),
                    row.get("n_selected_in_mixed_blocks"), row.get("n_excluded_in_mixed_blocks"))
        ci = row.get("block_bootstrap_ci95", [])
        if (observed != expected or row.get("block_bootstrap_replicates_valid") != 10000
                or len(ci) != 2 or ci[0] > 0 or ci[1] < 0):
            raise ValueError(f"adjusted 20/20 selection support or interval changed for {band}: {observed}")
    for path_text, digest in adjusted_selection.get("inputs_sha256", {}).items():
        path = Path(path_text)
        if not path.is_file() or sha256(path) != digest:
            raise ValueError(f"adjusted 20/20 selection input hash mismatch: {path}")
        provenance[str(path)] = digest
    provenance[str(adjusted_table)] = sha256(adjusted_table)

    forced_radial = payloads["icl_forced_position_radial"]
    forced_radial_dir = data_root / "icl_forced_position_radial_relations_2026-11-06_01"
    forced_radial_table = forced_radial_dir / "icl_forced_position_radial_relations.tsv"
    forced_radial_support = forced_radial_dir / "radial_support.tsv"
    forced_radial_max_t = forced_radial_dir / "maxT_null.tsv"
    forced_radial_runner = Path(__file__).resolve().parent / "analyze_icl_forced_position_radial_relations.py"
    forced_radial_rows = read_tsv(forced_radial_table)
    if (forced_radial.get("schema") != "astrafek-icl-forced-position-radial-relations-v1"
            or forced_radial.get("status") != "complete_exploratory"
            or forced_radial.get("science_measurement") is not False
            or forced_radial.get("hsc_network_used") is not False
            or forced_radial.get("image_pixels_read") is not False
            or forced_radial.get("n_icl_supported_pairs") != 224
            or forced_radial.get("n_relation_tests") != 9
            or forced_radial.get("replicates") != 10000
            or forced_radial.get("runner_sha256") != sha256(forced_radial_runner)
            or len(forced_radial_rows) != 9):
        raise ValueError("forced-position radial analysis has unexpected scope or status")
    for name, path in ((forced_radial_table.name, forced_radial_table),
                       (forced_radial_support.name, forced_radial_support),
                       (forced_radial_max_t.name, forced_radial_max_t)):
        if forced_radial.get("output_sha256", {}).get(name) != sha256(path):
            raise ValueError(f"forced-position radial output hash mismatch: {path}")
        provenance[str(path)] = sha256(path)
    provenance[str(forced_radial_runner)] = sha256(forced_radial_runner)
    expected_forced_support = {
        ("i", "20-30"): (16, 9), ("i", "30-40"): (18, 9),
        ("i", "40-50"): (25, 12), ("i", "50-75"): (34, 12),
        ("i2", "15-20"): (14, 12), ("i2", "20-30"): (27, 19),
        ("i2", "30-40"): (43, 29), ("i2", "40-50"): (53, 29),
        ("i2", "50-75"): (78, 35),
    }
    observed_forced_support = {
        (row["filter"], row["radial_bin_kpc"]):
            (int(row["n_clusters"]), int(row["n_connected_5deg_blocks"]))
        for row in forced_radial_rows
    }
    if observed_forced_support != expected_forced_support:
        raise ValueError("forced-position radial support counts changed")
    if any(float(row["p_maxT_across_supported_bins"]) < 0.05 for row in forced_radial_rows):
        raise ValueError("forced-position radial max-T result unexpectedly crosses 0.05")
    for path_text, digest in forced_radial.get("input_sha256", {}).items():
        path = Path(path_text)
        if not path.is_file() or sha256(path) != digest:
            raise ValueError(f"forced-position radial input hash mismatch: {path}")
        provenance[str(path)] = digest
    forced_radial["validated_results"] = forced_radial_rows

    outer_contrast = payloads["icl_radial_outer_contrast"]
    outer_relation = payloads["icl_satellite_outer_profile_relation"]
    outer_dir = data_root / "icl_radial_outer_baseline_contrast_2026-11-06_01"
    outer_table = outer_dir / "outer_baseline_contrast.tsv"
    outer_draws = outer_dir / "block_bootstrap_outer_contrast_draws.npz"
    outer_runner = Path(__file__).resolve().parent / "analyze_icl_radial_outer_baseline_contrast.py"
    outer_rows = read_tsv(outer_table)
    if (outer_contrast.get("schema") != "astrafek-icl-radial-outer-baseline-contrast-v1"
            or outer_contrast.get("status") != "complete_exploratory"
            or outer_contrast.get("science_measurement") is not False
            or outer_contrast.get("network_used") is not False
            or outer_contrast.get("image_pixels_read") is not False
            or outer_contrast.get("bootstrap", {}).get("replicas") != 20000
            or outer_contrast.get("support", {}).get("i", {}).get("n_pairs") != 20
            or outer_contrast.get("support", {}).get("i", {}).get("n_connected_blocks") != 10
            or outer_contrast.get("support", {}).get("i2", {}).get("n_pairs") != 48
            or outer_contrast.get("support", {}).get("i2", {}).get("n_connected_blocks") != 30
            or len(outer_rows) != 4
            or outer_contrast.get("result_files_sha256", {}).get(str(outer_table)) != sha256(outer_table)
            or outer_contrast.get("result_files_sha256", {}).get(str(outer_draws)) != sha256(outer_draws)):
        raise ValueError("radial outer-baseline contrast has unexpected scope or provenance")
    expected_outer = {
        ("i", "raw"): (20, 10), ("i", "corrected"): (20, 10),
        ("i2", "raw"): (48, 30), ("i2", "corrected"): (48, 30),
    }
    observed_outer = {(row["filter"], row["extraction"]):
                      (int(row["n_pairs"]), int(row["n_connected_blocks"])) for row in outer_rows}
    if observed_outer != expected_outer or any(
            row.get("interpretation") != "exploratory radial contrast; not a pure ICL estimate"
            for row in outer_rows):
        raise ValueError("radial outer-baseline result support or interpretation changed")
    for band in ("i", "i2"):
        for mode in ("raw", "corrected"):
            row = next(item for item in outer_rows if item["filter"] == band and item["extraction"] == mode)
            lo, hi = float(row["block_bootstrap_ci95_low"]), float(row["block_bootstrap_ci95_high"])
            if band == "i" and not (lo > 0 and hi > 0):
                raise ValueError("i-band radial contrast interval changed")
            if band == "i2" and not (lo < 0 < hi):
                raise ValueError("i2 radial contrast interval changed")
    for path_text, digest in outer_contrast.get("source_hashes", {}).items():
        path = Path(path_text)
        if not path.is_file() or sha256(path) != digest:
            raise ValueError(f"radial contrast source hash mismatch: {path}")
        provenance[str(path)] = digest

    outer_relation = payloads["icl_satellite_outer_profile_relation"]
    relation_dir = data_root / "icl_satellite_outer_profile_relation_2026-11-06_01"
    relation_table = relation_dir / "association_results.tsv"
    relation_pairs = relation_dir / "icl_satellite_common_support_pairs.tsv"
    relation_bootstrap = relation_dir / "block_bootstrap_partial_correlations.npz"
    relation_runner = Path(__file__).resolve().parent / "analyze_icl_satellite_outer_profile_relation.py"
    relation_rows = read_tsv(relation_table)
    expected_relation_support = {
        ("i", "raw"): (20, 10), ("i", "corrected"): (20, 10),
        ("i2", "raw"): (47, 30), ("i2", "corrected"): (47, 30),
    }
    observed_relation_support = {
        (row["filter"], row["icl_extraction"]):
        (int(row["n_clusters"]), int(row["n_connected_blocks"])) for row in relation_rows
    }
    if (outer_relation.get("schema") != "astrafek-icl-satellite-outer-profile-relation-v1"
            or outer_relation.get("status") != "complete_exploratory"
            or outer_relation.get("science_measurement") is not False
            or outer_relation.get("network_used") is not False
            or outer_relation.get("image_pixels_read") is not False
            or outer_relation.get("bootstrap", {}).get("replicas") != 20000
            or outer_relation.get("covariates") != ["cluster redshift", "log CAMIRA richness"]
            or observed_relation_support != expected_relation_support
            or outer_relation.get("output_sha256", {}).get(str(relation_table)) != sha256(relation_table)
            or outer_relation.get("output_sha256", {}).get(str(relation_pairs)) != sha256(relation_pairs)
            or outer_relation.get("output_sha256", {}).get(str(relation_bootstrap)) != sha256(relation_bootstrap)):
        raise ValueError("ICL/satellite outer-profile relation has unexpected scope or provenance")
    if any(float(row["block_bootstrap_ci95_low"]) >= 0
           or float(row["block_bootstrap_ci95_high"]) <= 0 for row in relation_rows):
        raise ValueError("ICL/satellite outer-profile intervals unexpectedly exclude zero")
    for path_text, digest in outer_relation.get("source_hashes", {}).items():
        path = Path(path_text)
        if not path.is_file() or sha256(path) != digest:
            raise ValueError(f"ICL/satellite outer-profile source hash mismatch: {path}")
        provenance[str(path)] = digest
    for path in (relation_table, relation_pairs, relation_bootstrap, relation_runner):
        provenance[str(path)] = sha256(path)
    outer_relation["validated_results"] = relation_rows
    provenance[str(outer_runner)] = sha256(outer_runner)
    provenance[str(outer_table)] = sha256(outer_table)
    provenance[str(outer_draws)] = sha256(outer_draws)
    outer_contrast["validated_results"] = outer_rows

    s23b_coverage = payloads["s23b_coverage_audit"]
    s23b_audit_dir = data_root / "icl_s23b_coverage_audit_2026-11-06_02"
    s23b_candidate_table = s23b_audit_dir / "per_candidate_coverage.tsv"
    s23b_coverage_runner = Path(__file__).resolve().parent / "audit_s23b_icl_coverage_metadata.py"
    if (s23b_coverage.get("schema") != "astrafek-s23b-icl-coverage-metadata-audit-v1"
            or s23b_coverage.get("status") != "complete_exploratory_audit"
            or s23b_coverage.get("science_measurement") is not False
            or s23b_coverage.get("hsc_network_used_by_audit") is not False
            or s23b_coverage.get("image_pixels_read") is not False
            or s23b_coverage.get("plan_sha256") != "62c28446f553bf45a940ea0126af1521d4ee0c1356c12cc73f42712df482c0c7"
            or s23b_coverage.get("n_candidates") != 1333
            or s23b_coverage.get("n_probes") != 11997
            or s23b_coverage.get("n_returned_patch_metadata_rows") != 13575
            or s23b_coverage.get("n_candidates_all_9_probes_covered") != 1043
            or s23b_coverage.get("n_candidates_all_9_probes_with_imeas_exists") != 1043
            or s23b_coverage.get("n_candidates_all_9_probes_with_idet_bkgd_exists") != 1043
            or s23b_coverage.get("runner_sha256") != sha256(s23b_coverage_runner)
            or s23b_coverage.get("output_sha256", {}).get(s23b_candidate_table.name) != sha256(s23b_candidate_table)):
        raise ValueError("S23B footprint-metadata audit has unexpected scope, status, or hashes")
    for path_text, digest in s23b_coverage.get("input_sha256", {}).items():
        path = Path(path_text)
        if not path.is_file() or sha256(path) != digest:
            raise ValueError(f"S23B footprint audit input hash mismatch: {path}")
        provenance[str(path)] = digest
    for path_text, digest in s23b_coverage.get("batch_sha256", {}).items():
        path = Path(path_text)
        if not path.is_file() or sha256(path) != digest:
            raise ValueError(f"S23B CAS result batch hash mismatch: {path}")
        provenance[str(path)] = digest
    provenance[str(s23b_candidate_table)] = sha256(s23b_candidate_table)
    provenance[str(s23b_coverage_runner)] = sha256(s23b_coverage_runner)
    s23b_coverage["validated_candidate_table"] = read_tsv(s23b_candidate_table)
    return payloads, provenance


def fmt(value: float, digits: int = 3) -> str:
    return f"{float(value):.{digits}f}"


def ci2(values: list[float]) -> str:
    return f"[{fmt(values[0])}, {fmt(values[-1])}]"


def ci5(values: list[float]) -> str:
    return f"[{fmt(values[0], 5)}, {fmt(values[-1], 5)}]"


def render(payloads: dict, provenance: dict, script_path: Path,
           query_audit: dict, s23b_query_audit: dict) -> str:
    support = payloads["icl_support94"]["results_by_filter"]
    pixel = payloads["icl_pixelnoise94"]["filters"]
    ipw = payloads["icl_ipw_pixelnoise"]["results_by_model"]["primary_support_covariates"]["association_by_filter"]
    member = {row["filter"]: row for row in payloads["icl_membercount"]["results"]}
    position = {row["filter"]: row for row in payloads["icl_positionmix"]["results"]}
    aperture = payloads["icl_aperture_covariance"]["filters"]
    cross = payloads["icl_injection_crosstalk"]
    grid = payloads["icl_transfer_grid"]
    shell = payloads["icl_shell_exclusion"]
    warp = payloads["icl_warp_variance"]
    warp94 = payloads["icl_warp_variance_supported94"]
    mask_assoc = payloads["icl_warp_mask_assoc94"]
    warp_annulus = next(
        row for row in warp["coadd_variance_comparison"]["radial_bins"]
        if row["inner_kpc"] == 200 and row["outer_kpc"] == 300)
    frame = payloads["frame_metadata"]
    frame_satellite = payloads["frame_metadata_satellite"]
    nonlinear = payloads["icl_nonlinear_adjustment"]
    spatial_cv = payloads["icl_spatial_cv"]["filter_results"]
    fixed_profile = payloads["icl_fixed_profile"]
    fixed_split = payloads["icl_fixed_profile_split"]
    wl = payloads["wl_stack"]
    wl_icl = payloads["wl_icl_relation"]
    regions = payloads["wl_region_influence"]
    boost = payloads["wl_boost"]
    random_placebo = payloads["wl_random_placebo"]
    annular = payloads["icl_annular_covariance"]
    full_annular = {
        patch_size: payloads[f"icl_annular_full783_covariance_{patch_size}px"]
        for patch_size in (64, 32)
    }
    full_support = {
        patch_size: payloads[f"icl_annular_full783_support_{patch_size}px"]
        for patch_size in (64, 32)
    }
    random_variance = payloads["icl_random_variance_audit"]
    random_variance_relation = payloads["icl_random_variance_satellite_relation"]
    triple_rows = payloads["icl_wl_triple_max_t"]["validated_metrics"]
    selection = payloads["icl_wl_overlap_selection"]
    forced_coverage = payloads["icl_forced_coverage"]
    block_selection = payloads["icl_forced_coverage_blocks"]
    adjusted_selection = payloads["icl_forced_coverage_adjusted"]
    forced_radial = payloads["icl_forced_position_radial"]
    outer_contrast = payloads["icl_radial_outer_contrast"]
    outer_relation = payloads["icl_satellite_outer_profile_relation"]
    s23b_coverage = payloads["s23b_coverage_audit"]
    triple_by_key = {(r["relationship"], r["component"], r["radial_bin"]): r for r in triple_rows}
    cross_satellite = triple_by_key[("wl_satellite", "cross", "2")]
    tangential_satellite = triple_by_key[("wl_satellite", "tangential", "6")]
    min_triple_p = min(float(r["wild_cluster_p_max_t_29"]) for r in triple_rows)
    attempts_by_host = {}
    for row in query_audit["attempts"]:
        attempts_by_host.setdefault(row["api_host"], []).append(row["http_status"])
    attempts_text = "; ".join(
        f"{host}: HTTP {status} on {statuses.count(status)} attempt(s)"
        for host, statuses in sorted(attempts_by_host.items())
        for status in sorted(set(statuses)))
    frame_satellite_pmin = min(
        result["maxT_adjusted_p"]
        for band in ("i", "i2")
        for result in frame_satellite["results"][band].values())
    frame_satellite_i2_r = [
        result["metadata_adjusted_partial_r"]
        for result in frame_satellite["results"]["i2"].values()]
    nonlinear_pmin = min(
        result["maxT_adjusted_p"]
        for band in ("i", "i2")
        for result in nonlinear["results"][band]["models"].values())
    nonlinear_i2_r = [result["adjusted_partial_r"]
                      for result in nonlinear["results"]["i2"]["models"].values()]

    lines = [
        "# ICL and weak-lensing evidence status (2026-11-06)",
        "",
        "This report consolidates current Astrafek diagnostics for the HSC cluster",
        "ICL–satellite relation and background-galaxy weak lensing. It reads saved",
        "JSON summaries and audit tables only. Building this report made no",
        "HSC request, downloaded no data, and read no image pixels.",
        "Every source summary has `science_measurement=false`; this report does not",
        "promote exploratory results to a science measurement.",
        "",
        "## ICL–satellite relation",
        "",
        "The parent is the selected 800-cluster HSC-I sample. Joint ICL 200–300 kpc",
        "and satellite 50–75 kpc support exists for 94 cluster/filter pairs",
        "(26 `i`, 68 `i2`). The primary statistic is the partial correlation",
        "controlling redshift and CAMIRA richness.",
        "",
        "| Filter | N pairs | Raw partial r (95% block interval) | Finite-lag corrected r (95% interval) | IPW + noise corrected r (95% interval) | Effective N |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for band in ("i", "i2"):
        raw = support[band]["support94"]
        noise = pixel[band]
        weighted = ipw[band]
        raw_ci = raw["block_bootstrap_ci95"]
        noise_ci = noise["noise_corrected_partial_r_bootstrap_2p5_50_97p5"]
        lines.append(
            f"| `{band}` | {raw['n_pairs']} | {fmt(raw['partial_r_z_richness'])} "
            f"{ci2(raw_ci)} | {fmt(noise['point_estimate']['pixel_noise_corrected_partial_r'])} "
            f"{ci2(noise_ci)} | {fmt(weighted['pixel_noise_corrected_partial_r'])} "
            f"{ci2(weighted['pixel_noise_corrected_ipw_partial_r_ci95'])} | "
            f"{fmt(weighted['effective_n'])} / {weighted['n_selected_rows']} |"
        )
    member_change = ci2(ast.literal_eval(member["i2"]["member_count_adjustment_delta_ci95"]))
    position_change = ci2(ast.literal_eval(position["i2"]["position_adjustment_delta_ci95"]))
    lines.extend([
        "",
        "The `i2` raw and finite-lag point estimates are positive, but the estimate",
        "moves to 0.238 after measured-covariate support weighting plus the same",
        "finite-lag noise term; its interval includes zero. The IPW result assumes",
        "missing satellite profiles are conditionally independent of their",
        "unobserved signal given modeled covariates, which is unverified.",
        "",
        "Two additional composition adjustments leave the result unresolved. In",
        f"`i2`, controlling contributing-member count changes r from {fmt(member['i2']['icl_vs_satellite_baseline_partial_r'])} "
        f"to {fmt(member['i2']['icl_vs_satellite_adjusted_for_member_count_partial_r'])} "
        f"(paired-change 95% interval {member_change}). Adding the",
        f"median host-centric radius changes r to {fmt(position['i2']['icl_vs_satellite_adjusted_for_position_mix_partial_r'])} "
        f"(paired-change 95% interval {position_change}). Both paired",
        "change intervals include zero. These are selected-sample sensitivities,",
        "not completeness or causal corrections.",
        "A separate post-hoc sensitivity lowers both target and rotated-control",
        "annular-coverage floors to 20% and adds forced measurements at unmatched",
        "CAMIRA positions. It has support for 50/71 i and 110/153 i2 ICL-supported",
        "pairs. For this sensitivity, partial r controls redshift and log richness.",
        f"Selected pairs have higher median 200–300 kpc ICL residual than excluded pairs (i: {fmt(forced_coverage['selection_support']['i']['selected_median_icl_delta_200_300'])} vs {fmt(forced_coverage['selection_support']['i']['excluded_median_icl_delta_200_300'])} nJy arcsec⁻²; i2: {fmt(forced_coverage['selection_support']['i2']['selected_median_icl_delta_200_300'])} vs {fmt(forced_coverage['selection_support']['i2']['excluded_median_icl_delta_200_300'])}). Support is related to the ICL outcome, so the selected subset cannot be treated as representative.",
        f"The full 20/20 partial r is {fmt(forced_coverage['association_by_filter']['i']['partial_r_20_20_full'])} in i and {fmt(forced_coverage['association_by_filter']['i2']['partial_r_20_20_full'])} in i2. The paired change from detected-only common support to the full 20/20 sample is {fmt(forced_coverage['association_by_filter']['i']['paired_change_total_full_minus_detected']['median'])} {ci2(forced_coverage['association_by_filter']['i']['paired_change_total_full_minus_detected']['ci95'])} in i and {fmt(forced_coverage['association_by_filter']['i2']['paired_change_total_full_minus_detected']['median'])} {ci2(forced_coverage['association_by_filter']['i2']['paired_change_total_full_minus_detected']['ci95'])} in i2; both intervals include zero.",
        "These forced CAMIRA positions are not confirmed isolated galaxies; this",
        "does not estimate a population correction or establish a physical relation.",
        f"Within connected 5-degree blocks, 9/18 i and 21/40 i2 blocks contain both selected and excluded pairs. The equal-block selected-minus-excluded ICL residual contrast is {fmt(block_selection['by_filter']['i']['mixed_block_selected_minus_excluded_icl_delta_mean'])} nJy arcsec⁻² in i (95% block-bootstrap interval {ci2(block_selection['by_filter']['i']['equal_block_mean_difference_bootstrap_ci95'])}) and {fmt(block_selection['by_filter']['i2']['mixed_block_selected_minus_excluded_icl_delta_mean'])} in i2 ({ci2(block_selection['by_filter']['i2']['equal_block_mean_difference_bootstrap_ci95'])}). The i2 contrast remains positive in this descriptive mixed-block analysis; i is unresolved. It is unadjusted for redshift/richness and does not validate transport to single-status blocks.",
        f"Adding redshift and log richness in a within-block linear model changes the selected-minus-excluded ICL residual to {fmt(adjusted_selection['by_filter']['i']['adjusted_selected_minus_excluded'])} nJy arcsec⁻² in i (95% interval {ci2(adjusted_selection['by_filter']['i']['block_bootstrap_ci95'])}) and {fmt(adjusted_selection['by_filter']['i2']['adjusted_selected_minus_excluded'])} in i2 ({ci2(adjusted_selection['by_filter']['i2']['block_bootstrap_ci95'])}). Both intervals include zero. This adjustment is limited to 9 i and 21 i2 mixed blocks and does not correct the parent sample.",
        "",
        "The 17×17 stationary finite-lag covariance contributes about "
        f"{100*aperture['i']['17x17_pixel_noise_sigma_to_spatial_interval_halfwidth']:.1f}% (`i`) and "
        f"{100*aperture['i2']['17x17_pixel_noise_sigma_to_spatial_interval_halfwidth']:.1f}% (`i2`) of the current",
        "spatial-block interval half-widths for the 200–300 kpc aperture. It omits",
        "source Poisson, covariance beyond eight pixels, sky transfer, and PSF wings.",
        "A separate direct coadd comparison measured covariance in 200–300 kpc",
        "and 350–500 kpc annuli. Every target/control and radial-difference interval",
        "includes zero, but strict fully contained masked patches leave only",
        f"{annular['paired_block_bootstrap']['i']['zero_lag']['icl_200_300']['n_pairs']} i and",
        f"{annular['paired_block_bootstrap']['i2']['zero_lag']['icl_200_300']['n_pairs']} i2 complete inner-annulus pairs; this is too little support to claim equivalence or",
        "transfer the outer covariance kernel to the ICL annulus. The outer ring",
        "may contain real cluster-envelope light.",
        "A broader delivered-coadd diagnostic then analyzed all 783 same-filter",
        "target/control pairs. With 64×64 patches, full annular support remained",
        f"low ({full_support[64]['support_by_filter'][0]['n_complete_pairs']}/302 i and",
        f"{full_support[64]['support_by_filter'][1]['n_complete_pairs']}/481 i2 pairs); no familywise-adjusted covariance metric",
        "resolved a difference. A 32×32-patch sensitivity, requiring at least",
        f"{full_annular[32]['method']['minimum_patches_per_center_annulus']} patches per center, increased support to",
        f"{full_support[32]['support_by_filter'][0]['n_complete_pairs']}/302 i and",
        f"{full_support[32]['support_by_filter'][1]['n_complete_pairs']}/481 i2 pairs, but support still shifts toward lower-redshift systems",
        f"(complete-minus-incomplete redshift SMD {full_support[32]['support_covariate_smd']['i']['redshift_complete_minus_incomplete_smd']:.2f} in i,",
        f"{full_support[32]['support_covariate_smd']['i2']['redshift_complete_minus_incomplete_smd']:.2f} in i2).",
        "In 32×32 results, within-filter familywise max-T intervals identify a",
        f"positive i-band ±8-pixel lag-sum difference ({ci5(full_support[32]['familywise_max_t_by_filter']['i']['metrics']['sum8']['max_t_simultaneous_ci95'])},",
        f"p={full_support[32]['familywise_max_t_by_filter']['i']['metrics']['sum8']['familywise_two_sided_max_t_p']:.3f}) and a negative i2 zero-lag difference",
        f"({ci5(full_support[32]['familywise_max_t_by_filter']['i2']['metrics']['zero_lag']['max_t_simultaneous_ci95'])},",
        f"p={full_support[32]['familywise_max_t_by_filter']['i2']['metrics']['zero_lag']['familywise_two_sided_max_t_p']:.3f}). Coverage also differs in both filters.",
        "These selected-support findings are sensitive to patch scale and mask",
        "coverage; lag sums can include diffuse astrophysical structure and are",
        "not pure detector-noise estimates. They diagnose the delivered coadds,",
        "do not justify a pixel-noise correction, and remain exploratory.",
        f"The separately approved HSC PDR3 random-point export returned {random_variance['n_returned_random_points']:,} points at",
        f"{random_variance['n_centers_with_points']:,}/{random_variance['n_planned_centers']:,} target/control centers. The per-center point-count distribution",
        f"has q16/median/q84 {random_variance['random_points_per_center_q16_median_q84'][0]}/",
        f"{random_variance['random_points_per_center_q16_median_q84'][1]}/",
        f"{random_variance['random_points_per_center_q16_median_q84'][2]}; {random_variance['matched_pairs_with_both_endpoint_statistics']['i_pix_variance']} of 798 matched pairs have a random-point variance statistic at both endpoints.",
        "This local sample measures sparse random-point variance and related",
        "metadata only; it is not a pixel variance/coadd-weight map and does not",
        "capture resampling covariance or source Poisson noise. No ICL correction",
        "is applied from these rows.",
        f"Adding target/control random-point variance, sky-RMS, input-count, and PSF-trace terms one at a time to the redshift/richness model on the same 94 ICL–satellite pairs changes the i2 partial r from {fmt(random_variance_relation['results']['i2']['delta_log_pix_variance']['baseline_partial_r'])} to {fmt(min(v['adjusted_partial_r'] for v in random_variance_relation['results']['i2'].values()))}–{fmt(max(v['adjusted_partial_r'] for v in random_variance_relation['results']['i2'].values()))}. All 16 cross-filter max-T change intervals include zero (minimum adjusted p={fmt(min(v['maxT_adjusted_p'] for band in ('i','i2') for v in random_variance_relation['results'][band].values()))}); no local-condition adjustment resolves the selected-sample relation. This underpowered sensitivity does not show that noise, masks, PSF wings, or selection have no effect, and no correction is applied.",
        "The approved patch-level frame metadata did not supply per-pixel coadd",
        "weights; the 10-test familywise audit found no resolved association with",
        "the coarse frame/visit/calibration predictors.",
        "A newly approved one-cluster HSC-I warp pilot compared delivered coadd",
        "variance with the inverse-variance sum reconstructed from nine visit",
        "VARIANCE planes. At 200–300 kpc the primary pixelwise coadd/prediction",
        f"ratio has median {fmt(warp_annulus['primary_coadd_to_inverse_variance_ratio']['q16_median_q84'][1], 5)} "
        f"(q16–q84 {fmt(warp_annulus['primary_coadd_to_inverse_variance_ratio']['q16_median_q84'][0], 5)}–"
        f"{fmt(warp_annulus['primary_coadd_to_inverse_variance_ratio']['q16_median_q84'][2], 5)}); "
        f"{100*warp_annulus['primary_coverage_fraction']:.1f}% of the annulus has at least three valid visits. A shared-mask",
        f"sensitivity retains {100*warp_annulus['strict_coverage_fraction']:.1f}% and has median ratio "
        f"{fmt(warp_annulus['strict_coadd_to_inverse_variance_ratio']['q16_median_q84'][1], 5)}. Two of nine",
        "visits have no finite image/variance pixels, and one has partial coverage.",
        "This single-position consistency check is not an independent noise",
        "calibration and does not measure source Poisson, effective gain,",
        "resampling covariance, sky-transfer covariance, or total ICL uncertainty.",
        "The approved extension covers all 94 matched target/control pairs (188 requested centers). Four centers returned repeat HTTP 404, while all expected visits were present at the remaining 184 centers. The endpoint returned 161 additional visit products within those same approved cutouts; the download totals 9.26 GiB, below its 12 GiB cap.",
        f"Across 184 analyzed centers, the median of center-level 200–300 kpc pixel-median variance ratios is {fmt(warp94['radial_200_300_primary_ratio_median_across_centers'], 5)}. Strict shared-mask coverage has median {100*warp94['radial_200_300_strict_coverage_median']:.1f}% of the full annulus. The paired target-minus-control ratios are {fmt(warp94['paired_target_minus_control_primary_median_ratio']['i']['median_target_minus_control'], 5)} in i ({ci5(warp94['paired_target_minus_control_primary_median_ratio']['i']['block_bootstrap_ci95'])}) and {fmt(warp94['paired_target_minus_control_primary_median_ratio']['i2']['median_target_minus_control'], 5)} in i2 ({ci5(warp94['paired_target_minus_control_primary_median_ratio']['i2']['block_bootstrap_ci95'])}); both sky-block intervals include zero.",
        "This remains an approximate variance-plane consistency check: nearest-neighbor WCS sampling does not reproduce the exact HSC resampling kernel or covariance, and the analysis does not measure source Poisson, sky transfer, PSF wings, or total ICL-estimator covariance.",
        "A kernel-aware secondary mapping then used separable Lanczos-5 coefficients to move the already-resampled visit-warp VARIANCE planes onto coadd pixels, propagating diagonal input variance with squared normalized coefficients. HSC PDR3 describes fifth-order Lanczos resampling in its [processing documentation](https://hsc-release.mtk.nao.ac.jp/doc/index.php/processing__pdr3/). Across 184 centers, median target/control coadd-to-prediction ratios are 1.0071/1.0076 in i and 1.0064/1.0061 in i2; strict common-mask coverage medians range from 69.0% to 71.9%. Relative to nearest-neighbor, the median center-level ratio changes by only -0.00024 in i and -0.00016 in i2. Paired target-control block intervals include zero in both filters. This is a mapping sensitivity only: PDR3's Lanczos resampling of CCD data already induces off-diagonal warp covariance, which these diagonal warp variances do not contain; exact coadd weights, clipping, and mask interpolation also remain unreproduced.",
        "A second, image-level sensitivity recomputed the ICL 200–300 kpc target-minus-control aperture under the strict warp-valid mask while retaining the original K=2.5 source masks. Of 90 pairs with both warp cutouts, nine fall below the original 40% target/control coverage floor; 81 remain (23 i, 58 i2). Comparing baseline and strict masks on those same retained pairs, partial r changes by +0.002 in i and -0.0048 in i2. The connected-block paired 95% intervals for these changes include zero (i: "
        f"{ci5(mask_assoc['association_by_filter']['i']['paired_delta_ci95'])}; i2: "
        f"{ci5(mask_assoc['association_by_filter']['i2']['paired_delta_ci95'])}). The retained-sample i2 baseline r is {fmt(mask_assoc['association_by_filter']['i2']['baseline_partial_r'], 3)} and strict-mask r is {fmt(mask_assoc['association_by_filter']['i2']['strict_partial_r'], 3)}; the change is small, but the 13-pair loss from the planned 94 and residual selection mean this does not resolve the parent-sample association.",
        "A second sensitivity added five paired CCD-frame metadata offsets one at",
        "a time to the supported ICL–satellite correlation. Across ten",
        f"filter-by-adjustment tests, the smallest max-T adjusted p is {fmt(frame_satellite_pmin, 3)}, and",
        "all simultaneous intervals include zero. In i2, adjusted r remains",
        f"{fmt(min(frame_satellite_i2_r))}–{fmt(max(frame_satellite_i2_r))} around the 0.293 baseline. This does not test pixel-specific",
        "depth or noise and leaves the association exploratory.",
        "A separate post-hoc model-form check adds redshift/richness curvature",
        "and their interaction to the same linear baseline. The i2 estimates stay",
        f"between {fmt(min(nonlinear_i2_r))} and {fmt(max(nonlinear_i2_r))}; the smallest",
        f"eight-test max-T adjusted p is {fmt(nonlinear_pmin, 3)}. No adjustment",
        "resolves the selected-sample association or validates a measurement.",
        "",
        "Target-only injection transfer is incomplete: valid matrices cover",
        f"{cross['n_pairs_with_valid_injection_transfer']} of {cross['n_observed_association_pairs']} supported pairs. "
        "The matrix inversion is unstable, generates unphysical amplitudes, and",
        "its spatial bootstrap conditions on fixed responses. Across the 0.5–2.0",
        "satellite-injection amplitude grid, the inferred i2 association changes",
        "substantially; this is a sensitivity diagnostic, not de-mixed science.",
        "Do not apply a median leakage coefficient or promote the inversion.",
        "An independent audit of the hybrid-PSF satellite-shell-exclusion injections",
        "verified all 48 saved target/filter frames against the fixed association",
        "sample. Removing supported satellite-shell pixels changes the diffuse-to-ICL",
        f"median response by {shell['filters'][0]['responses']['icl_response_to_diffuse_unit_injection']['paired_shell_excluded_minus_fixed']['median']:+.5f} (`i`) and",
        f"{shell['filters'][1]['responses']['icl_response_to_diffuse_unit_injection']['paired_shell_excluded_minus_fixed']['median']:+.5f} (`i2`), while paired median mask coverage falls by",
        f"{100*abs(shell['filters'][0]['mask_coverage_200_300_kpc']['paired_shell_excluded_minus_fixed']['median']):.2f} and",
        f"{100*abs(shell['filters'][1]['mask_coverage_200_300_kpc']['paired_shell_excluded_minus_fixed']['median']):.2f} percentage points. The sample ranges include no common scalar correction; the",
        "public HSC-I extended wing is not matched to each frame, and is only a proxy",
        "for i2. The audit is model sensitivity, not observed-flux validation.",
        "A leave-one-connected-block-out prediction check adds ICL to a",
        "redshift/richness model for satellite outer light. In i, held-out MSE",
        f"increases by {-100*spatial_cv['i']['relative_oof_mse_reduction_from_icl']:.1f}%; in i2 it decreases by",
        f"{100*spatial_cv['i2']['relative_oof_mse_reduction_from_icl']:.1f}%. Both block-bootstrap intervals for paired",
        "error gain include zero, so this selected sample does not show clear",
        "out-of-block predictive improvement from ICL.",
        "A complementary radial analysis fixed the cluster set across all seven",
        f"satellite bins (N={fixed_profile['fixed_support']['i']['n_fixed_pairs']} in i, "
        f"N={fixed_profile['fixed_support']['i2']['n_fixed_pairs']} in i2). The outer i2",
        "slope is positive pointwise, but its 14-test max-T p-value is 0.410;",
        "the minimum familywise max-T p-value is 0.237. No radial bin is",
        "resolved after the declared familywise correction.",
        "A complementary median ICL split uses the same fixed complete-profile",
        "clusters. Adjusted outer high-minus-low differences are +4.61 nJy",
        "arcsec^-2 in i and +1.35 in i2; the minimum max-T p across 14 tests",
        f"is {fmt(fixed_split['minimum_maxT_p_14'])}. The split profile likewise shows no resolved radial contrast.",
        "",
        "### Forced unmatched-position radial sensitivity",
        "Across the 224 ICL-supported pairs, forced 20/20%-coverage profiles provide",
        "nine supported filter/radius tests. The response is light at unmatched",
        "CAMIRA positions, not confirmed satellite-galaxy light. Models adjust for",
        "redshift and log richness; inference uses connected 5-degree blocks and",
        "10,000 restricted Webb wild-cluster draws. No max-T adjusted p-value is",
        "below 0.05. The i 50–75 kpc raw p-value is 0.030, but its max-T p-value",
        "is 0.059 and Holm p is 0.274. These post-hoc results remain exploratory:",
        "the 20/20 coverage sample is outcome-selected, inner bins have sparse",
        "support, and measurement/noise/PSF/sky-transfer uncertainty is not",
        "propagated. They do not establish or exclude an ICL–satellite relation.",
        "",
        "| Filter | Radius (kpc) | Clusters | 5° blocks | ICL slope per SD | CR1 95% CI | Webb p | max-T p |",
        "|---|---:|---:|---:|---:|---:|---:|---:|",
        *[f"| {row['filter']} | {row['radial_bin_kpc']} | {row['n_clusters']} | {row['n_connected_5deg_blocks']} | {fmt(row['icl_slope_per_sd'])} | [{fmt(row['cr1_ci95_low'])}, {fmt(row['cr1_ci95_high'])}] | {fmt(row['p_webb_raw'])} | {fmt(row['p_maxT_across_supported_bins'])} |"
          for row in forced_radial["validated_results"]],
        "",
        "### S23B candidate coverage screen",
        f"A metadata-only PDR3 audit covers {s23b_coverage['n_candidates']:,} S23B candidate centers and {s23b_coverage['n_probes']:,} center/boundary probes. {s23b_coverage['n_candidates_all_9_probes_covered']:,} candidates have HSC-I mosaic matches at the center and all eight 1-Mpc boundary probes; all of those also carry the image-product and background-product availability flags.",
        f"The export returned {s23b_coverage['n_returned_patch_metadata_rows']:,} rows across 27 batches. This is a nine-point footprint and QA-metadata screen only: it does not validate every pixel, full-aperture masks/WCS, matched controls, or actual cutout products. New image retrieval remains a separate approval-gated step.",
        "",
        "### Radial contrast with the 300–500 kpc profile",
        "The full 50–500 kpc complete-case stack requires 40% common coverage in all six annuli with the 2.5 Kron mask. The support is 20 i-band pairs in 10 connected blocks and 48 i2 pairs in 30 blocks.",
        "The exploratory contrast subtracts the area-weighted 300–500 kpc target-minus-control surface brightness from the 200–300 kpc value. The 300–500 kpc region can contain cluster light and is not blank sky. Source completeness, source-Poisson variance, exact coadd covariance, sky transfer, and extended PSF wings remain unresolved.",
        "| Filter | Extraction | Pairs | Blocks | Contrast (nJy arcsec⁻²) | 95% block-bootstrap interval |",
        "|---|---|---:|---:|---:|---:|",
        *[f"| {row['filter']} | {row['extraction']} | {row['n_pairs']} | {row['n_connected_blocks']} | {fmt(float(row['delta_200_300_minus_outer_njy_arcsec2']))} | [{fmt(float(row['block_bootstrap_ci95_low']))}, {fmt(float(row['block_bootstrap_ci95_high']))}] |"
          for row in outer_contrast["validated_results"]],
        "The raw and outer-plane-corrected i intervals exclude zero. Both i2 intervals include zero. The profile contrast is exploratory and does not isolate pure ICL or establish an ICL–satellite relation.",
        "",
        "### Relation to the satellite-centered outer profile",
        "A cluster-level partial-correlation analysis joins that radial ICL contrast to the existing detected-member satellite-centered 50–75 proper-kpc target-minus-azimuthal-control profile, adjusting for redshift and log richness. The sample contains 20 i-band pairs in 10 connected blocks and 47 i2 pairs in 30 blocks; one otherwise eligible i2 cluster has no supported satellite profile and is kept missing.",
        "| Filter | ICL extraction | Clusters | Blocks | Partial r | 95% block-bootstrap interval |",
        "|---|---|---:|---:|---:|---:|",
        *[f"| {row['filter']} | {row['icl_extraction']} | {row['n_clusters']} | {row['n_connected_blocks']} | {fmt(float(row['partial_correlation_control_z_logrichness']))} | [{fmt(float(row['block_bootstrap_ci95_low']))}, {fmt(float(row['block_bootstrap_ci95_high']))}] |"
          for row in outer_relation["validated_results"]],
        "All four intervals span zero. This selected-sample exploratory association is not resolved by the present data and does not establish that satellite outer light and ICL share a physical component. The predictor is detected-member light after azimuthal-control subtraction, not a complete satellite census. Measurement error, source-Poisson noise, coadd covariance, sky transfer, and extended PSF wings are not propagated.",
        "",
        "## Background-galaxy weak lensing",
        "",
        f"The paired cluster-minus-random stack uses {wl['input']['matched_sample']['n_common_shape_support_lenses']} lenses and "
        f"{wl['input']['matched_sample']['n_random_centers_in_common_support']} matched random centers, with the existing S19A",
        "shape calibration and HSC-Y3 redshift calibration. The tangential profile",
        "is positive and declining, but the footprint audit assigns lenses to only",
        f"{regions['n_occupied_stack_components']} occupied connected tract-center components. Deleting the largest",
        f"component removes {regions['region_deletions'][0]['deleted_lens_count']} of {regions['common_support_lenses']} lenses and shifts the",
        "profile by 28.6% at 0.1–0.2 Mpc and 26.8% at 0.4–0.7 Mpc. Three occupied",
        "regions are insufficient for a calibrated field covariance; the lensing",
        "amplitude and its formal significance remain exploratory.",
        f"The random-only placebo reuses {random_placebo['n_random_centers']} centers and {random_placebo['n_random_pairs_used']:,} selected source pairs. Its seven-bin tangential",
        f"profile is consistent with zero (χ²={fmt(random_placebo['global_null_tests']['tangential']['chi_square_vs_zero'], 2)}, rank {random_placebo['global_null_tests']['tangential']['effective_dof_covariance_rank']}, p={fmt(random_placebo['global_null_tests']['tangential']['p_value_chi_square_approx'], 3)}); the cross profile is also consistent (p={fmt(random_placebo['global_null_tests']['cross']['p_value_chi_square_approx'], 3)}). This checks for a large random-stack residual, but does not validate the cluster signal or the approximate block covariance.",
        "",
        f"The continuous ICL–lensing relation is measured only for {wl_icl['n_common_support_lenses']} i2 lenses. "
        "The nominal outer-bin slope has p=0.021, but the seven-bin multiple-test",
        "minimum adjusted p is 0.149; no radial bin is significant after that",
        "correction. A proposed estimator-weight boost changes the candidate",
        "profile by 12.5–16.3%, but it is not validated as a member-dilution",
        "correction and has not been applied.",
        "",
        "A familywise audit of all 29 ICL–satellite–WL relations in their selected",
        "27-lens i2 overlap finds no resolved association: the cross-component",
        f"satellite relation at 0.4–0.7 Mpc has max-T p={fmt(cross_satellite['wild_cluster_p_max_t_29'], 3)}, and",
        f"the tangential 2–3 Mpc relation has p={fmt(tangential_satellite['wild_cluster_p_max_t_29'], 3)}. The minimum adjusted p across the",
        f"family is {fmt(min_triple_p, 3)}. This low-power common-support result does not establish no physical relation.",
        f"The i2 overlap audit shows that only 27 of 68 ICL–satellite-supported systems have random-subtracted WL support. Among 33 connected sky blocks, only 9 contain both groups; those blocks include {selection['mixed_block_selected_n']} selected and {selection['mixed_block_excluded_n']} excluded lenses. The other 42 lenses occupy single-status blocks, so spatial positivity fails for transporting the 27-lens relation to all 68. Four-variable within-block max-T balance found no resolved redshift, richness, ICL, or satellite contrast (minimum adjusted p={fmt(selection['familywise_maxT_p_min'], 3)}), but only nine mixed blocks inform that test. Do not IPW-correct the overlap to the parent sample.",
        "",
        "## Current scientific assessment and next gates",
        "",
        "1. Treat the ICL–satellite i2 correlation as an exploratory, selected-sample",
        "   signal. The member-count and measured host-position checks do not remove",
        "   it; the support-weighted interval does include zero. The current injections",
        "   do not support a stable physical de-mixing.",
        "2. Treat the weak-lensing stack as an exploratory profile, not a calibrated",
        "   halo-mass result. Its spatial influence and full field covariance are not",
        "   controlled by the current equal-area jackknife.",
        "3. Treat the Lanczos-5 result only as a secondary mapping sensitivity, not exact covariance validation; recover the upstream warp covariance or exact coadd resampling/weight behavior, then resolve source-Poisson response,",
        "   extended PSF wings, upstream sky transfer, satellite-support selection,",
        "   and validated independent HSC field geometry before a science claim.",
        "4. Increase independent field support, not just cluster count. Existing",
        "   catalog expansion alone has not established image coverage or matched",
        "   controls. Attempts to submit the approved S20A coverage-metadata",
        "   plan returned these HTTP statuses:",
        f"   {attempts_text}.",
        "   No CAS job ID or result rows were recorded. Diagnose the rejection",
        "   state before any further submission.",
        f"   The separately approved S23B metadata plan ({s23b_query_audit['authorized_scope']['candidates']} candidates, {s23b_query_audit['authorized_scope']['probes']} probes, 14 metadata columns) completed after two pre-job rejections (HTTP 401 at the wrong host and HTTP 406 on unsupported geometry SQL); the corrected HSC CAS query returned {s23b_query_audit['final_output_rows']:,} footprint/QA rows in {s23b_query_audit['completed_batches']} batches. No images, source rows, shapes, or photo-z were returned.",
        "",
        "## Provenance",
        "",
    ])
    for path, digest in sorted(provenance.items()):
        lines.append(f"- `{path}` — SHA-256 `{digest}`")
    lines.extend([
        f"- Report builder `{script_path}` — SHA-256 `{sha256(script_path)}`",
        "",
    ])
    return "\n".join(lines)


def build(data_root: Path, output_dir: Path) -> None:
    if output_dir.exists():
        raise FileExistsError(f"refusing to overwrite {output_dir}")
    payloads, provenance = load_sources(data_root)
    query_audit_path = QUERY_AUDIT if data_root == DATA_ROOT else data_root / QUERY_AUDIT.relative_to(DATA_ROOT)
    query_audit = json.loads(query_audit_path.read_text(encoding="utf-8"))
    if (query_audit.get("schema") != "astrafek-icl-expansion-coverage-submission-attempt-audit-v1"
            or query_audit.get("science_measurement") is not False
            or query_audit.get("scope", {}).get("n_candidates") != 2596
            or query_audit.get("scope", {}).get("n_probes") != 23364
            or query_audit.get("scope", {}).get("n_batches") != 52):
        raise ValueError("S20A coverage attempt audit has unexpected scope or status")
    expected_plan_hashes = {
        "9da220c05f1f24c409a602fcc4ddffc8cbd5bc33f459d335cc8587f1cd235006",
        "974d4e3db450703481ac2b004d37860806586900097e626c5b456702a9e0346f",
    }
    if any(item.get("plan_sha256") not in expected_plan_hashes
           or item.get("job_id_returned") is not False
           or item.get("http_status") not in (401, 406)
           for item in query_audit.get("attempts", [])):
        raise ValueError("S20A coverage attempt audit contains an unexpected plan/job state")
    if query_audit.get("result_csv_count") != 0:
        raise ValueError("S20A coverage attempt audit unexpectedly records result files")
    for path_text, expected_hash in query_audit.get("checkpoint_sha256", {}).items():
        checkpoint = Path(path_text)
        if not checkpoint.is_file() or sha256(checkpoint) != expected_hash:
            raise ValueError(f"S20A coverage checkpoint hash mismatch: {checkpoint}")
    provenance[str(query_audit_path)] = sha256(query_audit_path)
    s23b_query_audit_path = S23B_QUERY_AUDIT if data_root == DATA_ROOT else data_root / S23B_QUERY_AUDIT.relative_to(DATA_ROOT)
    s23b_query_audit = json.loads(s23b_query_audit_path.read_text(encoding="utf-8"))
    s23b_plan_path = s23b_query_audit_path.parents[1] / "coverage_query_plan.json"
    s23b_checkpoint_path = s23b_query_audit_path.parent / "coverage_status.json"
    s23b_summary_path = s23b_query_audit_path.parent / "coverage_summary.json"
    s23b_jobs = s23b_query_audit.get("attempts", [])
    s23b_scope = s23b_query_audit.get("authorized_scope", {})
    if (s23b_query_audit.get("schema") != "astrafek-s23b-coverage-metadata-submission-audit-v2"
            or s23b_query_audit.get("final_plan_sha256") != "62c28446f553bf45a940ea0126af1521d4ee0c1356c12cc73f42712df482c0c7"
            or s23b_query_audit.get("plan_sha256_verified") != sha256(s23b_plan_path)
            or s23b_scope != {"candidates": 1333, "probes": 11997, "max_rows": 47988, "columns": 14,
                              "images": False, "object_rows": False, "shapes": False, "photoz_rows": False}
            or s23b_query_audit.get("final_status") != "complete"
            or s23b_query_audit.get("final_output_rows") != 13575
            or s23b_query_audit.get("completed_batches") != 27
            or len(s23b_jobs) != 3
            or [item.get("http_status") for item in s23b_jobs] != [401, 406, 200]
            or any(item.get("cas_job_created") is not False or item.get("result_rows") != 0
                   for item in s23b_jobs[:2])
            or s23b_jobs[2].get("cas_job_created") is not True
            or s23b_jobs[2].get("completed_jobs") != 27
            or s23b_jobs[2].get("result_rows") != 13575
            or s23b_query_audit.get("credentials_saved") is not False
            or s23b_query_audit.get("science_measurement") is not False
            or s23b_query_audit.get("checkpoint_sha256") != sha256(s23b_checkpoint_path)
            or s23b_query_audit.get("summary_sha256") != sha256(s23b_summary_path)):
        raise ValueError("S23B coverage submission audit has unexpected scope or status")
    provenance[str(s23b_query_audit_path)] = sha256(s23b_query_audit_path)
    provenance[str(s23b_plan_path)] = sha256(s23b_plan_path)
    provenance[str(s23b_checkpoint_path)] = sha256(s23b_checkpoint_path)
    provenance[str(s23b_summary_path)] = sha256(s23b_summary_path)
    script_path = Path(__file__).resolve()
    markdown = render(payloads, provenance, script_path, query_audit, s23b_query_audit)
    output_dir.mkdir(parents=True, exist_ok=False)
    report_path = output_dir / "icl_wl_current_evidence_report.md"
    report_path.write_text(markdown, encoding="utf-8")
    summary_path = output_dir / "summary.json"
    summary = {
        "schema": "astrafek-icl-wl-current-evidence-report-v1",
        "status": "complete_exploratory_synthesis",
        "science_measurement": False,
        "hsc_network_used": False,
        "image_pixels_read": False,
        "input_sha256": provenance,
        "coverage_query_attempts": query_audit["attempts"],
        "s23b_coverage_query_attempt": s23b_query_audit,
        "report_sha256": sha256(report_path),
        "report_builder_sha256": sha256(script_path),
    }
    summary_path.write_text(json.dumps(summary, ensure_ascii=False, indent=2, sort_keys=True) + "\n",
                            encoding="utf-8")
    print(json.dumps({"status": summary["status"], "science_measurement": False,
                      "report": str(report_path), "summary": str(summary_path),
                      "inputs": len(provenance)}, ensure_ascii=False, indent=2))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, default=DATA_ROOT)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    build(args.data_root.resolve(), args.output_dir.resolve())


if __name__ == "__main__":
    main()
