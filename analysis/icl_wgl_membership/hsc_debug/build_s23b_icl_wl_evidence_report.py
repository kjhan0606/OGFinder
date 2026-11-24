#!/usr/bin/env python3
"""Build a hash-pinned S23B ICL, satellite, PSF, and weak-lensing checkpoint.

This report builder reads saved JSON summaries and locally prepared plans only.
It makes no network requests, downloads no data, and does not open image pixels.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import shlex
import sys


DATA_ROOT = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
SOURCES = {
    "icl_stack": "s23b_icl_target_control_stack_2026-11-07_01/summary.json",
    "icl_contrasts": "s23b_icl_relative_radial_contrasts_2026-11-07_02/summary.json",
    "member_profiles": "s23b_camira_member_outer_profiles_full_2026-11-07_01/summary.json",
    "member_icl_relation": "s23b_camira_member_outer_icl_relation_2026-11-07_01/summary.json",
    "member_injection": "s23b_camira_member_outer_injection_recovery_2026-11-07_01/summary.json",
    "direct_member_relation": "s23b_icl_satellite_relation_2026-11-07_03/summary.json",
    "field8_sensitivity": "s23b_icl_satellite_field8_weight_sensitivity_2026-11-07_01/summary.json",
    "field8_weighted": "s23b_icl_satellite_field8_weighted_density_2026-11-07_02/summary.json",
    "icl_covariance": "s23b_icl_relative_covariance_propagation_2026-11-07_02/summary.json",
    "gaia_expansion": "s23b_psf_expansion_2026-11-07_01/gaia_audit_2026-11-07_01/summary.json",
    "psf_pair_comparison": "s23b_psf_target_control_profile_comparison_2026-11-07_03/summary.json",
    "wl_overlap": "wl_s23b_camira_overlap_audit_2026-11-07_01/summary.json",
    "wl_source_union": "wl_s23b_existing_source_union_2026-11-07_01/summary.json",
    "wl_count_plan": "wl_s23b_shape_count_plan_2026-11-07_01/shape_count_plan.json",
}
APPROVED_GAIA_PLAN_SHA256 = "c00eda5f896f8160022a7e4198b3f8ce69803f0e30adc8f560ce04c8beba0e7d"
GAIA_RESPONSE_SHA256 = "151f7c3f227ab7a16b1a6015563afd46fd53ccf4067d44d0264e27be75f665fb"
WL_COUNT_PLAN_SHA256 = "b713eebded489c280fd9ce1b5840bc8d7d3e42d029035f8dbd414b72fb31e984"

MEMBER_DIAGNOSTICS = "s23b_member_bin_diagnostics_full_2026-11-17_03"
MEMBER_PSF = "s23b_member_psf_injection_run_2026-11-17_01"
CURRENT_ICL_SOURCES = {
    "stack": (SOURCES["icl_stack"], "c8bac8dbdd0a886604c216b2410725eb97883b0cdc1a91fd47e7a7cfcaaaa20b"),
    "contrasts": (SOURCES["icl_contrasts"], "a47c6679ad318abcbee5ee5ab46e1a984976bd01c6f1dcbbe3d6de99b89306a8"),
    "covariance": (SOURCES["icl_covariance"], "9c1ef1740490f7cbc4811a9bd2573d69405f5482d7ce754a09642108507fc94e"),
    "candidate_density": (SOURCES["direct_member_relation"], "fab29d7aef464707c31e2616fe1db9fc3bad9d281ad6d5dfa6d63f2b484bbb39"),
    "measurement": (f"{MEMBER_DIAGNOSTICS}/measurement/summary.json", "23b376b26ea7355c94cf1d556f437f32525998561bd3a94b057afe45b44c3b4c"),
    "audit": (f"{MEMBER_DIAGNOSTICS}/audit/summary.json", "50a158970521616dc0c9ce29826b1c2e98f15f260e69d30a578fc8479805db28"),
    "coverage": (f"{MEMBER_DIAGNOSTICS}/coverage_followup/analysis/summary.json", "f530bf87a59e45205cbc85ff7fec9c7929058e3090ad59fdaa1e9fe2e3404ccf"),
    "unconvolved_injection": (SOURCES["member_injection"], "e9d4cced59f0b0f8ea3f189989886f28af74a339cc5a9b3421f8bca187defc77"),
    "core_injection": (f"{MEMBER_PSF}/injection/summary.json", "7cdd34b10f3e976a1ecd94518c338e286aaa59fd4b42c310c3607da4c27431c4"),
    "paired_comparison": (f"{MEMBER_PSF}/comparison/summary.json", "3ecb1ad6900ac6255b4143b461782a3978e020e3dbed9e22b9209e534bec19cf"),
    "normalization": ("s23b_member_psf_common_normalization_2026-11-17_01/summary.json", "cfb3d8fc5422fad350939ec762b5e58533efc731209eebedfa2af784069e27a8"),
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_sources(data_root: Path) -> tuple[dict[str, dict], dict[str, str]]:
    payloads: dict[str, dict] = {}
    provenance: dict[str, str] = {}
    for name, relative in SOURCES.items():
        path = data_root / relative
        if not path.is_file():
            raise FileNotFoundError(f"required evidence summary is missing: {path}")
        payload = json.loads(path.read_text(encoding="utf-8"))
        if name != "wl_count_plan" and payload.get("science_measurement") is not False:
            raise ValueError(f"{name} unexpectedly has science_measurement=true")
        payloads[name] = payload
        provenance[str(path)] = sha256(path)

    stack = payloads["icl_stack"]
    if (stack.get("status") != "complete_exploratory"
            or stack.get("n_targets_available") != 1043
            or stack.get("n_pairs_stacked") != 1022
            or stack.get("n_by_filter") != {"i": 241, "i2": 781}
            or stack.get("images_ready_for_science") is not False):
        raise ValueError("S23B ICL stack scope or science gate changed")

    gaia = payloads["gaia_expansion"]
    gaia_plan = data_root / "s23b_psf_expansion_2026-11-07_01/gaia_plan/plan.json"
    if (not gaia_plan.is_file()
            or sha256(gaia_plan) != APPROVED_GAIA_PLAN_SHA256
            or gaia.get("approved_plan_sha256") != APPROVED_GAIA_PLAN_SHA256
            or gaia.get("aggregate_response_sha256") != GAIA_RESPONSE_SHA256
            or gaia.get("n_candidate_positions") != 4299
            or gaia.get("n_completed_batches") != 18
            or gaia.get("n_gaia_match_rows") != 706
            or gaia.get("psf_approved") is not False):
        raise ValueError("completed Gaia query does not match its approved plan and audit")
    provenance[str(gaia_plan)] = sha256(gaia_plan)

    wl_plan = payloads["wl_count_plan"]
    if (wl_plan.get("query_plan_sha256") != WL_COUNT_PLAN_SHA256
            or wl_plan.get("status") != "prepared_not_submitted"
            or wl_plan.get("query_submitted") is not False
            or wl_plan.get("network_used") is not False
            or wl_plan.get("execution_requires_explicit_user_approval") is not True
            or wl_plan.get("n_lenses") != 427
            or wl_plan.get("output_columns") != ["lens_id", "n_shape_sources"]):
        raise ValueError("pending HSC count-only plan changed or is no longer unsubmitted")

    for name in ("icl_contrasts", "direct_member_relation", "member_icl_relation",
                 "member_injection", "field8_sensitivity", "field8_weighted",
                 "icl_covariance", "psf_pair_comparison", "wl_overlap", "wl_source_union"):
        if payloads[name].get("status", "").startswith("failed"):
            raise ValueError(f"{name} is not a completed evidence product")
    return payloads, provenance


def summarize(payloads: dict[str, dict], provenance: dict[str, str]) -> dict:
    stack = payloads["icl_stack"]
    contrasts = payloads["icl_contrasts"]
    contrast_rows = {
        row["filter"]: row for row in contrasts["primary_5deg_profiles"]
        if row.get("contrast") == "50-300_minus_300-800"
    }
    direct = payloads["direct_member_relation"]
    direct_primary = direct["primary_layout_results"]
    direct_sensitivity = direct["spatial_partition_sensitivity"]
    floor_rows = direct["common_support_floor_sensitivity_primary_layout"]
    field8 = payloads["field8_weighted"]
    field8_sensitivity = payloads["field8_sensitivity"]
    member_relation = payloads["member_icl_relation"]
    profiles = payloads["member_profiles"]
    injection = payloads["member_injection"]
    gaia = payloads["gaia_expansion"]
    psf = payloads["psf_pair_comparison"]
    wl = payloads["wl_overlap"]
    union = payloads["wl_source_union"]
    wl_plan = payloads["wl_count_plan"]

    profile_support = {
        row["filter"]: row["paired_frame_count_by_profile_bin"]
        for row in psf["filters_and_subsets"]
        if row.get("subset") == "all"
    }
    direct_rows = {}
    for band in ("i", "i2"):
        row = direct_primary[band]
        floors = floor_rows[band]
        ten_percent = next(
            (item for item in floors if item.get("minimum_common_support_fraction_both_annuli") == 0.1),
            None,
        )
        weighted_primary = field8["primary_layout_results"][band]
        direct_rows[band] = {
            "n_pairs": row["n_pairs"],
            "n_blocks": row["n_connected_blocks"],
            "unweighted_slope": row["density_ratio_slope_njy_arcsec2_per_dex"],
            "unweighted_ci95": [row["block_bootstrap_ci95_low"], row["block_bootstrap_ci95_high"]],
            "layouts_excluding_zero": direct_sensitivity[band]["n_layouts_ci_excluding_zero"],
            "n_layouts": direct_sensitivity[band]["n_layouts"],
            "ten_percent_support": ({
                "n_pairs": ten_percent["n_pairs"],
                "slope": ten_percent["density_ratio_slope_njy_arcsec2_per_dex"],
                "ci95": [ten_percent["block_bootstrap_ci95_low"], ten_percent["block_bootstrap_ci95_high"]],
            } if ten_percent else None),
            "field8_weighted_slope": weighted_primary["weighted_slope_field8_njy_arcsec2_per_dex"],
            "field8_weighted_ci95": [weighted_primary["weighted_ci95_low"], weighted_primary["weighted_ci95_high"]],
            "field8_layouts_excluding_zero": field8_sensitivity["spatial_partition_sensitivity"][band]["n_layouts_ci_excluding_zero"],
            "field8_n_layouts": field8_sensitivity["spatial_partition_sensitivity"][band]["n_layouts"],
        }

    member_outer_relation = {
        band: {
            "n_clusters": values["n_clusters"],
            "slope": values["slope"],
            "ci95": values["ci95"],
            "n_blocks": values["n_blocks"],
        }
        for band, values in member_relation["primary_results"].items()
    }
    injection_response = {}
    for band in ("i", "i2"):
        unit = injection["amplitude_response_by_filter"][band]["1.0"]
        injection_response[band] = {
            "n_supported_frames": unit["n_supported_frames"],
            "median_recovered_input_ratio": unit["median_recovered_response"],
            "fraction_detection_count_changed": unit["fraction_detection_count_changed"],
            "paired_amplitude_difference_ci_include_zero": all(
                interval[0] <= 0 <= interval[1]
                for key, interval in (
                    ("0.5_minus_1", injection["amplitude_response_by_filter"][band]["paired_0.5_minus_1"]["cluster_bootstrap_median_ci95"]),
                    ("2_minus_1", injection["amplitude_response_by_filter"][band]["paired_2_minus_1"]["cluster_bootstrap_median_ci95"]),
                )
            ),
        }

    covariance = payloads["icl_covariance"]
    covariance_at_primary = {
        row["filter"]: {
            "conditional_finite_lag_sigma": row["measured_finite_lag_noise_sigma_of_equal_pair_mean"],
            "pair_sem": row["pair_sem_njy_arcsec2"],
            "finite_lag_to_diagonal_sigma_ratio": row["finite_lag_to_diagonal_sigma_ratio"],
        }
        for row in covariance["filter_summary"]
        if row.get("contrast") == "50-300_minus_300-800"
    }

    return {
        "schema": "astrafek-s23b-icl-wl-evidence-report-v1",
        "status": "complete_exploratory_checkpoint",
        "scope": "S23B CAMIRA ICL, candidate satellite positions, PSF diagnostics, and their connection to the existing S19A weak-lensing sample",
        "this_report_network_used": False,
        "this_report_image_pixels_read": False,
        "this_report_data_downloaded": False,
        "science_measurement": False,
        "images_ready_for_science": False,
        "psf_approved": False,
        "provenance": provenance,
        "icl_stack": {
            "targets_available": stack["n_targets_available"],
            "matched_pairs_stacked": stack["n_pairs_stacked"],
            "pairs_by_filter": stack["n_by_filter"],
            "filter_mismatch_excluded": stack["n_pairs_excluded_filter_mismatch"],
            "zero_support_excluded": stack["n_pairs_excluded_zero_support"],
            "primary_relative_contrasts": contrast_rows,
            "conditional_pixel_covariance": covariance_at_primary,
        },
        "satellite_relation": {
            "unweighted_candidate_position_relation": direct_rows,
            "member_centered_outer_light_relation": member_outer_relation,
            "member_profile_sample": profiles["sample"],
            "injection_recovery_at_unit_amplitude": injection_response,
            "field8_semantics": field8["field_interpretation"],
        },
        "psf_diagnostics": {
            "gaia": {
                "approved_plan_sha256": gaia["approved_plan_sha256"],
                "response_sha256": gaia["aggregate_response_sha256"],
                "candidate_positions": gaia["n_candidate_positions"],
                "batches_completed": gaia["n_completed_batches"],
                "matched_sources": gaia["n_gaia_match_rows"],
                "strict_diagnostic_candidates": gaia["evidence_levels"]["strict_candidate"],
                "non_matches": gaia["evidence_levels"]["no_gaia_match"],
                "max_match_separation_arcsec": gaia["max_separation_arcsec"],
            },
            "paired_profile_support": profile_support,
            "paired_flux_difference_estimable": psf["paired_flux_difference_estimable"],
            "psf_approved": psf["psf_approved"],
        },
        "weak_lensing_overlap": {
            "existing_sample_counts": wl["sample_counts"],
            "nearest_s23b_to_s19a_lens_arcsec": wl["positional_match"]["nearest_separation_arcsec_quantiles"]["min"],
            "matches_within_60_arcsec": wl["positional_match"]["n_matches_within_60_arcsec_in_common_support_435"],
            "old_export_apertures_with_any_source_rows": union["targets_with_any_union_sources"],
            "old_export_apertures_without_source_rows": union["targets_with_no_union_sources"],
            "old_export_union_area_fraction": union["all_query_union_area_fraction"],
            "pending_count_only_plan": {
                "query_plan_sha256": wl_plan["query_plan_sha256"],
                "status": wl_plan["status"],
                "n_lenses": wl_plan["n_lenses"],
                "radius_mpc_proper": wl_plan["radius_mpc_proper"],
                "n_batches": wl_plan["n_query_batches"],
                "output_columns": wl_plan["output_columns"],
                "estimated_result_rows": wl_plan["estimated_result_rows"],
                "query_submitted": wl_plan["query_submitted"],
                "explicit_approval_required": wl_plan["execution_requires_explicit_user_approval"],
            },
        },
        "conclusion": {
            "icl": "A positive 50–300 kpc versus 300–800 kpc relative excess appears in both filters, but this is a candidate profile gradient, not a validated ICL measurement.",
            "satellite_association": "The unweighted i2 candidate-density slope is positive in the primary layout but weakens to an interval including zero at 10% common-area support and is not significant in all block layouts. The field-8-weighted result is a sensitivity only because the headerless field semantics are unverified.",
            "member_outer_light": "The adjusted outer member-light relation has intervals including zero in both filters; the current data do not establish an ICL–satellite outer-structure association.",
            "weak_lensing": "The existing S19A shape sample has no positional overlap with these S23B targets. Existing exported rows cover only a small, geometrically incomplete subset, so there is no same-cluster three-way WL–ICL–satellite result.",
            "next_gate": "The prepared 427-target HSC count-only plan remains unsubmitted and needs its own explicit approval. Any shape-row export and matched-random support query would be separate approvals.",
        },
        "limitations": [
            "The target-minus-control stack has sparse central common-pixel support; its outer annulus is a reference, not proven sky.",
            "The conditional covariance propagation assumes stationary outer-field kernels and omits full residual-sky, source-Poisson, PSF, extinction, mask, centering, and spatial covariance uncertainty.",
            "CAMIRA member positions are photometric candidates; the profile subset is coverage-selected and member-file field 8 is not authoritatively mapped.",
            "The injection test is conditional on previously detected members and omits PSF convolution and source Poisson noise.",
            "Gaia matches are astrometric diagnostics, not a calibrated stellar or PSF classifier; effective HSC coadd epochs are unknown.",
            "The older WL export did not apply the full source-redshift calibration, B-mode polarity, or shape-quality selection needed for a new S23B lensing measurement.",
        ],
    }


def render_markdown(report: dict) -> str:
    icl = report["icl_stack"]
    satellite = report["satellite_relation"]
    psf = report["psf_diagnostics"]
    wl = report["weak_lensing_overlap"]
    lines = [
        "# S23B ICL, satellite, PSF, and weak-lensing evidence checkpoint",
        "",
        "## Research status",
        "",
        "This hash-pinned report consolidates the saved S23B exploratory products. It reads summaries and plans only; it made no network request, downloaded no data, and read no image pixels. `science_measurement=false`, `images_ready_for_science=false`, and `psf_approved=false` remain in force.",
        "",
        "## ICL stack",
        "",
        f"The stack contains {icl['matched_pairs_stacked']:,} matched target/control pairs from {icl['targets_available']:,} S23B targets: {icl['pairs_by_filter']['i']:,} in HSC-i and {icl['pairs_by_filter']['i2']:,} in HSC-i2. Seventeen pairs were removed for filter mismatch and two for zero support.",
        "",
        "The primary 50–300 kpc minus 300–800 kpc relative contrasts are:",
        "",
        "| Filter | Pairs | Contrast (nJy arcsec⁻²) | 95% 5° spatial-block interval |",
        "|---|---:|---:|---:|",
    ]
    for band in ("i", "i2"):
        row = icl["primary_relative_contrasts"][band]
        lines.append(f"| {band} | {row['n_supported_pairs']} | {row['mean_contrast_njy_arcsec2']:.2f} | [{row['block_bootstrap_ci95_low']:.2f}, {row['block_bootstrap_ci95_high']:.2f}] |")
    lines += [
        "",
        "Both contrasts are positive, but common-mask support is sparse and the 300–800 kpc annulus is not established as blank sky. Conditional finite-lag pixel-noise σ is only a subset of the uncertainty budget; it does not correct the spatial-block intervals.",
        "",
        "## Satellite and member outer-light relations",
        "",
        "The unweighted relation uses CAMIRA candidate-position density contrast against the ICL relative contrast, adjusted for richness, redshift, extinction mismatch, and common support:",
        "",
        "| Filter | Pairs / blocks | Slope | 95% block interval | Layouts excluding zero | 10% support interval |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for band in ("i", "i2"):
        row = satellite["unweighted_candidate_position_relation"][band]
        ten = row["ten_percent_support"]
        lines.append(f"| {band} | {row['n_pairs']} / {row['n_blocks']} | {row['unweighted_slope']:.2f} | [{row['unweighted_ci95'][0]:.2f}, {row['unweighted_ci95'][1]:.2f}] | {row['layouts_excluding_zero']}/{row['n_layouts']} | [{ten['ci95'][0]:.2f}, {ten['ci95'][1]:.2f}] (N={ten['n_pairs']}) |")
    lines += [
        "",
        "The HSC-i2 coefficient is positive in the primary layout, but its interval includes zero after requiring at least 10% common support in both annuli, and only four of six spatial layouts exclude zero. A raw-field-8 weighted sensitivity has a positive HSC-i2 coefficient, but the headerless field's meaning is not verified, so it cannot establish robustness.",
        "",
        "A separate member-centered outer-light fit gives adjusted slopes and intervals:",
        "",
        "| Filter | Clusters | Slope | 95% block interval |",
        "|---|---:|---:|---:|",
    ]
    for band in ("i", "i2"):
        row = satellite["member_centered_outer_light_relation"][band]
        lines.append(f"| {band} | {row['n_clusters']} | {row['slope']:.3f} | [{row['ci95'][0]:.3f}, {row['ci95'][1]:.3f}] |")
    lines += [
        "",
        "Both intervals include zero. The injection recovery is conditional on already-detected members; it does not measure catalog completeness and omits PSF convolution and source Poisson noise.",
        "",
        "## PSF evidence",
        "",
        f"The already approved Gaia DR3 expansion is complete: {psf['gaia']['candidate_positions']:,} morphology candidates in {psf['gaia']['batches_completed']} batches yielded {psf['gaia']['matched_sources']} one-to-one matches; {psf['gaia']['strict_diagnostic_candidates']} pass a conservative astrometric screen. These are diagnostic candidates only. The paired target/control profile comparison has no radial bin with two valid matched frame pairs in either filter, so it cannot calibrate the PSF-wing difference. PSF approval remains false.",
        "",
        "## Weak-lensing connection",
        "",
        f"The existing S19A lensing sample has {wl['existing_sample_counts']['wl_lenses_in_verified_common_shape_random_support']} lenses with verified common source/random support, but none lies within 60 arcsec of an S23B target; the nearest center is {wl['nearest_s23b_to_s19a_lens_arcsec']:.1f} arcsec away. Existing lens and random source exports place rows in only {wl['old_export_apertures_with_any_source_rows']}/427 S23B apertures; {wl['old_export_apertures_without_source_rows']} have none, and estimated union-area coverage has median {wl['old_export_union_area_fraction']['median']:.3f}.",
        "",
        f"The next prepared gate is an HSC PDR3 count-only query for {wl['pending_count_only_plan']['n_lenses']} S23B clusters within {wl['pending_count_only_plan']['radius_mpc_proper']:.0f} proper Mpc. It returns only `lens_id` and an integer source count ({wl['pending_count_only_plan']['estimated_result_rows']} rows, {wl['pending_count_only_plan']['n_batches']} batches), has not been submitted, and still requires explicit approval. Shape rows and matched-random counts are separate future approvals.",
        "",
        "## Reproducibility",
        "",
        "The machine-readable `summary.json` records SHA-256 for every source summary and the approved Gaia plan. Rebuild with:",
        "",
        "```bash",
        "python3 analysis/icl_wgl_membership/hsc_debug/build_s23b_icl_wl_evidence_report.py \\",
        "  --output-dir /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_s23b_wl_evidence_2026-11-07_01",
        "```",
        "",
        "## Remaining scientific gates",
        "",
    ]
    lines.extend(f"- {item}" for item in report["limitations"])
    return "\n".join(lines) + "\n"


def current_icl_report(data_root: Path) -> dict:
    """Replay the current ICL evidence from exact saved local products."""
    products = {}
    provenance: dict[str, str] = {}

    def verify(path: Path, expected: str) -> None:
        actual = sha256(path)
        if actual != expected:
            raise ValueError(f"ICL evidence SHA-256 differs: {path}")
        provenance[str(path)] = actual

    for name, (relative, digest) in CURRENT_ICL_SOURCES.items():
        path = data_root / relative
        verify(path, digest)
        payload = json.loads(path.read_text(encoding="utf-8"))
        if payload.get("science_measurement") is not False:
            raise ValueError(f"Missing diagnostic science flag in {path}")
        products[name] = payload

    stack = products["stack"]
    contrasts = products["contrasts"]
    covariance = products["covariance"]
    measurement = products["measurement"]
    audit = products["audit"]
    coverage = products["coverage"]
    injection = products["core_injection"]
    norm = products["normalization"]
    source_hashes = {name: digest for name, (_, digest) in CURRENT_ICL_SOURCES.items()}
    contrast_inputs = [digest for name, digest in contrasts["inputs_sha256"].items()
                       if Path(name).parent.name == Path(SOURCES["icl_stack"]).parent.name
                       and Path(name).name == "summary.json"]
    if contrast_inputs != [source_hashes["stack"]]:
        raise ValueError("Relative-contrast stack lineage differs")
    if (coverage["measurement_summary_sha256"] != source_hashes["measurement"]
            or audit["measurement_summary_sha256"] != source_hashes["measurement"]
            or coverage["audit_summary_sha256"] != source_hashes["audit"]
            or not audit["historical_profile_rows_identical"]
            or not audit["historical_candidate_diagnostics_identical"]
            or audit["max_flux_reconstruction_difference"] != 0):
        raise ValueError("Member coverage results do not match the audited profile replay")
    if (norm["input_sha256"]["new_summary"] != source_hashes["core_injection"]
            or norm["input_sha256"]["old_summary"] != source_hashes["unconvolved_injection"]
            or norm["input_sha256"]["paired_summary"] != source_hashes["paired_comparison"]
            or injection["psf_convolution_applied"] is not True
            or injection["n_frames_failed"] != 0):
        raise ValueError("PSF normalization lineage or receipt differs")

    interval_checks = {}
    coverage_dir = data_root / Path(CURRENT_ICL_SOURCES["coverage"][0]).parent
    for filename in ("coverage_slope_layouts.tsv", "paired_slope_changes.tsv"):
        path = coverage_dir / filename
        verify(path, coverage["output_sha256"][filename])
        counts = {"supported": 0, "insufficient_sample": 0, "intervals_excluding_zero": 0}
        with path.open(encoding="utf-8", newline="") as stream:
            for row in csv.DictReader(stream, delimiter="\t"):
                if row["status"] == "insufficient_sample":
                    if row["ci_low"] or row["ci_high"] or row["slope"]:
                        raise ValueError(f"Unsupported row contains an estimate in {path}")
                    counts["insufficient_sample"] += 1
                elif row["status"] == "complete":
                    low, high = float(row["ci_low"]), float(row["ci_high"])
                    if not all(map(math.isfinite, (low, high))) or low > high:
                        raise ValueError(f"Invalid interval in {path}")
                    counts["supported"] += 1
                    counts["intervals_excluding_zero"] += int(not low <= 0 <= high)
                else:
                    raise ValueError(f"Unexpected sensitivity status in {path}")
        interval_checks[filename] = counts

    for label, filename, expected in (
        ("contrasts", "paired_relative_contrasts.tsv", contrasts["outputs"]["paired_relative_contrasts.tsv"]["sha256"]),
        ("normalization", "common_normalization_paired.tsv", norm["output_sha256"]["common_normalization_paired.tsv"]),
    ):
        verify(data_root / Path(CURRENT_ICL_SOURCES[label][0]).parent / filename, expected)

    return {
        "schema": "astrafek-s23b-icl-current-evidence-v1",
        "status": "complete_exploratory_evidence_replay",
        "science_measurement": False,
        "images_ready_for_science": False,
        "psf_approved": False,
        "this_report_network_used": False,
        "this_report_image_pixels_read": False,
        "this_report_data_downloaded": False,
        "runner": {"path": str(Path(__file__).resolve()), "sha256": sha256(Path(__file__))},
        "provenance": provenance,
        "icl_stack": {
            "targets": stack["n_targets_available"],
            "controls": stack["n_matched_controls_available"],
            "pairs_by_filter": stack["n_by_filter"],
            "relative_contrast": contrasts["estimand"],
            "primary_intervals": contrasts["primary_5deg_profiles"],
            "conditional_covariance": covariance["filter_summary"],
        },
        "satellite_candidate_density": products["candidate_density"]["primary_layout_results"],
        "member_coverage": {
            "frames": coverage["n_frames"],
            "measurements": measurement["sample"],
            "coverage_floors": coverage["coverage_floors"],
            "native_primary_layout": [row for row in coverage["primary_layout"] if row["mode"] == "native"],
            "interval_checks": interval_checks,
            "historical_flux_reconstruction_max_difference": audit["max_flux_reconstruction_difference"],
        },
        "psf_core_normalization": {
            "frames_completed": injection["n_frames_complete"],
            "frames_failed": injection["n_frames_failed"],
            "member_amplitude_rows": injection["n_member_amplitude_rows"],
            "paired_supported_rows": norm["n_paired_supported"],
            "unit_amplitude_results": [row for row in norm["results"] if row["amplitude"] == 1.0],
        },
        "open_science_gates": [
            "Extended PSF wings and residual sky transfer",
            "Source Poisson noise and effective coadd weights",
            "Long-range and shared-image covariance",
            "Mask-selection effects and photometric membership completeness/purity",
        ],
    }


def render_current_icl(report: dict, data_root: Path, output: Path) -> str:
    stack = report["icl_stack"]
    coverage = report["member_coverage"]
    psf = report["psf_core_normalization"]
    parts = [
        "# S23B ICL and satellite evidence replay", "",
        f"The local replay verified {len(report['provenance'])} saved summaries and tables. "
        "The work supports exploratory ICL and satellite diagnostics. `science_measurement=false`, "
        "`images_ready_for_science=false`, and `psf_approved=false`.", "",
        f"The stack uses {sum(stack['pairs_by_filter'].values()):,} matched target/control pairs "
        f"from {stack['targets']:,} targets and {stack['controls']:,} controls. The relative metric "
        "subtracts the 300–800 proper-kpc annulus, which may retain cluster light.", "",
        "| Filter | 50–300 minus 300–800 kpc (nJy arcsec⁻²) | 95% spatial-block interval |",
        "|---|---:|---:|",
    ]
    for row in stack["primary_intervals"]:
        if row["contrast"] == "50-300_minus_300-800":
            parts.append(f"| {row['filter']} | {row['mean_contrast_njy_arcsec2']:.2f} | "
                         f"[{row['block_bootstrap_ci95_low']:.2f}, {row['block_bootstrap_ci95_high']:.2f}] |")
    slope = coverage["interval_checks"]["coverage_slope_layouts.tsv"]
    change = coverage["interval_checks"]["paired_slope_changes.tsv"]
    parts += ["", f"The audited {coverage['frames']:,}-frame member profile replay reproduced the saved "
              "fluxes exactly. For member outer-light association, "
              f"{slope['intervals_excluding_zero']}/{slope['supported']} supported slope intervals "
              f"and {change['intervals_excluding_zero']}/{change['supported']} supported paired-change intervals "
              "exclude zero. The remaining configurations had insufficient support to estimate intervals.", "",
              "| Filter | Clusters at 40% coverage | 60% | 80% |",
              "|---|---:|---:|---:|"]
    for band in ("i", "i2"):
        values = {row["coverage_floor"]: row["n_target_clusters"] for row in coverage["native_primary_layout"]
                  if row["filter"] == band}
        counts = " | ".join(str(values[floor]) for floor in coverage["coverage_floors"])
        parts.append(f"| {band} | {counts} |")
    parts += ["", f"The PSF-core experiment completed {psf['frames_completed']} frames. "
              f"The normalization comparison uses {psf['paired_supported_rows']:,} member–amplitude rows "
              "supported in both variants. At 1 nJy arcsec⁻², holding the original nominal input fixed gives:", "",
              "| Filter | Supported pairs | Response change | 95% sky-cell interval |",
              "|---|---:|---:|---:|"]
    for row in psf["unit_amplitude_results"]:
        if row["metric"] == "common_denominator_response_change":
            parts.append(f"| {row['filter']} | {row['n_both']} | {row['median']:+.5f} | "
                         f"[{row['ci_low']:+.5f}, {row['ci_high']:+.5f}] |")
    parts += ["", "Changing the denominator offsets much of the measured recovery change in the "
              "original ratios. These conditional simulations address compact PSF cores. The observational "
              "calibration still requires:", ""]
    parts.extend(f"- {gate}" for gate in report["open_science_gates"])
    command = shlex.join([sys.executable, str(Path(__file__).resolve()), "--current-icl",
                          "--data-root", str(data_root), "--output-dir", str(output) + "_replay"])
    parts += ["", "Replay from the pinned summaries and tables into a fresh directory.", "", "```bash",
              command, "```", ""]
    return "\n".join(parts)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, default=DATA_ROOT)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--current-icl", action="store_true", help="replay current saved ICL and satellite evidence")
    args = parser.parse_args()
    output = args.output_dir.resolve()
    if output.exists():
        raise FileExistsError(f"refusing to overwrite evidence report: {output}")
    data_root = args.data_root.resolve()
    if args.current_icl:
        report = current_icl_report(data_root)
        assessment = render_current_icl(report, data_root, output)
    else:
        payloads, provenance = load_sources(data_root)
        report = summarize(payloads, provenance)
        assessment = render_markdown(report)
    output.mkdir(parents=True)
    (output / "summary.json").write_text(
        json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    (output / "assessment.md").write_text(assessment, encoding="utf-8")
    print(f"status={report['status']} science_measurement={report['science_measurement']}")
    print(f"summary={output / 'summary.json'} sha256={sha256(output / 'summary.json')}")
    print(f"assessment={output / 'assessment.md'} sha256={sha256(output / 'assessment.md')}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
