#!/usr/bin/env python3
"""Build a hash-provenance ICL--satellite robustness report from saved audits.

This is a local summary step: it reads existing JSON diagnostics only, makes
no HSC requests, and does not promote exploratory measurements to science
results.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path


DEFAULT_DATA_ROOT = Path("/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13")
SOURCES = {
    "overlap": "icl_satellite_overlap_sensitivity_2026-11-03_01/summary.json",
    "pixel_noise": "icl_satellite_pixel_noise_bias_2026-11-03_01/summary.json",
    "outer_sky": "icl_outer_sky_association_reblocked_2026-11-03_01/summary.json",
    "joint_sky": "icl_satellite_joint_sky_covariance_2026-11-03_01/summary.json",
    "poisson_audit": "icl_poisson_noise_input_audit_2026-11-04_01/summary.json",
    "supported_94": "icl_satellite_supported_94_association_2026-11-05_01/summary.json",
    "selection": "icl_satellite_association_selection_2026-10-20_03/summary.json",
    "sample_size": "icl_sample_size_plan_2026-10-21_01/sample_size_plan.json",
    "catalog_expansion": "icl_s23b_expansion_feasibility_2026-10-21_02/summary.json",
    "support_ipw": "icl_satellite_joint_support_ipw_2026-11-05_06/summary.json",
    "psf_sample": "icl_joint_injection_psf_sample_representativeness_2026-11-05_01/summary.json",
    "frame_count": "icl_frame_count_sensitivity_2026-11-05_01/summary.json",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def load_sources(data_root: Path) -> tuple[dict[str, dict], dict[str, str]]:
    payloads, hashes = {}, {}
    for name, relative in SOURCES.items():
        path = data_root / relative
        if not path.is_file():
            raise FileNotFoundError(f"required diagnostic summary is missing: {path}")
        payload = json.loads(path.read_text(encoding="utf-8"))
        if payload.get("science_measurement") is not False:
            raise ValueError(f"{name} diagnostic has an unexpected science gate")
        payloads[name] = payload
        hashes[str(path)] = sha256(path)
    return payloads, hashes


def fmt(value: float) -> str:
    return f"{float(value):.3f}"


def interval(values: list[float]) -> str:
    if len(values) == 2:
        return f"[{fmt(values[0])}, {fmt(values[1])}]"
    return f"[{fmt(values[0])}, {fmt(values[2])}] (median {fmt(values[1])})"


def build_report(data_root: Path, output: Path) -> None:
    audits, input_hashes = load_sources(data_root)
    overlap = audits["overlap"]
    pixel = audits["pixel_noise"]
    outer = audits["outer_sky"]
    joint = audits["joint_sky"]
    poisson = audits["poisson_audit"]
    supported = audits["supported_94"]
    selection = audits["selection"]
    sample_size = audits["sample_size"]
    expansion = audits["catalog_expansion"]
    ipw = audits["support_ipw"]
    psf_sample = audits["psf_sample"]
    frame_count = audits["frame_count"]

    lines = [
        "# ICL--satellite association: evidence status",
        "",
        "This report consolidates saved exploratory diagnostics for the fixed",
        "ICL--satellite association sample. It does not reread image pixels or",
        "submit HSC queries. All component summaries retain",
        "`science_measurement=false`; the numerical intervals below are",
        "diagnostic and do not establish a detection or causal relation.",
        "",
        "## Main robustness comparisons",
        "",
        "The first table uses the same pair sample and z/richness-adjusted",
        "partial-correlation estimator for the raw and finite-lag pixel-noise",
        "results. The shell-exclusion values are a paired aperture-overlap",
        "sensitivity: supported satellite 50--75 kpc annulus pixels are removed",
        "from the target ICL 200--300 kpc ring; matched-control ICL and satellite",
        "measurements remain fixed.",
        "",
        "| Filter | Pairs | Raw partial r (95% block interval) | Finite-lag pixel-noise corrected r (95% interval) | After shell-pixel exclusion (95% interval) | Paired change, excluded minus raw (95% interval) |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for band in ("i", "i2"):
        raw = pixel["filters"][band]
        shell = overlap["correlations"][band]
        lines.append(
            f"| {band} | {raw['n_pairs']} | "
            f"{fmt(raw['point_estimate']['raw_partial_r'])} "
            f"{interval(raw['raw_partial_r_bootstrap_2p5_50_97p5'])} | "
            f"{fmt(raw['point_estimate']['pixel_noise_corrected_partial_r'])} "
            f"{interval(raw['noise_corrected_partial_r_bootstrap_2p5_50_97p5'])} | "
            f"{fmt(shell['satellite_shells_masked_partial_r_z_richness_adjusted'])} "
            f"{interval(shell['block_bootstrap']['satellite_shells_masked_partial_r_2p5_50_97p5'])} | "
            f"{fmt(shell['point_change_masked_minus_original'])} "
            f"{interval(shell['block_bootstrap']['paired_change_masked_minus_original_2p5_50_97p5'])} |"
        )

    lines.extend([
        "",
        "## Support-selection and measurement-sample checks",
        "",
        "The 94 satellite-profile pairs are a subset of 224 ICL-supported",
        "cluster/filter pairs. An exploratory stabilized inverse-probability",
        "analysis reweighted the 94 rows using filter, redshift, log richness,",
        "ICL aperture coverage, and the central-mask flag. The model assumes",
        "selection is conditionally independent of the unobserved satellite",
        "signal; that assumption is unverified.",
        "",
        "| Filter | Supported / ICL-supported | Raw partial r | Selection-weighted partial r (95% block interval) | Effective N |",
        "|---|---:|---:|---:|---:|",
    ])
    for band in ("i", "i2"):
        row = ipw["association_by_filter"][band]
        lines.append(
            f"| {band} | {row['n_satellite_supported']} / {row['n_icl_supported']} | "
            f"{fmt(row['raw_partial_r_z_richness'])} | "
            f"{fmt(row['ipw_partial_r_z_richness'])} "
            f"{interval(ipw['block_bootstrap']['by_filter'][band]['ipw_partial_r_ci95'])} | "
            f"{fmt(row['effective_n'])} / {row['n_satellite_supported']} |"
        )
    lines.extend([
        "",
        "The i2 estimate falls from 0.293 to 0.227 after measured-covariate",
        "weighting; the i estimate moves from 0.079 to −0.011. The associated",
        "block intervals remain broad and include zero. Stable point estimates",
        "across the stored alternative propensity specifications do not address",
        "unmeasured reasons that satellite profiles are missing.",
        "",
        "A separate representativeness audit shows that the 57 PSF-convolved",
        "joint-injection pairs are concentrated in well-covered data: their",
        "200–300 kpc ICL coverage is higher than the 224-pair parent by 0.800 SD",
        "overall (0.736 SD in i2), and 54/57 pass the 60% satellite-support cut",
        "versus 98/224 parent pairs. Injection transfer from this subset cannot",
        "be generalized to the full supported sample.",
        "",
        "The count-only frame audit finds no clear relation between target-minus-",
        "control frame/visit/CCD count offsets and ICL contrast in 224 supported",
        "pairs. This rules out only a simple integer-count imbalance explanation;",
        "it does not test pixel depth, coadd weights, sky residuals, or variance.",
        "",
        "## Selection and sample-size limits",
        "",
        "The second table treats the measured 600--900 kpc target-minus-control",
        "plane intercept as a residual-sky diagnostic covariate. This annulus",
        "can contain real cluster light and is not subtracted from the ICL",
        "measurement.",
        "",
        "| Filter | Pairs | Baseline partial r (z, richness) | With outer-residual covariate | Paired change (95% block interval) |",
        "|---|---:|---:|---:|---:|",
    ])
    for band in ("i", "i2"):
        row = outer["association_outer_residual_sensitivity_by_filter"][band]
        lines.append(
            f"| {band} | {row['n_pairs']} | {fmt(row['baseline_partial_r_z_richness'])} | "
            f"{fmt(row['outer_residual_adjusted_partial_r_z_richness_sky'])} | "
            f"{interval(row['paired_delta_ci95'])} |"
        )

    lines.extend([
        "",
        "## Full existing joint-support sample sensitivity",
        "",
        "The earlier correlation tables use a 74-pair diagnostic subset. A",
        "separate Astrafek analysis recomputed the same raw ICL and satellite",
        "partial-correlation estimator on all 94 pairs with both ICL and",
        "satellite-profile support. The 94-pair subset remains selected from the",
        "800-cluster parent sample. The paired bootstrap compares the nested",
        "74-pair and 94-pair estimates under identical 5-degree connected-sky",
        "block resampling.",
        "",
        "| Filter | Support sample | Pairs | Partial r (95% connected-block interval) |",
        "|---|---|---:|---:|",
    ])
    for band in ("i", "i2"):
        result = supported["results_by_filter"][band]
        for sample_key, sample_label in (("previous74", "nested diagnostic"),
                                         ("support94", "full joint support")):
            row = result[sample_key]
            lines.append(
                f"| {band} | {sample_label} | {row['n_pairs']} | "
                f"{fmt(row['partial_r_z_richness'])} "
                f"{interval(row['block_bootstrap_ci95'])} |"
            )
        change = result["paired_scope_change"]
        lines.append(
            f"| {band} | paired 94 minus 74 change | — | "
            f"{fmt(change['scope_delta_full94_minus_narrow74'])} "
            f"{interval(change['paired_delta_block_bootstrap_ci95'])} |"
        )
    lines.extend([
        "",
        "Both paired scope-change intervals include zero. The i2 point estimate",
        "rises from 0.234 to 0.293, but the paired interval does not establish a",
        "change in the association. The i estimate shifts from 0.120 to 0.079",
        "with a paired interval that also includes zero. The expansion addresses",
        "sample accounting within existing support and does not correct support",
        "selection or image-level systematics.",
        "",
        "## Selection and sample-size limits",
        "",
        "A separate support-selection audit has 224 ICL-supported rows; 94 also",
        "have the required satellite-profile support (26 in i, 68 in i2). The",
        "finite-lag, shell-exclusion, and outer-sky tables use the narrower",
        "74-pair sample (21 in i, 53 in i2). The joint-support table reports",
        "both the nested 74-pair subset and all 94 supported pairs. The gates",
        "differ and their row counts should not be interchanged.",
        "",
        "The supported rows are selected: in i2, association rows differ from",
        "the non-associated ICL-supported rows by +1.055 SD in target image",
        "coverage, +0.318 SD in redshift, and -0.325 SD in richness. The audit",
        "describes selection; it does not correct it or make the sample",
        "representative of all 800 parent clusters.",
        "",
        "The following are design scenarios, not observed-effect claims. They",
        "assume 80% power, two covariates, and conservative Bonferroni control",
        "for 168 hypotheses. Parent-cluster projections hold the observed i2",
        "support yield fixed; design effects are sensitivity assumptions.",
        "",
        "| Assumed partial r | Design effect | Required i2 usable rows | Projected parent clusters |",
        "|---:|---:|---:|---:|",
    ])
    planned = sample_size["scenarios"]
    for effect in (0.20, 0.293):
        for design_effect in (1.0, 1.5):
            scenario = next(row for row in planned
                            if row["planning_partial_r"] == effect
                            and row["target_power"] == 0.8
                            and row["multiplicity_scenario"] == "bonferroni_168"
                            and row["design_effect"] == design_effect)
            lines.append(
                f"| {fmt(effect)} | {fmt(design_effect)} | "
                f"{scenario['inflated_usable_rows']:,} | "
                f"{scenario['projected_parent_lenses_i2']:,} |"
            )
    projected_i2 = expansion["projected_total_support_rows_if_yield_transfers"]["i2"]
    lines.extend([
        "",
        f"The offline S23B catalog audit gives at most "
        f"{expansion['offline_new_candidate_upper_bound']:,} additional catalog",
        f"candidates ({expansion['catalog_only_parent_upper_bound']:,} parent",
        f"clusters total), projecting about {fmt(projected_i2)} i2 support rows",
        "if the current yield transferred. This remains below every displayed",
        "family-wise planning target. It is a catalog-only upper bound: candidate",
        "PDR3 images, masks, depth, PSF, footprint, and matched controls are not",
        "verified, and no HSC query or new image download is represented here.",
        "The r=0.293 row is a post-selection exploratory estimate and is not a",
        "confirmatory effect-size assumption.",
    ])

    pairs = joint.get("n_pairs_by_filter", {})
    noise_groups = poisson.get("groups", {})
    lines.extend([
        "",
        "## Interpretation and limits",
        "",
        f"- The earlier robustness analysis has {joint.get('n_pairs', 'unknown')} pairs "
        f"({pairs.get('i', 'unknown')} in i and {pairs.get('i2', 'unknown')} in i2).",
        "  The expanded support analysis has 94 pairs, with 26 in i and 68 in i2.",
        "- The i2 point estimate is positive in these diagnostics, but the",
        "  supported-shell exclusion attenuates it; the paired change interval",
        "  is an exploratory sensitivity with only 28 spatial blocks. The i-band",
        "  estimate is much less precise (21 pairs; 13 blocks).",
        "- The finite 9x9 lag-kernel correction covers only the supplied",
        "  delivered-variance noise model. It is not a correction for large-scale",
        "  sky modes, PSF wings, selection, completeness, extinction, or calibration.",
        "- The shell-pixel test is not a deblending or PSF correction. It removes",
        "  supported pixels from the ICL aperture while keeping the satellite",
        "  estimator and matched-control ICL fixed.",
        "- The outer annulus is a residual diagnostic, not clean sky; no flux was",
        "  subtracted. Its paired change interval includes zero in both filters.",
        "- The coadd-header audit found no gain or exposure metadata in the 1,598",
        "  target/control images. A physically based injected-source Poisson",
        "  variance and completeness correction therefore remains unavailable",
        "  from the currently saved coadds and variance planes alone.",
        "- These results support a tentative i2 association worth following up,",
        "  not a confirmed ICL--satellite relation. Connected HSC fields, PSF",
        "  wings, sky transfer, source completeness, and noise covariance need",
        "  stronger validation before a science claim.",
        "",
        "## Provenance",
        "",
        "Bootstrap intervals use 5,000 replicates in their component analyses.",
        "Source-summary SHA-256 values:",
        "",
    ])
    for path, digest in sorted(input_hashes.items()):
        lines.append(f"- `{path}`: `{digest}`")
    script_path = Path(__file__).resolve()
    lines.extend([
        f"- report builder `{script_path}`: `{sha256(script_path)}`",
        "",
        "This consolidation reads saved JSON summaries only; it makes no network",
        "request, does not alter image data, and leaves the science gate false.",
        "",
    ])
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("x", encoding="utf-8") as stream:
        stream.write("\n".join(lines))
    print(json.dumps({
        "status": "complete_exploratory_synthesis",
        "science_measurement": False,
        "output": str(output.resolve()),
        "input_summaries": len(input_hashes),
    }, ensure_ascii=False))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, default=DEFAULT_DATA_ROOT)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    build_report(args.data_root.resolve(), args.output.resolve())


if __name__ == "__main__":
    main()
