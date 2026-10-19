#!/usr/bin/env python3
"""Plan independent spatial-block counts for ICL--satellite association tests.

This is a Fisher-z planning approximation, not a re-analysis of the HSC data.
The unit count is the connected spatial block, not the number of cluster rows.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from statistics import NormalDist


HERE = Path(__file__).resolve().parent
DEFAULT_EVIDENCE = HERE / "satellite_icl_spatial_block_assessment_2026-10-19.md"


def two_sided_power(n_units: int, effect_r: float, alpha: float,
                    n_controls: int = 2) -> float:
    """Approximate power for a partial correlation using Fisher's z transform."""
    r = abs(float(effect_r))
    alpha = float(alpha)
    if n_units <= n_controls + 3:
        raise ValueError("independent units must exceed n_controls + 3")
    if not math.isfinite(r) or not 0 < r < 1:
        raise ValueError("effect_r must be finite and strictly between 0 and 1")
    if not math.isfinite(alpha) or not 0 < alpha < 1:
        raise ValueError("alpha must be finite and between 0 and 1")
    if n_controls < 0:
        raise ValueError("n_controls must be non-negative")
    zcrit = NormalDist().inv_cdf(1.0 - alpha / 2.0)
    noncentrality = math.atanh(r) * math.sqrt(n_units - n_controls - 3)
    normal = NormalDist()
    return (1.0 - normal.cdf(zcrit - noncentrality)
            + normal.cdf(-zcrit - noncentrality))


def required_units(effect_r: float, alpha: float, target_power: float,
                   n_controls: int = 2) -> int:
    """Smallest independent-unit count that reaches target approximate power."""
    if not math.isfinite(target_power) or not 0 < target_power < 1:
        raise ValueError("target_power must be finite and between 0 and 1")
    n = n_controls + 4
    while n <= 1_000_000:
        if two_sided_power(n, effect_r, alpha, n_controls) >= target_power:
            return n
        n += 1
    raise ValueError("required sample size exceeds the supported search range")


def plan(observed_blocks: int, observed_r: float, n_comparisons: int,
         n_controls: int = 2, power_targets=(0.8, 0.9),
         effect_sizes=(0.2, 0.293, 0.4)):
    if observed_blocks <= n_controls + 3:
        raise ValueError("observed_blocks must exceed n_controls + 3")
    if n_comparisons < 1:
        raise ValueError("n_comparisons must be positive")
    alpha_nominal = 0.05
    alpha_adjusted = alpha_nominal / n_comparisons
    rows = []
    for effect in effect_sizes:
        rows.append({
            "effect_r_scenario": float(effect),
            "nominal_alpha": alpha_nominal,
            "multiplicity_adjusted_alpha_bonferroni": alpha_adjusted,
            "observed_independent_blocks": int(observed_blocks),
            "observed_power_nominal": two_sided_power(
                observed_blocks, effect, alpha_nominal, n_controls),
            "observed_power_adjusted": two_sided_power(
                observed_blocks, effect, alpha_adjusted, n_controls),
            "required_blocks": {
                f"power_{target:.1f}": {
                    "nominal_alpha": required_units(
                        effect, alpha_nominal, target, n_controls),
                    "adjusted_alpha": required_units(
                        effect, alpha_adjusted, target, n_controls),
                }
                for target in power_targets
            },
        })
    return {
        "schema": "astrafek-icl-spatial-block-power-plan-v1",
        "status": "planning_approximation",
        "science_measurement": False,
        "method": "two-sided Fisher-z approximation for partial correlation; variance uses n_blocks - n_controls - 3",
        "multiplicity_method": "Bonferroni upper bound for the existing family of eligible tests",
        "observed_blocks": int(observed_blocks),
        "observed_effect_r_scenario": float(observed_r),
        "n_controls": int(n_controls),
        "n_comparisons": int(n_comparisons),
        "power_targets": [float(x) for x in power_targets],
        "rows": rows,
        "limitations": [
            "The observed partial correlation is exploratory and is not treated as a confirmed effect.",
            "The Gaussian Fisher-z approximation does not replace the connected-block bootstrap or maxT test.",
            "Required counts refer to independent spatial blocks, not clusters within the current footprint.",
            "Selection, sky transfer, PSF wings, masks, shared pixels, and measurement covariance remain unresolved.",
        ],
    }


def render_markdown(result, evidence_path: Path, evidence_sha256: str) -> str:
    lines = [
        "# ICL–satellite spatial-block power plan",
        "",
        "**Status:** planning approximation; `science_measurement=false`.",
        "",
        "This plan uses the observed exploratory i2 partial-correlation scenario "
        f"r = {result['observed_effect_r_scenario']:.3f}, {result['observed_blocks']} "
        "connected target/control sky blocks, two adjustment variables, and 30 "
        "eligible tests from the existing familywise analysis. It counts spatial "
        "blocks as independent units, not cluster rows.",
        "",
        "The calculations use a two-sided Fisher-z approximation. Bonferroni "
        "adjustment is shown as a conservative planning bound; the actual "
        "analysis must retain the connected-block bootstrap and maxT procedure.",
        "",
        "| Assumed partial r | Current power, α=.05 | Current power, α=.05/30 | "
        "Blocks for 80% power, α=.05/30 | Blocks for 90% power, α=.05/30 |",
        "|---:|---:|---:|---:|---:|",
    ]
    for row in result["rows"]:
        req = row["required_blocks"]
        lines.append(
            f"| {row['effect_r_scenario']:.3f} "
            f"| {row['observed_power_nominal']:.3f} "
            f"| {row['observed_power_adjusted']:.3f} "
            f"| {req['power_0.8']['adjusted_alpha']} "
            f"| {req['power_0.9']['adjusted_alpha']} |"
        )
    lines.extend([
        "",
        "For the observed r = 0.293 scenario, the reported requirement is a count "
        "of independent connected sky blocks. Adding more clusters within the "
        "same blocks does not meet that requirement. These counts are design "
        "guidance only; they do not establish a signal or correct the current "
        "post-selection and measurement-systematic limitations.",
        "",
        f"Evidence report: `{evidence_path}`",
        f"Evidence SHA-256: `{evidence_sha256}`",
        "",
    ])
    return "\n".join(lines)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--evidence", type=Path, default=DEFAULT_EVIDENCE)
    parser.add_argument("--observed-blocks", type=int, default=33)
    parser.add_argument("--observed-r", type=float, default=0.293)
    parser.add_argument("--comparisons", type=int, default=30)
    parser.add_argument("--n-controls", type=int, default=2)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    evidence = args.evidence.resolve()
    if not evidence.is_file():
        parser.error(f"evidence report does not exist: {evidence}")
    evidence_sha = hashlib.sha256(evidence.read_bytes()).hexdigest()
    result = plan(args.observed_blocks, args.observed_r, args.comparisons,
                  args.n_controls)
    result["evidence_report"] = str(evidence)
    result["evidence_report_sha256"] = evidence_sha
    markdown = render_markdown(result, evidence, evidence_sha)
    if args.output:
        output = args.output.resolve()
        if output.exists():
            parser.error(f"refusing to overwrite existing output: {output}")
        output.mkdir(parents=True)
        (output / "summary.json").write_text(
            json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        (output / "power_plan.md").write_text(markdown, encoding="utf-8")
    else:
        print(markdown, end="")


if __name__ == "__main__":
    main()
