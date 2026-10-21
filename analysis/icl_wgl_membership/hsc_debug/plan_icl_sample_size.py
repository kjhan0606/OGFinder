#!/usr/bin/env python3
"""Plan ICL–satellite sample sizes from partial-correlation scenarios.

This is an exploratory design calculation, not a post-hoc power claim. It uses
the selected 800-cluster sample's current supported-row yield only to translate
required usable rows into a rough parent-catalogue size. It never queries HSC.
"""
from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
import sys

from scipy.stats import nct, t


def _power(n: int, effect: float, alpha: float, covariates: int) -> float:
    df = n - covariates - 2
    if df <= 0:
        return 0.0
    threshold = t.ppf(1.0 - alpha / 2.0, df)
    noncentrality = effect * math.sqrt(df / (1.0 - effect * effect))
    return float(nct.cdf(-threshold, df, noncentrality)
                 + nct.sf(threshold, df, noncentrality))


def required_independent_rows(effect: float, desired_power: float, alpha: float,
                              covariates: int = 2) -> int:
    """Smallest row count reaching power for a two-sided partial-correlation t test."""
    if not 0 < abs(effect) < 1:
        raise ValueError("effect must be a partial correlation strictly between -1 and 1")
    if not 0 < desired_power < 1 or not 0 < alpha < 1:
        raise ValueError("power and alpha must lie strictly between 0 and 1")
    if covariates < 0:
        raise ValueError("covariates must be non-negative")
    lo = max(4 + covariates, 5)
    hi = lo
    while _power(hi, abs(effect), alpha, covariates) < desired_power:
        hi *= 2
        if hi > 10_000_000:
            raise ValueError("required sample size exceeds the supported search range")
    while lo < hi:
        mid = (lo + hi) // 2
        if _power(mid, abs(effect), alpha, covariates) >= desired_power:
            hi = mid
        else:
            lo = mid + 1
    return lo


def build_plan(*, parent_lenses: int, usable_i: int, usable_i2: int,
               family_size: int = 168,
               effects=(0.10, 0.20, 0.293), powers=(0.80, 0.90),
               design_effects=(1.0, 1.25, 1.5, 2.0), covariates: int = 2):
    if parent_lenses <= 0 or family_size <= 0:
        raise ValueError("parent_lenses and family_size must be positive")
    if not 0 < usable_i <= parent_lenses or not 0 < usable_i2 <= parent_lenses:
        raise ValueError("usable row counts must be positive and no greater than parent_lenses")
    if any(not 0 < effect < 1 for effect in effects):
        raise ValueError("all planning effect sizes must lie between 0 and 1")
    if any(effect < 1 for effect in design_effects):
        raise ValueError("design effects must be at least 1")

    family_alpha = 0.05 / family_size
    yields = {"i": usable_i / parent_lenses, "i2": usable_i2 / parent_lenses}
    rows = []
    for effect in effects:
        for power in powers:
            for multiplicity, alpha in (("single_test_0.05", 0.05),
                                        (f"bonferroni_{family_size}", family_alpha)):
                independent = required_independent_rows(effect, power, alpha, covariates)
                for deff in design_effects:
                    usable_required = math.ceil(independent * deff)
                    rows.append({
                        "planning_partial_r": effect,
                        "target_power": power,
                        "multiplicity_scenario": multiplicity,
                        "alpha_two_sided": alpha,
                        "design_effect": deff,
                        "independent_usable_rows": independent,
                        "inflated_usable_rows": usable_required,
                        "projected_parent_lenses_i": math.ceil(usable_required / yields["i"]),
                        "projected_parent_lenses_i2": math.ceil(usable_required / yields["i2"]),
                    })
    return {
        "schema": "astrafek-icl-satellite-sample-size-plan-v1",
        "status": "exploratory_design_sensitivity",
        "science_measurement": False,
        "input": {
            "current_parent_lenses": parent_lenses,
            "current_joint_support_rows": {"i": usable_i, "i2": usable_i2},
            "observed_support_yield": yields,
            "family_size": family_size,
            "covariates": covariates,
            "effect_sizes_are_assumptions_not_confirmatory_estimates": True,
        },
        "method": {
            "test": "two-sided noncentral-t power for partial correlation, controlling the stated number of covariates",
            "multiplicity": "single-test alpha 0.05 and conservative Bonferroni alpha 0.05/family_size",
            "dependence": "design-effect multipliers are sensitivity scenarios, not estimated from these data",
            "parent_projection": "current supported-row yield assumed constant; i and i2 projected separately",
        },
        "scenarios": rows,
        "limitations": [
            "The endpoint correlation was selected during exploratory analysis; it is not an unbiased effect-size estimate.",
            "The current supported-row yield reflects image, mask, filter, and joint ICL/satellite selection and may not transfer to new clusters.",
            "A design effect cannot fix non-random missingness or PSF, sky, completeness, calibration, or coadd-covariance systematics.",
            "The calculation is a planning sensitivity, not evidence that a signal exists or that a larger sample would make it measurable.",
        ],
    }


def _floats(value: str):
    return tuple(float(item) for item in value.split(",") if item.strip())


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--parent-lenses", type=int, default=800)
    parser.add_argument("--usable-i", type=int, default=26)
    parser.add_argument("--usable-i2", type=int, default=68)
    parser.add_argument("--family-size", type=int, default=168)
    parser.add_argument("--effects", type=_floats, default=(0.10, 0.20, 0.293))
    parser.add_argument("--powers", type=_floats, default=(0.80, 0.90))
    parser.add_argument("--design-effects", type=_floats, default=(1.0, 1.25, 1.5, 2.0))
    parser.add_argument("--covariates", type=int, default=2)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args(argv)
    plan = build_plan(parent_lenses=args.parent_lenses, usable_i=args.usable_i,
                      usable_i2=args.usable_i2, family_size=args.family_size,
                      effects=args.effects, powers=args.powers,
                      design_effects=args.design_effects, covariates=args.covariates)
    rendered = json.dumps(plan, indent=2, allow_nan=False) + "\n"
    if args.output:
        if args.output.exists():
            parser.error(f"refusing to overwrite {args.output}")
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
        print(f"Wrote exploratory design plan: {args.output}")
    else:
        sys.stdout.write(rendered)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
