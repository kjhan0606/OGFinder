# 2026-11-06: covariate-adjusted ICL difference by forced-profile support

## Question and estimand

Within connected 5-degree blocks that contain both selected and excluded
cluster/filter pairs under the post-hoc 20/20 coverage rule, how much does the
ICL 200–300 kpc target-minus-control residual differ by profile support after
adjusting for redshift and log CAMIRA richness?

## Method

The analysis fits a pair-weighted within-block linear model with predictors
`selected_20_20`, redshift, and log richness. Demeaning all terms within each
connected target/control block is equivalent to including block fixed effects.
Uncertainty is summarized with 10,000 percentile-bootstrap resamples of the
mixed blocks. The estimand is limited to observed mixed-block support; it does
not transport to blocks with only selected or only excluded pairs.

## Results

| Filter | Mixed blocks | Pairs (selected / excluded) | Adjusted selected-minus-excluded ICL residual (nJy arcsec⁻²) | 95% block-bootstrap interval |
|---|---:|---:|---:|---:|
| i | 9 | 40 / 19 | +1.47 | [−2.10, +4.79] |
| i2 | 21 | 77 / 40 | +2.13 | [−0.64, +5.91] |

Both adjusted intervals include zero. The positive i2 contrast from the
unadjusted equal-block summary is attenuated after redshift and richness
adjustment, and its interval includes zero. The evidence therefore does not
resolve whether 20/20 support selects systems with different ICL residuals
after these covariates are controlled.

## Limitations

- The 20% threshold is post-hoc, and this analysis is not corrected for trying
  multiple support rules or model forms.
- Coverage selection is deterministic; the coefficient is descriptive, not a
  causal effect or a correction for missing profiles.
- Only 9 `i` and 21 `i2` mixed blocks contribute. The 95% intervals do not
  include ICL measurement uncertainty.
- Forced CAMIRA positions without unique SEP detections are not confirmed
  isolated galaxies.
- No correction or science promotion follows; `science_measurement=false`.

## Reproduction

Runner: [`estimate_icl_forced_coverage_selection_adjusted.py`](estimate_icl_forced_coverage_selection_adjusted.py), SHA-256 `5f81094c8a8f4a0b37567a6c05cc9e2b482fcd786588a308ae734b3bb2530cf0`.

Output: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_forced_coverage_selection_adjusted_2026-11-06_01/summary.json`
(SHA-256 `7e3c57cbd95b32986a7f52b707429f682a6d4e358e12724f9866f14f69b6f07c`);
mixed-block table SHA-256 `d7ed0601f1154d0ed359c39abbc70532d4322e77775abc33357686447da6de3c`.
The runner uses local saved tables and manifests only; no HSC query or image
pixel read occurred.
