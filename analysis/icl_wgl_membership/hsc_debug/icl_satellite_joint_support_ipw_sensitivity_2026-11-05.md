# ICL–satellite joint-support weighting sensitivity (2026-11-05)

## Aim

The existing 800-cluster analysis has 224 cluster/filter pairs with usable
200–300 kpc ICL support. A 50–75 kpc satellite profile is available for 94 of
those pairs (26 in `i`, 68 in `i2`). This offline sensitivity asks whether the
observed association within those 94 pairs changes after reweighting for
measured differences between profile-supported and unsupported pairs.

## Method

The selection outcome is whether the saved 50–75 kpc satellite profile exists.
For all 224 pairs, a pooled logistic model predicts that outcome from filter,
redshift, log richness, ICL 200–300 kpc coverage, and the existing central-mask
flag. Each predictor is available for both selected and unselected pairs. The
selected sample receives stabilized inverse-probability weights. The reported
partial correlation adjusts linearly for redshift and richness, matching the
existing 94-pair analysis. A log-richness adjustment is included as a secondary
check.

Uncertainty comes from 2,000 resamples of connected 5-degree target/control
sky blocks, refitting the selection model in every resample. The original
covariance-adjusted, unweighted estimates reproduce the saved result. The
weighting model assumes that profile availability is conditionally independent
of unmeasured satellite signal given the modeled covariates. That assumption
cannot be tested from these data, so this is a sensitivity analysis rather than
a completeness correction.

## Results

| Filter | Available / ICL-supported | Raw partial `r` | IPW partial `r` | IPW spatial-block 95% interval | IPW effective N |
|---|---:|---:|---:|---:|---:|
| `i` | 26 / 71 | 0.079 | −0.011 | [−0.566, 0.581] | 22.0 / 26 |
| `i2` | 68 / 153 | 0.293 | 0.227 | [−0.034, 0.416] | 59.9 / 68 |

The association weakens after weighting in both filters. In `i2`, it decreases
by 0.065 but the interval still reaches zero. The `i` estimate is near zero
with broad uncertainty. The two filters therefore do not provide a robust
selection-adjusted detection.

The weighted point estimates were stable across three prespecified model
alternatives: omitting central-mask state, adding a quadratic coverage term,
and allowing filter-specific slopes. Across these alternatives, `i` estimates
range from −0.014 to +0.021, and `i2` estimates range from +0.226 to +0.231.
Effective sample sizes remain 21.8–24.2 of 26 in `i` and 56.7–60.0 of 68 in
`i2`. This reduces concern that the point shift is driven only by the primary
logistic specification; it does not address unmeasured profile-selection
causes or change the broad uncertainty intervals.

The profile-supported subset is observably selected. Relative to unsupported
pairs, the 200–300 kpc ICL coverage standardized mean difference is +0.44 in `i`
and +0.70 in `i2`; richness differences are −0.29 and −0.34, and the redshift
difference is +0.05 and +0.32. Central-mask fractions are nearly identical.
The pooled model's odds ratio per one-standard-deviation change is 1.86 for ICL
coverage, 0.80 for log richness, and 1.27 for redshift. Predicted probabilities
span 0.103–0.870 across all pairs and 0.151–0.870 in the selected set. Maximum
stabilized weights are 2.31 (`i`) and 1.98 (`i2`); effective sample size remains
about 85% and 88% of the observed pairs. This suggests reasonable overlap for
these measured covariates, while leaving unmeasured selection unresolved.

The first runner version counted SciPy BFGS precision-loss messages as failed
fits even when the score gradient met the predeclared stationarity tolerance.
The corrected runner records optimizer status, gradient, and coefficient
diagnostics in the JSON output. All 216 precision-loss draws met the
`2×10⁻⁵` stationarity threshold; the maximum infinity-norm gradient was
`1.40×10⁻⁶` and the largest absolute coefficient was `3.63`. No selection
model draw failed stationarity. The corrected run produced 1,994 valid `i` and
2,000 valid `i2` draws; the six omitted `i` draws had fewer than six supported
pairs. Estimates are nearly unchanged and both intervals still include zero.

## Interpretation and next step

Measured support differences can account for part of the raw `i2` association,
but the weighted interval remains broad and includes zero. The weighting does
not establish that the other 130 ICL-supported pairs would have the same
satellite signal if their profiles were measurable. Do not use these estimates
as final science measurements, detections, or population completeness
corrections (`science_measurement=false`).

Next, validate the selection mechanism with injected satellites spanning the
224-pair sample, including the masked and low-coverage cases. Only after
recovery tests establish an appropriate missingness model should a corrected
population association be considered. Expanding the HSC sample remains
necessary for stronger precision.

## Reproduction

Run from the Astrafek OGFinder environment:

```bash
python3 analysis/icl_wgl_membership/hsc_debug/analyze_icl_satellite_joint_support_ipw_sensitivity.py --replicas 2000
```

Astrafek runner:
`analyze_icl_satellite_joint_support_ipw_sensitivity.py` (SHA-256
`00131be3d6f91246d7e1733ee2e06845525359dd183d8464b87db26d822deb85`).
Output directory:
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_joint_support_ipw_2026-11-05_06/`.
Summary SHA-256 `8dbfae2169361c9fb6ce5ccfcce684b486162dbf4820e531112cdfc4f39fdbc7`;
pair propensity table SHA-256
`912571de06313fa03ba08aea35c5d8edd51f15080115087a86c069f721def4e6`.
The analysis used no HSC network query or new download.
