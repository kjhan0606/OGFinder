# Variance terms in the ICL–satellite support model (2026-11-06)

## Question and inputs

Satellite-profile support is available for 94 of 224 cluster/filter pairs with
an ICL measurement. This sensitivity analysis checks whether the saved
target/control delivered VARIANCE-plane ratio and a coarse pair-depth proxy
predict that support, and whether adding them to the inverse-probability
weighting (IPW) model changes the ICL–satellite association.

The 224 ICL-supported pairs were joined at 200–300 proper kpc to the existing
full-sample variance table by cluster, filter, and matched-control ID. The
outcome association is calculated only on the 94 pairs that also have a
50–75 kpc satellite profile. Every model includes filter, redshift, log
richness, ICL annular coverage, and center-mask state. The second adds the
standardized log target/control variance ratio; the third also adds mean log
target/control variance as a coarse depth proxy. Each propensity model is
refit in 2,000 shared resamples of connected 5-degree target/control sky
blocks. The same selected pairs and block draws are used across models. No
image pixels were read and no HSC query was submitted.

## Results

| Support model | `i` IPW partial r (95% block interval) | `i2` IPW partial r (95% block interval) | Change vs baseline, `i2` (paired 95% interval) |
|---|---:|---:|---:|
| filter + z + log richness + ICL coverage + center mask | −0.011 (−0.511, 0.586) | 0.227 (−0.035, 0.419) | — |
| baseline + log variance ratio | −0.015 (−0.526, 0.579) | 0.266 (−0.026, 0.470) | +0.039 (−0.069, 0.162) |
| baseline + variance ratio + mean log variance | −0.083 (−0.582, 0.595) | 0.193 (−0.048, 0.401) | −0.034 (−0.163, 0.092) |

Both variance-augmented `i2` paired-change intervals contain zero. The ratio
alone moves the estimate upward; adding the coarse depth proxy moves it
downward relative to baseline. The data therefore do not identify a stable
variance-related correction to the support-adjusted association.

The richer model also reduces the `i` effective sample size from 22.0 to 13.5
of 26 selected pairs and raises its maximum stabilized weight from 2.31 to
5.08. This indicates weak overlap for the small `i` sample once the extra
variance covariates are included. For `i2`, the effective sample size is
59.9 at baseline, 57.7 with variance ratio, and 49.8 with both variance
terms, from 68 selected pairs. The `i` IPW intervals remain very broad.

## Interpretation limits

These are model-dependent sensitivity estimates, not a completeness
correction or causal result. The assumption that satellite-profile support
is independent of the unobserved satellite signal conditional on the measured
covariates remains unverified. The delivered variance plane omits the
effective-gain/source-Poisson contribution and is not a calibrated total-noise
model. Mean log variance is only a coarse depth proxy; neither variance term
represents ICL estimator covariance, coadd weights, sky transfer, extended
PSF wings, or detection completeness. Connected 5-degree groups are
approximate sky blocks. Keep `science_measurement=false`.

## Reproduction and hashes

Astrafek runner:
[`analyze_icl_satellite_ipw_variance_sensitivity.py`](analyze_icl_satellite_ipw_variance_sensitivity.py),
SHA-256 `932fcf912aeb532dd90bff85807f499b94f586b674fac89b8a9c6b6e2a297339`.

Results:
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_ipw_variance_sensitivity_2026-11-06_01/`.
`summary.json` SHA-256
`a1e2972d728bb88590dff3f3ff5197e1f10ea3bbeeec7a1aa0d93f6718fc6b2c`;
`support_variance_join.tsv` SHA-256
`4e5d8983b3963cefae36f2b62b959435838dda16afd456a6fc5c18b6ff471374`.
There were 49 approximate connected blocks (largest: 17 of 224 pairs); all
three propensity models passed their stationarity checks in the full sample.
Of 2,000 block replicates, `i` returned 1,997 valid estimates per model and
`i2` returned 2,000. No new HSC data was requested.
