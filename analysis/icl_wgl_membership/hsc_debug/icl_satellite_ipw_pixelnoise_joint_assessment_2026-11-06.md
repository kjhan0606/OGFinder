# Joint satellite-support and pixel-noise sensitivity (2026-11-06)

## Question and method

The 50–75 kpc satellite profile exists for 94 of 224 ICL-supported
cluster/filter pairs. Two earlier analyses addressed satellite-support
selection and finite-lag pixel-noise covariance separately. This analysis
combines them: it fits the support propensity among all 224 pairs, applies
stabilized inverse-probability weights to the selected 94, and uses the saved
finite-lag ICL/satellite covariance to correct the weighted residual
cross-products.

Three prespecified propensity models are compared: filter, redshift, log
richness, ICL coverage, and center mask; that baseline plus log delivered
target/control variance ratio; and the second model plus mean log target and
control variance. Each model is refit in every one of 2,000 resamples of the
same connected 5-degree target/control blocks. The pixel covariance support
auditor verified exact coverage of the same 94 selected pairs. No HSC query,
download, or image read occurred.

## Results

| Support propensity model | Filter | Noise-corrected IPW partial r | 95% connected-block interval | Effective N |
|---|---|---:|---:|---:|
| Baseline covariates | `i` | −0.011 | −0.579–0.602 | 22.0 / 26 |
| Baseline covariates | `i2` | 0.238 | −0.023–0.446 | 59.9 / 68 |
| + log variance ratio | `i` | −0.015 | −0.590–0.603 | 21.9 / 26 |
| + log variance ratio | `i2` | 0.278 | −0.018–0.495 | 57.7 / 68 |
| + variance ratio + mean log variance | `i` | −0.085 | −0.622–0.620 | 13.5 / 26 |
| + variance ratio + mean log variance | `i2` | 0.201 | −0.054–0.410 | 49.8 / 68 |

Every `i2` interval includes zero. The finite-lag correction alone changes
the baseline IPW point estimate from 0.227 to 0.238. Adding variance ratio
moves the corrected point estimate to 0.278; adding the coarse depth term
moves it to 0.201. Paired intervals for those changes from the baseline are
`[-0.079, 0.167]` and `[-0.174, 0.093]`, so neither change is resolved by the
sample. All `i` intervals are much broader. The fuller support model lowers
the `i` effective sample size to 13.5 pairs and is particularly unstable.

This combined estimate is more conservative about the support subset than
the unweighted pixel-noise-corrected `i2` value (0.306), but it is not a
population correction: the unmeasured-selection assumption remains
unverified.

## Limits and science status

The finite-lag adjustment treats the saved stationary covariance estimates
as fixed and does not propagate their estimation uncertainty. It excludes
source Poisson/effective gain, residual-sky transfer, extended PSF wings,
calibration, and astrophysical covariance. IPW requires satellite-profile
availability to be independent of the unobserved satellite signal given the
modeled covariates; this missing-at-random condition is not demonstrated.
The delivered variance plane is not total noise. The 5-degree connected
components are approximate blocks, not validated contiguous HSC regions.
These remain selected-sample exploratory associations; keep
`science_measurement=false`.

## Reproduction and integrity

Runner: [`analyze_icl_satellite_ipw_pixelnoise_joint.py`](analyze_icl_satellite_ipw_pixelnoise_joint.py),
SHA-256 `493d32b165bf1a2166623135261697d3082f1770a101bd11d9fbd651ff2a7d5d`.
The upstream support/variance propensity table SHA-256 is
`4e5d8983b3963cefae36f2b62b959435838dda16afd456a6fc5c18b6ff471374`; the
pixel covariance input SHA-256 is
`cc85ffd59fc94154a6597d3feca0e4d811682b9fcda7b79169e27ab778765167`.
The independent covariance-scope verification reports all 94 selected pairs
covered and no missing pairs.

Output directory:
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_ipw_pixelnoise_joint_2026-11-06_01/`.
`summary.json` SHA-256
`761e3448fd5d01b2937c62cc3668bc8a2f7e3af0bcdc74a3a2c821c193d3d01d`.
