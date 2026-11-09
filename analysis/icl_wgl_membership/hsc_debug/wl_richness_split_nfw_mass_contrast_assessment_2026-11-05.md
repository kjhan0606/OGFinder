# Richness-split NFW mass-contrast screen (2026-11-05)

Status: exploratory model comparison; `science_measurement=false`.

## Data and fit

The fit uses the matched-random, redshift-matched high-minus-low richness
`DeltaSigma` profile from 204 lens pairs and its shared 50-block spatial
jackknife covariance. It fits two group masses simultaneously to the four
tangential radial bins from 0.2 to 1.5 proper Mpc. For each group, the NFW
profile is averaged over the exact lens-source rows and the same HSC estimator
weights used by the stack. The model has one effective `M200c` per richness
group and a fixed concentration; concentrations 3, 4, and 5 are evaluated
separately. The covariance is the covariance of the observed high-minus-low
contrast, so no independent-profile covariance approximation is added.

The all-bin fit uses source bins 1–4. A second fit uses the same lens pairs,
random centers, and covariance procedure but excludes source bin 1 from both
target and random estimators. Both fits use weak shear. Conditional profile
intervals are one-parameter `Delta chi2=1` intervals after profiling over the
other group mass. The two-mass fit has two degrees of freedom.

## Results

Masses are effective single-halo stack scales, not arithmetic mean masses.
The ratio interval is a conditional profile interval at fixed concentration.

| Fixed c | All bins: low / high `M200c` (`10^14 Msun`) | All bins: high/low ratio [profile interval] | Excluding bin 1: low / high `M200c` (`10^14 Msun`) | Excluding bin 1: high/low ratio [profile interval] |
|---:|---:|---:|---:|---:|
| 3 | 1.79 / 3.39 | 1.89 [1.13, 8.63] | 0.70 / 1.98 | 2.85 [1.31, 28.03] |
| 4 | 5.23 / 7.30 | 1.40 [1.08, 2.86] | 2.26 / 3.90 | 1.73 [1.16, 9.86] |
| 5 | 12.36 / 14.95 | 1.21 [1.07, 1.83] | 4.98 / 6.92 | 1.39 [1.10, 3.42] |

The best-fit high-richness mass is above the low-richness mass for all tested
concentrations and both source-bin selections. However, the ratio changes
substantially with concentration and tomographic selection; the profile
interval upper bounds are very broad. Individual mass intervals span several
dex. Conditional equal-mass profile fits worsen by `Delta chi2` about 44 (all
bins) and 51 (source bin 1 excluded), while the best-fit two-mass profile has
`chi2/dof` about 1.6–1.8. These conditional statistics inherit the approximate
jackknife covariance and are not detection significances.

The data therefore support an exploratory higher-mass direction in this
selected sample, but do not measure a reliable richness-mass ratio. The
source-bin-1-excluded result preserves the direction while showing that the
mass scale shifts appreciably. The fixed-concentration and source-bin
sensitivities are material compared with the quoted conditional intervals.

## Interpretation limits

The single NFW mass per richness group does not represent the broad
richness-selected mass distribution. The fit omits concentration variation,
reduced shear, magnification, miscentering, two-halo structure, photo-z and
shear-calibration uncertainties, matched-aperture depth/mask/PSF/sky
differences, and covariance uncertainty. The fit spans only four radial bins,
and the 50 equal-area blocks are not validated contiguous HSC survey regions.
It does not establish a calibrated richness–mass relation or the ICL–satellite
association. Keep `science_measurement=false`.

## Reproducibility

Astrafek runner:
[`fit_wl_richness_split_nfw_mass_contrast.py`](fit_wl_richness_split_nfw_mass_contrast.py),
SHA-256 `25e9c57d59cd3c61568bbbf5478f9726bc31263f49ad442931c7909a5dd22f1d`.

All-source-bin summary:
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_richness_split_nfw_mass_contrast_2026-11-05_04/summary.json`,
SHA-256 `c0e4c92e5ef00fc77ffef7e26e918d3da312049e6d5ee2877edbe4385c32217b`.

Source-bin-1-excluded summary:
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_richness_split_nfw_mass_contrast_2026-11-05_03/summary.json`,
SHA-256 `6d293f41c12a2d44fc4dd7303be238c9bfbcca27570db180f0bf93e3224eaa94`.
The source-bin-excluded stack and covariance are in
`wl_random_subtracted_richness_split_zbin1_excluded_2026-11-05_01/summary.json`
(SHA-256 `814e9162aad25cc6c80c65a19ddfc52b88f1acc6f7d324723253395cea456df6`).
No HSC query or image read was made.
