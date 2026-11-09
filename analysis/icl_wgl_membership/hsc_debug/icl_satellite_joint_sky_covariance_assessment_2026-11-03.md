# Joint sky-residual and coadd-covariance sensitivity — 2026-11-03

## Question and scope

Does the existing ICL–satellite surface-brightness association change when both
the paired outer-residual indicator and the paired local coadd-correlation
proxy are included with redshift and richness? This is an offline sensitivity
analysis of the fixed 74-pair association sample (53 `i2`, 21 `i`). It uses the
existing 600–900 proper-kpc target-minus-control plane intercept and a local
approximate top-hat sigma-inflation factor from the variance-normalized 9×9
pixel-lag covariance grid. No images were reread or modified, and no HSC query
or download was made.

For each filter, partial Pearson correlation between the paired target-minus-
control ICL contrast (200–300 kpc) and satellite contrast (50–75 kpc) was
calculated after nuisance adjustment in four models: redshift and richness;
those plus outer residual; those plus the target-minus-control coadd sigma
proxy; and all four nuisance variables. Uncertainty uses 5,000 resamples of
the connected 5-degree blocks recomputed from that filter's association pairs.
All model-to-baseline changes use paired bootstrap replicates.

## Results

| Filter | Pairs / blocks | Redshift + richness | + outer residual | + coadd proxy | + both proxies |
|---|---:|---:|---:|---:|---:|
| `i2` | 53 / 28 | 0.239 [0.022, 0.423] | 0.252 [0.036, 0.456] | 0.211 [−0.006, 0.403] | 0.223 [0.009, 0.440] |
| `i` | 21 / 13 | 0.112 [−0.229, 0.514] | 0.107 [−0.273, 0.535] | 0.108 [−0.315, 0.671] | 0.098 [−0.374, 0.727] |

Entries are partial correlations with percentile 95% spatial-block bootstrap
intervals. For the combined model, the paired change from the redshift/richness
baseline is −0.016 in `i2` (95% interval −0.091 to +0.083) and −0.014 in `i`
(−0.300 to +0.272). The added terms therefore do not materially shift the
point estimate in this sample; the `i2` interval remains broad, while the `i`
sample is weakly constrained. These intervals are descriptive resampling
uncertainties for selected association pairs, not confirmatory detection tests.

## Interpretation and limits

The `i2` association direction persists under both diagnostic covariates, but
this does not establish a physical ICL–satellite relation. The `i2` sample is
small and selected; the `i` interval includes both substantial negative and
positive correlations. The coadd proxy is a finite-lag local residual statistic,
not the covariance matrix of the radial ICL estimator. The outer annulus may
contain cluster/ICL light, PSF wings, cirrus, and unmasked sources; it is not a
pure sky measurement and was not subtracted. The model does not calibrate
upstream HSC sky transfer, target-specific extended PSF wings, satellite
completeness, dust/passband effects, mask selection, or the full image-pixel
covariance.

The result remains exploratory with `science_measurement=false`. It is a joint
robustness check on the existing target/control stack, not a calibrated ICL
measurement.

## Reproduction and provenance

- Astrafek analysis: [`analyze_icl_satellite_joint_sky_covariance.py`](analyze_icl_satellite_joint_sky_covariance.py)
- Independent point-estimate recomputation: the saved 74-row table was refit by direct least squares; all four model estimates in both filters matched to `1e-12`.
- Output directory: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_joint_sky_covariance_2026-11-03_01/`
- Summary SHA-256: `947d9ecd8bb728e5a3bb1ec193f88687a893752403c9c7a4d409a3e6752888e2`
- Joined 74-row table SHA-256: `9ece7d3eb5ab88d5d8c7dbf298d33f90697a9414f33fbe8fe7847177bd423560`
- The run records the hashes of the association replay, outer-residual diagnostics, and independently audited coadd-covariance frame table in `summary.json`.
