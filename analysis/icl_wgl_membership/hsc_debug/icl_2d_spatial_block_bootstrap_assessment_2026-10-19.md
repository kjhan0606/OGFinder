# Spatial-block uncertainty for the HSC ICL 2-D stack — 2026-10-19

## Scope and method

To test whether the radial target-minus-matched-control profile depends on treating every cluster/control pair as independent, the three verified 2-D stacks (Kron-mask factors 2.5, 3.5, and 5.0) were resampled in connected sky blocks. Each pair was linked to both its cluster center and matched-control center. Pairs sharing any rectangular RA/Dec sky cell were joined into the same connected block, including links where one pair's target lies in another pair's control cell. The 0–50 kpc bin was left out of the profile covariance because central BCG light and sparse support dominate it.

The bootstrap sampled connected blocks with replacement and retained equal weighting per supported pair within each replicate. A pair bootstrap is reported beside it. Both are descriptive uncertainty sensitivities; neither supplies a survey-calibrated detection significance. Cell widths of 2.5°, 5°, 10°, and 15° were tested to expose dependence on the chosen sky scale.

## Results

The intervals below show the envelope of 95% spatial-block bootstrap intervals across the four cell widths. They are ranges over the lower and upper endpoints, not a combined confidence interval. Mean profiles are the same values reported in the [mask-radius assessment](icl_2d_mask_sensitivity_assessment_2026-10-19.md).

| Filter | Kron mask | Blocks across widths | 50–100 kpc interval envelope | 100–200 kpc interval envelope | 200–300 kpc interval envelope |
|---|---:|---:|---:|---:|---:|
| i | 2.5 | 9–33 | [37.99, 45.10] | [10.62, 16.28] | [−0.76, 5.84] |
| i | 3.5 | 9–33 | [30.02, 40.59] | [7.66, 12.82] | [−1.26, 4.33] |
| i | 5.0 | 9–33 | [16.68, 38.97] | [4.48, 9.75] | [−1.55, 5.96] |
| i2 | 2.5 | 16–93 | [39.79, 49.21] | [8.65, 14.04] | [1.66, 4.58] |
| i2 | 3.5 | 16–93 | [32.09, 43.76] | [7.04, 12.38] | [0.87, 3.88] |
| i2 | 5.0 | 16–93 | [29.50, 44.41] | [6.87, 12.51] | [0.73, 3.82] |

Across all tested mask choices and cell widths, the block-bootstrap intervals for the 50–100 and 100–200 kpc bins stay above zero in both filters. The amplitude declines as the masks grow, especially in the 50–100 kpc bin. This gives a useful indication that the intermediate-radius excess is not solely an independent-pair error-bar artifact, while also showing that its size depends on how source wings are masked.

For `i`, the 200–300 kpc interval crosses zero at every cell width and mask factor. For `i2`, it stays positive in these resampling runs, but the median common-area coverage falls to about 4% at K=5.0. At 300–500 kpc, `i2` remains compatible with zero. The `i` outer offset changes with mask size and has only a few percent common-area coverage under K=5.0, so it is not a stable diffuse-light constraint.

Spatial resampling did not uniformly enlarge the pair-bootstrap errors. The two intervals are similar for much of `i2`, while some `i` bins have narrower spatial-block scatter than pair-bootstrap scatter. This reflects the finite, unequal connected blocks and should not be taken as evidence that either estimate is conservative. Depending on filter and cell width, the profile is supported by 9–93 band-specific blocks; the largest connected block contains 62–151 pairs across the four widths.

## Scientific interpretation

The 50–200 kpc target-minus-control excess is a repeatable exploratory feature under these mask and resampling choices. It is not yet an ICL detection. The block boundaries are approximate rectangular RA/Dec cells; cells may combine disconnected HSC fields, and large blocks reduce the effective number of independent units. The intervals omit residual sky transfer, PSF wings, extinction/calibration, source and sky photon noise, pixel-level coadd covariance, and selection effects. The 2–4% common support at K=5.0 also makes the outer profile vulnerable to small-area residuals.

This analysis does not measure the relation between ICL and satellite outer light. The existing satellite association remains exploratory and inconsistent across filters/selections. The next ICL science gate is to establish the sky/PSF transfer and support completeness before combining cluster-centered ICL measurements with satellite-centered profiles.

## Reproducibility

The Astrafek utility is [`bootstrap_icl_2d_spatial_blocks.py`](bootstrap_icl_2d_spatial_blocks.py), with focused tests in [`test_bootstrap_icl_2d_spatial_blocks.py`](test_bootstrap_icl_2d_spatial_blocks.py). Four hash-linked run directories outside Git retain the 30-row tables, covariance matrices, bootstrap settings, figures, input hashes, and code hash:

- 2.5°: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_2d_spatial_block_bootstrap_w2p5_2026-10-19_03/` (110 connected blocks; table SHA-256 `eb6d636eb0b6daabe136bcc11fe73566dabcbef395121c525f9cd7e7c435c337`)
- 5°: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_2d_spatial_block_bootstrap_w5_2026-10-19_03/` (38 blocks; table SHA-256 `2733f76e4479a3ef9366c88efb5fc6642655b143247e6dda81eefa1e941f517f`)
- 10°: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_2d_spatial_block_bootstrap_w10_2026-10-19_03/` (20 blocks; table SHA-256 `42a87ec4a156c09edd27e3e907ee0c798ca0b09f8dc47fd71167dd5d5242b0e3`)
- 15°: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_2d_spatial_block_bootstrap_w15_2026-10-19_03/` (17 blocks; table SHA-256 `2e8c48cf4ba5f4cf2c979ec797285baa6312004703bf1a6e5cd32ecfac3be517`)

An independent consistency pass checked all 120 table rows, their output hashes, the per-bin pair counts, means, and pair SEMs against the three source evaluation tables. Two focused tests pass. All outputs retain `science_measurement=false`.
