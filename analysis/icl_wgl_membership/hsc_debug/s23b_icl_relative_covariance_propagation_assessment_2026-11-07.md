# S23B relative ICL contrast covariance propagation

Date: 2026-11-07
Status: completed conditional diagnostic; `science_measurement=false`

## Assessment

The measured short-lag covariance was propagated through the saved target-minus-control radial-contrast estimator for the 1,022 exploratory S23B matched pairs. The calculation replays the actual bilinear metric-grid warp, saved K=2.5 masks, pair-common pixels, and equal-pixel annular means. It checks that the replayed profile means agree with the saved stack and that the number of propagated pairs matches the existing 5-degree spatial-block support table.

The finite-lag image-noise standard deviation of the equal-pair mean is 8.9–9.8% larger than the independent-pixel VARIANCE-plane estimate across the six filter/annulus combinations. This component is smaller than the existing spatial-block interval half-widths, but it is not a total uncertainty estimate. It should not be added to profile error bars or used as a correction: the covariance kernel was estimated outside a 60-arcsec center exclusion and assumed spatially stationary when applied to cluster annuli. The result remains exploratory, and neither `science_measurement` nor the science-readiness gates are promoted.

## Method and validation

- The input stack contains 1,022 target/control pairs: 241 HSC `i` and 781 HSC `i2`.
- The three signed contrasts are the mean in 50–300, 100–300, or 200–300 proper kpc minus the mean in 300–800 proper kpc. Each target/control pair uses the same output-grid pixels valid on both sides.
- For each estimator, the output-grid coefficients are scattered back to native coadd pixels with the transpose of the production order-1 bilinear warp. The native coefficients are contracted with the delivered HSC VARIANCE plane, first diagonally, then with the measured zero-lag residual variance, and finally with each frame's normalized 9×9 lag kernel (lags through ±4 pixels).
- Profile means reproduced the saved stack within the runner's 2×10⁻⁵ nJy arcsec⁻² check. Propagated pair counts equal the saved 5-degree, phase-zero block-support counts for all six combinations.
- Support is 232, 232, and 231 pairs for HSC `i`; 756, 756, and 752 for HSC `i2`. The remaining contrast rows lack enough common pixels in the inner annulus and are marked `insufficient_common_support`.
- One `i2` target frame lacks a measured outer-field kernel. A per-filter/role median kernel was substituted in its three contrast rows; all three rows have insufficient common support and do not enter the reported summaries.
- The result contains 3,066 pair/contrast rows. Every row retains `science_measurement=false`.

## Results

All contrasts and noise values are in nJy arcsec⁻². The block intervals are the existing 5-degree, phase-zero spatial-block intervals; they are shown for scale and were not recalculated by this diagnostic.

| Filter | Inner minus outer annulus (kpc) | Pairs | Mean contrast | Existing block 95% interval | Diagonal σ | Zero-lag σ | Finite-lag σ | Finite-lag / diagonal σ |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| HSC `i` | 50–300 − 300–800 | 232 | 6.561 | [3.844, 8.667] | 0.423 | 0.452 | 0.463 | 1.094 |
| HSC `i` | 100–300 − 300–800 | 232 | 4.528 | [1.901, 6.537] | 0.430 | 0.460 | 0.470 | 1.094 |
| HSC `i` | 200–300 − 300–800 | 231 | 3.253 | [0.498, 5.441] | 0.562 | 0.602 | 0.616 | 1.095 |
| HSC `i2` | 50–300 − 300–800 | 756 | 7.611 | [5.722, 8.877] | 0.259 | 0.277 | 0.282 | 1.089 |
| HSC `i2` | 100–300 − 300–800 | 756 | 4.686 | [3.239, 5.707] | 0.353 | 0.377 | 0.388 | 1.098 |
| HSC `i2` | 200–300 − 300–800 | 752 | 1.908 | [0.401, 3.035] | 0.414 | 0.442 | 0.454 | 1.098 |

The finite-lag conditional noise is about 18–34% of the corresponding block-interval half-width. This is a scale comparison only: the two quantities describe different uncertainty components and must not be combined as if independent or complete. The cluster-to-cluster pair SEM is also larger than this conditional pixel-noise term in all six combinations.

## Scientific limits and use

- The outer 300–800 kpc annulus is a reference region, not demonstrated blank sky. Positive inner-minus-outer contrasts do not establish an absolute ICL surface-brightness profile.
- The measured 9×9 kernel comes from off-center outer-field patches, outside a 60-arcsec center exclusion. Applying it to the cluster annuli assumes spatial stationarity, which has not been tested.
- Only short-range covariance through four native coadd pixels is included. Long-range sky covariance, uncertainty in sky subtraction, source Poisson noise, PSF wings, extinction, centering, mask selection, and covariance-kernel uncertainty are absent.
- Target and control image-noise terms are treated as independent across clusters. Overlapping coadd footprints and large-scale covariance are not modeled.
- Keep the existing spatial-block uncertainty as the population-level reference. Do not use this conditional pixel term as a final error bar or as a correction to the stack.
- These radial contrasts do not by themselves establish the ICL–satellite association. Keep `science_measurement=false` and `images_ready_for_science=false`.

## Reproducibility record

- Runner: `analysis/icl_wgl_membership/hsc_debug/propagate_s23b_icl_relative_contrast_covariance.py`
- Runner SHA-256: `08d5d9997a9f8ffac922eacbb5f6cd4148ca52c82f4e8e68f2d3c5425bc8622a`
- Output directory: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_icl_relative_covariance_propagation_2026-11-07_02/`
- `summary.json` SHA-256: `9c1ef1740490f7cbc4811a9bd2573d69405f5482d7ce754a09642108507fc94e`
- `run_manifest.json` SHA-256: `aceb157f4a17e37adae1252af7f0ef2a8f6ce00251a89a91a2a26b04ef8fd747`
- `pair_measurements.tsv` SHA-256: `cc626c8c593c199fb296369da5e7fe695c40a7ee9d20c39642302b4a2408cd46`
- `filter_summary.tsv` SHA-256: `a69031bea620175ccbe08b63d0f3f4b5d7432eab5ee80bf7f52fbf94ab738ace`
- `relative_contrast_noise_budget.png` SHA-256: `2fcd9553cfdce60701f95c938ba0cae63db954652b71b88559986c60d9b3f6f3`
- Pinned stack profile SHA-256: `285906b03360e9d137727e948c15419767c356d3343e61f409562567ac4abc43`
- Pinned relative-contrast table SHA-256: `c5de5f32b51c6899c8d17fe2b3e4b0a3966d810afd8db5ee27ba649fc0642800`
- The run read existing local files only; it made no network request or data transfer.
