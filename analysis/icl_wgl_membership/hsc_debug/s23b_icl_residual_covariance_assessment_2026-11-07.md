# S23B coadd residual-covariance diagnostic (2026-11-07)

## Scope and decision

This CPU-only diagnostic measured short-range residual covariance in the exact
1,022 target/matched-control pairs used by the exploratory S23B flux stack:
241 HSC `i` pairs and 781 HSC `i2` pairs. It used existing coadds, their
delivered VARIANCE planes, and saved K=2.5 source masks. No new data or network
query was used. All 2,044 frames were read; 2,043 had enough outer-field
patches for a covariance estimate. One HSC `i2` target frame had fewer than
eight eligible patches.

The diagnostic does not correct the stack or promote any measurement. Keep
`science_measurement=false` and `images_ready_for_science=false`.

## Method

For each frame, the analysis selected up to 16 non-overlapping 64×64-pixel
patches whose centers lie at least 60 arcsec from the target or random-control
center. A patch needed at least 80% valid pixels after source masks, HSC
bad-pixel masks, and finite positive VARIANCE checks. The saved K=2.5 mask was
required to include the base exclusions.

An iterative, 4-sigma-clipped affine plane was removed from each patch. The
residual was divided by the square root of the delivered VARIANCE plane. The
mean product of normalized residuals was measured at every pixel displacement
from −4 to +4 pixels in both axes. Patch estimates were summarized by their
median for each frame. The reported approximate variance factor is the sum of
that finite 9×9 lag-covariance grid, which corresponds only approximately to a
square 9×9-pixel top-hat aperture.

## Results

| Filter and frames | Eligible frames (target/control) | Median zero-lag residual variance / VARIANCE (target/control) | Median approximate 9×9 variance factor (target/control) | Median sigma factor (target/control) | Paired target minus control factor |
|---|---:|---:|---:|---:|---:|
| HSC `i` | 241 / 241 | 1.143 / 1.143 | 1.586 / 1.580 | 1.259 / 1.257 | +0.0226; frame-pair 16th–84th spread −0.170 to +0.203 |
| HSC `i2` | 780 / 781 (one target frame insufficient) | 1.139 / 1.140 | 1.572 / 1.555 | 1.254 / 1.247 | +0.0291; frame-pair 16th–84th spread −0.160 to +0.205 |

The frame-level 16th–84th percentile ranges of the approximate variance factor
are 1.455–1.747 (`i` target), 1.433–1.713 (`i` control), 1.418–1.746 (`i2`
target), and 1.406–1.704 (`i2` control). These describe frame-to-frame spread;
they are not confidence intervals. The paired-difference ranges above are also
descriptive spreads, not confidence intervals.

The zero-lag ratios are about 1.14 in both target and control images. The
clipped, plane-subtracted outer-patch residual variance is therefore roughly
14% above the delivered VARIANCE prediction in this diagnostic. The
approximate 9×9 factors suggest that short-lag correlations could increase a
simple aperture variance by roughly 56–59% (about 25% in standard deviation)
relative to an independent-pixel calculation. Target and control medians are
close; paired differences are small compared with the frame-pair spread. This
points to a shared coadd-noise issue and does not establish a target-specific
excess.

## Limits and next analysis

These values must not be applied as a blanket error-bar multiplier. The
outer-field patches can still contain unresolved sources or diffuse cluster
light, so the residual covariance is not a pure instrumental-noise estimate.
The finite 9×9 kernel is not propagated through Lanczos reprojection or the
exact signed annular-profile weights. Long-range covariance, source Poisson
noise, sky-transfer uncertainty, PSF-wing uncertainty, and member-selection
uncertainty are not included. The frame-level percentiles do not quantify
uncertainty in the population median.

The next covariance step is to propagate image-level lag kernels through the
actual reprojection and annular estimator, then compare against connected
sky-block resampling and injection/recovery. The present result indicates
that independent-pixel errors are likely optimistic, but the stack is not yet
ready for science measurements.

## Reproducibility

- Runner: `analysis/icl_wgl_membership/hsc_debug/estimate_s23b_icl_coadd_noise_correlation.py`
- Helper: `analysis/icl_wgl_membership/hsc_debug/estimate_hsc_coadd_noise_correlation.py`
- Data directory: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_icl_residual_covariance_2026-11-07_01/`
- Output summary SHA-256: `8034a140deb4e52cb9c47fff7a2475dcaec19f0410305c5c0520995ab9eedbe1`
- Frame table SHA-256: `f6411a83c0d28b0a60eba2beaafc816e0d34671611feca518a534e664fe39e5d`
- Lag table SHA-256: `e5b13a6edff5d5d3a7f2475e29bb5235165ea92a78d6336f7552480ffb08609a`
- Plot SHA-256: `5ed405db3705e287302e90b5869af90540c64b2ea549515c1692bc1b9ff27b17`
- Runner SHA-256: `ba9f51e3c84d65ad6b38dc41d0fda97f6b07c072da01fa98f07a6feba3763513`
- Analysis-helper SHA-256: `bddebcf1839e5ac58fe9281321cc29d29221f88670eccacaeb7f12fc9cc33b9d`

The saved input hashes, all three output hashes, frame count, lag-row count,
and final state summary hash were independently checked after completion.
