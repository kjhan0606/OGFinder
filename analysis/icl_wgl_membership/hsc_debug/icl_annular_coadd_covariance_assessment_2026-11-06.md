# Annular coadd covariance comparison (2026-11-06)

## Result

Compared the existing HSC coadds in the 200–300 proper-kpc ICL science annulus
with an outer 350–500 proper-kpc annulus for the fixed, same-filter 94-pair
target/control sample. The outer region is not assumed to be blank sky; it may
contain cluster-envelope light. This is a covariance diagnostic, not an ICL
measurement or a correction.

Each frame was measured after the existing local affine-plane subtraction,
iterative clipping, and variance normalization. Non-overlapping 64×64-pixel
patches had to be fully contained in the annulus and at least 80% unmasked;
the saved target Kron-2.5 masks and matched-control SEP/Kron-2.5 masks were
used. The estimator records a normalized 17×17 lag kernel (±8 pixels). Paired
target-minus-control contrasts and radial difference-in-differences use 5°
connected-sky-block bootstrap intervals with 5,000 replicates.

## Paired results

Values below are medians of target-minus-control contrasts; intervals are
95% connected-block bootstrap intervals. The radial row is
`(target-control)_200–300 − (target-control)_350–500`.

| Filter | Metric | ICL 200–300 kpc (pairs / blocks) | Outer 350–500 kpc (pairs / blocks) | Radial difference-in-differences (pairs / blocks) |
|---|---|---:|---:|---:|
| i | zero-lag | +0.00073 [−0.01410, +0.01293] (8 / 7) | −0.00153 [−0.00851, +0.00284] (24 / 13) | +0.00924 [−0.01349, +0.03235] (8 / 7) |
| i | sum of lags ±8 | +0.234 [−0.283, +1.251] (8 / 7) | +0.048 [−0.073, +0.202] (24 / 13) | +0.179 [−0.800, +1.068] (8 / 7) |
| i | valid-pixel coverage | −0.087 [−0.208, +0.116] (8 / 7) | +0.004 [−0.015, +0.105] (24 / 13) | −0.044 [−0.154, +0.076] (8 / 7) |
| i2 | zero-lag | −0.00386 [−0.01663, +0.01095] (16 / 9) | −0.00244 [−0.00673, +0.00192] (65 / 31) | +0.00755 [−0.00406, +0.01736] (15 / 8) |
| i2 | sum of lags ±8 | −0.031 [−0.354, +0.220] (16 / 9) | +0.021 [−0.066, +0.182] (65 / 31) | −0.256 [−0.524, +0.471] (15 / 8) |
| i2 | valid-pixel coverage | +0.019 [−0.100, +0.067] (16 / 9) | +0.054 [−0.0004, +0.095] (65 / 31) | −0.055 [−0.136, +0.055] (15 / 8) |

Every reported interval includes zero. The usable inner-annulus sample is
small because both target and control must have at least four fully contained
valid patches: only 8 `i` and 16 `i2` pairs survive. Thus the comparison is
underpowered to establish equivalence or transfer the outer-annulus kernel to
the ICL annulus. The larger outer-ring counts do not remove that limitation.

## Reproducibility and limits

- Runner: `analyze_icl_annular_coadd_covariance.py`, SHA-256
  `16ad75b7f9ac3f2e0cd2e1759a84037f4c9519307e982a03e32174526b80f7d3`.
- Results: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_annular_coadd_covariance_2026-11-06_01/`.
- `summary.json` SHA-256: `c03b1c60fcc3415aee615c8685700b42633fa769df0ac4a8173013bafe2f93b0`.
- Output tables: `per_center_annular_covariance.tsv`,
  `annular_lag_covariance.tsv`, and `paired_annular_covariance_comparison.tsv`;
  their hashes are pinned in `summary.json`.
- The summary records input manifest, mask, supported-pair, and helper-script
  hashes. It reports 94 input pairs, 264 successful center-annulus estimates,
  no failed inputs, and no HSC network use.

This comparison does not model source-Poisson noise, exact HSC resampling
covariance, PSF wings, or radial sky transfer. The outer annulus can contain
real diffuse cluster light. No covariance correction was applied and
`science_measurement=false` remains.

The separately approved random-point variance export remains a distinct data
step. Its exact plan is already approved, but the prior secure-credential
submission session ended before any CAS job ID was recorded; its checkpoint
still has zero jobs and no query/data transfer. This analysis made no query or
download.
