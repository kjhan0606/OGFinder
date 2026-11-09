# Random-subtracted, redshift-matched richness split (2026-11-05)

Status: exploratory population contrast; `science_measurement=false`.

## Sample and method

This Astrafek analysis tests whether the earlier cluster-only richness contrast
remains after matched-random subtraction and lens-redshift matching. It reads
only the previously approved local HSC S19A target and random shape-row
exports, the existing CAMIRA parent manifest, and the HSC Y3 tomographic
`n(z)` file. It made no HSC request and read no images.

Richness labels follow the earlier deterministic split of the full 800-lens
parent, so the richness boundary was not selected from these lensing results.
Among 435 lenses in the saved target/random common source-row support, 415
have at least one estimator-usable target pair and at least one
estimator-usable matched random center after the exact B-mode, source quality,
tomographic, and radial cuts. This comprises 211 lower-richness and 204
upper-richness lenses. Maximum-cardinality 1:1 redshift matching within
`|delta-z| <= 0.02` yields 204 pairs (204 target lenses per group); mean
`|delta-z|=0.00650`, median `0.00341`, maximum `0.01998`. Median richness is
17.41 versus 26.63; mean lens redshift is 0.4333 versus 0.4289.

Each matched lens retains its five planned random centers (2,040 total); 1,911
centers have exported shape rows. The target and random stacks are computed
separately for each richness group. A shared 5-degree equal-area jackknife
deletes both lens and source endpoints across all four stacks, producing 50
common spatial blocks. The contrast is
`(high target - high random) - (low target - low random)`.

## Results

The profile and errors are in `10^12 M_sun Mpc^-2`.

| Proper R (Mpc) | High−low random-subtracted `DeltaSigma` | Jackknife error | High−low cross / error |
|---:|---:|---:|---:|
| 0.1–0.2 | −9.63 | 27.19 | +0.24 |
| 0.2–0.4 | +60.16 | 18.09 | +1.03 |
| 0.4–0.7 | +21.23 | 10.15 | +1.89 |
| 0.7–1.0 | +32.60 | 8.61 | −1.22 |
| 1.0–1.5 | +20.75 | 5.90 | +0.31 |
| 1.5–2.0 | +13.74 | 4.11 | +0.17 |
| 2.0–3.0 | +5.21 | 2.85 | +0.07 |

The higher-richness group's random-subtracted tangential profile is larger in
all bins from 0.2 to 3 Mpc. The seven-bin conditional tangential Wald statistic
is `chi2=67.84` (rank 7; nominal `p=4.04e-12`). This p-value is not a detection
significance: the covariance uses approximate equal-area cells rather than
validated contiguous HSC regions, and does not include the systematics below.
The cross-profile contrast has no bin above 1.9 jackknife standard errors.

The first 214-pair run (`..._01`) used source-row presence as support and left
10 lenses per group with no tomographic contribution to the target profile.
That output is retained as an audit trail but is superseded. The result here
filters on the exact estimator-usable target and random support before
redshift matching; do not quote the earlier 214-pair contrast.

## Interpretation limits

This result shows that the richness contrast survives random subtraction and
redshift matching in this selected HSC-supported subset. It does not establish
a calibrated richness–mass relation, a causal effect, or an ICL–satellite
association. The sample is conditioned on target and random source support and
then reduced from 415 usable lenses to 204 redshift-matched pairs. The matched
random apertures are not proven equivalent in pixel depth, masks, PSF, sky, or
coadd weights. The equal-area jackknife blocks are approximate. Fiducial
photo-z `n(z)` uncertainty, shear calibration, reduced shear, magnification,
concentration and mass distributions, and covariance uncertainty are not
propagated. Keep `science_measurement=false`.

## Reproducibility

Astrafek runner:
[`analyze_wl_random_subtracted_richness_split.py`](analyze_wl_random_subtracted_richness_split.py),
SHA-256 `6ecd77d5c1ce81b0595a4e04b8aa1b09b118929921d1e4926500adc25004cadb`.

Machine summary:
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_random_subtracted_richness_split_2026-11-05_02/summary.json`,
SHA-256 `b38a78f6005dc17d8a12ddd6117debd63a005d0b21ff8dbda05961bd31eb65b7`.
The summary records the target/random row batch hashes, count-plan hash,
common-support ID hash, parent manifest, n(z) SACC, jackknife replicates,
and the science gate.
