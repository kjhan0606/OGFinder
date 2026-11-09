# Weak-lensing source-density boost and DeltaSigma joint covariance (2026-11-05)

## Scope

This Astrafek diagnostic uses the existing approved HSC S19A shape-row exports
for 435 common-support cluster lenses and 2,175 matched random centers. Five
random centers accompany each lens. The source rows pass the same magnitude,
resolution, shape, B-mode, tomographic-efficiency, and radial cuts as the saved
weak-lensing stack. The analysis makes no HSC request, reads no images, and
does not modify the saved `DeltaSigma` profile.

The 2,175-center denominator retains all planned centers, including 129 random
centers with no exported shape rows. Paired delete-one equal-area sky cells
remove a lens/random center and source endpoint together. The analysis uses 52
endpoint cells at 5 degrees and 159 at 2.5 degrees. The five additional
2.5-degree cells and one additional 5-degree cell contain planned centers but
no lens-source pairs. The additional cells therefore change the density
denominator while leaving the shear estimator unchanged.

## Corrected estimator definition

The estimator-weight boost uses the actual inverse critical surface density,

\[
\Sigma_{\rm crit}^{-1}=\frac{4\pi G}{c^2}D_l\,\beta,
\qquad \beta=\left\langle\frac{D_{ls}}{D_s}\right\rangle_{n(z)}.
\]

The estimator weight per selected source is the HSC S19A response-adjusted
weight

\[
w_{\rm shear}\,\Sigma_{\rm crit}^{-2}
=w\,[2R(1+m)(1+m_{\rm sel})]^2\,\Sigma_{\rm crit}^{-2},
\]

with radial-bin responsivity $R$, per-object multiplicative correction $m$,
and HSC edge-selection correction $m_{\rm sel}$. The efficiency $\beta$ uses
the saved HSC Year 3 tomographic $n(z)$ and fixed cosmology
$H_0=70\ {\rm km\ s^{-1}\ Mpc^{-1}}$ and $\Omega_m=0.3$.

The previous spatial-boost outputs `_01` and `_02` used $\beta^2$ without the
lens-distance factor $D_l^2$ and therefore did not contain the exact
S19A estimator-weight ratio. The corrected `_03` output below includes
$\Sigma_{\rm crit}^{-2}$. The count ratios are unchanged. Do not use the
estimator-weight values from the earlier outputs.

## Results

The count ratio $B_N$ compares selected source pairs per planned lens and
random center. The response-adjusted weight ratio $B_W$ compares the exact
S19A estimator-weight sums per planned center. Values below show the full
sample point estimate and the 2.5-degree delete-one jackknife standard error.

| Proper radius (Mpc) | Count ratio $B_N$ | Exact S19A weight ratio $B_W$ |
|---:|---:|---:|
| 0.1–0.2 | $1.765\pm0.064$ | $1.139\pm0.048$ |
| 0.2–0.4 | $1.520\pm0.050$ | $1.142\pm0.039$ |
| 0.4–0.7 | $1.358\pm0.041$ | $1.163\pm0.034$ |
| 0.7–1.0 | $1.248\pm0.037$ | $1.132\pm0.034$ |
| 1.0–1.5 | $1.190\pm0.033$ | $1.125\pm0.033$ |
| 1.5–2.0 | $1.165\pm0.032$ | $1.134\pm0.032$ |
| 2.0–3.0 | $1.148\pm0.030$ | $1.126\pm0.031$ |

The 5-degree point ratios agree with the 2.5-degree ratios to stored
precision. The count excess declines from 76.5% in the innermost bin to 14.8%
at 2–3 Mpc. The exact estimator-weight excess remains between 12.5% and 16.3%.
The count ratio therefore overstates the weight change entering the current
shear estimator, especially at small radii.

The analysis stores a 42-by-42 jackknife covariance across four boost metrics
and fourteen random-subtracted tangential/cross `DeltaSigma` values. The
cross-correlation coefficients between boost parameters and `DeltaSigma`
parameters range from -0.346 to 0.325 at 2.5 degrees and from -0.393 to 0.465
at 5 degrees. The dimensionless correlation matrices are symmetric and
positive definite, but poorly conditioned. The minimum-to-maximum eigenvalue
ratios are about `2.94e-6` at 2.5 degrees and `6.89e-8` at 5 degrees. Any
future generalized least-squares fit needs a documented covariance
regularization and block-scale sensitivity. The complete covariance must be
used for any future profile-level propagation.

The sufficient-statistic reconstruction reproduces the saved full-sample
`DeltaSigma` and cross-shear profiles at both block widths. The reconstruction
also reproduces the saved `DeltaSigma` covariance when restricted to the
original 51 pair-endpoint blocks at 5 degrees and 154 at 2.5 degrees. The
current joint matrices then add the planned-center-only cells, retaining the
same jackknife deletion set used for the boost denominator.

## Interpretation and remaining gates

The corrected source population around cluster centers carries about 12–16%
more S19A estimator weight per center than the matched random population over
0.1–3 Mpc. The count excess is larger near the center. Source-density excess
alone does not identify cluster-member contamination. Magnification, depth,
mask, PSF, sky transfer, and source-selection differences can also change the
source density and estimator weight.

The spatial covariance quantifies variation across the chosen cells. The
covariance does not establish matched full-aperture depth, mask, PSF, or sky
response. The HSC Year 3 $n(z)$ uncertainty, dilution calibration, and reduced
shear remain unpropagated. Do not multiply the saved `DeltaSigma` profile by
the boost ratio. Keep `boost_correction_applied=false` and
`science_measurement=false` until the coverage and calibration gates pass.

## Reproduction and provenance

Run the Astrafek analysis from the repository root with:

```bash
python3 analysis/icl_wgl_membership/hsc_debug/analyze_wl_source_density_boost_blocks.py \
  --output /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_source_density_boost_spatial_2026-11-05_05
```

The corrected diagnostic is
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_source_density_boost_spatial_2026-11-05_04/`.
The summary records `status=complete_diagnostic`,
`hsc_network_used=false`, `images_read=false`,
`boost_correction_applied=false`, and `science_measurement=false`.

SHA-256 values:

- Analysis script: `76e29cc7c9c284c3f77beacdcb3b2d2624f86eb4bd1b72b62b628daf367cecdb`
- Corrected summary: `04fff9521b4b9f78eca88ea3c72620a3cc4eeadefdeb5fc3bec238cbfef66a70`
- Seven-bin table: `9a9ce71bde5ca994728008269b58c80797369d3de845de5b7f944a6d37abd834`
- Tomographic $n(z)$ input: `a28f9e2e088e92d96d2d87083a6eeeac4c0c84c6b48e033bd3d1dc1bf0d8958f`
