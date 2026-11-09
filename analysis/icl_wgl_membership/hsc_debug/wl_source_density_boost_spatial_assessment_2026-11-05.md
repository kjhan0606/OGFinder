# Spatially paired weak-lensing source-density boost audit (2026-11-05)

> Superseded for estimator-weight values by
> [`wl_source_density_boost_joint_covariance_assessment_2026-11-05.md`](wl_source_density_boost_joint_covariance_assessment_2026-11-05.md).
> This preliminary analysis used dimensionless `beta` without the lens-distance
> factor in `SigmaCritInv`. Its count ratios remain valid. Do not cite its
> estimator-weight ratios or covariance.

## Question and sample

The weak-lensing stack subtracts matched random-center shear, but the stack has
not applied a source-density boost correction. This audit compares the selected
HSC S19A source density and estimator weight around 435 cluster lenses with
2,175 matched random centers. Five random centers accompany each lens. The
denominator retains all planned centers, including 129 random centers with no
exported shape rows.

The audit reads the already approved local shape-row exports, the saved 5-degree
and 2.5-degree stack summaries, and the HSC Year 3 tomographic $n(z)$ data. The
audit makes no HSC request, reads no images, and leaves the saved
random-subtracted $\Delta\Sigma$ profile unchanged.

## Estimators and uncertainty

For each proper-radius bin, the count ratio is

\[
B_N(R)=\frac{N_{\rm lens}(R)/435}{N_{\rm random}(R)/2175}.
\]

The estimator-weight ratio uses the sum of the exact S19A response-adjusted
\(w\beta^2\) lensing estimator weights per planned center. The lensing
efficiency $\beta$ comes from the saved tomographic $n(z)$ and the fixed
cosmology $H_0=70\ {\rm km\ s^{-1}\ Mpc^{-1}}$ and $\Omega_m=0.3$. Both ratios
use the same magnitude, resolution, shape, B-mode, positive-efficiency, and
radial selection as the saved stack.

Uncertainties use paired delete-one equal-area cells in RA and $\sin(\mathrm{Dec})$.
Each replicate removes a pair when either the center endpoint or source endpoint
falls in the deleted cell. The replicate also removes center positions from the
density denominator, including positions without exported source rows. The
analysis uses 52 endpoint cells at 5 degrees and 159 endpoint cells at 2.5
degrees. Reported $\sigma_{\rm JK}$ values are delete-one jackknife standard
errors. The reported ranges are the minimum and maximum leave-one-cell ratios,
not confidence intervals. The analysis stores the full 28-by-28 covariance
across four boost metrics and seven radial bins for each block width. The
covariance diagonal reproduces the reported per-bin errors.

## Results

| Proper radius (Mpc) | Count ratio $B_N$ | Exact estimator-weight ratio $B_W$ |
|---:|---:|---:|
| 0.1–0.2 | $1.765\pm0.064$ | $1.157\pm0.052$ |
| 0.2–0.4 | $1.520\pm0.050$ | $1.145\pm0.041$ |
| 0.4–0.7 | $1.358\pm0.041$ | $1.167\pm0.036$ |
| 0.7–1.0 | $1.248\pm0.037$ | $1.132\pm0.036$ |
| 1.0–1.5 | $1.190\pm0.033$ | $1.126\pm0.036$ |
| 1.5–2.0 | $1.165\pm0.032$ | $1.134\pm0.035$ |
| 2.0–3.0 | $1.148\pm0.030$ | $1.127\pm0.033$ |

Each value gives the full-sample point ratio and its delete-one jackknife
standard error. The 5-degree and 2.5-degree point ratios agree to stored
precision. The count ratio falls from a 76.5% excess in the innermost bin to a
14.8% excess at 2–3 Mpc. The exact estimator-weight ratio stays between 12.6%
and 16.7% above unity. Source counts therefore overstate the radial change in
the weight entering the current shear estimator.

The analysis also computed raw shape-weight and $w\beta^2$ ratios. The inner-bin
raw shape-weight ratio is 1.951 while the $w\beta^2$ ratio is 1.143. The
difference shows why the unweighted count ratio cannot be used as a direct
multiplicative correction to $\Delta\Sigma$.

Adjacent radial bins in the exact estimator-weight boost are strongly
correlated. Their correlation ranges from 0.788 to 0.961 at 2.5-degree blocks
and from 0.843 to 0.979 at 5-degree blocks. The covariance matrices are
symmetric and positive definite in both saved diagnostics. A diagonal-only
boost error would discard substantial cross-bin information.

## Interpretation and remaining gates

The selected source population around cluster centers carries a higher
per-center estimator weight than the matched random population across all seven
radial bins. The estimator-weight excess is much smaller and flatter than the
raw count excess. The result motivates a source-dilution and selection audit,
but the ratio does not identify cluster-member contamination by itself.
Magnification, depth, mask, PSF, sky transfer, and source-selection differences
can also change the source density and weight.

The delete-one errors quantify variation over the chosen spatial cells. The
errors do not establish matched full-aperture depth, mask, PSF, or sky response.
The HSC Year 3 $n(z)$ uncertainty, dilution calibration, and reduced shear
remain unpropagated. The boost covariance is available, but the cross-covariance
between boost and random-subtracted $\Delta\Sigma$ remains uncomputed. Do not
multiply the saved $\Delta\Sigma$ profile by either ratio.
Keep `boost_correction_applied=false` and `science_measurement=false` until the
coverage and calibration gates pass.

## Reproduction and provenance

Run the Astrafek analysis from the repository root with:

```bash
python3 analysis/icl_wgl_membership/hsc_debug/analyze_wl_source_density_boost_blocks.py \
  --output /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_source_density_boost_spatial_2026-11-05_03
```

The analysis refuses to overwrite an existing output directory. The saved
diagnostic used for this assessment is
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_source_density_boost_spatial_2026-11-05_02/`.
The saved JSON reports `status=complete_diagnostic`,
`hsc_network_used=false`, `images_read=false`,
`boost_correction_applied=false`, and `science_measurement=false`.

SHA-256 values:

- Analysis script: `99a315fc6792e01482bae02c68760847997504886f976a21f15387368812327b`
- Diagnostic summary: `ed915f00e335253ffe01f9710bca6e0b8da4d4a0f5ab2f2a4a81d342e2cde18f`
- Seven-bin profile: `3b4b4358fd2028c2b8b31e5d617a40e875a615b108ec2b7fab570f17742c6f2b`
- Tomographic $n(z)$ input: `a28f9e2e088e92d96d2d87083a6eeeac4c0c84c6b48e033bd3d1dc1bf0d8958f`
