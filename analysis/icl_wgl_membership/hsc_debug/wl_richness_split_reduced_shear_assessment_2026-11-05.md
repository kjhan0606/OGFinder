# Redshift-matched richness contrast: reduced-shear sensitivity (2026-11-05)

## Aim and fixed inputs

This Astrafek model sensitivity asks whether the weak-shear NFW richness contrast changes when the forward model includes reduced shear integrated over the calibrated HSC Y3 source `n(z)`. Both fits retain the same 204 redshift-matched lens pairs in each richness group, the exact target shape-pair rows and estimator weights, the saved random-subtracted high-minus-low `DeltaSigma` profile, and its shared 50-block jackknife covariance. The fit range is 0.2–1.5 proper Mpc; concentration is fixed at 3, 4, or 5. A second fit excludes source tomographic bin 1 from both target groups.

The runner uses source-redshift efficiency moments through power 17 (geometric-series order 16), restricted to the nonzero support of each tomographic `n(z)`. It preserves the same single effective NFW mass per richness group used in the preceding weak-shear screen.

## Results

| Source bins | Fixed concentration | Weak-shear high/low mass ratio (conditional 68% profile interval) | Reduced-shear ratio (conditional 68% profile interval) | Reduced-shear fit `chi2/dof` | Maximum supported `kappa beta` | Conservative omitted-series bound |
|---|---:|---:|---:|---:|---:|---:|
| 1–4 | 3 | 1.894 [1.128, 8.629] | 1.168 [1.095, 2.399] | 2.40/2 | 0.691 | 0.60% |
| 1–4 | 4 | 1.396 [1.084, 2.856] | 1.200 [1.102, 2.033] | 2.91/2 | 0.555 | 0.010% |
| 1–4 | 5 | 1.209 [1.066, 1.833] | 1.214 [1.107, 1.887] | 3.58/2 | 0.485 | 0.0009% |
| 2–4 | 3 | 2.847 [1.307, 28.030] | 1.150 [1.095, 1.453] | 2.03/2 | 0.782 | 6.94% |
| 2–4 | 4 | 1.727 [1.160, 9.859] | 1.172 [1.102, 1.810] | 2.44/2 | 0.659 | 0.242% |
| 2–4 | 5 | 1.391 [1.099, 3.425] | 1.198 [1.109, 1.888] | 2.97/2 | 0.554 | 0.0097% |

The order-8 to order-16 mass-ratio changes are at most 0.3% across the six source-bin/concentration cases. The conservative *per-pair* omitted-series bound is larger in the most extreme `c=3`, source-bin-1-excluded case (6.94%); the fitted mass ratio remains numerically stable, but this local bound is retained as a model limitation.

## Interpretation

- Reduced shear pulls the `c=3` and `c=4` ratios toward approximately 1.2; the `c=5` ratio changes little. The previous weak-shear ratio was therefore substantially concentration-dependent, while the reduced-shear point ratios are more consistent across these three fixed concentrations.
- Excluding source bin 1 changes the reduced-shear point ratios by less than 2% at each concentration. This is a model sensitivity, not a source-photo-z calibration uncertainty propagation.
- The `Delta chi2 = 1` intervals are conditional profile intervals using the saved covariance. They are not detection significance: mass-distribution scatter, miscentering, magnification, photo-z calibration uncertainty, matched-random depth/mask/PSF/sky mismatch, and the validity of the approximate 50-block covariance remain unresolved.
- The observed sample is selected and the model assigns one effective mass per group. These values do not calibrate richness or measure a cluster-by-cluster ICL–satellite–mass relation.

Keep `science_measurement=false`. No HSC query, image read, or data download occurred.

## Reproduction and provenance

Runner: [`fit_wl_richness_split_reduced_shear.py`](fit_wl_richness_split_reduced_shear.py), SHA-256 `ec851c2d869a9680a22ab32379156df651306cf0ca6dae6c4425b4ee26c9f0a8`.

- Source bins 1–4 summary: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_richness_split_reduced_shear_2026-11-05_09/summary.json`, SHA-256 `1b1fbbafc7aaf90e57b1b651dc29979b955bb593cf2e8e5e99dbd89dd4c91be4`.
- Source bins 2–4 summary: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_richness_split_reduced_shear_2026-11-05_10/summary.json`, SHA-256 `c3885d73ec8dced7bb9a87fd96510a00b075de2dfdcd6d54156c11c601af4044`.

Each JSON binds the runner, the matched-richness split and source-bin-exclusion summaries, the source-row export integrity audit, the S19A lens-stack/n(z) provenance, and the exact row-export plan. Both outputs retain `science_measurement=false` and record no new data request.
