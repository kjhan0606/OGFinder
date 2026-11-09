# Random-subtracted weak-lensing reduced-shear sensitivity (2026-11-05)

Status: exploratory model sensitivity; `science_measurement=false`.

## Scope and calculation

The analysis refits the saved random-subtracted profile of the exact 435-lens
common sample. It reuses source bins 2–4, the HSC Y3 calibrated tomographic
`n(z)`, exact pair radii and lens redshifts, the stack's selection and shear
weights, the 0.2–1.5 proper-Mpc fit range, and the saved paired endpoint
jackknife covariance at 5° and 2.5°. No HSC request or image read was made.

For each source bin and lens redshift, the runner integrates
`mu_p = integral n(z_s) [D_ls/D_s]^p dz_s` through `p=9`. It evaluates the
source-redshift-averaged reduced shear using the geometric expansion through
order 8 in convergence, then applies the same estimator weights as the weak
shear fit. With `beta <= 1`, the maximum `kappa_infinity` over the fitted rows
is 0.35–0.37. The resulting conservative omitted-series tail bound is
0.012–0.022% across the six fits.

## Fits

Masses are effective single-halo `M200c` scales in solar masses. Intervals are
conditional profile intervals from the saved covariance; they do not include
the systematics listed below.

| Fixed concentration | 5° `M200c` | 5° `chi2/dof` | 2.5° `M200c` | 2.5° `chi2/dof` |
|---:|---:|---:|---:|---:|
| 3 | `1.750e14` | 1.39/3 | `1.730e14` | 2.11/3 |
| 4 | `1.604e14` | 12.58/3 | `1.525e14` | 8.69/3 |
| 5 | `1.450e14` | 32.54/3 | `1.374e14` | 19.03/3 |

Relative to the weak-shear model, the fitted mass shifts down by 1.7–5.5%,
depending on concentration and covariance block size. The `c=5` fits remain
poor; only `c=3` has acceptable conditional profile residuals for both block
widths. The fixed-concentration grid still changes the fitted mass by about
20%, larger than the reduced-shear shift. This ranking is conditional on the
selected one-halo model and four radial bins; it is not a concentration
measurement.

## Interpretation limits

This correction addresses only reduced shear within the single spherical NFW
model. The mass distribution of the richness-selected sample, concentration
distribution, magnification, miscentering, two-halo structure, photo-z
calibration uncertainty, random-aperture depth/mask/PSF/sky matching, and
covariance uncertainty remain unresolved. The first-order diagnostic stored
in the earlier `_01` output had a conservative omitted-term bound as large as
19%; it is superseded by the order-8 result here and should not be used for
interpretation. The `_02` output is an identical order-8 run before recording
the runner hash; `_03` is the provenance-complete result. None of these fits
is a detection or science-grade mass calibration. Keep
`science_measurement=false`.

## Reproducibility

Astrafek runner: [`fit_wl_random_stack_reduced_shear_sensitivity.py`](fit_wl_random_stack_reduced_shear_sensitivity.py),
SHA-256 `92394e9ba11bc49582c7aab23a79aeb01d8d7798d1656a6bf6a42d80a8a3c556`.

Machine summary:
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_random_stack_reduced_shear_sensitivity_2026-11-05_03/summary.json`,
SHA-256 `310a03e29069b7799239f21a9ae3c0d0ef0abda5e2d0da93c9e8bdc2b50bc286`.
It binds the weak-shear NFW screen, random-subtracted stack, row-integrity
audit, tomography profile, full stack, HSC Y3 `n(z)` SACC, and shape-row plan.
