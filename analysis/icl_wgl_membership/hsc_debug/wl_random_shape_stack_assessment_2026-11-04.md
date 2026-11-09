# HSC S19A random-subtracted weak-lensing stack (2026-11-04)

## Result status

The Astrafek CPU pipeline completed the approved 124-batch HSC random-center
shape-row export and formed a matched random-subtracted profile. The result is
**exploratory only**: `science_measurement=false` remains set in the machine
summary. This assessment does not promote it to a detection claim.

Machine result:
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_random_shape_stack_lens_reuse_2026-11-04_01/summary.json`

SHA-256: `83c5484454033b2652adcba5bd0699f28082fa6983eefa960613da75491bafc8`

Figure:
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_random_shape_stack_lens_reuse_2026-11-04_01/random_subtracted_profile.png`

Figure SHA-256: `3fe5cd60281a42935adb7a97093979288e5a273472b7d9d8e4dad4745a8478f3`.

## Export and common support

- The approved plan SHA-256 is `fc883444d2b23663d5f6034d773b76db948add29f1cda0c23384d17cfed2996a`.
- All 124 batches completed: 3,100 random centers associated with 620 matched
  cluster IDs, 9,010,324 center-source rows, and 2,962,522,918 CSV bytes.
- The original random-center manifest covered five potential centers per
  matched lens. The shape count preflight found 150 of the 620 cluster lenses
  with no target-side shape rows, and 166 with no source rows at any associated
  random center. These groups overlap in 131 lenses, leaving 435 lenses with
  common target and random support.
- The matched estimator therefore uses those same 435 lenses on both sides,
  with 2,175 associated random centers. Of these, 2,046 centers have positive
  source counts and 129 have zero; zero-count centers contribute no pairs.
  The paired random input contains 8,905,892 center-source rows.
- This support correction is essential: using all 620 target IDs with random
  rows would compare different effective lens samples.

## Radial profile

Errors below are the diagonal errors from the 51-block paired delete-one
jackknife. Marginal signal-to-error ratios are descriptive, not independent
bin detections.

| Proper radius (Mpc) | Random-subtracted ΔΣ (10¹³ M☉ Mpc⁻²) | Cross component (10¹³ M☉ Mpc⁻²) |
|---:|---:|---:|
| 0.1–0.2 | 16.01 ± 1.52 | −0.58 ± 1.53 |
| 0.2–0.4 | 11.17 ± 0.76 | 0.67 ± 0.77 |
| 0.4–0.7 | 7.33 ± 0.60 | −0.39 ± 0.45 |
| 0.7–1.0 | 5.17 ± 0.41 | −0.05 ± 0.37 |
| 1.0–1.5 | 3.94 ± 0.31 | −0.32 ± 0.20 |
| 1.5–2.0 | 2.68 ± 0.31 | 0.41 ± 0.21 |
| 2.0–3.0 | 1.68 ± 0.17 | 0.07 ± 0.15 |

The tangential signal is positive and declines with radius. The cross component
has one bin near 2σ at 1.5–2 Mpc; its seven-bin covariance-weighted null check is
χ²=10.59 for 7 degrees of freedom (p=0.157) under this estimated covariance.
These checks are conditional on the current sample and covariance assumptions.

## Numerical and pipeline checks

- Reproduced the saved 800-lens target profile from its verified source rows
  before restricting the comparison to common support.
- Random and cluster profile bins use the same SACC n(z), B-mode polarity,
  cosmology, and radial edges. The random-subtracted and jackknife profile
  point estimates agree exactly in the saved output.
- The 14×14 paired covariance is finite and positive definite; its maximum
  asymmetry is only floating-point roundoff relative to its matrix scale.
- 51 equal-area sky cells were occupied. The covariance deletes lens and
  source endpoints from both matched samples in the same cell.
- The stack summary records 9,010,324 exported rows, 3,100 centers, 435 common
  lenses, random subtraction, and spatial block covariance.

A later 2.5-degree versus 5-degree block-width comparison found unchanged point estimates but 0.820–1.067 error ratios; see [`wl_random_shape_jackknife_sensitivity_assessment_2026-11-05.md`](wl_random_shape_jackknife_sensitivity_assessment_2026-11-05.md).

## Limitations and next gates

- `science_measurement=false` remains authoritative. Random-center validity
  still depends on full-aperture depth, mask, PSF, and sky matching; those gates
  are not established by this row export.
- The 5-degree equal-area cells may cross disconnected HSC footprints. Repeat
  with alternative block widths and validated contiguous survey regions before
  using the covariance for a science claim.
- HSC Y3 tomographic n(z) is applied, but its calibration uncertainty is not
  propagated. Reduced-shear correction is also absent.
- The cross component and source-selection response need follow-up systematics
  checks. The 1.5–2 Mpc cross residual is not a detection, but should remain in
  the audit record.
- This weak-lensing result does not settle the separate ICL–satellite
  association. That branch remains exploratory and needs its own image, sky,
  PSF, source-injection, and common-support gates.

## Reproduction

The result was produced by
`run_wl_random_shape_stack.py` from the hash-complete export and existing HSC
Y3 n(z), cluster stack, and cluster source-row products. The post-export
watcher, submitter, stacker, checkpoint, row summary, and run logs remain in
Astrafek or the documented `/gpfs/kjhan/Astrafek/data` result tree. The figure
is reproducible with `plot_wl_random_shape_stack.py`.
