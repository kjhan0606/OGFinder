# Joint source-density and weak-lensing diagnostic

This figure brings the corrected HSC S19A source-density diagnostics and the
random-subtracted weak-lensing profile into one view, alongside their spatial
delete-one correlation matrix. It is descriptive only: the boost correction
has not been applied and `science_measurement=false` remains in force.

![Joint source-density boost and weak-lensing covariance](/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_source_density_boost_spatial_2026-11-05_04/joint_boost_delta_sigma_covariance.png)

## Reading the figure

- The inner-bin count-density ratio is `1.7648 ± 0.0638`, while the corrected
  S19A response-adjusted estimator-weight ratio is `1.1387 ± 0.0483`. The
  count ratio alone is therefore not the appropriate multiplicative factor
  for this calibrated shear estimator.
- Across 0.1–3 Mpc, the estimator-weight ratio stays between 1.125 and 1.163
  at the seven saved radii. The count and shape-weight ratios vary more
  strongly with radius. No one of these curves is applied as a correction.
- The tangential profile is positive over the displayed range and declines
  with radius; the cross component is consistent with zero at the shown
  delete-one errors. This is a descriptive stack, not a detection statement.
- The joint matrix contains all 42 saved parameters: four boost definitions
  over seven radial bins, plus tangential and cross `DeltaSigma` over seven
  bins each. It makes clear that adjacent radial estimates are strongly
  correlated and that boost parameters share substantial covariance with one
  another. The full covariance is poorly conditioned, so direct GLS inference
  needs a documented regularization and block-scale sensitivity study.

## Reproduction and provenance

Run the Astrafek plotting script with the archived corrected diagnostic:

```bash
python3 analysis/icl_wgl_membership/hsc_debug/plot_wl_boost_delta_sigma_joint.py \
  --summary /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_source_density_boost_spatial_2026-11-05_04/summary.json \
  --output /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_source_density_boost_spatial_2026-11-05_04/joint_boost_delta_sigma_covariance.png \
  --block-width-deg 2.5
```

- Summary SHA-256: `04fff9521b4b9f78eca88ea3c72620a3cc4eeadefdeb5fc3bec238cbfef66a70`
- Figure SHA-256: `2ecfcae2a312386a09504463251c6d328598c287857624a2063c7bfb2d100d2e`
- Plot script: `plot_wl_boost_delta_sigma_joint.py`
- Sample: 435 common-support lenses, 2,175 random centers, 159 spatial blocks
  at 2.5°; saved-stack point and pair-endpoint covariance reproduction gates
  both pass.

The plot reads the saved JSON only and makes no HSC request. It does not
change the science gate, perform a fit, or alter the saved stack.
