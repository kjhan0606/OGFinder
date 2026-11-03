# Supported i2 hybrid-PSF grid assessment — 2026-11-03

## Scope and validation

The completed Slurm job 434556 ran 12 local-core/public-HSC-I-wing splice models on
the same 44 PSF-supported `i2` cluster frames. The models vary splice radius
(1.8, 2.4, 3.0 arcsec) and wing amplitude (0.50, 0.65, 0.80, 1.00), with a fixed
0.6-arcsec blend. No observed science coadd was modified, and no HSC service or new
data request was used.

Astrafek's independent summarizer recomputed the saved response arithmetic for all
528 frame-model records. All 12 configurations contain the identical 44 cluster IDs,
all 44 records per configuration validate, and the largest stored-versus-recomputed
arithmetic difference is zero. The reproduced result is
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_core_wing_hybrid_i2_supported44_2026-10-22_01/grid_audit_repro_2026-11-03.json`
(SHA-256 `eaf08ddc0cbaec3243f3e13935bdd5c960b1acdf450f682ca94f44fcbce13712`).

## Quantitative findings

- For a smooth diffuse unit injection measured with the original masks fixed, the
  median diffuse-to-ICL response changes only from 0.225702 to 0.225874 across the
  12 splice models (0.08% total envelope). Redetection-branch model medians are
  0.223556–0.224220, but the representative `r2.4_a0.65` frame distribution has
  q16–q84 of 0.1710–0.2533 and spans 0.00464–0.31791. The median's model stability
  therefore does not imply stable frame-level recovery after redetection.
- For satellite-envelope injection measured in the ICL annulus, fixed-mask medians
  are 0.58897–0.58939. Redetection medians are lower, 0.48460–0.49332; the median
  redetection-minus-fixed response is −0.10954 to −0.08790 across splice models.
  This is a detection/mask response difference, not an astrophysical correction.
- Satellite-envelope recovery in the satellite measurement is broad. In the central
  `r2.4_a0.65` model, fixed-mask q16–q84 is 0.010–0.898, while redetection is
  −0.613–1.187 and the full redetection frame range is −5.92 to 10.78. The
  diffuse-to-satellite redetection response also has a broad q16–q84 interval
  (−0.184 to +0.179) around a median 0.0206. Neither response supports a scalar
  transfer correction for the satellite observable.
- In the 200–300 kpc annulus, the median saved-mask coverage is 0.620 (q16–q84
  0.512–0.756). A separate noiseless 1 nJy point-source profile sends a median
  0.000724–0.001411 of its input flux into the unmasked annulus across splice
  models. The generic HSC-I wing is not matched to these `i2` frame positions or
  epochs, so this is a proxy sensitivity, not a measurement of actual coadd leakage.

## Interpretation and next gate

The fixed-mask transfer of smooth diffuse light to the ICL statistic is numerically
stable across this particular model grid. Re-detection and satellite-related
responses remain highly frame dependent, including negative and greater-than-unity
values. The 44-frame subset is PSF-supported, not a random draw of all 68 `i2`
association rows. The q16–q84 ranges describe this selected subset and are not
population confidence intervals.

This grid does not include source Poisson noise, coadd covariance, realistic sky
transfer, matched-control residuals, or full-sample completeness. It does not validate
the public HSC-I wing as an `i2` PSF, provide a PSF correction, or test the observed
ICL–satellite association itself. Continue to mark `science_measurement=false`,
`psf_approved=false`, `observational_validation_passed=false`, and
`flux_covariance_available=false`.

The next useful ICL validation is a noise-aware, matched target/control injection
test through the same mask and sky estimator, followed by block-level uncertainty
propagation. Preserve the current observed-stack products and keep the synthetic
transfer result separate from measured ICL profiles.
