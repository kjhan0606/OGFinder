# 2-D ICL stack mask-radius sensitivity — 2026-10-19

## Result and scope

The matched target-minus-control 2-D stack was repeated with circular source-mask radii of 2.5, 3.5, and 5.0 Kron radii, holding the selected 800-cluster manifest, image data, pair matching, pixel grid, and equal-pair estimator fixed. Each run stacked 783 eligible target/control pairs (302 in `i`, 481 in `i2`). All three runs passed the independent stack verifier and radial-profile evaluation. The comparison contains 12 matched filter/radius rows and is recorded in `icl_2d_mask_sensitivity_assessment_2026-10-19_01/summary.json`; the table checksum is `dd07ffa209492a21df5e4c19ec192c868cd22027be4abe87e04a9b02836de601`.

This is an exploratory mask-sensitivity check. It does not promote the stack to a science measurement: all outputs retain `science_measurement=false`.

## Radial comparison

Values are the equal-pair target-minus-control mean in nJy/arcsec². `n` is the number of pairs with annular support; coverage is median common-area coverage among those pairs.

| Filter | Radius (kpc) | K=2.5: Δ (n, coverage) | K=3.5: Δ (n, coverage) | K=5.0: Δ (n, coverage) |
|---|---:|---:|---:|---:|
| i | 50–100 | 41.39 (235, 0.150) | 34.86 (215, 0.074) | 28.59 (153, 0.020) |
| i | 100–200 | 12.92 (277, 0.145) | 10.17 (272, 0.098) | 7.27 (260, 0.031) |
| i | 200–300 | 2.08 (283, 0.117) | 1.13 (281, 0.082) | 2.64 (275, 0.034) |
| i | 300–400 | 4.36 (295, 0.103) | 3.70 (295, 0.070) | 0.78 (292, 0.027) |
| i | 400–500 | 4.34 (299, 0.101) | 1.35 (298, 0.071) | 1.65 (295, 0.029) |
| i2 | 50–100 | 44.59 (378, 0.190) | 38.08 (357, 0.084) | 36.85 (239, 0.031) |
| i2 | 100–200 | 10.75 (447, 0.177) | 9.16 (444, 0.110) | 9.46 (435, 0.037) |
| i2 | 200–300 | 3.18 (467, 0.154) | 2.43 (465, 0.106) | 2.17 (457, 0.043) |
| i2 | 300–400 | 0.14 (476, 0.130) | −0.10 (475, 0.090) | −0.11 (472, 0.038) |
| i2 | 400–500 | −0.30 (480, 0.127) | −0.66 (480, 0.089) | −1.29 (477, 0.040) |

The innermost 0–50 kpc bins are dominated by the central source and mask geometry and should not be read as diffuse-light constraints. Their pair counts are only 46–141 in `i` and 80–240 in `i2`; the mean changes non-monotonically as the masks grow.

## Interpretation

The positive 50–200 kpc target-minus-control excess persists in both filters under all three mask sizes. Its amplitude falls as the masks grow: in `i`, 50–100 kpc declines from 41.39 to 28.59 and 100–200 kpc from 12.92 to 7.27 nJy/arcsec² between K=2.5 and K=5.0. In `i2`, the corresponding values decline from 44.59 to 36.85 and from 10.75 to 9.46. This shows material sensitivity to source-mask choice even where the sign is stable.

The available area also contracts sharply. At K=5.0 the median common coverage is only 2–4% in the 50–500 kpc bins, compared with 10–19% at K=2.5. The annulus-support pair count remains high at larger radii, but it does not compensate for the very small common pixel area. Those outer means therefore remain vulnerable to local sky structure, correlated coadd noise, residual source wings, and control mismatch.

At 300–500 kpc, `i2` is consistent with a near-zero offset at all mask sizes (−1.29 to +0.14 nJy/arcsec²). The `i` outer offset is positive for each mask choice, but falls from about 4.3 at K=2.5 to 0.8–1.7 at K=5.0; this band-dependent behavior, together with low coverage, does not establish an ICL detection.

The reported pair SEM is descriptive pair-to-pair scatter only. Mask changes alter both the pixels and the set of pairs contributing to each radial bin. These comparisons are not independent tests, do not account for full pixel covariance, and do not include an observationally validated sky/PSF transfer model. They cannot be interpreted as detection significance or as an ICL–satellite covariance measurement.

## Reproducibility

The underlying run directories are:

- K=2.5: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_2d_matched_stack_2026-10-19_02/`
- K=3.5: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_2d_matched_stack_2026-10-19_k3p5/`
- K=5.0: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_2d_matched_stack_2026-10-19_k5p0/`

K=5.0 and its independent verification/evaluation completed in Slurm job 428633 (39 min 8 sec; peak RSS 574,056 KiB). The reproducible comparison command is `python3 analysis/icl_wgl_membership/hsc_debug/compare_icl_2d_mask_sensitivity.py`. It refuses incomplete or unverified inputs and writes the comparison TSV, two diagnostic plots, and a hash-linked summary beneath the data directory. The figures are `mask_sensitivity_i.png` and `mask_sensitivity_i2.png` in the assessment output directory.
