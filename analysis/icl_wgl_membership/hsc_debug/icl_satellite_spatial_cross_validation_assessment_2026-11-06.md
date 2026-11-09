# Spatial cross-validation of the ICL–satellite relation (2026-11-06)

## Question and result

Does adding the 200–300 kpc ICL contrast improve prediction of the saved
50–75 kpc satellite-light excess for spatial blocks left entirely out of model
fitting, after controlling redshift and CAMIRA richness?

Astrafek compared two ordinary least-squares models within each filter:
redshift plus log richness, and the same predictors plus ICL contrast. Each
connected target/control 5-degree block was held out in turn. Predictors were
standardized using only the training fold. The summary reports a descriptive
block-bootstrap interval for the paired change in held-out squared error.

| Filter | Pairs / blocks | OOF MSE: nuisance only | OOF MSE: plus ICL | Relative MSE reduction from ICL | Mean squared-error gain (95% block interval) |
|---|---:|---:|---:|---:|---:|
| `i` | 26 / 13 | 99.33 | 105.27 | −6.0% | −5.94 [−18.53, 1.67] |
| `i2` | 68 / 33 | 45.13 | 43.57 | +3.5% | +1.57 [−2.98, 6.81] |

The direction matches the earlier in-sample pattern only weakly: ICL adds a
small held-out improvement in `i2`, while it worsens held-out error in `i`.
Both intervals include zero. These selected pairs therefore show no clear
evidence that the ICL measurement improves prediction in unseen spatial
blocks. This is a robustness check, not evidence against a physical relation;
the test has low sample size and inherits the profile-selection and measurement
limitations.

## Reproducibility and limitations

- Astrafek runner: [`analyze_icl_satellite_spatial_cross_validation.py`](analyze_icl_satellite_spatial_cross_validation.py)
- Runner SHA-256: `a068c2ac0330a012964596e7bafc11a3c3d2be7aca293866a0e90b4024b8b510`
- Output summary: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_spatial_cv_2026-11-06_01/summary.json`
- Summary SHA-256: `7059c7b196467afb581f0edf99ce69a65b9e132979f7a3284be8182c6cc1ceb7`
- Fold table SHA-256: `a2606fdc08e490428f8c53feb0cc063790b1631ba7a825faaddc7490db3d0eb8`
- Bootstrap: 10,000 resamples of spatial blocks; seeds 20261109 (`i`) and 20261110 (`i2`).
- No HSC network access or image-pixel reads occurred. `science_measurement=false` remains.

Connected equal-area sky blocks approximate survey fields. The input is a
post-hoc intersection of 94 supported cluster/filter pairs (26 `i`, 68 `i2`),
not the full 800-cluster parent. The result is predictive rather than causal
and does not remove uncertainty from source Poisson noise, coadd weights, sky
transfer, extended PSF wings, or missing satellite profiles. The bootstrap
interval is descriptive and does not account for the full history of
exploratory analyses.
