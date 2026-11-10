# S23B fixed-52 random-center weighting covariance (2026-11-11)

Status: exploratory estimator sensitivity; `science_measurement=false`.

## Scope and estimators

The offline Astrafek run used the previously approved 20-file S19A shape-row export, pinned HSC Y3 SACC `n(z)`, and existing radial-support/QA products. It made no CAS request, network request, download, or image-pixel read. The source-bin selection is 2–4 over seven 0.1–3 proper-Mpc annuli.

The fixed sample contains 52 lens sets, each with target support and at least four of five matched control centers supported in every annulus. For each set, the target profile is recalibrated separately. The pooled estimator combines supported random source pairs by their lensing weights. The equal-center estimator calibrates each random center separately and gives each positive center equal weight within its annulus. Both estimators then use an equal outer mean over the same 52 sets.

For each 2.5° and 5° equal-area RA/sin(Dec) deletion, a source-center pair is removed if either endpoint lies in the deleted cell. Target and control calibration is recomputed. Both estimators use a common lens-set intersection across all seven bins and keep the four-control floor. The equal-minus-pooled covariance uses these same paired deletions.

## Tangential profile

Values and errors are in units of `10^13 Msun Mpc^-2`. The pooled point estimate exactly reproduces the pinned fixed-52 profile. The paired-difference error is the uncertainty of the estimator contrast under the same approximate deletion scheme.

| R proper (Mpc) | Pooled ΔΣt (5°) | Pooled error (2.5° / 5°) | Equal-center ΔΣt (5°) | Equal-center error (2.5° / 5°) | Equal−pooled (5°) | Contrast error (2.5° / 5°) |
|---:|---:|---:|---:|---:|---:|---:|
| 0.1–0.2 | 23.5538 | 5.2856 / 4.4101 | 24.2044 | 5.3995 / 4.2079 | +0.6506 | 1.1416 / 1.1015 |
| 0.2–0.4 | 6.0478 | 2.3783 / 2.0390 | 5.9004 | 2.3630 / 1.8968 | −0.1474 | 0.5283 / 0.4594 |
| 0.4–0.7 | 4.0030 | 1.4803 / 1.3399 | 4.0850 | 1.5545 / 1.3870 | +0.0821 | 0.3479 / 0.2617 |
| 0.7–1.0 | 3.5211 | 1.5141 / 1.2548 | 3.3786 | 1.4789 / 1.2341 | −0.1425 | 0.1541 / 0.1317 |
| 1.0–1.5 | 2.5206 | 1.0271 / 0.9687 | 2.5513 | 1.0076 / 0.9436 | +0.0307 | 0.1159 / 0.0890 |
| 1.5–2.0 | 0.6750 | 0.6897 / 0.5796 | 0.7726 | 0.6828 / 0.5758 | +0.0976 | 0.0809 / 0.0723 |
| 2.0–3.0 | 1.2551 | 0.5563 / 0.5352 | 1.3007 | 0.5590 / 0.5415 | +0.0456 | 0.0612 / 0.0384 |

The largest absolute estimator shift is about 1.35 paired jackknife errors across these bins. This is not an equivalence test: the covariance is conditional on approximate sky cells, and the sample is selected for radial support. Neither weighting rule is adopted as final from this comparison.

## Cross-component diagnostics

These nominal quadratic-form diagnostics use the cross-component covariance submatrix. Their block-width sensitivity prevents a stable null assessment.

| Cell width | Pooled cross profile χ² / dof (p) | Equal-center cross profile χ² / dof (p) | Equal−pooled cross contrast χ² / dof (p) |
|---:|---:|---:|---:|
| 2.5° | 6.349 / 7 (0.4997) | 6.834 / 7 (0.4464) | 7.427 / 7 (0.3859) |
| 5.0° | 13.344 / 7 (0.0642) | 14.065 / 7 (0.0500) | 13.673 / 7 (0.0573) |

These are conditional diagnostics, not confirmatory tests or detection probabilities. The numerical covariance gate passed for pooled, equal-center, and estimator-difference profiles at both widths: all six matrices are 14×14, full rank, and PSD. The 2.5° grid has 64 occupied cells and a minimum of 49/52 common sets after deletion; the 5° grid has 28 cells and a minimum of 47/52. These equal-area cells are not validated contiguous HSC survey regions, so the numerical gates do not establish survey-footprint covariance.

## Assessment and remaining gates

The pooled and equal-center point profiles are similar relative to their paired contrast errors under this partition. The cross-component zero-profile diagnostic changes from p≈0.45–0.50 at 2.5° to p≈0.05–0.06 at 5°, so the cross null depends on the spatial partition. The weighting comparison therefore does not resolve the weak-lensing interpretation.

Target/control aperture depth, mask, PSF and sky matching remain unverified. Full `n(z)` uncertainty, photo-z dilution, reduced shear, magnification, cluster-member dilution/boost, selection-response uncertainty, and an empirical contiguous-footprint covariance also remain open. There was no PSF convolution or completeness correction in this run. Keep `science_measurement=false`; do not infer a halo mass or detection from these profiles.

## Reproducibility

- Runner: [`analyze_s23b_wl_fixed_equal_control_spatial_covariance.py`](analyze_s23b_wl_fixed_equal_control_spatial_covariance.py)
- Runner SHA-256: `f31f7f30fc9d244d20dd6efaccaa51a9b5b833cf76bbb8c05b310aecb7039fc7`
- Fixed lens-set IDs SHA-256: `18312f27139b3750abdbf71543f22d9b7b18b2d3b24cbd29ed686a377524bb7d`
- Approved source-row plan SHA-256: `7bfb998149f30a50b5ee7a537bb703ed302b7ea269a3296c35204a836e75aabc`
- Output directory: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_s23b_equal_random_center_spatial_covariance_2026-11-11_02/`
- Summary SHA-256: `00c8427ef6810b78f0900da1d515e3fca2ed77233bd012ebbeddbf5bc5e5c033`
- Generated assessment SHA-256: `63272539aa6d2beff5eec1c9b2b2d86b84d1db4c9bbc4075db1ce7f99215bd74`
- Cross diagnostics TSV SHA-256: `88e1628443cd93100068bcaf4463902e928f2316368985c38c01bb40c7ff5d2a`
- Estimator profiles TSV SHA-256: `d194d838f6505bf1211202fddf8ce68512827fca74d0ebe241d15cd3e026604f`
- Covariance matrices TSV SHA-256: `64f20a7d69877b738ac6b9734b8d3e08d737c99ee416da9c10a7a1fe12c9ad09`
