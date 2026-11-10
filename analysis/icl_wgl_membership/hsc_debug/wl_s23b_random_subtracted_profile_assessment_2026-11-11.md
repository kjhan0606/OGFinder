# S23B matched-sample weak-lensing profile — 2026-11-11

## Scope and result

This exploratory Astrafek analysis formed an HSC S19A shape-source `DeltaSigma` profile around the approved S23B CAMIRA lenses and matched random centers. It read the verified 80-lens/400-control source-row export and the pinned HSC Y3 SACC tomographic `n(z)`. The calculation ran on `lageunha` CPU. It made no CAS query, network request, additional download, or image-pixel read. It did not fit halo masses. The result remains `science_measurement=false`.

The primary policy retained 74 lens sets with nonzero calibrated lens support and at least four positive-support controls out of five (367 control centers). It used source bins 2–4, the existing HSC per-source `m/c` calibration, radial-bin responsivity and edge-selection response, and the verified B-mode-mask convention (raw 0 excluded). The estimator uses the weak-shear approximation. The random profile pools the selected positive-support control centers; lens sets with five supported controls therefore contribute one more control aperture than sets with four. This weighting choice needs an equal-per-lens sensitivity before final inference.

## Primary profile

The errors below are paired endpoint equal-area jackknife standard errors for 5-degree RA/sin(Dec) cells. A pair is deleted when either its lens/random center or source lies in the deleted cell. The cells are not validated contiguous HSC survey regions.

| Projected radius (proper Mpc) | Random-subtracted tangential `DeltaSigma` (`10^13 Msun Mpc^-2`) | Cross (`10^13 Msun Mpc^-2`) | Lens/control pair rows |
|---:|---:|---:|---:|
| 0.1–0.2 | 23.845 ± 3.432 | 5.391 ± 4.288 | 649 / 3,158 |
| 0.2–0.4 | 6.835 ± 1.803 | -2.183 ± 2.634 | 2,579 / 11,874 |
| 0.4–0.7 | 3.214 ± 1.141 | 0.153 ± 1.328 | 6,577 / 32,574 |
| 0.7–1.0 | 4.111 ± 1.102 | -1.168 ± 0.831 | 10,221 / 49,451 |
| 1.0–1.5 | 2.630 ± 0.688 | 0.645 ± 0.560 | 24,390 / 118,910 |
| 1.5–2.0 | 1.503 ± 0.382 | 0.955 ± 0.487 | 33,924 / 162,249 |
| 2.0–3.0 | 0.762 ± 0.400 | 0.455 ± 0.401 | 95,724 / 462,270 |

The tangential profile is positive and generally declines with radius. Its conditional zero-profile Wald diagnostic is `chi2=221.85` for rank 7 (`p=2.67e-44`), but this is not a detection probability: the covariance uses approximate equal-area cells and does not include the calibration and survey-selection systematics below. Pair rows are center/source pairs, not globally unique galaxies.

## Sensitivity and null checks

| Case | Lens sets / controls | Occupied cells | First / last tangential bin (`10^13 Msun Mpc^-2`) | Nominal cross-null `p` |
|---|---:|---:|---:|---:|
| At least 4/5 controls; bins 2–4; 5° | 74 / 367 | 36 | 23.845 / 0.762 | 0.0741 |
| At least 4/5 controls; bins 1–4; 5° | 74 / 367 | 36 | 22.565 / 0.824 | 0.1254 |
| 5/5 controls; bins 2–4; 5° | 71 / 355 | 36 | 23.262 / 0.845 | 0.0341 |
| At least 4/5 controls; bins 2–4; 2.5° | 74 / 367 | 80 | 23.845 / 0.762 | 0.3182 |

The stricter 5/5-control case has a nominal cross-null discrepancy, while the 2.5-degree case does not. The null behavior is therefore sensitive to control selection and block width; it is not a clean, stable systematic check. All quoted `p` values are nominal conditional diagnostics with no correction for this sensitivity grid. Restoring source bin 1 shifts the innermost tangential point from 23.845 to 22.565 (about 5.4%); the low-bin `n(z)` and estimator-weight concerns from earlier HSC audits remain.

## Assessment and next gates

This is useful as a reproducible pilot profile for the 80-lens source-row export, not as a calibrated mass or ICL–satellite relation. Before making physical claims:

1. Reweight the random stack so each matched lens set contributes equally, and compare it with the pooled-control result.
2. Compare the fiducial SACC `n(z)` with the available S19A-specific DESI calibration curves/envelopes. Those files do not contain joint posterior draws, so they cannot provide a full uncertainty propagation; photo-z dilution uncertainty is currently unpropagated.
3. Replace or validate the equal-area cells against contiguous HSC survey regions, and quantify mask, depth, PSF, and sky matching between lens and random apertures.
4. Assess reduced shear, magnification, cluster-member dilution/boost, and selection-response uncertainty before any mass fit or ICL–satellite interpretation.

Until these gates pass, retain `science_measurement=false`. No additional remote data are needed for the first two sensitivity analyses; later access to new data would require its own scoped approval.

## Reproducibility

- Runner: [`analyze_s23b_wl_random_subtracted_profile.py`](analyze_s23b_wl_random_subtracted_profile.py)
- Runner SHA-256: `7e7ca33190f59fb9daeee91994231127c58facfa9d452554758a43fccf046c65`
- Approved row-plan SHA-256: `7bfb998149f30a50b5ee7a537bb703ed302b7ea269a3296c35204a836e75aabc`
- HSC Y3 SACC `n(z)` SHA-256: `a28f9e2e088e92d96d2d87083a6eeeac4c0c84c6b48e033bd3d1dc1bf0d8958f`
- Radial/tomographic QA summary SHA-256: `eb73aa0365813d880d3edf60b12e4039f2575f8a97eab8868fc918121a34d445`
- Per-center support table SHA-256: `e630442cbaeb054763897a46d68d3b4baea7c0884e80ac22f46ece906583f48e`
- Output directory: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_s23b_random_subtracted_profile_2026-11-11_01/`
- `summary.json` SHA-256: `2c899b0090dd94a34f1fd6b0cb111e3d1ada696da2178f98e9fd5b659e699edf`
- `assessment.md` SHA-256: `d152c22f21a244e9848f453523c0a346120e30f4897b9a427d65776aa804ba2a`
- `profile.tsv` SHA-256: `3dd1ef2b7f4ac7aa33d6b046ef9b528ae6c8818d1f3992037b6b70aa6e5548b2`
- `profile.png` SHA-256: `4c3593236d5a46dd504d0f54b3f2396e0a925460c27f3f879bc6fef7bef92ce9`
