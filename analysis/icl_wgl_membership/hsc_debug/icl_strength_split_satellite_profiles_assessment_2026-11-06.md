# ICL-strength split satellite profiles (2026-11-06)

## Question and data

This exploratory analysis asks whether the cluster-mean satellite excess
profile differs between high- and low-ICL systems. It uses the fixed set of 94
supported target/control pairs (26 in `i`, 68 in `i2`) and the saved HSC
satellite profiles from the 800-cutout analysis. Within each filter, the ICL
signal at 200–300 kpc defines a median split. Each radial-bin profile is
joined at cluster level; the number of usable clusters changes with radius,
and no missing profile is imputed. The estimand is the equal-cluster mean
target-minus-azimuthal-control satellite excess in nJy arcsec⁻².

The high-minus-low profile contrast is estimated with OLS adjusted for
standardized redshift and log richness. Uncertainty uses a CR1 cluster-robust
covariance and restricted HC2 Webb six-point wild-cluster bootstrap-t. Spatial
blocks are recomputed jointly from target/control sky coordinates across both
filters so a shared physical region receives the same bootstrap weight in all
14 filter-by-radius comparisons. Holm and single-step max-T adjustments cover
those 14 contrasts.

## Results

| Filter | Radius (kpc) | N | Adjusted high−low (nJy arcsec⁻²) | CR1 95% interval | max-T p |
|---|---:|---:|---:|---:|---:|
| i | 5–10 | 26 | −338.4 | [−730.6, 53.8] | 0.520 |
| i | 10–15 | 25 | −78.0 | [−226.5, 70.5] | 0.961 |
| i | 15–20 | 24 | −33.8 | [−122.4, 54.9] | 0.999 |
| i | 20–30 | 24 | +1.3 | [−46.6, 49.2] | 1.000 |
| i | 30–40 | 23 | −4.9 | [−30.7, 21.0] | 1.000 |
| i | 40–50 | 22 | −2.3 | [−17.3, 12.7] | 1.000 |
| i | 50–75 | 26 | +4.7 | [−1.0, 10.5] | 0.580 |
| i2 | 5–10 | 66 | +208.2 | [−54.8, 471.2] | 0.704 |
| i2 | 10–15 | 66 | +60.2 | [−24.6, 145.1] | 0.821 |
| i2 | 15–20 | 65 | +29.4 | [−3.7, 62.4] | 0.565 |
| i2 | 20–30 | 63 | +18.1 | [−2.1, 38.2] | 0.559 |
| i2 | 30–40 | 59 | +4.4 | [−0.7, 9.6] | 0.602 |
| i2 | 40–50 | 61 | +10.4 | [−5.4, 26.3] | 0.877 |
| i2 | 50–75 | 68 | +1.3 | [−1.9, 4.5] | 1.000 |

All 14 CR1 intervals include zero; all Holm-adjusted p-values are 1.0, and the
smallest max-T adjusted p-value is 0.520. The i2 point estimates are positive
at all radii, but this selected-sample pattern is not resolved against the
spatially clustered uncertainty. The i-band estimates change sign with
radius and are substantially less precise. These results do not establish an
ICL–satellite relation or a physical difference between high- and low-ICL
clusters.

## Audit trail and limitations

The runner is
[`analyze_icl_strength_split_satellite_profiles.py`](analyze_icl_strength_split_satellite_profiles.py)
(SHA-256 `2f6208e7f34504f3dc7a5c76d2d202b4c34ae2bdf2ed05ec9e655d8e6fac512b`).
The plotter is
[`plot_icl_strength_split_satellite_profiles.py`](plot_icl_strength_split_satellite_profiles.py).
Its SHA-256 is `dfed877d3140fabffa5c3c64dad537293248f98a71c723bd7f3effe0f338f036`;
the rendered figure is
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_icl_strength_stacked_profiles_2026-11-06_01/icl_strength_split_satellite_profiles.png`
(SHA-256 `43f9843839a276bf3464c50b73346bce7f6ee598399168844bb924aff9ca32d`).
Input and output hashes, seed (`20261106`), 10,000 bootstrap replicas, methods,
and the `science_measurement=false` flag are recorded in
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_icl_strength_stacked_profiles_2026-11-06_01/summary.json`.
The profile table and bootstrap max-T null draws are saved alongside it.

This is a post-hoc split on the same selected 94-pair support intersection;
it is not an independent validation and cannot undo previous analysis
selection. The 5-degree connected blocks are approximate sky partitions.
Intervals omit remaining photometric covariance, source Poisson noise,
coadd-weight uncertainty, sky-transfer uncertainty, PSF-wing residuals,
satellite completeness, and injection-transfer uncertainty. No new HSC data
were requested or downloaded for this step. `science_measurement=false`.
