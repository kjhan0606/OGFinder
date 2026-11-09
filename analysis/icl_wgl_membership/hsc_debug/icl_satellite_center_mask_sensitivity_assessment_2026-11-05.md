# ICL–satellite association by central-mask eligibility (2026-11-05)

## Scope and method

This offline sensitivity uses the 94 existing cluster/filter pairs with both
ICL 200–300 kpc and satellite 50–75 kpc profile support. Within each filter,
the pairs are separated by the saved central-mask eligibility flag. The
estimand is the partial Pearson correlation between raw target-minus-control
ICL surface brightness and raw satellite contrast, controlling linearly for
cluster redshift and richness. Each subgroup receives its own percentile
bootstrap over connected 5-degree target/control sky blocks (5,000
replicates). A paired comparison resamples the full filter sample's connected
blocks once per replicate and applies the same selected blocks to both
eligibility groups. This is a conditional subgroup diagnostic; it does not
model profile availability or correct fluxes.

## Results

| Filter | Group | Pairs | Connected blocks | Partial `r` | 95% block interval |
|---|---|---:|---:|---:|---:|
| `i` | All supported | 26 | 13 | 0.079 | [−0.183, 0.477] |
| `i` | Center mask eligible | 8 | 7 | −0.421 | [−1.000, 1.000] |
| `i` | Center mask ineligible | 18 | 10 | 0.234 | [0.014, 0.802] |
| `i2` | All supported | 68 | 33 | 0.293 | [0.039, 0.500] |
| `i2` | Center mask eligible | 15 | 11 | 0.307 | [−0.069, 0.792] |
| `i2` | Center mask ineligible | 53 | 32 | 0.312 | [−0.070, 0.551] |

| Filter | Eligible minus ineligible `r` | Paired 95% block interval | Valid paired draws |
|---|---:|---:|---:|
| `i` | −0.655 | [−1.665, 0.632] | 4,306 / 5,000 |
| `i2` | −0.005 | [−0.516, 0.693] | 4,977 / 5,000 |

Only 8 `i` and 15 `i2` pairs pass the saved central-mask eligibility flag.
The `i2` point estimates are nearly unchanged between eligible and ineligible
groups, but both subgroup intervals are broad and include zero. The `i`
eligible estimate changes sign, with an interval spanning nearly the full
correlation range. The paired difference intervals include zero in both
filters. Only 86% of paired `i` resamples meet the minimum subgroup-size rule,
so its difference interval is especially weak. No group-difference p-value is
reported.

## Interpretation

The current data cannot establish that central-mask eligibility explains or
does not explain the ICL–satellite association. The `i2` point estimate is
similar in the two groups, while the small clean-center sample lacks precision.
The `i` result is especially unstable. This confirms that central-mask
eligibility is a real sample-size constraint and that the present association
is not yet robust enough for a science claim. It does not justify relaxing the
mask or pooling the subgroups as if they were exchangeable.

The analysis does not correct central light, PSF wings, sky residuals,
extinction, pixel covariance, satellite completeness, or selection. The 94
pairs are a selected support intersection rather than a random subset of the
800-lens parent. All machine products retain `science_measurement=false`.

## Reproduction and provenance

```bash
python3 analysis/icl_wgl_membership/hsc_debug/analyze_icl_satellite_center_mask_sensitivity.py
```

Astrafek runner SHA-256:
`305323e2499a7524a45e92de2ef206882942276e965cef6a6232a25acdf1f3c0`.
Output directory:
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_center_mask_sensitivity_2026-11-05_02/`.
Summary SHA-256:
`8b5ad1230b78816a68b8d8ea810f64df32740f1ec199caac86f6be65c7e0e9ee`.
The analysis used existing pair summaries and coordinate manifests only; no
HSC query or image read occurred.
