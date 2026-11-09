# Redshift-matched richness split in the HSC cluster-only weak-lensing stack — 2026-11-03

## Scope and method

The earlier 400/400 CAMIRA richness split had slightly different parent-sample
lens-redshift distributions (median `z=0.442` and `0.427`). This follow-up asks
whether its cluster-only `DeltaSigma` contrast persists after matching lens
redshift. It reuses only the already-approved HSC S19A shape-row export and HSC
Y3 tomographic `n(z)`; no HSC query, new data transfer, or source-row export was
performed.

The 800 parent lenses were first split by richness, sorting on `(richness,
lens_id)` and assigning 400 to each group. Before reading the shear profiles,
the script formed a maximum-cardinality 1:1 minimum-`|delta z|` match between
the two groups, with `|delta z| <= 0.02`. This produced 392 candidate pairs.
The verified shape-row export has no rows for some parent lenses, so the stack
kept only pairs where both lenses have exported source-row support. The final
balanced sample is 251 lens pairs (251 lenses per richness group); 141 candidate
matches were removed by this support requirement. Mean/median/maximum matched
`|delta z|` in the final sample are 0.00720/0.00253/0.01999. Mean lens redshift
is 0.4259 vs 0.4198. Median CAMIRA richness remains separated (17.314 vs
26.558), as intended.

Both stacks use the same HSC Y3 four-bin SACC `n(z)`, the existing raw B-mode
exclusion value 0, seven proper-radius bins, and 5-degree equal-area sky cells.
Each delete-one replicate removes pairs whose lens or source falls in the
deleted cell and recomputes all tomographic and shear calibrations. The matched
profile difference uses paired spatial-cell jackknife covariance across the
two stacks.

## Results

`DeltaSigma` and its paired jackknife error are in units of
`10^12 M_sun Mpc^-2`.

| Proper R (Mpc) | 400/400 split upper−lower | Redshift-matched upper−lower |
|---:|---:|---:|
| 0.1–0.2 | 15.13 ± 21.71 | 35.28 ± 24.90 |
| 0.2–0.4 | 41.58 ± 12.26 | 40.73 ± 14.38 |
| 0.4–0.7 | 35.38 ± 8.71 | 32.79 ± 8.74 |
| 0.7–1.0 | 27.92 ± 6.76 | 24.67 ± 6.63 |
| 1.0–1.5 | 15.32 ± 4.88 | 15.31 ± 5.56 |
| 1.5–2.0 | 17.23 ± 3.56 | 17.61 ± 3.57 |
| 2.0–3.0 | 9.25 ± 2.12 | 9.59 ± 2.28 |

The radial contrast remains similar from 0.2 to 3 Mpc after the redshift match.
The innermost-bin contrast grows, with a broad paired error. A conditional
seven-bin tangential-profile Wald diagnostic is `chi2=62.66`, rank 7, nominal
`p=4.44e-11`; this is not a detection significance because the profile is still
cluster-only and the covariance is based on approximate equal-area cells.

## Interpretation and open gates

Matching lens redshift alone does not remove the richness-dependent profile
contrast in this selected export-supported subset. This strengthens the case
for checking a richness dependence after the random-lens step, while the loss
of 141/392 candidate matches to source-row support shows that the analyzed
subset is strongly conditioned on catalog support. It cannot establish a
causal richness–mass relation or a weak-lensing detection.

The stack still lacks matched random-lens subtraction and matched-field
variance. The cells are not validated contiguous HSC survey regions. Lens
selection, source-selection differences, shear-calibration uncertainty,
uncertainty in `n(z)`, reduced shear, and richness–mass calibration are not
propagated. The separate random-center aggregate-count plan remains
`prepared_not_submitted`; this analysis did not submit it. Keep
`science_measurement=false`.

## Reproduction and provenance

- Analysis: [`audit_wl_redshift_matched_richness_split.py`](audit_wl_redshift_matched_richness_split.py)
- Machine-readable output: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_redshift_matched_richness_2026-11-03_01/summary.json`
- Summary SHA-256: `f050c2c90236a4215d951cb939bb65c145b98d6526dea90a88007beff2f152fb`
- Analysis SHA-256: `97f948fce84dfe22c292d287157e3e0b5c388f8283dc518f772258b942ce4bdd`
- Input row export: 3,069,544 verified pairs, 800 parent lenses; existing plan, batch, count-query, sample-manifest, and SACC hashes are embedded in `summary.json`.
- Execution: lageunha CPU, 49 occupied cells; no HSC network access.
