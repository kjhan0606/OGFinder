# WL source-density boost count audit (2026-11-05)

## Purpose

The exploratory HSC S19A weak-lensing stack subtracts the matched-random
`DeltaSigma` profile, but it has not applied a source-density boost correction.
This offline audit asks whether eligible source-pair counts around the 435
common-support cluster lenses exceed those around the 2,175 matched random
centers (five per lens). It reuses the saved 5° and 2.5° stack summaries; it
does not reread images or source rows, and it makes no HSC request.

## Estimator

For each proper-radius bin, the count-based candidate is

\[
B_{\rm count}(R)=\frac{N_{\rm lens}(R)/435}{N_{\rm random}(R)/2175}.
\]

The counts come from the paired cluster-minus-random profile and use the same
shape, B-mode, tomographic positive-efficiency, magnitude, resolution, and
radial cuts as the saved lensing estimator. All 2,175 planned random centers
remain in the denominator, including 129 centers with no exported shape rows.
The two block-width summaries must have identical samples and radial pair
counts; the resulting point ratios are therefore block-width independent.

## Result

| Proper radius (Mpc) | Lens source pairs | Random source pairs | Count-based boost candidate | Count excess |
|---:|---:|---:|---:|---:|
| 0.1–0.2 | 7,680 | 21,759 | 1.765 | 76.5% |
| 0.2–0.4 | 25,699 | 84,530 | 1.520 | 52.0% |
| 0.4–0.7 | 62,878 | 231,552 | 1.358 | 35.8% |
| 0.7–1.0 | 88,842 | 355,842 | 1.248 | 24.8% |
| 1.0–1.5 | 204,909 | 861,144 | 1.190 | 19.0% |
| 1.5–2.0 | 277,768 | 1,191,783 | 1.165 | 16.5% |
| 2.0–3.0 | 775,269 | 3,376,139 | 1.148 | 14.8% |

The count excess is strongest at small radii and remains positive through 3 Mpc.
This is large enough that source-density dilution could materially affect the
lensing profile if the count ratio traces cluster-member contamination. It can
also reflect residual differences in aperture coverage, depth, masks, PSF, or
sky between cluster and random positions.

## Interpretation and next gate

These are point-only counts: no uncertainty interval is available because the
saved summary retains total counts but not radial-by-spatial-block counts. Do
not multiply the published exploratory `DeltaSigma` profile by these ratios.
The random-center full-aperture matching gates and a paired spatial covariance
for the density ratio must first be verified. A follow-up Astrafek analysis
should recompute per-bin source counts and lensing-weight sums by matched
spatial block from the existing approved local pair rows; it should report
both count-weighted and lensing-weighted boosts, and propagate their covariance
jointly with random subtraction. Photo-z dilution calibration and reduced
shear remain separate open corrections. Keep `science_measurement=false`.

## Reproduction and provenance

```bash
python3 analysis/icl_wgl_membership/hsc_debug/audit_wl_source_density_boost.py \
  --output /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_source_density_boost_count_audit_2026-11-05_01
```

The machine-readable output is in
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_source_density_boost_count_audit_2026-11-05_01/`.
The summary SHA-256 is `25b935584ed114f59c2ce22363f60a60e2a7536a093ca8e322b40ba84adbabc1`,
and the seven-bin TSV SHA-256 is
`4c95f293e7f86f5a29533b49a9365ebf05f8335972dd6b5d3267e809b3dca81a`. The
Astrafek script SHA-256 is
`84697e3a7379591d05170c5a7b9f86d43ce7890c3d723857ea219b9f77994549`.
Input summary hashes are recorded in `summary.json` and match the documented
5° and 2.5° products.
