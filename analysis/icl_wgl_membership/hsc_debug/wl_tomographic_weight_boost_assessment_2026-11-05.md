# Tomographic weak-lensing weight boost and source-bin sensitivity (2026-11-05)

## Scope and validation

This Astrafek analysis reuses the approved local HSC S19A lens/source and
matched-random source-row exports. It includes the same 435 common cluster
lenses and 2,175 random centers (five per lens), with four HSC Y3 source bins
and paired endpoint jackknife covariance at 5-degree and 2.5-degree widths.
No HSC request or image read occurred. The source-bin pair counts sum exactly
to the saved full-stack counts in every radial bin, and the four-bin sum
reproduces the saved full-sample `DeltaSigma` profile. The output remains
diagnostic with `science_measurement=false`; no boost correction is applied.

## Redshift dependence of estimator-weight excess

The exact S19A response-adjusted estimator-weight ratio (cluster centers over
matched random centers, per planned center) differs strongly by source
tomographic bin. Values below use the 2.5-degree paired jackknife; 5-degree
point estimates are identical and give the same qualitative result.

| Radius (Mpc) | z-bin 1 | z-bin 2 | z-bin 3 | z-bin 4 |
|---|---:|---:|---:|---:|
| 0.1–0.2 | 4.605 ± 0.274 | 1.149 ± 0.055 | 0.973 ± 0.056 | 1.042 ± 0.072 |
| 0.4–0.7 | 2.231 ± 0.114 | 1.157 ± 0.040 | 1.108 ± 0.035 | 1.149 ± 0.038 |
| 1.5–2.0 | 1.298 ± 0.058 | 1.160 ± 0.037 | 1.129 ± 0.033 | 1.093 ± 0.033 |
| 2.0–3.0 | 1.260 ± 0.049 | 1.150 ± 0.035 | 1.121 ± 0.031 | 1.092 ± 0.031 |

In every radial bin, z-bin 1 differs from each of bins 2–4 in paired
delete-one-cell contrasts. All 21 contrasts remain below Holm-adjusted
`p=9.4×10⁻⁴` at 2.5 degrees and `p=1.23×10⁻³` at 5 degrees, conditional on
the selected cells and a normal approximation to the jackknife contrasts.
The 5-degree 56-parameter cross-tomographic covariance has rank 51 and is
singular as expected from its number of delete-one replicates; it was not
inverted. Individual contrast variances are computed directly from that
covariance. The 2.5-degree covariance is full rank, with minimum/maximum
eigenvalue ratio `2.16×10⁻⁶`.

The pronounced low-bin excess is consistent with enhanced source selection or
photo-z leakage where the low-redshift source population overlaps the lens
population. It does not identify the cause: member contamination,
magnification, depth, mask, PSF, sky response, and redshift-calibration errors
can all affect these ratios. Do not treat any of the ratios as a dilution
correction.

## Excluding the lowest source bin

As a sample-sensitivity check, the `DeltaSigma` stack was recomputed after
excluding HSC Y3 source bin 1, retaining bins 2–4 and the same 435 lenses,
2,175 random centers, and matched spatial deletions. The full four-bin sum
first reproduced the saved stack profile. The table reports tangential
`DeltaSigma` in `10^13 M_sun Mpc^-2`; differences are (bins 2–4) minus (all
bins), and the uncertainty is from paired jackknife replicas.

| Radius (Mpc) | All bins | Bins 2–4 | Paired change ± jackknife σ |
|---|---:|---:|---:|
| 0.1–0.2 | 16.01 | 18.03 | +2.025 ± 0.595 |
| 0.2–0.4 | 11.17 | 12.27 | +1.095 ± 0.200 |
| 0.4–0.7 | 7.33 | 7.55 | +0.225 ± 0.091 |
| 0.7–1.0 | 5.17 | 5.20 | +0.027 ± 0.069 |
| 1.0–1.5 | 3.94 | 3.91 | −0.025 ± 0.046 |
| 1.5–2.0 | 2.68 | 2.66 | −0.016 ± 0.031 |
| 2.0–3.0 | 1.68 | 1.71 | +0.035 ± 0.026 |

The seven-bin tangential profile remains positive after excluding bin 1:
`χ²=653.38/7` (2.5-degree blocks) and `485.88/7` (5-degree blocks) against
zero under the estimated covariance. The bin-1-exclusion cross profile is
consistent with zero (`p=0.264` and `p=0.054`, respectively). The paired
change in the tangential profile is nonzero under the same jackknife
approximation (`χ²=57.43/7`, `p=4.90×10⁻¹⁰` at 2.5 degrees;
`χ²=45.59/7`, `p=1.05×10⁻⁷` at 5 degrees), while the cross-profile change is
null-consistent (`p=0.372` and `p=0.441`). These are sensitivity statistics,
not corrected detection significances: n(z) calibration uncertainty,
reduced shear, contiguous-footprint covariance, and full-aperture
cluster/random depth, mask, PSF, and sky matching remain unpropagated.

The result shows that the tangential profile is not driven solely by the
lowest source bin, but its inner amplitude changes materially when that
overlapping bin is included. Report bins 2–4 as a robustness result alongside
the fiducial profile; do not select one as the corrected measurement until
the source-redshift and random-field systematics are calibrated.

## Reproduction and provenance

Run on `lageunha` CPU from the repository root:

```bash
python3 analysis/icl_wgl_membership/hsc_debug/analyze_wl_tomographic_weight_boost.py
```

The run uses the already approved local row exports and makes no HSC network
request. Astrafek runner SHA-256:
`497b1ac799f261beae91dd24bb8a2a61da65817f2943ba5411d95173d3e5d1673`.
Machine summary:
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_tomographic_weight_boost_2026-11-05_04/summary.json`
(SHA-256 `47419ac229ffec9aab2a001f629252627126a456812b01cf11908ad44122fe14`).
The summary contains all four-bin profiles, cross-bin covariance,
multiple-comparison-adjusted contrasts, source-bin-exclusion profiles, and
paired covariance matrices.
