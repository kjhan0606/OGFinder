# Nonlinear redshift/richness adjustment of the ICL--satellite relation (2026-11-06)

## Question and method

The existing 94-pair analysis linearly controls cluster redshift and CAMIRA
richness. This post-hoc sensitivity asks whether that model form drives the
ICL 200--300 kpc versus satellite 50--75 kpc partial correlation. On the same
26 `i` and 68 `i2` pairs, four models add redshift squared, richness squared,
their interaction, or all three terms to the linear baseline. Terms are
centered and scaled before constructing the polynomial basis. Each model uses
the same 5,000 connected 5-degree target/control block resamples within its
filter. A single-step max-T adjustment covers eight model changes (four per
filter).

## Results

| Filter | Baseline r | Added terms | Adjusted r | Change (pointwise 95% interval) | max-T p |
|---|---:|---|---:|---:|---:|
| `i` | 0.079 | redshift² | 0.005 | −0.073 [−0.299, +0.215] | 0.917 |
| `i` | 0.079 | richness² | 0.094 | +0.015 [−0.101, +0.158] | 0.999 |
| `i` | 0.079 | redshift × richness | 0.107 | +0.028 [−0.096, +0.168] | 0.989 |
| `i` | 0.079 | all quadratic terms | 0.036 | −0.043 [−0.375, +0.299] | 0.998 |
| `i2` | 0.293 | redshift² | 0.294 | +0.001 [−0.040, +0.030] | 1.000 |
| `i2` | 0.293 | richness² | 0.296 | +0.003 [−0.033, +0.033] | 0.999 |
| `i2` | 0.293 | redshift × richness | 0.298 | +0.005 [−0.033, +0.038] | 0.996 |
| `i2` | 0.293 | all quadratic terms | 0.305 | +0.012 [−0.088, +0.057] | 0.993 |

All paired-change intervals include zero. The `i2` point estimate is stable
across these four smooth model-form changes. The `i` estimate moves more, but
its 26-pair intervals are broad. This rules out neither unmeasured confounding
nor selection effects; it only shows that these particular quadratic controls
do not resolve the existing ambiguity.

## Interpretation limits

This is a post-hoc sensitivity family on a selected, small sample. The block
partition is approximate and is not a validated contiguous-field covariance.
The analysis uses no image pixels and does not address source-Poisson noise,
coadd weights, sky transfer, PSF wings, or profile completeness. No profile
correction is applied, and `science_measurement=false` remains.

## Reproducibility

- Astrafek runner: `analyze_icl_satellite_nonlinear_adjustment.py`, SHA-256
  `3c9e24ea50bc400a3527e6be8b7e534e85aa93d87cb2c52853233903e75bc49b`.
- Result directory:
  `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_nonlinear_adjustment_2026-11-06_01/`.
- Summary SHA-256:
  `4777d93629348a8849a74af8aa3217b4553ff27142d7802724b8c6091a5c332c`.
- Joined-pair table SHA-256:
  `2f97e22a45e1613cd0602041ab02d8d55ebd4777ef032be626b39100c035f697`.
- Input support table SHA-256:
  `79f7c4c6b18b29ad837e91da4a4359433b09a564cad62efaf60b4cce655685da`.
- No HSC query, download, or pixel read occurred.
