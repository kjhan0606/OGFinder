# Supported-sample warp variance: Lanczos-5 mapping sensitivity

## Scope and reproducibility

This offline sensitivity uses the same 94 matched ICL–satellite target/control
pairs (188 planned HSC-I centers) and the already downloaded visit-warps. It
does not contact HSC or read any new remote data. Four centers remain unavailable
after repeat HTTP 404 responses; the analysis includes 184 centers, giving 24
complete `i` pairs and 64 complete `i2` pairs for the paired comparisons.

The HSC PDR3 processing documentation describes fifth-order Lanczos resampling
for CCD-to-warp processing:
[HSC PDR3 processing documentation](https://hsc-release.mtk.nao.ac.jp/doc/index.php/processing__pdr3/).
This analysis applies a separable Lanczos-5 kernel only to map downloaded warp
variance samples onto coadd pixel centers. It propagates the diagonal input
variance as the sum of squared normalized kernel coefficients times the input
variances. This does not undo or model the off-diagonal covariance already
introduced when CCD pixels were resampled into the warp.

## Results: 200–300 proper-kpc annulus

| Filter | Role | Centers | Median coadd/predicted variance ratio | Median strict common-mask coverage |
|---|---:|---:|---:|---:|
| `i` | target | 26 | 1.00707 | 69.9% |
| `i` | control | 24 | 1.00760 | 70.1% |
| `i2` | target | 67 | 1.00639 | 71.9% |
| `i2` | control | 65 | 1.00611 | 69.0% |

Using connected 5-degree sky-block bootstrap intervals, the paired median
target-minus-control ratio is `+0.0000322` in `i` (24 pairs, 12 blocks; 95%
interval `[-0.001928, +0.001551]`) and `+0.000458` in `i2` (64 pairs, 33
blocks; `[-0.000110, +0.001831]`). Both intervals include zero.

Compared center by center with the previous nearest-neighbor mapping, the
median change in the primary variance ratio is `-0.000245` in `i` and
`-0.000161` in `i2`. Primary coverage is lower by a median 0.27 and 0.24
percentage points, respectively, because the full Lanczos support must fit
inside valid warp pixels. Strict shared-mask coverage is substantially lower
than primary coverage and varies by center.

## Interpretation limits

The Lanczos-5 mapping changes the variance ratios only slightly relative to
nearest-neighbor mapping, and the paired target/control intervals still include
zero. This is a resampling sensitivity, not an exact HSC covariance
reconstruction. The calculation treats warp-pixel errors as diagonal and does
not reproduce coadd weights, flux-scale factors, clipping, mask interpolation,
or the CCD-to-warp covariance. It does not measure source-Poisson response,
residual-sky covariance, extended PSF wings, or total ICL-estimator uncertainty.
Do not use these values as a noise correction or promote the ICL profile;
`science_measurement=false` remains.

## Output provenance

- Result directory: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_supported94_warp_noise_lanczos5_2026-11-06_02/`
- Summary SHA-256: `16c8afbcc8287444860bacdcfd01038105bd51d5e82d61e44d8b82d8853e9f74`
- Per-center table SHA-256: `a0322cc88e63768d172ccc8da52ce5148bef1a4407232b53e93e00bc9ae0f06f`
- Runner SHA-256: `0df14f6c6be1d858a38437c3ae7300adce81fba859047c8d29d8fe4d07232db9`
- GPU used: no. HSC network used: no. New data requested: no.
