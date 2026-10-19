# ICL–satellite association under independent pixel support — 2026-10-20

## Result

The fixed 94-row cluster/filter sample was recomputed with the 200–300 kpc ICL
target and control means formed from each frame's own valid pixels, using four
25 kpc radial bins weighted by their valid-pixel counts. Each independent
support profile passed its per-run verifier, and all 94 baseline rows joined
one-to-one at each SEP/Kron mask factor (`K=2.5, 3.5, 5.0`). The association
statistic is partial Pearson correlation between the ICL target-minus-control
signal at 200–300 kpc and the satellite target-minus-control signal at
50–75 kpc, residualized against redshift and richness. Its paired uncertainty
uses 5,000 shared resamples of connected 5-degree sky blocks.

For all supported pairs, the i2 estimate remains near 0.29 under all three
mask factors; the i-band estimate remains weak and near zero. Every paired
95% block-bootstrap interval for the change in partial correlation includes
zero in both filters. Thus the estimator change does not materially alter the
full supported-pair i2 association in this sensitivity test, while the i-band
data do not show a clear association.

| Filter | K | N | Common-support partial r | Independent-support partial r | Paired Δr (95% block interval) |
|---|---:|---:|---:|---:|---:|
| i2 | 2.5 | 68 | 0.293 | 0.304 | +0.012 [−0.066, +0.070] |
| i2 | 3.5 | 68 | 0.293 | 0.279 | −0.014 [−0.061, +0.051] |
| i2 | 5.0 | 68 | 0.293 | 0.291 | −0.001 [−0.100, +0.093] |
| i | 2.5 | 26 | 0.079 | −0.021 | −0.067 [−0.213, +0.059] |
| i | 3.5 | 26 | 0.079 | −0.021 | −0.100 [−0.235, +0.032] |
| i | 5.0 | 26 | 0.079 | 0.029 | −0.050 [−0.207, +0.112] |

The smaller center-eligible subset has 15 i2 and 8 i rows. Its estimated
independent-support changes vary substantially by mask and have broad block
intervals; the K=2.5 i2 row alone has a negative interval, so this subset does
not provide a stable cross-mask result. These small-subset estimates are
exploratory and should not be treated as a separate detection.

## Interpretation limits

This is estimator sensitivity within the original sample selected by the
common-support K=2.5 gate, not a re-selection or a population-level analysis.
The independent estimator uses different valid target and control pixels and
assumes azimuthal symmetry and matched populations. The block interval
captures the change between estimators within this selected sample; it does
not include mask-selection, extended-PSF, sky-transfer, extinction, or full
coadd-covariance uncertainty. No ICL or weak-lensing science measurement is
claimed (`science_measurement=false`).

The reproducible Astrafek analysis is
`analyze_satellite_icl_independent_support.py`, submitted by
`run_satellite_icl_independent_support.sbatch` as Slurm job 428768. The
machine-readable output is
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_independent_support_association_2026-10-20_01/`.
The summary SHA-256 is
`4c6ba99f4050d72163dfa00c37a3347978d12a99910cf9c55ed94fdb14ae356f`, and the
association table SHA-256 is
`83aeb0ea08476236df9abdec20abf0312ddf1c498d71b854a3edc97f0c7d34ff`.
