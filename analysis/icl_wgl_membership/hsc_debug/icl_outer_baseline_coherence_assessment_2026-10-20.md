# ICL inner-to-outer baseline coherence — 2026-10-20

## Result

Within the existing 5%-support-selected target/control pairs, clusters with a
larger 400–500 kpc target-minus-control offset also tend to have a larger inner
target-minus-control profile at 100–300 kpc. This pattern remains after linear
adjustment for cluster redshift, CAMIRA richness, and the inner and outer common
area fractions. In 50–100 kpc, the i-band interval is broad and includes zero;
the i2 estimate is weak. The six annulus/filter comparisons are exploratory,
not multiplicity-adjusted tests.

| Filter | Inner annulus (kpc) | Pairs / connected blocks | Adjusted partial r (95% block interval) | Adjusted slope (95% block interval) |
|---|---:|---:|---:|---:|
| i | 50–100 | 157 / 18 | 0.124 (−0.075, 0.274) | 0.510 (−0.405, 1.170) |
| i | 100–200 | 177 / 18 | 0.338 (0.091, 0.464) | 0.494 (0.123, 0.762) |
| i | 200–300 | 189 / 17 | 0.455 (0.283, 0.549) | 0.477 (0.314, 0.595) |
| i2 | 50–100 | 279 / 44 | 0.144 (0.005, 0.263) | 0.536 (0.020, 1.023) |
| i2 | 100–200 | 338 / 44 | 0.355 (0.211, 0.483) | 0.463 (0.295, 0.630) |
| i2 | 200–300 | 368 / 44 | 0.421 (0.287, 0.541) | 0.583 (0.329, 0.955) |

The adjusted slope is the change in inner surface brightness per unit change
in the pair's 400–500 kpc difference after removing the listed covariates.
Values use equal pair weights. Intervals are percentile ranges from 5,000
resamples of connected 5-degree target/control spatial blocks.

## Interpretation and limits

The outer baseline and inner profile are not independent across pairs, especially
for the 100–300 kpc annuli. Subtracting a per-pair outer annulus therefore
changes the profile in a way that is associated with the inner signal in this
selected sample. This coherence can arise from residual sky structure, genuine
radially extended ICL, cluster-to-cluster amplitude variation, or a combination;
the present data do not distinguish these causes. The 400–500 kpc annulus is not
certified blank sky.

The input rows inherit the 5% common-support floor, and the selected pair set
differs by annulus. The block bootstrap represents only the chosen connected
5-degree target/control grouping; it omits full pixel/coadd covariance, PSF
wings, and residual-sky transfer uncertainty. These diagnostics do not correct
the ICL profile or establish an ICL detection. `science_measurement=false`.

Reproduce with
[`analyze_icl_outer_baseline_coherence.py`](analyze_icl_outer_baseline_coherence.py).
The verified result is
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_outer_baseline_coherence_2026-10-20_02/summary.json`
(SHA-256 `1b243b33da40cb7f29abc3fbe2fa65785f8120cdd9925841a3dd253d8a319722`);
its compact six-row output is `coherence_metrics.jsonl` (SHA-256
`5b9d7080f8f7a7d8d85b0ceecc136460b62249abe6543628cadcd826321fdd00`). The JSON
stores all input hashes, the code hash, and the 5,000 bootstrap results. No HSC
query or download was made.
