# Random-subtracted weak-lensing n(z) mass sensitivity — 2026-11-06

## Scope

This offline analysis propagates the existing HSC Y3 source-redshift calibration
shift scenarios through the exact approved 435-lens common-support stack. For
each scenario it recalculates both the cluster and five-matched-random-center
profiles using the same HSC shape cuts and estimator, subtracts the random
profile, and refits the reduced-shear NFW model at fixed concentrations
`c=3,4,5`. It uses existing source shape rows, random-center shape rows, and
the local HSC Y3 SACC `n(z)` file only. No query, image read, or new data export
was made.

The scenarios are the local SACC fiducial curves, the prior published central
higher-redshift translations for source bins 3 and 4 (`+0.043`, `+0.050`),
and the four corners of their marginal 68% endpoint intervals
(`bin 3: +0.021/+0.063`; `bin 4: +0.038/+0.062`). They are constant
translations of two curves, not the published calibrated spline posteriors.

Runner: `hsc_debug/analyze_wl_random_subtracted_nz_mass_sensitivity.py`
(SHA-256 `8265a52dc71a35f9b6354f5550b4b4e31e7ffe1a5dc10238cfa204abc6f88c89`).
Machine output:
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_random_subtracted_nz_mass_sensitivity_2026-11-06_01/summary.json`
(SHA-256 `e6d6956607ffd0702678c475ab8e22188f579701ce29263157695b82fc13af2e`).

## Results

The zero-shift pair-level re-stack exactly reproduces the saved 5-degree
source-bin-1-excluded tomography profile: maximum fractional difference
`8.0e-16`. The sample contains 435 cluster lenses and 2,175 matched random
centers; source bin 1 is excluded. The baseline spatial-jackknife covariance
is held fixed for every scenario so this isolates profile and model response
to the `n(z)` shifts.

For the fixed-concentration `c=3` reduced-shear model, the local SACC fiducial
mass is `1.750e14 M☉` with conditional Δχ²=1 interval
`[1.626,1.877]e14 M☉`. The published central shift gives `1.706e14 M☉`,
2.60% below the fiducial. The four marginal-endpoint corners span
`1.689–1.724e14 M☉`, which is −0.96% to +1.07% relative to the shifted
central case. The endpoint spread is much smaller than the conditional
statistical interval for this fixed model.

Masses depend strongly on the assumed fixed concentration: the fiducial fits
are `1.750e14`, `1.604e14`, and `1.450e14 M☉` for `c=3,4,5`; their χ² values
for 3 degrees of freedom are 1.39, 12.58, and 32.54, respectively. Thus the
`c=4,5` profiles fit poorly under the saved covariance, and these results do
not constitute a mass calibration. The concentration/model dependence exceeds
the endpoint n(z) uncertainty in this grid.

## Limits

- The calibration grid is a rectangular set of marginal shifts; it has no
  posterior covariance between source bins and does not include the bin-1
  calibration shift.
- The covariance is fixed to the saved 5-degree, source-bin-1-excluded
  jackknife matrix. Changed estimator weights can change the covariance, which
  was not recomputed for each scenario.
- Fixed-concentration single-halo NFW fits do not represent the mass
  distribution, miscentering, magnification, two-halo term, or full
  concentration uncertainty. Conditional intervals omit these effects.
- The result does not validate random-aperture depth matching or the
  contiguous HSC field covariance, and it does not change the scientific gate.

All output keeps `science_measurement=false`. The weak-lensing signal remains
exploratory until the unresolved survey and calibration systematics are
resolved.
