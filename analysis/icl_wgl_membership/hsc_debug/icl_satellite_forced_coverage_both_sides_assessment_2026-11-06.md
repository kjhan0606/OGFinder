# 2026-11-06: 20% target/control coverage and ICL–satellite selection sensitivity

## Question

Does adding forced CAMIRA-position measurements at a 20% annular-coverage floor
change the selected ICL–satellite partial association? Is the comparison
affected by which ICL-supported pairs pass that floor?

## Inputs and method

The analysis uses only existing HSC cutouts and saved tables. It combines the
20% target / 20% rotated-control forced-profile run with the 200–300 proper-kpc
ICL annulus table, retaining every ICL-supported cluster/filter pair in the
selection table, whether or not it has a satellite profile. It compares the
full 20/20 sample, the common subset that also has a detected-only profile, and
that same common subset using detected-only measurements. Partial correlations
control redshift and log CAMIRA richness. Paired percentile intervals use
10,000 resamples of connected 5-degree spatial blocks defined on the supported
20/20 sample.

The 20% threshold is a post-hoc sensitivity choice. Forced CAMIRA coordinates
without SEP detections are not confirmed isolated galaxies. The measurements
are not a calibrated satellite stellar-envelope estimator.

## Results

| Filter | ICL-supported pairs | 20/20 selected | Excluded | Detected-only common | Full 20/20 partial r | Detected-only common r |
|---|---:|---:|---:|---:|---:|---:|
| i | 71 | 50 | 21 | 26 | −0.088 | +0.082 |
| i2 | 153 | 110 | 43 | 68 | +0.102 | +0.285 |

Selected and excluded pairs have similar median redshift and richness, but
selected pairs have substantially higher median ICL 200–300 kpc residual:
3.060 versus 1.116 nJy arcsec⁻² in `i`, and 4.446 versus 1.853 nJy arcsec⁻²
in `i2`. This is direct evidence that the low-coverage selection is related
to the ICL outcome; the full-sample correlation cannot be transported by
assuming selection is ignorable.

The paired decomposition is:

| Filter | Change component | Median Δr | 95% block-bootstrap interval |
|---|---|---:|---:|
| i | Sample expansion: full 20/20 minus common 20/20 | −0.007 | [−0.396, +0.432] |
| i | Forced measurements: common 20/20 minus detected-only | −0.175 | [−0.451, +0.113] |
| i | Total: full 20/20 minus detected-only common | −0.184 | [−0.471, +0.127] |
| i2 | Sample expansion: full 20/20 minus common 20/20 | −0.164 | [−0.339, +0.012] |
| i2 | Forced measurements: common 20/20 minus detected-only | −0.003 | [−0.161, +0.166] |
| i2 | Total: full 20/20 minus detected-only common | −0.164 | [−0.391, +0.065] |

All paired intervals include zero. The `i2` point-estimate shift mainly comes
from adding newly supported pairs, while the added forced measurements on the
common subset barely change its point estimate. In `i`, adding forced values
shifts the point estimate downward, but uncertainty is broad. These results
neither establish an ICL–satellite association nor show that none exists.

## Limits and disposition

- Coverage selection depends on image support and is associated with ICL
  residual, so this is not a missing-at-random correction or a representative
  parent-sample estimate.
- The 20% threshold was selected after prior sensitivities; no multiplicity
  correction across thresholds/specifications is included.
- Spatial-block intervals are exploratory and do not propagate measurement
  uncertainty from the ICL or forced satellite profiles.
- PSF wings, source-Poisson variance, total coadd covariance, sky transfer,
  blend contamination, and completeness remain unresolved.
- No correction is applied; `science_measurement=false` remains.

## Reproduction and artifacts

Runner: [`analyze_icl_satellite_forced_coverage_both_sides.py`](analyze_icl_satellite_forced_coverage_both_sides.py), SHA-256 `7e9ec1672e0ac000e5e28cf9810c9f98b1e538e60465664853b1dc9a0fe5609f`.

Machine output: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_forced_coverage_both_sides_2026-11-06_03/summary.json` (SHA-256 `b8b79a65e1c838e757fc4ee6c1e1fdf61340d72e989bf035ab29bc0a7a6a314a`); paired table SHA-256 `b613b53457dc9a85aa1be2b3ed0de0026a47d86fa32897af392859c34f5fabd2`.

The inputs are pinned in `summary.json`. No HSC query or download was made for
this analysis. The previously approved frame-metadata export (plan SHA-256
`49aa4b7601dd42243209d7b3b394a4dd3b4cbf00d87d28cbb9467553c29f4ce3`) had
already completed all 16 batches; this turn sent no duplicate request.
