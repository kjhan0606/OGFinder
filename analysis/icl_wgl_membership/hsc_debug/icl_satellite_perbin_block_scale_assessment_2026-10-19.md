# Full per-annulus ICL–satellite block-scale sensitivity — 2026-10-19

## Scope and provenance

This is the block-width sensitivity test for the primary 94-row per-annulus
matched sample: 26 i and 68 i2 cluster/filter pairs. It reruns the
redshift/richness-adjusted partial-Pearson bootstrap with 2.5°, 5°, and 10°
connected target/control sky cells, retaining the fixed rows and seed rule.
No sample reselection, HSC query, or download was performed.

The Astrafek implementation is
[`audit_satellite_icl_perbin_block_scale.py`](audit_satellite_icl_perbin_block_scale.py).
Its loader verifies the association, ICL-annulus, and matched-control
coordinate inputs before analysis. At 5°, it reproduces the archived spatial
bootstrap to within 1e-12 for 28 numeric comparisons. The output is
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/satellite_icl_perbin_block_scale_sensitivity_2026-10-19_01/`.
The 94-row association SHA-256 is
`9a6befa87a81257f8faf1d34fa7c019a91d64484fbfa64136bd93496dfc7b58c`; the
result TSV SHA-256 is
`f6ced42df5b578cbf71d91e4426ab3a8810eba9b66f7f10a581c31720bfe63c5`, and
the summary SHA-256 is
`32bfe7067c62d4a865c47e955b54522b44a12b09656bed8ee01019b7c37c9527`.

## Results

| Filter and selection | Cell width | Pairs | Blocks | Partial r | Spatial-block 95% interval |
|---|---:|---:|---:|---:|---:|
| i2, all support | 2.5° | 68 | 47 | 0.293 | [−0.003, 0.504] |
| i2, all support | 5° | 68 | 33 | 0.293 | [0.042, 0.497] |
| i2, all support | 10° | 68 | 23 | 0.293 | [0.057, 0.491] |
| i2, center eligible | 2.5° | 15 | 14 | 0.307 | [−0.150, 0.702] |
| i2, center eligible | 5° | 15 | 11 | 0.307 | [−0.059, 0.812] |
| i2, center eligible | 10° | 15 | 10 | 0.307 | [−0.071, 0.805] |

The point estimate is invariant because each run uses the same 68 i2 pairs.
The interval for all supported i2 pairs narrowly crosses zero at 2.5° and
excludes zero at 5° and 10°. Thus interval-based interpretation depends on the
assumed spatial dependence scale. The center-eligible subset is too small to
constrain the relation at any tested width. The i-band results are also
preserved in the machine-readable summary.

These are unadjusted bootstrap intervals from post-hoc sensitivity widths, not
confirmatory significance tests. The separate familywise block-wild analysis
found no adjusted significant association. Pixel/coadd covariance, PSF wings,
sky transfer, extinction, membership completeness, and selection uncertainty
remain unresolved. `science_measurement=false` remains in force.

The 5-test suite passed. The output records exact input hashes, each width's
pair/block counts, valid and invalid bootstrap draws, and the baseline
reproduction check.
