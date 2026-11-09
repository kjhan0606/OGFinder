# Radial ICL–satellite relations (2026-11-06)

## Analysis

This analysis uses the existing 94-pair supported ICL/satellite sample and the
cluster-level target-minus-azimuthal-control satellite excess profiles in seven
radial bins. It evaluates two complementary predictors at each radius:

1. Continuous ICL excess at 200–300 kpc, scaled to one within-filter standard
   deviation.
2. High-minus-low ICL median split within each filter.

Both models adjust for standardized redshift and log richness. The 5-degree
connected target/control pair blocks were recomputed jointly over both filters.
Restricted HC2 Webb six-point wild-cluster bootstrap-t used the same random
weight per physical block across all 28 tests. Holm and single-step max-T
corrections cover the complete family of both relation types, filters, and
radii. Profile support ranges from 22 to 68 clusters per bin.

## Results

The smallest raw bootstrap p-value is 0.0369 for the continuous i2 relation at
50–75 kpc, with a slope of +1.96 nJy arcsec⁻² per 1σ ICL contrast. Its
max-T-adjusted p-value is 0.770. The smallest max-T p-value in the 28-test
family is 0.190 (continuous i at 5–10 kpc); its slope is −240 nJy arcsec⁻² per
1σ ICL contrast. All 28 Holm-adjusted p-values are 1.0. Thus neither apparent
outer positive i2 trend nor inner negative i trend is supported after the
joint spatial bootstrap and multiplicity correction.

Three continuous i inner-bin CR1 intervals (5–10, 10–15, and 15–20 kpc) lie
below zero pointwise. The corresponding restricted-bootstrap raw p-values are
0.052, 0.068, and 0.136; the max-T values are 0.190, 0.246, and 0.315. Treat these
pointwise intervals as descriptive; they do not establish an inner
anti-correlation. The i2 continuous slopes are positive across the seven bins,
but their intervals and familywise tests do not support a detection. The
median-split profile likewise yields no familywise-significant contrast.

These are selected-sample exploratory results from 26 i and 68 i2 pairs. They
do not establish causality or generalize to all 800 cutouts. No images or new
HSC data were requested in this calculation. `science_measurement=false`.

## Reproduction and limitations

Runner:
[`analyze_icl_satellite_radial_relations.py`](analyze_icl_satellite_radial_relations.py),
SHA-256 `40170376c87b69ff07abf9972f70020a798a33dcdd17639cec14de944582ca3b`.
Plotter:
[`plot_icl_satellite_radial_relations.py`](plot_icl_satellite_radial_relations.py),
SHA-256 `53b3f1267f13a7d14b72057349856a6fee9a8dc26b81ff8c573f7026a7662dd5`.
The figure SHA-256 is `933fb1aa4f27c2bca76335e5756c346e458133328d2aa41d8050169f5b569ca9`.
Inputs, exact method definitions, seed `20261107`, 10,000 replicates,
per-bin estimates, p-values, and output checksums are recorded in
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_radial_relations_joint_2026-11-06_01/summary.json`.
The table and max-T null draws are stored alongside the summary.

This post-hoc 94-pair support intersection is affected by profile availability,
mask, control, and earlier analysis selection. The connected 5-degree blocks
are approximate, not HSC field jackknife regions. Intervals omit uncertainty
from photometric covariance, source Poisson noise, coadd weights, sky transfer,
PSF wings, satellite completeness, and injection response. The 28-test
adjustment controls only this declared profile family and cannot correct for
the full prior analysis history.
