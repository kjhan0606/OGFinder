# S23B CAMIRA selection scenarios (2026-11-05)

## Scope and provenance

This is an offline catalog-pool comparison using the already received HSC
S23B CAMIRA Wide star-mask catalog, local S20A CAMIRA Wide catalog, and the
current 800-lens manifest. It submitted no HSC query and downloaded no new
data. The candidates are not established as image-supported or usable ICL / weak-lensing lenses; `science_measurement=false`.

The audit excludes S23B rows within 1 arcsec of S20A, then applies the existing
richness-first greedy selector with a 6-arcmin minimum separation from the
current 800-lens sample and previously selected additions. The baseline
reproduced the exact 1,333 candidate IDs in the pinned prior audit.

| Lens selection | Eligible S23B rows | After S20A dedup | Greedy additions | New 5° cells | Total 5° cells |
|---|---:|---:|---:|---:|---:|
| `0.30 ≤ z_cl < 0.55`, BCG spectroscopic consistency required | 2,288 | 1,440 | 1,333 | 22 | 100 |
| `0.25 ≤ z_cl < 0.65`, BCG spectroscopic consistency required | 3,462 | 2,124 | 1,983 | 25 | 103 |
| `0.20 ≤ z_cl < 0.70`, BCG spectroscopic consistency required | 4,119 | 2,515 | 2,339 | 26 | 104 |
| `0.30 ≤ z_cl < 0.55`, cluster redshift only | 3,574 | 2,412 | 2,184 | 24 | 102 |
| `0.25 ≤ z_cl < 0.65`, cluster redshift only | 5,751 | 3,756 | 3,420 | 28 | 106 |
| `0.20 ≤ z_cl < 0.70`, cluster redshift only | 7,033 | 4,579 | 4,137 | 31 | 109 |

## Interpretation and next sample decision

Relaxing the redshift range while retaining the BCG spectroscopic-consistency
cut increases the catalog-only additions from 1,333 to 1,983 or 2,339, but
adds only 3 or 4 new 5-degree cells. The resulting population spans a wider
redshift interval and cannot be pooled with the baseline without explicit
redshift calibration, matched selection, and redshift-stratified results.

Dropping BCG spectroscopic consistency yields larger candidate counts, but it
also changes the lens-redshift validation. Those scenarios should not be used
as the primary sample until lens-redshift calibration and the intended
selection function are established. The 5-degree cells are coarse footprint
bins, not proven independent covariance blocks. Even the largest scenario has
31 new cells and 109 total cells; it does not by itself meet the 92-block
single-endpoint or 180-block multiple-comparison design targets recorded in
the power assessment.

Keep the present selection as the replication baseline. If sample expansion
is pursued, the broader spectroscopically consistent selections are secondary
catalog candidates for a prespecified, redshift-stratified extension. Before
any image download, first audit PDR3 image and mask availability, depth, PSF,
full-aperture coverage, member completeness, and matched-control feasibility.
Those checks and the actual per-target/control support determine the
downloadable sample; catalog counts do not.

## Reproduction

Run on `lageunha` CPU:

```bash
python3 analysis/icl_wgl_membership/hsc_debug/audit_s23b_camira_selection_scenarios.py
```

The runner verifies the S23B catalog, S20A catalog, current manifest, and
prior baseline hashes before analysis. It stores six candidate ledgers and a
machine summary at
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_s23b_camira_selection_scenarios_2026-11-05_02/`.
The summary SHA-256 is
`e0fddda7ecc973327842796278dd01b8e2ae3de33f2d07244590cfb085507ff6`; the
Astrafek runner SHA-256 is
`dcbe3a106ea472facea0b317779a4f86ccee182d2fec8b098c030a9339c282c5`.

All conclusions are catalog-selection ceilings. No new HSC service query,
image download, shape-row export, or science measurement occurred.
