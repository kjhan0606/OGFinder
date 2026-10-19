# CAMIRA sample expansion pool assessment — 2026-10-19

## Scope and provenance

This offline diagnostic tests whether the existing local CAMIRA parent catalog
can substantially enlarge the 800-cluster sample under the current redshift
and BCG-spectroscopy cuts, or by relaxing them. It makes no CAS request, exports
no coordinates, downloads no images, and performs no science measurement.

The reproducible tool is
[`audit_camira_sample_expansion.py`](audit_camira_sample_expansion.py), with
unit tests in
[`test_audit_camira_sample_expansion.py`](test_audit_camira_sample_expansion.py).
It verifies the full-catalog checksum against `sample_manifest.json`, excludes
the current 800 IDs, applies a six-arcminute minimum center separation in
richness-priority order, and records newly occupied 5-degree sky cells. Results
are in
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/camira_expansion_pool_2026-10-19_01/summary.json`.
An independent repeat is in `_02/summary.json`; both result files have SHA-256
`27ada4545891397f0da0bd8ed4f410a32cda0dfeea197b088a3dcc204523c2ed` and use
catalog SHA-256
`d929985d7f6ee38c4435fc12099e9df1e5dbfd9ce9603cea319d5704129d22aa`.

## Candidate-pool results

| Selection scenario | Eligible candidates before separation | Greedy additions | New 5° cells | Total occupied cells |
|---|---:|---:|---:|---:|
| Existing `0.30 ≤ z < 0.55` plus consistent BCG spectroscopy | 69 | 57 | 0 | 78 |
| `0.25 ≤ z < 0.65` plus consistent BCG spectroscopy | 575 | 544 | 2 | 80 |
| `0.20 ≤ z < 0.70` plus consistent BCG spectroscopy | 856 | 813 | 2 | 80 |
| Existing redshift range, cluster redshift only | 843 | 792 | 2 | 80 |
| `0.25 ≤ z < 0.65`, cluster redshift only | 2,086 | 1,949 | 6 | 84 |
| `0.20 ≤ z < 0.70`, cluster redshift only | 2,802 | 2,596 | 6 | 84 |

The current selection is nearly exhausted: it offers only 57 further separated
catalog candidates and no new 5-degree cell. Broadening the redshift range while
keeping consistent BCG spectroscopy adds candidates, but only two additional
coarse sky cells. Removing the spectroscopy requirement increases the pool
substantially but adds only six such cells and changes the redshift-selection
and lensing-calibration problem. The cell count is a footprint diagnostic; it
does not establish independent covariance blocks or usable HSC image coverage.

## Recommendation and limits

Do not request more images under the unchanged selection expecting a decisive
increase in independent statistical power. The existing exploratory power
assessment requires about 180 independent connected spatial blocks for 80%
power at its conservative 30-test adjustment, while the current analysis has
33 i2 blocks. This catalog-only exercise does not show that the available
additional candidates supply those blocks. A wider sky footprint or a
pre-registered, scientifically justified redshift design is needed before
proposing another imaging batch.

The broader cluster-redshift scenarios are candidate-pool counts only. They
need a separate photo-z, source-selection, and weak-lensing calibration plan;
they must not be treated as approved targets. No new candidate has been
selected, queried, or approved for download here. The greedy six-arcminute
algorithm is deterministic but not a globally optimal maximum independent set.
Image availability, masks, depth, PSF, member completeness, and spatial-block
covariance for new candidates were not evaluated. `science_measurement=false`.
