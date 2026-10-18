# ICL sky-estimator association contrast — 2026-10-19

This Astrafek analysis compares the existing pixel-SFD ICL estimator with the
plane-first estimator on the same matched cluster/control pairs. The reusable
implementation is
[`audit_icl_sky_estimator_delta.py`](audit_icl_sky_estimator_delta.py), with
regression checks in
[`test_audit_icl_sky_estimator_delta.py`](test_audit_icl_sky_estimator_delta.py).

## Fixed comparison

The input table contains 74 unique matched pairs: 53 in `i2` and 21 in `i`.
The preregistered primary band is `i2`; `i` is descriptive. Within each band,
the calculation residualizes both ICL and satellite surface-brightness deltas
against redshift and richness, then computes their partial Pearson
correlation. It holds the center-SFD-corrected satellite delta and pair sample
fixed while comparing `icl_delta_pixel_sfd_njy_arcsec2` with
`icl_delta_pixel_sfd_after_plane_njy_arcsec2`. Uncertainty comes from 5,000
resamples of connected joint target/control 5-degree sky blocks. Each pair of
estimators receives the identical block draw in every replicate.

## Results

| Band | Pairs | Connected blocks | Pixel-SFD partial r (95% block CI) | Plane-first partial r (95% block CI) | Paired difference, plane-first minus pixel-SFD (95% CI) |
|---|---:|---:|---:|---:|---:|
| i2 (primary) | 53 | 28 | 0.239 (0.013, 0.425) | 0.124 (−0.168, 0.375) | −0.115 (−0.329, 0.063) |
| i (descriptive) | 21 | 13 | 0.112 (−0.217, 0.527) | 0.114 (−0.207, 0.592) | 0.002 (−0.146, 0.223) |

For the primary `i2` sample, the point estimate falls by about 0.115 under the
plane-first sky estimator. The paired interval includes zero, so this dataset
does not establish that the two estimator correlations differ. The broad
intervals also do not establish a stable ICL–satellite association under the
plane-first estimator. This is an estimator-sensitivity result, not evidence
for or against a physical association.

## Reproduction and provenance

Run from the OGFinder repository on `syntax` with the approved existing local
GPFS input; the calculation is small and CPU-only:

```bash
python3 analysis/icl_wgl_membership/hsc_debug/audit_icl_sky_estimator_delta.py \
  --work /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_sky_estimator_delta_2026-10-19_01
```

Machine-readable results and provenance are in
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_sky_estimator_delta_2026-10-19_01/summary.json`.
The manifest records the input-table and program SHA-256 hashes, seed, method,
and `science_measurement: false`. The input is the existing
`icl_pixel_extinction_2026-10-19_03/pixel_extinction_satellite_association.tsv`;
this run retrieved no HSC data.

The result conditions on the selected 74 pairs and chosen spatial blocks. It
does not model upstream sky-transfer errors, PSF wings, mask leakage,
selection/completeness, injection/recovery, or full image-level covariance.
The ICL science gate remains closed pending those checks.
