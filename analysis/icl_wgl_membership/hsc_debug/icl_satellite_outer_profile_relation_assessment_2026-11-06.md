# ICL radial contrast and satellite-centered outer light

## Result

This exploratory analysis joins the existing ICL radial contrast
(`200–300 kpc` minus the area-weighted `300–500 kpc` target-minus-control
profile) to each cluster's measured satellite-centered `50–75 proper-kpc`
target-minus-azimuthal-control residual profile. It computes the partial
Pearson correlation after linear adjustment for cluster redshift and log CAMIRA
richness. Uncertainty intervals use 20,000 resamples of connected paired
target/control 5-degree sky blocks.

| Filter | ICL extraction | Clusters / blocks | Partial r | 95% block-bootstrap interval |
|---|---|---:|---:|---:|
| `i` | raw | 20 / 10 | −0.140 | [−0.542, 0.630] |
| `i` | outer-plane corrected | 20 / 10 | −0.126 | [−0.541, 0.633] |
| `i2` | raw | 47 / 30 | +0.417 | [−0.154, 0.714] |
| `i2` | outer-plane corrected | 47 / 30 | +0.391 | [−0.195, 0.698] |

One cluster with an otherwise eligible ICL profile lacks a supported satellite
profile in `i2`; it remains missing and is not replaced by zero. The `i` sample
has no such exclusions. The raw and corrected estimates reuse the same images,
masks, clusters, and sky fits; they are sensitivity variants, not independent
replications.

## Interpretation and limits

All four block-bootstrap intervals include zero. The point estimates differ by
filter and are imprecise, so this selected sample does not resolve a relation
between the ICL radial contrast and satellite-centered outer light. It neither
establishes nor rules out a physical connection.

The satellite profile is based on detected CAMIRA-member matches and its
original 40% annular-coverage selection. It is not a complete satellite census
or a clean measurement of satellite stellar envelopes: asymmetric diffuse
cluster light and extended PSF wings can contribute after azimuthal-control
subtraction. The ICL contrast also uses the 300–500 kpc profile as a relative
baseline, which can contain real cluster light. The block intervals do not
propagate measurement errors, source-Poisson variance, exact coadd covariance,
upstream sky transfer, or PSF-wing uncertainty. No p-values or causal claims are
reported. `science_measurement=false` remains in force.

## Reproducibility

Runner: [`analyze_icl_satellite_outer_profile_relation.py`](analyze_icl_satellite_outer_profile_relation.py).

Output directory:
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_outer_profile_relation_2026-11-06_01/`.

The runner reads saved profile archives and compact tables only. It makes no
HSC request and reads no raw image pixels. The summary records source and
output hashes, support, bootstrap settings, and the science gate.
