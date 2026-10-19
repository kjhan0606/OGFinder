# 800-cluster ICL–satellite sample retention assessment — 2026-10-19

## Scope and reproducibility

This assessment traces how the existing 800 HSC cluster frames become the
usable ICL–satellite measurements. It uses only the already prepared mask
records, the authoritative corrected full-sample satellite profile run, the
per-annulus satellite stack, and the spatial-block ICL–satellite association.
No CAS query or data download was made. Astrafek's
[`audit_icl_satellite_sample_retention.py`](audit_icl_satellite_sample_retention.py)
verifies each saved input summary/table hash and writes a machine-readable
result under
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_sample_retention_2026-10-19_01/`.
The output summary SHA-256 is
`3bd1dca509aef518db4ae686cfb666ee8922cd3f3415e2f6683286d2f7d367d6`.

The parent sample contains 309 i-band and 491 i2-band frames; all 800 mask
preparations succeeded. This is an image-complete sample, not an 800-object
measurement sample. The mask audit found a unique automatic center candidate
for 258 frames (94 i, 164 i2). With the K2.5 mask, only 129 also have an
unmasked center (51 i, 78 i2). Requiring at least 50% open area in the ICL
200–300 kpc annulus leaves 57 center-eligible frames (24 i, 33 i2); requiring
50% area in all six 50–500 kpc annuli leaves 16 (7 i, 9 i2).

## Satellite profile retention

In the corrected full-profile run `_03`, 9,468 CAMIRA positions lie in the
50–500 kpc candidate range. The target-coverage and matched-control gates
retain 90 satellite profiles from 60 cluster/filter frames:

| Filter | Candidate positions | Target coverage passed | Retained after control gate | Retention |
|---|---:|---:|---:|---:|
| i | 3,689 | 130 | 12 | 0.33% |
| i2 | 5,779 | 358 | 78 | 1.35% |

Coverage failure dominates for the six-annulus requirement. Among the
positions that pass target coverage, 398 fail the minimum of three valid
control profiles (118 i, 280 i2). The preceding `_02` run retained 655 profiles
but mishandled SEP flag 32; `_03` is the corrected run and is the one used
here.

For the preselected outer-satellite annulus, Astrafek also measures each
radial bin independently, without requiring every other bin to pass. The
50–75 kpc target-minus-azimuthal-control stack then contains 42 i-band clusters
and 140 i2-band clusters, with 70 and 308 satellites, respectively. Its
cluster-bootstrap means are `2.51 [−0.30, 5.50]` and
`1.19 [−0.24, 2.49]` nJy/arcsec². Both intervals include zero; nominal pixel
variance does not include correlated coadd noise.

## ICL association and power

Joining the 200–300 kpc ICL measurement to the independent 50–75 kpc satellite
bin yields 26 i-band and 68 i2-band cluster/filter pairs. The i2 partial
Pearson coefficient adjusted for redshift and richness is 0.293. Its 95%
connected spatial-block interval is `[0.042, 0.497]`, while the ordinary-pair
interval is `[-0.007, 0.526]` and the Spearman spatial-block interval is
`[-0.088, 0.301]`. The sample-size calculation estimates only 6.1% power at
that effect under the 30-comparison Bonferroni threshold with 33 i2 spatial
blocks; it estimates 35.9% nominal power. Approximate targets are 92 blocks
for 80% nominal power and 180 for 80% multiplicity-adjusted power. Those are
independent blocks, not additional clusters inside the same footprint.

This is therefore an exploratory hint, not a confirmed ICL–satellite
association. The positive i2 spatial-block interval alone is not decisive
because the ordinary-pair and rank-based intervals include zero, the analysis
is post hoc, and the relevant sky, mask, PSF, and coadd-noise systematics remain
unresolved. `science_measurement=false` remains in all outputs.

## Next analysis gates

Keep the radial-bin-specific estimator for the 50–75 kpc primary satellite
measurement; all-bin completeness unnecessarily discards valid outer-annulus
data. Before changing the mask or control criteria, inspect why the 398
target-qualified positions lack three controls and retain the current
selection as the reference. The 800-frame sample cannot provide the planned
confirmatory power for a moderate effect under the current multiplicity and
spatial-block assumptions. A larger independent cluster footprint or a
pre-registered lower-comparison analysis will be needed; further cluster data
requests require a separate exact plan and approval.

The weak-lensing stack remains exploratory until B-mode polarity is grounded
in catalog semantics and random-lens subtraction plus spatial-block covariance
are available. No weak-lensing data were altered for this assessment.
