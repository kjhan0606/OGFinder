# Radial selection in the forced CAMIRA-position satellite stack (2026-11-06)

## Scope

This descriptive audit uses the 3,197 CAMIRA positions without a nearby SEP
detection and the existing forced-profile run with a 20% annular-coverage
floor, unmasked neighboring unmatched CAMIRA positions, and at least three
rotated controls. It reads no new HSC data and computes no new fluxes. It
describes which parent positions survive the stack support rules; it does not
test an ICL–satellite relation or a missing-at-random assumption.

## Radial support is strongly selected

At 50–75 kpc from each candidate position, 174 `i` and 292 `i2` candidates
meet the target-only 20% area threshold. Only 25 `i` and 40 `i2` remain after
the rotated-control requirement. The retained position profiles are strongly
shifted toward the outskirts of their host clusters:

| Filter | Parent candidate positions | Parent median clustercentric radius (kpc) | Retained profiles | Retained median radius (kpc) | Retained at <100 kpc | Retained at 300–500 kpc |
|---|---:|---:|---:|---:|---:|---:|
| `i` | 1,292 | 200 | 25 | 388 | 0 | 15 |
| `i2` | 1,905 | 195 | 40 | 339 | 0 | 26 |

Thus 41/65 (63%) retained profiles lie at 300–500 kpc, while no retained
50–75 kpc forced profile lies inside 100 kpc. The threshold passes at least
one profile for only 18/282 candidate-bearing `i` clusters and 32/439
candidate-bearing `i2` clusters. Retention rates vary with host-centric radius;
the parent-position median radius and retained-profile median radius differ
by 188 kpc in `i` and 144 kpc in `i2`.

This pattern means the 20% forced-position sample cannot be treated as a
representative supplement to the detected satellite sample. In particular,
there is no radial overlap for an inverse-probability correction inside
100 kpc, so positivity fails for recovering the full 50–500 kpc parent
distribution. The combined cluster-level ICL–satellite correlation may also
change because the added forced positions disproportionately sample outer
cluster environments.

## Scientific interpretation

The 20% threshold increased the profile count, but it did not solve the
coverage-selection problem. The 65 profiles are light around CAMIRA catalog
positions without nearby SEP detections; individual satellite identities are
unconfirmed. Low target coverage and the control-placement requirement
select a different radial population. Therefore the earlier 20% combined
correlation (`i2` `r=0.257`) must remain a sensitivity result; it cannot be
used as a corrected ICL–satellite relation. Keep the detection-based and
forced-position stacks separate in the primary science analysis. A useful
completeness correction requires validated overlap across cluster-centric
radius and injection-calibrated recovery under the actual masks. The random
point variance export can diagnose local depth but cannot restore this missing
radial support. `science_measurement=false` remains.

## Reproduction and provenance

- Astrafek runner:
  [`audit_icl_satellite_forced_position_coverage_selection.py`](audit_icl_satellite_forced_position_coverage_selection.py)
- Candidate-level radius table:
  `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_forced_position_coverage_selection_2026-11-06_01/coverage_selection_by_radius.tsv`
- Cluster support table:
  `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_forced_position_coverage_selection_2026-11-06_01/cluster_retention_50_75.tsv`
- Summary:
  `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_forced_position_coverage_selection_2026-11-06_01/summary.json`
- Runner SHA-256:
  `1413a0084ff840c17b003a97a90c51197654d28ec2aa8a52b8cd6c1bc68d9fc0`
- Summary SHA-256:
  `6d6650a0ff81d203d11ad012ba82867b9150904af38f511212d1e25931ce9e28`
- No HSC query or download; output `science_measurement=false`.
