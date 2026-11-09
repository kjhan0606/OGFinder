# ICL and weak-lensing sample overlap (2026-11-05)

## Result

The local Astrafek audit [`audit_icl_wl_sample_overlap.py`](audit_icl_wl_sample_overlap.py)
reconstructed the exact common-lens set used by the approved random-subtracted
weak-lensing stack from the verified target shape rows, random-center count
exports, and the selected 620-lens plan. The resulting 435 lens IDs reproduce
the stack's stored SHA-256 (`82c0e737ebdbab8027f94267c3f7c5d77e6ed9f3535c6b20b8e8ddf31394eed1`).
The analysis is local; it submitted no HSC query and returned no new data.

| ICL filter | ICL profile support | In 620-lens WL plan | Target shape rows | Positive random support | Random-subtracted common support | 50–75 kpc satellite support in common WL sample |
|---|---:|---:|---:|---:|---:|---:|
| `i` | 71 | 48 | 69 | 38 | 38 | 12 of 26 |
| `i2` | 153 | 130 | 92 | 79 | 70 | 27 of 68 |

The joint ICL-plus-satellite subset is therefore small for filter-specific
weak-lensing inference: only 12 `i` and 27 `i2` clusters have both a saved
50–75 kpc satellite profile and verified cluster/random shape support. These
counts do not justify a confirmatory mass–satellite correlation. The broader
ICL-supported sets (38 `i`, 70 `i2`) can support a separate exploratory
random-subtracted lensing stack, but that estimates lensing for ICL-supported
clusters and does not establish a satellite-structure relation.

The target shape export is also the main common-sample limit: only 35 of the 68
`i2` joint-support clusters and 24 of the 26 `i` clusters have target source
rows. The random count exports show that 31 and 12 respectively have positive
random support; after requiring both sides, 27 and 12 remain. Eight of the 68
`i2` clusters are outside the 620-lens matched-random plan. Do not interpret
missing source rows as a physical absence of lensing signal.

## Next analysis

1. Run an exploratory random-subtracted weak-lensing stack for the 70 `i2` and
   38 `i` ICL-supported common lenses, with the same radial estimator, matched
   random centers, and spatial-block covariance as the existing 435-lens
   product. Keep the filters separate and `science_measurement=false`.
2. Report the 27 `i2` and 12 `i` joint ICL/satellite/WL subset as a small-sample
   diagnostic only. Do not pool filter-specific rows as independent clusters.
3. Before any confirmatory relation, diagnose the target-shape zero-support
   lenses against the saved HSC coadd footprint and review whether the shape
   catalog selection or common-sample definition can be improved. Any new HSC
   extraction remains subject to a separate exact plan and approval.

## Reproducibility

The machine-readable audit and 94-pair eligibility table are under
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_wl_sample_overlap_audit_2026-11-05_03/`.
Summary SHA-256:
`98d94834b895f1275337f494a801aab7e021048ef502d4f7074635bb26b25e12`.
The audit verifies all 32 target-row batch hashes, all 124 random-count batch
hashes, and the reconstructed common-lens ID hash. Its output explicitly keeps
`science_measurement=false`.
