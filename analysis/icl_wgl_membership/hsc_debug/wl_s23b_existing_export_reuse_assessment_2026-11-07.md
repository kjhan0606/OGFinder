# S23B Reuse of the Existing S19A Shape Export (2026-11-07)

## Question

Can the already downloaded HSC S19A weak-lensing source rows for 800 lenses
replace a new S19A source-support query or source-row export for the 427 S23B
CAMIRA targets in the member outer-profile analysis?

## Audit performed

The local audit verified all 32 source-row CSV batches against the saved row
plan, run receipts, byte counts, SHA-256 hashes, and count-only expected row
totals. It checked the 22-column schema, uniqueness of each (lens_id,
object_id) pair, consistency of source-level fields for repeated object_id
values, binary raw B-mode values, and valid non-negative tomographic bins.
Every downloaded source coordinate lies within its originating 3 proper-Mpc
lens query cone; the maximum distance-minus-radius value is
−9.26×10⁻⁵ arcsec.

The export contains 3,069,544 lens-source pair rows and 2,929,303 unique
source IDs. Repeated source rows across lenses are consistent; no source-level
metadata conflict was found. The existing source rows were then searched
locally around the exact 427 S23B target centers using 3 proper-Mpc apertures
in a flat cosmology with H0=70 km s⁻¹ Mpc⁻¹ and Ωm=0.3.

Geometric reuse was estimated independently by sampling 4,096 deterministic
Sobol points uniformly in solid angle within each S23B aperture and testing
them against the union of the 800 old S19A query cones. The reported fractions
are estimates of geometric query-cone overlap, not proof that the HSC shape
catalog has uniform coverage in those regions.

## Results

- 43/427 S23B apertures intersect at least one old S19A query cone.
- Only 21/427 apertures contain any source ID from the existing export.
- No S23B aperture is fully contained by one old query cone.
- The estimated old-cone union coverage has median 0, 90th percentile
  0.00713, and maximum 0.604. No aperture reaches 95% coverage.
- No S23B target center is within 60 arcsec of an old S19A lens center.
- Existing source-row counts range from zero for most apertures to 2,995.
  Because the old export samples only its original cones, these counts are
  incomplete wherever the old-query union does not cover the full S23B
  aperture.

The existing lens-centered S19A export therefore cannot provide complete
427-target shape support. Reusing the partial intersections would create a
spatially incomplete sample and cannot replace target-specific queries.

## Extension: add the approved random-center source rows

To test whether previously approved random-center rows improve reuse, a
second local audit verified the existing 3,100-center, 124-batch S19A export
and combined its source IDs with the 800-lens export. All 9,010,324
random-center/source pair rows were checked against their saved batch hashes,
byte and row receipts, exact 23-column schema, requested centers, unique
center/source pairs, and per-batch count-only totals. Repeated `object_id`
values have consistent source fields across both exports. No source rows were
downloaded or sent to a remote service during this audit.

The combined exports contain source rows inside 32/427 S23B apertures: 21 have
rows from the lens export, 16 from the random-center export, and 11 of the 32
have random-center rows but no lens-export rows. The remaining 395 apertures
have no source rows in either saved export. A 4,096-point solid-angle Sobol
estimate of the combined query-cone union has median coverage 0, 90th
percentile 0.310, and maximum 0.994. Only 2 apertures reach 95%, one reaches
99%, and none reaches 99.9%; no aperture is fully contained by any one old
query cone. These are geometric upper bounds, not shape-catalog completeness
or uniform-depth measurements.

Adding the random-center export improves partial coverage but still cannot
support a common 427-target S23B sample. The separate 427-target count-only
plan remains the next data gate; even that count query will not authorize or
replace a later source-row export.

## Combined-export reproducibility

- Runner: `analysis/icl_wgl_membership/hsc_debug/audit_wl_s23b_reuse_existing_s19a_sources.py`
- Results: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_s23b_existing_source_union_2026-11-07_01/`
- Per-target coverage SHA-256: `c6cc5467d031e76d6f867e222f6a1a40b0b76e86b36129f5fd7901853261cede`
- Batch-integrity table SHA-256: `8c0869a8489a348b870196bf51a3e7a86b3631445862a929b0fdf0ecafe4223b`
- Summary SHA-256: `876bed400f93a55c391c89128aa9e9214ed489e664ea0e618ea3df6920f91ed7`
- Random-center export plan SHA-256: `fc883444d2b23663d5f6034d773b76db948add29f1cda0c23384d17cfed2996a`
- No network request, data transfer, commit, or push was made.

The count and selection fields in the per-target table are diagnostics only.
The B-mode mask polarity is not decoded; no lensing-efficiency or n(z)
calibration was applied; photo-z PDFs were not included. Keep
`science_measurement=false`.

## Next data gate

A separate, already prepared HSC PDR3 CAS count-only plan covers exactly the
427 S23B member-profile targets in 18 batches. It returns one namespaced lens
ID and one integer source count per target, with 122,198 bytes of pinned SQL.
Its canonical plan SHA-256 is
b713eebded489c280fd9ce1b5840bc8d7d3e42d029035f8dbd414b72fb31e984.
The plan is in prepared_not_submitted state and has passed its local no-network
plan check. It is distinct from the completed Gaia query and has not been
submitted. A source-row export would be a later, separate approval.

The count query alone will not establish matched random-center support or
produce a calibrated lensing measurement. HSC Y3 tomographic n(z), B-mode
mask polarity, spatial support, and estimator validation remain necessary.
Keep science_measurement=false.

## Reproducibility

- Runner: analysis/icl_wgl_membership/hsc_debug/audit_wl_s23b_reuse_existing_s19a_shape_rows.py
- Results: /gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_s23b_existing_export_overlap_2026-11-07_01/
- Summary SHA-256: 61026a5d82d47303724f576922415c451589f2516f39337adee3abe021c2c54a
- Network use: none
- New HSC rows or images: none
