# Cross-release CAMIRA expansion feasibility — 2026-10-21

## Question and method

The current 800-cluster parent sample leaves only 26 `i` and 68 `i2` rows in
the joint ICL 200–300 kpc / satellite 50–75 kpc support cell. I audited whether
the already-local S20A and S23B CAMIRA catalogs, combined under the current
selection, could materially improve this support. The offline audit applies
`0.30 <= z_cl < 0.55`, requires a consistent BCG spectroscopic redshift with
`|z_BCG-z_cl|/(1+z_cl) < 0.01`, excludes candidates within 1 arcsec of an
existing target, then greedily enforces 6 arcmin between all selected
coordinates. Candidates are ranked by richness, with a deterministic
release-prefixed catalog-ID tie break. This is an upper-bound pool audit, not
an image-supported sample.

## Results

The pinned source catalogs contain 7,319 S20A and 12,162 S23B rows. The common
selection yields 869 and 2,288 candidates, respectively. After current-target
matching, the union contains 1,798 rows. There are 607 eligible cross-release
pairs within 1 arcsec; the final 6-arcmin greedy selection produces 1,620
additions (24 S20A and 1,596 S23B), occupying 22 new coarse 5-degree RA/Dec
cells. This gives a catalog-only parent upper bound of 2,420 clusters and 100
occupied coarse cells. These cells are not independent covariance blocks.

If the current selected-sample joint-support yields transferred unchanged,
1,620 additions would contribute an expected 53 `i` and 138 `i2` joint-support
rows, for projected totals of 79 and 206. Those projections remain below all
four illustrative 80%-power, 168-family Bonferroni scenarios:

| Assumed partial r | Design effect | Required i2 rows | Projected total | Fraction |
|---:|---:|---:|---:|---:|
| 0.20 | 1.0 | 536 | 206 | 38% |
| 0.20 | 1.5 | 804 | 206 | 26% |
| 0.293 | 1.0 | 323 | 206 | 64% |
| 0.293 | 1.5 | 485 | 206 | 42% |

The `r=0.293` scenario is a post-selection exploratory estimate, not a
confirmatory effect size. The calculation indicates that combining these two
catalog releases alone is unlikely to provide enough jointly supported rows
for a multiplicity-aware ICL–satellite test at the current support yield.

## Limits and next evidence gate

This result uses catalog positions, redshifts, richness, and BCG spectroscopy
only. It does not establish that the 1,620 candidates have HSC-I coadds,
full-aperture mask/depth/PSF support, matched controls, or independent spatial
covariance. S20A and S23B can assign different properties to the same nearby
cluster; the 6-arcmin selection reduces close duplicates but does not prove
cross-release identity. The current support yield may also fail to transfer to
the expanded selection. No CAS query, images, shapes, or photo-z rows were
requested, and `science_measurement=false` remains in force.

The approved PDR3 nine-point coverage query remains the next catalog-to-image
gate for the S20A expansion candidates. Once its metadata screen and image
support are available, candidate selection must be revisited before any
additional image request. Reaching the illustrative sample sizes likely
requires a wider CAMIRA pool and, more importantly, better joint ICL/satellite
support efficiency; simply adding catalog rows is insufficient.

## Reproduction and provenance

- Astrafek implementation: [`audit_icl_cross_release_expansion.py`](audit_icl_cross_release_expansion.py).
- Tests: [`test_audit_icl_cross_release_expansion.py`](test_audit_icl_cross_release_expansion.py), two passed.
- Output directory: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_cross_release_expansion_2026-10-21_01/`.
- Summary SHA-256: `ac62ab540bae2bede06b04367aac964d1153a2a3b7248ed15882b013ed6b792f`.
- Candidate table SHA-256: `4231d9e3cd74a0108f9bb808d19e188a2b8c719e375fab2f1c9f81f834385625`.
- Analysis script SHA-256: `b75feba5f3305805e6455c5562fdddd4fe6abb46a7e693e0e6b1172d441082ec`.
- Inputs were checked against pinned S20A/S23B catalog, sample-manifest, and sample-size-plan hashes recorded in `summary.json`.
