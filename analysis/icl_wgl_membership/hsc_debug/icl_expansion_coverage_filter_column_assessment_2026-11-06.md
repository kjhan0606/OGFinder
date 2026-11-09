# PDR3 expansion-coverage query field audit — 2026-11-06

## Finding

The two earlier expansion-coverage query generators selected and filtered on
`pdr3_wide.mosaic.band`. The HSC PDR3 stored-function documentation's mosaic
example names the filter column `filter01`, not `band`. The HSC PDR3 database
documentation confirms that `mosaic` is an image-metadata table. This gives a
specific likely cause for the earlier HTTP 406 responses, but the CAS response
body was suppressed, so the cause is not confirmed by the server.

The generators now emit `m.filter01 AS band` and constrain the same metadata
join with `m.filter01 = 'HSC-I'`. The output schema, candidate identities,
probe positions, geometry, and requested metadata columns remain unchanged.
The result label `band` is an output alias, not a database column.

## Prepared plans

The corrected S20A plan covers 2,596 candidates and 23,364 center/1-Mpc
boundary probes in 52 sequential batches. Its plan SHA-256 is
`96ad3b9078beff082d3dce4de1f9bb7608791915ffc7a46a656d9e6cf30e1149` and its
SQL payload is 1,349,895 bytes. It is at
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_expansion_coverage_plan_2026-11-06_02/`.

The corrected S23B plan covers 1,333 candidates and 11,997 probes in 27
batches. At most 47,988 metadata rows are expected from the four-patch-per-
probe ceiling; SQL payload is 751,434 bytes. Its plan SHA-256 is
`ee6d2b11893dab0e9ceb0c0d2a3c80c2dee27b600683839d2be294441578ef05` and it is
at
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_s23b_coverage_plan_2026-11-06_02/`.
S23B IDs remain release-prefixed and are never joined to S20A IDs by numeric
identifier.

Both plans passed the Astrafek submitter's offline integrity and scope checks.
No CAS request or network access was made while preparing or checking these
plans. They request only 14 image-footprint/depth/seeing metadata columns; they
exclude images, object/source rows, shapes, and photo-z. The S23B plan remains
pending approval of its exact hash. A successful metadata screen would still
not establish full-aperture image/mask support or matched-control availability;
those gates precede any separate image request. `science_measurement=false`.

## Reproducibility

- S20A/S23B query generator: `prepare_icl_expansion_coverage_plan.py`, SHA-256
  `564014aaf6c3143b741b385064d368d6bc0c66d8e071d68916e59975a4ee0cee`.
- S23B plan generator: `prepare_s23b_icl_coverage_plan.py`, SHA-256
  `30b7041eea5830c94c6ae04383494fd600e5db5b7403e057fb6d62852777e60d`.
- Approval-gated submitter: `submit_icl_expansion_coverage_metadata.py`,
  SHA-256 `22dd1f486648babcca5d818544e9994d38bce3becc33db2fd29d47e5ec444108`.
- HSC reference: [PDR3 stored-function examples](https://hsc-release.mtk.nao.ac.jp/rsrc/pdr3/stored_functions.html)
  and [PDR3 database documentation](https://hsc-release.mtk.nao.ac.jp/doc/index.php/database__pdr3/).
