# S23B HSC-I footprint metadata audit (2026-11-06)

## Scope and retrieval history

The user approved a metadata-only screen for 1,333 S23B CAMIRA candidate
centers. The final plan contains 11,997 probes (each center plus eight
great-circle points at 1,000 proper kpc), 27 SQL batches, and 14 HSC-I mosaic
and patch-QA metadata columns. It explicitly excludes images, object/source
catalog rows, shape rows, and photo-z. Final plan SHA-256:
`62c28446f553bf45a940ea0126af1521d4ee0c1356c12cc73f42712df482c0c7`.

Two attempts failed before CAS job creation: HTTP 401 at the `hscdata` host and
HTTP 406 for the query using `public.skymap.patch_area`. The corrected query
uses `pdr3_wide.mosaic.band` with the documented `mosaic.areacube` and
`radec_to_coord()` footprint operator. The official [HSC PDR3 stored-function
guide](https://hsc-release.mtk.nao.ac.jp/rsrc/pdr3/stored_functions.html)
documents the coordinate conversion and `areacube` spatial search. The
corrected plan completed all 27 batches and returned 13,575 rows, with no image
or source data. Credentials were not saved.

## Audit findings

Astrafek independently checked the exact plan hash, candidate/probe inventory,
all 27 CSV checksums and schemas, output column contract, HSC-I band labels,
patch identifiers, metadata value ranges, and boolean availability flags.

- 1,048/1,333 candidate centers have at least one HSC-I mosaic metadata match.
- 9,437/11,997 center/boundary probes have a mosaic match.
- 1,043 candidates have matches at the center and all eight 1-Mpc boundary
  probes.
- All 1,043 fully screened candidates have `imeas_exists=true` and
  `idet_bkgd_exists=true` at every probe.
- The export contains 13,575 patch rows. Multiple rows may map one probe to
  overlapping patches.

These results nominate up to 1,043 candidates for subsequent image-level
review. Nine probes do not establish full-disk coverage, usable per-pixel masks,
image WCS integrity, actual cutout availability, or matched-control
availability. Depth and seeing fields are QA metadata, not a validated image
quality cut. No image download was initiated by this audit and
`science_measurement=false` remains.

## Reproducibility

Plan and results:

- Plan: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_s23b_coverage_plan_2026-11-06_04/coverage_query_plan.json`.
- CAS state and result files: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_s23b_coverage_plan_2026-11-06_04/cas_coverage_metadata/`.
- Audit output: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_s23b_coverage_audit_2026-11-06_02/`.

Astrafek tooling:

- Generator: `prepare_s23b_icl_coverage_plan.py` (uses the shared S20A/S23B SQL renderer).
- Approval-gated submitter: `submit_icl_expansion_coverage_metadata.py`.
- Independent audit: `audit_s23b_icl_coverage_metadata.py`.

Integrated report version 31 validates the plan, CAS checkpoint, submission
history, audit summary, all batch files, and per-candidate result table.
`science_measurement=false`.
