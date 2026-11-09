# S20A HSC-I coverage metadata submission audit (2026-11-06)

## Approved scope and preflight

The user-approved S20A metadata-only scope was revalidated on `lageunha`:
2,596 candidates, 23,364 center/1-Mpc-boundary probes, 52 SQL batches, and the
14 approved mosaic/patch-QA metadata columns. Both the plan and the submitter
exclude images, object/source rows, shapes, photo-z, and science measurements.
The exact approval record and scope are in
[`icl_expansion_coverage_query_plan_2026-10-21.md`](icl_expansion_coverage_query_plan_2026-10-21.md).

The first submission of plan SHA-256
`9da220c05f1f24c409a602fcc4ddffc8cbd5bc33f459d335cc8587f1cd235006` was
rejected by HSC CAS with HTTP 406. Inspection found an inner `JOIN` to
`pdr3_wide.mosaic` with no `ON`/`USING` condition. The SQL generator was
corrected to constrain that join with `m.band = 'HSC-I'`, preserving the same
candidate and probe files byte-for-byte and the same output columns and data
scope.

The regenerated plan has SHA-256
`974d4e3db450703481ac2b004d37860806586900097e626c5b456702a9e0346f` and passed
the submitter's offline `--check`. Its first CAS submission also returned HTTP
406. Both checkpoints show `coverage_mosaic_000: pending` with no CAS job ID;
no subsequent batch was submitted and no result CSV exists. The response body
was deliberately suppressed, so the remaining cause of HTTP 406 is not yet
known. Do not retry until the CAS rejection reason is diagnosed and any
possible server-side job state is ruled out.

## Current status

- No coverage metadata results have been received.
- No image, object/source catalog, shape, or photo-z data were requested or
  downloaded by these attempts.
- `science_measurement=false` remains.
- Existing HSC coadds remain unchanged; this query would only prioritize
  candidates for a later, separately authorized image-support evaluation.

## Reproducibility

- Generator: `prepare_icl_expansion_coverage_plan.py`, SHA-256
  `f980d2178270728b8fb7b78949c3ee65ff4c3b40d2acbbd450428162b94c9c9c`.
- Submitter: `submit_icl_expansion_coverage_metadata.py`, SHA-256
  `4307ecd3aeba624bbcfb550fcbaca60918ad7e647fd50f1e9c469098a0602c78`.
- Corrected plan: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_expansion_coverage_plan_2026-11-06_01/coverage_query_plan.json`.
- Corrected plan SHA-256: `974d4e3db450703481ac2b004d37860806586900097e626c5b456702a9e0346f`.
- Corrected candidate centers and probes exactly match the approved v1 files.
