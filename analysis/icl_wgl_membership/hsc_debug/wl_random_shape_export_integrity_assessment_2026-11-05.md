# Random-center shape-row export integrity audit (2026-11-05)

## Scope and result

Astrafek's [`audit_wl_random_shape_export_integrity.py`](audit_wl_random_shape_export_integrity.py)
reconciles the approved random-center HSC S19A source-shape export against its
immutable plan, checkpoint, completion summary, and the two saved random-
subtracted stack runs. It makes no HSC/CAS request, reads no credentials, and
does not recompute the scientific estimator.

The audit completed with zero integrity discrepancies:

- All 124/124 batch jobs are marked complete, and all 124 CSV files match the
  checkpoint SHA-256 and byte counts.
- The verified total is 9,010,324 source-pair rows in 2,962,522,918 CSV bytes
  (about 2.76 GiB).
- Both saved stacks reference the same random-row plan SHA-256
  `fc883444d2b23663d5f6034d773b76db948add29f1cda0c23384d17cfed2996a` and
  reproduce the current SHA-256 for every batch.
- The `KeyboardInterrupt` text remains in the top-level checkpoint despite
  `status=complete`, a completed row summary, and 124 complete jobs. The audit
  records this stale field; it does not invalidate the matched file hashes or
  stack provenance. The checkpoint was not edited.

## Scientific scope

The stacks use 620 selected cluster IDs and five matched random centers per
selected cluster. Only 435 lenses have both target and positive random-source
support, so the random-subtracted estimator uses 435 identical lenses and
2,175 random centers. The paired block widths are 5.0 and 2.5 degrees. These
remain exploratory diagnostics: upstream random apertures still lack full
depth, mask, PSF, and sky matching, and the covariance cells do not establish
contiguous survey regions. The audit does not validate those scientific
assumptions; `science_measurement=false` remains in force.

## Reproducibility

- Audit script SHA-256: `1fb7ccf2122a118d0be373ec0fbfb70f486413cd125eb13bbb74c0897d51d60a`.
- Machine summary:
  `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_random_shape_export_integrity_2026-11-05_01/summary.json`
- Summary SHA-256: `b7a371a6ff5e4c8d546cbc88ec81a2289f57c09816fd5874eef135e2cc7da27d`.
- Output retains all per-batch job IDs, row/byte counts, and hashes, plus the
  two stack summary hashes.

No new query, download, or stack was run during this audit.
