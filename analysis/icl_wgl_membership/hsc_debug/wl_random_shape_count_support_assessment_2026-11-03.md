# Random-center shape-source support: interim audit

## Scope and status

This is an interim audit of the separately approved HSC S19A count-only query for
3,100 random centers (620 matched lenses, five centers per lens), within each
center's matched lens redshift and proper 3 Mpc aperture. At the audit snapshot,
68/124 batches were complete. The audit uses only batches marked complete in
the checkpoint and verifies each CSV SHA-256, exact requested center/lens IDs,
and one non-negative integer count per center. No source rows, coordinates,
photo-z values, or images were requested or returned. This is not a lensing
measurement.

## Interim result

The first 1,700 verified centers cover 340 matched lenses. Their summed
per-center source count is 4,749,572. Of these centers, 600 (35.3%) have zero
shape sources. The complete-lens support distribution is:

| Positive centers among five | Lenses |
|---:|---:|
| 0 | 107 |
| 1 | 8 |
| 2 | 8 |
| 3 | 3 |
| 4 | 3 |
| 5 | 211 |

The first 340 lenses follow the current manifest order and are not a random
subset of the full 620-lens sample. Therefore these fractions are diagnostic
only and must not be presented as final sample support. The result does show
that random-center source availability needs an explicit paired-support and
footprint audit before any random source-row export or random-subtracted stack.
Selecting centers only from positive source counts could alter the angular
selection and source-density distribution; assess that against the lens sample
and HSC mask/depth before defining eligibility.

## Reproduction

The Astrafek audit utility is
[`audit_wl_random_shape_count_support.py`](audit_wl_random_shape_count_support.py).
It is offline, hash-checks completed batches, reports partial/full lens support,
and refuses to overwrite its output. The snapshot JSON is stored outside Git at
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_random_shape_count_support_audit_2026-11-03_01/partial_audit.json`.
It records plan SHA-256
`18aeb1e831d847818b2a671dfe5f75c3122c38637c91b9abbd8fabff5faae625`,
68 verified batches, 1,700 centers, and `science_measurement=false`.

At this initial snapshot the count-only submission was still active. It has now
completed and the full support pattern is recorded below. A separate random
shape-row export still requires its own explicit approval.

## Complete paired target-lens comparison

A complete offline audit verified all 124 batches (3,100 centers, 620 complete
matched lenses) against the already approved and independently hash-verified
800-lens count-only preflight. It did not query HSC or read new source rows. The
preflight has 168 zero-count lenses among all 800; 150 of those are in the
complete 620-lens matched subset. Of these 150 matched zero-count lenses, 131
also have zero sources at all five random centers. Among the 470 matched lenses
with positive lens-position counts, 385 have positive counts at all five random
centers, while 35 have zero counts at all five. Thus many zero-support random
sets are paired with target lenses that themselves have no source-count support.

For the 470 matched lenses with positive lens-position counts, the median ratio
of mean count across five random centers to the lens-position count is 0.878;
the 16th–84th percentile range is 0.419–1.110. This checks only raw candidate
availability; it is not a correction and does not establish that selected
source density, shape weights, or redshift distributions match.

Do not interpret the raw zero-center fraction alone as random-control failure.
The completed count query now permits target and random support to be assessed
jointly, but identical quality, tomography, B-mode, and radial source selection
has not yet been applied to both endpoints. These preflight counts do not give a
final retention fraction or lensing detection.

Across all 3,100 random centers, 1,011 (32.6%) have count zero and the verified
count sum is 9,010,324. The paired snapshot is
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_random_shape_count_support_audit_2026-11-03_01/lens_paired_support_audit_ratio_v25.json`,
SHA-256 `76c7b87f5c099ad2770b56b39aeac0e324ea293e3e644b7495d244fa4bf0c1c4`.

## Next gate: random shape-row export plan

The offline export plan is prepared and passed its immutable-plan check; it has
not been submitted to CAS. It contains 124 batches for the same 3,100 random
centers and predicts 9,010,324 `(random_center_id, object_id)` rows from the
verified count-only results. The planning estimate for uncompressed CSV is
2.15–4.03 GiB (256–480 bytes per row; actual returned size may differ).

The proposed 23 columns include source coordinates, ellipticities, shape weights
and calibration, aperture magnitude, resolution, PSF shape, raw B-mode mask,
and the HSC tomographic-bin label. They omit individual photo-z values/PDFs and
images. A source may recur at different random centers. These are shape source
rows and need a separate explicit user approval; the prior approval covered
count-only rows only. No shape-row query has been submitted.

Plan directory:
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_random_shape_rows_lens_reuse_2026-11-03_01/`.
Manifest SHA-256 is `fc883444d2b23663d5f6034d773b76db948add29f1cda0c23384d17cfed2996a`.
The submission utility's `--check` returned `prepared_not_submitted`, `network_used=false`,
and `query_submitted=false`.
