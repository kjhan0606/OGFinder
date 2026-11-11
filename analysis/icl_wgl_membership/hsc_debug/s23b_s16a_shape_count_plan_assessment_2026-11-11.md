# S23B target coverage in the S16A shape catalog: count-only plan

## Purpose and current state

The existing HSC S19A source-row export supports only a selected 52-lens joint
ICL/WL sample. This offline plan checks whether the public S16A weak-lensing
shape catalog has source coverage around the full set of 1,043 S23B CAMIRA
centers. The exact plan was approved for count-only submission and has
completed. The immutable plan remains pinned to its reviewed hash; query
progress and results are stored outside Git.

## Exact proposed query

- Catalog: S16A HSC wide shape catalog in the `pdr1` CAS release, using the
  `s16a_wide` schema.
- Centers: all 1,043 IDs in the pinned S23B ICL target manifest
  (`e487582ba40cd3e58da67704242db4d7219e381fb08f4b33fed47ce4ff268b44`).
- Aperture: 3 proper Mpc at each center's saved cluster redshift and cosmology.
- Submission shape: 42 SQL batches, at most 25 centers per batch.
- Returned values: exactly two columns per center, `lens_id` and
  `n_shape_sources`, where the count is `COUNT(DISTINCT object_id)` after joins
  to the S16A re-Gaussianization shape table and positional measurements.
- Estimated response: 1,043 rows; at most about 31.4 kB of CSV if each 64-bit
  count uses 19 decimal digits. The 42 SQL files total 307,104 bytes.
- Excluded: source identifiers, source positions, ellipticities, weights,
  calibration values, photo-z values, images, and PSF products.

The count is a catalog-availability ceiling before source-redshift, radial,
mask, estimator, and image-matching cuts. It does not establish usable source
support or a lensing measurement. S16A and S19A have different processing and
calibration; their source rows must not be concatenated without a separate
release and calibration audit. Any future shape-row export is a new scope and
requires its own explicit approval.

## Plan integrity and references

- Planner: [`prepare_s23b_s16a_shape_count_plan.py`](prepare_s23b_s16a_shape_count_plan.py)
- Prepared files (outside Git):
  `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_s23b_s16a_shape_count_plan_2026-11-11_01/`
- Plan SHA-256: `c426129b150d440c22fca97d1cf2c722414a7f708132a6267563d77ccda997e0`
- Exact plan approved: yes; SHA-256
  `c426129b150d440c22fca97d1cf2c722414a7f708132a6267563d77ccda997e0`.
- Execution complete (2026-11-11): all 42 batches returned 1,043 center counts.
  The sum of the 3-proper-Mpc aperture counts is 765,197; 177 centers have
  positive counts and 866 have zero. Among the 177 positive centers, the
  median count is 4,622 (range 4–10,138). These are pre-cut catalog counts,
  not usable weak-lensing source counts.
- Independent local QA: all 42 batch CSV hashes matched the checkpoint; every
  batch contained the expected IDs and row count; the combined table has 1,043
  unique IDs exactly matching the pinned target manifest. Nonnegative integer
  values and both output columns were verified. The fixed S19A 52-lens sample
  overlaps positive S16A counts at 21 lenses; 31 have zero.
- Result directory (outside Git):
  `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_s23b_s16a_shape_count_plan_2026-11-11_01/results/`
- Combined count TSV SHA-256: `0945329119d1100b783211fc128cd0a329b449f5acf6648313a1484ac55a7b12`
- Count summary JSON SHA-256: `0cf90929c2a69eda0e62f3a9350359fd430b37bf869b2ead53246129ac51bb8b`
- Checkpoint (outside Git):
  `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_s23b_s16a_shape_count_plan_2026-11-11_01/results/count_status.json`
- No source rows, images, or photo-z were requested or returned. Credentials
  were not saved. `science_measurement=false`.
- Official schema/query examples: [S16A Shape Catalog](https://hsc-release.mtk.nao.ac.jp/doc/index.php/s16a-shape-catalog-pdr2/), [HSC summary-table SQL usage](https://hsc-release.mtk.nao.ac.jp/doc/index.php/database/summary_table_usage-2/), and [HSC FAQ on the PDR1 schema](https://hsc-release.mtk.nao.ac.jp/doc/index.php/faq__pdr3/).

The count-only result is an availability ceiling, not a science-ready WL sample.
Only 177/1,043 centers have any S16A shape-catalog objects in the aperture,
and only 21/52 centers in the existing fixed S19A overlap are nonzero. The
count-only result does not support a 1,043-cluster S23B stack; a selected
177-center subset may have catalog entries but still needs source cuts, matched
controls, and spatial-selection checks. S16A shape-row transfer and any
S16A/S19A merge remain separate approval gates requiring release/calibration
validation.

A later local audit found that among the same 427 targets already counted in
S19A, only 3 are S16A-positive/S19A-zero; 66 are S19A-positive/S16A-zero. A
further 142 S16A-positive centers (126 with raw count ≥1,000) are outside that
427-center S19A audit and have unknown S19A count status. The fixed 52-lens
sample is also imbalanced by ICL band (S16A-positive for 19/25 i lenses but
2/27 i2 lenses). See [the spatial and release-selection audit](s23b_s16a_shape_count_selection_assessment_2026-11-11.md).
