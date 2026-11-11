# S23B target coverage in the S16A shape catalog: count-only plan

## Purpose and current state

The existing HSC S19A source-row export supports only a selected 52-lens joint
ICL/WL sample. This offline plan checks whether the public S16A weak-lensing
shape catalog has source coverage around the full set of 1,043 S23B CAMIRA
centers. The exact plan was approved for count-only submission. The immutable
plan file stays pinned to its reviewed hash; query progress and results are
stored separately.

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
- Latest verified execution snapshot (2026-11-11): 31/42 batches complete,
  775 center counts, 100 centers with positive counts, 675 with zero sources,
  and count sum 416,948. Batches 25–30 returned zero for every center.
  Completed rows match the approved center IDs and batch checkpoint hashes;
  batch 31 is submitted and still running.
- Checkpoint (outside Git):
  `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_s23b_s16a_shape_count_plan_2026-11-11_01/results/count_status.json`
- No source rows, images, or photo-z were requested or returned; credentials
  were not saved. `science_measurement=false`.
- Official schema/query examples: [S16A Shape Catalog](https://hsc-release.mtk.nao.ac.jp/doc/index.php/s16a-shape-catalog-pdr2/), [HSC summary-table SQL usage](https://hsc-release.mtk.nao.ac.jp/doc/index.php/database/summary_table_usage-2/), and [HSC FAQ on the PDR1 schema](https://hsc-release.mtk.nao.ac.jp/doc/index.php/faq__pdr3/).

The plan hash and scope are the approval boundary for this count-only query.
