# S23B S16A weak-lensing random-center count plan (2026-11-11)

## Local selection

The previous approved S16A query measured raw catalog counts at the 116 saved
one-to-one ICL image-control positions. Those positions are not the five
random centers per lens required by the weak-lensing subtraction. A local-only
selection therefore reused the pinned S23B candidate-center audit and the
existing deterministic global-uniqueness matcher.

For each of the 116 high-raw-count ICL-paired lenses, candidate positions had
to pass the saved center-flag and metadata-match checks, have no CAMIRA cluster
within 2 proper Mpc, and lie more than 6 proper Mpc from their matched lens.
The 23,200 candidate rows yielded 2,799 positions after these cuts. Five
globally unique positions were available for 114 lenses: 96 `i` and 18 `i2`.
Two `i` lenses were excluded because only four (`S23B_10380`) and three
(`S23B_2473`) candidates passed. The 570 selected positions are unique by both
random-object ID and rounded sky coordinate.

The upstream candidate audit still marks all controls as not ready for science.
The selected positions have not passed image, PSF, or mask QA. This is a
provisional query sample, not a validated weak-lensing random-control sample.
The saved 3,100-center random pool for the separate 620-lens S19A design covers
only four of these 116 candidate lenses; it cannot stand in for the 570
positions selected here.

## Prepared S16A count-only plan

The local planner prepared a 3 proper-Mpc S16A PDR1 cone-count query at each
selected random center, using the matched lens redshift. It has 23 SQL batches
of at most 25 centers, 186,051 SQL bytes, and a 25,113-byte conservative result
bound. If submitted, the coordinates and cone radii are query predicates sent
to HSC CAS; the returned CSV contains only `random_center_id` and the integer
`COUNT(DISTINCT weaklensing_hsm_regauss.object_id)` (570 rows). It requests no
source ID, returned coordinates, shape, photo-z, image, or PSF data.

The plan was validated locally against its pinned selection artifacts and SQL
batch hashes. It has **not** been submitted. The exact plan SHA-256 is
`dfd2fc75650fc08ce4acfe1eef6d0ed476a7bdf8c84b5d083b607f4009678de0`.
Submitting this count-only query requires separate explicit user approval.
Even if the counts are favorable, shape-row transfer remains a separate
approval and requires image/PSF/mask and common-support gates. S16A and S19A
rows must remain separate until release and calibration compatibility is
validated. `science_measurement=false`.

## Reproducibility

- Selector: [`prepare_s23b_s16a_candidate_wl_random_centers.py`](prepare_s23b_s16a_candidate_wl_random_centers.py)
- Count-plan generator: [`prepare_s23b_s16a_wl_random_count_plan.py`](prepare_s23b_s16a_wl_random_count_plan.py)
- Approval-gated submitter: [`submit_s23b_s16a_wl_random_counts.py`](submit_s23b_s16a_wl_random_counts.py)
- Selection output: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_s16a_expansion_wl_random_centers_2026-11-11_02/`
- Selection summary SHA-256: `197dd377f54dd645a0acdaa3059a50c5c6a6a8d87e494faa75f60970981bce33`
- Selected random-center TSV SHA-256: `d63bc851903882e8ad21bc85de522ac2722a49d7b8533f7fc1ea8335fe9cac24`
- Per-lens selection ledger SHA-256: `904744b0d6efd0c7ded8570fe1bc941b9e97367aaaae5282df5448151da5a2f6`
- Prepared plan: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_s23b_s16a_expansion_random_count_plan_2026-11-11_01/`
- Plan SHA-256: `dfd2fc75650fc08ce4acfe1eef6d0ed476a7bdf8c84b5d083b607f4009678de0`
- `submit_s23b_s16a_wl_random_counts.py --plan` passed; it made no network or CAS request.
