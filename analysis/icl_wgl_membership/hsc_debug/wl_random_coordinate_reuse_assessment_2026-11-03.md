# Weak-lensing random-coordinate reuse assessment — 2026-11-03

## Result

The approved random-coordinate export completed all 32 batches: 780,191 lens–coordinate
associations for 800 lenses, with only `lens_id`, `ra_deg`, and `dec_deg`. Slurm job
448283 completed the offline CAMIRA geometry screen in 31 seconds (peak RSS 669,240
KiB). It found 95,847 associations whose random center has no projected overlap with
either the matched lens's proper 3 Mpc aperture or a CAMIRA cluster's proper 3 Mpc
aperture at the corresponding redshift.

The original screen also required each coordinate to occur exactly once across the
entire 780,191-row table. That stricter rule leaves 68,352 candidates across 544
matched lenses, with at least five candidates for 483 lenses. It discards a clean
association whenever the same sky position appears for another lens, including cases
where that other lens association fails its own overlap test.

An offline sensitivity audit instead applies the no-overlap test to each
lens–coordinate pair, forbids duplicate positions within one lens, and permits the
same clean sky position to be paired with different lenses. This retains 95,847
candidates across 667 lenses; 620 lenses have at least five candidates and 577 have at
least ten. The clean candidate associations contain 82,366 distinct sky positions;
11,514 positions are shared by multiple matched lenses (24,995 associations, or
13,481 assignments beyond one per position). A shared position is used by at most five
lenses.

This comparison supports keeping lens-specific CAMIRA-overlap checks while treating
cross-lens coordinate reuse as a covariance relationship to track, rather than
automatically discarding the association. It does not establish that shared controls
are statistically independent or suitable for a final random-lens stack. Spatial
block covariance must include the shared positions, and the full-aperture footprint,
mask, depth, PSF, and sky must still be checked.

## Prepared count-only feasibility gate

For a balanced preflight, a deterministic SHA-256 ranking selected five centers for
each of the 620 lenses with at least five lens-scoped candidates: 3,100
lens–center associations covering 3,020 distinct positions. Seventy-seven selected
positions are shared across matched lenses. The plan contains 124 SQL batches of 25
centers each. Each returns only
`random_center_id`, `matched_lens_id`, and the integer `n_shape_sources`, counting
distinct HSC S19A shape-catalog object IDs in a 3 proper Mpc aperture before shape,
redshift, B-mode, or radial cuts.

The plan is stored at
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_random_shape_count_plan_lens_reuse_2026-11-03_01/`
with status `prepared_not_submitted`. It has not contacted HSC CAS. It requests no
source rows, shapes, photo-z values, images, or additional coordinates. Because this
is a new HSC aggregate-count export scope, submit it only after explicit user approval.
Even if approved and completed, these counts are a candidate-support check, not a
weak-lensing measurement or a shape-row export approval.

## Reproduction and provenance

- Geometry screen: `run_wl_random_coordinate_geometry.sbatch`; Slurm job 448283.
- Geometry table: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_random_point_geometry_2026-10-19_01/candidate_geometry_audit.tsv`; SHA-256 `b0eb810934639733c2978abc21c8bb79a7e8e0486ec90d97e825cc1fe17f3ac6`.
- Reuse audit: `audit_wl_random_coordinate_reuse.py`; focused tests cover cross-lens reuse and same-lens duplicate rejection.
- Verification: 30 focused tests across the reuse audit, center selector, count-plan builder, and guarded submitter passed; both modified Python modules compile, and `git diff --check` passes.
- Reuse audit output: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_random_coordinate_reuse_audit_2026-11-03_01/reuse_sensitivity.json`; SHA-256 `e464a875f7c97734b797c76e55ce34e6422dd2b587fae59debe5675ce6287bc0`.
- Selected-center manifest: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_random_centers_lens_reuse_2026-11-03_01/random_centers.json`; SHA-256 `728c8e7ae25b961a6d68f4dfff1ac37913535ef20ee074b48b999b6a9bb2f133`.
- Count plan: `random_shape_count_plan.json`; SHA-256 `18aeb1e831d847818b2a671dfe5f75c3122c38637c91b9abbd8fabff5faae625`.
- Count plan state: 3,100 centers, 620 matched lenses, 124 batches, 25 centers per batch; `network_used=false`, `query_submitted=false`, `science_measurement=false`.
