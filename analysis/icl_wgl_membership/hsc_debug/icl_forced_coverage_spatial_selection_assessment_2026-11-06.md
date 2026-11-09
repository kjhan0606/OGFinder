# 2026-11-06: spatial overlap of forced-profile coverage selection

## Question

Does the outcome-related 20/20 forced-profile coverage selection occur only
between separated sky regions, or is selected/excluded support present within
the same connected target/control sky blocks?

## Method

The analysis joins all 224 rows from the hash-pinned 20/20 comparison to the
existing target and rotated-control download manifests. Per filter, it builds
connected paired-sky blocks by linking pairs that share either endpoint's
5-degree sky cell. For blocks containing both selected and excluded pairs, it
computes the within-block selected-minus-excluded mean ICL 200–300 kpc residual.
The reported interval is a percentile bootstrap of these block-level contrasts,
with equal weight per mixed block. It is descriptive and does not adjust for
redshift or richness.

No image pixels or HSC endpoints were accessed.

## Results

| Filter | Connected blocks | Mixed blocks | Selected / excluded in mixed blocks | Selected / excluded in single-status blocks | Equal-block mean ICL difference, selected minus excluded (nJy arcsec⁻²) | 95% block-bootstrap interval |
|---|---:|---:|---:|---:|---:|---:|
| i | 18 | 9 | 40 / 19 | 10 / 2 | +1.80 | [−3.26, +6.66] |
| i2 | 40 | 21 | 77 / 40 | 33 / 3 | +3.73 | [+0.55, +6.93] |

Thus the higher ICL residual among selected `i2` pairs persists within mixed
spatial blocks in this equal-block descriptive summary. The `i` comparison is
too imprecise to resolve its sign. Many selected pairs still occupy
single-status blocks, particularly in `i2`; block overlap therefore does not
make a parent-sample transport correction valid.

## Interpretation and limits

- Coverage support is outcome-related within the observed spatial overlap,
  especially in `i2`; between-region composition alone does not account for
  the descriptive difference.
- The 20% threshold is post-hoc. The bootstrap resamples only 9 `i` or 21 `i2`
  mixed blocks and does not propagate ICL measurement uncertainty.
- Connected 5-degree blocks approximate shared target/control dependence.
  Block contrasts are unadjusted for redshift and richness and are not a
  randomized or causal comparison.
- Forced CAMIRA positions without unique SEP detections are not confirmed
  isolated galaxies. No selection correction or science promotion follows;
  `science_measurement=false` remains.

## Reproduction

Runner: [`audit_icl_forced_coverage_selection_blocks.py`](audit_icl_forced_coverage_selection_blocks.py), SHA-256 `35888cca6115619e77b23a465ce769aa61e6af495bdd4c6bacd45382a73d3ae8`.

Output directory: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_forced_coverage_selection_blocks_2026-11-06_02/`.
Summary SHA-256 `3eb63f02507fa83cf1eea377ebc3574bd2b66e86f204bb76face995dcf5f5742`;
block table SHA-256 `fc0df26e806db6a1c16adab0ec2ee2a3faae455d493062578c829d0efdc24a98`;
mixed-block contrast table SHA-256 `0a50f8658552739eb318837d4d58809e0894caba42db525e58939d7ff180609b`.
