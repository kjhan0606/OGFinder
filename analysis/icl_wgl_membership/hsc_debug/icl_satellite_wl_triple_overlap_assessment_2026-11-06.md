# ICL–satellite–weak-lensing triple overlap (2026-11-06)

## Scope and method

`analyze_icl_satellite_wl_triple_overlap.py` intersects the already audited i2
ICL–satellite support table with the continuous random-subtracted weak-lensing
per-lens profiles. It verifies the three pinned input hashes and the 27 IDs in
the saved sample-overlap audit, then checks that the two ICL measurements,
redshift, and richness agree lens by lens. No HSC request or network access is
made. The result is exploratory and keeps `science_measurement=false`.

The direct ICL–satellite partial correlation controls redshift and log
richness. In each of seven proper-radius bins, the WL profile is the saved
target-minus-pooled-random `DeltaSigma`. A lens is retained only when both
target and random source-pair counts are positive and all model quantities are
finite. The conditional WL–ICL relation controls satellite ICL, redshift, and
log richness; the WL–satellite relation controls ICL, redshift, and log
richness. Each correlation interval resamples connected 5-degree target and
control sky blocks 5,000 times. The seven radial bins span 0.1–0.2, 0.2–0.4,
0.4–0.7, 0.7–1.0, 1.0–1.5, 1.5–2.0, and 2.0–3.0 proper Mpc.

## Findings

The common sample contains 27 i2 lenses, from 68 satellite-supported systems
and 70 continuous-WL lenses. The direct ICL–satellite partial correlation is
`r=0.067`, with a pointwise 95% connected-block bootstrap interval
`[-0.144, 0.293]`; this selected sample does not resolve a direct association.

For the WL relations, per-bin lens counts are 23, 24, 25, 25, 26, 26, and 27.
All 14 pointwise intervals for WL versus continuous ICL include zero. Two of
the 14 WL-versus-satellite pointwise intervals exclude zero: the cross
component at 0.4–0.7 Mpc (`r=-0.517`, interval `[-0.767, -0.136]`) and the
tangential component at 2.0–3.0 Mpc (`r=0.412`, interval `[0.018, 0.700]`).
These are unadjusted intervals among 29 reported relationships (including the
direct ICL–satellite statistic), with a selected 27-lens sample and estimated
sky blocks. They are post-hoc exploratory patterns, not detections; their
physical or instrumental origin is unresolved.

There are 18–20 recomputed connected blocks per WL fit, with no block larger
than two lenses. These equal-area cells are approximate spatial resampling
units, not validated HSC survey fields. The inherited WL profiles retain
unresolved source photo-z, shear calibration, magnification, masking, PSF,
depth, and covariance systematics. This analysis does not establish a causal
ICL–satellite–halo-mass connection and does not supersede the larger-sample WL
analysis.

## Reproduction and artifacts

Runner:
`analysis/icl_wgl_membership/hsc_debug/analyze_icl_satellite_wl_triple_overlap.py`.

Inputs (SHA-256):

- `supported_pairs.tsv`: `79f7c4c6b18b29ad837e91da4a4359433b09a564cad62efaf60b4cce655685da`
- `per_lens_profiles.tsv`: `d3e07ea875511027739a1a519a7beac765bd00a4ccd3c6b40e23dcc3beb789d8`
- overlap audit `summary.json`: `98d94834b895f1275337f494a801aab7e021048ef502d4f7074635bb26b25e12`

Result directory:
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_wl_triple_overlap_2026-11-06_01/`.
`summary.json` records the full ID list, inputs, estimator, and output hashes;
`metrics.tsv` contains point estimates and block-bootstrap intervals;
`triple_overlap_lenses.tsv` contains the joined per-lens sample. The metrics
and joined-table SHA-256 values are respectively
`afe6fc396119bcb2b6c23af1afda3ae79d8080ecea5b75cc71800c1c6df34b8d` and
`2eedc5d349c4ad2ecc4b95c32c18925c8dab18380ef6e77d3d08bbfe34c94c25`.
