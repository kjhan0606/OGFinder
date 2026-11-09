# Forced CAMIRA-position stack: coverage-cut sensitivity (2026-11-06)

## Question

The previous forced-position analysis measured 3,197 in-frame CAMIRA member
positions without a SEP detection within 3 arcsec. It retained only nine
positions at 50–75 kpc (four `i`, five `i2`) because each target annulus and at
least three rotated control annuli had to retain 60% of their geometric area.
The target-coverage 95th percentile at 50–75 kpc was only 27.6%. This
diagnostic asks whether the low support is caused by masking other unmatched
CAMIRA positions or by the per-position coverage cut.

All runs use the existing 800 target cutouts and masks; no HSC request,
download, or image product was added. Every output remains exploratory and
keeps `science_measurement=false`.

## Mask-setting comparison at 60% coverage

At fixed 60% coverage, disabling the fallback mask around other unmatched
CAMIRA positions retained the same 20 position–radial-bin measurements.
Across those measurements, the largest absolute target-coverage change was
zero, the largest control-coverage change was 0.00184, and the largest
surface-brightness contrast change was 0.0275 nJy arcsec⁻². Across all 3,197
candidates, the maximum change in the radial coverage vector was 0.0251. The
mask choice therefore does not explain the loss of support at the original
60% floor. The main bottleneck is the high per-position area threshold after
the HSC quality and SEP-source masks are applied.

## Lower-coverage sensitivity

A separate run lowered the minimum target and control annular coverage to
20%, retained the requirement for at least three rotated controls, and left
other unmatched CAMIRA positions unmasked in both target and controls.
Support increased from nine to 65 forced positions at 50–75 kpc: 25 in `i`
and 40 in `i2`, spanning 18 and 32 clusters. Across all radii, 123 of the
3,197 candidates contributed at least one annulus and 263 position-bin
measurements were retained.

The forced-only 50–75 kpc `i2` sample has 21 cluster/filter pairs in 17
connected 5-degree blocks. Its redshift/richness-adjusted partial correlation
is `r=0.147` with a broad spatial-block interval `[-0.485, 0.425]`. The `i`
forced-only sample has only eight pairs and remains below the minimum support
for correlation analysis.

Combining the forced positions with the detected-position profiles gives 29
`i` and 74 `i2` cluster/filter pairs, compared with 26 and 68 at the original
60% cut. In `i2`, the 20% full-sample point estimate is `r=0.257`
(`[0.092, 0.435]` in the unpaired percentile bootstrap). That interval is
conditional on this post-hoc coverage rule and is not multiplicity-adjusted.
The paired comparison uses a shared connected-block partition and shows:

| `i2` comparison | Partial-correlation change | Paired 95% interval |
|---|---:|---:|
| 20% vs 60%, same 68 pairs | +0.003 | [−0.070, 0.126] |
| Add the six new pairs at 20% | −0.025 | [−0.127, 0.079] |
| 20% full 74 pairs vs 60% 68 pairs | −0.019 | [−0.149, 0.147] |

The cut change does not produce a resolved shift on the common sample, and the
six-pair sample expansion is also consistent with no change. The apparent
positive full-sample `i2` interval is not a robust independent detection:
the forced-only estimate is imprecise, support is selected by coverage, and
the threshold was not pre-registered across the broader family of analyses.
No ICL–satellite association is established by this sensitivity run.

## Interpretation and remaining gates

Relaxing the per-object area floor recovers more forced positions, but those
positions are still catalog-selected locations without confirmed SEP galaxy
detections. Missing pixels may be non-random because the SEP masks trace real
nearby sources. Rotated controls match cluster-centric radius and retain three
valid placements, but do not provide exact pixel-by-pixel matched masks. The
candidate mean light therefore cannot yet be assigned to satellite stellar
envelopes. The ICL annulus remains BCG plus diffuse light, and both components
still lack complete PSF-wing, sky-transfer, source-Poisson, and total
covariance calibration. Keep `science_measurement=false`.

The 20% setting is an explicit sensitivity option in the Astrafek forced
profile runner; the default remains 60%. Do not replace the primary science
sample with the 20% result without validating the mask-selection model and
completeness using representative injections.

## Reproduction and provenance

- Forced-profile runner:
  [`measure_forced_satellite_profiles.py`](measure_forced_satellite_profiles.py)
- Paired sensitivity analysis:
  [`analyze_icl_satellite_forced_coverage_sensitivity.py`](analyze_icl_satellite_forced_coverage_sensitivity.py)
- Sensitivity runner SHA-256:
  `901301a8fb492f73bdfa227e7672fb97e9192ab5fe557ef714adf64ba944bb53`
- Sensitivity summary:
  `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_forced_coverage_sensitivity_2026-11-06_02/summary.json`
- Sensitivity summary SHA-256:
  `14367972f95c554c989509f7f1152f2d7fffe09e28c8bcb45dd448cf595ca7a2`
- 60% unmasked-CAMIRA run:
  `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/satellite_forced_outer_profile_unmatched_camira_unmasked_2026-11-06_01/`
- 20% unmasked-CAMIRA run:
  `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/satellite_forced_outer_profile_unmatched_camira_unmasked_cov20_2026-11-06_01/`
- 20% run summary SHA-256:
  `8de336e6d22da5ef1679d60a73ffd15e299f1c6f5303ced6bd6d5a1545e5852c`
- 20% run completed all 721 eligible target frames, 10,000 paired comparison
  replicates, and 5,000 within-run association replicates.

The initial sensitivity summary directory `_01` used a combined maximum across
quantities with different units and is superseded by the separated, unit-aware
summary in `_02`.


A follow-up position-support audit found that the 20%-coverage forced profiles
are selected toward large host-centric radii: none of 65 retained 50–75 kpc
profiles is inside 100 kpc, and 41 are at 300–500 kpc. This precludes
reweighting the forced sample to the full candidate population. Keep the
forced-position result separate from the detected-satellite sample. See
[`radial-selection assessment`](icl_satellite_forced_position_selection_assessment_2026-11-06.md).
