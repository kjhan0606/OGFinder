# S23B member outer-light coverage sensitivity — completed 2026-11-17

## Completed scope and exact reproduction

The checkpointed `_03` run completed all 979 selected frames (226 i, 753 i2)
and its automatic historical audit. The gated sensitivity and figure jobs
also completed. All four final stages' recorded output hashes were checked
again, and the exported figure was visually inspected. Neither original
runner nor follow-up remains live; the terminal products establish completion
of these local steps. This does not complete the broader ICL/WL science goal.

The diagnostic measurement contains 16,503 candidates and 115,521 candidate
by radial-bin records. All 13,269 retained profile rows and all original
candidate diagnostics are exactly identical to the historical full-sample
products. Reconstructing retained target-minus-rotated-control fluxes from
the new raw diagnostic records gives maximum difference **0.0**. These checks
cover the full selected sample, extending the earlier two-frame verification.
The recorded diagnostic audit retains a generic pilot limitation string;
its explicit 979-frame scope and exact comparisons are the authoritative
scope for this completed full run.

Only previously approved images and tables were used. No new query, download,
credential request, GPU or Slurm job was involved. All scientific products
retain `science_measurement=false`.

## Full outer-bin failure census

At 50–75 proper kpc around each photometric member candidate, with 60%
target/control coverage and at least three accepted control profiles:

| Exclusive diagnostic combination | i | i2 |
| --- | ---: | ---: |
| No unique SEP match | 1,281 | 4,269 |
| Target area fails alone | 271 | 1,556 |
| Control-bin area fails alone | 166 | 839 |
| Both target and control-bin area fail | 1,625 | 5,138 |
| Target area fails and too few accepted control centers | 1 | 0 |
| Retained | 160 | 1,197 |
| **Candidate total** | **3,504** | **12,999** |

Joint target/control masked-area failure is the largest matched-candidate
loss category. A shortage of accepted rotation centers is rare in these saved
records. This does not identify individual mask components or rotation
rejection causes, or establish that simply lowering coverage cuts is safe.

## Support and native-sample adjusted association

| Coverage floor | i: clusters / candidates | i2: clusters / candidates |
| ---: | ---: | ---: |
| 0.4 | 133 / 511 | 609 / 3,024 |
| 0.6 | 62 / 160 | 365 / 1,197 |
| 0.8 | 0 / 0 | 60 / 113 |

These counts are conditional retention, not completeness for all satellites.
Native-sample slopes of candidate outer light against relative ICL contrast,
adjusted by the established model, are shown for the 5° phase-zero connected
sky-block layout (5,000 bootstrap draws):

| Filter | Floor | Slope | Marginal 95% spatial-block interval |
| --- | ---: | ---: | ---: |
| i | 0.4 | 0.0490 | [−0.0342, 0.1364] |
| i | 0.6 | 0.0857 | [−0.1000, 0.3455] |
| i | 0.8 | unavailable | no supported candidates |
| i2 | 0.4 | −0.0383 | [−0.1008, 0.0243] |
| i2 | 0.6 | −0.0409 | [−0.1449, 0.0906] |
| i2 | 0.8 | 0.0291 | [−0.1118, 0.2122] |

The 0.6 point estimates reproduce the earlier association. The intervals
use the newly recorded seed, so their endpoints need not equal the old run's
bootstrap endpoints. Native samples and their coverage nuisance covariates
change with the floor; the i2 point-estimate sign change alone cannot diagnose
a photometric threshold effect.

## Fixed support

No i candidate passes 0.8, so an intersection over all three floors is empty
for i. Its common-cluster and common-candidate fits are correctly unavailable.
This is a support limitation, not a zero slope or evidence of no association.

For i2, all three floors share 60 clusters. Their native candidate counts
within that same cluster set are 584, 389 and 113. The strict common-candidate
set is **113 candidates in 60 clusters**. Coverage nuisance covariates are
frozen at the 0.6 reference in fixed modes, giving identical design matrices
for paired slope-change fits. On common individual candidates:

| Floor | Slope | Marginal 95% spatial-block interval |
| ---: | ---: | ---: |
| 0.4 | 0.0242 | [−0.1502, 0.1876] |
| 0.6 | 0.0273 | [−0.1431, 0.1888] |
| 0.8 | 0.0439 | [−0.1324, 0.2113] |

The paired slope changes relative to 0.6 are −0.00315 [−0.01890, 0.00689]
for 0.4 and +0.01657 [−0.04442, 0.09211] for 0.8. Each includes zero. Fixed
clusters but varying member lists also give zero-compatible paired changes:
+0.07599 [−0.01457, 0.18373] and +0.08569 [−0.07969, 0.29559]. These numerical
differences do not establish which mechanism dominates; their uncertainty
and selection restrictions must accompany them.

Across all six sky layouts, **66/108** slope configurations are supported;
42 have insufficient samples, all from the empty i intersection or i at 0.8.
Every supported marginal slope interval includes zero. Of 48 paired-change
configurations, 24 are supported (i2), and every interval includes zero.
No configuration failed rank or bootstrap-validity checks. Layout-dependent
intervals are sensitivity summaries, not 66 independent tests or a proof
that the physical association is absent.

## Decision and next work

The outer-light association remains unresolved. The current data do not
support a positive or negative physical claim, and the broad intervals do
not exclude meaningful associations. Strong coverage selection remains,
especially for i. The strict i2 common sample is small and selected by the
highest coverage threshold. The existing membership, source-Poisson,
extended-PSF, sky transfer, coadd covariance and centering gates remain open.
Do not promote the result to a calibrated science measurement.

Next local steps are to examine the joint target/third-best-control coverage
ceiling and individual mask causes, and to compare common i support between
0.4 and 0.6 as a separately labeled two-floor sensitivity. The latter cannot
use the empty three-floor intersection and must not be selected to maximize
significance. No additional image download is justified by this result alone.
The separate S16A WL random-center query remains approval-gated and unsubmitted.

## Products and provenance

Root:
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_member_bin_diagnostics_full_2026-11-17_03/`.

- `measurement/`: full v2 diagnostic measurement and original profile tables.
- `audit/`: exact historical comparison, failure combinations and threshold support.
- `coverage_followup/analysis/`: cluster profiles, six-layout slopes and paired changes.
- `coverage_followup/figures/`: PNG/PDF sensitivity figure and provenance.

Summary SHA-256 values:

- Measurement: `23b376b26ea7355c94cf1d556f437f32525998561bd3a94b057afe45b44c3b4c`
- Audit: `50a158970521616dc0c9ce29826b1c2e98f15f260e69d30a578fc8479805db28`
- Sensitivity: `f530bf87a59e45205cbc85ff7fec9c7929058e3090ad59fdaa1e9fe2e3404ccf`
- Figures: `faee11db9885b4856a9f6e04d11cf9ebf8e46b918de869a5523a30a93f2493db`

The method and reproduction commands are in
`s23b_member_coverage_sensitivity_methods_2026-11-17.md`; execution history
and bounded-resume instructions are in
`s23b_member_bin_diagnostics_assessment_2026-11-17.md`. Preserve interrupted
`_01`/`_02` artifacts and the successful `_03` checkpoints. No commit/push was
performed as part of this completed local analysis.
