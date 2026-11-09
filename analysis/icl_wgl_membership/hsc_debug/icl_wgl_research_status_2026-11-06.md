# HSC ICL, satellite, and weak-lensing research status

## Current sample and ICL stack

Astrafek holds 800 HSC-I cluster cutouts and 798 matched control cutouts.
The 200–300 proper-kpc ICL annulus has target-minus-control profiles for 224
cluster/filter pairs. The supported means are 3.23 nJy arcsec⁻² in `i`
(71 pairs, 95% connected-block interval 2.15–4.23) and 3.79 nJy arcsec⁻² in
`i2` (153 pairs, 2.68–4.80). The positive contrasts describe the existing
coverage-selected sample. The current photometric measurement has not
separated BCG light from diffuse ICL and therefore is not a pure ICL flux.
Mask, sky, PSF, completeness, source-Poisson, and total covariance corrections
remain incomplete.

## ICL and satellite outer light

The saved satellite 50–75 kpc profile overlaps 94 of 224 ICL-supported pairs,
with 26 `i` and 68 `i2` pairs. After redshift and richness adjustment, the raw
partial correlation is 0.079 in `i` (95% interval −0.183–0.477) and 0.293 in
`i2` (0.039–0.500). The `i2` estimate weakens to 0.238 (−0.023–0.446) after
finite-lag pixel-noise adjustment and satellite-support weighting. Adding the
delivered variance ratio changes the estimate to 0.278 (−0.018–0.495). Adding
the variance ratio and coarse depth proxy changes the estimate to 0.201
(−0.054–0.410). Every adjusted interval includes zero. The association remains
conditional on a selected profile sample and an unverified missing-at-random
assumption.

The existing satellite profiles match CAMIRA member positions to SEP
detections within 1 arcsec, then measure light at 5–75 kpc from each satellite
center. Up to eight azimuthally rotated controls at the same cluster-centric
radius estimate the local sky and diffuse-light baseline. The reported
satellite contrast subtracts the mean control profile from the satellite
profile. The subtraction reduces the average host ICL contribution but does
not remove asymmetric cluster substructure, residual sources, or PSF wings.
Astrafek has not separated a calibrated satellite stellar envelope from those
residuals, so the measured profile does not yet identify tidal truncation or a
satellite outer-halo component. Satellite and cluster support also change
between radial bins, so the radial profile does not follow a fixed satellite
sample.

## Weak lensing

The random-subtracted HSC S19A stack has common target and random support for
435 lenses. The tangential `DeltaSigma` profile declines from
`16.01 ± 1.52` to `1.68 ± 0.17` × 10¹³ M☉ Mpc⁻² between 0.1–0.2 and
2–3 proper Mpc under the saved 5° block jackknife. The seven-bin cross-component
null statistic is χ²=10.59 for seven degrees of freedom, with nominal p=0.157.
Alternative block-width tests retain the point estimates but change the error
scale. The covariance blocks are approximate and do not represent validated
contiguous HSC fields.

A random-only placebo reran the same calibrated estimator over the 2,175 matched
random centers for those lenses. The full tangential profile is consistent with
zero (5-degree delete-one-block χ²=6.90, rank 7, approximate p=0.439); the cross
profile is also consistent with zero (p=0.773). This finds no large residual in
the random-center stack, but it does not validate the cluster signal or the
regional covariance. Assessment:
[`random-only placebo`](wl_random_only_placebo_assessment_2026-11-06.md).

The continuous ICL–lensing analysis uses 70 common-support `i2` lenses. The
outer 2–3 Mpc tangential slope has nominal p=0.021, but the seven-bin
multiple-comparison-adjusted minimum is 0.149. The radial slopes change sign,
and all cross-component slopes have nominal p≥0.20. The median-split `i` sample
has only ten matched pairs and a strong cross-component null failure. The
triple overlap of ICL, satellite, and lensing support contains 27 `i2` lenses.
The direct ICL–satellite partial correlation is 0.067 (95% interval
−0.144–0.293). The small common sample and unresolved lensing calibration
systematics do not establish an ICL–halo-mass relation.

## Remaining noise and sample gates

The approved frame-metadata checkpoint contains 67,071 rows for 800 target and
798 control centers. Frame exposure, gain, zeropoint, sky-RMS, and MJD fields
have zero missing values. Frame membership does not supply pixel-to-frame
weights or an amplifier map, so the rows cannot determine local coadd variance
or source-Poisson correction.

An additional fixed-sample sensitivity joined these audited frame summaries to
the 94 ICL/satellite pairs and added five target/control metadata offsets,
one at a time, to the existing redshift- and linear-richness-adjusted
correlation. Across the ten filter-by-metadata comparisons, all max-T adjusted
p-values are at least 0.895 and all simultaneous intervals include zero. In
`i2`, adjusted r remains 0.286--0.309 around the baseline 0.293. This does not
show that coarse metadata explain the association, nor does it constrain
subpixel variance, Poisson noise, sky transfer, PSF wings, or selection. The
profiles remain uncorrected and `science_measurement=false`. See the
[`frame-metadata sensitivity assessment`](icl_frame_metadata_satellite_relation_assessment_2026-11-06.md)
and result summary
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_frame_metadata_satellite_relation_2026-11-06_03/summary.json`.

A separate post-hoc sensitivity adds quadratic redshift and richness terms and
their interaction to the same fixed sample. The `i2` estimate ranges from
0.294 to 0.305 around 0.293, and all eight max-T adjusted model-change tests
remain unresolved. The result is insensitive to these smooth model-form
changes, but selection, measurement covariance, and the small common support
still prevent a science claim. See
[`nonlinear adjustment assessment`](icl_satellite_nonlinear_adjustment_assessment_2026-11-06.md)
and the consolidated report
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_wl_current_evidence_2026-11-06_08/icl_wl_current_evidence_report.md`.

The previously saved random-point sample contains 200 points per target in a
30′ cone. Only 109 points fall in the target 200–300 kpc ICL annulus across all
800 targets, and 98 targets have at least one point. Only three points fall
within 50–75 kpc. The saved random-point sample has no matched-control entries.
The sample cannot measure paired local noise in either target or control
annuli.

Astrafek prepared an approval-gated PDR3 random-point export for the 200–300
kpc ICL annulus at all 1,598 target/control centers. The plan requests local
`i_pix_variance`, sky statistics, input counts, PSF moments, and mask flags.
The plan estimates 46,424 rows and approximately 10 MB in 16 SQL batches.
Plan SHA-256 is
`eeef98ef2210337f714e36d8a4afcd352fe312a536c90e8e236ff9f843f53e52`.
The plan is at
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_random_pixel_variance_plan_2026-11-06_03/plan.json`.
The user approved this exact plan on 2026-11-06. Submission attempts reached
the hidden-credential prompt but the execution session exited before recording
a job ID; the checkpoint remains `pending` with zero jobs,
`query_submitted=false`, and `network_used=false`. No HSC request or data
transfer occurred. The approval remains valid for this pinned scope. A
random-point variance export will improve the local noise audit but will not
replace a complete coadd weight map or correlated-noise and source-Poisson
modeling.

The main science flag remains `science_measurement=false`. The next data step
is to submit the already-approved pinned random-point plan once secure hidden
credential input works, audit target/control annular coverage and variance, and
rerun the ICL contrast sensitivity. The
satellite-envelope interpretation then requires a model that separates
satellite-centered stellar light from the local diffuse ICL. Weak-lensing
inference requires validated survey-region covariance and propagation of
photo-z, shear, boost, reduced-shear, magnification, depth, mask, and PSF
uncertainties.

## Reproducibility

- Science plan: [`science_research_plan.md`](../science_research_plan.md).
- ICL/satellite IPW and finite-lag sensitivity:
  [`icl_satellite_ipw_pixelnoise_joint_assessment_2026-11-06.md`](icl_satellite_ipw_pixelnoise_joint_assessment_2026-11-06.md).
- Continuous ICL–lensing relation:
  [`wl_continuous_icl_relation_assessment_2026-11-06.md`](wl_continuous_icl_relation_assessment_2026-11-06.md).
- ICL/satellite/lensing triple overlap:
  [`icl_satellite_wl_triple_overlap_assessment_2026-11-06.md`](icl_satellite_wl_triple_overlap_assessment_2026-11-06.md).
- Existing random-point support audit:
  [`icl_random_variance_support_assessment_2026-11-06.md`](icl_random_variance_support_assessment_2026-11-06.md).
- HSC random-point field definitions:
  [PDR3 random-point documentation](https://hsc-release.mtk.nao.ac.jp/doc/index.php/random-points__pdr3/).

## Latest satellite-support and random-point variance status

An offline forced-position profile run tested 3,197 CAMIRA member positions
without a nearby SEP detection in the existing 800 target cutouts. Only 14
positions passed any annulus and 9 passed 50–75 kpc, so this method does not
recover the detection-selected satellite sample. Forced-only ICL association
support is insufficient. Combining the few measurements with the existing
detected profiles leaves the exploratory i2 partial correlation near its
previous value (`r=0.287`, 95% approximate block interval 0.037–0.490); this
does not remove selection or noise systematics. Full result and hashes:
[`icl_satellite_forced_position_assessment_2026-11-06.md`](icl_satellite_forced_position_assessment_2026-11-06.md).

The approved paired PDR3 random-point variance export has completed. All 16 SQL batches returned 40,048 rows (11,138,863 bytes) for the 200–300 kpc annulus at 800 targets and 798 controls. Points are present at 1,587/1,598 centers; 787/798 matched pairs have valid `i_pix_variance` summaries at both ends. The post-download auditor checked row schemas, hashes, center/role mapping, and annulus boundaries. Median center-level random-point count is 22 (q16–q84: 11–39); median per-center median `i_pix_variance` is 0.003649 catalog units. This sparse local statistic is not a coadd variance map and does not measure resampling or source-Poisson covariance. It has not yet been used to adjust the ICL–satellite relation. `science_measurement=false`. See `icl_random_pixel_variance_export_assessment_2026-11-06.md` and audit summary `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_random_pixel_variance_audit_2026-11-06_01/summary.json`.

## Random-subtracted lensing n(z) sensitivity

Using the same 435 common lenses and 2,175 matched random centers, Astrafek
re-estimated the source-bin-1-excluded HSC Y3 stack for the fiducial SACC
curves, published central source-bin 3/4 shifts, and four marginal 68% shift
corners. The fiducial pair-level stack reproduces the saved jackknife profile
to relative precision `8.0e-16`. For fixed `c=3`, the effective NFW mass is
`1.750e14 M☉` fiducially and `1.706e14 M☉` under the central shift (−2.60%);
the four endpoint corners vary by about ±1% around that central-shift mass.
Changing fixed concentration to `c=4` or `5` gives poor conditional profile
fits, so model/concentration assumptions dominate this limited n(z) grid.
The saved covariance was held fixed and the curves were translated rather
than using posterior spline samples. This is an exploratory calibration
sensitivity, not a mass measurement; `science_measurement=false` remains.
Details: [`wl_random_subtracted_nz_mass_sensitivity_assessment_2026-11-06.md`](wl_random_subtracted_nz_mass_sensitivity_assessment_2026-11-06.md).

## Supported-sample long-lag covariance

For the 224 pairs supporting the 200–300 kpc ICL contrast (71 `i`, 153 `i2`),
the local residual covariance kernel was extended from ±4 to ±8 pixels using
32 deterministic patches per image and paired 5-degree sky-block bootstrap.
The approximate square-kernel target−control median factor is 0.025
(95% interval −0.002–0.047) in `i` and 0.059 (0.026–0.090) in `i2` at 4 px;
at 8 px it is 0.081 (−0.001–0.132) and 0.077 (0.030–0.129), respectively.
The cumulative factors are non-monotonic because some lag shells contribute
negative covariance, so ±4 is not shown to be converged. The `i2` filter-wise
interval at 8 px excludes zero, while the `i` interval includes zero; interpret
both as exploratory and not multiplicity-adjusted. Spatially varying and
longer-range covariance, source Poisson, pixel coadd weights, sky transfer,
and PSF wings remain unmodeled. No science gate changes. See
[`long-lag covariance assessment`](icl_supported_long_lag_noise_assessment_2026-11-06.md).

The ±8 kernels were then propagated through the actual annulus masks for the
same 224 pairs. Annular mean and coverage matched the saved profiles for all
pairs. The independent-frame finite-lag pixel-noise sigma is 0.0793 nJy
arcsec⁻² in `i` and 0.0745 in `i2`, only 0.8% and 1.3% above the previous 9×9
estimate; these are 7.6% and 7.0% of the saved spatial-block interval
half-widths. This makes finite-lag pixel noise a small component of the
observed stack scatter under the stationarity and independent-frame
assumptions, but it is not the total error: source Poisson, sky-transfer,
PSF, >8-pixel covariance and field covariance remain. No ICL signal correction
was applied. See
[`17×17 aperture propagation`](icl_long_lag_aperture_covariance_assessment_2026-11-06.md).

## Target-only ICL–satellite injection cross-talk

An offline join of the observed 94-pair association sample to audited
PSF-convolved injections found valid target-side transfer responses for 58 pairs
(14 `i`, 44 `i2`). The median ICL-annulus response to an injected satellite
profile was 0.357 in `i` and 0.490 in `i2`; the corresponding diffuse-injection
responses were 0.211 and 0.224. The satellite-profile response distributions
were broad and crossed zero. Target-only injection measures the incremental
signal response of the target-minus-control ICL contrast when matched random
controls remain fixed.

A median-ratio leakage stress test reversed the observed partial correlation,
but the calculation is not a correction and does not establish a physical
anti-correlation. Injection support is selected, 36/94 observed pairs lack a
valid response, and the transfer matrix omits source Poisson, coadd weights,
extended PSF wings, and upstream sky transfer. A pairwise 2×2 matrix inversion
on all 58 available injection pairs changed the i2 partial correlation from
0.399 to −0.381. The conditional block interval was [−0.971, −0.281], but 21/44
diffuse amplitudes and 18/44 satellite amplitudes were negative, with maximum
absolute values of 1,840 and 2,211 nJy arcsec⁻². The bootstrap held the transfer
matrices fixed and omitted response uncertainty, so the interval is not a
physical anti-correlation. The extreme amplitudes show that the current
zero-intercept linear inversion is unstable for the selected sample. The next
calibration must estimate response uncertainty, expand representative support,
and model measurement covariance before interpreting a de-mixed relation. See
the [pairwise inversion assessment](icl_satellite_transfer_matrix_deconvolution_assessment_2026-11-06.md) and
the
[`target-only cross-talk assessment`](icl_satellite_injection_cross_talk_assessment_2026-11-06.md).
The science gate remains `science_measurement=false`.

The same 137 PSF-supported target frames were also processed at satellite
injection amplitudes 0.5, 1.0, and 2.0 nJy arcsec⁻². The 1.0-amplitude response
reproduces the cross-talk table exactly. Redetection support leaves a common
54-pair sample across all amplitudes (12 `i`, 42 `i2`). The pairwise i2 matrix
inversion changes from partial `r=−0.922` at amplitude 0.5 to `r=−0.075` at
amplitude 2.0. Maximum matrix condition numbers reach 2,685 and 962, and
inferred amplitudes include large negative and extreme positive values. The
bootstrap holds measured matrices fixed. The negative intervals are not a
physical anti-correlation because pairwise inversion is unstable and response
uncertainty is omitted. See the
[`amplitude-grid matrix assessment`](icl_satellite_amplitude_grid_transfer_matrix_assessment_2026-11-06.md).
`science_measurement=false` remains.

## Random-subtracted lensing and connected HSC-region influence

An offline leave-one-region-out analysis reconstructed the existing paired
lensing stack from 2,075,062 cluster-lens and 8,905,892 random-center source
rows. The reconstructed seven-bin tangential and cross profiles match the
saved profile exactly. The HSC patch-QA tract-center graph contains four
connected components, three of which overlap the 435-lens common-support stack.
Removing the dominant component removes 261 lenses and shifts ΔΣ by −28.6% at
0.1–0.2 Mpc and −26.8% at 0.4–0.7 Mpc. The shifts reach 3.29 times the
corresponding 5-degree jackknife diagonal error. The graph uses approximate
nearest tract-center assignment and has only three occupied components, so
the analysis measures regional influence rather than a field covariance. The
positive declining profile remains exploratory and `science_measurement=false`.
See [`connected HSC-region influence assessment`](wl_random_tract_region_influence_assessment_2026-11-06.md),
runner `analyze_wl_random_tract_region_influence.py`, and output
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_random_tract_region_influence_2026-11-06_01/summary.json`.

The tract-center connected-component partition remains unchanged at adjacency
thresholds 1.6°, 2.0°, 2.4°, and 3.0°. The same three occupied components
contain 261, 141, and 33 common-support lenses. The threshold audit strengthens
the influence result against that graph choice but does not validate tract
polygons or turn three regions into a covariance estimate. See
`audit_wl_random_tract_partition_sensitivity.py` and
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_random_tract_partition_sensitivity_2026-11-06_01/summary.json`.

## S23B sample-expansion feasibility

The saved 1,333-object S23B candidate pool adds almost no new connected HSC
region support to the 800-cluster discovery sample. A nearest-tract-center
audit places 1,330 candidates in the three already occupied components and
only three candidates in the fourth component. The candidate pool remains an
unmeasured catalog upper bound. Downloading the full pool could add exploratory
measurements but would not yield a useful independent field validation sample.
The audit used the saved candidate-center ledger and patch-QA catalogue only.
See `analysis/icl_wgl_membership/hsc_debug/audit_s23b_candidate_tract_region_support.py`
and `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_s23b_tract_region_support_2026-11-06_01/summary.json`.

## ICL-strength split satellite profiles

An existing-data re-stack compared the 94 supported ICL/satellite pairs in
seven radial bins, split at the within-filter median ICL signal and adjusted
for redshift and richness. All 14 spatial-CR1 intervals include zero; the
10,000-replicate shared-block Webb bootstrap gives Holm p=1.0 for all profiles
and minimum max-T p=0.520. Positive i2 point estimates do not survive the
spatially clustered uncertainty, and the i profile is noisy. This selected,
post-hoc comparison does not establish the ICL–satellite relation.
See `icl_strength_split_satellite_profiles_assessment_2026-11-06.md` and
`analyze_icl_strength_split_satellite_profiles.py`. The result remains
`science_measurement=false`; no HSC request was made.

## Joint continuous and median-split radial ICL–satellite relation

Added a 28-test family spanning continuous ICL slopes and median-split
contrasts across both filters and seven satellite radii, with redshift/log-
richness adjustment and shared spatial-block Webb bootstrap weights. The
smallest raw p is 0.0369 (continuous i2 at 50–75 kpc), but max-T p=0.770; the
minimum max-T p is 0.190 and all Holm p-values are 1.0. Three continuous i
inner-bin pointwise CR1 intervals lie below zero, but the restricted-bootstrap
and familywise tests do not support a detection. The 94-pair support subset is
selected; this does not establish a physical relation. See
`icl_satellite_radial_relations_assessment_2026-11-06.md` and
`analyze_icl_satellite_radial_relations.py`. No new HSC query; the gate remains
`science_measurement=false`.

## Radial ICL profile contrast against the outer annuli (2026-11-06)

Reconstructed the saved 2.5-Kron, 40%-coverage complete-case profile across
50–500 proper kpc. The common six-annulus support is 20 `i` pairs in 10
connected 5-degree blocks and 48 `i2` pairs in 30 blocks. The 200–300 kpc
target-minus-control surface brightness exceeds the area-weighted 300–500 kpc
profile by 2.780 nJy arcsec⁻² in raw `i` ([0.825, 5.115]) and 2.494 in
outer-plane-corrected `i` ([0.640, 4.774]). The corresponding `i2` contrasts
are 2.359 ([−0.110, 4.219]) and 2.159 ([−0.296, 4.047]). Intervals are
20,000-resample percentile intervals from connected paired-sky block
resampling. The outer annuli can contain cluster light and are not blank sky.
The cut is exploratory, the pairs are coverage-selected, and source-Poisson,
exact coadd covariance, upstream sky transfer, and extended PSF wings remain
unresolved. No pure-ICL or ICL–satellite claim follows. Assessment:
`icl_radial_outer_baseline_contrast_assessment_2026-11-06.md`.

The consolidated report was refreshed to version 32 with 933 hash-verified
inputs. Report SHA-256 `6a76bec4c8549d7cd887a5d2b115a3559d1ec8e03e10b719391813ebde83d06a`;
summary SHA-256 `ff8f03cc0593ff3072d0d2242eff20b5cbcfa4275e65c5dabafa44cb1a46a87c`.
Independent input-hash audit found zero mismatches. `science_measurement=false`.


## Weak-lensing boost sensitivity (2026-11-06)

A paired nonlinear jackknife sensitivity calculation multiplies the saved
random-subtracted DeltaSigma profile by the exact S19A estimator-weight boost
for each radius. Across 0.1–3 Mpc, the candidate profile is 12.5–16.3% higher;
2.5- and 5-degree central values agree. This does not establish that the
excess weight is cluster-member dilution, and the correction remains
unvalidated. The saved profile is unchanged; full-aperture depth, mask, PSF,
sky, source selection, magnification, source-redshift uncertainty, reduced
shear, source Poisson noise, and field covariance remain open. Keep
`boost_correction_applied=false` and `science_measurement=false`. See
[`wl_boost_profile_sensitivity_assessment_2026-11-06.md`](wl_boost_profile_sensitivity_assessment_2026-11-06.md)
and result directory
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_source_density_boost_spatial_2026-11-05_06/`.


## Forced-position coverage sensitivity (2026-11-06)

At the original 60% area floor, disabling fallback masks around other
unmatched CAMIRA members leaves the nine 50–75 kpc forced positions
unchanged. At 20% coverage with at least three rotated controls, support rises
to 65 positions (25 `i`, 40 `i2`), and combined support rises to 74 `i2`
cluster/filter pairs. The forced-only `i2` relation is weakly constrained
(`r=0.147`, connected-block 95% interval −0.485–0.425). The combined
`i2` estimate is `r=0.257`, but paired 20%-versus-60% changes include zero
for both the common sample (+0.003, −0.070–0.126) and the six-pair expansion
(−0.025, −0.127–0.079). The lower threshold is post-hoc; mask missingness is
not random by construction, and unmatched CAMIRA positions are not confirmed
satellites. No science association is established. Assessment:
[`icl_satellite_forced_coverage_sensitivity_assessment_2026-11-06.md`](icl_satellite_forced_coverage_sensitivity_assessment_2026-11-06.md).
Results: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_forced_coverage_sensitivity_2026-11-06_02/`.
The separately prepared HSC random-point variance plan has exact-plan approval,
but remains unsubmitted while secure hidden-credential input is unavailable.


## Forced-position sample radial-selection audit

The retained 20%-coverage 50–75 kpc forced profiles are concentrated at
large clustercentric radii: zero lie within 100 kpc, while 41/65 lie at
300–500 kpc. Candidate-radius medians are 200 kpc in `i` and 195 kpc in `i2`;
retained medians are 388 and 339 kpc. Only 18/282 candidate-bearing `i`
clusters and 32/439 `i2` clusters contribute any profile. This lacks radial
positivity for correction to the parent sample; keep forced profiles separate
and do not interpret the combined 20%-cut correlation as a corrected relation.
Assessment: [`icl_satellite_forced_position_selection_assessment_2026-11-06.md`](icl_satellite_forced_position_selection_assessment_2026-11-06.md).

## Spatial out-of-block prediction check (2026-11-06)

Leave-one-connected-block-out validation compared satellite 50–75 kpc excess
predictions from redshift and richness with predictions that also use ICL
200–300 kpc contrast. In `i`, adding ICL increased out-of-block MSE by 6.0%;
in `i2`, it reduced MSE by 3.5%. The paired spatial-block bootstrap intervals
for mean squared-error gain include zero in both filters. The selected 94-pair
sample therefore provides no clear evidence that ICL improves prediction in
unseen blocks. This predictive check does not test causality or remove
photometric systematics. Reproduction and details:
[`icl_satellite_spatial_cross_validation_assessment_2026-11-06.md`](icl_satellite_spatial_cross_validation_assessment_2026-11-06.md).
Runner hash `a068c2ac0330a012964596e7bafc11a3c3d2be7aca293866a0e90b4024b8b510`;
summary hash `7059c7b196467afb581f0edf99ce69a65b9e132979f7a3284be8182c6cc1ceb7`.

## Fixed-support radial ICL–satellite test (2026-11-06)

A complementary seven-bin analysis fixed the cluster sample within each filter
before estimating any radial slope. The complete-profile sample contains 21
`i` pairs/13 connected blocks and 53 `i2` pairs/28 blocks. The outer `i2`
50–75 kpc pointwise slope is positive, but its familywise max-T p-value is
0.410; the minimum p-value across all 14 filter/radius tests is 0.237. The
negative `i` inner slopes likewise fail familywise correction. This avoids
changing cluster support from bin to bin but does not fix member-detection or
mask changes. The association remains exploratory. See
[`icl_satellite_fixed_profile_assessment_2026-11-06.md`](icl_satellite_fixed_profile_assessment_2026-11-06.md).
Result summary hash `08e60ec574a624f4c177c7fbbc66c9ee93a57a4072cbea1755acf23aaeda113d`.

## Fixed-sample high/low ICL profile contrast (2026-11-06)

Using the same 21 `i` and 53 `i2` complete-profile clusters, a median ICL split
compared adjusted high-minus-low satellite excess across the same seven radial
bins. At 50–75 kpc the contrasts are +4.61 nJy arcsec⁻² (`i`) and +1.35
nJy arcsec⁻² (`i2`); their 14-test max-T p-values are 0.896 and 0.967. The
smallest familywise p-value is 0.452. No radial bin is resolved. This is a
complementary estimand on the same selected support, not an independent
replication. Assessment:
[`icl_satellite_fixed_profile_split_assessment_2026-11-06.md`](icl_satellite_fixed_profile_split_assessment_2026-11-06.md).
Summary hash `49403d38eb267f61bee9812c7b2a24694a2fa45b5463fb2db60c13cbd4d1f909`.


### 2026-11-06 continuation: 16-pixel ICL covariance sensitivity

Added an Astrafek runner extending the 224 supported target/control pairs from
±8 to ±16-pixel variance-normalized residual covariance. All 448 frames
completed; the 129,472 lag cells overlapping the previous ±8 result are exactly
identical. Contracting the 33×33 kernels with the exact saved 200–300 kpc
annulus masks gives finite-lag stack sigma 0.06311 nJy arcsec⁻² in `i` and
0.06290 in `i2`, versus 0.07931/0.07449 for 17×17. Paired 5° block-bootstrap
5,000-replicate intervals for the sigma change are −21.3% to −19.1% (`i`) and
−21.5% to −9.4% (`i2`). These conditional intervals do not include uncertainty
in the lag estimates.

The 16-pixel square-kernel sum is below one because negative long-lag shells
dominate after affine-plane removal and separate per-lag patch medians; this
assembled median kernel is not guaranteed positive-semidefinite. Treat the
change as estimator sensitivity, not a physical noise correction. The finite-
lag term is only 6.1%/5.9% of the saved spatial 95% interval half-width in
`i`/`i2`. Profiles remain unchanged and `science_measurement=false`.
Assessment `analysis/icl_wgl_membership/hsc_debug/icl_covariance_16px_aperture_assessment_2026-11-06.md`.
Kernel, aperture, and block-bootstrap summaries are under `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_supported_covariance_to_16px_2026-11-06_01/`,
`icl_aperture_covariance_200_300_16px_2026-11-06_01/`, and
`icl_aperture_covariance_16px_block_bootstrap_2026-11-06_01/`. No HSC access,
additional download, GPU, Slurm, commit, or push.


### 2026-11-06 continuation: affine-plane covariance-null calibration

Added `analysis/icl_wgl_membership/hsc_debug/analyze_icl_plane_subtraction_covariance_null.py`
to quantify lag anticorrelation induced by the local residual-plane fit. For an
unmasked, homoscedastic 64×64 white-noise patch, the exact OLS projection
kernel yields square-lag sums 0.8087 at radius 8 and 0.3629 at radius 16. A
100-realization simulation through the actual clipped estimator gives means
0.848 and 0.358, consistent within its broad finite-simulation spread. The
observed HSC target/control medians are 1.56–1.66 at radius 8 and 0.74–0.83 at
radius 16. The ideal null confirms that plane removal creates negative tails,
but its predicted 8-to-16 decrease is smaller than the observed ≈−0.8 shell;
it does not explain the entire HSC result. Masks, heteroscedasticity, sky
structure, source masking, and the per-lag median remain entangled. No ICL
correction or science promotion follows; `science_measurement=false`.
Assessment: `analysis/icl_wgl_membership/hsc_debug/icl_plane_subtraction_covariance_null_assessment_2026-11-06.md`.
Runner SHA-256 `150fedeb2a8e00ec011dea525f0bf3c71c31fe9d66f31dc1bc0463724f9a9aee`;
summary SHA-256 `bf93d7c0b9189782c9a510cffee6964a45ab87b49bb7d5630ecc6822419b2760`.
No HSC request or download.

## Same-visit split-stack covariance pilot (2026-11-06)

Added `analyze_icl_visit_split_covariance_pilot.py` and used the already
downloaded warps for supported `i2` cluster 101 and its matched control. The
same seven visit IDs were divided into alternating sorted halves at both
centers, then inverse-variance stacked on the coadd grid through the existing
affine-WCS/nearest-neighbor mapping. One visit per center had no valid pixels,
leaving usable split counts of 3/2 for the target and 4/2 for the control.
Clean common annulus coverage is 67.5% and 58.3%, respectively. The normalized
split-difference RMS before plane removal is 1.093/1.166; after local plane fit
and clipping, zero-lag covariance is 1.182 in both centers. The ±8-pixel
cumulative lag sums are 0.948 (target) and 1.187 (control). With one pair,
unequal split depth, approximate nearest-neighbor mapping, and no matched
white-noise propagation yet, this is a diagnostic only. No correction or
association claim is justified; `science_measurement=false`. Assessment:
`icl_visit_split_covariance_pilot_assessment_2026-11-06.md`. Runner SHA-256
`1f28f66ed9477abbef03bb53188f0cb14cd9354498a94eaf43989802be68c929`; summary
SHA-256 `09e557de604d4161cd0281065a59f32669744ab480b180029d192f79ba192359`.

## Supported-sample visit-split covariance audit (2026-11-06)

Scaled the pilot to all complete pairs in the approved 94-pair warp download:
90 pairs (24 `i`, 66 `i2`), with four centers missing from the acquisition.
The first aggregate output incorrectly treated zero-support lag cells as zero;
the hash-pinned `_02` rerun requires valid pixel-pair support at every one of
the 289 cells before using the ±8 cumulative kernel. Use only the `_02` output.
The corrected common complete-case sample is 21 `i` pairs/11 sky blocks and 64
`i2` pairs/33 blocks. Four-metric within-filter max-T intervals show no resolved
target/control contrast in normalized split RMS or ±8 cumulative covariance.
The i2 zero-lag difference is +0.0060 (simultaneous 95% interval +0.0010 to
+0.0110) and clean annulus coverage is +0.074 (+0.006 to +0.141), both
indicating a small residual matching/coverage difference. Split-RMS medians are
1.082/1.082 (`i`) and 1.085/1.079 (`i2`) for target/control; these also include
sky/source residuals and are not a detector-noise calibration. Visit-ID Jaccard
overlap medians are 0.80/0.845 in `i`/`i2`. Nearest-neighbor mapping, unequal
visit support, residual sky/source terms and the estimator remain unresolved.
No correction or science promotion; `science_measurement=false`. Assessment:
`icl_visit_split_covariance_sample_assessment_2026-11-06.md`. Sample-runner
SHA-256 `bc2c001d4f4715d74a808fc893d8d53227c3dc12598fbcfe57454fdb871b37f7`;
joint-bootstrap runner SHA-256
`64c36f4724d71030936f7b765cc64438a0ba8e96f5e4b35c2f8afdb9727db903`;
sample summary SHA-256
`82f7ebaa199dca4d278152969fe2a37bf5e8b4a3380bc999a74e6a276ef4f35c`;
joint-bootstrap summary SHA-256
`973d2def87d97e8c5ef5a1527ddd2e66a76fb8bdeb0dc52ac3c6e7fca0a68f26`.

## Heteroscedastic white-noise null for the visit-split estimator (2026-11-06)

Added `simulate_icl_visit_split_covariance_null.py` and ran 100 realizations
for target/control cluster 101. Independent Gaussian noise was generated on
native warp pixels from each supplied variance plane, then passed through the
observed masks, visit split, nearest-neighbor mapping, inverse-variance stacks,
local affine-plane fit, iterative clipping, and the ±8 lag estimator. The null
produces raw normalized RMS ≈1 and fitted zero-lag covariance ≈0.997 at both
centers, but the plane estimator lowers the cumulative ±8 sum to ≈0.715. The
observed values (target/control) are RMS 1.093/1.166, zero lag 1.182/1.181,
and ±8 sum 0.948/1.187; all exceed the simulated central 68% interval. Thus
independent pixel white noise plus supplied variance and this approximate map
does not explain the observed split residuals. Correlated warp noise, residual
sky/source structure and variance calibration remain unresolved; no correction
is applied, and `science_measurement=false`. Assessment:
`icl_visit_split_covariance_null_assessment_2026-11-06.md`. Runner SHA-256
`7c4dcb02f7d6b0b31cd3a1c85f47923c2e6e9ffbef75520f0687143b64d57410`; summary
SHA-256 `d9dc7090b1cee14fdb88a020adc5164a6da610bcaa76d9d190142b8a2ad86890`.


## Annular covariance comparison (2026-11-06)

The 94-pair coadd audit compared the 200–300 kpc ICL annulus to 350–500 kpc. All paired target-control and radial difference-in-differences 95% block-bootstrap intervals include zero. The complete inner-annulus samples are only 8 `i` and 16 `i2` pairs after strict mask/patch cuts; this is low power and does not establish equivalent covariance or support radial transfer. The outer ring may include cluster-envelope light. No correction was applied; `science_measurement=false`. Assessment: `icl_annular_coadd_covariance_assessment_2026-11-06.md`; output `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_annular_coadd_covariance_2026-11-06_01/`; summary SHA-256 `c03b1c60fcc3415aee615c8685700b42633fa769df0ac4a8173013bafe2f93b0`.

The exact random-point variance export plan remains previously approved, but its prior credential-gated submission attempt ended without a CAS job ID or data transfer. Its checkpoint remains zero-job/pending; no duplicate request was submitted during this annular audit.


## Familywise ICL–satellite–WL association audit (2026-11-06)

The saved 27-lens i2 triple-overlap analysis reported 29 conditional relationships, including two pointwise patterns. Added a shared-block restricted Webb wild-cluster max-T test (10,000 draws): the smallest familywise p is 0.672 for WL cross versus satellite contrast at 0.4–0.7 Mpc; the tangential 2–3 Mpc pattern has adjusted p=0.923. No association survives familywise adjustment. Point estimates reproduce the prior table exactly. The small selected sample and approximate 20 connected sky blocks mean this is low-power exploratory evidence, not evidence of no physical relation. `science_measurement=false`. Assessment: `icl_satellite_wl_triple_overlap_max_t_assessment_2026-11-06.md`; output `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_wl_triple_overlap_maxt_2026-11-06_01/`; summary SHA-256 `d5e2d3aa05be32b07be5fa62b8de467aea5314411daba49cb733444850126b0e`. No HSC query or download.


## Updated consolidated evidence report (2026-11-06)

Extended `build_icl_wl_current_evidence_report.py` to validate and include the annular coadd covariance table and the 29-way ICL–satellite–WL max-T metrics. New report: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_wl_current_evidence_2026-11-06_19/icl_wl_current_evidence_report.md`; report SHA-256 `5a4d2f7517e3fa78625c19ec96e51dbe8920c0b3c1bcf6ef3293edb6a7bf5426`; summary SHA-256 `ff2d720b09d6f0112670c956dad9053d440b02508f6ce2f4484d6a28b12d941f`; builder SHA-256 `af028a910fe9e5308d4b382029544f9968818e9dd77c97ca08ba1cbd3298418a`. The report pins 28 input artifacts and remains `science_measurement=false`. The already-approved random-point variance export is still pending with zero CAS jobs; this session confirmed no `HSC_USER`/`HSC_PASSWORD` environment variables and no live submission process, so no request was sent.


## i2 ICL–satellite–WL overlap selection (2026-11-06)

The 68-lens i2 ICL–satellite sample attrits to 27 verified random-subtracted WL lenses (60 in the matched 620 parent, 35 with target shape rows, 31 with positive random counts). Across 33 connected sky blocks, only nine contain both selected and excluded lenses; within them there are 10 selected and 16 excluded, leaving 17 selected and 25 excluded in single-status blocks. Spatial positivity therefore fails for transporting the 27-lens relation to the full 68, regardless of a fitted propensity model. Four-variable restricted-block max-T finds no adjusted imbalance in redshift, richness, ICL, or satellite contrast, but the nine mixed blocks make the audit low power. Target shape-pair count is a direct eligibility gate (median 4,342 selected versus 0 excluded), not an independent physical result. Do not IPW-correct to the parent sample. `science_measurement=false`. Assessment: `icl_wl_overlap_selection_i2_assessment_2026-11-06.md`; output `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_wl_overlap_selection_i2_2026-11-06_03/`; summary SHA-256 `fef9cdbfdd595f770d1a1a933d29618443524896a318298266a6a1134cb92a98`.

The 31-input consolidated report now includes this audit and the preceding covariance/max-T results: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_wl_current_evidence_2026-11-06_20/icl_wl_current_evidence_report.md`; report SHA-256 `dae6062a10ca70f0b2eeb6f2fdfbcae76e52665a1bf2ffac31eea9cfa9a7826a`, summary SHA-256 `c665590ba37c5ecbbd9ba08d291434d5f98fdc7d5d6dfd1d8c3fef5c4cea6c2a`; builder SHA-256 `507c0d07e6d545c5afe8f34f6f59d3ab5260129e94aadb2be7245017dd5bba80`.


## i2 ICL–satellite–WL overlap selection audit

The 68 i2 ICL–satellite-supported systems reduce to 27 verified random-subtracted WL lenses. Only nine of 33 connected 5-degree blocks include both selection groups; these contain 10 selected and 16 excluded lenses, while 17 selected and 25 excluded lie in single-status blocks. Thus spatial positivity is absent for transporting the 27-lens relationship to the 68-lens parent. Four-variable within-block max-T shows no adjusted redshift, richness, ICL, or satellite-contrast difference, but nine mixed blocks make the test low-power. Target shape-pair count is a direct eligibility gate (median 4,342 vs 0) and is not an independent physical imbalance. Do not apply IPW to represent the parent sample. Keep `science_measurement=false`. Assessment `icl_wl_overlap_selection_i2_assessment_2026-11-06.md`; result `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_wl_overlap_selection_i2_2026-11-06_03/`; summary SHA-256 `fef9cdbfdd595f770d1a1a933d29618443524896a318298266a6a1134cb92a98`.

Consolidated report version 20 now hash-validates this audit, annular covariance, and the ICL–satellite–WL max-T results. Report `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_wl_current_evidence_2026-11-06_20/icl_wl_current_evidence_report.md`; report SHA-256 `dae6062a10ca70f0b2eeb6f2fdfbcae76e52665a1bf2ffac31eea9cfa9a7826a`; summary SHA-256 `c665590ba37c5ecbbd9ba08d291434d5f98fdc7d5d6dfd1d8c3fef5c4cea6c2a`.


## Lanczos-5 warp-to-coadd variance sensitivity (2026-11-06)

Applied a separable Lanczos-5 mapping to the already downloaded warp variance planes for the same 94-pair sample. Of 188 requested centers, 184 were analyzed; four repeated HTTP 404 gaps remain. Median 200–300 proper-kpc coadd/predicted variance ratios are 1.00707/1.00760 for i target/control and 1.00639/1.00611 for i2. Paired block-bootstrap target-minus-control intervals include zero in both filters. Median per-center ratio shifts from nearest-neighbor mapping are −0.000245 in i and −0.000161 in i2. Strict shared-mask coverage medians range from 69.0% to 71.9%. HSC PDR3 documents fifth-order Lanczos resampling for CCD-to-warp processing, but this secondary mapping only propagates diagonal warp variances; it does not recover upstream correlated warp noise or exact coadd weights, clipping, and mask interpolation. This is not an exact covariance estimate or noise correction; `science_measurement=false`. Assessment `icl_supported94_warp_noise_lanczos5_assessment_2026-11-06.md`; output `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_supported94_warp_noise_lanczos5_2026-11-06_02/`; summary SHA-256 `16c8afbcc8287444860bacdcfd01038105bd51d5e82d61e44d8b82d8853e9f74`; table SHA-256 `a0322cc88e63768d172ccc8da52ce5148bef1a4407232b53e93e00bc9ae0f06f`.

The consolidated evidence report was rebuilt to include the Lanczos-5 sensitivity with input hash checks: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_wl_current_evidence_2026-11-06_22/icl_wl_current_evidence_report.md` (SHA-256 `2ffcb4383e9923f4a904f4b6b684a5c9c03d3e3fc52d5bafe23c8e1f2b0a6389`); summary SHA-256 `51ac24ee822eedd5b8af107ef67ab9b3b893686d422a5604b59411c6482b40cf`; builder SHA-256 `cdd70be6227f272adefc1343c1c19c5d0911fb3f26718763fd3c7c8b0c57eaef`. It records 34 pinned input artifacts, no HSC network or image-pixel access during report generation, and keeps `science_measurement=false`.

## 2026-11-06: full-sample annular coadd covariance scale sensitivity

The delivered-coadd covariance diagnostic now covers all 783 same-filter target/control pairs. Complete 64×64 patch support was only 17/302 i and 48/481 i2 pairs, with strong low-redshift selection; familywise intervals included zero. A 32×32 sensitivity (up to 64 patches, minimum 8 per center, 80% valid, lag radius 8) increased support to 215/302 i and 402/481 i2, but redshift SMD remains −0.65 and −0.39. Its within-filter familywise intervals identify an i lag-sum difference and an i2 zero-lag difference, while coverage differs in both filters. Since estimates change with patch scale and support/coverage, these are exploratory delivered-coadd diagnostics, not a pure detector-noise calibration or correction. `science_measurement=false`. See `icl_full783_annular_covariance_32px_assessment_2026-11-06.md`; the consolidated report will be rebuilt after the approved random-point variance export completes.

## 2026-11-06: approved random-point variance export completed

The pinned HSC PDR3 random-point query completed all 16/16 batches: 40,048 rows
(11,138,863 CSV bytes), versus a data-dependent estimate of 46,424. Points were
returned for 1,587/1,598 target/control centers, including 791/800 targets and
796/798 controls; both endpoints have `i_pix_variance` for 787 matched pairs.
The independent audit verified row contracts, all batch hashes, center/role
mapping, and the 200–300 proper-kpc annulus. This is sparse random-point local
variance/sky/PSF metadata, not a pixel weight map or complete covariance model;
no ICL correction follows and `science_measurement=false`. Assessment:
`icl_random_pixel_variance_export_assessment_2026-11-06.md`; audit summary
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_random_pixel_variance_audit_2026-11-06_01/summary.json`
(SHA-256 `96245d565375635a994c93ec3294fc581d007161c4479232d5091bbd3cf52389`).

The post-download consolidated report, including the 64/32-pixel 783-pair
covariance audits, variance-export integrity and center-support audit, is
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_wl_current_evidence_2026-11-06_23/icl_wl_current_evidence_report.md`
(SHA-256 `14ebf8e0965e12f9dfc5aabad4bc359d8e670428de2919c643ca3830f2d1936d`);
summary SHA-256 `5fe134a256ea35e4a7780a796f9b6c38d6cc460432dcddad956c314b318a5cfc`.

## 2026-11-06: random-point conditions versus ICL–satellite association

Joined the audited 40,048-row random-point export to the existing 94-pair
ICL/satellite support table (26 `i`, 68 `i2`). Eight one-at-a-time terms
(target/control log ratio and pair-mean log value for pixel variance, sky RMS,
input count, and PSF trace) were tested beyond redshift and richness with 5,000
connected-block draws and one 16-test max-T family. The baseline partial r is
0.079 in `i` and 0.293 in `i2`; the i2 adjusted values range 0.260–0.295. All
change intervals include zero, with minimum max-T adjusted p=0.951. The data
show no resolved relation change from these local-condition proxies, but the
sample is selected, `i` has only 26 pairs, and random-point median uncertainty
is not propagated. No correction or science promotion; `science_measurement=false`.
Assessment `icl_random_point_variance_satellite_relation_assessment_2026-11-06.md`;
summary `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_random_point_variance_satellite_relation_2026-11-06_01/summary.json`
(SHA-256 `0c75520ad6b04dd500ba66ae10982e30d6fa9758379615a4372e7aa129569f0f`).

The current 69-input integrated evidence report is
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_wl_current_evidence_2026-11-06_24/icl_wl_current_evidence_report.md`
(SHA-256 `bb45e9d3ed817222debbe0d5aff2bad54582eb3b8d71b974234e58e933a0d122`);
summary SHA-256 `87a7bed6f574ea12a1946f16d2b9da7654b258f02c310c1e532612cc30ddb78f`.

## 20/20 forced-profile support and ICL outcome selection (2026-11-06)

The corrected paired analysis retains all 224 ICL-supported cluster/filter pairs
in the selection ledger: 71 `i` and 153 `i2`. The 20/20 target/control
coverage floor supports 50 `i` and 110 `i2` pairs. Selected pairs have much
higher median ICL 200–300 kpc residual than excluded pairs (3.060 vs 1.116
nJy arcsec⁻² in `i`; 4.446 vs 1.853 in `i2`), while median redshift and
richness are similar. Support is therefore outcome-related, and the selected
subset is not representative by assumption.

The full 20/20 partial correlations (redshift and log richness controlled) are
−0.088 (`i`, 50 pairs) and +0.102 (`i2`, 110 pairs), compared with +0.082 and
+0.285 in the respective common detected-only subsets. Paired connected-5°
block-bootstrap intervals for all decomposed changes include zero. The i2
point shift is mostly associated with sample expansion; the common-sample
forced-measurement point change is near zero. For `i`, the forced-measurement
point shift is negative but uncertain. The threshold is post-hoc, intervals
are exploratory, and profile/ICL measurement uncertainty is not propagated.
No correction or detection follows; `science_measurement=false` remains.
Assessment `icl_satellite_forced_coverage_both_sides_assessment_2026-11-06.md`;
output `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_forced_coverage_both_sides_2026-11-06_03/`
(summary SHA-256 `b8b79a65e1c838e757fc4ee6c1e1fdf61340d72e989bf035ab29bc0a7a6a314a`; runner SHA-256 `7e9ec1672e0ac000e5e28cf9810c9f98b1e538e60465664853b1dc9a0fe5609f`).

## Consolidated report refresh after forced-coverage sensitivity

The latest evidence synthesis is version 27 at
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_wl_current_evidence_2026-11-06_27/icl_wl_current_evidence_report.md`
(report SHA-256 `af88c214eadad945c5177e701f4d7496f95222b297ad84e0be2ba91f595f8b0f`,
summary SHA-256 `42780af71c7392f90f80d53057a8efefe5f9062a84e8f08c595dbbc045eb89b5`).
It hash-pins 77 artifacts and incorporates the forced-coverage support-selection
result. All input hashes were verified; no HSC access or image-pixel read
occurred during report generation. The synthesis remains exploratory with
`science_measurement=false`.

## Spatial overlap of 20/20 satellite-profile support

Coordinate-linked connected paired-sky blocks contain both selected and
excluded pairs in 9/18 `i` blocks and 21/40 `i2` blocks. The equal-block mean
selected-minus-excluded ICL residual is +1.80 nJy arcsec⁻² (`i`, 95% block
interval −3.26 to +6.66) and +3.73 (`i2`, +0.55 to +6.93). Thus the positive
i2 support/ICL difference persists within mixed blocks in this descriptive
analysis; the i result is imprecise. Redshift and richness are not adjusted
within blocks, and single-status blocks prevent parent-sample transport. This
post-hoc diagnostic does not supply a correction and leaves
`science_measurement=false`. Assessment
`icl_forced_coverage_spatial_selection_assessment_2026-11-06.md`; summary
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_forced_coverage_selection_blocks_2026-11-06_02/summary.json`
(SHA-256 `3eb63f02507fa83cf1eea377ebc3574bd2b66e86f204bb76face995dcf5f5742`).

The latest integrated 80-artifact report is
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_wl_current_evidence_2026-11-06_28/icl_wl_current_evidence_report.md`
(report SHA-256 `33e41f660a70d3a81114b11801b47793fbdad161b8087da87c3a2458d1417fc3`;
summary SHA-256 `bc97fe8ba720c1edb90538d09cbdb918cfb9ac976f5b8cf8a4001a75028351f8`).

## Covariate-adjusted forced-profile selection contrast

Within mixed connected 5-degree blocks, the selected-minus-excluded ICL residual
coefficient adjusted for redshift and log richness is +1.47 nJy arcsec⁻² in
`i` (95% block interval −2.10 to +4.79; 9 blocks) and +2.13 in `i2` (−0.64 to
+5.91; 21 blocks). Both intervals include zero, so the unadjusted i2 difference
is not resolved after covariate adjustment. This selected-sample analysis does
not correct for missing profile support or transport to single-status blocks.
It remains exploratory and `science_measurement=false`. Assessment
`icl_forced_coverage_adjusted_selection_assessment_2026-11-06.md`; output
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_forced_coverage_selection_adjusted_2026-11-06_01/summary.json`
(SHA-256 `7e3c57cbd95b32986a7f52b707429f682a6d4e358e12724f9866f14f69b6f07c`).

Latest consolidated evidence report is version 29 with 82 hash-verified inputs:
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_wl_current_evidence_2026-11-06_29/icl_wl_current_evidence_report.md`
(SHA-256 `8350d58462b4eabd228920d970de9852d5b66272e58b629a52577a280dff8d10`).

## 2026-11-06: forced unmatched-position radial analysis

Added `analyze_icl_forced_position_radial_relations.py` and its assessment to
measure whether the 200–300 proper-kpc ICL residual covaries with light at
forced unmatched CAMIRA positions across radial bins. The analysis uses the
existing post-hoc 20/20%-coverage sensitivity profiles, adjusts for redshift
and log richness, clusters covariance over connected paired 5-degree blocks,
and uses 10,000 shared restricted Webb wild-cluster bootstrap-t draws with
Holm and max-T correction. Of nine supported tests, none has max-T p<0.05. The
i 50–75 kpc nominal Webb p=0.030 has Holm p=0.274 and max-T p=0.059. These
locations lack unique SEP identities, so the outcome is unmatched-position
light rather than confirmed satellite light. Selection and omitted
measurement/noise/PSF/sky terms remain material. No correction or promotion;
`science_measurement=false`.

The integrated report builder now validates the radial runner hash, summary,
three result tables, five source hashes, and the S23B submission audit. Version
30 has 91 verified inputs at
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_wl_current_evidence_2026-11-06_30/icl_wl_current_evidence_report.md`
(report SHA-256 `62a41b863e1ee2d2150c8647829626972fb3a99d5b123cd993289cc130a65515`;
summary SHA-256 `8d166119df6d21a1d51d3f26a9514e38c6bda8aee727b33d7b4edb1d6e1d2b63`;
builder SHA-256 `0028c1420146454be8b8ce765e957ea133bf41b6e47695c882d7a23330cec713`).
All input hashes were rechecked. Report generation made no HSC request or pixel
read; the S23B request got HTTP 401 before any CAS job or result row.

## 2026-11-06: S23B HSC-I metadata coverage completed

The previously approved S23B metadata scope completed after correcting the CAS
host and using the PDR3 mosaic `areacube` footprint operator documented by HSC.
The exact scope remained 1,333 candidate centers, 11,997 center/boundary
probes, 14 metadata columns, 27 batches, and no image/source/shape/photo-z
data. All batches passed row-contract and checksum audits; 13,575 metadata rows
were returned. Center metadata exists for 1,048 candidates; 1,043 have matches
at all nine probes, with image-product and background-product flags true at
each point. These are metadata flags only and do not verify image WCS, masks,
full-area support, or controls. The separate image step remains approval-gated.
Assessment `icl_s23b_coverage_metadata_assessment_2026-11-06.md`.

The current 124-input report is version 31 at
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_wl_current_evidence_2026-11-06_31/icl_wl_current_evidence_report.md`
(report SHA-256 `3d6aebd01cb2727ca21d2617ef2de1fd17a85c526199315510710f167a1acf15`;
summary SHA-256 `2a0b4edfa9a6536e1f19491365cc5c397a49106b6e1a81d547895b9d9a2f5f9e`;
builder SHA-256 `c7fd66ed582e948b002528bd0b72661608dc5cc2ecf56b209219cf53656ba47f`).
Every pinned input hash was independently checked; the synthesis made no HSC
request and read no image pixels. `science_measurement=false`.

## S23B matched-control candidate metadata plan (2026-11-06)

To advance beyond the existing 800-cluster image sample, prepared an offline
HSC PDR3 random-point query for the 1,043 S23B candidates that passed the
nine-probe HSC-I metadata screen. The query returns at most 200 random-position
candidates per target within 30 arcmin, with 29 location/ID and scalar I-band
quality/depth fields; maximum 208,600 rows across 105 batches, estimated 80 MB.
It excludes images, detected-source catalog rows, shapes, and photo-z. The plan
and SQL hashes were checked offline before submission. Plan SHA-256:
`143a1bd0a3a3c341576ce420f77ba0089fbfe80fdf58fbfd7afe30ea8e3f6714`.
Generator: `prepare_s23b_icl_control_candidate_plan.py`; plan directory:
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_icl_control_candidate_plan_2026-11-06_01/`.
This candidate metadata enables local depth/seeing and cluster-contamination
matching, but does not itself provide matched controls. The exact plan was
submitted only after user approval. Its submitter is pinned to the exact
SHA-256 and validates every returned row against the 29-column contract. The
subsequent image download requires separate approval. Keep
`science_measurement=false`.

## 2026-11-06 follow-up: satellite outer profile and approved S23B query

The existing 800-cluster profile archives were used to relate the ICL radial
contrast (200–300 kpc minus area-weighted 300–500 kpc) to the detected-member
satellite-centered 50–75 kpc target-minus-azimuthal-control profile, controlling
redshift and log richness. The partial correlations are −0.140/−0.126 in `i`
(20 pairs, 10 connected blocks) and +0.417/+0.391 in `i2` (47 pairs, 30
blocks). All four 95% connected-block intervals include zero. One unsupported
i2 satellite profile remains missing. This is selected-sample exploratory
evidence; source-Poisson, measurement, exact coadd covariance, sky transfer,
and extended PSF-wing uncertainty are not propagated. Assessment:
`icl_satellite_outer_profile_relation_assessment_2026-11-06.md`.

The integrated report is version 33 with 940 independently hash-verified
inputs at
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_wl_current_evidence_2026-11-06_33/icl_wl_current_evidence_report.md`
(SHA-256 `1ef9acc5bff699fed8e1d75131192cc03e88d48832e8ea485918b70c94832742`).
Its science gate remains false.

After explicit user approval, the exact S23B candidate metadata plan above was
submitted from LagEunha. Follow-up on 2026-11-07: batches 000–079 returned and
passed validation (160,000 rows total); batch 080 was submitted in live session
3123. The checkpoint is
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_icl_control_candidate_results_2026-11-06_01/status.json`.
Checkpoint flags confirm that credentials were not saved and that images,
detected-source rows, shapes, and photo-z were not returned. The export remains
limited to the approved random-position and scalar I-band QA/depth fields.
Target-cutout image receipt was subsequently approved and is recorded below.

The local-only postprocessor `audit_s23b_icl_random_control_candidates.py` is
prepared to independently verify the completed batch files, screen position
flags, depth/seeing/SFD extinction compatibility, and exclude positions near
the union of the local S20A and approved S23B CAMIRA catalogs. Its output is
provisional: selected 1-Mpc cutouts still need image-level WCS, mask, sky,
variance, and PSF checks. The audit has not been run because the metadata export
is still in progress; image retrieval is outside its scope.

## 2026-11-07: user-approved S23B target image receipt

`prepare_s23b_icl_target_image_manifest.py` prepared 1,043 screened S23B CAMIRA
centers for HSC PDR3 `pdr3_wide` HSC-I-family `coadd/bg` cutouts with IMAGE,
MASK, and VARIANCE planes and a 1,000 proper-kpc half-width. The target manifest
SHA-256 is
`e487582ba40cd3e58da67704242db4d7219e381fb08f4b33fed47ce4ff268b44`. Estimated
payload is 58,311,064,944 bytes (54.31 GiB); the required free-space margin is
75,804,384,428 bytes. `download_pilot.py` was extended to accept up to 5,000
sample frames while retaining per-file size limits and FITS/HDU/WCS checks.

Download session 5633 is active in
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_icl_target_image_plan_2026-11-07_01/target_download_01/`.
The latest checkpoint records 121/1,043 images and 6,826,320,000 bytes validated;
credentials are not saved and `science_measurement=false`. This target-only
scope uses CAMIRA centroids, not verified BCG centroids; no member catalog is
included, and no control images are included in this plan. The matched-control
metadata query remains active separately in session 3123.

`audit_s23b_icl_target_images.py` is prepared for the completed download. It
will recheck all 1,043 FITS hashes and IMAGE/MASK/VARIANCE layout, HSC-I filter,
WCS alignment, and cutout footprint before stacking. It is not run against the
partial transfer.

Live checkpoint follow-up (2026-11-07): session 5633 remains active; its
checkpoint records 373/1,043 target cutouts (21,290,526,720 bytes) validated.
Session 3123 completed all 105 metadata batches. No control images have been
requested or included. The completed query confirms no saved credentials,
image, detected-source, shape, or photo-z returns, and `science_measurement`
remains false. Its audit and control-image plan follow below.

The random-control candidate export subsequently completed in session 3123:
105/105 batches, 208,600 rows. The local audit completed with 164,984 unique
random positions, 148,682 center-QA passes, and 132,886 metadata matches. At the
1 proper-Mpc CAMIRA exclusion, 1,041 non-overlapping provisional controls are
available for 1,041 targets; sensitivity is 1,041 at 0.5 Mpc and 1,037 at
2 Mpc. Audit files and hashes are recorded in `HANDOVER.md` and under
`icl_s23b_control_candidate_audit_2026-11-07_01/`.

The offline control image plan was prepared and specifically approved:
1,041 HSC PDR3 `pdr3_wide` HSC-I `coadd/bg` IMAGE/MASK/VARIANCE cutouts with
1,000 proper-kpc half-width; maximum estimate 58,209,418,068 bytes (54.21 GiB),
requiring 75,672,243,489 bytes free with a 30% margin. Manifest SHA-256 is
`4e182a3160055f1dfec4171035bceb1f4a8b574d10f0b12660f67a39dd037c35`.
Control-image session 81219 stopped after 14/1,041 validated files (876,274,560
bytes), when `CTRL_S23B_456` returned HTTP 404 at tract 9718/patch 507. No
partial file remains. The first non-TTY attempt stopped before network access;
the later TTY `--resume` preserved its checkpoint. Credentials are not saved.

`replan_s23b_icl_control_manifest_after_404.py` created a revised manifest
offline, replacing only the unavailable `S23B_456` random point with ID
`42740795815504337` at RA 8.842668812326, Dec 2.786488591937 (tract 9718/patch
407). It passed the same 1-Mpc CAMIRA exclusion, metadata match, and global
footprint non-overlap screens. V2 manifest SHA-256 is
`ae9eb406cfc179d9980031bb67f85f8ee71c4580ba59e97713278823db5c06df`; estimate
and 30%-margin storage requirement are unchanged. The user approved this exact
v2 hash. The original 14-file v1 failure checkpoint is preserved as
`download_status_v1_404_checkpoint.json`; the resume checkpoint now binds to
v2. All 14 existing files match unchanged v2 rows and were revalidated on
resume. Control-download session 56526 is active; its latest checkpoint has
86/1,041 files and 5,233,340,160 bytes, including the replacement position.
Credentials are not saved.

Target-image session 5633 remains active with 850/1,043 images (47,507,889,600
bytes) validated at the latest checkpoint. The partial filter count is 685 i2
and 165 i; this split may change before the transfer completes.
`df -B1 /gpfs/kjhan/Astrafek/data` reports 157,141,126,086,656 bytes available,
above the 30%-margin storage requirements for the approved target and control
image plans.

The separate CPU-only control-image verifier
`audit_s23b_icl_control_images.py` is now prepared. It pins the approved v2
manifest, verifies all 1,041 control files against the download checkpoint,
rechecks FITS/HDU/filter/WCS evidence, and audits the unique target/control map.
It has not been run against the partial transfer. The existing 2D stacker
assumes the earlier 800-target/798-control input contract; the new S23B images
will require their own audited pipeline. No stack or new pixel-quality claim
was made in this continuation.

Centering interpretation correction (2026-11-07): the S23B paper states that
CAMIRA selects a bright galaxy as the BCG and recomputes richness and cluster
redshift around that position, but the public seven-column table documents only
generic RA/Dec and gives no separate BCG coordinates. It is not established
whether the published RA/Dec fields are exactly the final algorithm-selected
BCG coordinates. The v1 image manifest's phrase “not a BCG position” therefore
overstates the evidence; interpret it as “not separately verified against a
BCG position.” The approved manifest and ongoing download are unchanged. The
generator `prepare_s23b_icl_target_image_manifest.py` now uses the more
cautious wording for any future plan. Sources: S23B method paper
https://arxiv.org/abs/2512.13954 and author catalog directory
https://github.com/oguri/cluster_catalogs/tree/main/hsc_s23b_camira.

The local S23B CAMIRA wide catalog contains 12,162 cluster rows and only seven
published fields (ID/name, RA/Dec, cluster photometric redshift, richness, and
BCG spectroscopic redshift); it has no satellite/member positions. File
`/gpfs/kjhan/Astrafek/data/catalogue_s23b_camira_2026-10-20_01/camira_s23b_wide_sm.tbl`
is 827,339 bytes, SHA-256
`5af057028feab7976cd3ed389a7246cbe73add8bbb96891101f0f692b1de4d71`. Member
data and validation of CAMIRA-centroid versus BCG-centroid alignment are needed
before relating this new ICL sample to satellites; neither was included in the
approved image scope.
The local Astrafek audit `audit_s23b_camira_catalog_content.py` verified the
file hash, all 12,162 rows, and the seven declared fields without network or
image access. Its summary is
`/gpfs/kjhan/Astrafek/data/s23b_camira_catalog_content_audit_2026-11-07_01/summary.json`
(SHA-256 `28f7836a0520f6b59a32241531a8e3ce5fc2e0bb7b71248dfad05f8870743ac8`).
The public `oguri/cluster_catalogs` S23B directory lists the cluster
catalog files and points separately to photometric LRG catalogs; it does not
list a CAMIRA member table. The linked LRG product has not been validated as a
satellite membership catalog.
The associated paper reports over 600,000 probable CAMIRA members, while
describing the photometric LRG catalog as independently constructed. The
public GitHub directory does not identify a downloadable member table; keep
this source unresolved and do not use the LRG list as confirmed members. See
https://arxiv.org/abs/2512.13954.
