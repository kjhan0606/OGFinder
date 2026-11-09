# S23B Gaia-supported empirical PSF-wing pilot (2026-11-07)

## Scope and decision

The user approved one Gaia DR3 query for the 492 morphology-only compact-source
positions in the existing 50-frame S23B target pilot (plan SHA-256
`394083a8cdaf0770376959e71c8bd46bd004c762dcb3e77a6256bac656ee52f0`). The
query used two batches (250 and 242 positions), a 2-arcsec match radius, and
the 18 approved astrometry, quality, variability, and classification fields.
It requested no images, HSC shape rows, or photo-z values. The response
completed on 2026-11-07 at 05:10:42 UTC with 77 candidate-source rows. Its
SHA-256 is
`1f1f8b947a3b23e1dad4675447247137d66d5f2586d6aea3a73989df3cb88e6b`.

The query is an identification aid, not PSF approval. An offline audit found
77 unique matched positions among 492 (15.7%); 415 positions had no Gaia
match. No position had multiple Gaia sources in the cone, and no Gaia source
was reused for two HSC candidates. A conservative diagnostic screen found 48
strict candidates (20 in `i`, 28 in `i2`): one Gaia source within 2 arcsec,
RUWE < 1.4, no duplicate/QSO/galaxy flag, significant parallax or proper
motion, and combined-star probability >= 0.5. These thresholds are not a
calibrated stellar classifier. Gaia non-matches do not establish that an HSC
source is non-stellar.

The header-only audit found no DATE, MJD, or EPOCH key in any of the 50 HSC
coadd headers. The effective HSC epoch therefore remains unknown, and the
Gaia positions were not propagated to an assumed epoch.

## Wide-stamp and profile method

The 48 strict candidates were checked against the existing HSC target images.
Thirty-six 401x401-pixel stamps were complete (16 `i`, 20 `i2`) across 22
frames; 12 candidates were blocked because the required stamp crossed an
image edge. Stamps retain the candidate star and mask detected neighbors,
known CAMIRA member positions, and bad pixels. They are neither
background-subtracted nor subpixel-recentered. No matched-control star sample
was used.

Four frame-weighted diagnostic stacks were formed, with 5,000 bootstrap draws
resampling frames rather than individual stars. Each per-star profile is the
annular mean minus a local background, divided by the masked core flux per
pixel. The two tested sky annuli are 120–150 and 160–190 pixels (20.2–25.2 and
26.9–31.9 arcsec); radial profiles extend to 96 pixels (16.1 arcsec). Two
neighbor-mask expansions (Kron 2.5 and 3.5) were evaluated. The `outer500`
ledger additionally requires the star to lie at least 500 proper kpc from its
CAMIRA cluster center.

Only 82–159 of the 528–880 profile rows per stack passed the geometric and
sky-coverage gates (15.5–22.6%). The resulting outer-bin estimates usually
depend on only 2–5 frames. The bootstrap intervals below are the 16th–84th
percentiles of frame-resampled means; they do not include sky-transfer,
coadd-covariance, selection, or PSF-model uncertainty. Values are normalized
mean flux per pixel, in units of 1e-6 core flux.

| Filter, star subset | 32–48 px, sky 120–150 / 160–190 | 64–96 px, sky 120–150 / 160–190 |
|---|---:|---:|
| `i`, all | −1.09 [−2.16, −0.02] / +0.13 [−0.89, +1.16] | +0.08 [−0.40, +0.57] / −1.42 [−1.93, −0.90] |
| `i`, >=500 kpc | −1.48 [−4.05, +1.10] / +0.53 [−0.50, +1.56] | −0.37 [−1.30, +0.57] / −1.42 [−1.93, −0.90] |
| `i2`, all | +3.09 [+0.23, +5.94] / +2.39 [+0.71, +4.06] | +1.68 [−0.20, +3.55] / +0.71 [−0.35, +1.66] |
| `i2`, >=500 kpc | +4.94 [+1.69, +8.19] / +2.39 [+0.61, +4.06] | +2.74 [+0.47, +5.01] / +0.71 [−0.35, +1.66] |

The frame sets differ between annuli and bins, so these are sensitivity
comparisons, not paired significance tests. `i2` shows a positive extended
component in this small pilot, but its amplitude declines by factors of about
2 at 32–48 pixels and 3.9 at 64–96 pixels when the farther background annulus
is used. The latter interval also crosses zero. In `i`, outer profiles change
sign with radius and background-annulus choice. The apparent `i2` wings are
therefore not yet separable from sky-annulus transfer at useful precision.

## 801-pixel follow-up

A geometry screen of the same 48 Gaia candidates found enough detector margin
for 801×801 stamps for 8 `i` and 14 `i2` candidates. The local CPU runner
produced 22 stamps across 13 frames; 26 candidates were blocked by image
edges. In the >=500-proper-kpc subset there are 4 `i` stamps from 4 frames
and 10 `i2` stamps from 8 frames.

Profiles now extend to 200 pixels (33.6 arcsec), with more distant background
annuli at 260–300 and 340–390 pixels (43.7–50.4 and 57.1–65.5 arcsec). In
`i2`, bins beyond 96 pixels are supported by only 1–2 frames. For the
>=500-kpc subset, the 96–128-pixel bin is −2.49 [−5.84, +0.87] and −3.00
[−8.46, +2.45] ×10⁻⁶ core-flux per pixel for the two annuli; the 160–200-pixel
bin is −2.33 [−6.59, +1.94] and −2.84 [−9.21, +3.52] ×10⁻⁶. The 128–160-pixel
bin has only one contributing frame in each annulus and no frame-bootstrap
interval. Only 96/560 profile rows (17.1%) pass geometry and sky-coverage
gates in `i2`; only 13/224 (5.8%) pass in `i`.

This smaller wide-stamp subset does not establish the positive inner `i2` wing
seen in the 401-pixel stack beyond 16 arcsec. Because the stamp subset and
background annuli differ, the two results are not a paired comparison.
Negative outer residuals may reflect over-subtraction if the sky annuli
contain PSF light; they do not establish a negative true PSF wing. The result
is that this sparse target-only sample cannot separate the extended wing from
its local background. No ICL correction follows.

The expanded stamp output is
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_psf_wing_stamps_2026-11-07_02/`;
the four profile outputs are
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_psf_wing_stack_{i_all,i_outer500,i2_all,i2_outer500}_2026-11-07_03/`.

## Scientific assessment

This pilot demonstrates that HSC target images contain Gaia-supported stars
with measurable profiles beyond the PSF core. It does **not** establish a
target/control-matched PSF, an effective coadd PSF, a stable extended-wing
correction, or an ICL surface-brightness correction. The sample is sparse,
limited to a stratified 50-frame subset, and has no matched-control stars.
The selected stars are not a representative sample of all 1,043 S23B targets.
Local sky annuli are close to the stamp edge and may include PSF wings,
residual cluster light, unresolved sources, or sky structure.

Keep `science_measurement=false` and `psf_approved=false`. Do not subtract
these profiles from ICL stacks or use them as a calibrated PSF correction.

## Next gates

1. Expand the matched-control PSF pilot beyond the 49-frame subset. The
   separately approved control Gaia query yielded 38 strict diagnostics, but
   no radial bin has two valid paired frame profiles. Any additional Gaia
   query needs its own exact approval.
2. Recover the effective HSC coadd epoch and propagate Gaia astrometry before
   assessing core-centering effects.
3. Compare the empirical profiles with available HSC model PSFs in the core,
   then run member-preserving ICL injection/recovery and sky/covariance
   transfer checks before reopening the science-measurement gate.

## Reproducibility

- Gaia audit runner:
  `analysis/icl_wgl_membership/hsc_debug/audit_s23b_psf_gaia_matches.py`
  (SHA-256 `4097d1fa06b105b5214a2f32ede3fc13c49e9a63634d13009e12b068dfc7c403`).
- Stamp runner:
  `analysis/icl_wgl_membership/hsc_debug/prepare_s23b_psf_wing_stamps.py`
  (401-pixel run SHA-256 `6200b88ab95b6347c38a5015e793ca822b7372c7f1db80c67210276e5f98c5b7`; expanded-stamp run SHA-256 `bd2d5708dbfdbb6b39654bb8a59ac36312a1e13775206e1c1180e596477e819c`).
- Stack runner:
  `analysis/icl_wgl_membership/hsc_debug/stack_empirical_psf_star_profiles.py`
  (401-pixel run SHA-256 `7a0898f2ddaf14f275a45d0b54fb9a568d7ba96eb693f4fc022794347ed5e238`; expanded-stack run SHA-256 `24a381edc1d5fe1a8702871c7635126a8b3d4c2c42c94dcc77fd268acdbf682f`).
- Gaia audit output:
  `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_psf_gaia_audit_2026-11-07_01/`.
- Wide-stamp output:
  `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_psf_wing_stamps_2026-11-07_01/`.
- Four stack outputs:
  `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_psf_wing_stack_{i_all,i_outer500,i2_all,i2_outer500}_2026-11-07_02/`.
- The `i2`, >=500-kpc figure is
  `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_psf_wing_stack_i2_outer500_2026-11-07_02/frame_weighted_psf_wings.png`.
- Expanded outputs use 801-pixel stamps, radial edges through 200 pixels,
  sky annuli 260–300 and 340–390 pixels, and 5,000 frame-bootstrap draws.

## Matched-control candidate, Gaia, and profile comparison (2026-11-07)

The same local morphology-candidate preparation was run on the already
received, seeing-matched random-sky control images. All 49 matched image pairs
were prepared locally (25 `i`, 24 `i2`); 48 control fields yielded at least one
candidate (24 per filter). The morphology screen found 458 control positions
versus 492 target positions. Candidate counts were 180 versus 229 in `i` and
278 versus 263 in `i2`. Median per-field counts were 7 versus 8 in `i` and 11
versus 10.5 in `i2`. The median detection-moment FWHM was 5.312 versus 5.229
pixels in `i` and 5.359 versus 5.434 pixels in `i2` (control versus target).
These descriptive morphology measurements do not establish stellar identity
or PSF suitability.

The nearest control candidate is 58.12 arcsec from any of the 492 target
candidates; no target/control pair lies within 2 arcsec. No two control
candidates lie within 2 arcsec of one another, so the proposed Gaia cones are
not duplicated at the adopted match radius.

A separate Gaia DR3 plan covered exactly the 458 control
positions in two query batches (250 and 208 positions; 2-arcsec cones). It
requests 16 Gaia source fields: source ID, position/reference epoch,
parallax-significance and proper-motion fields, RUWE/duplication indicators,
G magnitude/variability, QSO/galaxy flags, and combined star probability.
Candidate indices and exact separations are attached locally; the final local
table has 18 columns. The rough estimate, scaled from the completed target
query, is 72 matched rows and 14,326 response bytes. No images, HSC shape rows,
photo-z values, or CAS data were requested. The user separately approved
this exact plan (SHA-256
`cbf9bb94c73738b06cd5f616ba1c190f4f0e11ec71a5d914333b4b291d0e7daf`). The
two synchronous Gaia TAP batches completed and returned 58 match rows, versus
the rough estimate of 72. The response SHA-256 is
`37a8d515b415990a3444bd4654856b0fe10e2d2c2aa2e10c0c517d699eccb4b7`.

The local audit classified 38 candidates as strict diagnostic matches, 20 as
matched without clean astrometric evidence, and 400 as Gaia non-matches. The
38 strict candidates comprise 16 `i` and 22 `i2`; the target sample has 20 `i`
and 28 `i2` under the same diagnostic thresholds. The paired control-minus-
target differences in strict-candidate counts per image are −0.16 [−0.84,
0.52] for `i` and −0.25 [−0.88, 0.38] for `i2` (95% matched-frame bootstrap
intervals). Both include zero. Gaia non-matches are not evidence that HSC
candidates are non-stellar, and the thresholds do not calibrate a PSF-star
classifier.

The first two draft plans were not submitted. The final plan stored and hashed
both coordinate-specific batch queries, the synchronous Gaia TAP endpoint,
request parameters, returned fields, and candidate list. The execution guard
checked the approved plan and generated SQL before sending it. Batch SQL
SHA-256 values are `d00cb316520d992c990a72411b010283ec3dda5c7b25feea4649c2ab0590bba6`
and `fc4208c01dd5d2f3489a7f649469efb0167db3c84e001f5f0ec186db386bf852`.

Plan directory:
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_psf_control_gaia_plan_2026-11-07_03/`.
Control morphology audit:
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_psf_control_candidate_audit_2026-11-07_01/`.

Using those 38 control Gaia candidates, the local 401×401-pixel runner made
26 stamps (11 `i`, 15 `i2`); 12 were blocked by image edges. Four control
stacks used the same radial grid, local sky annuli (120–150 and 160–190 px),
and 5,000-draw frame bootstrap as the target pilot. A paired comparison first
averaged valid stars within each frame and then matched each target frame to
its control. Across all filters, `all` and center-distance-≥500-kpc subsets,
no radial bin had two valid matched frame pairs; the maximum was one. Thus no
paired control-minus-target PSF profile or uncertainty interval is estimable.
The separate target and control stacks remain diagnostics, not a calibrated
PSF or ICL correction. For controls, center distance is from the random sky
position, not a cluster center.

The target/control count comparison output is
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_psf_target_control_gaia_comparison_2026-11-07_01/`.
Control stamp output is
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_psf_control_wing_stamps_2026-11-07_01/`.
Control stack outputs are
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_psf_control_wing_stack_{i_all,i_outer500,i2_all,i2_outer500}_2026-11-07_01/`.
Paired profile-support result is
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_psf_target_control_profile_comparison_2026-11-07_03/`.

Keep `science_measurement=false` and `psf_approved=false`. The local matched
sample has since been expanded using the already received images; its Gaia
follow-up is documented below and remains unsubmitted pending a separate
approval. Coadd-epoch recovery and PSF/sky/covariance calibration remain
required before any ICL correction.

No commit or push has been made.

## 2026-11-07 additional seeing-stratified matched-pair screen

The original 49-pair comparison had at most one valid matched target/control
frame pair in any radial bin. To increase support, a local plan selected 100
new matched pairs in each actual HSC-I subfilter, `i` and `i2`, after excluding
the 49 pilot pairs. Seventeen other eligible pairs were excluded because the
target and control FITS headers reported different subfilters. Within each
filter, 100 pairs were drawn across five equal-count target-seeing strata.
The selection came from the full already-received pool of 1,043 target and
1,041 control images, all of which had passed file/FITS/WCS audits.

The approved CAMIRA member-position import supplied 7,505 member coordinates
inside these selected target images. They were used only to screen nearby
compact detections. Target centers are CAMIRA centroids and control centers
are random sky positions; neither is declared a BCG. No HSC or Gaia network
query was made in the local image preparation.

Preparation succeeded for 199 target frames and all 200 controls. Target
`S23B_11786` failed because fewer than 100 valid image pixels remained, and
three other prepared target frames produced no compact candidates. The
candidate totals are 2,069 target and 2,230 control positions. Per-filter
counts are target/control 1,016/1,000 in `i` and 1,053/1,230 in `i2`. The
candidate list remains morphology-only.

The user later approved the exact Gaia DR3 plan for 4,299 positions with SHA-256
`c00eda5f896f8160022a7e4198b3f8ce69803f0e30adc8f560ce04c8beba0e7d`. The
completed query and local audit appear in the following section. The plan
requested 2-arcsec cones and 16 Gaia fields. No HSC image, shape row, photo-z,
or CAS row was requested.

Plan and candidate audit directory
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_psf_expansion_2026-11-07_01/gaia_plan/`.
Selected-pair plan and manifests
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_psf_expansion_2026-11-07_01/plan/`.
Target and control preparation summaries
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_psf_expansion_2026-11-07_01/target_preparation/summary.json`
and
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_psf_expansion_2026-11-07_01/control_preparation/summary.json`.

The Gaia match serves as an identity and astrometry diagnostic. A Gaia match
alone cannot resolve coadd-epoch propagation, PSF-wing completeness,
sky-annulus transfer, coadd covariance, or member-preserving injection
recovery. Keep `science_measurement=false` and `psf_approved=false`.

## Matched 200-pair Gaia and PSF-wing expansion (2026-11-07)

### Gaia query and audit

The user approved the exact Gaia DR3 query for 4,299 morphology-only candidate
positions. The approved plan SHA-256 is
`c00eda5f896f8160022a7e4198b3f8ce69803f0e30adc8f560ce04c8beba0e7d`. The
query used 18 pinned batches, 2-arcsec cones, and 16 Gaia source fields. The
approved scope excluded images, CAS rows, HSC shape rows, and photo-z values.

All 18 batches completed. The aggregate response contains 706 candidate-source
rows and has SHA-256
`151f7c3f227ab7a16b1a6015563afd46fd53ccf4067d44d0264e27be75f665fb`. The
response size is 138,892 bytes, above the rough 111,955-byte estimate. A local
audit regenerated every candidate-source association from the saved Gaia
cone-union responses. The audit verified each planned SQL hash, raw-response
hash, batch match table, aggregate row, and candidate-to-image-detection link.
The 706 matches map to 706 unique candidates and 706 unique Gaia sources. No
candidate cone has multiple Gaia sources, and no Gaia source is reused for a
second HSC candidate. The maximum separation is 0.2876 arcsec.

A conservative diagnostic screen found 448 strict candidates. The screen
requires one Gaia source within 2 arcsec, RUWE below 1.4, no duplicate, QSO,
or galaxy flag, significant parallax or proper motion, and combined-star
probability of at least 0.5. The screen does not provide a calibrated
stellar or PSF classification. Gaia non-matches do not establish that an HSC
candidate is non-stellar.

| Population and filter | Morphology candidates | Gaia matches | Strict diagnostics |
| --- | ---: | ---: | ---: |
| Target `i` | 1,016 | 160 | 93 |
| Target `i2` | 1,053 | 186 | 132 |
| Control `i` | 1,000 | 129 | 70 |
| Control `i2` | 1,230 | 231 | 153 |

The remaining evidence counts are 3,593 Gaia non-matches, 255 matches without
clean astrometry, and 3 matches with clean astrometric evidence but without the
combined-star threshold.

### Stamps and radial profiles

The local CPU pass used existing approved HSC images and made 168 target and
166 control stamps with 401×401 pixels. The target pass blocked 57 candidates
at image edges. The control pass blocked 56 edge-crossing candidates and one
known variable. The stamps preserve raw image values and the fractional
SExtractor center. The profile code measures radius from the stored fractional
center without shifting stamp pixels. Detected-neighbor, CAMIRA-member, and
bad-pixel masks remain active.

The profile code subtracts a local sky estimate measured in either 120–150 or
160–190 pixels, then divides each annular mean by the masked core flux. Radial
bins extend to 96 pixels, or 16.1 arcsec. Eight separate frame-weighted stacks
cover both populations, filters `i` and `i2`, and all-star or outer-500-kpc
subsets. The bootstrap resamples matched target-control frame pairs 5,000
times. A control outer-500-kpc cut measures radius from the random sky center,
not from a cluster center.

| Filter and subset | Bins with at least two matched pairs | Maximum pairs per bin | Target valid profile rows | Control valid profile rows |
| --- | ---: | ---: | ---: | ---: |
| `i`, all | 0/44 | 1 | 480/2,992 | 504/2,508 |
| `i`, outer 500 kpc | 0/44 | 1 | 330/2,156 | 302/1,584 |
| `i2`, all | 43/44 | 14 | 855/4,400 | 1,885/4,796 |
| `i2`, outer 500 kpc | 30/44 | 8 | 598/3,168 | 1,436/3,388 |

The `i2` all-star frame-weighted means are negative from 32 to 96 pixels for
targets and controls under both sky annuli. The paired control-minus-target
mean in the 24–32-pixel bin changes sign when the sky annulus changes.

| Sky annulus | Matched pairs | Control minus target, in 10⁻⁶ core flux per pixel | 16th–84th percentile interval |
| --- | ---: | ---: | ---: |
| 120–150 pixels | 7 | −6.73 | [−13.32, 0.00] |
| 160–190 pixels | 5 | +1.67 | [−8.04, +12.02] |

Two outer-500-kpc `i2` bins have nominal 95% bootstrap intervals excluding
zero under the farther sky annulus. Each bin has only three matched pairs
among 44 radial bins. The low support and annulus-dependent sign changes do
not establish a target-control PSF-wing difference.

### Assessment and next gate

The 200-pair expansion supplies a larger empirical PSF diagnostic. The profile
measurements do not calibrate the extended PSF wing, the effective coadd PSF,
or an ICL surface-brightness correction. The ICL stack must not use the
profiles as a PSF correction. Keep `science_measurement=false` and
`psf_approved=false`.

The effective coadd epoch remains unknown, and the Gaia positions were not
propagated to an assumed epoch. The previously approved 67,071-row CCD export
covers 800 S19A lens centers and 798 controls. The prior coordinate audit found
no match within 60 arcsec between the S19A lenses and S23B CAMIRA centers. The
saved MJD rows cannot assign coadd epochs to the 200 S23B matched pairs. Any
new metadata request needs separate exact approval. Sky-annulus transfer,
coadd covariance, and injection recovery also remain open.

### Reproducibility

- Gaia expansion audit runner
  `analysis/icl_wgl_membership/hsc_debug/audit_s23b_psf_expansion_gaia_matches.py`
- Gaia match runner
  `analysis/icl_wgl_membership/hsc_debug/match_psf_gaia_batched.py`
- Stamp runner
  `analysis/icl_wgl_membership/hsc_debug/prepare_s23b_psf_wing_stamps.py`
- Frame-weighted stack runner
  `analysis/icl_wgl_membership/hsc_debug/stack_empirical_psf_star_profiles.py`
- Paired comparison runner
  `analysis/icl_wgl_membership/hsc_debug/compare_s23b_target_control_psf_profiles.py`
- Gaia audit output
  `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_psf_expansion_2026-11-07_01/gaia_audit_2026-11-07_01/`
- Target and control stamp outputs
  `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_psf_expansion_2026-11-07_01/target_wing_stamps_2026-11-07_01/`
  and
  `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_psf_expansion_2026-11-07_01/control_wing_stamps_2026-11-07_01/`
- Eight target/control filter/subset stacks
  `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_psf_expansion_2026-11-07_01/`
- Paired profile table and figures
  `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_psf_expansion_2026-11-07_01/target_control_profile_comparison_2026-11-07_01/`
