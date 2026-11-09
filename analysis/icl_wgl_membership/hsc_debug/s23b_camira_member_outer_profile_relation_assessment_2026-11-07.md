# S23B CAMIRA member outer profiles and ICL relation assessment

Date: 2026-11-07
Status: exploratory methods assessment; `science_measurement=false` and
`images_ready_for_science=false`.

## Scope and data use

This assessment uses the already approved S23B CAMIRA member catalog, the
already approved HSC-I coadd images and masks, saved target/control ICL
profiles, and existing metadata. No files were downloaded and no network
requests were made in these runs. The CAMIRA member file was not fetched a
second time.

The source-aware stack selected all 979 eligible target/filter frames with a
saved ICL contrast, a successful image mask, and at least four in-frame CAMIRA
positions in the 50–500 proper-kpc selection range: 226 HSC-i and 753 HSC-i2.
There were 16,503 CAMIRA candidate positions in that range. The stack retained
13,269 candidate-by-radial-bin measurements after source matching and paired
coverage cuts. A 50–75 kpc cluster-level outcome was available for 62 HSC-i
and 365 HSC-i2 frames.

## Measurement

CAMIRA candidates are photometric member candidates, not spectroscopically
confirmed satellites. For each candidate, the pipeline attempts a unique SEP
detection match within 1 arcsec and preserves that detection's segment while
masking unrelated detections, bad/bright pixels, and fallback circles around
unmatched candidates. It estimates a local baseline from up to eight
same-radius azimuthal rotations, rejecting controls close to detections or
other CAMIRA positions. A candidate profile bin is retained only when both
target and control annuli have at least 60% coverage and at least three valid
rotated controls.

The radial bins are 5–10, 10–15, 15–20, 20–30, 30–40, 40–50, and 50–75
proper kpc. The median cluster-level candidate-minus-control surface
brightness decreases strongly with radius. In the outer 50–75 kpc bin it is
0.68 nJy arcsec⁻² for HSC-i (cluster q16–q84: −3.77 to 5.64; 62 clusters) and
1.46 nJy arcsec⁻² for HSC-i2 (−4.23 to 8.94; 365 clusters). These are
descriptive cluster distributions, not confidence intervals or detections.

## ICL–outer-member-light relation

The association model relates the paired target-minus-control ICL contrast
(50–300 minus 300–800 proper kpc) to the cluster mean candidate-minus-rotated-
control light at 50–75 proper kpc. It adjusts for log richness, redshift,
minimum target/control common-pixel support, and median target/control member
profile coverage. Uncertainty is estimated from 5,000 connected spatial-block
bootstrap replicates for each of six RA-Dec block layouts (2.5, 5, and 10
degrees, with two grid phases).

| Filter | Clusters / connected blocks (5° phase 0) | Adjusted slope | 95% spatial-block interval | Layouts excluding zero |
| --- | ---: | ---: | ---: | ---: |
| HSC-i | 62 / 19 | 0.086 | [−0.098, 0.334] | 0/6 |
| HSC-i2 | 365 / 43 | −0.041 | [−0.145, 0.090] | 0/6 |

The slope is the change in outer candidate light per 1 nJy arcsec⁻² change in
the relative ICL contrast; its numerical units cancel. Every primary-layout
interval includes zero. Leave-one-block-out estimates change sign in 0 HSC-i
and 1 HSC-i2 deletion. The association is therefore unresolved in these data;
the results do not support a positive or negative physical claim.

The conclusion is stable to the examined ICL common-support floors in the
sense that no available block-bootstrap interval excludes zero. At floors of
5%, 10%, and 20%, the HSC-i samples are 56, 47, and 23 clusters; the HSC-i2
samples are 348, 302, and 155. HSC-i has only three valid layouts at the 20%
floor; the other layouts are insufficient in connected blocks. The fitted
point-estimate sign changes as the support floor changes, further cautioning
against interpreting a single fit as a detection.

## Limitations and decision

This is a source-aware exploratory profile, not a calibrated measurement of
diffuse satellite structure. Important remaining effects include uncertain
CAMIRA member purity/completeness and unverified member-file fields, residual
sky transfer, extended PSF wings, Galactic extinction, source Poisson noise,
coadd pixel covariance, centering and mask selection. Rotated controls share
the cluster image and can contain cluster light. The outer ICL reference
annulus can also contain cluster light. SEP matching and the coverage cut
select a subset of the CAMIRA candidates, and the sample is much smaller in
HSC-i than HSC-i2 at 50–75 kpc. The analysis does not model PSF wings or
recover missed/faint candidates.

**Decision:** the profiles are useful for method development and for planning
injection/recovery and mask/covariance work. They are not ready for a science
claim or a manuscript result. Keep `science_measurement=false` and
`images_ready_for_science=false`.

Recommended next gates:

1. Verify the CAMIRA member-file column schema and membership-weight meaning
   against an authoritative catalog description.
2. Calibrate member detection and masking with injections that preserve each
   target member while handling its PSF and unrelated sources.
3. Measure residual-sky transfer, extended-PSF contamination, source Poisson
   noise, and coadd spatial covariance at these radii.
4. Re-evaluate centering, mask/coverage selection, and matched-control
   construction; then repeat the block-level association analysis.

## Reproducibility

Analysis runners:

- `analyze_s23b_member_outer_icl_relation.py` — spatial-block adjusted relation.
- `measure_s23b_camira_member_outer_profiles.py` — full eligible source-aware
  profile stack.

The adjusted relation summary is outside Git at
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_camira_member_outer_icl_relation_2026-11-07_01/summary.json`
(SHA-256 `584082c4907be8f9d37cd361a5bdb95887cc5163a4cf18b928fa86786149a884`).
The full profile summary is at
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_camira_member_outer_profiles_full_2026-11-07_01/summary.json`
(SHA-256 `e6b91fd7fe07caf89a6012053a46bba38be17fdb5298eb265b931ccc5e0c4cdd`).
The relation runner SHA-256 is
`ce0c56c0d234f805071880144e70a3ac5238d93c4664e9998374ceb20a02bc8f`; the
profile runner SHA-256 is
`cc91779a92b3e6cfa86a564d9b9072f17f519292f95fca17c7186f1cbe3f4c87`.

Figures:

- Full radial profiles and raw descriptive relation:
  `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_camira_member_outer_profiles_full_2026-11-07_01/member_outer_profiles_and_icl_relation.png`
  (SHA-256 `32243dbb514ea0f61c4594e78ae802345636c7ad5aa5ef82edf671d0fcd9fbb1`).
- Adjusted spatial-block relation:
  `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_camira_member_outer_icl_relation_2026-11-07_01/adjusted_outer_relation.png`
  (SHA-256 `634d22ff128abe041482d7d7cac4f13f6506fd69a0d18d182ecbbbb31620c1fa`).

No commit or push has been made.
