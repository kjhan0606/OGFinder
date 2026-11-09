# S23B CAMIRA member outer-profile raw-field-8 sensitivity

Date: 2026-11-07

Status: exploratory sensitivity only; `science_measurement=false` and
`images_ready_for_science=false`.

## Why this check was run

The source-aware profile stack treated each retained photometric CAMIRA member
candidate equally. The imported member file has an eighth raw field with
values in [0.0010, 1.0], but its column name and meaning are not documented in
the local headerless file. The S23B catalog paper describes the richness as a
membership-weighted quantity and refers to a member weight factor
`w_mem`; this establishes that such a CAMIRA quantity exists, but does not
map this file's eighth field to it. The public S23B repository README lists
cluster catalogs and links to separate photometric-LRG data, but it does not
provide the downloaded member-file column layout. See [Oguri et al. (2026)](https://arxiv.org/abs/2512.13954)
and the [official S23B catalog README](https://raw.githubusercontent.com/oguri/cluster_catalogs/main/hsc_s23b_camira/README.md).

For that reason this analysis labels the quantity only as **raw field 8**.
It does not call it `w_mem`, a membership probability, or use it to promote a
science measurement.

## Calculation

The runner joined the imported raw field to every saved candidate profile bin
by target ID and source row. All 13,269 profile-bin rows joined uniquely; the
input included 45,634 member-position records. For each target/filter/radial
bin it computed the raw-field-8 weighted mean of the already-saved
candidate-minus-rotated-control profiles. The unweighted mean was reproduced
for every target/filter outer bin before fitting the sensitivity.

The ICL relation reuses the original 50–75 proper-kpc outcome and the same
covariates: log richness, redshift, minimum ICL common support, and median
target/control profile coverage. Each layout uses 5,000 connected target and
control sky-block draws. The weighted and unweighted fits use the same draws,
allowing a paired interval for their slope difference. No science-image pixels
were read and no data was downloaded in this analysis.

## Results

Raw field 8 has median 0.373, minimum 0.0010, and maximum 1.0 in the imported
positions. Across the 427 target/filter clusters with an outer profile, the
weighted 50–75 kpc candidate-minus-control median is 1.03 nJy arcsec⁻² in
HSC-i (62 clusters; q16–q84 −3.87 to 6.36) and 1.99 nJy arcsec⁻² in HSC-i2
(365 clusters; −4.43 to 10.89). These descriptive cluster spreads still
include zero and are not confidence intervals.

| Filter | Clusters / blocks (5° phase 0) | Unweighted slope (95% block interval) | Raw-field-8 weighted slope (95% block interval) | Paired weighted-minus-unweighted slope (95% interval) |
| --- | ---: | ---: | ---: | ---: |
| HSC-i | 62 / 19 | 0.086 [−0.098, 0.334] | 0.017 [−0.247, 0.334] | −0.069 [−0.184, 0.057] |
| HSC-i2 | 365 / 43 | −0.041 [−0.145, 0.090] | −0.013 [−0.118, 0.129] | 0.028 [0.0003, 0.080] |

The weighted association interval includes zero in all six spatial layouts
for both filters. The paired HSC-i2 slope change is marginally positive in
3/6 layouts and includes zero in the other 3/6; HSC-i changes include zero in
all 6/6 layouts. No familywise adjustment is applied to this exploratory
layout sensitivity. At 5%, 10%, and 20% ICL-support floors, every available
weighted slope interval still includes zero. HSC-i has only three valid
connected-block layouts at the 20% floor.

**Interpretation:** weighting by this undocumented field changes the
HSC-i2 coefficient slightly toward zero in the primary layout, but the
weighted relation itself remains unresolved, and the change is spatial-layout
sensitive. This does not establish that field 8 is a membership weight and
does not support a physical ICL–satellite association.

## Limits and decision

The field's exact semantics remain unverified. This weighting does not correct
candidate incompleteness or contamination, SEP matching selection, annular
coverage selection, residual sky transfer, extended PSF wings, Galactic
extinction, source Poisson noise, coadd covariance, or possible cluster light
in rotated controls and the outer ICL reference. Keep
`science_measurement=false` and `images_ready_for_science=false`.

Next gates are to obtain an authoritative file-specific column definition,
calibrate member detection and profile recovery with member-preserving
injections, and repeat the relation after PSF/sky/covariance and center/mask
sensitivities. Do not relabel raw field 8 before the file layout is verified.

## Reproducibility

Runner: `analyze_s23b_member_outer_raw_field8_sensitivity.py`
Runner SHA-256: `5a2608904f9631965e9fb712adbb561d073ef86966b2c7ec0c0afe6eb5336330`.

Output directory (outside Git):
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_camira_member_outer_field8_sensitivity_2026-11-07_01/`

- `summary.json` SHA-256: `7156397a478932a9d1fb133b88a758a814bfc0543708719cd36f2a21c1b14875`
- `raw_field8_weight_sensitivity.png` SHA-256: `13dae58995cbbac004462b5fabedc0512c26ab7436b701a08050f4e112e3d717`
- `raw_field8_relation_layout_sensitivity.tsv` SHA-256: `874093ad09fd22eebce1e2dc3df71a2cbf47e0abde07b914f33aecef75e6bf98`
- `raw_field8_support_floor_sensitivity.tsv` SHA-256: `eee2c2a25dd4608e2ed8fde1f941f8bf42c35b4679babf33954d7ed37935cd9d`
- `raw_field8_weighted_radial_profiles.tsv` SHA-256: `79f77a94b639e75c721164bd5369b08a182c922a5c12e3619530fe064790e78e`

No commit or push has been made.

## Follow-up on field-8 identity (2026-11-07)

The subsequent all-catalog audit found 637,146 member rows with raw field 8
entirely in [0.001000069, 1.0], all meeting the published `w_mem >= 10^-3`
member-catalog threshold, across all 12,162 S23B Wide star-mask clusters.
This strongly supports field 8 as probable `w_mem`, although the public README
does not document the headerless member-file column ordering. The sensitivity
above remains exploratory; its intervals still include zero in every layout
and support-floor fit. Keep the stored field name as `field8_raw` pending
file-specific schema confirmation.

See [Oguri et al. 2026](https://arxiv.org/html/2512.13954v1) and the [official
S23B catalog README](https://raw.githubusercontent.com/oguri/cluster_catalogs/main/hsc_s23b_camira/README.md).
The independent schema assessment is
`s23b_camira_member_raw_field_schema_assessment_2026-11-07.md`.
