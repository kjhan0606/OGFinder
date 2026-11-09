# S23B ICL–satellite candidate relation: raw member-field-8 sensitivity — 2026-11-07

## Purpose and scope

This analysis asks whether the exploratory relation between S23B relative ICL
contrast and CAMIRA candidate-member density changes when each row is weighted
by the raw eighth numeric field in the headerless member file. The file schema
has not been mapped well enough to identify that field as CAMIRA membership
probability, so this is explicitly a raw-field sensitivity analysis. No claim
that the field is w_mem is made.

The response is the common-support target-minus-control surface-brightness
contrast at 50–300 proper kpc minus that at 300–800 proper kpc. The predictor
is the base-10 logarithm of the ratio of inner to outer surface densities,
where each candidate contributes its raw eighth-field value and a +0.5
continuity correction is added in raw-field units. The model adjusts for log
richness, redshift, target-control SFD E(B−V) mismatch, and minimum common
pixel support. Uncertainty is estimated from 5,000 connected RA-Dec block
bootstrap draws for each of six spatial layouts (2.5, 5, and 10 degree
blocks, with two phases).

This is an offline calculation on saved profile/member-position tables. It
made no HSC requests and read no science-image pixels. The supported samples
contain 232 HSC-i and 756 HSC-i2 target/control pairs.

## Raw-field audit

All 45,634 target-member rows have finite eighth-field values in [0, 1]:
minimum 0.001001, median 0.373253, maximum 1.000000. Per-target sums of this
field divided by catalog richness have q10/median/q90 of 0.830/0.945/1.071.
Those properties are consistent with a probability-like field, but they do
not establish the column mapping or calibration. The S23B paper discusses a
CAMIRA membership weight, while the published catalog README does not document
the local headerless member-file columns; neither source by itself verifies
that this particular raw field is that weight.

References checked:

- Oguri et al., S23B CAMIRA cluster catalog paper: <https://arxiv.org/html/2512.13954>
- Public S23B catalog repository README: <https://raw.githubusercontent.com/oguri/cluster_catalogs/main/hsc_s23b_camira/README.md>

## Results

At the primary 5-degree, phase-zero layout, the adjusted slope per dex of
field-8-weighted candidate-density ratio is:

| Filter | Pairs / connected blocks | Slope (nJy arcsec⁻² dex⁻¹) | 95% connected-block interval |
| --- | ---: | ---: | ---: |
| HSC-i | 232 / 26 | 5.13 | [−0.80, 13.05] |
| HSC-i2 | 756 / 42 | 5.01 | [2.44, 7.62] |

The i2 interval excludes zero in all six block layouts; the i interval does so
in two of six. This is spatial-partition stability for the permissive
common-support sample, not evidence that the underlying relation has been
measured robustly. The intervals remain wide, especially in i.

The common-support floor sensitivity is less favorable:

| Minimum common support in each annulus | Filter | Pairs | Slope | 95% interval |
| ---: | --- | ---: | ---: | ---: |
| 5% | HSC-i | 169 | 3.38 | [−3.42, 12.84] |
| 10% | HSC-i | 103 | 3.93 | [−1.90, 10.33] |
| 20% | HSC-i | 28 | 17.47 | [−17.96, 33.39] |
| 5% | HSC-i2 | 659 | 4.68 | [1.98, 8.11] |
| 10% | HSC-i2 | 525 | 3.29 | [−0.24, 6.76] |
| 20% | HSC-i2 | 197 | 3.54 | [−0.51, 7.56] |

At the 10% support floor, both filters' intervals include zero. The 20% i
estimate is based on only 28 pairs and is highly imprecise. These changes show
that the apparent i2 precision depends on retaining low-common-support pairs.

## Assessment

The weighted sensitivity does not change the cautious interpretation of the
unweighted analysis. It suggests a positive i2 trend under the most permissive
support cut and across all six spatial layouts. The trend is not stable to the
10% support requirement, and the raw field has not been identified. HSC-i
remains unresolved. Large target-level residuals also remain, so a small
number of unusual clusters can matter even when leave-one-block-out slopes
stay positive.

Do not interpret this as a detection, a causal effect, or a calibrated
membership-weight result. The 300–800 kpc reference annulus may contain
cluster light; member candidates are photometric; support is sparse; and the
pixel-level PSF wings, sky transfer, mask selection, and coadd covariance are
not included in this regression uncertainty. Keep
science_measurement=false and images_ready_for_science=false.

## Next analysis gates

1. Resolve the raw member-file column definitions from a primary data product,
   producer documentation, or an independently verified catalog/code release;
   until then, keep field 8 labeled as raw.
2. Audit member positions against usable image pixels and masks, not just FITS
   footprint geometry.
3. Measure member-centered outer-light profiles against matched control
   positions and quantify selection with injection/recovery.
4. Propagate residual-sky, PSF-wing, and spatially correlated coadd-noise
   uncertainty before evaluating a physical ICL–satellite relation.

## Reproducibility

Runner:
analysis/icl_wgl_membership/hsc_debug/analyze_s23b_icl_satellite_field8_weight_sensitivity.py
(SHA-256 `6ce39b603f5fdbd75d616200c2cc24a2f302b5dd304e42c1f7939152bfa81fd2`).

Output directory:
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_icl_satellite_field8_weight_sensitivity_2026-11-07_01/`

Summary SHA-256:
`ec2f6a325ac0a6e349b6f0f52cd106c544d3d811594c35568f5712e3273b2500`.

Additional product hashes:

- paired relation TSV: `f06a35f744339b655786fe9375ccdfe776bcdbe8cf4bd2c5b0a476fc66b7fce`
- six-layout sensitivity TSV: `ce16e32a2cea8f91a15b090f3c292d6a5516a08e2d59cd322e22ccb8b4909e15`
- leave-one-block-out TSV: `09dacbc9fca170a83a7a2a3e82c239160e59395e06483140150170c1f7190481`
- support-floor sensitivity TSV: `ffe12d91576fd26ac8fefbe4ddbf90b009e7eaa529388e0a204b3b410b68cac3`
- plot: `fd804bb1a43cb7c7ad830b76ded365491af4d4cdafefb5c7ee15a4b049eb7f53`

All products are exploratory and retain `science_measurement=false`.

## Follow-up on field-8 identity (2026-11-07)

An offline audit of all 637,146 rows subsequently found that raw field 8 ranges
from 0.001000069 to 1.0 and every row satisfies the published S23B CAMIRA
member-catalog threshold `w_mem >= 10^-3`. The rows map to all 12,162 Wide
star-mask clusters. This is strong evidence that field 8 is probably
`w_mem`. The public README still does not document the headerless member-file
column order, so the analysis products remain labeled as raw field 8 and the
weighted relation remains a sensitivity result.

The published richness is not the unqualified sum of those weights: CAMIRA
applies stellar-mass and spatial filters and a compensated background filter.
Thus the imperfect field-8-sum/richness match does not disprove the likely
`w_mem` mapping. See [Oguri et al. 2026](https://arxiv.org/html/2512.13954v1)
and the [official S23B catalog README](https://raw.githubusercontent.com/oguri/cluster_catalogs/main/hsc_s23b_camira/README.md).

The schema-consistency assessment and all-cluster audit products are recorded
in `s23b_camira_member_raw_field_schema_assessment_2026-11-07.md` and
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_camira_member_schema_audit_2026-11-07_04/`.
