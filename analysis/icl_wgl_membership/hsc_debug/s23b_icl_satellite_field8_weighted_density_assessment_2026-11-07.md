# S23B ICL–satellite raw-field-8 weighted-density sensitivity

Date: 2026-11-07
Status: exploratory sensitivity; `science_measurement=false`.

## Scope and field interpretation

This offline analysis uses the already received, user-approved
`camira_s23b_wide_sm_member.dat.gz` file, its pinned cluster/WCS join, the
saved S23B target/control radial profiles, and their sky-block assignments.
It makes no network or CAS request, downloads no data, and reads no image
pixels. The compressed member file is 10,626,590 bytes, SHA-256
`66529e02c6ae11ebdccc3fb1d4756b639364f372f9d0f09a721786a19906945a`.

The member file has no header. Published CAMIRA sources describe a member
catalog containing stellar mass and a membership weight; the S23B paper
counts members above `w_mem >= 10^-3`, and Oguri (2014) defines `w_mem` from
the galaxy likelihood, stellar-mass filter, and spatial filter. The imported
eighth field lies in `[0.001000069, 1.0]`, which strongly supports it as
`w_mem`; the seventh field's numerical range is consistent with log stellar
mass. However, the public S23B catalog README and paper do not state the
column order for this specific headerless file. The analysis therefore keeps
the field named `member_field8_raw` and treats its weighting as a sensitivity,
not a verified membership-probability measurement.

References:

- [Oguri et al. 2026, S23B CAMIRA catalogs](https://arxiv.org/html/2512.13954v1), especially Sections 2.3 and 3.1.
- [Oguri 2014, CAMIRA definition](https://academic.oup.com/mnras/article/444/1/147/1017165), Equation 20.
- [Toba et al. 2023, CAMIRA member-catalog use](https://academic.oup.com/pasj/article/75/6/1246/7286233), Sections 2.1 and 4.1.1.
- [Official S23B catalog README](https://raw.githubusercontent.com/oguri/cluster_catalogs/main/hsc_s23b_camira/README.md), which does not document the member-file columns.

## Measurement

The predictor is the log ratio of the sums of raw field 8 for in-image
candidate members in 50–300 and 300–800 proper-kpc annuli, divided by their
geometric areas. A 0.5 pseudo-count is added to each annular sum. The outcome
is the existing common-pixel target-minus-control ICL contrast between those
same annuli. The regression retains the existing adjustments for log richness,
cluster redshift, target/control SFD mismatch, ICL common-pixel support, and
member image-footprint fraction. The weighted and unweighted fits use exactly
the same targets and paired spatial-block bootstrap draws (5,000 per layout).

The primary layout is the 5-degree, phase-zero grid:

| Filter | Pairs / connected blocks | Raw-field-8 weighted slope, 95% block interval | Unweighted slope, 95% block interval | Paired slope change, 95% interval |
| --- | ---: | ---: | ---: | ---: |
| HSC-i | 232 / 26 | 5.13 [−0.80, 13.05] | 5.72 [−0.06, 13.63] | −0.59 [−3.26, 1.60] |
| HSC-i2 | 756 / 42 | 5.01 [2.44, 7.62] | 4.24 [0.99, 7.80] | +0.77 [−1.05, 2.46] |

The fitted change from field-8 weighting is compatible with zero in all six
sky-grid layouts for both filters. In the full supported sample, the i2
weighted slope interval excludes zero in all six layouts, while i is
layout-dependent. When requiring at least 10% common image pixels in both ICL
annuli, the primary i2 weighted interval becomes [−0.21, 6.57]; only two of
six layouts exclude zero. The i interval includes zero in every layout. This
support sensitivity prevents treating the full-sample i2 interval as robust.

## Assessment

Raw field-8 weighting does not materially change the fitted relation relative
to unweighted candidate counts: the paired slope differences are small
compared with their block-bootstrap intervals. A positive i2 point estimate
persists in the full sample, but the interval crosses zero under the primary
layout after the 10% support requirement, and the i filter is unresolved.
This is not a robust ICL–satellite detection. Photometric-member purity and
completeness, source-detection selection, the assumed column mapping, the
outer reference light, residual sky, extended PSF wings, source Poisson noise,
and coadd covariance remain unresolved. Keep `science_measurement=false`.

The result does not change the weak-lensing sample conclusion: the existing
S19A lens sample and S23B member-profile sample have no positional overlap.
The separate S19A count-only plan for 427 S23B profile targets remains
unsubmitted; a member-catalog download approval does not submit that CAS
query. A same-cluster three-way weak-lensing/ICL/satellite measurement still
requires that independently approved count query and later validation of
matched random support before any source-shape export.

## Reproducibility

Runner:
`analysis/icl_wgl_membership/hsc_debug/analyze_s23b_icl_satellite_field8_weighted_density.py`
(SHA-256 `e33390cfaec625ff0c70c9f2b56497c1809ec0083af8e2826a1195f6f96a87e3`).

Output directory:
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_icl_satellite_field8_weighted_density_2026-11-07_02/`

- Summary SHA-256: `174ec016278ca69dcdbb5497ac9da79449911033f91a3c73b03e5c1e2bc33ed3`
- Paired measurements: `paired_field8_weighted_density.tsv`
- Six-layout comparison: `paired_layout_sensitivity.tsv`
- Support-floor sensitivity: `common_support_floor_sensitivity.tsv`
- Plot: `field8_weighted_density_relation.png`

No tests, CAS submissions, or image downloads were run for this analysis.
