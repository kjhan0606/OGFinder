# S23B ICL–CAMIRA member-position relation — 2026-11-07

## Member-file receipt and validation

The user approved the single S23B CAMIRA Wide star-mask member file, file ID
1DEnjGj9P1ilUooEm2U5_A-qXUrxeiJdE, at 10,626,590 compressed bytes. No other
file was downloaded, and no credentials were saved.

The gzip stream is intact. The headerless table has 637,146 rows and eight
numeric fields. Its first four fields match the pinned S23B cluster catalog
(cluster RA, Dec, richness, and photometric redshift) for all 12,162 distinct
cluster groups, with no unmatched or ambiguous group. All 1,043 target IDs
join uniquely and provide 45,634 member-position rows (16/39/151 minimum,
median, and maximum rows per target). Fields 5 and 6 are retained as member
coordinates and checked against the target FITS WCS. Trailing fields 7 and 8
are retained as raw values; this import does not assign them meanings or use
them as weights.

The WCS footprint check places 39,166 of the 45,634 catalog positions inside
the 1,043 target images (85.83%). The remaining 6,468 positions fall outside
the downloaded square cutouts. In the profile-supported relation sample,
all CAMIRA positions in both the 50–300 and 300–800 kpc annuli are within their
target image pixel-center bounds. This geometric result does not establish
that those positions land on usable, unmasked pixels.

Import runner:
analysis/icl_wgl_membership/hsc_debug/import_s23b_camira_member_catalog.py
(SHA-256
46e4d4d00a0e828f23a90d644e3e8d1fd522b69ab59fddf3d1f920d031089817).
Import summary:
/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_camira_member_catalog_2026-11-07_01/import_2026-11-07_02/summary.json
(SHA-256
3d1442420aeae051f5a970c6b717a93c69fdaa57cebb096063c46d8afe5f4d4e).
Compressed member-file SHA-256:
66529e02c6ae11ebdccc3fb1d4756b639364f372f9d0f09a721786a19906945a.

## Direct exploratory association

The analysis uses the 1,022 filter-matched target/control image pairs and
joins catalog positions to the saved pair radial profiles. The outcome is the
common-support target-minus-control mean at 50–300 kpc minus that at
300–800 kpc. The predictor is the base-10 logarithm of the CAMIRA candidate
member surface-density ratio between those annuli, with a 0.5-count
continuity correction. Candidate counts are unweighted. The fit adjusts for
log CAMIRA richness, cluster redshift, target-minus-control SFD E(B−V)
mismatch, and the minimum common-pixel fraction in the two annuli. The
position-footprint fraction is constant at 100% in the supported annuli and
therefore is not a fitted covariate.

The saved annular estimates were reproduced within 1e-9 for all 988
supported pair/filter records. Only 232 i and 756 i2 pairs have common
pixels in both ICL annuli. Median common-pixel fractions are only 15.0% and
11.6% (i, inner and outer), and 20.1% and 16.8% (i2). Median in-image
candidate-member counts are 8/13 (i) and 9/15 (i2), inner/outer.

At the primary 5° phase-zero connected-block layout, adjusted slopes are:

| Filter | Pairs / blocks | Slope (nJy arcsec⁻² dex⁻¹) | 95% block interval |
|---|---:|---:|---:|
| HSC-i | 232 / 26 | 5.72 | [−0.06, 13.63] |
| HSC-i2 | 756 / 42 | 4.24 | [0.99, 7.80] |

The point estimates are positive, and deleting any single 5° block leaves
both filter slopes positive (HSC-i 4.52–8.04; HSC-i2 3.50–5.16). The intervals
are less stable across the six block-grid layouts: they exclude zero in 3/6
layouts for HSC-i and 4/6 for HSC-i2.

The support-floor sensitivity is more restrictive:

| Minimum common support in both annuli | Filter | Pairs / blocks | Slope | 95% block interval |
|---:|---|---:|---:|---:|
| 0% (any positive support) | HSC-i | 232 / 26 | 5.72 | [−0.06, 13.63] |
| 0% (any positive support) | HSC-i2 | 756 / 42 | 4.24 | [0.99, 7.80] |
| 5% | HSC-i | 169 / 23 | 4.14 | [−3.47, 13.75] |
| 5% | HSC-i2 | 659 / 44 | 4.49 | [0.94, 7.97] |
| 10% | HSC-i | 103 / 17 | 3.84 | [−2.07, 9.70] |
| 10% | HSC-i2 | 525 / 44 | 3.33 | [−0.91, 7.32] |
| 20% | HSC-i | 28 / 10 | 17.61 | [−8.46, 31.83] |
| 20% | HSC-i2 | 197 / 39 | 4.21 | [−0.13, 8.27] |

The support threshold sensitivity is at 5° phase zero. It shows that the
HSC-i2 interval includes zero once both annuli require at least 10% common
support; the HSC-i result remains unresolved throughout. The 20% HSC-i
estimate is especially imprecise at only 28 pairs.

## Interpretation and remaining gates

This is the first direct S23B association with CAMIRA member positions, but it
is exploratory, not a detection. The positive HSC-i2 slope depends on a
permissive support sample and is not stable under stricter common-area
requirements. Sparse common-pixel coverage, large cluster-level residuals,
photometric-member contamination/completeness, and unassigned member-file
fields limit the inference. The 300–800 kpc reference annulus may contain
cluster light, so the outcome is a relative ICL contrast rather than an
absolute ICL surface brightness.

Residual sky transfer, extinction correction, extended PSF wings, complete
coadd covariance, center/mask sensitivity, and member photometric completeness
remain unresolved. Do not promote either the imported catalog or this
association to science_measurement=true.

## Reproduction and outputs

Runner:
analysis/icl_wgl_membership/hsc_debug/analyze_s23b_icl_satellite_relation.py
(SHA-256
d961b373e52289d374c203c3c5a1887cf731c21b6dc8d6d562ecb6a5af3bfd53).
It reads saved profile tables and catalog positions, uses six connected
RA-Dec block layouts with 5,000 bootstrap draws each, and makes no network
requests or image-pixel reads.

Output directory:
/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_icl_satellite_relation_2026-11-07_03/
(summary SHA-256
fab29d7aef464707c31e2616fe1db9fc3bad9d281ad6d5dfa6d63f2b484bbb39).
It contains the pair-level joined table, six-layout slope sensitivity,
5° leave-one-block-out diagnostics, common-support floor sensitivity, and
the adjusted relation plot. All products retain science_measurement=false.
