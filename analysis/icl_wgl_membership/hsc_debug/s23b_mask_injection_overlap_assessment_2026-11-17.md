# Member-mask and earlier injection-recovery overlap

## Scope

This exploratory local-only audit joins the existing 96-cluster spatial mask
sample to the prior 427-frame, 0.5/1/2 nJy arcsec⁻² member-halo injections.
The overlap is 13 i and 21 i2 clusters, each with one selected injected member
and all three amplitudes. The join is unique for 102 cluster/amplitude rows;
input and output hashes were checked. No images were reread and no network
requests or data transfer were made.

The mask predictor is the difference between target and control mask fractions
pooled over **all measured members in each cluster**. It is not the mask at the
single injected member. The overlap and its injection-supported subset were
selected after the original studies and are not population weighted.

## Existing recovery outcomes

| Filter | Signal at 50 kpc (nJy arcsec⁻²) | Supported / overlapping clusters | Median recovered/input |
| --- | ---: | ---: | ---: |
| i | 0.5 | 13/13 | 0.590 |
| i | 1.0 | 12/13 | 0.518 |
| i | 2.0 | 12/13 | 0.699 |
| i2 | 0.5 | 20/21 | 0.564 |
| i2 | 1.0 | 20/21 | 0.606 |
| i2 | 2.0 | 20/21 | 0.656 |

At each unsupported row, the injected target annulus fails the coverage
criterion. These medians use only supported rows; they are conditional
measurement-transfer diagnostics, not completeness, unbiased corrections or
estimates for all members. The three amplitude samples need not be identical,
and the small overlap was selected post hoc.

At the 1.0 reference amplitude, the supported subset contains 12 i clusters
in 9 sky cells and 20 i2 clusters in 20 cells. Spearman correlations of
recovered/input with cluster-pooled control-minus-target quality-bit or
unmasked-area differences have marginal 5,000-draw sky-cell intervals that
all include zero:

| Filter | Mask difference | Spearman rho | Marginal 95% sky-cell interval |
| --- | --- | ---: | ---: |
| i | Quality bits | 0.224 | [−0.312, 0.804] |
| i | Unmasked area | 0.266 | [−0.385, 0.656] |
| i2 | Quality bits | 0.146 | [−0.404, 0.582] |
| i2 | Unmasked area | −0.060 | [−0.466, 0.403] |

These exploratory, unadjusted and noncausal intervals are broad, especially
with only nine i cells. They do not establish absence of a mask effect, and
cluster-pooled masks do not measure the injected member's exact mask.

## Remaining calibration and provenance

The prior synthetic halos were **not PSF-convolved**, did not add source
Poisson, and held the delivered coadd variance fixed. Extended PSF wings,
coadd weights, sky transfer and covariance are unresolved. This overlap
cannot establish an ICL–satellite association or be used to correct a science
profile. Keep `science_measurement=false`.

The prepared 427-position member PSF plan remains a separate explicit
user-approval request. After approval and receipt, validate the exact files,
then rerun the position-specific PSF injection; the PSF cores alone do not
supply extended wings or an effective coadd gain. Do not substitute a nearby
PSF as if it were position-specific.

- Runner: `analysis/icl_wgl_membership/hsc_debug/audit_s23b_mask_injection_overlap.py`.
- Output: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_member_mask_injection_overlap_2026-11-17_01/`.
- Summary SHA-256: `79f52d8f8b3a3052925083117635039ba442c648f0e46f29ca5a169dc9d690b6`.
- Unique 102-row ledger SHA-256: `290a200bccb22ff67f7ab6070906c3f9594e986bc71ea1df2a29167227e8d76c`.
- Both input hashes and output hash were verified, and frame/amplitude
  uniqueness was independently counted. Shared bootstrap and tied-rank
  helper behavior was verified in the preceding paired-mask analysis.

No commit or push was performed.
