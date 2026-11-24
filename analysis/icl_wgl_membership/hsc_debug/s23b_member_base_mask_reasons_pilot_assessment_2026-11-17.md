# S23B member base-mask reasons: 16-frame pilot

## Method and checks

A deterministic selection (seed 20261107, eight eligible frames per filter)
was remeasured on the existing images. The pilot records 1,636 target and
rotated-control stamps at 50–75 proper kpc. It separates the actual base mask
into the selected quality-bit union, the bright-object bit, nonfinite image
pixels, and nonfinite or nonpositive-variance pixels. These components overlap.
Individual selected bit-plane footprints are recorded using each FITS header's
plane names. It also counts the annulus on the complete discrete pixel grid,
including coordinates outside the image, so image-boundary loss is not inferred
from continuous area alone.

For every stamp, the split components exactly reconstruct the original base
mask, and all four broad mask components reproduce the original source-mask
union. The discrete annulus counts reconcile inside/outside image bounds. All
179 retained candidate-bin flux/coverage rows reproduce the saved full-sample
measurement within absolute 1e-9. The saved header records for all 1,043
downloaded target frames have the same mask-plane positions as the current
measurement code's fixed assumptions; this verifies the sample, not arbitrary
future FITS files. No network or new download was used.

## Descriptive result

Each value below is the median of eight cluster-level fractions. Within a
cluster, pixels are pooled over its target or rotated-control stamps. Repeated
stamps from one image are correlated; the eight clusters are the descriptive
units. Percentages use in-image annulus pixels as denominator, except boundary
loss which uses the complete discrete annulus. Footprints overlap and cannot
be added. Ranges and all individual bit counts can be reproduced from the
saved table. The sample is small and selected from the existing eligible set.

| Filter / role | Stamps | Selected quality bits | Bright-object bit | `INEXACT_PSF` bit | `REJECTED` bit |
| --- | ---: | ---: | ---: | ---: | ---: |
| i target | 72 | 10.63% | 7.40% | 10.23% | 5.68% |
| i controls | 575 | 10.75% | 34.02% | 10.09% | 5.65% |
| i2 target | 110 | 8.55% | 35.28% | 7.60% | 5.23% |
| i2 controls | 879 | 7.59% | 24.77% | 6.21% | 5.02% |

The bit names are the FITS header labels; this audit does not determine the
physical origin of individual flags or whether removing one is scientifically
acceptable. Bright-object footprints varied substantially among these few
clusters. All 1,636 measured annuli had **zero discrete pixels outside the
image**, and there were zero nonfinite image or invalid-variance pixels in
those annuli. These exact zeros apply only to the pilot's measured stamps,
not all candidate positions or images. They do not imply all coverage loss
comes from the base mask: source shapes and segmentation remain additional,
overlapping components.

The bright-object and selected quality-bit footprints motivate a broader,
cluster-stratified audit, but the present output does not justify relaxing
quality cuts. The 16 frames cannot establish population-wide proportions or
a causal explanation of the empty 80% i support. ICL–satellite association,
PSF, sky, Poisson, covariance and membership gates remain unresolved;
`science_measurement=false`.

## Reproduction

Runner: `analysis/icl_wgl_membership/hsc_debug/audit_s23b_member_base_mask_reasons.py`.
Run with `--output` pointing to a new directory; `--n-per-filter 8` is the
explicit default for this deterministic pilot. The selected IDs are in its
summary. Existing pilot outputs are not overwritten.

Output:
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_member_base_mask_reasons_pilot_2026-11-17_01/`.

- Summary SHA-256: `d37955ec4b6698fce76c86b9e5221e6aefa0c2bbf21cb2962499f41f6b23f7e8`.
- Stamp table SHA-256: `adba04ffade665a0315d33ee7accc51af268426c0715ce6acac9c5049e3ceac1`.
- Runner SHA-256: `3d81324181bf3075a68ecba66b0e06ed7c7f5f6cb482dbf3ea0ba13f48823a57`.

The output table hash, exact pixel counts, mask unions and historical flux
comparison were checked. The full measurement products were not changed.
