# S23B member mask-area audit: predeclared spatial sample

## Selection and verification

Before reading image pixels, a deterministic plan selected 48 eligible i and
48 eligible i2 cluster images from the existing 979-frame measurement sample.
Using 5° RA/Dec cells and a fixed seed, selection visited occupied cells in a
randomized round robin: i covers all 33 occupied i cells; i2 covers 48 of 76
occupied i2 cells. The fixed 48-per-filter design and cell balancing are not
population-proportional sampling. Values below describe this selected sample;
they must not be reported as unbiased full-survey fractions.

The audit used already approved images on LagEunha CPU. It measured 9,349
member-centered or rotated-control annuli at 50–75 proper kpc. For each, the
quality/bright/invalid base components exactly reconstruct the original base;
the broader mask components exactly reconstruct the original full mask; and
the discrete full annulus equals its inside plus outside pixels. The saved
full-sample flux and coverage were reproduced to 1e-9 for all 953 retained
candidate/bin rows from the selected frames. The sample plan, audit table and
summary hashes were independently checked. No new network request, image
transfer, GPU or Slurm was used.

## Equal-cluster descriptive results

For each cluster and role, stamp pixels are pooled before calculating an
in-image fraction; the 48 resulting cluster fractions receive equal weight.
The interval is a 5,000-replicate bootstrap of the **sampled 5° sky cells**,
carrying all clusters in a selected cell together. There are 33 sampled cells
for i and 48 for i2. These intervals describe variability within the sampled
cells; they do not correct the sky-stratified selection or establish a
population-level confidence interval. Target and control stamps from a given
image are correlated and should not be treated as independent observations.
Mask components overlap, so these percentages cannot be added.

| Filter / role | Quality-bit footprint | Bright-object footprint | Neighbor-shape footprint |
| --- | ---: | ---: | ---: |
| i target | 13.05% [10.45, 16.19] | 44.53% [36.11, 52.91] | 46.21% [42.58, 50.07] |
| i controls | 22.60% [17.24, 28.52] | 43.76% [36.21, 51.41] | 41.27% [37.44, 45.13] |
| i2 target | 9.53% [7.59, 11.65] | 35.32% [28.71, 41.96] | 43.21% [39.24, 47.32] |
| i2 controls | 12.70% [9.39, 16.78] | 37.40% [31.24, 43.46] | 39.82% [36.12, 43.58] |

Brackets are marginal 95% sky-cell bootstrap intervals in percentage points.
For example, the i target quality-bit union includes `INEXACT_PSF` 12.34%
[9.74, 15.46] and `REJECTED` 7.56% [6.19, 9.22]; these flags overlap each
other and other masks. Their names are FITS header labels, not proof of the
physical source of a particular masked pixel. Exclusive-component counts and
all per-cluster measurements are in the saved tables. The figures above do
not establish a causal improvement from removing any bit or mask component.

Across these 9,349 sampled annuli, the exact discrete outside-image pixel
count and nonfinite-image count were zero. One i control stamp had 90 invalid
variance pixels. These counts concern **measured stamps only**; an all-zero
bootstrap interval is not an upper bound on missed or unmeasured populations.

## Decision and remaining work

Bright-object and neighbor-source masks cover substantial, overlapping area
in this deliberately spatially balanced sample. Quality-bit footprints are
also material. The result identifies candidate contributors to limited
coverage, but it does not distinguish the physical sources of the flags,
validate a relaxed mask, establish satellite completeness, or resolve the
ICL–outer-light association. Keep all existing masks and the PSF, source
Poisson, sky transfer, coadd covariance and membership gates;
`science_measurement=false`.

A useful next analysis is a **paired cluster-level comparison** of target and
rotated-control footprints, with the same sky-cell bootstrap draw per pair,
then comparison against coverage and relative ICL quantities. Such selection
associations are exploratory and require the same sampling limitations.

## Reproducibility and outputs

- Predeclared plan: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_member_mask_spatial_sample_2026-11-17_01/plan.json`, SHA-256 `c32033e8108fa2d09b419e4aec785793882f5145a029ef3efa9f5af1a49d69ed`.
- Image audit root: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_member_mask_spatial_audit_2026-11-17_01/`, summary SHA-256 `3e66d7d026f121322a8b4a76722dd0dfe7f09b3a9c7765a954d7a74c2845fe85`.
- Cluster/sky-cell summary root: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_member_mask_spatial_summary_2026-11-17_01/`, summary SHA-256 `b80bfbbe0b0982566be98bada466b32680b4395189fa46d2b49955894e5f5138`.
- Runners: `prepare_s23b_mask_spatial_sample.py`, `audit_s23b_member_mask_spatial_sample.py`, `summarize_s23b_member_mask_spatial_sample.py` under `analysis/icl_wgl_membership/hsc_debug/`.
- A synthetic two-cell check confirmed that the bootstrap resamples all
  clusters of each drawn cell together. The stamp-table, sample-set and
  output hashes were verified, as were all exact reconstruction gates above.

No commit or push occurred during this local analysis. The S16A weak-lensing
random-center count query is a separate, approval-gated operation.
