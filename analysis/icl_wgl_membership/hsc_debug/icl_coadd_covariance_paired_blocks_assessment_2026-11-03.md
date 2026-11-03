# Full-sample paired coadd-covariance comparison — 2026-11-03

## Scope and method

Astrafek compared the local coadd residual-covariance measurements for 783
same-filter HSC target and matched-control pairs from the existing 800-cluster
sample. The parent audit measured a 9×9 pixel-lag grid in at most 96 unmasked
64×64 patches per frame. The paired comparison used the per-frame sum of the
variance-normalized lag grid as an approximate square top-hat variance factor.
The comparison used the target K=2.5 mask, SEP/Kron-2.5 masks for controls, and
the saved acquisition FILTER values.

Astrafek formed connected 5-degree sky blocks separately within each filter.
The i sample contains 14 blocks and 302 pairs. The i2 sample contains 45 blocks
and 481 pairs. Astrafek resampled blocks with replacement for 5,000 replicas
and calculated paired target-minus-control medians. The calculation used only
existing images and acquisition records. No HSC query or download occurred.

## Results

| Filter | Target median variance factor | Control median variance factor | Paired median difference | 95% block-bootstrap interval |
|---|---:|---:|---:|---:|
| i | 1.631 | 1.595 | +0.035 | [+0.021, +0.047] |
| i2 | 1.626 | 1.589 | +0.031 | [+0.019, +0.035] |

The approximate sigma-inflation medians are 1.277 for i targets and 1.263 for i
controls. The i2 medians are 1.275 for targets and 1.260 for controls. The
paired median differences are +0.014 for i and +0.012 for i2. The zero-lag
variance-plane ratios are nearly identical between target and control frames.
Their paired median differences are +0.00008 for i and −0.00089 for i2.

The lag-grid statistic is slightly higher in target frames in both filters.
The result supports a small difference in local pixel-correlation structure
between the selected target and control samples. The near-equal zero-lag ratios
show that the difference comes from non-zero-lag covariance rather than a
measurable shift in the variance-plane diagonal term. The bootstrap intervals
describe spatial-block stability and are not a detection test.

## Limits and status

The 9×9 finite-lag sum is not the exact covariance of the radial ICL estimator.
The per-frame kernel is a median over sampled patches, rather than a position-
dependent covariance model across the ICL aperture. The calculation omits
source Poisson noise, HSC upstream sky transfer, PSF wings, and mask-growth
response. The variance factor must not be applied directly to the ICL profile
or to the spatial-jackknife errors.

The ICL–satellite result remains exploratory and retains
`science_measurement=false`. The analysis does not establish a physical
ICL–satellite association.

## Reproducibility

- Analysis: [`analyze_icl_coadd_covariance_paired_blocks.py`](analyze_icl_coadd_covariance_paired_blocks.py)
- Parent covariance runner: [`audit_icl_full_sample_coadd_covariance.py`](audit_icl_full_sample_coadd_covariance.py)
- Parent independent audit: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_coadd_noise_correlation_full_2026-11-03_01/independent_audit.json`
- Paired-block output: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_coadd_covariance_paired_blocks_2026-11-03_01/summary.json`
- Output SHA-256: `00025ec509ece86a1ca7720c259f4ced2e3b214cdc289d0c1973e23fcfd40fa8`
- Bootstrap replicas: 5,000 per filter, seed 20261103 for i and 20271103 for i2.
- New HSC requests: none.
