# Visit-split covariance audit across the supported ICL sample

## Scope

This offline analysis processes the previously approved and downloaded HSC-I
warp images, masks, and variance planes. It covers the supported 94-pair ICL
sample. Four center requests from the approved download returned HTTP 404,
leaving 90 complete target/control pairs: 24 `i` and 66 `i2` (180 centers).
No additional HSC query or download was made.

## Method

Each center's expected visit IDs are sorted and alternated between two
inverse-variance half-stacks. The target and matched control use the same split
rule, but their visit-ID sets do not always match exactly. Warp pixels are
affine-mapped to the delivered coadd grid and nearest-neighbor sampled, matching
the earlier approximate warp audit rather than HSC's production resampler. The
normalized split difference is `(A-B)/sqrt(VAR_A+VAR_B)`.

The estimator measures local covariance over ±8 pixels in 64×64 patches within
the 200–300 proper-kpc annulus after a local affine-plane fit and iterative 4σ
clipping. The covariance kernels and pair counts are saved for later aperture
contraction. Four paired diagnostics are compared: normalized split RMS before
plane removal, fitted zero-lag covariance, cumulative covariance sum through
±8 pixels, and clean annulus coverage. Paired uncertainty resamples connected
5-degree sky blocks. A studentized max-|T| interval adjusts the four metrics
within each filter; the two filters are not combined into one multiplicity
family.

## Results

The split-difference RMS medians are 1.082 target versus 1.082 control in `i`,
and 1.085 versus 1.079 in `i2`. These are about 8% above unit RMS, but the
split difference includes residual sky and source terms in addition to detector
noise, so this is not a direct variance-plane calibration.

For the common complete-case paired sample (21 `i` pairs in 11 sky blocks; 64
`i2` pairs in 33 blocks), the max-T simultaneous 95% intervals for target minus
control are:

| Filter | Diagnostic | Median difference | Simultaneous 95% interval |
|---|---|---:|---:|
| `i` | normalized split RMS | +0.0006 | −0.0285 to +0.0296 |
| `i` | zero-lag covariance | +0.0049 | −0.0011 to +0.0109 |
| `i` | cumulative covariance, ±8 px | +0.1305 | −0.1726 to +0.4336 |
| `i` | clean annulus coverage | +0.0488 | −0.0806 to +0.1782 |
| `i2` | normalized split RMS | +0.0044 | −0.0095 to +0.0182 |
| `i2` | zero-lag covariance | +0.0060 | +0.0010 to +0.0110 |
| `i2` | cumulative covariance, ±8 px | +0.0069 | −0.1075 to +0.1212 |
| `i2` | clean annulus coverage | +0.0736 | +0.0063 to +0.1409 |

The i2 zero-lag and coverage intervals remain above zero within the four-metric
filter family. Their small size and correlation with unequal clean support
make them a warning that target/control noise conditions are not perfectly
matched, not a correction to the ICL signal. The median target/control visit-ID
Jaccard overlap is 0.80 in `i` and 0.845 in `i2`; some pairs have substantially
less overlap. Five centers have no clean 200–300 kpc split support, so the
covariance complete-case count is below the coverage count.

## Interpretation and limits

The paired split RMS and long-lag covariance differences do not show a resolved
target/control contrast under the within-filter simultaneous intervals. This
does not prove equal total noise. The approximate nearest-neighbor mapping,
unequal visit support, residual sky, faint sources, visit-dependent PSF and
masks, local plane fit, and non-linear clipping remain mixed. The first
aggregate pass treated empty lag cells as zero; it was superseded by this
hash-pinned rerun, which requires observed support at all 289 lags before using
the ±8 sum. The superseded output is retained for audit and must not be cited.

No covariance correction was applied to ICL or weak-lensing profiles, no
association claim is made, and `science_measurement=false` remains. The next
noise gate is to propagate heteroscedastic white-noise simulations through the
same masks, visit split, WCS mapping, and plane estimator, then compare the
predicted and observed kernels before contracting them with the exact annulus
apertures.

## Reproducibility

- Sample runner: `analysis/icl_wgl_membership/hsc_debug/analyze_icl_visit_split_covariance_sample.py`
- Sample runner SHA-256: `bc2c001d4f4715d74a808fc893d8d53227c3dc12598fbcfe57454fdb871b37f7`
- Joint bootstrap runner: `analysis/icl_wgl_membership/hsc_debug/analyze_icl_split_covariance_joint_bootstrap.py`
- Joint bootstrap runner SHA-256: `64c36f4724d71030936f7b765cc64438a0ba8e96f5e4b35c2f8afdb9727db903`
- Sample summary: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_visit_split_covariance_sample_2026-11-06_02/summary.json`
- Sample summary SHA-256: `82f7ebaa199dca4d278152969fe2a37bf5e8b4a3380bc999a74e6a276ef4f35c`
- Paired comparison table SHA-256: `b41d428bf423692da0d849b8b590bf1c106fbfd0945c93cf5c3a837466638443`
- Covariance kernel table SHA-256: `ae1aa6219138628b814f5294d0c1fb76fd68cf95b04696e68e0d3e0bae60e32f`
- Joint bootstrap summary: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_visit_split_covariance_joint_bootstrap_2026-11-06_01/summary.json`
- Joint bootstrap summary SHA-256: `973d2def87d97e8c5ef5a1527ddd2e66a76fb8bdeb0dc52ac3c6e7fca0a68f26`
- Approved parent warp plan SHA-256: `bd936a94bcdbff83a607acbfca66616c396db526ada50cc89c9e3e312d37007a`
- Approved parent download summary SHA-256: `0039321462f8158d1c006a1df64bb3e96ae988744ed780c845f74a60c79b95b4`
