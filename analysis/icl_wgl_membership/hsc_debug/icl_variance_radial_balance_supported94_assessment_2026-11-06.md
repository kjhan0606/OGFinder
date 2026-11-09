# Radial coadd-variance balance in the 94 supported ICL–satellite pairs (2026-11-06)

## Question and method

This offline analysis asks whether the delivered HSC coadd VARIANCE plane
differs radially between each target and its matched control, and whether the
target/control variance ratio tracks the saved ICL contrast. It reads only
the existing 94 target/control FITS pairs (26 `i`, 68 `i2`), verifies each
file against its download SHA-256, checks the saved HSC IMAGE/MASK/VARIANCE
mapping and center WCS, and uses the matched-pair redshift to convert angular
distance to proper kpc.

For ten annuli from 0 to 900 kpc, the statistic is the median delivered
variance after converting native units to nJy² arcsec⁻⁴. Non-finite and
non-positive variance pixels and reviewed bad/bright-star mask bits are
excluded. No source-detection masks are applied because this measures the
instrumental variance plane; valid-pixel fractions are reported. A lens-bin
median is omitted if fewer than 50 valid pixels remain. Target/control
differences and the partial correlation with ICL are summarized with 5,000
resamples of connected 5-degree target/control sky blocks per filter and
annulus. The ICL correlation controls redshift and log richness.

## Results

The median target/control variance ratios are close to one across radius. Three
of the 20 unadjusted pointwise 95% block-bootstrap intervals for the paired
log-ratio exclude zero:

| Filter | Radius (kpc) | Median target/control variance ratio | 95% interval for ratio |
|---|---:|---:|---:|
| `i` | 0–50 | 1.058 | 1.009–1.265 |
| `i2` | 50–75 | 0.974 | 0.931–0.995 |
| `i2` | 75–100 | 0.952 | 0.923–0.993 |

The remaining 17 pointwise intervals include equality. These are 20 radial
comparisons without a multiplicity correction, from a selected support sample;
they are weak, filter-dependent diagnostics, not evidence for a calibrated
noise difference. In all 20 annuli, the ICL partial-correlation bootstrap
interval includes zero. Its point estimates range from −0.065 to +0.209 in
`i` and −0.118 to +0.140 in `i2`; this selected sample does not show a clear
trend between ICL contrast and target/control variance ratio.

Median valid-pixel fractions decline toward the outskirts: about 0.53 in
600–900 kpc for both targets and controls, with target/control support
differences varying by filter and radius. Thus outer-bin medians describe the
available variance pixels rather than complete annular coverage. The 0–50 kpc
`i` comparison has 25 pairs; all other reported bins retain 26 `i` or 64–68
`i2` pairs under the 50-pixel floor.

## Interpretation and limits

The delivered coadd variance captures a background/instrumental variance
model, not the missing effective-gain term for adding synthetic-source
Poisson noise. This radial comparison does not estimate the covariance of the
masked ICL estimator, does not correct the ICL profile, and does not include
source shot noise, sky-transfer uncertainty, PSF wings, or astrophysical
covariance. The 5-degree blocks are approximate equal-area connected cells,
not validated HSC survey fields. Keep `science_measurement=false`.

## Reproduction and artifacts

Runner:
[`analyze_icl_variance_radial_balance.py`](analyze_icl_variance_radial_balance.py),
SHA-256 `56b7de553c689fe4986f7ce09cddf62ee848173a956d8d32317aeabaf01eefd8`.

Results:
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_variance_radial_balance_supported94_2026-11-06_01/`.
`summary.json` records input and output hashes and each filter/bin result;
`target_control_radial_comparison.tsv` has 930 supported pair/bin rows,
`per_pair_radial_variance.tsv` has 1,880 target/control/bin rows, and
`image_hashes.tsv` records the 94 verified image pairs. The output hashes are
recorded in the summary. No HSC query or download occurred.
