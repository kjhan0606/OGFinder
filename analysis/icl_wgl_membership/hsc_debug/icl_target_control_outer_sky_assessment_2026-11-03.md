# Full-sample target/control outer-residual assessment — 2026-11-03

## Scope and processing

Astrafek measured 600–900 proper-kpc outer-annulus tile residuals in the existing
HSC PDR3 global-sky `coadd/bg` cutouts for all 783 same-filter target/matched-
control pairs (302 `i`, 481 `i2`). The run processed 1,566 frames with 64-pixel
tiles, a saved target K=2.5 mask, and the same SEP/Kron-2.5 mask family for
matched controls. Tile intensities and variances were converted to nJy/arcsec².
No input pixels were modified, subtracted, or downloaded. The existing acquisition
SHA values and file sizes were cross-referenced; image content was not rehashed in
this pass.

Slurm job 448403 completed on `syntax` in 26:08, with 610,944 KiB maximum RSS and
exit 0. The independent table audit verified 1,566 frame rows, 783 unique
same-filter pairs, the pairwise tile/intercept arithmetic, and 38 connected
5-degree blocks. Its result is
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_target_control_outer_sky_2026-11-03_01/independent_output_audit.json`
(SHA-256 `990209f909dd4ea57479c9f466cb8b55ca391cc5f9e56837e6b82e980593b099`).

## Outer residual distribution

The fitted 600–900 kpc target-minus-control plane intercepts have these
per-pair distributions. The quantiles below describe the observed pair
distribution; they are not confidence intervals.

| Filter | Pairs | Median intercept difference | 16th–84th percentile | Median target/control residual-MAD to nominal-sigma |
|---|---:|---:|---:|---:|
| i | 302 | −0.055 nJy/arcsec² | [−3.942, +3.816] | 3.16 / 3.10 |
| i2 | 481 | −0.370 nJy/arcsec² | [−4.863, +3.596] | 3.10 / 3.04 |

The median outer residual is near zero in both filters, while individual target-
control differences span several nJy/arcsec². The tile-plane residual scatter
is also about three times the nominal independent-pixel tile-mean sigma. Neither
statistic is a calibrated sky-error multiplier: the target annulus can contain
real cluster/ICL light, PSF wings, cirrus, undetected sources, and mask-dependent
structure. No residual was subtracted from the science profile.

## ICL–satellite association sensitivity

For the existing 74-pair association sample, the pixel-SFD 200–300 kpc ICL
contrast and the fixed center-SFD satellite 50–75 kpc contrast were residualized
against redshift and richness. We then added the paired target-minus-control
outer-plane intercept as a third linear covariate. Each filter used 5,000
paired resamples of connected 5-degree blocks formed from that filter's
association sample alone.

| Filter | Pairs / blocks | Baseline partial r [95% block interval] | With outer-residual covariate [95% interval] | Adjusted − baseline [95% paired interval] |
|---|---:|---:|---:|---:|
| i2 (primary) | 53 / 28 | 0.239 [0.022, 0.423] | 0.252 [0.036, 0.456] | +0.013 [−0.034, +0.111] |
| i (descriptive) | 21 / 13 | 0.112 [−0.229, +0.514] | 0.107 [−0.273, +0.535] | −0.005 [−0.137, +0.106] |

The point estimates change little, and both paired intervals include zero. This
does not show that sky residuals are negligible: the outer indicator may contain
physical cluster light and shares the same selected imaging sample. It is a
covariate sensitivity check, not a correction or a confirmatory ICL–satellite
detection.

### Block-assignment correction

The first 783-pair run assigned the 74 association rows their labels from the
full 783-pair connected-component graph. Unrelated clusters could therefore
bridge separate association regions. That historical run reported only 20 i2
and 8 i connected blocks. The image measurements and pair table remain valid,
but its association-bootstrap intervals are superseded and are not used above.
The corrected replay recomputed blocks inside the fixed association sample for
each filter, restoring 28 i2 and 13 i blocks. This produced wider, more
appropriate intervals. Both the historical table arithmetic and the corrected
replay were independently audited; the latter audit is
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_outer_sky_association_reblocked_2026-11-03_01/independent_audit.json`
(SHA-256 `fbabe476139839ea6287109a649d8a4a249c7b7fda40c33f3bcf0408d87c84cd`).

## Reproducibility and status

- Full-frame runner: [`audit_icl_target_control_outer_sky.py`](audit_icl_target_control_outer_sky.py)
- Association-sample corrected replay: [`replay_icl_outer_sky_association.py`](replay_icl_outer_sky_association.py)
- Full-sample output audit: [`verify_icl_target_control_outer_sky.py`](verify_icl_target_control_outer_sky.py)
- Corrected association replay audit: [`verify_icl_outer_sky_association_replay.py`](verify_icl_outer_sky_association_replay.py)
- Full pair table SHA-256: `94691b193244209960dfebe55047b3270455b745818bdc3e05f2fa55a20b6e0f`
- Full run summary SHA-256: `581b173f37dc3a7f3e60bd71a676e69a89e38f6fb3d76b194f33b5eb19968808`
- Reblocked association summary SHA-256: `89c0f39d72f76a8f4864a68f8ca4a09fa05fe0d047f937c888055c81b5c8c983`

The research status remains `science_measurement=false`. Upstream HSC sky
transfer, extended PSF wings, member/satellite completeness, extinction and
passband systematics, and full image-pixel covariance remain uncalibrated. This
diagnostic does not promote the observed profile or settle the physical
ICL–satellite association.
