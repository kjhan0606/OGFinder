# Full-sample ICL annular covariance: 32-pixel sensitivity (2026-11-06)

## Scope and method

This sensitivity extends the local delivered-coadd lag-covariance audit from
the 94 pairs with visit-level warp products to all 783 same-filter target and
matched-control pairs in the already downloaded HSC sample. It reads saved
coadd image and mask FITS files only; it makes no HSC request and downloads no
data. The tested annulus is 200–300 proper kpc.

The primary full-sample run used non-overlapping 32×32-pixel patches, up to 64
patches per center, at least 8 accepted patches per center, at least 80% valid
pixels per patch, and a square lag kernel through ±8 pixels. The estimator
removes a local affine plane, clips, and computes variance-normalized lag
covariance. Target-minus-control uncertainty uses 10,000 connected 5-degree
block-bootstrap draws per filter with a shared max-T family across zero lag,
the cumulative ±8-pixel lag sum, and annulus coverage.

## Support

| Filter | Same-filter pairs | Complete target/control annuli | Support | Connected blocks |
|---|---:|---:|---:|---:|
| i | 302 | 215 | 71.2% | 13 |
| i2 | 481 | 402 | 83.6% | 35 |

This is much broader than the 64×64 sensitivity (17/302 i and 48/481 i2), but
the supported systems remain selected. Complete-minus-incomplete redshift SMD
is −0.65 in i and −0.39 in i2; richness SMD is +0.08 and −0.24. The 32-pixel
results therefore do not represent the full parent sample without further
selection analysis.

## Familywise-adjusted target–control contrasts

| Filter | Metric | Median difference | Simultaneous 95% interval | max-T p |
|---|---|---:|---:|---:|
| i | zero lag | −0.00021 | [−0.00297, +0.00255] | 0.985 |
| i | cumulative lag sum ±8 px | +0.0976 | [+0.0203, +0.1748] | 0.0265 |
| i | annulus coverage | −0.0591 | [−0.0983, −0.0198] | 0.0133 |
| i2 | zero lag | −0.00194 | [−0.00384, −0.00004] | 0.0480 |
| i2 | cumulative lag sum ±8 px | +0.0395 | [−0.0107, +0.0898] | 0.164 |
| i2 | annulus coverage | −0.0446 | [−0.0690, −0.0202] | 0.0002 |

The i-band cumulative lag sum and the i2 zero-lag ratio have small adjusted
intervals that exclude zero in this 32-pixel specification. Coverage differs
in both bands. The 64-pixel specification had no familywise interval excluding
zero for its three metrics. This scale sensitivity, together with redshift
selection and unequal valid coverage, means these are diagnostics of the
delivered coadds and masks rather than a calibrated noise correction. The
cumulative lag sum can also contain diffuse cluster/ICL residual structure;
zero-lag variance is not a complete aperture covariance model.

The input run's numeric provenance records `patch_size_pixels=32`. Its
human-readable `method.patches` description was accidentally hard-coded to
64×64. The runner and support auditor have been corrected to generate this
description from the numeric configuration; the original input tables and
summary remain preserved. The corrected support audit below records the actual
32×32 configuration explicitly.

No correction is applied to the ICL or weak-lensing measurements, and
`science_measurement=false` remains in force. Remaining requirements include a
representative support model, exact coadd-weight/resampling covariance,
source-Poisson response, and sky-transfer validation.

## Reproducibility

- Analysis runner: `analyze_icl_full783_annular_covariance.py`, executed SHA-256
  `abd1912aaf052146f31b99e60b8387cd41b7ec8a30ba33b4ef055aab08279b1d`.
- Analysis summary:
  `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_annular_coadd_covariance_full783_32px_2026-11-06_01/summary.json`, SHA-256
  `6e9c786f10116de6973bdf70dda7c0d335d4cc71df717d5ab207d8d42ab25bd2`.
- Paired covariance table SHA-256
  `355c7943cc78a597e95367bf671e902fad40b9b4f859d3734780f21a3b6d1675`;
  per-center table SHA-256
  `59227ae694ca7fd8021f25ec9dbd2e97361304bd4e0c01fc3b6b1f3798ea79fe`.
- Corrected support auditor: `audit_icl_full783_annular_covariance_support.py`,
  SHA-256 `10372f3b00c6903506ac6ab0c2d4fe3a8d62dd5143d5d8ca42db2983e38f2a65`.
- Corrected support summary:
  `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_annular_covariance_full783_32px_support_audit_2026-11-06_02/summary.json`, SHA-256
  `8ad6bcba35b9ce9e7e87a014b97608b4f095b85cebab1e4639b76e0e0357343a`.
- HSC network used: false; image pixels read by the covariance analysis: true;
  science measurement: false.
