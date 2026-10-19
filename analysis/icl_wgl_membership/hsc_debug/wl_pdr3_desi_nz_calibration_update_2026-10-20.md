# HSC PDR3 weak-lensing redshift calibration update — 2026-10-20

## Finding

The current Astrafek weak-lensing adapter uses the official HSC Y3 four-bin
SACC fiducial `n(z)` with the S19A PDR3 shape export. HSC documents
`hsc_y3_zbin` as the mapping to its four Y3 bins and provides the corresponding
shape calibration fields, so the existing bin mapping is supported. The
remaining limitation is the redshift-distribution calibration and uncertainty,
not a missing bin correspondence.

A newer study, *Full calibration of the tomographic redshift distribution from
the HSC PDR3 Shape Catalog with DESI*, applies DESI DR1/DR2 clustering-redshift
calibration directly to the S19A HSC PDR3 shape catalog. The v3 manuscript was
revised 2026-08-25 and is published as JCAP 05 (2026) 073. It adds high-redshift
ELG and QSO tracers above `z=1.2`, which the earlier HSC Y3 calibration did not
have. The paper reports a small lower-redshift shift for bin 1, agreement with
photometric calibration for bin 2, and higher-redshift shifts for bins 3 and 4
relative to Rau+2022. In the paper's convention, the latter are
`Δz₃=−0.043^{+0.022}_{−0.020}` and `Δz₄=−0.050±0.012`.

This new result is directly relevant to the bin-dependent profile disagreement
in [`wl_tomographic_consistency_assessment_2026-10-20.md`](wl_tomographic_consistency_assessment_2026-10-20.md).
It motivates an `n(z)` calibration sensitivity run; it does not establish that
redshift calibration caused the discrepancy or validate a weak-lensing signal.

## Data and next gate

The paper states that its calibration code and figure data are public in the
authors' [`desi-y3-hsc` repository](https://github.com/JeanCHDJdev/desi-y3-hsc/),
with result files under `results/distributions/`. No files from that directory
have been downloaded into Astrafek or the project data tree. Before importing
any result, inspect the exact file names, format, provenance, and license; keep
only the calibrated `n(z)` products needed for the four-bin lensing calculation
and their posterior samples if available. Do not retrieve additional HSC image,
shape, or per-object photo-z data as part of this calibration update.

Once separately authorized, add a manifest with source URL, revision, file
hashes, and schema checks. Extend the Astrafek SACC loader only after the
published file schema is confirmed. Compare the new distributions to the current
SACC baseline on identical cluster/source pairs, with the same random-lens
subtraction and paired spatial-block jackknife. Propagate redshift-distribution
posterior uncertainty through the pair-wise `Sigma_crit^{-1}` calculation.
Keep `science_measurement=false` until the random controls, source-overlap
covariance, redshift calibration, and remaining shear-systematic gates pass.

## Sources

- Choppin de Janvry et al., [arXiv:2511.18133v3](https://arxiv.org/abs/2511.18133),
  *Full calibration of the tomographic redshift distribution from the HSC PDR3
  Shape Catalog with DESI*.
- Authors' public [calibration repository](https://github.com/JeanCHDJdev/desi-y3-hsc/).
- Official HSC [S19A Shape Catalog (PDR3) documentation](https://hsc-release.mtk.nao.ac.jp/doc/index.php/s19a-shape-catalog-pdr3/).
