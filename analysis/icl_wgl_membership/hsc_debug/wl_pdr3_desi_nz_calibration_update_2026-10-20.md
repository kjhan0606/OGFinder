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
with result files under `results/distributions/`. A metadata-only GitHub API
inspection pinned the current `main` revision to
`abc70af29910e49fce17a63ef3083160d515fb6c` (2026-08-28). The repository API
does not declare an SPDX license; the arXiv paper is marked CC BY 4.0. Retain
the paper citation and check file-level terms when the artifacts arrive.

The authors' README marks `_v_1p1` as the corrected release and describes the
`distributions_*.npz` files as the calibrated distributions. Metadata-only
inspection found two scale-cut variants and their small sidecar metadata files:

| File | Size |
|---|---:|
| `results/distributions/0.3_3_v_1p1/distributions_0.3_3_v_1p1.npz` | 272,416 bytes |
| `results/distributions/0.3_3_v_1p1/nz_res_metadata_0.3_3_v_1p1.json` | 2,494 bytes |
| `results/distributions/1_5_v_1p1/distributions_1_5_v_1p1.npz` | 272,416 bytes |
| `results/distributions/1_5_v_1p1/nz_res_metadata_1_5_v_1p1.json` | 2,492 bytes |
| **Total** | **549,818 bytes** |

These four files are the requested import scope; no files have been downloaded
into Astrafek or the project data tree. The approval request is limited to these
four author-published distribution products for comparing the fiducial 0.3–3
Mpc/h cut with the 1–5 Mpc/h sensitivity cut. Do not retrieve HSC images, shape
catalog rows, per-object photo-z data, or the repository's other derived files
as part of this step. After approval, verify the NPZ schema and whether it
contains posterior samples before implementing an importer; the GitHub README
does not establish the internal array schema.

After authorization, store the four files outside Git and add a manifest with
source URL, pinned revision, file hashes, license notes, and schema checks. Extend
the Astrafek loader only after confirming the published NPZ schema and whether it
provides posterior samples. Compare the two scale cuts and the current SACC
baseline on identical cluster/source pairs, with the same random-lens subtraction
and paired spatial-block jackknife. Propagate redshift-distribution uncertainty
through the pair-wise `Sigma_crit^{-1}` calculation. Keep
`science_measurement=false` until the random controls, source-overlap covariance,
redshift calibration, and remaining shear-systematic gates pass.

## Sources

- Choppin de Janvry et al., [arXiv:2511.18133v3](https://arxiv.org/abs/2511.18133),
  *Full calibration of the tomographic redshift distribution from the HSC PDR3
  Shape Catalog with DESI*.
- Authors' public [calibration repository](https://github.com/JeanCHDJdev/desi-y3-hsc/).
- Official HSC [S19A Shape Catalog (PDR3) documentation](https://hsc-release.mtk.nao.ac.jp/doc/index.php/s19a-shape-catalog-pdr3/).
