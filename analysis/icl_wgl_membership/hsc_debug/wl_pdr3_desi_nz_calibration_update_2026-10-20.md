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

The user approved exactly these four files on 2026-10-20. They were downloaded
from the pinned commit into
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/desi_y3_hsc_nz_abc70af2_2026-10-20/`.
The external-data manifest records the pinned URLs, exact byte counts, and
SHA-256 values. No HSC images, shape rows, per-object photo-z data, or other
repository files were fetched. Schema inspection found four bins on a shared
1,000-point `z=0..3` grid, with mean/median and pointwise lower/upper curves,
including magnification-corrected variants. The NPZ files do not contain
posterior samples.

The approved curves were compared to the SACC baseline with the same verified
800-lens/source-row export and shared 5-degree equal-area spatial jackknife.
Both DESI scale cuts were tested for mean and magnification-corrected mean; the
published lower/upper curves were also recomputed as pointwise-envelope
sensitivity profiles. The envelope curves are not joint posterior draws and do
not supply a covariance. The full assessment and validation are documented in
[`wl_desi_nz_distribution_sensitivity_assessment_2026-10-20.md`](wl_desi_nz_distribution_sensitivity_assessment_2026-10-20.md).
The central curves change median `DeltaSigma` by +12.5–14.8% in source bin 1
and −1–3.5% in bins 2–4. The 0.1–0.7 Mpc source-bin inconsistency is reduced
but remains. The random-center shape export is unavailable, so this stage is
cluster-only and cannot pass the science gate. Keep
`science_measurement=false` until the random controls, source-overlap covariance,
redshift calibration, and remaining shear-systematic gates pass.

## Sources

- Choppin de Janvry et al., [arXiv:2511.18133v3](https://arxiv.org/abs/2511.18133),
  *Full calibration of the tomographic redshift distribution from the HSC PDR3
  Shape Catalog with DESI*.
- Authors' public [calibration repository](https://github.com/JeanCHDJdev/desi-y3-hsc/).
- Official HSC [S19A Shape Catalog (PDR3) documentation](https://hsc-release.mtk.nao.ac.jp/doc/index.php/s19a-shape-catalog-pdr3/).
