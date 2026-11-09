# Affine-plane subtraction null for the ICL covariance kernel

## Question

The 33×33 residual-lag kernel built from the 224 supported ICL pairs had large
negative 9–16 pixel shells. This diagnostic estimates how much negative
covariance the local affine-plane subtraction can create by itself.

## Null model and result

For a complete 64×64 patch with independent unit-variance noise, no mask, and
an unweighted ordinary least-squares fit of a constant plus x and y slopes,
the exact residual covariance is (I-H). Summing its expected lag means over
square lag radii gives factors 0.9993 at radius 0, 0.8087 at radius 8, and
0.3629 at radius 16. A 100-realization Gaussian simulation using Astrafek's
actual clipping and lag estimator gives means 0.848 at radius 8 and 0.358 at
radius 16; the finite simulation spread is broad (q16–q84 0.510–1.214 and
−0.159–0.884, respectively), and contains the analytic values.

In HSC data, the median per-frame square-kernel sums at radii 8 and 16 are:

| Filter | Role | Radius 8 | Radius 16 | Shell 9–16 |
|---|---|---:|---:|---:|
| `i` | target | 1.655 | 0.825 | −0.803 |
| `i` | control | 1.559 | 0.755 | −0.819 |
| `i2` | target | 1.661 | 0.826 | −0.804 |
| `i2` | control | 1.611 | 0.740 | −0.824 |

The idealized plane-fit null therefore confirms that affine detrending can
produce a negative long-lag tail. Its expected radius-8-to-16 change (−0.446)
is smaller than the observed median shell (about −0.8), so the null does not
explain the full HSC tail. Masks, heteroscedasticity, residual sky structure,
source masking, and the per-lag median operation remain entangled.

## Interpretation

This demonstrates that the estimated long-lag kernel depends strongly on the
plane-subtraction estimator. It does not identify which remaining effects
cause the additional HSC anticorrelation, and it does not establish the true
coadd covariance. The 33×33 lag medians still need an estimator that avoids
the local projection bias or is calibrated against independent split-coadd
noise. The current finite-lag noise is not applied to ICL fluxes;
`science_measurement=false` remains.

## Reproduction

- Runner: `analyze_icl_plane_subtraction_covariance_null.py` (SHA-256
  `150fedeb2a8e00ec011dea525f0bf3c71c31fe9d66f31dc1bc0463724f9a9aee`).
- Output summary:
  `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_plane_subtraction_covariance_null_2026-11-06_02/summary.json`
  (SHA-256 `bf93d7c0b9189782c9a510cffee6964a45ab87b49bb7d5630ecc6822419b2760`).
- Inputs: the checksum-validated 33×33 lag table and 224-pair exact-aperture
  covariance summary from existing local HSC images. No HSC network request,
  new download, GPU, or Slurm job was used.
