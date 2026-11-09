# Random-only weak-lensing placebo stack (2026-11-06)

## Question and method

The 435-lens cluster-minus-random analysis uses five matched random centers per
lens. This diagnostic asks whether those random centers alone produce a
residual tangential or cross profile. It revalidates the approved 3,100-center
random shape-row export and 800-lens source-row export, retains exactly the
same 435 common-support lenses and their 2,175 random centers, and recomputes
the HSC Y3 tomographic estimator with the saved redshift calibration, source
cuts, B-mode mask, radial bins, and cosmology. The full random-only profile
reproduces the saved matched-random component of the cluster-minus-random
stack.

Uncertainty is estimated from 50 delete-one 5-degree equal-area spatial
blocks. Each replicate removes a random lens-source pair if its center or
source falls in the deleted cell, then reruns the complete calibrated stack.
Separate covariance tests compare the seven-bin tangential and cross profiles
with zero.

## Results

| Radius (proper Mpc) | Random-only tangential ΔΣ (10¹² M☉ Mpc⁻²) | Random-only cross ΔΣ (10¹² M☉ Mpc⁻²) |
|---|---:|---:|
| 0.1–0.2 | −4.03 ± 6.31 | 12.37 ± 7.47 |
| 0.2–0.4 | 1.81 ± 3.03 | −1.89 ± 3.00 |
| 0.4–0.7 | 1.53 ± 2.11 | 1.15 ± 2.60 |
| 0.7–1.0 | 2.23 ± 1.52 | 0.18 ± 1.61 |
| 1.0–1.5 | 1.10 ± 1.05 | −0.01 ± 0.98 |
| 1.5–2.0 | 1.19 ± 1.33 | −0.14 ± 1.25 |
| 2.0–3.0 | −0.38 ± 1.00 | −0.74 ± 0.71 |

The covariance-conditional global tangential null is χ²=6.90 for rank 7
(approximate p=0.439). The cross-component null is χ²=4.06 for rank 7
(approximate p=0.773). Neither random-only vector has a resolved departure
from zero under this jackknife. The first-bin cross estimate is nonzero at
about 1.7 standard errors, but the full cross profile is consistent with the
null and the bins are correlated.

## Interpretation and limits

The random-only null is a useful estimator/selection placebo: it finds no
large residual profile in the matched random sample. It does **not** validate
the cluster tangential signal, which remains sensitive to the footprint's
three occupied connected tract regions and to calibration systematics. The
5-degree cells are approximate equal-area sky blocks, not validated
contiguous HSC survey regions. The global chi-square p-values are
covariance-conditional approximations; photo-z calibration, shear response,
reduced shear, magnification, source selection, and full regional covariance
are not propagated. `science_measurement=false` remains.

## Reproducibility

- Astrafek runner: `analyze_wl_random_only_placebo.py`, SHA-256
  `44fb39796180f452dc6283f443cb5b4c0de21cc6fd6fcdcb7bcf3034ac35606b`.
- Result directory:
  `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_random_only_placebo_2026-11-06_01/`.
- Summary SHA-256:
  `9c2b6ad3f1525f18e26c8bf3ec947e40de49abc40271dd75313dad32b2b50fa2`.
- Random-only covariance SHA-256:
  `c312159b7efc861b299f7b513348fdbac1bebbcb2033e99d19ea3f6793d4ed99`.
- Astrafek plotter: `plot_wl_random_only_placebo.py`, SHA-256
  `0baf2ccb1ff376cae335c29a637883f958eb2de3477db356568cfa13d3d77fd7`.
- Figure: `random_subtracted_vs_random_only.png` (SHA-256
  `85d796d51338b74dc26412ba1e9c58d2cdf372b3f28a5430aeac13141b14c374`);
  vector PDF SHA-256 `5c3af5b555956d8153f03fc25d72b1745a1dc534c6c06f0f15d8b6933c25e772`.
- Inputs are the saved paired random stack, verified random/lens shape-row
  exports, and the existing HSC Y3 SACC file. No HSC query or image download
  occurred.
