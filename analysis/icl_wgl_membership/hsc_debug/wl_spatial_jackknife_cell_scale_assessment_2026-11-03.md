# Weak-lensing spatial-jackknife cell-scale sensitivity — 2026-11-03

## Scope and method

Astrafek reran the existing cluster-only HSC Y3 tomographic spatial jackknife at
equal-area sky-cell widths of 2.5 and 7.5 degrees and compared both with the
saved 5-degree baseline. The same 800-lens shape-row export, published source
cuts, raw B-mode exclusion value 0, fiducial HSC Y3 SACC `n(z)`, cosmology, and
proper-radius bins were used in all three runs. A replicate removes every pair
whose lens or source falls in the deleted cell and recomputes the calibrated
profile. Only the spatial grid width changes.

The two new CPU-only jobs ran on `syntax`: Slurm 448324 (2.5 degrees, 28:54,
2.47 GiB max RSS) and 448323 (7.5 degrees, 7:09, 2.33 GiB max RSS). Both exited
0. They read existing approved data; no HSC query, download, or source-row
export occurred.

## Results

| Cell width | Occupied cells | Used lens–source pairs | Cross-component chi-square / dof | Diagnostic p-value |
|---:|---:|---:|---:|---:|
| 2.5° | 156 | 2,153,004 | 7.10 / 7 | 0.418 |
| 5.0° | 52 | 2,153,004 | 9.30 / 7 | 0.232 |
| 7.5° | 36 | 2,153,004 | 6.21 / 7 | 0.515 |

The seven tangential `DeltaSigma` central values are identical across the runs,
as expected because the full-sample stack is unchanged. Across radial bins, the
jackknife errors at 2.5 degrees range from 0.829 to 1.167 times the 5-degree
errors; at 7.5 degrees they range from 0.878 to 1.177 times the 5-degree
errors. Thus the reported diagonal errors shift by up to about 18% over these
three grid choices. All three cross-component chi-square diagnostics are
consistent with zero under their stated Gaussian approximation, but their
p-values depend on the grid and are not detections or validation of the
tangential profile.

## Interpretation and status

The decreasing cluster-centered profile and its full-sample normalization are
stable under this grid-width perturbation. The jackknife uncertainty is
moderately sensitive to the arbitrary partition. These cells are not HSC
tracts or guaranteed independent fields; equal-area bins can combine
disconnected survey regions and do not replace a paired cluster-minus-random
covariance. The exercise also does not add random-lens subtraction, source
selection or shear-calibration systematics, uncertainty in `n(z)`, or
reduced-shear corrections. Raw B-mode value 0 remains the working exclusion
assumption, not a newly decoded catalog convention.

This is a cluster-only covariance sensitivity diagnostic. Keep
`science_measurement=false`; do not interpret the displayed `DeltaSigma` as a
final mass profile. The prepared random-center shape-count feasibility query
remains unsubmitted pending explicit approval for that separate aggregate HSC
export.

## Reproducibility

- Runner: [`run_wl_cluster_spatial_jackknife.py`](run_wl_cluster_spatial_jackknife.py)
- CPU-only Slurm wrapper: [`run_wl_cluster_jackknife_scale.sbatch`](run_wl_cluster_jackknife_scale.sbatch)
- 2.5-degree output:
  `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_cluster_jackknife_scale_2p5deg_2026-11-03_01/summary.json`
  SHA-256 `714a399f57aaf7935e38214d6d059f9241f1137f34b6641559fc74d33e35060b`
- 5-degree baseline:
  `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_cluster_spatial_jackknife_2026-10-19_01/summary.json`
  SHA-256 `16579be08d02ee29f5ad912f9274f64b728bafc408ddca0746de3ad69494a6c3`
- 7.5-degree output:
  `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_cluster_jackknife_scale_7p5deg_2026-11-03_01/summary.json`
  SHA-256 `8bb4aae4f046ebbcecd5a104d5d1bb6d131c1aed5af947bec6bb6ef3d0b5722c`

All three outputs record `science_measurement=false` and retain the full
covariance and block-deletion diagnostics.
