# HSC Y3 published mean-shift sensitivity — 2026-10-20

## Method and result

With no new data download, Astrafek translated the local SACC `n(z)` curves for
source bins 3 and 4 toward higher redshift by `+0.043` and `+0.050`, respectively,
then renormalized them on the same redshift grid. The values are the central
mean shifts reported by the DESI DR1/DR2 calibration paper for HSC PDR3; its
paper convention gives `Delta z_3=-0.043` and `Delta z_4=-0.050` while describing
the distributions as shifted to higher redshift ([v3 paper](https://arxiv.org/abs/2511.18133v3)).
The translated SACC curves retain the reported mean shifts to `2.4×10⁻⁶` or
better after grid truncation.

The same 800-lens row export, HSC cuts, four source bins, 52 spatial cells,
radial bins, and joint delete-one-cell jackknife were used as in the baseline.
Of 3,069,544 approved lens–source rows, 2,254,646 are in source bins 1–4 and
2,153,004 pairs pass the shape, B-mode, distance, and radial selections. The
jackknife covariance is 56×56, rank 51, finite, and symmetric to relative
precision `2.6×10⁻¹⁷`.

| R (Mpc) | Baseline χ²/3 | Baseline Holm p (7) | Shifted χ²/3 | Shifted Holm p (7) |
|---|---:|---:|---:|---:|
| 0.1–0.2 | 25.80 | 6.29×10⁻⁵ | 25.27 | 8.12×10⁻⁵ |
| 0.2–0.4 | 45.15 | 6.02×10⁻⁹ | 43.60 | 1.28×10⁻⁸ |
| 0.4–0.7 | 21.22 | 4.74×10⁻⁴ | 19.91 | 8.86×10⁻⁴ |
| 0.7–1.0 | 2.18 | 1.000 | 1.65 | 1.000 |
| 1.0–1.5 | 0.51 | 1.000 | 0.43 | 1.000 |
| 1.5–2.0 | 6.26 | 0.398 | 6.12 | 0.423 |
| 2.0–3.0 | 3.03 | 1.000 | 2.77 | 1.000 |

The inner `0.1–0.7 Mpc` source-bin inconsistency remains significant after the
published high-bin mean shifts. At `0.2–0.4 Mpc`, bins 3 and 4 move by only
`−3.86×10¹²` and `−2.25×10¹² M☉ Mpc⁻²`; bins 1 and 2 remain unchanged. These
central shifts alone therefore do not explain the mismatch.

## Limits

This is a constant-translation approximation of two fiducial SACC curves. It is
not the authors' calibrated spline distributions or posterior samples. It omits
the small bin-1 shift to lower redshift, the reported shift uncertainties, and
other changes to distribution shape. The DESI NPZ and metadata files were not
downloaded; their import still requires the user's explicit approval. The stack
also has no random-lens subtraction or connected-footprint covariance. Keep
`science_measurement=false`.

## Reproducibility

- Astrafek runner: [`audit_wl_published_nz_shift_sensitivity.py`](audit_wl_published_nz_shift_sensitivity.py)
- Focused tests: [`test_audit_wl_published_nz_shift_sensitivity.py`](test_audit_wl_published_nz_shift_sensitivity.py)
- Slurm job 428861: completed in 12m45s; maximum RSS 1,968,092 KiB.
- Summary: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_published_nz_shift_sensitivity_2026-10-20_01/summary.json`
- Summary SHA-256: `5ecee0f5c41697cbda827f7a4ccb706fb3927b00c0228c7688e15df52baaf3f4`
- Baseline summary SHA-256: `39b967e91d27ba3ebe1f05f24779245afff2fa197077469c8227d4d5efcfa643`

The output records the row, SACC, baseline-summary, and runner hashes. All
numerical products remain outside Git.
