# HSC Y3 tomography split by lens redshift — 2026-10-20

## Result

Astrafek split the unique lenses represented in the approved shape-pair table
into two equal-count cohorts after sorting by `(lens_z, lens_id)`. There are
632 represented lenses, 316 per cohort; the boundary is between `z=0.42950`
and `z=0.42955`. Of the original 3,069,544 lens–source rows, 2,254,646 have
HSC Y3 source-bin assignments 1–4, and 2,153,004 pairs pass the published
shape, B-mode, source-distance, and 0.1–3 Mpc cuts. The full selected sample
contains 609 unique contributing lenses. Each split has 52 occupied 5-degree
equal-area cells and uses the same delete-one-cell jackknife for all source bins.

The tangential `DeltaSigma` consistency test compares four source-bin estimates
of one common profile within each lens-redshift cohort. The table gives the
four-bin chi-square with 3 degrees of freedom and its Holm-adjusted p-value
across all 14 (7 radii × 2 cohorts) tests. Values are diagnostic and are not
lensing detections.

| Lens cohort | R (Mpc) | χ²/dof | Holm p (14 tests) | Result |
|---|---:|---:|---:|---|
| lower z | 0.1–0.2 | 14.98/3 | 0.0220 | inconsistent |
| lower z | 0.2–0.4 | 28.48/3 | 3.75×10⁻⁵ | inconsistent |
| lower z | 0.4–0.7 | 11.74/3 | 0.0916 | not significant after correction |
| lower z | 0.7–1.0 | 1.25/3 | 1.000 | consistent in this test |
| lower z | 1.0–1.5 | 1.34/3 | 1.000 | consistent in this test |
| lower z | 1.5–2.0 | 7.46/3 | 0.585 | consistent in this test |
| lower z | 2.0–3.0 | 1.01/3 | 1.000 | consistent in this test |
| upper z | 0.1–0.2 | 6.60/3 | 0.772 | consistent in this test |
| upper z | 0.2–0.4 | 28.82/3 | 3.42×10⁻⁵ | inconsistent |
| upper z | 0.4–0.7 | 6.17/3 | 0.827 | consistent in this test |
| upper z | 0.7–1.0 | 2.14/3 | 1.000 | consistent in this test |
| upper z | 1.0–1.5 | 3.33/3 | 1.000 | consistent in this test |
| upper z | 1.5–2.0 | 1.50/3 | 1.000 | consistent in this test |
| upper z | 2.0–3.0 | 4.41/3 | 1.000 | consistent in this test |

The inner `0.2–0.4 Mpc` disagreement occurs in both lens-redshift cohorts,
so it is not confined to the low-z half. The additional low-z-only result at
`0.1–0.2 Mpc` is consistent with greater sensitivity to foreground dilution
and the poorly separated first source bin, but this split alone cannot identify
the cause. The plotted profiles show source bin 1 near zero or negative in the
innermost bins while bins 2–4 are positive; the upper-z bin-1 errors are broad.

## Interpretation and next test

This is a cluster-only stack using the fiducial HSC Y3 `n(z)`. It has no
random-lens subtraction, no posterior/systematic `n(z)` propagation, no
connected-footprint covariance, no source-overlap correction, and no
reduced-shear correction. It therefore diagnoses a calibration/selection
consistency issue rather than measuring a cluster mass signal. Keep
`science_measurement=false`.

The first run computed the 112×112 joint covariance in memory but did not
persist the matrix in its summary. The per-cohort consistency statistics above
are from that covariance, but a direct test of whether relative source-bin
offsets change with lens redshift cannot be recovered from the saved diagonal
errors. The corrected Astrafek runner now writes the full covariance. Slurm job
428856 is rerunning the same approved data to a new directory; its saved
covariance will feed a difference-in-differences interaction test that permits
each cohort its own common lensing amplitude. No additional catalog data are
requested for that test.

## Reproducibility

- Analysis: [`audit_wl_tomographic_lens_redshift.py`](audit_wl_tomographic_lens_redshift.py)
- Plotter: [`plot_wl_tomographic_lens_redshift.py`](plot_wl_tomographic_lens_redshift.py)
- Interaction test: [`audit_wl_tomographic_redshift_interaction.py`](audit_wl_tomographic_redshift_interaction.py)
- Initial Slurm job 428849: completed in 12m27s; maximum RSS 2.045 GiB.
- Initial machine-readable output: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_tomographic_lens_redshift_2026-10-20_01/summary.json`
- Initial summary SHA-256: `e80d85b031b6184511cfcb475d3c19da41c685a8867859034dba333d4aace43c`
- Profile figure: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_tomographic_lens_redshift_2026-10-20_01/tomographic_profiles_by_lens_redshift_v2.png`
- Figure SHA-256: `cacff2c9ed94c100015452a464b32ce42aee2ed68e0878ff0a8b6ccacf2ededb`

The CSV row hashes and SACC checksum are retained in the machine-readable
summary. The figure and all numerical outputs remain outside Git.
