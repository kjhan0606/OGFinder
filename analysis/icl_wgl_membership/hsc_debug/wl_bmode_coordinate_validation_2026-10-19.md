# HSC Y3 B-mode mask coordinate validation and cluster-only stack — 2026-10-19

## Mask polarity validation

The official HSC S19A shape-catalog page says `b_mode_mask` removes the region
with large B-mode, but its prose does not assign a raw integer value. A
published HSC Y3 systematics analysis identifies the excluded GAMA09H rectangle
as RA 132.5–140 degrees and Dec 1.6–5 degrees, with area about 20 deg².
Astrafek tested the approved shape rows against that external sky definition,
without using the shear profile to choose a value.

After the published shape-quality and tomographic-bin cuts, the coordinate
match is exact:

| Source position | raw `b_mode_mask=0` pairs | raw `b_mode_mask=1` pairs | Unique sources (0 / 1) |
|---|---:|---:|---:|
| Inside published rectangle | 77,021 | 0 | 75,684 / 0 |
| Outside published rectangle | 0 | 2,157,489 | 0 / 2,056,176 |

The same separation holds for all approved lens-source pairs at the pair level;
the table's pair counts include repeated appearances of sources in overlapping
lens apertures. This independently decodes **raw value 0 as the excluded
region** for this S19A export. The machine-readable audit records
`b_mode_polarity_decoded=true`, `b_mode_masked_raw_value=0`, and
`raw_value_meaning_inferred_from_shear=false`.

The audit is
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_bmode_spatial_diagnostic_2026-10-19_04/summary.json`
(SHA-256 `1d1311809b502b0f2fd95dfdbb16a4c57b3546dba9855550b1c716fc316d5428`).
Slurm job 428546 completed successfully. The coordinate classifier and exact
two-way decoder are in
[`audit_wl_bmode_spatial_distribution.py`](audit_wl_bmode_spatial_distribution.py).

## Cluster-only weak-lensing stack

Using the now-resolved mask value, Slurm job 428553 ran
[`run_wl_shape_stack.py`](run_wl_shape_stack.py) on the existing 800-lens row
export and official four-bin HSC Y3 SACC `n(z)` file. The runner applies HSC
shape cuts, tomographic `n(z)`, per-object shear calibration, and the catalog's
selection-bias correction. It made no HSC request or download. The first
attempt, Slurm 428549, stopped on a CLI argument-name mismatch before reading
the input; a parser regression test now covers the option and the successful
rerun used the same unused output path.

The completed output is
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_cluster_stack_bmode0_2026-10-19_01/summary.json`
(SHA-256 `bbcf5cb5ef4a03a2a542f21f09d4b90fecdb56e9ac860a0ac7d03435606754c5`).
It retained 2,153,004 lens-source pairs across 609 unique lenses. The table
shows weak-shear `DeltaSigma` and the adapter's independent-pair shape-noise
uncertainty; these are not a total covariance.

| Proper radius (Mpc) | Pairs | Lenses | Tangential `DeltaSigma` ± shape noise (`M_sun pc^-2`) | Cross `DeltaSigma` ± shape noise (`M_sun pc^-2`) |
|---|---:|---:|---:|---:|
| 0.1–0.2 | 11,251 | 541 | 161.5 ± 11.5 | 3.0 ± 11.5 |
| 0.2–0.4 | 37,896 | 564 | 116.3 ± 5.8 | −2.0 ± 5.8 |
| 0.4–0.7 | 92,829 | 577 | 72.4 ± 3.5 | −2.7 ± 3.5 |
| 0.7–1.0 | 131,813 | 586 | 54.6 ± 2.9 | −3.3 ± 2.9 |
| 1.0–1.5 | 304,179 | 594 | 38.5 ± 1.8 | −1.7 ± 1.8 |
| 1.5–2.0 | 414,174 | 599 | 27.9 ± 1.6 | 3.1 ± 1.6 |
| 2.0–3.0 | 1,160,862 | 609 | 15.3 ± 0.9 | 0.05 ± 0.93 |

The tangential profile is positive and declines with radius in this
cluster-only stack. This is a preliminary lensing signal, not a final
measurement: there is no random-lens subtraction, source-overlap correction,
photo-z uncertainty propagation, or spatial-block covariance, and the
weak-shear approximation is used. The 1.5–2.0 Mpc cross component is about two
shape-noise standard errors from zero and warrants a null-test follow-up. The
output keeps `science_measurement=false`.

The exact rectangular mask test is supported by the [official HSC S19A shape
catalog documentation](https://hsc-release.mtk.nao.ac.jp/doc/index.php/s19a-shape-catalog-pdr3/)
and the [HSC Y3 PSF-systematics paper](https://academic.oup.com/mnras/article/525/2/2441/7202336).
