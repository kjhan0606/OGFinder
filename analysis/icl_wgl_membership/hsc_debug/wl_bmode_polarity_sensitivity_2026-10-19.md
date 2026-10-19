# HSC Y3 cluster-lensing stack: B-mode bit-polarity sensitivity

**Status:** completed two-hypothesis sensitivity calculation; `science_measurement=false`.

## Input and method

Astrafek's `run_wl_bmode_polarity_sensitivity.py` checksum-verified all 32
previously approved shape-row batches for 800 CAMIRA lenses (3,069,544
lens–source pairs), loaded the approved HSC Y3 SACC file, and computed a
cluster-centered tomographic stack twice. Both runs use the same HSC Y3
selection cuts, per-object multiplicative/additive terms, radial selection-bias
correction, fiducial bin `n(z)`, and flat cosmology (`H0=70 km/s/Mpc`,
`Omega_m=0.3`). The calculation used the seven proper-radius bins from 0.1 to
3 Mpc used in the prior raw-profile diagnostic.

The HSC documentation says `b_mode_mask` removes regions with large B-mode,
but it does not identify whether raw 0 or 1 denotes the masked region. The two
cases below are therefore explicitly competing hypotheses. No polarity was
chosen from the observed shear profile or bit frequency.

## Results

`DeltaSigma` is the weak-shear approximation in `Msun/pc^2`; the uncertainty
shown is the adapter's independent-pair shape-noise term only. It is not a
survey or total covariance.

| Proper R (Mpc) | Pairs if raw 0 is masked | `DeltaSigma` ± shape noise | Pairs if raw 1 is masked | `DeltaSigma` ± shape noise |
|---|---:|---:|---:|---:|
| 0.1–0.2 | 11,251 | 161.5 ± 11.5 | 430 | 58.1 ± 61.1 |
| 0.2–0.4 | 37,896 | 116.3 ± 5.8 | 1,437 | 122.1 ± 30.2 |
| 0.4–0.7 | 92,829 | 72.4 ± 3.5 | 3,468 | 60.1 ± 18.2 |
| 0.7–1.0 | 131,813 | 54.6 ± 2.9 | 5,045 | 57.6 ± 14.6 |
| 1.0–1.5 | 304,179 | 38.5 ± 1.8 | 11,020 | 31.9 ± 9.8 |
| 1.5–2.0 | 414,174 | 27.9 ± 1.6 | 14,048 | 28.6 ± 8.5 |
| 2.0–3.0 | 1,160,862 | 15.3 ± 0.9 | 41,416 | 12.8 ± 4.9 |

Across 0.1–3 Mpc, the raw-zero-masked hypothesis retains 2,153,004 pairs and
609 lenses; the raw-one-masked hypothesis retains 76,864 pairs and 26 lenses.
The raw-zero-masked profile has positive tangential signal decreasing with
radius. Its cross component is generally near zero, although individual bins
are comparable to or somewhat larger than the independent-pair shape-noise
term. The raw-one-masked hypothesis has larger cross-component excursions and
shape-noise errors comparable to the tangential values because so few lenses
remain. Those differences do not determine the bit meaning: choosing whichever
case looks more lens-like would be circular.

## Quality gates and limits

The run completed on Slurm job 427885 in 71 seconds, with 1.6 GiB peak batch
memory. The output JSON is
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_bmode_polarity_sensitivity_2026-10-19_01/stack.json`
and has SHA-256
`67adb9a2a6abb67a6dac6c43d09c2c8699b8378ba0c782464cd4ad25a1e4163c`.
It records the source-row batch hashes, row-plan/count provenance, and SACC
checksum (`a28f9e2e088e92d96d2d87083a6eeeac4c0c84c6b48e033bd3d1dc1bf0d8958f`).

The SACC distributions are fiducial population-level `n(z)` and their
posterior/systematic uncertainty is not propagated. This cluster-only stack
has no matched random-lens subtraction, source-overlap correction, or spatial
block covariance; it also uses the weak-shear approximation and does not
include a reduced-shear fit. Consequently these profiles are exploratory, not
a mass measurement or detection, and `science_measurement=false` remains in
the output.

## Spatial polarity diagnostic prepared

The HSC release documentation describes `b_mode_mask` as the cut for the
large-B-mode region, but does not define whether raw zero or one is the
excluded value. Astrafek now has
[`audit_wl_bmode_spatial_distribution.py`](audit_wl_bmode_spatial_distribution.py),
which verifies the 32 approved row-batch hashes, reapplies the published
shape/tomographic cuts, deduplicates repeated sources by `object_id`, and
counts both raw values on a one-degree sky grid. It deliberately does not use
the measured shear profile to select a polarity. The Slurm runner is
[`run_wl_bmode_spatial_distribution.sbatch`](run_wl_bmode_spatial_distribution.sbatch).

The first submitted job (Slurm 428534) failed before producing an output
because the approved CSV export prefixes its first header field with `# `.
The parser now strips that marker only from the first field and retains the
exact expected-column check. The successful rerun was Slurm 428536; its output
is `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_bmode_spatial_diagnostic_2026-10-19_02/summary.json`
(SHA-256 `8260ad1cba3b926af322a5407751762b1741cbe7f75a5cd9c7c08fa16fa90448`).
All 3,069,544 exported lens-source rows were read. After the published cuts,
2,234,510 pair rows remained, representing 2,131,860 unique sources. Raw value
0 accounts for 77,021 pairs (3.45%) and 75,684 unique sources; raw value 1
accounts for 2,157,489 pairs and 2,056,176 unique sources.

At the one-degree coordinate-grid resolution, raw value 0 appears in 26 of
495 occupied cells, spanning RA 132–140 degrees and Dec 1–5 degrees. Nineteen
cells contain only raw 0 and seven contain both values; all other occupied
cells contain only raw 1. This is strong spatial localization in the approved
800-lens export and is consistent with one value tracing a survey footprint or
mask. The grid is not an area-corrected density map, and at that stage this correspondence
did not establish which bit the catalog excludes. The output explicitly
retains `b_mode_polarity_decoded=false`, `raw_value_meaning_inferred_from_shear=false`,
and `science_measurement=false`.

The HSC PDR3 documentation says public random points are available in its
database, but our current worktree has no approved random-lens export. The next
data step is to inspect that product's schema and propose a count-only matched
random-position query before requesting approval for any row export.

References: [HSC S19A shape-catalog documentation](https://hsc-release.mtk.nao.ac.jp/doc/index.php/s19a-shape-catalog-pdr3/);
[HSC PDR3 available data](https://hsc-release.mtk.nao.ac.jp/doc/index.php/available-data__pdr3/).

## Historical status update

The text above describes the state before the coordinate-to-footprint match. The subsequent
[`wl_bmode_coordinate_validation_2026-10-19.md`](wl_bmode_coordinate_validation_2026-10-19.md)
exact-matched raw value 0 to the published HSC Y3 excluded GAMA09H rectangle and raw value
1 to outside positions, independently of the measured shear. A follow-up cluster-only
calibrated stack using raw 0 as the masked value is documented there. That stack remains
exploratory because random-lens subtraction, source-overlap and spatial covariance, and
photo-z uncertainty are still missing.
