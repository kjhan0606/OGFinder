# HSC S19A shape-row quality and published source cuts (2026-10-19)

## Scope

This streaming audit checks the already approved and downloaded HSC S19A lens-source
export for pair counts, finite values, sky-coordinate ranges, shape-weight signs,
tomographic assignments, B-mode raw-value counts, and source-selection support. It
does not submit CAS queries, download data, infer B-mode polarity, calculate shear,
or promote the rows to a science measurement.

The export contains 800 lens IDs and 3,069,544 lens-source pairs in 32 SHA-256
verified CSV batches (909,833,598 bytes). Repeated sources across different lens
cones remain separate lens-source pairs. The per-lens median is 4,438 candidate rows;
168 lenses have no rows in the 3 Mpc query cone.

## Selection-rule correction

The HSC Y3 shape-catalog paper's weak-lensing sample requires ellipticity amplitude
`|e|<2`, measurement uncertainty `0<=sigma_e<=0.4`, resolution `R2>=0.3`, and
1-arcsec aperture magnitude `mag_A<=25.5`. It defines `e1/e2` as re-Gaussianization
ellipticity components; those components are noisy measurements, so the adapter must
apply the published sample cut rather than impose the noiseless geometric bound
`|e|<1` on every row. The official paper also describes the first two quality
thresholds and the resolution/magnitude boundaries ([Li et al. 2022, HSC Y3 shape
catalog](https://arxiv.org/abs/2107.00136), Table 2).

The adapter previously rejected any input containing `|e|>=1` before it could apply
source selection. In this export, 88,460 rows have `1<|e|<2`, while none reaches
`|e|>=2`. The prior validator would therefore have stopped the whole stack on
otherwise allowed catalog rows. The adapter now requires `sigma_e`, accepts the
documented shape range, and applies the published ellipticity, uncertainty,
resolution, and magnitude cuts per pair. A regression case explicitly verifies that
`|e|=1.2` remains eligible and `|e|=2` is rejected.

## Observed support

| Diagnostic | Count |
|---|---:|
| All exported lens-source pairs | 3,069,544 |
| `source_zbin=0` | 814,898 |
| `source_zbin=1..4` | 2,254,646 |
| `|e|>1` but `<2` | 88,460 |
| `|e|>=2` | 0 |
| `sigma_e` outside `[0, 0.4]` | 39,828 |
| Shape cuts + `source_zbin=1..4`, before B-mode/radius cuts | 2,234,510 |

After all four published shape cuts and the HSC Y3 tomographic-bin requirement,
the pair counts are 604,274, 712,085, 590,222, and 327,929 for bins 1–4. These
remain cone-query support counts: radial selection and B-mode handling have not been
applied. The rows have no invalid sky coordinates, all exported weights are positive,
and all z-bin assignments are integers from 0 through 4.

The raw `b_mode_mask` counts are 105,692 zeros and 2,963,852 ones. This frequency
imbalance does not establish which stored value is masked; polarity remains undecoded.
No B-mode-filtered count is reported.

## Reproducibility and status

The quality output is
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/wl_shape_row_quality_2026-10-19_03/summary.json`.
It records hashes for the export plan, row checkpoint, row summary, count summary,
the audit script, and each of the 32 source-row batches. It reports
`complete_descriptive`, `science_measurement=false`, and
`credentials_saved=false`. The 17 S19A adapter tests pass after the selection-rule
correction.

The next calibrated stack still needs the official HSC Y3 SACC `n(z)` input,
release-verified B-mode polarity, redshift-dilution uncertainty, matched random-lens
subtraction, and a paired spatial-block covariance that accounts for source overlap.
No calibrated lensing profile has been calculated from these counts.
