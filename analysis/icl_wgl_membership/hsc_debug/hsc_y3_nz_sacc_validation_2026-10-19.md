# HSC Y3 SACC n(z) input validation — 2026-10-19

The user-approved retrieval was limited to the official HSC Y3 SACC file. No
additional images, shape rows, or per-object photo-z PDFs were downloaded.

## Retrieved artifact

- Official HSC resource: `hsc_y3_fourier_space_data_vector.sacc`
- Local path: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/hsc_y3_nz_sacc_2026-10-19_01/`
- Size: 22,340,160 bytes
- SHA-256: `a28f9e2e088e92d96d2d87083a6eeeac4c0c84c6b48e033bd3d1dc1bf0d8958f`
- Retrieval manifest records the pinned HSC URL, final URL, UTC time, byte count,
  checksum, and `credential_values_saved=false`.

The downloaded file is a valid FITS/SACC artifact and the local loader opened it
successfully using the Astropy fallback (`sacc` is not installed on the current
analysis host). It yielded exactly four ordered tracers, `wl_0` through `wl_3`,
each with 161 samples on the same strictly increasing redshift grid from 0 to 4.
All values are finite and nonnegative. Each unnormalized table integrates to
0.025 on this grid; the stack adapter normalizes each distribution when
computing lensing efficiency. The file contains two copies of each tracer table;
the duplicate z and n(z) arrays match exactly for all four tracers. The loader
deduplicates exact copies and rejects conflicting duplicates.

## Source-bin correspondence

The official HSC S19A shape-catalog documentation describes `hsc_y3_zbin` as
four tomographic source bins numbered 1–4, with bin 0 excluded, and identifies
the SACC `wl_0`–`wl_3` distributions as the corresponding four source-bin n(z)
tracers. Thus the analysis adapter maps source bin `k` to SACC row `wl_{k-1}`.
The reference is the [official HSC S19A shape-catalog page](https://hsc-release.mtk.nao.ac.jp/doc/index.php/s19a-shape-catalog-pdr3/).

The saved 800-lens shape-row QA summary reports:

| HSC source bin | Raw exported rows | Published shape cuts + bin, before B-mode/radius selection |
|---:|---:|---:|
| 0 (unassigned/excluded) | 814,898 | excluded |
| 1 → `wl_0` | 604,638 | 604,274 |
| 2 → `wl_1` | 717,886 | 712,085 |
| 3 → `wl_2` | 598,433 | 590,222 |
| 4 → `wl_3` | 333,689 | 327,929 |
| **Bins 1–4 total** | **2,254,646** | **2,234,510** |

These counts are consistent with the documented four-bin mapping and the
SACC-backed adapter. They are a catalog-definition/count cross-check; no
calibrated shear profile was calculated in this validation step.

## Interpretation and remaining gates

The official HSC page notes that point/photo-z-based lensing estimates require
dilution and selection-bias corrections. The current adapter can integrate the
fiducial bin n(z), but this file does not provide the uncertainty propagation
needed for a final measurement. `science_measurement=false` remains in effect.
Before interpreting a stack, still resolve the stored B-mode mask polarity,
random-lens subtraction and null tests, source reuse/overlap, spatial covariance,
and calibration uncertainty. The raw row export remains unchanged.
