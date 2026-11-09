# Injection-amplitude sensitivity of the ICL–satellite transfer matrix (2026-11-06)

## Scope and matched-input audit

Astrafek compared the existing PSF-convolved, redetection-enabled batches at
satellite amplitudes 0.5, 1.0, and 2.0 nJy arcsec⁻² per member at 50 kpc. Each
batch contains the same 137 target/filter frames. The image and ancillary input
hashes match across amplitudes, and the no-injection baseline redetection is
identical. Diffuse injection remains 1.0 nJy arcsec⁻² at 100 kpc. Source
detection and masks are regenerated for each injected scenario.

The 1.0-amplitude response coefficients reproduce the corrected 58-pair
cross-talk table exactly. The 0.5-amplitude batch has valid satellite response
for 57 of the 58 observed association pairs. The 2.0-amplitude batch has valid
response for 54 pairs. Four pairs lose support at the brighter amplitude, and
one pair loses support at the fainter amplitude. The strict common support
therefore contains 54 pairs, 12 in `i` and 42 in `i2`, across eight and 22
connected 5-degree blocks.

## Pairwise matrix-inversion sensitivity

| Filter | Injection amplitude | Observed partial *r* | Inverted partial *r* | Conditional 95% block interval | Median / maximum matrix condition number |
|---|---:|---:|---:|---:|---:|
| `i` | 0.5 | 0.483 | −0.443 | [−0.979, 0.771] | 4.54 / 43.55 |
| `i` | 1.0 | 0.483 | −0.946 | [−0.999, 0.943] | 4.10 / 16.32 |
| `i` | 2.0 | 0.483 | −0.841 | [−0.994, 0.563] | 4.10 / 11.26 |
| `i2` | 0.5 | 0.417 | −0.922 | [−0.999, −0.318] | 6.45 / 2,684.74 |
| `i2` | 1.0 | 0.417 | −0.381 | [−0.971, −0.285] | 6.79 / 126.72 |
| `i2` | 2.0 | 0.417 | −0.075 | [−0.745, −0.032] | 6.52 / 962.43 |

The medians of the i2 satellite-to-ICL response terms change little across
amplitudes. The pairwise matrix inversion changes substantially. The inferred
i2 partial correlation shifts from −0.922 at amplitude 0.5 to −0.075 at
amplitude 2.0. The per-pair matrix condition numbers reach 2,685 and 962 at
those amplitudes. Among 42 i2 pairs, 16–20 inferred diffuse amplitudes and
16–20 inferred satellite amplitudes are negative across the three amplitudes.
The maximum absolute inferred i2 diffuse amplitude ranges from 1,840 to
16,127 nJy arcsec⁻². The maximum absolute satellite amplitude ranges from 353
to 2,212 nJy arcsec⁻².

## Interpretation

The amplitude sweep shows that stable ensemble-median response coefficients do
not ensure stable pairwise matrix inversion. Redetection and masking alter
pair-level transfer matrices, and a small number of nearly singular matrices
produce extreme inferred amplitudes. The connected-block intervals condition
on the measured matrices. The intervals omit injection-response uncertainty,
measurement covariance, source-Poisson noise, coadd weights, extended PSF wings,
and upstream sky transfer. The support also excludes five of the original
58-pair overlap at one or more amplitudes.

The negative i2 intervals are not physical anti-correlation detections. The
zero-intercept linear inversion fails the physical-amplitude and stability
checks. A forward model must propagate injection-response uncertainty,
measurement covariance, and support selection without unstable pairwise
division. Until the forward model and image-level noise gates pass,
`science_measurement=false` remains.

## Reproducibility

- Astrafek runner: `analyze_icl_satellite_amplitude_grid_transfer_matrix.py`
- Runner SHA-256: `aefebb159ea64c062a4a8e9951eac80ee70b3d9134d2a51ff9f9baa7f46fc891`
- Output directory: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_amplitude_grid_transfer_matrix_2026-11-06_01/`
- Summary SHA-256: `d7f3dd636b8efe3201c0bba9f79e690149722567b37f2b7c01e1c839ee6013e3`
- Output-hash file SHA-256: `245df97f723a497b0e286ca1fcbc4f9173b466308ff19e77a0d4a15aed1bb1e1`
- Input batch-summary hashes: 0.5 `c61584ba81f7a21305c50f3bb3dd3738a8c6fc48a0bf5b593cd3aad86e99c656`, 1.0 `be3c2e8cb22185c06768d9022ae3e6192072f180fa54d55e21b85920e5572750`, and 2.0 `55454c092f7fcace28d55bacab7f17d523e294d3ec4fbfef04ea75745218fbcb`
- HSC queries and downloads: none
