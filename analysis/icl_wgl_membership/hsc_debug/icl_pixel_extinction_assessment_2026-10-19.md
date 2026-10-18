# Pixel-level SFD correction of the HSC ICL annulus (2026-10-19)

## Scope and reproducibility

This audit asks whether evaluating Milky-Way extinction at each usable image pixel changes the existing 200–300 kpc HSC-I target-minus-control profile. It uses the already downloaded SFD north/south maps, the existing 800 target images and masks, and the 798 extracted random-control profiles. It makes no network requests and does not alter source images or masks.

The analysis applies the saved filter coefficient to each pixel as `10**(0.4 * A_lambda)`, where `A_lambda = coefficient * E(B-V)`. A 128×128 WCS grid is interpolated to the image plane; 512 random image pixels per frame are compared with direct native-map WCS interpolation. The largest observed absolute E(B−V) interpolation difference across the 448 selected frames was `3.05e-5`. The native SFD map resolution remains the limiting angular resolution.

The script re-measures the fixed K2.5 / common coverage ≥40% / 200–300 kpc selection, verifies every image, mask, control-baseline, dust-map and table hash, and requires the raw target and control profiles and coverages to reproduce the saved values. It also computes a separate sensitivity after subtracting the local plane fit from the 600–900 kpc sky annulus before applying extinction. That plane-subtracted path is diagnostic; it is not assumed to be a superior estimate because it can remove real extended light.

Final output:

`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_pixel_extinction_2026-10-19_03/`

The 224-row profile table and 74-row association table have verified SHA-256 hashes, unique cluster/filter keys, and reproduce the original raw target/control table to the set 2×10⁻⁶ nJy arcsec⁻² tolerance. The recorded source-hash map contains 912 input entries. Output status is `complete_exploratory`; `science_measurement=false`.

## ICL profile result

| Filter | Pairs | Raw mean ΔSB | Pixel-SFD corrected mean ΔSB | Center-SFD corrected mean ΔSB | Pixel minus center |
|---|---:|---:|---:|---:|---:|
| i | 71 | 3.226 | 3.379 `[2.189, 4.484]` | 3.380 | −0.00049 |
| i2 | 153 | 3.787 | 4.000 `[2.836, 5.133]` | 4.001 | −0.00014 |

Units are nJy arcsec⁻². Intervals are central 95% percentile intervals from 5,000 resamples of connected joint target/control 5-degree sky blocks (18 i blocks and 40 i2 blocks). The per-pixel extinction correction is effectively indistinguishable from applying the field-center correction to the raw target-minus-control profile. This confirms that the earlier center-based correction is adequate for this annulus at the current SFD resolution.

The plane-first pixel-corrected means were 3.410 `[1.900, 4.610]` for i and 4.521 `[3.417, 5.648]` for i2. These differ from the raw-profile correction by +0.031 and +0.521 nJy arcsec⁻², respectively. This is sensitivity to the local sky-plane treatment, not evidence that plane subtraction gives the true ICL profile.

## Satellite association sensitivity

The association reuses the fixed 74 cluster/filter pairs (i 21, i2 53) and the existing satellite 50–75 kpc profile. It residualizes both measurements against redshift and CAMIRA richness, then resamples connected 5-degree target/control blocks 5,000 times. The satellite profile receives the field-center SFD factor; only the cluster ICL annulus receives the pixel-by-pixel correction.

| Filter | Raw ICL r (95% block interval) | Pixel-SFD raw-profile r (95% interval) | Plane-first pixel-SFD ICL r (95% interval) |
|---|---:|---:|---:|
| i | 0.116 `[-0.225, 0.518]` | 0.112 `[-0.224, 0.515]` | 0.114 `[-0.202, 0.579]` |
| i2 | 0.235 `[0.013, 0.422]` | 0.239 `[0.013, 0.427]` | 0.124 `[-0.165, 0.373]` |

The ICL association barely changes under pixel-level dust correction. The i2 point estimate and interval change materially when a local plane is removed, and the interval then includes zero. Thus the apparent i2 association is more sensitive to sky-profile handling than to the spatial variation of foreground extinction. Neither path is confirmatory: the sample was selected before these corrections, the relation is exploratory, and the prior familywise audit found no group surviving multiplicity correction.

## Remaining limits and next step

This analysis holds source masks and the pair selection fixed. The 50–75 kpc satellite flux is not corrected pixel-by-pixel, and the target-minus-control data still lack validated upstream sky transfer, PSF-wing subtraction, mask-leakage calibration, and a full joint ICL/satellite covariance. The plane-first result is not adopted as a replacement profile. The next ICL gate is therefore a pre-registered comparison of sky estimators with realistic diffuse-light injections and matched random fields, tracking how each estimator transmits a known 200–300 kpc signal and how it changes the satellite association.

The weak-lensing source rows are already exported, but the calibrated HSC Y3 source `n(z)` input and documented `b_mode_mask` polarity are still required before interpreting a lensing stack. The HSC documentation describes `b_mode_mask` as the cut for large-B-mode regions without defining the stored 0/1 polarity, so the raw bit remains undecoded ([official S19A shape-catalog page](https://hsc-release.mtk.nao.ac.jp/doc/index.php/s19a-shape-catalog-pdr3/)).
