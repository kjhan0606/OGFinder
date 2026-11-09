# Satellite-injection amplitude response assessment

## Scope

Ran the same position-dependent PSF and detection/mask-reselection pipeline on
the 57-frame PSF-supported subset (13 `i`, 44 `i2`) at satellite surface
brightness amplitudes 0.5, 1, and 2 nJy arcsec⁻² per galaxy at 50 kpc. The
diffuse ICL-like injection remained 1 nJy arcsec⁻² at 100 kpc. Native images,
PSF file content, satellite profile shape, and all no-injection baseline masks
were identical across the three runs. Detection and masks were regenerated
for the injected scenarios.

## Results

| Filter | Satellite amplitude | Unit satellite recovery median (q16–q84) | Paired median change from 1.0 (95% cluster-bootstrap interval) |
|---|---:|---:|---:|
| i | 0.5 | −0.273 (−3.146, 0.641) | −0.444 (−0.991, 0.310) |
| i | 1.0 | −0.144 (−1.807, 0.348) | reference |
| i | 2.0 | 0.412 (−0.989, 0.876) | +0.533 (−0.526, 1.810) |
| i2 | 0.5 | 0.593 (−0.258, 1.446) | −0.006 (−0.401, 0.162) |
| i2 | 1.0 | 0.572 (−0.516, 1.311) | reference |
| i2 | 2.0 | 0.560 (−0.087, 1.023) | +0.000 (−0.221, 0.230) |

The unit response is normalized to input amplitude. For `i2` it is stable
across this 4× amplitude range at the ensemble median; both paired intervals
contain zero. Individual-frame spread is large, and the 13-frame `i` result is
especially imprecise. Median 50–75 kpc detections/masks do change with
amplitude: the number of frames whose satellite-only detection count differs
from baseline is 41/44, 41/44, and 41/44 for `i2` at amplitudes 0.5, 1, and 2;
the median increase in masked pixels is about 230, 631, and 1,014. Thus the
stable median transfer coefficient does not mean selection masks are
amplitude-independent.

The ICL-aperture response to satellite envelopes is also broadly stable in
`i2` (median unit response 0.480, 0.490, 0.485 across the three amplitudes).
In `i`, the 2.0-versus-1.0 paired interval is narrowly below zero
(`−0.025`, 95% interval `[-0.088,-0.0003]`), but it is one of several
exploratory metrics from only 13 frames and is not treated as evidence of a
physical effect. All component files retain `science_measurement=false`.

## Interpretation

This check supports approximate amplitude scaling of the median `i2`
satellite-profile transfer over the simulated range. The mask footprint
changes with brightness, while frame-level recovery has broad scatter. It is a
pipeline transfer and selection diagnostic on a PSF-selected subset; it does
not estimate the real satellite luminosity distribution, completeness, or an
ICL–satellite detection. HSC coadd correlated noise, the physical Poisson
variance of the injected light, and extended PSF wings are not calibrated by
these runs.

## Reproducibility

- Analysis script: [`analyze_icl_satellite_injection_amplitude_grid.py`](analyze_icl_satellite_injection_amplitude_grid.py).
- Analysis script SHA-256 `fa64fd9d554e9b44b8306859bc1ce3a1a2a7679db84898fa8eda079498561852`.
- Output directory: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_injection_amplitude_grid_2026-11-05_02`.
- Paired table SHA-256 `dbd5489eb3c799bb6a7fb53f9251a800923dd2fd0cf30f679589a5b421cab5e5`.
- Summary SHA-256 `733bcb90a7471be061b0280259e6cc9bc87fb05f1761f500d25cd91727d22d63`.
- Assessment SHA-256 `be79677ea2a43308c99c2b77e8055d524c9a6ee21d4a9f4cebdb60abe0c78a59`.
- Paired cluster-ID bootstrap: 10,000 resamples, fixed seed `20261105`.
