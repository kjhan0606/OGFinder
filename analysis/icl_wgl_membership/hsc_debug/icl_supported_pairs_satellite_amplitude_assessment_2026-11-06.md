# Satellite-injection amplitude response on the 137-pair ICL sample — 2026-11-06

## Design and validation

The satellite injection amplitude grid was extended from the original 57-frame
subset to all 137 frames selected by the measured 200–300 kpc ICL coverage
threshold (44 `i`, 93 `i2`). Three runs use 0.5, 1.0, and 2.0 nJy arcsec⁻² per
satellite at 50 kpc; the diffuse component stays at 1.0 nJy arcsec⁻² at 100 kpc.
SEP detections and masks are regenerated for the injected scenarios.

Each run completed 137/137 frames with zero errors. The position-dependent HSC
PSF-bank auditor verified every frame and file-role hash for all three runs.
The paired grid analysis confirms identical frame IDs, cluster images, local
PSF contents, and no-injection baseline detections at the three amplitudes.
Missing satellite-profile measurements are excluded from a metric's valid
summary and from its paired cluster bootstrap; they are never filled with zero.

## Results

| Filter | Input amplitude | Satellite-profile unit response median (q16–q84) | Valid frames |
|---|---:|---:|---:|
| `i` | 0.5 | 0.285 (−0.374–1.462) | 15 |
| `i` | 1.0 | 0.349 (−0.306–0.955) | 14 |
| `i` | 2.0 | 0.564 (0.300–1.396) | 12 |
| `i2` | 0.5 | 0.612 (−0.418–1.232) | 44 |
| `i2` | 1.0 | 0.409 (−0.384–1.263) | 44 |
| `i2` | 2.0 | 0.411 (−0.404–1.032) | 43 |

| Filter | ICL annulus response to satellite envelopes, amplitude 0.5 / 1 / 2 | Paired median Δ (95% cluster-bootstrap interval), 0.5−1 / 2−1 |
|---|---|---|
| `i` | 0.347 / 0.342 / 0.284 | −0.008 [−0.041, 0.037] / −0.026 [−0.077, −0.003] |
| `i2` | 0.451 / 0.449 / 0.437 | +0.007 [−0.010, 0.038] / −0.007 [−0.023, 0.014] |

Median added masked pixels per frame grow as injection amplitude rises:

| Filter | Injection | Satellite-only scenario | ICL + satellite scenario |
|---|---:|---:|---:|
| `i` | 0.5 | 349 | 363 |
| `i` | 1.0 | 942 | 842 |
| `i` | 2.0 | 1,423 | 1,327 |
| `i2` | 0.5 | 151 | 318 |
| `i2` | 1.0 | 542 | 445 |
| `i2` | 2.0 | 871 | 1,104 |

At least one satellite-profile support gap remains in 78/137, 79/137, and
82/137 frames for amplitudes 0.5, 1.0, and 2.0, respectively. The gap counts
are scenario-level profile-support outcomes, not failed frame executions.

## Interpretation

The larger injected galaxies lead the regenerated source masks to cover more
pixels, so the measurement pipeline's selection footprint is amplitude
dependent. The satellite-profile unit-response distributions remain very broad.
For `i2`, paired median response shifts relative to amplitude 1 are +0.070
with interval [−0.239, 0.213] at amplitude 0.5, and −0.232
[−0.345, 0.106] at amplitude 2. For `i`, only 14 and 12 frames contribute to
those comparisons. The ICL-annulus contribution from satellite envelopes is
fairly stable in `i2`; the `i` 2-versus-1 paired interval narrowly excludes
zero, one of several exploratory metrics in a smaller 44-frame sample.

These runs do not support a universal completeness correction. They show that
mask growth is material while sample-level recovery is imprecise and affected
by substantial satellite-support gaps. They remain selected-sample pipeline
sensitivities, not an astrophysical detection or an ICL–satellite covariance
measurement. No injected-source Poisson variance, HSC upstream sky transfer,
extended PSF-wing calibration, or validated flux covariance is included.
Keep `science_measurement=false`.

## Reproduction and provenance

- Analysis tool: [`analyze_icl_satellite_injection_amplitude_grid.py`](analyze_icl_satellite_injection_amplitude_grid.py), SHA-256 `2291744adf04d752c3f9ffa957541d2ed78c9ca6196ca28925a3928e11529cbb`.
- Base 1.0 run: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/joint_icl_satellite_psf_redetection_icl_supported224_2026-11-06_01/summary.json`, SHA-256 `be3c2e8cb22185c06768d9022ae3e6192072f180fa54d55e21b85920e5572750`.
- 0.5 run: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/joint_icl_satellite_psf_redetection_icl_supported224_amp0p5_2026-11-06_01/summary.json`, SHA-256 `c61584ba81f7a21305c50f3bb3dd3738a8c6fc48a0bf5b593cd3aad86e99c656`; audit SHA-256 `2b2bfe36fb5348a5aebfa3aa86c8b8556e146943a2cf6a7bb0a47d0c444efb9e`.
- 2.0 run: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/joint_icl_satellite_psf_redetection_icl_supported224_amp2p0_2026-11-06_01/summary.json`, SHA-256 `55454c092f7fcace28d55bacab7f17d523e294d3ec4fbfef04ea75745218fbcb`; audit SHA-256 `e8c7cd409c0100d18a397cbfaca105f1d2bce4284c7f783316bca0cde34f094e`.
- Paired grid summary: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_satellite_injection_amplitude_grid_supported137_2026-11-06_01/summary.json`, SHA-256 `1191414120941dc2022001ed260f277601a4af038367482383f7226cf55707f3`.
- Paired frame table SHA-256 `3b0c2432820d412b7ac3895bca5dadcccc810ee97f8efe055fca71736e773198`; machine-generated assessment SHA-256 `c16b96ecf2cff7fcca21b8255293d40bc795b166db3ab191f2e3034173091851`.

The analysis supports both the original 57-frame input and an explicitly
selected larger paired set via `--expected-frames`; the default remains 57.
