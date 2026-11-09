# S23B CAMIRA member outer-profile injection recovery

## Scope and method

This CPU-only exploratory measurement-transfer check uses the existing S23B
CAMIRA member outer-profile sample. It selects one seeded-random member that
already has a supported 50–75 kpc profile in each target/filter frame, then
adds a circular exponential halo to the target coadd at 50 kpc with
amplitudes 0.5, 1, and 2 nJy arcsec⁻² and scale length 15 kpc. For every
amplitude it repeats SEP detection, one-to-one member matching, source masks,
and the same-radius rotated-control measurement. No injected images are
written.

The recovered response is the injected-minus-baseline target-minus-control
profile increment divided by the nominal, unmasked model mean in 50–75 kpc.
Intervals are 10,000-draw cluster-bootstrap intervals for the median. A
separate model-through-mask calculation records the effect of the revised
segmentation masks and rotated controls.

## Results

All 427 selected frames reproduced their saved baseline profile exactly
(maximum absolute contrast residual 0). There are 62 HSC `i` frames and 365
HSC `i2` frames, giving 1,281 frame-amplitude rows. Of these, 1,221 have
adequate target and control support; the other 60 fail because of insufficient
target annulus coverage (30), too few rotated controls (16), or failure to
recover a unique member match (14).

| Filter | Injected SB at 50 kpc (nJy arcsec⁻²) | Supported frames | Median recovered / input | Cluster-bootstrap 95% interval | Median mask-model transfer |
|---|---:|---:|---:|---:|---:|
| HSC `i` | 0.5 | 61 | 0.896 | [0.733, 1.043] | 0.929 |
| HSC `i` | 1.0 | 59 | 0.838 | [0.597, 0.974] | 0.921 |
| HSC `i` | 2.0 | 52 | 0.879 | [0.706, 1.035] | 0.914 |
| HSC `i2` | 0.5 | 355 | 0.808 | [0.650, 0.893] | 0.777 |
| HSC `i2` | 1.0 | 352 | 0.807 | [0.707, 0.893] | 0.803 |
| HSC `i2` | 2.0 | 342 | 0.819 | [0.734, 0.871] | 0.818 |

Paired response changes across amplitudes are small compared with their
cluster-bootstrap intervals: for `i`, the median differences are +0.001
(0.5−1) and +0.088 (2−1), with both intervals including zero. For `i2`, they
are +0.015 and +0.056, also with both intervals including zero. This supports
approximate amplitude stability over this injected grid, conditional on the
selected frames and measurement pipeline; it does not establish unbiased
recovery.

The source-selection masks respond substantially to the injection. Detection
counts change in 62%, 64%, and 85% of supported `i` frames and 57%, 68%, and
74% of supported `i2` frames at increasing amplitudes. Median segmentation
footprint increases by 17, 39, and 92 pixels in `i`, and 12, 30, and 65.5
pixels in `i2`. The median response remains below unity, especially in `i2`;
mask and control transfer therefore matter even though the paired amplitude
differences are not resolved from zero.

## Interpretation and limits

This result calibrates the pipeline response to a specific smooth, circular
model around already-detected CAMIRA members. It does not measure candidate
completeness, establish a physical satellite-halo profile, or demonstrate an
ICL–satellite association. The member sample and 50–75 kpc profiles are
conditional on prior detection, coverage, and rotated-control support.

The injected halo is not convolved with a position-dependent PSF. No source
Poisson noise is added and the delivered VARIANCE plane is held fixed. Only the
target image is injected; rotated controls remain on the original image and
can contain cluster light or subtract some injected wings. Consequently these
intervals omit important image and selection uncertainties. Keep
`science_measurement=false` and `images_ready_for_science=false`; do not use
the ratios as a correction to the science profile.

## Reproducibility

- Runner: `run_s23b_member_outer_profile_injection_recovery.py`
- Runner SHA-256: `ef4835b4491c1aa318d14f843eb95d9801c8d2238a9ca8a5842e1292a058f0d9`
- Result directory: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/s23b_camira_member_outer_injection_recovery_2026-11-07_01/`
- Summary SHA-256: `e9d4cced59f0b0f8ea3f189989886f28af74a339cc5a9b3421f8bca187defc77`
- Response table SHA-256: `b7cc4deeff36dc95dd0fcbd41cef929b49e6524761f7b561586364d4288bfa30`
- Plot: `member_halo_injection_recovery.png` in the result directory
- Network access and data downloads: none
