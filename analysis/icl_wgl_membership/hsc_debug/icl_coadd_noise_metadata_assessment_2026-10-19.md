# ICL observed-coadd noise metadata audit — 2026-10-19

## Scope and reproducibility

The ICL/satellite injection pipeline adds synthetic components to copies of
the actual downloaded HSC `coadd/bg` science arrays, reruns source detection and
masks, and measures the result against each image's delivered variance plane.
This audit checks whether the target coadd headers contain enough calibration
metadata to add the source-photon variance consistently. It reads FITS headers
only and does not modify or download images.

Astrafek's reproducible tool is
[`audit_icl_coadd_noise_metadata.py`](audit_icl_coadd_noise_metadata.py).
It follows the 57 cluster/filter frames in the checksum-audited
position-dependent PSF injection run and verifies their image paths against
the saved mask records and each frame's recorded image hash. The output is
`/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_coadd_noise_metadata_2026-10-19_01/summary.json`
(SHA-256 `c210ee9dd604e07ebe820883faab040d6d3f7e762a8724e6111c7244d6e37c7b`).
Its input PSF injection audit is `joint_icl_satellite_psf_redetection_2026-10-19_02/psf_redetection_audit.json`
(SHA-256 `ce377060cfa2143936b36a400f5bfdb33fb355956bf3513d7a2fa4964595ac1e`).

## Findings

| FITS metadata | Present frames | Missing frames |
|---|---:|---:|
| `FLUXMAG0` | 57/57 | 0 |
| science, mask, and variance image planes with matching dimensions | 57/57 | 0 |
| `GAIN` | 0/57 | 57 |
| `EXPTIME` | 0/57 | 57 |
| `DARKTIME`, `BUNIT`, `SKYSIGMA`, `NCOMBINE` | 0/57 each | 57 each |

All audited frames share `FLUXMAG0=63095734448.01944`, which establishes the
photometric conversion used by the surface-brightness measurements. The
variance plane supplies the delivered per-pixel variance. However, the headers
do not provide gain or exposure metadata for computing the additional photon
variance associated with synthetic source light. No gain was guessed, and the
existing variance plane was left unchanged during the injection tests.

HSC describes its coadds as resampled and averaged CCD images and distinguishes
global-sky coadds from local-sky-subtracted measurement images. The available
`coadd/bg` frames remain downstream of HSC's exposure-level and sky processing;
an injection into them cannot measure the full upstream transfer function. See
the [HSC PDR3 processing description](https://hsc-release.mtk.nao.ac.jp/doc/index.php/processing__pdr3/).

## Interpretation and next gate

The earlier reports' statement that photon noise was not added remains
correct, but the injections do operate on observed coadd pixels rather than
isolated synthetic arrays. The missing input is a defensible source-shot-noise
model, not access to the observed science arrays. To add that term, use
exposure-level calibrated inputs with their gain/weight metadata or establish
and validate the effective source-variance coefficient from the coadd
construction. Correlated coadd noise and upstream sky-subtraction transfer
remain separate systematics. The current injection and ICL–satellite results
remain exploratory: `science_measurement=false`,
`observational_validation_passed=false`, and
`flux_covariance_available=false`.
