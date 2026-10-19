# HSC PDR3 PSF Picker single-position test — 2026-10-19

## Scope and provenance

After user approval, one coadd-model PSF was requested for CAMIRA cluster 10
(HSC-I2 image center; RA 0.339748 deg, Dec +1.030287 deg). The HSC PDR3 PSF
Picker utility was used with `rerun=pdr3_wide`, `type=coadd`, and
`filter=HSC-I`. PDR3 coaddition combines the `i` and `i2` exposure variants;
the Picker's supported broad-band name is `HSC-I`.

The official HSC utility source was fetched from
[`ssp-software/data-access-tools/pdr3/downloadPsf`](https://hsc-gitlab.mtk.nao.ac.jp/ssp-software/data-access-tools/-/tree/master/pdr3/downloadPsf/).
The exact downloaded source had SHA-256
`17cee33a86de87a77f3f800049b3fdb056f2e0c50b96c2b8cf7e80c1f62f2298`.
The HSC documentation describes PSF Picker as a tool for downloading model
PSF images and states that coordinate-list requests return FITS models at the
requested positions ([data access](https://hsc-release.mtk.nao.ac.jp/doc/index.php/data-access__pdr3/),
[PDR3 FAQ](https://hsc-release.mtk.nao.ac.jp/doc/index.php/faq__pdr3/)).

## Received file and checks

- File: `/gpfs/kjhan/Astrafek/data/sample_800_2026-10-13/icl_psf_picker_test_2026-10-19_01/cluster10_i2_center_psf.fits`
- Size: 20,160 bytes; SHA-256: `1276d46113d389116ded6f63bae4c9865bb42d24244b22247a84d2c6929cb202`
- FITS: one primary image, 43 × 43 pixels, float64; all values finite.
- Pixel sum: `0.9999999999999993`; minimum `-5.045225748178976e-05`, maximum `0.06360799166085342`; 498 pixels are negative, so the original model values must be preserved during convolution.
- Positive-flux centroid: `(21.0273, 21.0122)` pixels, close to the geometric center `(21, 21)`; positive-flux radial second-moment scale: 2.72 pixels.
- The FITS header contains no sky coordinate, filter, tract, patch, or rerun metadata. The request coordinates and HSC product identity are therefore recorded here and must be retained alongside the file. Future batch requests should use an output-name template that includes the returned tract, patch, and filter metadata.

## Interpretation and next gate

At 0.168 arcsec/pixel, the 43-pixel stamp spans about 7.2 arcsec (3.6 arcsec
radius). It characterizes the local PSF core, not the extended wings needed to
validate ICL profiles over tens to hundreds of kpc. The already-downloaded
4001 × 4001 extended-PSF reference remains necessary for wing sensitivity;
this single model PSF does not establish a target-specific extended PSF.

The current redetection-injection sample contains 57 cluster/filter frames,
with 688 satellite injection positions. A full position-dependent core-PSF
set would therefore contain up to 745 position/filter requests (57 frame
centers plus 688 satellite centers), including the one test position above.
At the measured test-file size, the remaining 744 models would be about
15.0 MB before archive/metadata overhead. This estimate is provisional because
only one location has been tested. No additional PSFs, images, shape rows, or
photo-z values were requested in this run; the remaining PSF set requires a
separate user approval.

`science_measurement=false` remains in force. The result is a PSF-input
feasibility test, not an ICL or weak-lensing measurement.
